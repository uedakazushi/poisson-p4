#!/usr/bin/env python3
"""A universal quadratic polynomial section of S61 -> W16.

The degree-two part of the q52 constant-kernel condition is a constant
50x45 linear system. Its unique solution has186 quadratic monomials; all134
Poisson equations vanish identically after substitution. This realizes every
first jet but does not exhaust the other S61 components.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
from s61_model import S61,X,coefficients,bivector_from_coefficients,J as J5
from elliptic_weight1_dominance import projective_affine
from weight1_structure import J
ROOT=Path(__file__).resolve().parents[1]


def derive_section(model):
    weights=model.reduced['parameter_weights'];one=[i for i,w in enumerate(weights) if w==1];two=[i for i,w in enumerate(weights) if w==2]
    u=s.symbols('u1:62');raw=coefficients(X[0]**2*J5)
    for i,v in zip(model.certificate['raw_free_indices'],u,strict=True):raw[i]=v
    # Only parameter weights0,1,2 contribute to the affine two-jet. Omitting
    # weight3 avoids constructing the cubic tail before it is needed.
    for entry in model.data['eliminated_raw_variables']:
        degree={sum(weights[i] for i in mon) for mon,_,_ in entry['polynomial']}
        assert len(degree)<=1
        raw[entry['raw_index']]=model.evaluate(entry['polynomial'],u) if not degree or max(degree)<=2 else 0
    y=s.symbols('y0:4');p=projective_affine(bivector_from_coefficients(raw),X,y)
    assert all(s.Poly(p[i,j],*y).total_degree()<=2 for i in range(4) for j in range(4))
    zero=dict.fromkeys(y,0);linear=sum((p.diff(yi).subs(zero)*yi for yi in y),s.zeros(4))
    beta=(-J*linear*J).applyfunc(s.expand)
    n=s.Matrix([[0,p[2,3],-p[1,3],p[1,2]],[-p[2,3],0,p[0,3],-p[0,2]],
                [p[1,3],-p[0,3],0,p[0,1]],[-p[1,2],p[0,2],-p[0,1],0]])
    r=s.symbols('r0:5');z=s.Matrix([[0,r[0],r[1],r[2]],[-r[0],0,r[3],r[4]],
                                  [-r[1],-r[3],0,-r[0]],[-r[2],-r[4],r[0],0]])
    ell=s.Matrix([sum(s.diff(beta[i,j],y[k])*z[i,j] for i in range(4) for j in range(i+1,4)) for k in range(4)])
    w=-J*ell;yy=s.Matrix(y);kernel=z+yy*w.T-w*yy.T
    residual=s.Poly(s.expand(sum(n[i,j]*kernel[i,j] for i in range(4) for j in range(i+1,4))),*y,*r)
    selected=[(mon,c) for mon,c in residual.terms() if sum(mon[:4])==2]
    equations=[c for _,c in selected];q=[u[i] for i in two];a,b=s.linear_eq_to_matrix(equations,q)
    assert a.shape==(50,45) and not a.free_symbols
    rows=list(a.T.rref()[1]);assert len(rows)==45
    solution=(a[rows,:].inv()*b[rows,0]).applyfunc(s.expand)
    assert (a*solution-b).applyfunc(s.expand)==s.zeros(50,1)
    variables=[u[i] for i in one]
    assert all(s.Poly(f,*variables).total_degree()==2 for f in solution if f)
    assert sum(len(s.Poly(f,*variables).terms()) for f in solution)==186
    return u,one,two,a,b,rows,solution,[list(mon) for mon,_ in selected]


def verify_equations(model,u,one,two,solution):
    variables=[u[i] for i in one];base,*x=ring(','.join(map(str,variables)),QQ)
    values=[None]*61
    for i,v in zip(one,x,strict=True):values[i]=v
    for i,f in zip(two,solution,strict=True):values[i]=base.from_expr(f)
    for equation in model.reduced['equations']:
        value=base.zero
        for mon,n,d in equation:
            term=base(QQ(n,d))
            for i in mon:term*=values[i]
            value+=term
        assert not value
    return values



def verify_full_graph(model,u,two,solution):
    values=list(u)
    for i,f in zip(two,solution,strict=True):values[i]=f
    y=s.symbols('y0:4');p=projective_affine(model.reconstruct(values),X,y)
    zero=dict.fromkeys(y,0);linear=sum((p.diff(yi).subs(zero)*yi for yi in y),s.zeros(4));beta=-J*linear*J
    n=s.Matrix([[0,p[2,3],-p[1,3],p[1,2]],[-p[2,3],0,p[0,3],-p[0,2]],
                [p[1,3],-p[0,3],0,p[0,1]],[-p[1,2],p[0,2],-p[0,1],0]])
    r=s.symbols('r0:5');z=s.Matrix([[0,r[0],r[1],r[2]],[-r[0],0,r[3],r[4]],
                                  [-r[1],-r[3],0,-r[0]],[-r[2],-r[4],r[0],0]])
    ell=s.Matrix([sum(s.diff(beta[i,j],y[k])*z[i,j] for i in range(4) for j in range(i+1,4)) for k in range(4)])
    w=-J*ell;yy=s.Matrix(y);kernel=z+yy*w.T-w*yy.T
    assert s.expand(sum(n[i,j]*kernel[i,j] for i in range(4) for j in range(i+1,4)))==0
    affine_degree=max(s.Poly(p[i,j],*y).total_degree() for i in range(4) for j in range(i+1,4))
    assert affine_degree==3
    return affine_degree


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    model=S61();u,one,two,a,b,rows,solution,monomials=derive_section(model)
    verify_equations(model,u,one,two,solution)
    affine_degree=verify_full_graph(model,u,two,solution)
    point_source=ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json'
    record=next(p for p in json.loads(point_source.read_text())['points'] if p['component_id']=='T01')
    point=list(map(s.Rational,record['S61_coordinates_u1_through_u61']))
    sample=list(point);substitution={u[i]:point[i] for i in one}
    for i,f in zip(two,solution,strict=True):sample[i]=f.subs(substitution)
    assert all(model.evaluate(eq,sample)==0 for eq in model.reduced['equations'])
    assert sample!=point
    report={'id':'Q52_UNIVERSAL_SECTION_20261101','field':'QQ[u_weight1]','all_passed':True,
        'weight1_u_indices_zero_based':one,'weight2_u_indices_zero_based':two,
        'constant_linear_matrix':[[str(v) for v in row] for row in a.tolist()],
        'quadratic_rhs':list(map(str,b)),'degree2_y_and_r_monomials':monomials,
        'pivot_rows':rows,'pivot_determinant':str(a[rows,:].det()),
        'quadratic_section_q':list(map(str,solution)),'total_quadratic_monomials':186,
        'universal134_identities_checked':True,'full_graph_annihilation_checked':True,
        'full_affine_bivector_degree':affine_degree,
        'sample_at_saved_toric_beta':list(map(str,sample)),
        'claim':'Every first jet admits an explicit polynomial Poisson lift. The reduced T40 closure is this graph, using its established dominant constant-kernel characterization.',
        'scope':['No assertion that other components or other lifts are exhausted.','No classification of embedded structure at component intersections.'],
        'proof_note':'research/notes/Q52_UNIVERSAL_SECTION_20261101.md',
        'source_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in (Path(__file__),ROOT/'code/s61_model.py',point_source,ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json',ROOT/'legacy/v02/data/symplectic_slice_equations.json')}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'constant_linear_rank':45,'quadratic_terms':186,'every_jet_lifts':True,'output':str(args.output)}))

if __name__=='__main__':main()

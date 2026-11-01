#!/usr/bin/env python3
"""Exact sample isomorphisms from the toric cubic to the canonical E[2] cubic.

Two independent rational first jets, all96 associativity relations, invertible
trace forms, and exact primitive-element maps are verified. These samples
support a structural conjecture; they do not prove a generic identification.
"""
from __future__ import annotations
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from sympy.polys.numberfields import field_isomorphism
from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS
from toric_algebra_fiber import TRIPLES,multiplication_from_jet,associativity_quadrics,trace_form
from elliptic_curve_invariants_20261101 import curve_invariants_from_beta
ROOT=Path(__file__).resolve().parents[1]


def sample(vector,linear,beta_values):
    c=s.symbols('c0:20');t=s.Symbol('t');m=multiplication_from_jet(vector,dict(zip(TRIPLES,c,strict=True)))
    form,equations=associativity_quadrics(m);rref,pivots=linear.rref();assert len(pivots)==18
    free=[j for j in range(20) if j not in pivots];assert free==[18,19];x,y=[c[j] for j in free];solution=list(c)
    for row,index in enumerate(pivots):solution[index]=-sum(rref[row,j]*c[j] for j in free)
    substitution=dict(zip(c,solution,strict=True));residual=[s.expand(f.subs(substitution,simultaneous=True)) for f in equations]
    monomials=[x*x,x*y,y*y];leading=s.Matrix([[s.Poly(f,x,y).coeff_monomial(mon) for mon in monomials] for f in residual])
    rows=list(leading.T.rref()[1]);assert len(rows)==3
    nf=(-leading[rows,:].inv()*(s.Matrix([residual[i] for i in rows])-leading[rows,:]*s.Matrix(monomials))).applyfunc(s.expand)
    assert all(s.Poly(f,x,y).total_degree()<=1 for f in nf)
    conics=[s.expand(mon-f) for mon,f in zip(monomials,nf,strict=True)]
    assert all(s.expand(f-sum(leading[i,j]*conics[j] for j in range(3)))==0 for i,f in enumerate(residual))
    def nv(f):
        polynomial=s.Poly(f,x,y);return s.Matrix([polynomial.coeff_monomial(mon) for mon in (1,x,y)])
    mx=s.Matrix.hstack(s.Matrix([0,1,0]),nv(nf[0]),nv(nf[1]));my=s.Matrix.hstack(s.Matrix([0,0,1]),nv(nf[1]),nv(nf[2]));assert mx*my==my*mx
    unit=s.Matrix([1,0,0]);cyclic=s.Matrix.hstack(unit,my*unit,my**2*unit);assert cyclic.det()
    xp=cyclic.inv()*s.Matrix([0,1,0]);xpoly=sum(xp[i]*t**i for i in range(3));cvalues=[s.expand(f.subs({x:xpoly,y:t})) for f in solution]
    cubic=my.charpoly(t).as_expr();assert s.Poly(cubic,t).is_irreducible
    reduce=lambda f:s.rem(s.expand(f),cubic,t)
    trace=s.Matrix([[reduce(f.subs(dict(zip(c,cvalues,strict=True)),simultaneous=True)) for f in row] for row in form.tolist()])
    determinant=reduce(trace.det(method='domain-ge'));bezout=s.gcdex(determinant,cubic,t);assert bezout[2]==1
    assert s.expand(bezout[0]*determinant+bezout[1]*cubic)==1
    c4,c6,delta,j=curve_invariants_from_beta(vector);assert delta
    e2=t**3-27*c4*t-54*c6;isomorphism=None;root_index=None
    for index in range(3):
        possible=field_isomorphism(s.CRootOf(e2,index),s.CRootOf(cubic,0),fast=True)
        if possible:isomorphism=possible;root_index=index;break
    assert isomorphism
    X=sum(value*t**(len(isomorphism)-1-k) for k,value in enumerate(isomorphism))
    assert reduce(e2.subs(t,X))==0
    change=s.Matrix([[s.Poly(reduce(X**j),t).coeff_monomial(t**i) for j in range(3)] for i in range(3)]);assert change.det()
    # All20 tensor coordinates are residue expressions of the expected weights.
    inverse=change.inv();U=[];V=[]
    for value in cvalues:
        polynomial=s.Poly(reduce(value*(3*X*X-27*c4)),t)
        coefficient=inverse*s.Matrix([polynomial.coeff_monomial(t**i) for i in range(3)])
        assert coefficient[2]==0;V.append(coefficient[0]);U.append(coefficient[1])
    assert linear*s.Matrix(U)==linear*s.Matrix(V)==s.zeros(51,1)
    assert all(reduce(value*(3*X*X-27*c4)-a*X-b)==0 for value,a,b in zip(cvalues,U,V,strict=True))
    # Check the three96-vector identities that a universal C5,C9 formula would
    # need. These are still checked only at the two rational beta values here.
    C=s.Matrix([f.subs(dict.fromkeys(c,0)) for f in equations]);B=s.Matrix(equations).jacobian(c).subs(dict.fromkeys(c,0))
    H=s.Matrix(equations)-B*s.Matrix(c)-C
    hu=H.subs(dict(zip(c,U,strict=True))).applyfunc(s.expand);hv=H.subs(dict(zip(c,V,strict=True))).applyfunc(s.expand)
    huv=(H.subs(dict(zip(c,[a+b for a,b in zip(U,V,strict=True)],strict=True)))-hu-hv).applyfunc(s.expand)/2
    assert hu+3*B*s.Matrix(V)+81*c4*C==s.zeros(96,1)
    assert 2*huv+54*c4*B*s.Matrix(U)+486*c6*C==s.zeros(96,1)
    assert hv+162*c6*B*s.Matrix(U)-27*c4*B*s.Matrix(V)+729*c4*c4*C==s.zeros(96,1)
    ratio=s.cancel(s.discriminant(cubic,t)/delta);square=s.sqrt(ratio);assert square.is_Rational
    return {'beta_u_values':beta_values,'primitive_beta':list(map(str,vector)),
        'toric_free_c_indices':free,'toric_conics':list(map(str,conics)),'toric_cubic':str(cubic),'c_in_t':list(map(str,cvalues)),
        'all96_reduced_equations_span_three_conics':True,'toric_cubic_irreducible':True,
        'trace_form_determinant_mod_cubic':str(determinant),'trace_form_Bezout':list(map(str,bezout[:2])),
        'canonical_c4':str(c4),'canonical_c6':str(c6),'canonical_Delta':str(delta),'canonical_j':str(j),
        'E2_cubic':str(e2),'E2_root_index_used_for_search':root_index,'E2_X_in_t':str(X),'isomorphism_cyclic_determinant':str(change.det()),
        'exact_E2_polynomial_remainder_zero':True,'discriminant_ratio_square_root':str(square),
        'C5_residue':list(map(str,U)),'C9_residue':list(map(str,V)),
        'all20_residue_identities_checked':True,'three96_vector_identities_checked':True}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();start=time.time()
    model=S61();z,u,beta=primitive_beta(model);vector=s.Matrix([s.diff(beta[i,j],z[k]) for (i,j),k in COORDS])
    plane_path=ROOT/'reports/TORIC_UNIVERSAL_THREE_PLANE_20261101.json';plane=json.loads(plane_path.read_text());L=s.Matrix(plane['all51_homogeneous_linear_consequences']);results=[]
    for values in [[i+1 for i in range(16)],[(-1)**i*(i*i+3*i+5) for i in range(16)]]:
        substitution=dict(zip(u,values,strict=True));results.append(sample(vector.subs(substitution),L.subs(substitution),values))
    report={'id':'TORIC_TWO_TORSION_SAMPLES_20261101','all_passed':True,'field':'QQ','beta_u_order':list(map(str,u)),
        'samples':results,'sample_count':2,'proposed_generic_formula':'c=(C5(beta)*X+C9(beta))/(3X²−27c4(beta)), X³−27c4(beta)X−54c6(beta)=0',
        'scope':'Exact rational sample algebra isomorphisms and residue identities only. No generic identification with E[2] is proved by these samples.',
        'proof_note':'research/notes/TORIC_TWO_TORSION_SAMPLES_20261101.md'}
    sources=[Path(__file__),ROOT/'code/toric_algebra_fiber.py',ROOT/'code/elliptic_curve_invariants_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',plane_path]
    report['source_sha256']={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as out:json.dump(report,out,indent=2);out.write('\n')
    print(json.dumps({'all_passed':True,'exact_sample_isomorphisms':2,'generic_identification_proved':False,'seconds':round(time.time()-start,2),'output':str(args.output)}))

if __name__=='__main__':main()

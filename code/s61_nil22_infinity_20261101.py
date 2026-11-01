#!/usr/bin/env python3
"""Generic-W exclusion for every nonzero modular Jordan-(2,2) infinity stratum.

Uses the archived17-branch fixed-trace exhaustion as an explicit dependency.
All rank bounds and centralizer normalizations used here are recomputed.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from s61_infinity_flags import leading_equations_matrix
from s61_model import S61,X,J as J5
from weight1_structure import J,sp_basis

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(combinations(range(4),2))


def sm(matrix):return [[str(v) for v in row] for row in matrix.tolist()]


def tensor(values):
    Q=s.zeros(4)
    for (i,j),value in zip(PAIRS,values,strict=True):Q[i,j]=value;Q[j,i]=-value
    return Q


def beta_matrix(model,Q):
    B=X[0]**2*J5;B[1:,1:]+=Q
    point=model.normalize(B)['parameters'];matrix=leading_equations_matrix(model,point)
    assert all(s.cancel(v)==0 for v in matrix[:,0])
    return matrix[:,1:]


def polynomial_minor(M,values,variables,field=s.QQ):
    zero=M.subs(values)
    _,cols=DomainMatrix.from_Matrix(zero).convert_to(field).rref()
    rows=DomainMatrix.from_Matrix(zero[:,list(cols)].T).convert_to(field).rref()[1]
    domain=field.poly_ring(*variables)
    determinant=s.factor(domain.to_sympy(DomainMatrix.from_Matrix(M.extract(list(rows),list(cols))).convert_to(domain).det()))
    return list(rows),list(cols),determinant


def run():
    model=S61();x,y,z,w=X[1:];xyz=s.Matrix([x,y,z,w]);p=s.symbols('p0:3')
    summary_path=ROOT/'legacy/v07/data/nil22_classification_summary.json'
    fiber_path=ROOT/'legacy/v06/data/nilpotent22_certificate.json'
    summary=json.loads(summary_path.read_text());fiber=json.loads(fiber_path.read_text())
    local=dict(zip(['x','y','z','w'],list(xyz),strict=True));fixed=[]
    for label in ['Cplus','Cminus','Dplus','U0','U1','A0','A1','AB']:
        values=[]
        for i,j in PAIRS:
            value=s.sympify(summary[label]['Q'][str((i,j))],locals=local)
            if (i,j)==(1,3):value+=p[0]*x*x+p[1]*x*z+p[2]*z*z-(x*x+2*x*z+3*z*z)
            values.append(value)
        Q=tensor(values);M=beta_matrix(model,Q)
        field=s.QQ.algebraic_field(s.sqrt(3)) if label=='Dplus' else s.QQ
        rows,cols,det=polynomial_minor(M,dict.fromkeys(p,0),p,field)
        assert det!=0 and not det.has(*p)
        assert len(cols)>=(9 if label=='U0' else 10)
        entry={'label':label,'Q':sm(Q),'beta_matrix':sm(M),'rows':rows,'columns':cols,
               'uniform_rank_lower_bound':len(cols),'constant_minor':str(det)}
        if label=='U0':
            witness=M.subs(dict(zip(p,[1,2,3],strict=True)))
            ww,cc,dd=polynomial_minor(M,dict(zip(p,[1,2,3],strict=True)),p)
            assert len(cc)==10 and dd.subs(dict(zip(p,[1,2,3],strict=True)))!=0
            entry['rank10_witness']={'p':[1,2,3],'rows':ww,'columns':cc,'minor':str(dd),
                                    'evaluated_minor':str(dd.subs(dict(zip(p,[1,2,3],strict=True))))}
        fixed.append(entry)

    # The two closed lines retain their parameter and all magnetic terms.
    d=s.Symbol('d');lines=[]
    for label,c in [('L0',d),('L1',d-1)]:
        Q=tensor([(1-d)*x*z,0,(c-3*d+2)*x*x/2,-(c+d)*z*z/2,
                  p[0]*x*x+p[1]*x*z+p[2]*z*z+(2*d-c-1)*x*y+c*z*w,d*x*z])
        M=beta_matrix(model,Q);zero=M.subs(dict.fromkeys(p,0))
        nz=[i for i in range(M.rows) if any(zero[i,j]!=0 for j in range(16))]
        minors=[];gcd=s.S.Zero
        for shift in range(4):
            co=list(range(shift,16))+list(range(shift))
            offset=(3*shift)%len(nz);ro=nz[offset:]+nz[:offset]
            A=zero.extract(ro,co);_,pp=A.to_DM().rref();cc=[co[i] for i in pp]
            rr=A[:,list(pp)].T.to_DM().rref()[1];rr=[ro[i] for i in rr]
            determinant=s.factor(M.extract(rr,cc).to_DM().det().as_expr())
            assert len(cc)==12 and not determinant.has(*p)
            gcd=s.factor(s.gcd(gcd,determinant))
            minors.append({'rows':rr,'columns':cc,'determinant':str(determinant)})
        factors=s.factor_list(gcd,d)[1];assert all(s.degree(f,d)==1 for f,e in factors)
        roots=sorted({s.solve(f,d)[0] for f,e in factors})
        exceptions=[]
        for value in roots:
            special=M.subs(d,value)
            rr,cc,dd=polynomial_minor(special,dict.fromkeys(p,0),p)
            assert dd!=0 and not dd.has(*p)
            assert len(cc)>=(9 if label=='L0' and value in (0,1) else 10)
            exceptions.append({'d':str(value),'rows':rr,'columns':cc,'rank_bound':len(cc),'constant_minor':str(dd)})
        lines.append({'label':label,'Q':sm(Q),'beta_matrix':sm(M),'rank12_minors':minors,
                      'gcd':str(gcd),'exceptional_values':exceptions})

    # Rebuild the full open chart from the archived coordinate expression.
    ts=s.symbols('t1:12');base_loc={str(a):a for a in ts}|dict(zip(['x1','x2','x3','x4'],list(xyz),strict=True))
    raw=tensor([s.sympify(fiber['normal_form_brackets'][str(i+1)+str(j+1)],locals=base_loc) for i,j in PAIRS])
    chart={s.Symbol(key):s.sympify(value,locals=base_loc) for key,value in fiber['open_A_chart'].items()}
    A,b,v,dd=s.symbols('A b v d')
    open_Q=raw.subs(chart,simultaneous=True).subs({ts[8]:A,ts[1]:b,ts[6]:v,ts[7]:dd},simultaneous=True).applyfunc(s.cancel)
    rr,hh,rho=s.symbols('r h rho',nonzero=True)
    T=s.eye(4);T[1,0]=rr;T[3,2]=hh
    def pullback(Q,S):return (S.inv()*Q.subs(dict(zip(xyz,S*xyz,strict=True)),simultaneous=True)*S.inv().T).applyfunc(s.cancel)
    assert T.T*J*T==J
    transformed=pullback(open_Q,T)
    expected=open_Q.subs({b:b+4*A*rr,dd:dd+A*hh},simultaneous=True)
    assert (transformed-expected).applyfunc(s.cancel)==s.zeros(4)
    D=s.diag(rho,1/rho,1/rho,rho);assert D.T*J*D==J
    assert (pullback(open_Q,D)-open_Q.subs({A:A*rho**2,b:b*rho**4,v:v*rho**2},simultaneous=True)).applyfunc(s.cancel)==s.zeros(4)
    normalized=open_Q.subs({A:1,b:0,dd:0});M=beta_matrix(model,normalized)
    oo,cc,det=polynomial_minor(M,{v:0},[v])
    assert len(cc)==12 and det==s.Rational(2116013625,512)

    # The same torus normalizes each nonzero closed core parameter.
    a,bb,u,vv,c,ddd=s.symbols('a bb u vv c ddd')
    closed=raw.subs(dict(zip(ts,[a,bb,*p,c,vv,ddd,0,0,u],strict=True))).applyfunc(s.expand)
    closed_sub={a:a/rho**4,bb:bb*rho**4,u:u/rho**2,vv:vv*rho**2,
                p[0]:p[0]*rho**2,p[2]:p[2]/rho**2}
    assert (pullback(closed,D)-closed.subs(closed_sub,simultaneous=True)).applyfunc(s.cancel)==s.zeros(4)
    S=s.Matrix([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]])
    swap={s.Symbol(key):s.sympify(value,locals=base_loc) for key,value in fiber['centralizer_swap'].items()}
    assert (pullback(raw,S)-raw.subs(swap,simultaneous=True)).applyfunc(s.cancel)==s.zeros(4)
    modular=s.zeros(4);modular[1,2]=modular[3,0]=-1
    action=s.Matrix.hstack(*[(Z*modular-modular*Z).reshape(16,1) for Z in sp_basis()])
    assert action.rank()==6
    inputs=[Path(__file__),summary_path,fiber_path,ROOT/'code/s61_model.py',ROOT/'code/s61_infinity_flags.py']
    return {'id':'S61_NIL22_INFINITY_20261101','all_passed':True,'field':'QQ, and QQ(sqrt(3)) for Dplus',
        'fixed_cores':fixed,'closed_lines':lines,
        'open_chart':{'Q':sm(open_Q),'unipotent':sm(T),'parameter_action':['b->b+4*A*r','d->d+A*h'],
            'torus':sm(D),'torus_parameter_action':['A->A*rho^2','b->b*rho^4','v->v*rho^2','d->d'],
            'normalized_Q':sm(normalized),'beta_matrix':sm(M),'minor_rows':oo,'minor_columns':cc,
            'constant_rank12_minor':str(det)},
        'closed_torus':{'parameter_action':{str(aa):str(bb) for aa,bb in closed_sub.items()},
                        'all_coordinate_identities':True},
        'centralizer_swap_identity':True,'modular_orbit_dimension':6,
        'dimension_bounds':{'full_affine_nil22_Q':10,'projective_nil22_Q':9,
            'rank9_locus_affine':9,'rank9_locus_projective':8,'beta_incidence':15},
        'dependency':'legacy/v06 fixed-trace17-branch reduced exhaustion and its centralizer permutation.',
        'scope':'For sufficiently general beta, no nonzero modular Jordan-(2,2) point survives '
                'the original134 unsaturated leading equations. Not an all-Delta-open statement '
                'and not an affine-fiber or remaining-modular-strata exhaustion.',
        'proof_note':'research/notes/S61_NIL22_INFINITY_20261101.md',
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=run();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'fixed_cores':8,'closed_lines':2,'open_rank':12,
                      'beta_incidence_bound':15,'output':str(args.output)}))


if __name__=='__main__':main()

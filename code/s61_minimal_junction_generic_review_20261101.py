#!/usr/bin/env python3
"""Independent intrinsic-evaluation and dominance audit for minimal junction strata."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
ROOT=Path(__file__).resolve().parents[1]
def run():
    start=time.time();m=S61();z,u,B=primitive_beta(m);v=s.Matrix([0,1,0,0]);a,b,c,t=s.symbols('a b c t');w=s.Matrix([0,a,b,c]);row=v.T*B.subs(dict(zip(z,v)));omega=v.T*J[1:,1:];alpha=-row[0];tau=s.Matrix([row[2],row[3]]);tau_matrix=tau.jacobian(u);assert tau_matrix.rank()==2
    q=(w.T*B.subs(dict(zip(z,w)))*v)[0].expand();shift=w+t*v
    assert s.expand((shift.T*B.subs(dict(zip(z,shift)))*v)[0]-q-t*(w.T*B.subs(dict(zip(z,v)))*v)[0])==0
    K=s.Matrix.hstack(*tau_matrix.nullspace());r=s.symbols('r0:14');sub=dict(zip(u,K*s.Matrix(r)))
    assert (row-alpha*omega).subs(sub,simultaneous=True).applyfunc(s.expand)==s.zeros(1,4)
    qp=s.Poly(q.subs(sub,simultaneous=True),a,b,c);assert qp.degree(a)==0
    evaluations=s.Matrix([qp.coeff_monomial(b*b),qp.coeff_monomial(b*c),qp.coeff_monomial(c*c),alpha.subs(sub,simultaneous=True)])
    E=evaluations.jacobian(r);assert E.rank()==4
    reversed_cols=list(reversed(range(14)));ER=E[:,reversed_cols];piv=list(ER.to_DM().rref()[1]);minor=ER[:,piv].det();assert minor!=0
    qa,qb,qc,al,sigma=s.symbols('qa qb qc al sigma');f=s.Rational(25,9)*al**2-sigma**2*(qb*qb-4*qa*qc);det=s.hessian(f,[qa,qb,qc,al]).det();assert s.expand(det-s.Rational(1600,9)*sigma**6)==0
    # All nine Hamiltonian coefficients, not just the three baseline terms.
    x,y,zz,ww=X[1:];g=s.symbols('g0:9');G=sum(gg*mm for gg,mm in zip(g,[x*x*zz,x*x*ww,x*zz*zz,x*zz*ww,x*ww*ww,zz**3,zz*zz*ww,zz*ww*ww,ww**3]))
    for kind in ['M0','M1']:
        Q=s.zeros(4);Q[0,1]=x*x if kind=='M0' else x*x/2;Q[2,3]=0 if kind=='M0' else x*x/2;Q[1,2]=-s.diff(G,ww)-(x*zz if kind=='M1' else 0);Q[1,3]=s.diff(G,zz);Q-=Q.T
        divergence=s.Matrix([sum(s.diff(Q[i,j],X[j+1]) for j in range(4)) for i in range(4)])
        assert divergence==s.Matrix([0,-2*x,0,0])
    path=ROOT/'reports/S61_MINIMAL_JUNCTION_GENERIC_STRATA_20261101.json';main=json.loads(path.read_text());branches=json.loads((ROOT/'reports/S61_MINIMAL_JUNCTION_BRANCHES_20261101.json').read_text());review=json.loads((ROOT/'reports/S61_MINIMAL_JUNCTION_REVIEW_20261101.json').read_text());assert main['all_passed'] and review['all_passed'];out=[]
    allu=s.symbols('u1:62')
    for rec,known in zip(branches['cases'],main['cases']):
        kk=s.Rational(rec['k']);sig=2*kk+1;assert sig!=0 and known['sigma_squared']==str(sig**2)
        point=dict(zip(allu,map(s.Rational,rec['point61'])));qpoint=s.Poly(q.subs(point),b,c);A=qpoint.coeff_monomial(b*b);BB=qpoint.coeff_monomial(b*c);C=qpoint.coeff_monomial(c*c);ap=alpha.subs(point)
        assert BB*BB-4*A*C==1 and s.Rational(25,9)*ap**2==sig**2
        D=s.Matrix(rec['beta_differential']).applyfunc(s.Rational)[:,1:];columns=list(reversed(range(20)));DR=D[:,columns];pp=list(DR.to_DM().rref()[1]);dm=DR[:,pp].det(method='domain-ge');assert len(pp)==16 and dm!=0
        out.append(dict(kind=rec['kind'],k=rec['k'],sigma_squared=str(sig**2),alternative_beta_minor=str(dm),image_dimension=16,generic_multiplicity=rec['length_lower']))
    paths=[Path(__file__),path,ROOT/'code/s61_minimal_junction_generic_strata_20261101.py',ROOT/'research/notes/S61_MINIMAL_JUNCTION_GENERIC_STRATA_20261101.md',ROOT/'reports/S61_MINIMAL_JUNCTION_BRANCHES_20261101.json',ROOT/'reports/S61_MINIMAL_JUNCTION_SATURATION_20261101.json',ROOT/'reports/S61_MINIMAL_JUNCTION_REVIEW_20261101.json',ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json',ROOT/'reports/S61_MINIMAL_INFINITY_FORMAL_RANKS_20261101.json',ROOT/'reports/S61_MINIMAL_BASELINE_THICKENINGS_20261101.json',ROOT/'reports/S61_MINIMAL_CUBIC_TAIL_20261101.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py']
    return dict(id='S61_MINIMAL_JUNCTION_GENERIC_REVIEW_20261101',all_passed=True,vertex_rank=2,quotient_evaluation_rank=4,alternative_evaluation_minor=str(minor),quadric_hessian_determinant=str(det),all_nine_hamiltonian_modular_image_verified=True,cases=out,accepted_generic_entire_stratum_removal=True,scope='The eight main-junction image closures and Weyl-paired presentations; not all seventeen-dimensional curve/tail families.',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'cases':len(r['cases']),'seconds':r['seconds']}))

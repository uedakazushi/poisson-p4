#!/usr/bin/env python3
"""Intrinsic sixteen-dimensional vertex divisors supporting minimal junctions."""
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
ROOT=Path(__file__).resolve().parents[1]
def run():
    m=S61();z,u,B=primitive_beta(m);allu=s.symbols('u1:62');v=s.Matrix([0,1,0,0]);a,b,c=s.symbols('a b c');w=s.Matrix([0,a,b,c]);vertex={allu[3]:0,allu[11]:0};alpha=-3*allu[48]/5
    row=(v.T*B.subs(dict(zip(z,v)))).subs(vertex);assert row==alpha*(v.T*J[1:,1:])
    q=s.expand((w.T*B.subs(dict(zip(z,w)))*v)[0].subs(vertex));expected=allu[12]*b*b+(2*allu[13]-allu[48])*b*c-allu[4]*c*c;assert s.expand(q-expected)==0
    disc=s.discriminant(q.subs(c,1),b);assert s.expand(disc-(2*allu[13]-allu[48])**2-4*allu[12]*allu[4])==0
    evaluation=s.Matrix([-allu[11],allu[3],allu[12],2*allu[13]-allu[48],-allu[4],alpha]).jacobian(u);assert evaluation.rank()==6;ecols=list(evaluation.rref()[1]);edet=evaluation[:,ecols].det();assert edet!=0
    qa,qb,qc,al,sigma=s.symbols('qa qb qc al sigma');quad=s.Rational(25,9)*al**2-sigma**2*(qb**2-4*qa*qc);H=s.hessian(quad,[qa,qb,qc,al]);hdet=s.factor(H.det());assert hdet!=0 and s.factor(hdet/sigma**6).is_Rational
    branches_path=ROOT/'reports/S61_MINIMAL_JUNCTION_BRANCHES_20261101.json';sat_path=ROOT/'reports/S61_MINIMAL_JUNCTION_SATURATION_20261101.json';branches=json.loads(branches_path.read_text());sat=json.loads(sat_path.read_text());out=[]
    for rec in branches['cases']:
        k=s.Rational(rec['k']);sig=2*k+1;assert sig!=0;assert all(k!=a for a in [-s.Rational(7,9),-s.Rational(1,2),-s.Rational(2,9)])
        point=dict(zip(allu,list(map(s.Rational,rec['point61']))));assert s.expand(s.Rational(25,9)*alpha**2-sig**2*disc).subs(point)==0;assert disc.subs(point)==1
        D=s.Matrix(rec['beta_differential']).applyfunc(s.Rational)[:,list(range(1,21))];minor=rec['restricted_k_beta_minor'];det=D.extract(minor['rows_zero_based'],minor['columns_zero_based']).det(method='domain-ge');assert det==s.Rational(minor['determinant'])!=0
        upper=next(r for r in sat['cases'] if r['kind']==rec['kind'] and r['k']==rec['k']);assert upper['generic_rank']==rec['length_lower'] and upper['t_nilpotent']
        out.append(dict(kind=rec['kind'],k=rec['k'],sigma_squared=str(sig**2),vertex_divisor_dimension=16,junction_image_dimension=16,beta_dominant_minor=str(det),generic_multiplicity=rec['length_lower'],generic_beta_entire_stratum_removed=True))
    # The nonzero modular image of each main representative recovers [v].
    x,y,z0,w0=X[1:];aa,bb,cc=s.symbols('aa bb cc');G=aa*x*x*z0+bb*x*x*w0+cc*x*z0*w0
    for kind in ['M0','M1']:
        Q=s.zeros(4);Q[0,1]=x*x if kind=='M0' else x*x/2;Q[2,3]=0 if kind=='M0' else x*x/2;Q[1,2]=-s.diff(G,w0)-(x*z0 if kind=='M1' else 0);Q[1,3]=s.diff(G,z0);Q-=Q.T;mod=s.Matrix([sum(s.diff(Q[i,j],X[j+1]) for j in range(4)) for i in range(4)]);assert mod==s.Matrix([0,-2*x,0,0])
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',branches_path,sat_path,ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json',ROOT/'research/notes/S61_R0_GENERIC_EXHAUSTION_20261101.md']
    return dict(id='S61_MINIMAL_JUNCTION_GENERIC_STRATA_20261101',all_passed=True,field='QQ',q_v=str(q),alpha=str(alpha),discriminant=str(disc),evaluation_rank=6,evaluation_minor_columns=ecols,evaluation_minor=str(edet),quadric_hessian_determinant=str(hdet),weyl_pair='k maps to -1-k; sigma changes sign',cases=out,proof_note='research/notes/S61_MINIMAL_JUNCTION_GENERIC_STRATA_20261101.md',scope='Generic-beta removal of the eight main-junction image closures and their Weyl-paired presentations; not removal of entire seventeen-dimensional tail/curve families.',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'strata':len(r['cases'])}))
if __name__=='__main__':main()

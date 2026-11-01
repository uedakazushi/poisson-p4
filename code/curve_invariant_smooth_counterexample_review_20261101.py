#!/usr/bin/env python3
"""Independent global/pointwise review of the smooth-curve counterexample.

The polynomial pushforward, projective extension, Jacobi, symplectic frame and
primitive jet are independently reconstructed. The previously validated Fisher
invariant utility is shared; this does not recompute the QQ(t) formulas.
"""
import argparse, hashlib
from datetime import datetime, timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
import sys,json,sympy as s
from itertools import combinations
from elliptic_curve_invariants_20261101 import curve_invariants_from_beta
Y=s.Matrix(s.symbols('y0:4'));z=s.Matrix(s.symbols('z0:4'));lam=s.zeros(4);pairs=list(combinations(range(4),2));J=s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
for (i,j),v in zip(pairs,[2,8,12,3,5,7]):lam[i,j]=v;lam[j,i]=-v
phi=s.Matrix([z[0]+z[1]**2,z[1],z[2],z[3]]);push=phi.jacobian(z)*s.diag(*z)*lam*s.diag(*z)*phi.jacobian(z).T;P=push.subs(dict(zip(z,[Y[0]-Y[1]**2,Y[1],Y[2],Y[3]])),simultaneous=True).applyfunc(s.expand)
B=s.zeros(5);X=s.symbols('x0:5');B[0,1]=-2*X[2]**2;B[1,0]=-B[0,1]
for (i,j),v in zip(pairs,[2,8,12,3,5,7]):B[i+1,j+1]=v*X[i+1]*X[j+1];B[j+1,i+1]=-B[i+1,j+1]
for i,j in pairs:assert s.expand(P[i,j]-(B[i+1,j+1]-X[i+1]*B[0,j+1]+X[j+1]*B[0,i+1]).subs(dict(zip(X,[1,*Y]))))==0
assert all(s.expand(sum(P[i,h]*s.diff(P[j,k],Y[h])+P[j,h]*s.diff(P[k,i],Y[h])+P[k,h]*s.diff(P[i,j],Y[h]) for h in range(4)))==0 for i,j,k in combinations(range(4),3))
report=json.load(open(ROOT/'reports/CURVE_INVARIANT_SMOOTH_COUNTEREXAMPLE_20261101.json'));samples=[]
for t in [2,3]:
 p={Y[0]:t,Y[1]:1,Y[2]:1,Y[3]:1};K=P.subs(p);E=s.eye(4);r1=E[0,:];r2=E[1,:]/K[0,1]
 def orth(v):return v-(v*K*r2.T)[0]*r1+(v*K*r1.T)[0]*r2
 r3=orth(E[2,:]);r4=orth(E[3,:]);r4=r4/(r3*K*r4.T)[0];F=s.Matrix.vstack(r1,r2,r3,r4);assert F*K*F.T==J
 inv=F.inv();L=sum((Y[k]*F*sum((inv[h,k]*P.diff(Y[h]).subs(p) for h in range(4)),s.zeros(4))*F.T for k in range(4)),s.zeros(4));a=s.Matrix([s.diff(L[0,1]+L[2,3],v)/5 for v in Y]);v=J*a;Lp=L-2*a.dot(Y)*J+Y*v.T-v*Y.T;beta=-J*Lp*J
 vec=s.Matrix([s.expand(beta[i,j]).coeff(v) for i,j in pairs for v in Y]);actual=curve_invariants_from_beta(vec);orig=report['exact_rational_samples'][t-2];expected=tuple(s.Rational(orig[k]) for k in ['c4','c6','Delta','j']);assert actual==expected
 samples.append({'t':t,'frame':[[str(v) for v in row] for row in F.tolist()],'invariants':list(map(str,actual))});print(t,list(map(str,actual)),flush=True)
assert samples[0]['invariants'][3]!=samples[1]['invariants'][3]
paths=[Path(__file__),ROOT/'code/elliptic_curve_invariants_20261101.py',ROOT/'reports/CURVE_INVARIANT_SMOOTH_COUNTEREXAMPLE_20261101.json']
result={'id':'CURVE_INVARIANT_SMOOTH_COUNTEREXAMPLE_REVIEW_20261101','generated_utc':datetime.now(timezone.utc).isoformat(),'global_projective_extension':True,'independent_polynomial_pushforward':True,'direct_affine_Jacobi':True,'independent_normalization_samples':samples,'all_passed':True,'proof_note':'research/notes/CURVE_INVARIANT_SMOOTH_COUNTEREXAMPLE_REVIEW_20261101.md','input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'scope':'Independent global Jacobi and two pointwise invariants; Fisher contraction implementation shared, QQ(t) invariant formulas not independently repeated.'}
args.output.parent.mkdir(parents=True,exist_ok=True)
with args.output.open('x') as handle:
 json.dump(result,handle,indent=2);handle.write('\n')

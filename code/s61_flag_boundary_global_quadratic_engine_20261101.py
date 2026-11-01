#!/usr/bin/env python3
"""Exact ordinary-flag fourth-order core cover: quadratic."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
WD=Path(sys.argv[1])
import json,pickle,time
from pathlib import Path
from itertools import combinations_with_replacement as cwr
import sympy as s
from sympy.polys.matrices import DomainMatrix as DM
st=time.time();r=s.Symbol('b2');P=s.QQ.poly_ring(r);cover=json.load(open(ROOT/'reports/S61_FLAG_BOUNDARY_QUADRATIC_COVER_20261101.json'));mu=s.Matrix(list(map(s.sympify,cover['socle_multiplication_coefficients'])));maxdeg=max(s.Poly(p,r).degree()for p in mu)
for bound in range(5):
 n=10*(bound+1);A=s.zeros(maxdeg+bound+1,n+1);A[0,n]=1
 for i,p in enumerate(mu):
  for(e,),c in s.Poly(p,r).terms():
   for k in range(bound+1):A[e+k,i*(bound+1)+k]=c
 R,piv=A.to_DM().rref();print('socle bound',bound,'rank',len(piv),'consistent',all(p<n for p in piv),time.time()-st,flush=True)
 if all(p<n for p in piv):break
else:raise RuntimeError('no low socle')
z=s.zeros(1,10);RR=R.to_Matrix()
for j,p in enumerate(piv):i,k=divmod(p,bound+1);z[0,i]+=RR[j,n]*r**k
assert s.expand((z*mu)[0])==1
cache=pickle.load(open(WD/'flag-b1-universal-quadratic-cache.pkl','rb'));q=json.load(open(WD/'flag-b1-universal-quadratic.json'));mons=list(cwr(range(9),2));cols=[j for j,m in enumerate(mons)if all(i in {4,5,7,8}for i in m)];rawrows=q['pivot_rows'];contents=list(map(s.sympify,cover['row_contents']))
def dec(poly):return sum(s.Rational(c)*r**e[0]for e,c in poly)
B=s.Matrix([[s.cancel(dec(cache['A2'][i][j])/contents[k])for j in cols]for k,i in enumerate(rawrows)]);assert (B*mu).applyfunc(s.expand)==s.zeros(9,1)
Q=s.eye(10)-mu*z;piv=list(B.to_DM().rref()[1]);K=(Q[:,piv].to_DM().convert_to(s.QQ.frac_field(r))*B[:,piv].to_DM().to_field().inv()).to_Matrix();K=K.applyfunc(s.cancel)
assert (K*B-Q).applyfunc(s.cancel)==s.zeros(10);assert (B*K-s.eye(9)).applyfunc(s.cancel)==s.zeros(9)
enc=lambda M:[[str(v)for v in row]for row in M.tolist()]
Path(WD/'flag-global-low-quad.json').write_text(json.dumps(dict(mu=enc(mu),socle=enc(z),basis=enc(B),right_inverse=enc(K),row_contents=list(map(str,contents)),selected_raw_rows=rawrows,socle_degree=bound,basis_degree=max(s.Poly(p,r).degree()for p in B),inverse_denominator=str(s.factor(s.lcm([s.denom(p)for p in K]))),seconds=time.time()-st),indent=2)+'\n');print('LOW QUAD GLOBAL SPLIT PASS',time.time()-st,flush=True)

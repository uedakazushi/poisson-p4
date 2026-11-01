#!/usr/bin/env python3
"""T7 flag flat4 certificate: alternate coordinate lattice stage. Rebuild in supplied directory."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
WD=Path(sys.argv[1]);WD.mkdir(parents=True,exist_ok=True)
import sys,json,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
st=time.time();data=json.load(open(WD/'flag-T7-universal-transport.json'));b=list(map(s.Symbol,data['variables']));e=s.Symbol('e');C=[s.Matrix(data['coefficients'][str(i)]).applyfunc(s.sympify)for i in [1,2,3]];R,*_=ring(b+[e],s.QQ);ds=[s.lcm([s.denom(v)for v in M])for M in C];dp=[R.from_expr(v)for v in ds];P=[[R.from_expr(s.cancel(C[j][i]*ds[j]))for j in range(3)]for i in range(45)]
indices=[9,0,10];B=[P[i]for i in indices]
# Cyclic cofactor indices already have positive permutation sign.
adj=[[B[(j+1)%3][(i+1)%3]*B[(j+2)%3][(i+2)%3]-B[(j+1)%3][(i+2)%3]*B[(j+2)%3][(i+1)%3]for j in range(3)]for i in range(3)]
det=sum((B[0][j]*adj[j][0]for j in range(3)),R.zero);assert det
assert all(sum((B[i][k]*adj[k][j]for k in range(3)),R.zero)==(det if i==j else R.zero)for i in range(3)for j in range(3))
def val(p):return min(m[-1]for m in p)if p else 10000
def divide_e(p,k):return R.from_dict({m[:-1]+(m[-1]-k,):c for m,c in p.items()})
vd=val(det);D=divide_e(det,vd);print('det terms',len(det),'valuation',vd,'regular det terms',len(D),'seconds',time.time()-st,flush=True)
nums=[]
for row in P:
 out=[sum((row[k]*adj[k][j]for k in range(3)),R.zero)for j in range(3)];assert all(val(v)>=vd for v in out);nums.append([divide_e(v,vd)for v in out])
print('coordinate regularity PASS',time.time()-st,flush=True)
# Each selected coordinate function has pure-nu coefficients B[i,k]/d[k].
# Multiply two functions and express in the selected coordinate basis.
mult=[]
for i in range(3):
 for j in range(i,3):
  # Product rescaled by d[k] has entries [0, B_i0 B_j0*d1/d0^2,
  # (B_i0 B_j1+B_i1 B_j0)*d2/(d0*d1)]. The ratios are Laurent monomials.
  expr1=s.cancel(ds[1]/ds[0]**2);expr2=s.cancel(ds[2]/(ds[0]*ds[1]));denprod=s.lcm([s.denom(expr1),s.denom(expr2)]);pp=R.from_expr(denprod)
  pr=[R.zero,B[i][0]*B[j][0]*R.from_expr(expr1*denprod),(B[i][0]*B[j][1]+B[i][1]*B[j][0])*R.from_expr(expr2*denprod)]
  numer=[sum((pr[k]*adj[k][j]for k in range(3)),R.zero)for j in range(3)];den=pp*det;vden=val(den);assert all(val(v)>=vden for v in numer)
  mult.append(dict(numerators=[str(divide_e(v,vden).as_expr())for v in numer],denominator=str(divide_e(den,vden).as_expr())))
print('multiplication regularity PASS',time.time()-st,flush=True)
enc=lambda A:[[str(v.as_expr())for v in row]for row in A]
out=dict(variables=data['variables'],selected_q_rows=indices,coordinate_numerators=enc(nums),coordinate_denominator=str(D.as_expr()),multiplication=mult,column_denominators=list(map(str,ds)),seconds=time.time()-st)
(WD/'flag-T7-universal-lattice-ring.json').write_text(json.dumps(out,indent=2)+'\n');print('UNIVERSAL LATTICE RING PASS',time.time()-st,flush=True)

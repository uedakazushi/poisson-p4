#!/usr/bin/env python3
"""Exact rational eight-variable flag calculation: inverse stage.

Internal engine of s61_flag_boundary_eight_cubic_20261101.py.
All temporary products are rebuilt in the supplied working directory.
"""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
WD=Path(sys.argv[1]); WD.mkdir(parents=True,exist_ok=True)
import sys,json,time
from pathlib import Path
sys.path.insert(0,str(ROOT/'code'))
import sympy as s
from sympy.polys.rings import ring
start=time.time();r=json.load(open(WD / 'flag-b1-normalized-frame.json'));J=s.Matrix(r['jacobian']).applyfunc(s.sympify);rows=r['pivot_rows'];cols=r['pivot_columns'];A=J[rows,cols];variables=list(map(s.Symbol,r['variables']));core=s.Symbol('b2');zero={x:0 for x in variables if x!=core};A0=A.subs(zero);iv0=A0.inv(method='DM');den=s.factor(s.lcm([s.denom(v) for v in iv0]));print('core inverse den',den,'sec',time.time()-start,flush=True)
R,*_=ring(variables,s.QQ);poly=lambda v:R.from_expr(s.cancel(v));II=[[poly(den*v) for v in row] for row in iv0.tolist()];DD=[[poly(v) for v in row] for row in (A-A0).tolist()];D=poly(den)
def mult(U,V):
 W=[[R.zero for _ in range(len(V[0]))] for i in range(len(U))]
 for i,row in enumerate(U):
  for j,c in enumerate(row):
   if c:
    for k,v in enumerate(V[j]):
     if v:W[i][k]+=c*v
 return W
N=mult(II,DD);edges={i:{j for j,v in enumerate(row) if v} for i,row in enumerate(N)}
remaining=set(range(36));order=[]
while remaining:
 leaves=[i for i in remaining if not(edges[i]&remaining)]
 if not leaves:print('cycle',remaining,flush=True);break
 order+=leaves;remaining-=set(leaves)
print('acyclic',not remaining,'edges',sum(len(v) for v in edges.values()),'sec',time.time()-start,flush=True)
assert not remaining
# Solve (den*Id+N)X=II by the acyclic dependency order. Fractions have only core denominators.
X=[[R.zero for _ in range(36)] for i in range(36)];powers=[0]*36
for i in order:
 power=1+max([powers[j] for j in edges[i]]+[0]);powers[i]=power
 for k in range(36):
  v=II[i][k]*D**(power-1)
  for j in edges[i]:v-=N[i][j]*X[j][k]*D**(power-1-powers[j])
  X[i][k]=v
 for j in range(power):
  divs=[v.div(D) for v in X[i]]
  if any(rem for q,rem in divs):break
  X[i]=[q for q,rem in divs];powers[i]-=1
iv=s.Matrix([[s.cancel(v.as_expr()/den**powers[i]) for v in row] for i,row in enumerate(X)])
assert (A*iv-s.eye(36)).applyfunc(s.cancel)==s.zeros(36)
Path(WD / 'flag-b1-triangular-inverse.json').write_text(json.dumps(dict(inverse=[[str(v) for v in row] for row in iv.tolist()],core_denominator=str(den),dependency_order=order,seconds=time.time()-start),indent=2)+'\n');print('TRIANGULAR INVERSE PASS',time.time()-start,flush=True)

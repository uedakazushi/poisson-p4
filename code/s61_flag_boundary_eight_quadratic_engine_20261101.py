#!/usr/bin/env python3
"""Exact rational eight-variable flag calculation: quadratic stage.

Internal engine of s61_flag_boundary_eight_cubic_20261101.py.
All temporary products are rebuilt in the supplied working directory.
"""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
WD=Path(sys.argv[1]); WD.mkdir(parents=True,exist_ok=True)
import sys,json,time
from pathlib import Path
from itertools import combinations_with_replacement
sys.path.insert(0,str(ROOT/'code'))
import sympy as s
from sympy.polys.rings import ring
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
start=time.time();m=S61();r=json.load(open(WD / 'flag-b1-normalized-frame.json'));case=json.load(open(ROOT / 'reports' / 'S61_FLAG_UNIFORM_OPEN_20261101.json'))['flag'];H=s.Matrix(case['beta_kernel']).applyfunc(s.Rational);b=s.symbols('b:15');sub={s.Symbol(k):s.sympify(v) for k,v in r['substitutions'].items()};bet=(H*s.Matrix(b)).subs(sub);one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1]
_,F=specialize(m,dict(zip(one,bet)),True);J=s.Matrix(r['jacobian']).applyfunc(s.sympify);T=s.Matrix(r['tangent']).applyfunc(s.sympify);rows=r['pivot_rows'];cols=r['pivot_columns'];U=s.zeros(45,36);iv=s.Matrix(json.load(open(WD / 'flag-b1-triangular-inverse.json'))['inverse']).applyfunc(s.sympify)
for i,c in enumerate(cols):U[c,:]=iv[i,:]
ud=s.factor(s.lcm([s.denom(v) for v in U]));td=s.factor(s.lcm([s.denom(v) for v in T]));print('inverse denominators',ud,td,'sec',time.time()-start,flush=True)
variables=list(map(s.Symbol,r['variables']));R,*_=ring(variables,s.QQ);zero=R.zero
poly=lambda v:R.from_expr(s.cancel(v))
UU=[[poly(ud*v) for v in row] for row in U.tolist()];TT=[[poly(td*v) for v in row] for row in T.tolist()];JJ=[[poly(v) for v in row] for row in J.tolist()];UD=poly(ud)
two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];q=[s.Rational(case['point61'][i]) for i in two];fixed=next(i for i,v in enumerate(q) if v);chart=[i for i in range(45) if i!=fixed];mons=list(combinations_with_replacement(range(9),2));Q=[[zero for _ in mons] for i in range(135)]
for row,f in enumerate(F):
 for mon,c in f.items():
  if fixed in mon:continue
  i,j=[chart.index(k) for k in mon]+[44]*(2-len(mon));co=poly(c)
  for k,(a,bb) in enumerate(mons):
   value=TT[i][a]*TT[j][bb]
   if a!=bb:value+=TT[i][bb]*TT[j][a]
   if value:Q[row][k]+=co*value
print('raw quadratic sec',time.time()-start,flush=True)
W=[[zero for _ in mons] for i in range(45)]
for i,row in enumerate(UU):
 for j,co in enumerate(row):
  if co:
   for k,v in enumerate(Q[rows[j]]):
    if v:W[i][k]-=co*v
A=[[UD*v for v in row] for row in Q]
for i,row in enumerate(JJ):
 for j,co in enumerate(row):
  if co:
   for k,v in enumerate(W[j]):
    if v:A[i][k]+=co*v
active=sorted({k for row in A for k,v in enumerate(row) if v});print('active monomials',[mons[k] for k in active],flush=True)
import pickle
enc=lambda M:[[[[list(e),str(c)] for e,c in v.items()] for v in row] for row in M]
pickle.dump(dict(A2=enc(A),W2=enc(W),T=enc(TT),U=enc(UU),J=enc(JJ),F=[[[list(mon),str(c)] for mon,c in f.items()] for f in F],chart=chart,fixed=fixed,ud=str(ud),td=str(td),rows=rows,variables=list(map(str,variables))),open(WD / 'flag-b1-universal-quadratic-cache.pkl','wb'))
AA=s.Matrix([[v.as_expr() for v in row] for row in A]);sample={x:0 for x in variables};sample[s.Symbol('b2')]=1;chosen=list(AA.subs(sample).T.to_DM().rref()[1]);assert len(chosen)==9;R9,piv=AA[chosen,:].to_DM().rref();R9=R9.to_Matrix()[:9,:];print('quad rref symbols',sorted(set().union(*(v.free_symbols for v in R9)),key=str),'sec',time.time()-start,flush=True)
assert (AA-AA[:,piv]*R9).applyfunc(s.cancel)==s.zeros(135,45)
Path(WD / 'flag-b1-universal-quadratic.json').write_text(json.dumps(dict(variables=r['variables'],right_inverse=[[str(v) for v in row] for row in U.tolist()],right_inverse_denominator=str(ud),tangent_denominator=str(td),quadratic_rref=[[str(v) for v in row] for row in R9.tolist()],pivot_rows=chosen,pivot_columns=list(piv),monomials=mons,seconds=time.time()-start),indent=2)+'\n');print('UNIVERSAL QUAD RANK9 PASS',time.time()-start,flush=True)

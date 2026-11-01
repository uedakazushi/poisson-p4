#!/usr/bin/env python3
"""Exact h=0 ordinary-flag fourth-order core cover: project4."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
WD=Path(sys.argv[1])
import json,pickle,time
from pathlib import Path
from itertools import combinations_with_replacement as cwr
import sympy as s
from sympy.polys.rings import ring
st=time.time();c=pickle.load(open(WD/'flag-eight-flint-projection-input.pkl','rb'));v=list(map(s.Symbol,c['variables']));R,*_=ring(v,s.QQ);r=s.Symbol('b2');b8=s.Symbol('b8');active={4,5,7,8};rad=[0,1,2,3,6];ms={d:list(cwr(range(9),d))for d in [2,3,4]};q=json.load(open(WD/'flag-global-low-quad.json'));mu={m:s.sympify(q['mu'][i][0])for i,m in enumerate(cwr(sorted(active),2))};G=(r-2)**2
quot=[m for m in ms[4]if sum(i in active for i in m)<=1]+[(0,0,8,8),(0,1,8,8),(1,1,8,8)];assert len(quot)==213;ix={m:i for i,m in enumerate(quot)};sub={0:(0,1),1:(1,1),2:(1,0),3:(1,0),6:(1,0)}
proj=[]
for m in ms[4]:
 a=tuple(i for i in m if i in active);radm=tuple(i for i in m if i not in active)
 if len(a)>=3:proj.append(None)
 elif len(a)<=1:proj.append((ix[m],R.from_expr(G)))
 else:
  rr=tuple(sorted(sub[i][0]for i in radm));co=s.cancel(G*mu[a]*s.prod(sub[i][1]for i in radm));assert s.denom(co).is_Number;proj.append((ix[tuple(sorted(rr+(8,8)))],R.from_expr(co)))
print('projection ready',time.time()-st,flush=True)
out=[]
for i,row in enumerate(c['A4_surviving']):
 result=[R.zero for _ in quot]
 for j,poly in zip(c['A4_columns'],row):
  if not poly:continue
  dst,factor=proj[j];p=R.from_dict({tuple(e):s.QQ.convert(s.Rational(cc))for e,cc in poly});result[dst]+=p*factor
 out.append(result)
 if i%20==0:print('row',i,'terms',sum(len(x)for x in result),'seconds',time.time()-st,flush=True)
enc=lambda M:[[[[list(e),str(cc)]for e,cc in p.items()]for p in row]for row in M]
tout=[R.zero for _ in quot]
for j,poly in zip(c['A4_columns'],c['T4_surviving'][0]):
 if poly:
  dst,factor=proj[j];tout[dst]+=factor*R.from_dict({tuple(e):s.QQ.convert(s.Rational(cc))for e,cc in poly})
c['T4_projected']=enc([tout]);del c['T4_surviving']
c['A4_projected']=enc(out);c['projection_multiplier']=str(G);c['quotient_monomials']=quot;c['projection']=[None if p is None else [p[0],str(p[1].as_expr())]for p in proj]
del c['A4_surviving'];pickle.dump(c,open(WD/'flag-eight-flint-projected4.pkl','wb'));print('PROJECTION DONE',sum(len(p)for row in out for p in row),(WD/'flag-eight-flint-projected4.pkl').stat().st_size,time.time()-st,flush=True)

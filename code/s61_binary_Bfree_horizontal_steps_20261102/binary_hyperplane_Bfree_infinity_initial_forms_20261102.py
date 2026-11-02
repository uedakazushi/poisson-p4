from pathlib import Path
import json,time
import sympy as s
from sympy.polys.rings import ring
st=time.time();tw=[1,0,0,0,1,0,1,1];R,a,e=ring('a,e',s.QQ);rows=[]
with Path('/tmp/poisson-hyperplane-Bfree-cubic-20261102/rows').open() as f:
 for _ in range(int(f.readline())):
  idx,w,ct=map(int,f.readline().split());rr=[]
  for _ in range(ct):
   m=tuple(map(int,f.readline().split()));z=f.readline().split()
   for k in range(int(z[0])):
    c=s.QQ(z[1+4*k]);ex=tuple(map(int,z[2+4*k:5+4*k]));rr.append((m,ex,c))
  maxdegree=max(ex[2]-sum(i*j for i,j in zip(m,tw)) for m,ex,c in rr);q={}
  for m,ex,c in rr:
   if ex[2]-sum(i*j for i,j in zip(m,tw))==maxdegree:q.setdefault(m,{})[ex[:2]]=c
  q={m:R.from_dict(c) for m,c in q.items()};g=R.zero
  for c in q.values():g=g.gcd(c)
  q={m:c.exquo(g) for m,c in q.items()};rows.append(dict(source=idx,weight=w,d_shift=maxdegree,content=str(g.as_expr()),terms=[(list(m),str(c.as_expr())) for m,c in q.items()]))
print('forms',[(r['source'],r['weight'],len(r['terms'])) for r in rows],flush=True);print('coreterms',sum(len(R.from_expr(s.sympify(c))) for r in rows for m,c in r['terms']),flush=True);Path(__file__).with_suffix('.json').write_text(json.dumps(dict(weights=tw,rows=rows,seconds=time.time()-st),indent=2)+'\n');print('seconds',time.time()-st)

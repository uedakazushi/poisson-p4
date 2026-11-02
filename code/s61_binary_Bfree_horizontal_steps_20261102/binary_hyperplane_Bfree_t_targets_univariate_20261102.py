from pathlib import Path
import json,itertools,functools,time
p=32003;a0,e0=1,2;wd=Path('/tmp/poisson-hyperplane-Bfree-cubic-20261102');start=time.time();nw=[1,1,1,1,1,2,1,1]
def mons(w,i=0,m=()):
 if i==8:
  if w==0:yield m
  return
 for k in range(w//nw[i]+1):yield from mons(w-k*nw[i],i+1,m+(k,))
@functools.lru_cache(None)
def mod(c):
 if '/' in c:
  a,b=map(int,c.split('/'));assert b%p;return a%p*pow(b%p,-1,p)%p
 return int(c)%p
def poly(z,dim=3):
 q={}
 for k in range(int(z[0])):
  v=mod(z[1+(dim+1)*k])*pow(a0,int(z[2+(dim+1)*k]),p)*pow(e0,int(z[3+(dim+1)*k]),p)%p;d=int(z[4+(dim+1)*k]);q[d]=(q.get(d,0)+v)%p
 return {d:v for d,v in q.items() if v}
def addto(q,r,scale=1):
 for d,v in r.items():q[d]=(q.get(d,0)+scale*v)%p
 return {d:v for d,v in q.items() if v}
def mul(q,r):
 z={}
 for i,a in q.items():
  for j,b in r.items():z[i+j]=(z.get(i+j,0)+a*b)%p
 return {i:a for i,a in z.items() if a}
rows=[]
with (wd/'rows').open() as f:
 for i in range(int(f.readline())):
  id,w,ct=map(int,f.readline().split());q={}
  for j in range(ct):m=tuple(map(int,f.readline().split()));q[m]=poly(f.readline().split())
  for shift in mons(4-w):rows.append({tuple(x+y for x,y in zip(m,shift)):c for m,c in q.items() if c})
ms=[m for m in mons(4) if m[5]!=2];ix={m:i for i,m in enumerate(ms)};assert len(ms)==238 and len(rows)==321
with (wd/'coeff').open() as f:
 n=int(f.readline());z=[str(n)]
 for _ in range(n):z+=f.readline().split()
 ud=poly(z)
t={}
with (wd/'raw').open() as f:
 assert list(map(int,f.readline().split()))==[135,3,8,3]
 for degree in [2,3]:
  for i in range(180):
   ct=int(f.readline())
   for _ in range(ct):
    z=f.readline().split()
    if i!=44:continue
    m=tuple(map(int,z[4:]));assert sum(x*y for x,y in zip(m,nw))==3;v=mod(z[0])*pow(a0,int(z[1]),p)*pow(e0,int(z[2]),p)%p;d=int(z[3]);q={d:v}
    if degree==2:q=mul(q,ud)
    t[m]=addto(t.get(m,{}),q)
target=[]
for i in [0,1,2,3,4,6,7]:target.append({tuple(x+(j==i) for j,x in enumerate(m)):q for m,q in t.items() if q})
path=Path('/tmp/hyperplane-Bfree-t-targets-univariate-20261102.in')
with path.open('w') as f:
 f.write(f'{p} {len(rows)} {len(ms)} {len(target)}\n80 102\n')
 for row in rows+target:
  f.write(str(len(row))+'\n')
  for m,q in row.items():
   assert m in ix;f.write(str(ix[m])+' '+str(max(q)+1)+' '+' '.join(str(q.get(i,0)) for i in range(max(q)+1))+'\n')
print('READY',len(rows),len(ms),'tterms',len(t),'seconds',time.time()-start,flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps(dict(point=[a0,e0],rows=len(rows),columns=len(ms),target_count=7,target_normal_terms=len(t),seconds=time.time()-start),indent=2)+'\n')

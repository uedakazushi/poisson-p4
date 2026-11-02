from pathlib import Path
import json,itertools,functools,time
p=32003;a0,e0=1,2;wd=Path('/tmp/poisson-hyperplane-Bfree-cubic-20261102');start=time.time();nw=[1,1,1,1,1,2,1,1]
def mons(w,i=0,m=()):
 if i==8:
  if w==0:yield m
  return
 for k in range(w//nw[i]+1):yield from mons(w-k*nw[i],i+1,m+(k,))
import sympy as sy
from sympy.polys.rings import ring
R,a,e,d=ring('a,e,d',sy.QQ)
def poly(z,dim=3):
 return R.from_dict({tuple(map(int,z[2+4*k:5+4*k])):sy.QQ(z[1+4*k]) for k in range(int(z[0]))})
def addto(q,r,scale=1):return q+scale*r
def mul(q,r):return q*r
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
    m=tuple(map(int,z[4:]));assert sum(x*y for x,y in zip(m,nw))==3;q=R.from_dict({tuple(map(int,z[1:4])):sy.QQ(z[0])})
    if degree==2:q=mul(q,ud)
    t[m]=addto(t.get(m,R.zero),q)
target=[]
for i in [0,1,2,3,4,6,7]:target.append({tuple(x+(j==i) for j,x in enumerate(m)):q for m,q in t.items() if q})

path=Path('/tmp/hyperplane-Bfree-t-targets-exact-core-20261102.in')
with path.open('w') as f:
 f.write(f'{len(rows)+len(target)} {len(ms)} {len(rows)}\n')
 for row in rows+target:
  for m in ms:
   q=row.get(m,R.zero);f.write(str(len(q)))
   for mm,c in q.items():f.write(' '+str(c)+' '+' '.join(map(str,mm)))
   f.write('\n')
print('READY',len(rows),len(ms),'terms',sum(len(v) for r in rows for v in r.values()),'seconds',time.time()-start)

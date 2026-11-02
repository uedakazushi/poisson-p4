from pathlib import Path
import json,time,functools,sys
import sympy as s,numpy as np
from sympy.polys.rings import ring
from sympy.polys.matrices import DomainMatrix
sys.path.insert(0,'code');from s61_binary_deeper_uniform_quadratic_20261101 import elimination_rank,determinant_mod
st=time.time();R,a,e=ring('a,e',s.QQ);K=R.to_field();weights=[1,1,1,1,1,2,1,1];light=[0,1,2,3,4,6,7]
normal=json.loads(Path('work/generic_Bfree_infinity_constant_form_20261102.json').read_text());initial=json.loads(Path('work/binary_hyperplane_Bfree_infinity_initial_forms_20261102.json').read_text());p4=json.loads(Path('work/binary_hyperplane_infinity_constant_quadratic_projection_20261102.json').read_text());p3=json.loads(Path('work/binary_hyperplane_infinity_constant_cubic_projection_20261102.json').read_text())
N=[[K.from_expr(s.sympify(c))for c in row]for row in normal['coordinate_inverse8']];D=R.one
for i in light:
 for c in N[i]:D=D.lcm(c.denom)
print('D',D,flush=True);NN=[[c.numer*D.exquo(c.denom)for c in N[i]]for i in light];Nnum={i:NN[k]for k,i in enumerate(light)}
def mons(w,i=0,m=()):
 if i==8:
  if w==0:yield m
  return
 for k in range(w//weights[i]+1):yield from mons(w-k*weights[i],i+1,m+(k,))
ms={j:list(mons(j))for j in [1,2,3]};ms[4]=[m for m in mons(4)if m[5]!=2];ix={j:{m:i for i,m in enumerate(mons)}for j,mons in ms.items()}
qs=[{tuple(m):s.QQ(c)for m,c in row}for row in p4['quadrics']]
A=DomainMatrix({i:{ix[2][m]:c for m,c in q.items()}for i,q in enumerate(qs)},(6,29),s.QQ);rr,piv=A.rref();fr=[i for i in range(29)if i not in piv];fi={i:j for j,i in enumerate(fr)};raw=rr.to_dod();B2={i:{j:s.QQ(1)}for i,j in fi.items()}
for j,i in enumerate(piv):B2[i]={fi[k]:-c for k,c in raw[j].items()if k in fi}
proj={2:[B2.get(i,{})for i in range(29)],3:[{j:s.QQ(c)for j,c in row}for row in p3['projection_sparse']],4:[{j:s.QQ(c)for j,c in row}for row in p4['projection_sparse']]};free={2:fr,3:p3['free_columns'],4:p4['free_columns']};reps={w:[ms[w][j]for j in free[w]]for w in [2,3,4]}
def addto(out,j,p):
 if p:
  z=out.get(j,R.zero)+p
  if z:out[j]=z
  elif j in out:del out[j]
def multlinear(f,degree,i):
 out={}
 for k,c in f.items():
  m=reps[degree][k]
  for j in light:
   if not Nnum[i][j]:continue
   n=tuple(v+int(h==j)for h,v in enumerate(m));z=c*Nnum[i][j]
   for col,scale in proj[degree+1][ix[degree+1][n]].items():addto(out,col,z*scale)
 return out
@functools.lru_cache(None)
def square(i,j):
 out={}
 for k in light:
  for h in light:
   if Nnum[i][k] and Nnum[j][h]:
    mon=tuple(int(v==k)+int(v==h)for v in range(8));z=Nnum[i][k]*Nnum[j][h]
    for col,c in proj[2][ix[2][mon]].items():addto(out,col,z*c)
 return out
@functools.lru_cache(None)
def cubic(mon):
 if mon[5]:
  i=next(i for i in light if mon[i]);out={}
  for j in light:
   n=tuple(int(k==j)+int(k==5)for k in range(8))
   for col,c in proj[3][ix[3][n]].items():addto(out,col,Nnum[i][j]*c)
  return out
 ids=[i for i in light for _ in range(mon[i])];return multlinear(square(ids[0],ids[1]),2,ids[2])
@functools.lru_cache(None)
def parse(c):return R.from_expr(s.sympify(c))
# Independently check the original six leading quadrics reduce to zero after this change.
for row in initial['rows']:
 if row['weight']!=2:continue
 out={}
 for m,c in row['terms']:
  ids=[i for i in light for _ in range(m[i])];assert len(ids)==2
  for col,z in square(*ids).items():addto(out,col,parse(c)*z)
 assert not out,row['source']
print('Q2_PASS',time.time()-st,flush=True)
cr=[];contents=[];sources=[]
for row in initial['rows']:
 if row['weight']!=3:continue
 out={}
 for m,c in row['terms']:
  for col,z in cubic(tuple(m)).items():addto(out,col,parse(c)*z)
 g=R.zero
 for z in out.values():g=g.gcd(z)
 if g:out={j:z.exquo(g)for j,z in out.items()}
 cr.append(out);contents.append(str(g.as_expr()));sources.append(row['source']);print('CUBIC',row['source'],len(out),sum(len(z)for z in out.values()),time.time()-st,flush=True)
rows=[]
for f in cr:
 for i in light:
  out={}
  for k,c in f.items():
   m=reps[3][k];n=tuple(v+int(h==i)for h,v in enumerate(m))
   for col,scale in proj[4][ix[4][n]].items():addto(out,col,c*scale)
  g=R.zero
  for z in out.values():g=g.gcd(z)
  if g:out={j:z.exquo(g)for j,z in out.items()}
  rows.append(out)
prime=32003
@functools.lru_cache(None)
def evaluate(z):
 ans=0
 for (i,j),c in z.items():ans=(ans+int(c.numerator)*pow(int(c.denominator),-1,prime)*pow(2,j,prime))%prime
 return ans
num=np.array([[evaluate(row.get(j,R.zero))for j in range(106)]for row in rows],dtype=np.int64);rank,cols=elimination_rank(num,prime);_,rr=elimination_rank(num[:,cols].T,prime);assert rank==31;fr=[j for j in range(106)if j not in cols];minor=determinant_mod(num[np.ix_(rr,cols)],prime)
path=Path('/tmp/hyperplane-Bfree-infinity-reduced-kernel-20261102.in')
with path.open('w')as f:
 f.write(f'{len(rows)} 106 {rank} {len(fr)}\n'+' '.join(map(str,rr))+'\n'+' '.join(map(str,cols))+'\n'+' '.join(map(str,fr))+'\n')
 for row in rows:
  for j in range(106):
   z=row.get(j,R.zero);f.write(str(len(z))+''.join(' '+str(c)+' '+str(m[0])+' '+str(m[1])for m,c in z.items())+'\n')
rep=dict(all_passed=True,status='Candidate reduction over QQ(a,e), no degree225 or full polynomial module claim',D=str(D.as_expr()),normal_num=[[str(z.as_expr())for z in Nnum[i]]for i in light],light_indices=light,rows=len(rows),columns=106,rank_witness=rank,minor=minor,cubic_sources=sources,cubic_contents=contents,cubic_rows_sparse=[[[j,[[list(m),str(c)]for m,c in z.items()]]for j,z in row.items()]for row in cr],cubic_free_monomials=list(map(list,reps[3])),quartic_free_monomials=list(map(list,reps[4])),matrix_terms=sum(len(z)for row in rows for z in row.values()),input=str(path),seconds=time.time()-st)
Path(__file__).with_suffix('.json').write_text(json.dumps(rep,indent=2)+'\n');print('REDUCED',len(rows),106,rank,rep['matrix_terms'],'SECONDS',time.time()-st,flush=True)

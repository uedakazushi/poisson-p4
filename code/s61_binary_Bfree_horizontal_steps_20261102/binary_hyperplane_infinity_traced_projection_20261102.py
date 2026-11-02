from pathlib import Path
from fractions import Fraction as Q
import json,time
st=time.time();ref=json.loads(Path('work/binary_hyperplane_infinity_constant_quadratic_projection_20261102.json').read_text());weights=ref['normal_weights'];quad=[{tuple(m):Q(c)for m,c in row}for row in ref['quadrics']]
def mons(w,i=0,m=()):
 if i==8:
  if not w:yield m
  return
 for k in range(w//weights[i]+1):yield from mons(w-k*weights[i],i+1,m+(k,))
def addto(r,other,c):
 for j,v in other.items():
  z=r.get(j,Q(0))+c*v
  if z:r[j]=z
  elif j in r:del r[j]
results={}
for w in [2,3,4]:
 ms=[m for m in mons(w)if not(w==4 and m[5]==2)];ix={m:j for j,m in enumerate(ms)};shift=list(mons(w-2));original=[];labels=[]
 for k,q in enumerate(quad):
  for sh in shift:original.append({ix[tuple(a+b for a,b in zip(m,sh))]:v for m,v in q.items()});labels.append([k,list(sh)])
 n=len(ms);nr=len(original);A=[r|{n+i:Q(1)}for i,r in enumerate(original)];piv=[];rank=0
 for j in range(n):
  candidates=[i for i in range(rank,nr)if A[i].get(j)]
  if not candidates:continue
  i=min(candidates,key=lambda i:len(A[i]));A[rank],A[i]=A[i],A[rank];c=A[rank][j];A[rank]={k:v/c for k,v in A[rank].items()}
  for i in range(nr):
   if i!=rank and A[i].get(j):addto(A[i],A[rank],-A[i][j])
  piv.append(j);rank+=1
  if rank==nr:break
 for r in A:
  prod={}
  for j,c in r.items():
   if j>=n:addto(prod,original[j-n],c)
  assert prod=={j:c for j,c in r.items()if j<n}
 assert rank=={2:6,3:37,4:132}[w]
 results[str(w)]=dict(monomials=list(map(list,ms)),row_labels=labels,rank=rank,pivot_columns=piv,free_columns=[j for j in range(n)if j not in piv],relations=[[[j,str(c)]for j,c in r.items()if j<n]for r in A[:rank]],relation_lifts=[[[j-n,str(c)]for j,c in r.items()if j>=n]for r in A[:rank]],left_syzygies=[[[j-n,str(c)]for j,c in r.items()if j>=n]for r in A[rank:]])
print('traces',[ (w,r['rank'],len(r['left_syzygies']))for w,r in results.items()],'seconds',time.time()-st);out=Path(__file__).with_suffix('.json');assert not out.exists();out.write_text(json.dumps(dict(all_passed=True,projections=results,seconds=time.time()-st),indent=2)+'\n')

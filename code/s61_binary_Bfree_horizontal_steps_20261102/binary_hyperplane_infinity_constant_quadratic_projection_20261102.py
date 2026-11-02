from pathlib import Path
import json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
start=time.time();weights=[1,1,1,1,1,2,1,1]
def mons(w,i=0,m=()):
 if i==8:
  if w==0:yield m
  return
 for k in range(w//weights[i]+1):yield from mons(w-k*weights[i],i+1,m+(k,))
def m(**kw):return tuple(kw.get(str(i),0)for i in range(8))
def exp(*items):return tuple(items.count(i)for i in range(8))
qs=[{exp(0,0):1},{exp(0,1):1},{exp(0,4):1,exp(0,2):-1},{exp(1,1):1,exp(1,2):-1,exp(0,3):-1},{exp(3,4):1},{exp(4,2):1}]
ms=[t for t in mons(4)if t[5]!=2];ix={t:i for i,t in enumerate(ms)};rows=[]
for q in qs:
 for sh in mons(2):rows.append({ix[tuple(a+b for a,b in zip(t,sh))]:s.QQ(c)for t,c in q.items()})
A=DomainMatrix({i:r for i,r in enumerate(rows)},(len(rows),len(ms)),s.QQ);RR,piv=A.rref();assert len(rows)==174 and len(ms)==238 and len(piv)==132;free=[i for i in range(238)if i not in piv];fi={j:i for i,j in enumerate(free)};raw=RR.to_dod();B={j:{k:s.QQ(1)}for j,k in fi.items()}
for i,j in enumerate(piv):B[j]={fi[k]:-c for k,c in raw[i].items()if k in fi and c}
BM=DomainMatrix(B,(238,106),s.QQ);assert A.matmul(BM).is_zero_matrix
out=dict(all_passed=True,status='Constant quadratic quotient for the newly found infinity normal form; transformation into original model is separate.',normal_variables=['L','W','C','P','E','H','free1','free2'],normal_weights=weights,quadrics=[[[list(t),str(c)]for t,c in q.items()]for q in qs],quartic_monomials=list(map(list,ms)),matrix_shape=[174,238],rank=132,free_columns=free,pivot_columns=list(piv),projection_sparse=[[[j,str(c)]for j,c in sorted(B.get(i,{}).items())]for i in range(238)],all_quadratic_relations_killed=True,seconds=time.time()-start)
p=Path(__file__).with_suffix('.json');assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n');print('PASS132',len(free),'terms',sum(len(r)for r in B.values()),'seconds',time.time()-start)

from pathlib import Path
import json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
start=time.time();ref=json.loads(Path('work/binary_hyperplane_infinity_constant_quadratic_projection_20261102.json').read_text());weights=ref['normal_weights'];qs=[{tuple(t):s.QQ(c)for t,c in row}for row in ref['quadrics']]
def mons(w,i=0,m=()):
 if i==8:
  if w==0:yield m
  return
 for k in range(w//weights[i]+1):yield from mons(w-k*weights[i],i+1,m+(k,))
ms=list(mons(3));ix={t:i for i,t in enumerate(ms)};rows=[]
for q in qs:
 for sh in mons(1):rows.append({ix[tuple(a+b for a,b in zip(t,sh))]:c for t,c in q.items()})
A=DomainMatrix({i:r for i,r in enumerate(rows)},(len(rows),len(ms)),s.QQ);RR,piv=A.rref();assert len(rows)==42 and len(ms)==91 and len(piv)==37;free=[i for i in range(91)if i not in piv];fi={j:i for i,j in enumerate(free)};raw=RR.to_dod();B={j:{k:s.QQ(1)}for j,k in fi.items()}
for i,j in enumerate(piv):B[j]={fi[k]:-c for k,c in raw[i].items()if k in fi and c}
BM=DomainMatrix(B,(91,54),s.QQ);assert A.matmul(BM).is_zero_matrix
out=dict(all_passed=True,status='Constant cubic quotient for infinity quadratic normal form; original transformation separate.',normal_variables=ref['normal_variables'],normal_weights=weights,cubic_monomials=list(map(list,ms)),matrix_shape=[42,91],rank=37,free_columns=free,pivot_columns=list(piv),projection_sparse=[[[j,str(c)]for j,c in sorted(B.get(i,{}).items())]for i in range(91)],all_quadratic_relations_killed=True,seconds=time.time()-start)
p=Path(__file__).with_suffix('.json');assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n');print('PASS37',len(free),'terms',sum(len(r)for r in B.values()),'seconds',time.time()-start)

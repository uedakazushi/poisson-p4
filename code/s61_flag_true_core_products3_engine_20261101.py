import json,time,sys
from pathlib import Path
from itertools import combinations_with_replacement as cwr
import sympy as s
WD=Path(sys.argv[1]);st=time.time();q=json.load(open('reports/S61_FLAG_TRUE_CORE_QUADRATIC_20261101.json'));raw=json.load(open(WD/'flag-b1-universal-quadratic.json'));Q=s.Matrix(q['polynomial_quadratic_basis']).applyfunc(s.sympify);RR=s.Matrix(raw['quadratic_rref']).applyfunc(s.sympify);C=RR[[i for i in range(15)if not any(x.free_symbols for x in RR[i,:])],:];m2=list(cwr(range(10),2));m3=list(cwr(range(10),3));ix={m:i for i,m in enumerate(m3)}
def mult(A):
 B=s.zeros(A.rows*10,220)
 for i in range(A.rows):
  for j,m in enumerate(m2):
   if A[i,j]:
    for k in range(10):B[i*10+k,ix[tuple(sorted(m+(k,)))]]+=A[i,j]
 return B
BC=mult(C);RC,pc=BC.to_DM().rref();RC=RC.to_Matrix()[:len(pc),:];free=[i for i in range(220)if i not in pc];print('constant rows',C.rows,'product rank',len(pc),'free',len(free),'seconds',time.time()-st,flush=True)
B=mult(Q);E=(B[:,free].to_DM()-B[:,list(pc)].to_DM()*RC[:,free].to_DM()).to_Matrix();nz=[i for i in range(E.rows)if any(E[i,:])];E=E[nz,:];print('remaining product matrix',E.shape,'symbols',set().union(*(x.free_symbols for x in E)),'seconds',time.time()-st,flush=True)

E0=E.subs({x:0 for x in E.free_symbols});sel=list(E0.T.to_DM().rref()[1]);pe=list(E0[sel,:].to_DM().rref()[1]);M=E[sel,pe].to_DM().to_field();det=M.domain.to_sympy(M.det());assert det.is_Rational and det
RE=(M.inv()*E[sel,:].to_DM().to_field()).to_Matrix().applyfunc(s.cancel)
assert all(s.denom(x).is_Number for x in RE)
assert(E-E[:,pe]*RE).applyfunc(s.expand)==s.zeros(*E.shape)
enc=lambda M:[[str(x)for x in row]for row in M.tolist()]
(WD/'flag-true-core-cubic-products-global.json').write_text(json.dumps(dict(constant_rank=len(pc),constant_rref=enc(RC),constant_pivots=list(pc),free=free,remaining_rows=nz,remaining_rref=enc(RE),remaining_pivots=list(pe),total_rank=len(pc)+len(pe),minor=dict(rows=sel,columns=pe,value=str(det)),seconds=time.time()-st))+'\n')
print('GLOBAL PRODUCTS PASS',time.time()-st,flush=True)

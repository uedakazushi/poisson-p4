"""Polynomial degree-four basis of the established quadratic and cubic ideal."""
import json,sys,time
from pathlib import Path
from itertools import combinations_with_replacement as cwr
import sympy as s
from sympy.polys.matrices import DomainMatrix as DM

st=time.time();wd=Path(sys.argv[1])
q=json.load(open('reports/S61_FLAG_TRUE_CORE_QUADRATIC_20261101.json'))
proof3=json.load(open('reports/S61_FLAG_TRUE_CORE_CUBIC_20261101.json'))
p,t=proof3['quadratic_product_space'],proof3['cubic_space']
raw=json.load(open(wd/'flag-b1-universal-quadratic.json'))
Q=s.Matrix(q['polynomial_quadratic_basis']).applyfunc(s.sympify)
RR=s.Matrix(raw['quadratic_rref']).applyfunc(s.sympify)
C=RR[[i for i in range(15)if not any(x.free_symbols for x in RR[i,:])],:]
assert (C-C[:,q['constant_minor']['columns']]*Q).applyfunc(s.expand)==s.zeros(13,55)
ms={d:list(cwr(range(10),d))for d in [1,2,3,4]}
free3=[p['free'][i]for i in range(120)if i not in p['remaining_pivots']]
T=s.zeros(5,220);TC=s.Matrix(t['rref']).applyfunc(s.sympify)
for k,j in enumerate(free3):T[:,j]=TC[:,k]
def mult(A,d,e):
 ix={m:j for j,m in enumerate(ms[d+e])};out={}
 for i in range(A.rows):
  for j in range(A.cols):
   if A[i,j]:
    for k,m in enumerate(ms[e]):out[i*len(ms[e])+k,ix[tuple(sorted(ms[d][j]+m))]]=A[i,j]
 return s.MutableSparseMatrix(A.rows*len(ms[e]),len(ix),out)
BC=mult(C,2,2);RC,pc=BC.to_DM().rref();RC=RC.to_Matrix()[:len(pc),:]
assert len(pc)==440;free=[j for j in range(715)if j not in pc]
B=mult(Q,2,2).col_join(mult(T,3,1))
E=(B[:,free].to_DM()-B[:,list(pc)].to_DM()*RC[:,free].to_DM()).to_Matrix()
nz=[i for i in range(E.rows)if any(E[i,:])];E=E[nz,:]
E0=E.subs({x:0 for x in E.free_symbols});sel=list(E0.T.to_DM().rref()[1])
cols=list(E0[sel,:].to_DM().rref()[1]);assert len(sel)==64
M=E[sel,cols];vs=sorted(E.free_symbols,key=str);K=s.QQ.poly_ring(*vs)
M0=M.subs({x:0 for x in vs}).to_DM();M0i=M0.inv()
N=M0i.convert_to(K)*M.to_DM().convert_to(K)-DM.eye(64,K)
P=DM.eye(64,K);Inv=P;power_nnz=[]
for i in range(1,5):
 P=P*N;power_nnz.append(len(P.to_dok()))
 if P.is_zero_matrix:break
 Inv=Inv+P.scalarmul(K((-1)**i))
assert P.is_zero_matrix and i==4
MI=Inv*M0i.convert_to(K);assert MI*M.to_DM().convert_to(K)==DM.eye(64,K)
RE=MI*E[sel,:].to_DM().convert_to(K)
assert(E.to_DM().convert_to(K)-E[:,cols].to_DM().convert_to(K)*RE).is_zero_matrix
enc=lambda A:[[str(x)for x in row]for row in (A.to_Matrix()if isinstance(A,DM)else A).tolist()]
result=dict(constant_rref=enc(RC),constant_pivots=list(pc),free=free,
            remaining_rref=enc(RE),remaining_pivots=cols,remaining_rows=nz,
            selected_rows=sel,remaining_inverse=enc(MI),
            minor_value=str(M0.domain.to_sympy(M0.det())),
            nilpotence_index=4,power_nnz=power_nnz,total_rank=504,seconds=time.time()-st)
(wd/'flag-true-core-fourth-products-global.json').write_text(json.dumps(result)+'\n')
print('ALL PARAMETER FOURTH LOWER SPACE PASS',time.time()-st,flush=True)

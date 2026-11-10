#!/usr/bin/env python3
"""Exact 5-row witness for an ambient P4 Poisson quadratic obstruction.

All 150 quadratic bivector coefficients are tested, and the full 350x150
projective Jacobi linearization is ranked over QQ. No previous code imports.
"""
import argparse,itertools,json
from collections import defaultdict
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=Path('work/ambient_obstruction.json'));a=p.parse_args()
assert s.__version__=='1.14.0'
x=s.symbols('u0 u1 z w0 w1');u0,u1,z,w0,w1=x;E=s.Matrix(x)
q=u0**2+2*u1*w0+2*z*w1;G=s.hessian(q,E)/2
pairs=list(itertools.combinations(range(5),2));trip=list(itertools.combinations(range(5),3));quad=list(itertools.combinations(range(5),4))
def sgn(t):return (-1)**sum(t[i]>t[j]for i in range(len(t))for j in range(i+1,len(t)))
xi=(u0,u1,-2*z);gamma={(1,2):xi[0],(0,2):-xi[1],(0,1):xi[2]}
P=s.zeros(5)
for i,j in pairs:
 rest=tuple(k for k in range(5)if k not in (i,j))
 P[i,j]=s.expand(sgn((i,j,*rest))*sum((-1)**r*s.diff(q,x[k])*gamma.get(tuple(rest[t]for t in range(3)if t!=r),0)for r,k in enumerate(rest)));P[j,i]=-P[i,j]
C=s.Matrix([0,u1,-z,-w0,w1]);A=C.jacobian(E)
assert A.T*G+G*A==s.zeros(5)
K=(q*A*G.inv()/2+(E*C.T-C*E.T)/2).applyfunc(s.expand)
def jac(B):return {ijk:s.expand(sum(B[i,l]*s.diff(B[j,k],x[l])+B[j,l]*s.diff(B[k,i],x[l])+B[k,l]*s.diff(B[i,j],x[l])for l in range(5)))for ijk in trip for i,j,k in [ijk]}
J0=jac(P);JK=jac(K);J1=jac(P+K)
assert all(v==0 for v in J0.values()) and all(J1[t]==JK[t]for t in trip)
assert any(v for v in JK.values())
def mons(degree):
 out=[]
 for t in itertools.combinations_with_replacement(range(5),degree):
  ex=[0]*5
  for i in t:ex[i]+=1
  out.append(tuple(ex))
 return out
m2,m3,m4=(mons(d)for d in (2,3,4))
rows=[(t,e)for t in trip for e in m3];ri={v:i for i,v in enumerate(rows)}
prows=[(t,e)for t in quad for e in m4];rj={v:i for i,v in enumerate(prows)}
D={(i,i):{}for i in range(5)}
for i,j in pairs:
 terms={k:v for k,v in s.Poly(P[i,j],*x).terms()if v}
 D[i,j]=terms;D[j,i]={k:-v for k,v in terms.items()}
def der(poly,i):return {tuple(v-1 if j==i else v for j,v in enumerate(k)):a*k[i]for k,a in poly.items()if k[i]}
def mul(a,b):
 out=defaultdict(lambda:s.S.Zero)
 for ea,va in a.items():
  for eb,vb in b.items():out[tuple(p+q for p,q in zip(ea,eb))]+=va*vb
 return out
L=s.MutableSparseMatrix(350,150,{})
for col,(ij,m) in enumerate(itertools.product(pairs,m2)):
 i,j=ij;R={(u,v):{}for u in range(5)for v in range(5)}
 R[i,j]={m:s.S.One};R[j,i]={m:-s.S.One}
 for t in trip:
  i,j,k=t
  for u,v,w in ((i,j,k),(j,k,i),(k,i,j)):
   for l in range(5):
    for e,c in mul(D[u,l],der(R[v,w],l)).items():L[ri[t,e],col]+=c
    for e,c in mul(R[u,l],der(D[v,w],l)).items():L[ri[t,e],col]+=c
W=s.MutableSparseMatrix(350,350,{})
for i,(t,e)in enumerate(rows):
 for k in range(5):
  four=tuple(sorted(set(t)|{k}))
  if len(four)==4:
   v=list(e);v[k]+=1;W[rj[four,tuple(v)],i]+=(-1)**four.index(k)
M=W*L
b=-W*s.Matrix([s.Poly(JK[t],*x).coeff_monomial(s.prod(y**j for y,j in zip(x,e)))for t,e in rows])
rank=DomainMatrix.from_Matrix(M).convert_to(s.QQ).rank()
aug=DomainMatrix.from_Matrix(M.row_join(b)).convert_to(s.QQ).rank()
assert (rank,aug)==(95,96)
witness=[((0,1,2,3),(3,0,1,0,0),2),((0,1,2,3),(1,1,1,1,0),1),((0,1,2,3),(1,0,2,0,1),1),((0,1,2,4),(3,1,0,0,0),4),((0,1,2,4),(1,2,0,1,0),10)]
lam=s.zeros(1,350)
for t,e,c in witness:lam[0,rj[t,e]]=c
assert lam*M==s.zeros(1,150)
assert (lam*b)[0]==9
sel=[rj[t,e]for t,e,c in witness]
nonzero=[i for i in range(150)if any(M[j,i]for j in sel)]
assert nonzero==[4,18,60,67,71,88]
res={'date':'2026-11-10','version':'v8','status':'PASS','sympy':s.__version__,'q':str(q),
'raw_variables':150,'projective_jacobi_equations':350,'linearization_rank':rank,
'augmented_rank':aug,'projective_tangent_dimension':150-rank-24-1,
'nonzero_basis_columns_in_five_rows':nonzero,
'five_row_six_column_matrix':[[int(v) for v in row] for row in M.extract(sel,nonzero).tolist()],
'left_annihilator':[c for t,e,c in witness],
'obstruction_on_negative_Jacobi':'9',
'witness_rows':[{'indices':list(t),'monomial_exponents':list(e)}for t,e,c in witness],
'JK_nonzero':{''.join(map(str,t)):str(s.factor(v))for t,v in JK.items()if v},
'claim':'The embedded c-direction is first-order Poisson but is obstructed at order two in the full projective P4 Poisson scheme. Ambient smoothness fails; ambient nonreducedness is not established.'}
a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(res,indent=2)+'\n')
print('PASS full ambient projective obstruction; ranks',rank,aug,'witness',str((lam*b)[0]))

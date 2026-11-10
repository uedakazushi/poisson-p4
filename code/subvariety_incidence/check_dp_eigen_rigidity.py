#!/usr/bin/env python3
"""Universal exact minor at the five discriminant eigen-hyperplanes of a smooth (2,2) del Pezzo.

Coordinates x0..x4, q0=sum x_i^2, q1=sum lambda_i*x_i^2,
lambda0=0 and h=x0. The ten-by-ten mixed-Schouten minor is a
nonzero product of eigenvalue differences. The Japanese TeX proves why
injectivity here implies all twenty extension directions are rigid.
"""
import argparse
import itertools as it
import json
from pathlib import Path
import sympy as s

parser=argparse.ArgumentParser()
parser.add_argument("--output",type=Path,default=Path("work/dp_eigen_rigidity.json"))
parser.add_argument("--check",type=Path)
args=parser.parse_args()
assert s.__version__=="1.14.0"
x=s.symbols("x0:5")
ell=s.symbols("l1:5")
lam=(s.S.Zero,*ell)
pairs=list(it.combinations(range(5),2))
triples=list(it.combinations(range(5),3))
def signature(t):
 return (-1)**sum(t[i]>t[j]for i in range(len(t))for j in range(i+1,len(t)))
P=s.zeros(5)
for i,j in pairs:
 if i==0:continue
 k,l=sorted(set(range(1,5))-{i,j})
 P[i,j]=signature((i,j,0,k,l))*(lam[l]-lam[k])*x[k]*x[l]
 P[j,i]=-P[i,j]
assert all(s.expand(sum(s.diff(P[i,j],x[i])for i in range(5)))==0 for j in range(5))
q0=sum(v*v for v in x)
q1=sum(lam[i]*x[i]**2 for i in range(5))
assert all(s.expand(sum(P[i,j]*s.diff(q,x[j])for j in range(5)))==0
           for q in (q0,q1)for i in range(5))
rows=[]
desc=[]
for I in triples:
 i,j,k=I
 for m in range(5):
  v=[]
  for a,b in pairs:
   T=s.zeros(5);T[a,b]=1;T[b,a]=-1
   J=s.expand(sum(T[vv,l]*s.diff(P[ww,zz],x[l])
                  for vv,ww,zz in ((i,j,k),(j,k,i),(k,i,j))
                  for l in range(5)))
   v.append(J.coeff(x[m]))
  rows.append(v)
  desc.append((I,m))
M=s.Matrix(rows)
minor_rows=(3,4,9,19,31,32,33,36,37,41)
minor=s.factor(M.extract(list(minor_rows),list(range(10))).det())
expected=-(ell[0]-ell[1])*(ell[0]-ell[2])*(ell[0]-ell[3])**2*(ell[1]-ell[2])*(ell[1]-ell[3])**2*(ell[2]-ell[3])**3
assert s.expand(minor-expected)==0
specialized=dict(zip(ell,(1,2,3,4)))
assert minor.subs(specialized)==-72
assert M.subs(specialized).rank()==10
report={
 "date":"2026-11-11","version":"v19","status":"PASS","sympy":s.__version__,
 "pencil":"q0=sum_i x_i^2; q1=sum_i lambda_i*x_i^2 (lambda0=0)",
 "anticanonical_section":"h=x0 (one of the five discriminant eigen-hyperplanes)",
 "hypothesis":"lambda1,..,lambda4 distinct from each other and from lambda0",
 "strict_poisson_casimir_check":"Jacobian J(q0,q1,x0)/4 is divergence-free and q0,q1 Casimirs",
 "linear_map":"constant bivector T -> [J(q0,q1,x0)/4,T], linear triple-vector coefficients",
 "matrix_dimensions":[50,10],
 "minor_row_indices_zero_based":list(minor_rows),
 "minor_rows":[[list(desc[i][0]),desc[i][1]]for i in minor_rows],
 "minor_column_pairs":[list(pair)for pair in pairs],
 "minor_factorization":str(s.factor(minor)),
 "minor_at_lambda_1_2_3_4":"-72",
 "universal_rank_of_constant_bivector_map":10,
 "projective_Jacobi_full_extension_linearization_rank":20,
 "scheme_conclusion":"The Jacobian extension is a reduced isolated point of the full 20-parameter Poisson extension scheme for all smooth diagonal pencils and any eigen-hyperplane; other disconnected solutions are not ruled out.",
 "scope":"The implication from constant-bivector injectivity to all 20 extension directions uses q0,q1 coprime and divergence-free Koszul exactness, proved in accompanying Japanese TeX."
}
if args.check is not None:
 assert report==json.loads(args.check.read_text(encoding="utf-8")),"certificate mismatch"
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
print("PASS: universal discriminant minor and local rigidity at all five eigen-hyperplanes")

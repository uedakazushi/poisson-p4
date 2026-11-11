#!/usr/bin/env python3
"""Exact 20 universal coefficient identities forcing constant compatibility to vanish.

For any smooth (2,2) del Pezzo, q0,q1 have five pairwise distinct diagonal
eigenvalues. h=sum_i c_i*x_i is arbitrary nonzero and A arbitrary skew.
The proof of A=0 and the open-closed universal incidence split is in TeX.
"""
import argparse,itertools as it,json
from pathlib import Path
import sympy as s
parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,default=Path('work/dp_universal_jacobian_rigidity.json'))
parser.add_argument('--check',type=Path)
args=parser.parse_args()
assert s.__version__=='1.14.0'
x=s.symbols('x0:5');c=s.symbols('c0:5');lam=s.symbols('l0:5');aa=s.symbols('a0:10')
pairs=list(it.combinations(range(5),2))
A=s.zeros(5)
for a,(i,j)in zip(aa,pairs):A[i,j]=a;A[j,i]=-a
L=s.diag(*lam);M=A*L-L*A;d=A*s.Matrix(c);X=s.Matrix(x)
q0=sum(z*z for z in x);q1=sum(lam[i]*x[i]**2 for i in range(5))
F=s.expand(2*(X.T*M*X)[0])
U=s.expand(2*(X.T*d)[0])
V=s.expand(2*(X.T*L*d)[0])
# Directly evaluate the 5 contractions with dq0 ^ dq1 ^ dh.
for k in range(5):
 actual=0
 for i,j in pairs:
  co=s.det(s.Matrix([[s.diff(q0,x[t])for t in(i,j,k)],
                     [s.diff(q1,x[t])for t in(i,j,k)],
                     [c[t]for t in(i,j,k)]]))
  actual+=A[i,j]*co
 expected=F*c[k]-U*s.diff(q1,x[k])+V*s.diff(q0,x[k])
 assert s.expand(actual-expected)==0
beta=s.Matrix([s.expand(F*c[k]-2*U*lam[k]*x[k]+2*V*x[k])for k in range(5)])
for i,j in pairs:
 db=s.expand(s.diff(beta[j],x[i])-s.diff(beta[i],x[j]))
 expected_i=-4*(lam[j]-lam[i])*(A[i,j]*c[i]+d[j])
 expected_j=4*(lam[j]-lam[i])*(A[i,j]*c[j]-d[i])
 assert s.expand(db.coeff(x[i])-expected_i)==0,(i,j,'xi')
 assert s.expand(db.coeff(x[j])-expected_j)==0,(i,j,'xj')
assert all(s.expand(sum(A[j,i]*c[i]for i in range(5)if i!=j)-d[j])==0 for j in range(5))
report={
 'date':'2026-11-11','version':'v20','status':'PASS','sympy':s.__version__,
 'base_field':'QQ[lambda0:lambda4,c0:c4,Aij,xi]',
 'arbitrary_h':'h=sum_i c_i*x_i, any nonzero c',
 'pencil':'q0=sum x_i^2, q1=sum lambda_i*x_i^2, pairwise distinct lambda_i',
 'constant_bivector_dimension':10,
 'contraction_coefficients_checked':5,
 'closedness_linear_coefficients_checked':20,
 'coeff_xi':'-4*(lambda_j-lambda_i)*(A_ij*c_i+(A*c)_j)',
 'coeff_xj':'4*(lambda_j-lambda_i)*(A_ij*c_j-(A*c)_i)',
 'deduction':'for fixed j, (A*c)_j=sum_(i!=j) A_ji*c_i=4*(A*c)_j; hence A*c=0 and all A_ij*c_i=0; any c_k!=0 then forces all A_ij=0',
 'tangent_result':'every Jacobian extension is reduced isolated, universally in smooth (2,2) S and nonzero h',
 'universal_incidence':'canonical section is open and closed, by unramified neighborhood + separatedness',
 'full_extension_uniqueness':'NOT PROVED; disconnected rank-four solutions may exist',
 'scope':'All coefficient identities are polynomial and checked. The deduction and global incidence splitting use ordinary proofs in current Japanese TeX.'}
if args.check:
 assert report==json.loads(args.check.read_text(encoding='utf-8')),'certificate mismatch'
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print('PASS universal 5 contractions and 20 Jacobian compatibility identities')

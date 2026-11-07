#!/usr/bin/env python3
"""One complete fixed-coframe T01 limit test, with all six magnetic RHS.

This decides only the prescribed hyperplane confluence, never T01 membership.
Periods through order2 are sufficient for its inverse-frame pole of order1.
"""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]

def main():
 started=time.monotonic();source_path=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_STABLE_CUBIC_20261107.json'
 source=json.loads(source_path.read_text())['cases'][1]
 h,q1,q2,x,y,e=s.symbols('h x1 x2 x3 x4 e');xx=(q1,q2,x,y);mm=s.symbols('m0:6');env={str(z):z for z in(h,)+xx+mm}
 ell=s.Matrix([1+e*q1,x+3+e*q2,y-3,x-y-3]);B=ell.jacobian(xx);inverse=B.inv()
 assert s.expand(B.det()+e*e)==0
 assert all(s.Poly(s.expand(e*z),e).degree()<=1 for z in inverse)
 pairs=list(itertools.combinations(range(4),2));aa=s.symbols('a0:18');A=s.zeros(4)
 for k,(i,j) in enumerate(pairs):A[i,j]=sum(aa[6*r+k]*e**r for r in range(3))*ell[i]*ell[j];A[j,i]=-A[i,j]
 candidate=(inverse*A*inverse.T).applyfunc(s.expand)
 target=s.zeros(4)
 for ij,z in source['full_six_magnetic_source'].items():
  i,j=map(int,ij.split(','));assert i>0
  target[i-1,j-1]=s.sympify(z,locals=env).subs(h,1);target[j-1,i-1]=-target[i-1,j-1]
 equations=[];labels=[]
 for i,j in pairs:
  expr=s.Poly(s.expand(e*e*(candidate[i,j]-target[i,j])),e,*xx)
  for ex,co in expr.terms():
   if ex[0]<=2:equations.append(co);labels.append([i,j,ex[0]-2,*ex[1:]])
 matrix,rhs=s.linear_eq_to_matrix(equations,aa)
 rhs_columns=rhs.subs(dict.fromkeys(mm,0)).row_join(s.Matrix.hstack(*[rhs.diff(m)for m in mm]))
 dm=lambda Z:DomainMatrix.from_Matrix(Z).convert_to(s.QQ)
 rank=lambda Z:len(dm(Z).rref()[1])
 ranks=[rank(matrix.row_join(rhs_columns[:,i]))for i in range(7)];r=rank(matrix)
 left=dm(matrix.T).nullspace().to_Matrix();obstruction=left*rhs_columns
 assert left*matrix==s.zeros(left.rows,18)
 witnesses=[]
 for j in range(7):
  rows=[i for i in range(left.rows)if obstruction[i,j]!=0]
  if rows:
   i=rows[0];witnesses.append(dict(rhs='constant'if j==0 else str(mm[j-1]),nonzero_pairing=str(obstruction[i,j]),left_functional={str(k):str(a)for k,a in enumerate(left.row(i))if a},equation_labels={str(k):labels[k]for k,a in enumerate(left.row(i))if a}))
 all_consistent=all(z==r for z in ranks)
 result=dict(task_id='E_multiple_hyperplane',all_passed=True,input_sha256={str(source_path.relative_to(ROOT)):hashlib.sha256(source_path.read_bytes()).hexdigest()},producer_bytes=Path(__file__).read_text(),
  coframe=list(map(str,ell)),homogeneous_coframe=['h','h+e*x1','3*h+x3+e*x2','-3*h+x4','-3*h+x3-x4'],coframe_determinant=str(B.det()),inverse_pole_order=1,
  period_count=6,period_Taylor_max_degree=2,unknown_count=18,coefficient_equation_count=len(equations),coefficient_matrix_rank=r,constant_and_six_magnetic_augmented_ranks=ranks,
  full_six_magnetic_target=source['full_six_magnetic_source'],whole_source_regular_limit_in_this_fixed_coframe=all_consistent,
  exact_left_obstructions=witnesses,all_left_functionals_annihilate_complete_coefficient_matrix=True,
  finite_order_completeness='At an off-divisor affine point, bij={ell_i,ell_j}/(ell_i ell_j) is regular in e for every regular bivector family. Thus six periods are regular. InverseB has pole at most1 and its wedge at most2; negative powers and the constant limit depend only on period Taylor coefficients0,1,2. Higher terms cannot change these equations.',
  scope='A failed test rules out this one fixed hyperplane coframe with arbitrary regular periods. It does not exclude T01 membership, other hyperplane paths, MPS membership or fullsource componenthood.',seconds=time.monotonic()-started)
 dest=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_TORIC_CONFLUENCE_20261107.json'
 if dest.exists():raise FileExistsError(str(dest))
 dest.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:result[k]for k in ['coefficient_matrix_rank','constant_and_six_magnetic_augmented_ranks','whole_source_regular_limit_in_this_fixed_coframe','seconds']}))
if __name__=='__main__':main()

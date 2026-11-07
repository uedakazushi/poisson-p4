#!/usr/bin/env python3
"""Last finite highest-cross entrance for the fixedN degree2 unipotent limit.
This tests the sole degree2 case left by the complete nilpotence-range audit.
"""
from pathlib import Path
import hashlib,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS
from classification_three_hyperplanes_quadric_independent import ambient_matrices
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();path=ROOT/'reports/CLASSIFICATION_GA9_UNIPOTENT_LIMITS_20261107.json';packet=json.loads(path.read_text());case=next(c for c in packet['cases']if c['name']=='LINEAR_W');xx=s.symbols('h x y z w');E=s.Matrix(xx);N=s.Matrix(packet['N_matrix']);NX=N*E;P=s.zeros(5)
 for ij,f in case['canonical_point'].items():
  i,j=map(int,ij.split(','));P[i,j]=s.sympify(f,locals=dict(zip(map(str,xx),xx)));P[j,i]=-P[i,j]
 def vec(A):return s.Matrix([s.Poly(s.expand(A[ij]),*xx).coeff_monomial(ex)for ij,ex in COLUMNS])
 LL=[];CC=[]
 for ij,ex in COLUMNS:
  B=s.zeros(5);B[ij]=s.prod(x**a for x,a in zip(xx,ex));B[ij[1],ij[0]]=-B[ij];LL.append(vec(N*B+B*N.T-sum((NX[i]*B.diff(xx[i])for i in range(5)),s.zeros(5))));div=s.Matrix([sum(s.diff(B[i,j],xx[i])for i in range(5))for j in range(5)]);CC.append(vec(B-(E*div.T-div*E.T)/5))
 L=s.Matrix.hstack(*LL);C=s.Matrix.hstack(*CC);D,V,G,F=ambient_matrices(P,xx,[s.zeros(5)]);A=(L**2*C).col_join(D*L*C);rhs=vec(P).col_join(s.zeros(350,1));dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ);rank=lambda A:len(dm(A).rref()[1]);r,ar=rank(A),rank(A.row_join(rhs));result=dict(task_id='E_reduced_nonnormal',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [Path(__file__),path,ROOT/'code/classification_rank3_mixed_component_review.py',ROOT/'code/classification_three_hyperplanes_quadric_independent.py']},producer_bytes=Path(__file__).read_text(),matrix_shape=list(A.shape),rank=r,augmented_rank=ar,necessary_linear_conditions_consistent=r==ar,scope='Onlyhighestcross necessarylinearconditions forfixedN degree2 actualconjugation limit. NoPoissonQ/membership fromconsistententrance. Togetherwithfullrange imageaudit andsimplepencilobstructions thisclosesfixedNmodel onlyifinconsistent.',seconds=time.monotonic()-started)
 if ar>r:
  left=dm(A.T).nullspace().to_Matrix();i=next(i for i in range(left.rows)if(left*rhs)[i]!=0);f=left.row(i);assert f*A==s.zeros(1,150);result.update(left_functional={str(i):str(a)for i,a in enumerate(f)if a},nonzero_pairing=str((f*rhs)[0]))
 else:result['next_route']='StopthisfixedNroutewithoutsolvinglargeresidualquadraticQideal; useotherwhole-sourcegeometry.'
 out=ROOT/'reports/CLASSIFICATION_GA9_UNIPOTENT_QUADRATIC_LIMIT_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in ['rank','augmented_rank','necessary_linear_conditions_consistent','seconds']}))
if __name__=='__main__':main()

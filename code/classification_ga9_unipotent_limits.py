#!/usr/bin/env python3
"""Finite actual unipotent-orbit-limit entrance at three fullGa9c0 sources.

No jets: epsilon^k exp(N/epsilon)Q is a finite polynomial of bounded degree,
and it is Poisson for all nonzero epsilon exactly when Q is Poisson.
"""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,minor_certificate
from classification_three_hyperplanes_quadric_independent import ambient_matrices
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();xx=s.symbols('h x y z w');env=dict(zip(map(str,xx),xx));N=s.zeros(5);N[1,2]=1;N[3,0]=1;E=s.Matrix(xx);NX=N*E
 def vector(P):return s.Matrix([s.Poly(s.expand(P[ij]),*xx).coeff_monomial(ex)for ij,ex in COLUMNS])
 def canonical(P):
  div=s.Matrix([sum(s.diff(P[i,j],xx[i])for i in range(5))for j in range(5)]);return(P-(E*div.T-div*E.T)/5).applyfunc(s.expand)
 columns=[];projection=[]
 for ij,ex in COLUMNS:
  P=s.zeros(5);P[ij]=s.prod(x**a for x,a in zip(xx,ex));P[ij[1],ij[0]]=-P[ij];columns.append(vector((N*P+P*N.T-sum((NX[i]*P.diff(xx[i])for i in range(5)),s.zeros(5))).applyfunc(s.expand)));projection.append(vector(canonical(P)))
 L=s.Matrix.hstack(*columns);C=s.Matrix.hstack(*projection);dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ);rank=lambda A:len(dm(A).rref()[1]);assert L*C==C*L;assert C*C==C and rank(C)==126
 powers=[];T=L.copy()
 for k in range(1,12):
  A=T*C;r=rank(A);powers.append((k,A,r));
  if r==0:break
  T=L*T
 assert powers[-1][2]==0;records=[];inputs=[Path(__file__),ROOT/'code/classification_rank3_mixed_component_review.py',ROOT/'code/classification_three_hyperplanes_quadric_independent.py']
 for label in ['QUADRATIC_W_POWER1','QUADRATIC_W_POWER3','LINEAR_W']:
  path=ROOT/('reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL211_ZERO_C_'+label+'_NORMALFORM_20261107.json');inputs.append(path);packet=json.loads(path.read_text())['source_magnetic_normalform'];P=s.zeros(5)
  for ij,f in packet['canonical_brackets'].items():
   i,j=map(int,ij.split(','));P[i,j]=s.sympify(f,locals=env);P[j,i]=-P[i,j]
  P=canonical(P);den=s.ilcm(1,*[a.q for f in P for a in s.Poly(f,*xx).coeffs()]);P*=den;v=vector(P);assert L*v==s.zeros(150,1);data=[]
  for k,A,r in powers:
   ar=rank(A.row_join(v));item=dict(k=k,image_rank=r,augmented_rank=ar,point_in_image=ar==r)
   if ar!=r:
    left=dm(A.T).nullspace().to_Matrix();which=next(i for i in range(left.rows)if(left*v)[i]!=0);f=left.row(which);assert f*A==s.zeros(1,150);item.update(left_functional={str(i):str(a)for i,a in enumerate(f)if a},nonzero_pairing=str((f*v)[0]))
   data.append(item)
  D,V,G,F=ambient_matrices(P,xx,[s.zeros(5)]);A=(L*C).col_join(D*C);rhs=v.col_join(s.zeros(350,1));r,ar=rank(A),rank(A.row_join(rhs));pencil=dict(linear_system_shape=list(A.shape),rank=r,augmented_rank=ar,actual_simple_unipotent_pencil_linear_entrance=ar==r)
  if ar!=r:
   left=dm(A.T).nullspace().to_Matrix();which=next(i for i in range(left.rows)if(left*rhs)[i]!=0);f=left.row(which);pencil.update(left_functional={str(i):str(a)for i,a in enumerate(f)if a},nonzero_pairing=str((f*rhs)[0]))
  else:
   sol,pars=A.gauss_jordan_solve(rhs);image=C*sol;actualvars=sorted(set().union(*(a.free_symbols for a in image)),key=str);pencil.update(canonical_linear_solution={str(i):str(a)for i,a in enumerate(image)if a},free_actual_parameters=list(map(str,actualvars)),remaining_Poisson_problem='All350 Jacobi ofQ muststill be solved; no actualPoissonfamily or membership claimed.')
  record=dict(name=label,canonical_point={str(i)+','+str(j):str(P[i,j])for i,j in PAIRS if P[i,j]!=0},powers=data,simple_pencil=pencil);records.append(record);print(label,[(q['k'],q['point_in_image'])for q in data],pencil['rank'],pencil['augmented_rank'],len(pencil.get('free_actual_parameters',[])),flush=True)
 result=dict(task_id='E_reduced_nonnormal',all_passed=True,scope='Finite representation-theoretic necessary entrance for actual unipotent orbit degeneration, no generic componentmembership unless a fullPoissonQ is found. N-invariantpoints tested for the FULL finite nilpotentpower range.',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),N_matrix=[[int(a)for a in row]for row in N.tolist()],canonical_projection_rank=126,nilpotence_index_on_projective_bivectors=powers[-1][0],power_ranks=[dict(k=k,rank=r)for k,A,r in powers],cases=records,seconds=time.monotonic()-started)
 out=ROOT/'reports/CLASSIFICATION_GA9_UNIPOTENT_LIMITS_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print('PASS',result['seconds'])
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Whole sixmag doubleC fixedconic T02 confluence, complete pole2 linearclass."""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,jacobi
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();xx=s.symbols('h q1 q2 x y');h,q1,q2,x,y=xx;e=s.Symbol('epsilon');mag=s.symbols('m0:6');sourceparams=s.symbols('p0:4')+(s.Symbol('t0'),);forms=s.Matrix([x,x+3*h+y,h+e*q1,h+e*q2,x+h]);T=forms.jacobian(s.Matrix(xx));assert s.factor(T.det()) in [e*e,-e*e];inverse=T.inv();E=s.Matrix(xx)
 def canonical(P):
  div=s.Matrix([sum(s.diff(P[i,j],xx[i])for i in range(5))for j in range(5)]);return(P-(E*div.T-div*E.T)/5).applyfunc(s.expand)
 def vector(P):return s.Matrix([s.Poly(s.expand(P[ij]),*xx).coeff_monomial(ex)for ij,ex in COLUMNS])
 toricpath=ROOT/'reports/MPS_S61_RAW_CATALOGUE_20261103_v62/T02.json';sourcepath=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json';toric=json.loads(toricpath.read_text());assert toric['projective_Jacobi_identically_zero'];env=dict(zip(['x'+str(i)for i in range(5)],forms));env.update({str(p):p for p in sourceparams});P=s.zeros(5)
 for ij,text in toric['homogeneous_bivector'].items():
  i,j=map(int,ij.split(','));P[i,j]=s.sympify(text,locals=env);P[j,i]=-P[i,j]
 pulled=canonical(inverse*P*inverse.T);columns=[];labels=[]
 for param in sourceparams:
  basis=vector(pulled.diff(param))
  for degree in range(3):
   values=[s.Poly(s.cancel(e**(2+degree)*co),e)for co in basis];columns.append(s.Matrix([f.coeff_monomial(e**power)for power in range(3)for f in values]));labels.append(dict(parameter=str(param),Taylor_degree=degree))
 A=s.Matrix.hstack(*columns);source=json.loads(sourcepath.read_text())['cases'][2];original=s.zeros(5);env=dict(zip(['h','x1','x2','x3','x4'],xx));env.update({str(m):m for m in mag})
 for ij,text in source['full_six_magnetic_source'].items():
  i,j=map(int,ij.split(','));original[i,j]=s.sympify(text,locals=env);original[j,i]=-original[i,j]
 target=vector(canonical(original));rhs=s.zeros(300,7).col_join(s.Matrix.hstack(target.subs(dict.fromkeys(mag,0)),*[target.diff(m)for m in mag]));dm=lambda M:DomainMatrix.from_Matrix(M).convert_to(s.QQ);r=len(dm(A).rref()[1]);ar=len(dm(A.row_join(rhs)).rref()[1]);print('full450x15',r,ar,round(time.monotonic()-started,3),flush=True)
 result=dict(task_id='E_multiple_hyperplane',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [Path(__file__),toricpath,sourcepath,ROOT/'code/classification_rank3_mixed_component_review.py']},producer_bytes=Path(__file__).read_text(),whole_six_magnetic_target=True,source_component='T02',source_linear_forms=list(map(str,forms)),source_frame_determinant=str(s.factor(T.det())),parameter_Taylor_degree=2,priori_pole_bound=2,complete_fixedframe_regularlimit_class=True,completeness_proof='AtgenericfixedoffDpoint allfive sourcecoordinateforms at0 nonzero. Sourceperiods are recoveredfromnonedge canonicalbrackets bydividingnonzerocoordinateproducts; edge recoveredfrom01 afterperiods. Henceanyregularpullback hasregularp,t. Inverseframepole1,wedgepole2 meansonlyTaylor0..2 influencepoles/endpoint. Thisproof isforTHISfixedframe, notallPGLarcs.',linear_shape=list(A.shape),rank=r,all_six_rhs_augmented_rank=ar,linear_parameter_labels=labels,complete_matrix=[[str(a)for a in row]for row in A.tolist()],complete_rhs=[[str(a)for a in row]for row in rhs.tolist()],scope='OnecompletefixedconicT02collisionclass, fullEulerquotient150 andall6magnetic. Negativeexcludesonlythisframe, notT02/MPSmembership.',seconds=time.monotonic()-started)
 if ar>r:
  left=dm(A.T).nullspace(divide_last=True).to_Matrix();i,j=next((i,j)for i in range(left.rows)for j in range(rhs.cols)if(left.row(i)*rhs[:,j])[0]!=0);row=left.row(i);assert row*A==s.zeros(1,15);result.update(fixedframe_collision_impossible=True,left_functional={str(k):str(a)for k,a in enumerate(row)if a!=0},target_rhs_column=j,nonzero_pairing=str((row*rhs[:,j])[0]))
 else:
  solutions=[]
  for j in range(rhs.cols):
   sol,free=A.gauss_jordan_solve(rhs[:,j]);solutions.append(sol.subs(dict.fromkeys(list(free),0)))
  solution=s.Matrix.hstack(*solutions)*s.Matrix([1,*mag]);bind={param:s.expand(sum(solution[3*i+d]*e**d for d in range(3)))for i,param in enumerate(sourceparams)};family=pulled.subs(bind,simultaneous=True).applyfunc(s.cancel);assert all(not s.denom(f).has(e)for f in family);assert(family.subs(e,0)-canonical(original)).applyfunc(s.expand)==s.zeros(5);sparse={ij:dict(s.Poly(s.expand(family[ij]),*xx).terms())for ij in PAIRS if family[ij]!=0};assert jacobi(sparse).applyfunc(s.expand)==s.zeros(350,1)
  result.update(actual_whole_source_T02_family=True,full_parameter_substitutions={str(p):str(f)for p,f in bind.items()},whole_homogeneous_family={str(i)+','+str(j):str(family[i,j])for i,j in PAIRS if family[i,j]!=0},all350_Jacobi_zero=True,all_source_endpoint_exact_modEuler=True,entire_double_C_T02_containment=True)
 result['seconds']=time.monotonic()-started;out=ROOT/'reports/CLASSIFICATION_NIL22_DOUBLE_C_T02_CONFLUENCE_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print('PASS',result['seconds'])
if __name__=='__main__':main()

#!/usr/bin/env python3
"""The entire double_D source is a regular limit in the accepted CN31_25."""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from classification_rank3_mixed_component_review import PAIRS,jacobi
ROOT=Path(__file__).resolve().parents[1]
def main():
 start=time.monotonic();h,q1,q2,x,y,e=s.symbols('h q1 q2 x y e');xx=(h,q1,q2,x,y);aff=xx[1:];mm=s.symbols('m0:6');a,u=s.symbols('a u');params=mm+(a,u);env={str(z):z for z in xx+params};env.update(dict(zip(['x1','x2','x3','x4'],aff)))
 source_path=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json';old_path=ROOT/'reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json'
 source=json.loads(source_path.read_text())['cases'][3];target=s.zeros(4)
 for ij,f in source['full_six_magnetic_source'].items():
  i,j=map(int,ij.split(','));target[i-1,j-1]=s.sympify(f,locals=env).subs(h,1);target[j-1,i-1]=-target[i-1,j-1]
 X,Y,Z,W=old=s.symbols('X Y Z W');ns=s.symbols('n0:15');local={str(z):z for z in ns};local.update(dict(zip(['x1','x2','x3','x4'],old)));P=s.zeros(4)
 for (i,j),f in zip(itertools.combinations(range(4),2),json.loads(old_path.read_text())['cases'][0]['full_six_brackets']):P[i,j]=s.sympify(f,locals=local);P[j,i]=-P[i,j]
 M0,M1,M2,M3,M4,M5=mm
 periods={ns[0]:-M0/e-u/e**2,ns[1]:-M3/e,ns[2]:-M4/e,ns[3]:-M5/e-(5*a/2+1)/e**2,ns[4]:-M1/e,ns[5]:-u/e,ns[6]:0,ns[7]:-2/e,ns[8]:-M2/e-(3+(a+s.Rational(1,2))*u)/e**2,ns[9]:0,ns[10]:a,ns[11]:0,ns[12]:-(2*a+1)*u/2,ns[13]:0,ns[14]:-1-2*a}
 # The full n7-open Jacobi branch, including all its compatibility equations.
 assert s.cancel(periods[ns[12]]+(2*a+1)*periods[ns[5]]/periods[ns[7]])==0
 assert s.cancel(periods[ns[13]]+2*(3*a+1)*periods[ns[6]]/(3*periods[ns[7]]))==0
 B=s.Matrix([[0,1,0,-1/e],[e,0,0,0],[0,0,1,0],[0,0,0,1]]);assert B.det()==-e
 family=(B*P*B.T).subs(dict(zip(old,[q2/e,q1+y/e,x,y])),simultaneous=True).subs(periods,simultaneous=True).applyfunc(s.cancel)
 assert (family.subs(e,0)-target).applyfunc(s.expand)==s.zeros(4)
 for z in family:s.Poly(z,*params,*aff,e);assert s.Poly(z,*aff).total_degree()<=2
 strict=[s.expand(sum(family[i,l]*s.diff(family[j,k],aff[l])+family[j,l]*s.diff(family[k,i],aff[l])+family[k,l]*s.diff(family[i,j],aff[l])for l in range(4)))for i,j,k in itertools.combinations(range(4),3)];assert strict==[0]*4
 pf=s.factor(family[0,1]*family[2,3]-family[0,2]*family[1,3]+family[0,3]*family[1,2]);assert s.expand(pf+x*x*(y+u/2))==0
 N=s.Matrix([s.expand(sum(s.diff(family[i,j],aff[i])for i in range(4)))for j in range(4)]);assert N==s.Matrix([x,y+e*q1,0,0])
 Q=s.zeros(5)
 for i,j in itertools.combinations(range(4),2):Q[i+1,j+1]=s.expand(h*h*family[i,j].subs({z:z/h for z in aff},simultaneous=True));Q[j+1,i+1]=-Q[i+1,j+1]
 sparse={ij:dict(s.Poly(s.expand(Q[ij]),*xx).terms())for ij in PAIRS if Q[ij]!=0};assert jacobi(sparse).applyfunc(s.expand)==s.zeros(350,1)
 encode=lambda A:{str(i)+','+str(j):str(s.factor(A[i,j]))for i,j in itertools.combinations(range(A.rows),2)if A[i,j]!=0}
 inputs=[source_path,old_path,Path(__file__),ROOT/'code/classification_rank3_mixed_component_review.py']
 result=dict(task_id='E_multiple_hyperplane',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),all_eight_original_parameters=list(map(str,params)),source_case_index=3,
  actual_GL_matrix=[[str(z)for z in row]for row in B.tolist()],GL_determinant=str(B.det()),whole_old_source_Laurent_parameters={str(z):str(s.factor(f))for z,f in periods.items()},
  complete_n7_open_compatibilities_checked=True,n7_minus_two_over_e_nonzero=True,whole_regular_affine_family=encode(family),whole_regular_homogeneous_family=encode(Q),original_endpoint_equal=True,complete_strict_and_projective350_Jacobi_zero=True,
  actual_affine_Pfaffian=str(pf),whole_modular=list(map(str,N)),whole_source_in_accepted_CN31_25=True,new_component_count=0,
  proof='For every fixed eightparameter choice, each e!=0 tensor is an actualGL transform of the complete accepted nil31-image n7-open source. n7=-2/e never vanishes, and its Pfaffian is not identically0 for any u. The tensor is polynomial in e and its exact endpoint is the complete double_D source. Closedness of CN31_25 therefore covers all parameters, with no exceptional a/u/magnetic locus omitted.',independent_review_required=True,seconds=time.monotonic()-start)
 out=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_DOUBLE_D_LIMIT_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_passed=True,whole_double_D_covered=True,seconds=result['seconds'])))
if __name__=='__main__':main()

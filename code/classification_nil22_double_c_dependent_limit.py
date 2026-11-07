#!/usr/bin/env python3
"""Entire doubleC sixmag source has an actual dependent3H+Q confluence."""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from classification_rank3_mixed_component_review import PAIRS,jacobi as projective_jacobi
from classification_three_hyperplanes_quadric_dependent import projective_source
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();xx=s.symbols('h q1 q2 x y');h,q1,q2,x,y=xx;aff=xx[1:];ep=s.Symbol('epsilon');mag=s.symbols('m0:6');R=sum(a*b for a,b in zip(mag,[1,x,y,x*x,x*y,y*y]));Q=x*(y+1)-1;L=x+1;K=1+ep*L;A=-(1+ep)*y-(2+3*ep)-ep*x
 U=s.Matrix([0,-Q]);V=s.Matrix([L*K,-(1+ep)*(y+1)**2-(1+2*ep)*(y+1)-ep]);P=s.zeros(4);P[0,1]=R+A*q1;P[0,2],P[0,3]=U;P[1,2],P[1,3]=V
 for i,j in itertools.combinations(range(4),2):P[j,i]=-P[i,j]
 bracket=V.jacobian(s.Matrix([x,y]))*U-U.jacobian(s.Matrix([x,y]))*V;assert(bracket-A*U).applyfunc(s.expand)==s.zeros(2,1);assert s.expand(s.det(s.Matrix.vstack(U.T,V.T))-L*K*Q)==0
 strict=[s.expand(sum(P[i,l]*s.diff(P[j,k],aff[l])+P[j,l]*s.diff(P[k,i],aff[l])+P[k,l]*s.diff(P[i,j],aff[l])for l in range(4)))for i,j,k in itertools.combinations(range(4),3)];assert strict==[0]*4
 sourcepath=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json';source=json.loads(sourcepath.read_text())['cases'][2];assert source['name']=='nil22_image_double_C';env=dict(zip(['h','x1','x2','x3','x4'],xx));env.update({str(m):m for m in mag});target=s.zeros(5)
 for ij,text in source['full_six_magnetic_source'].items():
  i,j=map(int,ij.split(','));target[i,j]=s.sympify(text,locals=env);target[j,i]=-target[i,j]
 Ph=s.zeros(5)
 for i,j in itertools.combinations(range(4),2):
  Ph[i+1,j+1]=sum(c*h**(2-sum(ex))*s.prod(z**n for z,n in zip(aff,ex))for ex,c in s.Poly(P[i,j],*aff).terms());Ph[j+1,i+1]=-Ph[i+1,j+1]
 assert(Ph.subs(ep,0)-target).applyfunc(s.expand)==s.zeros(5);sparse={ij:dict(s.Poly(s.expand(Ph[ij]),*xx).terms())for ij in PAIRS if Ph[ij]!=0};assert projective_jacobi(sparse).applyfunc(s.expand)==s.zeros(350,1)
 pf=[]
 for missing in range(5):
  i,j,k,l=[a for a in range(5)if a!=missing];pf.append((-1)**missing*(Ph[i,j]*Ph[k,l]-Ph[i,k]*Ph[j,l]+Ph[i,l]*Ph[j,k]))
 D=s.factor(sum(xx[i]*pf[i]for i in range(5)));Qh=x*y+h*x-h*h;Lh=x+h;Kh=(1+ep)*h+ep*x;assert s.expand(D+h*Lh*Kh*Qh)==0
 # Literal all-eight-source recognition, ratherthan factor resemblance.
 X,Y,Z,W=s.symbols('X Y Z W');new=s.Matrix([-ep*(x+1),y+2+ep*(x+1)/(ep+1),q2,-ep*(ep+1)*q1]);B=new.jacobian(s.Matrix(aff));shift=new-B*s.Matrix(aff);assert s.factor(B.det())==-ep**2*(ep+1);back=(B.inv()*(s.Matrix([X,Y,Z,W])-shift)).applyfunc(s.cancel);moved=(B*P*B.T).subs(dict(zip(aff,back)),simultaneous=True).applyfunc(s.cancel)
 known,knownxx,params,f,knownq=projective_source();rh=s.expand(-ep*(ep+1)*R.subs(dict(zip(aff,back)),simultaneous=True));rvalues=[s.Poly(rh,X,Y).coeff_monomial(m)for m in [1,X,Y,X*X,X*Y,Y*Y]];assert s.cancel(rh-sum(v*m for v,m in zip(rvalues,[1,X,Y,X*X,X*Y,Y*Y])))==0;binding=dict(zip(params,[ep*(ep+1),ep+1,*rvalues]));expected=known[1:,1:].subs(knownxx[0],1).subs(binding,simultaneous=True).subs(dict(zip(knownxx[1:],[X,Y,Z,W])),simultaneous=True).applyfunc(s.cancel);assert(moved-expected).applyfunc(s.cancel)==s.zeros(4)
 magmap=s.Matrix(rvalues).jacobian(mag);assert s.factor(magmap.det())==ep**2*(ep+1)**6
 inputs=[Path(__file__),sourcepath,ROOT/'code/classification_rank3_mixed_component_review.py',ROOT/'code/classification_three_hyperplanes_quadric_dependent.py',ROOT/'research/classification/THREE_HYPERPLANES_QUADRIC_DEPENDENT.md',ROOT/'reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_DEPENDENT_REVIEW_20261107.json',ROOT/'reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_DEPENDENT_SOURCE_20261107.json'];result=dict(task_id='E_multiple_hyperplane',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),all_six_magnetic_parameters=list(map(str,mag)),source_family='nil22_image_double_C',source_case_index=2,whole_affine_family={str(i+1)+','+str(j+1):str(s.expand(P[i,j]))for i,j in itertools.combinations(range(4),2)if P[i,j]!=0},whole_homogeneous_family={str(i)+','+str(j):str(s.expand(Ph[i,j]))for i,j in PAIRS if Ph[i,j]!=0},source_endpoint_all6magnetic_exact=True,all4_strict_and350_projective_Jacobi_zero=True,actual_anticanonical_quintic=str(D),generic_three_hyperplanes=[str(h),str(Lh),str(Kh)],quadric=str(Qh),quadric_rank3_determinant=str(s.factor((s.hessian(Qh,(h,x,y))/2).det())),generic_hyperplane_span_rank=2,three_distinct_hyperplanes_condition='epsilon!=0',entire_family_polynomial_no_denominator=True,
  accepted_full_eight_source_coordinate_map=list(map(str,new)),coordinate_jacobian=[[str(a)for a in row]for row in B.tolist()],coordinate_determinant=str(s.factor(B.det())),new_coordinates=['X','Y','Z','W'],old_coordinates=list(map(str,aff)),full_inverse_coordinate_map=list(map(str,back)),accepted_full_source_parameter_substitution={str(p):str(s.factor(v))for p,v in binding.items()},all_six_parameter_magnetic_map=[[str(z)for z in row]for row in magmap.tolist()],magnetic_map_determinant=str(s.factor(magmap.det())),full_eight_source_all6_brackets_exactly_matched=True,source_a0_times_a0_plus_a1=str(s.factor(binding[params[0]]*(binding[params[0]]+binding[params[1]]))),generic_source_chart='epsilon!=0,-1',source_parameter_localizations=[],
  whole_double_C_in_known_finite_union=True,whole_double_C_specific_component='R4_DEPENDENT_THREE_HYPERPLANES_QUADRIC_24',containment_proof='Theactualfamily ispolynomialinall6magandep; endpointisENTIREdoubleC. For ep!=0,-1 theexplicitinvertibleaffinemap andall8parameterbinding identifyitwiththeacceptedENTIRE rank3/fiberHessian0dependent3H+Q source. Hence liesintheaccepted24component. Closednesscovers ep0 foreverymagparameter, noboundaryorperioddropped.',new_component_count=0,independent_review_required=True,seconds=time.monotonic()-started)
 out=ROOT/'reports/CLASSIFICATION_NIL22_DOUBLE_C_DEPENDENT_LIMIT_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(whole_double_C_known_containment=True,exact_full8source_map=True,seconds=result['seconds'])))
if __name__=='__main__':main()

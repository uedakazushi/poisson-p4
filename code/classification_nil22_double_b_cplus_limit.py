#!/usr/bin/env python3
"""Entire doubleB sixmag source is a regular CPLUS_EXT confluence.

No chosen magneticperiod, no jet: all six arbitrary coefficients are constant
along the exact rational Poisson curve, regular on epsilon!=2 including0.
"""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from classification_rank3_mixed_component_review import PAIRS,jacobi as projective_jacobi
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();xx=s.symbols('h q1 q2 x y');h,q1,q2,x,y=xx;coords=xx[1:];eps=s.Symbol('epsilon');mag=s.symbols('m0:6');R=sum(a*b for a,b in zip(mag,[1,x,y,x*x,x*y,y*y]));a=(eps+2)/(eps-2)
 U=s.Matrix([-x*x/4+2/(eps-2),-x*y/2-2*a*x]);V=s.Matrix([-x*y/4-a*x,8-y*y/2+8*eps*x*x/(eps-2)]);qa=y-4+eps*x*x;qb=y+4+2*x*x;co=-y/4-a
 assert (U.jacobian(s.Matrix([x,y]))*V-V.jacobian(s.Matrix([x,y]))*U+co*U).applyfunc(s.cancel)==s.zeros(2,1)
 assert s.cancel(s.det(s.Matrix.vstack(U.T,V.T))-qa*qb/(2-eps))==0
 P=s.zeros(4);P[0,1]=R+co*q1;P[0,2],P[0,3]=U;P[1,2],P[1,3]=V
 for i,j in itertools.combinations(range(4),2):P[j,i]=-P[i,j]
 strict=[s.cancel(sum(P[i,l]*s.diff(P[j,k],coords[l])+P[j,l]*s.diff(P[k,i],coords[l])+P[k,l]*s.diff(P[i,j],coords[l])for l in range(4)))for i,j,k in itertools.combinations(range(4),3)];assert not any(strict)
 source_path=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json';source=json.loads(source_path.read_text())['cases'][1];assert source['name']=='nil22_image_double_B';original=s.zeros(5);env=dict(zip(['h','x1','x2','x3','x4'],xx));env.update({str(m):m for m in mag})
 for ij,text in source['full_six_magnetic_source'].items():
  i,j=map(int,ij.split(','));original[i,j]=s.sympify(text,locals=env);original[j,i]=-original[i,j]
 homogeneous=s.zeros(5)
 for i,j in itertools.combinations(range(4),2):
  homogeneous[i+1,j+1]=sum(c*h**(2-sum(ex))*s.prod(z**n for z,n in zip(coords,ex))for ex,c in s.Poly(P[i,j],*coords).terms());homogeneous[j+1,i+1]=-homogeneous[i+1,j+1]
 assert (homogeneous.subs(eps,0)-original).applyfunc(s.expand)==s.zeros(5)
 sparse={ij:dict(s.Poly(s.expand(homogeneous[ij]),*xx).terms())for ij in PAIRS if homogeneous[ij]!=0};all350=projective_jacobi(sparse);assert all350.applyfunc(s.cancel)==s.zeros(350,1)
 pf=[]
 for missing in range(5):
  i,j,k,l=[a for a in range(5)if a!=missing];pf.append((-1)**missing*(homogeneous[i,j]*homogeneous[k,l]-homogeneous[i,k]*homogeneous[j,l]+homogeneous[i,l]*homogeneous[j,k]))
 D=s.factor(sum(xx[i]*pf[i]for i in range(5)));QA=h*y-4*h*h+eps*x*x;QB=h*y+4*h*h+2*x*x;assert s.cancel(D+h*QA*QB/(2-eps))==0
 # All generic fibers have exactly shared_tangent_generic configuration.
 X,Y,Z=s.symbols('X Y Z');change={h:X,x:Y,y:-2*Z-4*X};standard_q=s.expand(QB.subs(change,simultaneous=True)/-2);standard_other=s.expand(QA.subs(change,simultaneous=True)/-eps);assert standard_q==X*Z-Y*Y;assert s.cancel(standard_other-(2*X*Z/eps-Y*Y+8*X*X/eps))==0
 normals=[]
 for Q in [QA,QB]:
  A=s.hessian(Q,(h,x,y))/2;normals.append(dict(quadratic_matrix=[[str(a)for a in row]for row in A.tolist()],determinant=str(s.factor(A.det())),restriction_at_H=str(s.expand(Q.subs(h,0)))))
 assert all(s.factor(s.hessian(Q,(h,x,y)).det())!=0 for Q in [QA,QB])
 inputs=[Path(__file__),source_path,ROOT/'code/classification_rank3_mixed_component_review.py',ROOT/'research/classification/TWO_QUADRICS_PLANE_CONIC_CONFIGURATIONS.md',ROOT/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_SOURCES_REVIEW_20261107.json'];result=dict(task_id='E_multiple_hyperplane',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),all_six_magnetic_parameters=list(map(str,mag)),source_family='nil22_image_double_B',source_case_index=1,whole_affine_brackets={str(i+1)+','+str(j+1):str(s.factor(P[i,j]))for i,j in itertools.combinations(range(4),2)if P[i,j]!=0},whole_homogeneous_brackets={str(i)+','+str(j):str(s.factor(homogeneous[i,j]))for i,j in PAIRS if homogeneous[i,j]!=0},actual_anticanonical_quintic=str(D),quadratic_factors=[str(QA),str(QB)],quadratic_rank3_certificates=normals,all4_strict_and350_projective_Jacobi_identically_zero=True,source_endpoint_all6magnetic_parameters_exact=True,family_regular_at0=True,only_family_parameter_inverted='epsilon-2',generic_family_locus='epsilon!=0,2',generic_divisor_reduced_rank4=True,generic_accepted_plane_configuration='shared_tangent_generic',configuration_linear_map='(h,x,y)=(X,Y,-2Z-4X)',configuration_original_quadric_scalar='QB=-2*(X*Z-Y^2), QA=-epsilon*((2/epsilon)*X*Z-Y^2+(8/epsilon)*X^2)',configuration_parameters=dict(u='2/epsilon',nonzero_xsquare_coefficient='8/epsilon'),fixed_H_common_tangent_point='[h:x:y]=[0:0:1]',whole_source_component_containment='R4_CPLUS_EXT',containment_proof='For everyfixed6magnetic R, P_e is regularat0, allJacobi andactualD exact. On densee!=0,2 the two rank3conics have the accepted sharedtangentgeneric configuration (u!=0,1,d!=0). The accepted ENTIREsixmagnetic free-plane theorem places P_e in CPLUS_EXT24. Closedness places P0 there, for every R without localization. Alloriginalsourceparameters andtheirboundaries retained.',new_component_added=False,independent_review_required=True,seconds=time.monotonic()-started)
 out=ROOT/'reports/CLASSIFICATION_NIL22_DOUBLE_B_CPLUS_LIMIT_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(whole_double_B_CPLUS_containment=True,seconds=result['seconds'])))
if __name__=='__main__':main()

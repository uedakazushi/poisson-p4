#!/usr/bin/env python3
"""Exact whole-magnetic lower-root sources of the nil22-image repeated-H stratum."""
from pathlib import Path
import hashlib
ROOT=Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0,str(ROOT/'code'))
import sympy as s,itertools,time,json
from sympy.polys.matrices import DomainMatrix
from classification_three_hyperplanes_quadric_independent import ambient_matrices
from classification_rank3_mixed_component_review import minor_certificate
start=time.monotonic();h,q1,q2,x,y=s.symbols('h x1:5');xx=(h,q1,q2,x,y);vs=(q1,q2,x,y);p=s.Matrix([x,y]);mm=s.symbols('m0:6');F0=sum(a*b for a,b in zip(mm,[1,x,y,x*x,x*y,y*y]));a,u,n,c,d,lam=s.symbols('a u n c d lam')
cases=[
 ('double_A',x*x*y,s.Matrix([[s.Rational(8,3),1],[1,-s.Rational(3,64)]]),s.Matrix([-s.Rational(1,8),-s.Rational(5,3)]),s.diag(-s.Rational(1,4),-s.Rational(1,2)),[]),
 ('double_B',x*x*y,s.diag(-8,1),s.Matrix([0,1]),s.diag(-s.Rational(1,4),-s.Rational(1,2)),[]),
 ('double_C',x*x*y,s.Matrix([[2,1],[1,0]]),s.Matrix([0,-2]),s.diag(0,-1),[]),
 ('double_D',x*x*y,s.diag(u,0),s.Matrix([0,-1-(a+s.Rational(1,2))*u]),s.diag(a,-1-2*a),[a,u]),
 ('triple_A',x**3,s.diag(-18,-6),s.Matrix([1,0]),s.diag(-s.Rational(1,3),-s.Rational(1,2)),[]),
 ('triple_B',x**3,s.Matrix([[0,1],[1,0]]),s.zeros(2,1),s.diag(-s.Rational(1,3),-s.Rational(2,3)),[]),
 ('triple_C',x**3,s.zeros(2),s.Matrix([0,n]),s.Matrix([[-s.Rational(1,3),0],[c,d]]),[n,c,d]),
 ('quadratic_A',0,s.Matrix([[0,1],[1,0]]),s.zeros(2,1),s.diag(lam-s.Rational(1,2),-lam-s.Rational(1,2)),[lam]),
 ('quadratic_B',0,s.diag(1,0),s.Matrix([0,n]),s.Matrix([[-s.Rational(1,2),0],[c,d]]),[n,c,d])]
records=[]
encode=lambda Q:{str(i)+','+str(j):str(s.factor(Q[i,j]))for i,j in itertools.combinations(range(5),2)if Q[i,j]!=0}
for name,D3,S,V0,V1,extra in cases:
 D3=s.sympify(D3);D=S.det()+(p.T*S*p)[0]/2+D3;H=S+s.hessian(D3,(x,y))/2;V=V0+V1*p;M=(-H.adjugate()+p*V.T).applyfunc(s.expand);assert s.expand(M.det()-D)==0,(name,s.factor(M.det()-D))
 div=[sum(s.diff(M[i,j],(x,y)[j])for j in range(2))for i in range(2)];F=F0+(y+div[1])*q1+(-x-div[0])*q2;P=s.zeros(4);P[0,1]=F;P[1,0]=-F
 for i in range(2):
  for j in range(2):P[i,2+j]=M[i,j];P[2+j,i]=-M[i,j]
 J=[s.expand(sum(P[i,l]*s.diff(P[j,k],vs[l])+P[j,l]*s.diff(P[k,i],vs[l])+P[k,l]*s.diff(P[i,j],vs[l])for l in range(4)))for i,j,k in itertools.combinations(range(4),3)];assert J==[0]*4,(name,J)
 N=s.Matrix([sum(s.diff(P[i,j],vs[i])for i in range(4))for j in range(4)]).applyfunc(s.expand);assert N==s.Matrix([x,y,0,0]),(name,N)
 Q=s.zeros(5)
 for i,j in itertools.combinations(range(4),2):Q[i+1,j+1]=s.expand(h*h*P[i,j].subs({v:v/h for v in vs},simultaneous=True));Q[j+1,i+1]=-Q[i+1,j+1]
 params=list(mm)+extra;scale=s.ilcm(1,*[s.denom(z)for z in Q for z in s.Poly(z,*xx,*params).coeffs()]);values=dict(zip(mm,(1,2,3,5,7,11)));values.update({a:3,u:2,n:2,c:3,d:5,lam:7});point=scale*Q.subs(values);dirs=[scale*Q.diff(z).subs(values)for z in params];DD,VV,GG,FF=ambient_matrices(point,xx,dirs);rank=lambda Z:len(DomainMatrix.from_Matrix(Z).convert_to(s.QQ).rref()[1]);ranks=list(map(rank,[DD,GG,FF]));assert DD*FF==s.zeros(350,FF.cols)
 certs={key:minor_certificate(Z,rank(Z))for key,Z in [('D',DD),('vertical',VV),('orbit',GG),('source',FF)]}
 records.append(dict(name='nil22_image_'+name,free_parameters=list(map(str,params)),normal_D3=str(D3),S=[[str(a)for a in row]for row in S.tolist()],V0=list(map(str,V0)),V1=[[str(a)for a in row]for row in V1.tolist()],M=[[str(a)for a in row]for row in M.tolist()],affine_D=str(s.factor(D)),
  full_six_magnetic_source=encode(Q),universal_brackets=encode(scale*Q),witness_brackets=encode(point),source_integer_scale=str(scale),full_source_specification=dict(free_parameters=list(map(str,params)),witness={str(z):str(values[z])for z in params}),
  all_universal_Jacobi_modular_and_Pf_checked=True,all_D_times_full_actual_source_zero=True,highest_Pf_zero=True,ranks_QQ=dict(D=ranks[0],orbit=ranks[1],source=ranks[2]),certificates=certs,projective_source_dimension_lower=ranks[2]-25,full_tangent_upper=125-ranks[0],true_normal_count=150-ranks[0]-ranks[2],componenthood_or_membership_proved=False))
 print(name,'ranks',ranks,'D',s.factor(D),'params',len(params),'scale',scale,'seconds',time.monotonic()-start,flush=True)
here=Path(__file__);inputs=[here,ROOT/'research/classification/MODULAR_NIL_RANK2_FIRSTJET_BOUND.md',ROOT/'reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json',ROOT/'code/classification_three_hyperplanes_quadric_independent.py',ROOT/'code/classification_rank3_mixed_component_review.py']
result=dict(task_id='E_multiple_hyperplane',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=here.read_text(),cases=records,scope='Nine exact whole6-magnetic sources for lower binary-root strata; exhaustion requires the separate written coefficient branching and GL2 normalization. No componenthood follows from these rank bounds alone.',seconds=time.monotonic()-start)
out=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json'
if out.exists():raise FileExistsError(str(out))
out.write_text(json.dumps(result,indent=2)+'\n')

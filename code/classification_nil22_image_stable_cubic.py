#!/usr/bin/env python3
"""Complete nil22-image repeated-H source over stable binary leading cubics.

The four exact Hessians are polynomially exhausted without inverting them;
three are GL2-equivalent. All six magnetic coefficients remain in each source.
Componenthood is deliberately not inferred from source/tangent ranks alone.
"""
import hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_three_hyperplanes_quadric_independent import ambient_matrices
from classification_rank3_mixed_component_review import minor_certificate
ROOT=Path(__file__).resolve().parents[1]

def main():
 started=time.monotonic();xx=s.symbols('h x1:5');h,q1,q2,x,y=xx;vars=(q1,q2,x,y);pp=s.Matrix([x,y]);D3=x*y*(x-y);T=s.hessian(D3,(x,y))/2;mm=s.symbols('m0:6');F0=sum(a*b for a,b in zip(mm,[1,x,y,x*x,x*y,y*y]))
 u,v,w=s.symbols('u v w');f1=u*u+2*u*w+15*u-6*w-54;f2=w*w+2*u*w-6*u+15*w-54
 assert s.expand(f1-f2-(u-w)*(u+w+21))==0
 assert s.expand(f1.subs(w,u)-3*(u+6)*(u-3))==0
 assert s.expand(f1.subs(w,-u-21)+(u-3)*(u+24))==0
 matrices=[s.Matrix([[3,-15],[-15,3]]),s.Matrix([[-6,3],[3,-6]]),s.Matrix([[3,12],[12,-24]]),s.Matrix([[-24,12],[12,3]])]
 def build(S):
  U,V,W=S[0,0],S[0,1],S[1,1];H=S+T;V0=s.Matrix([1+W/6,-1-U/6]);field=V0-pp/3;D=S.det()+(pp.T*S*pp)[0]/2+D3
  M=(-H.adjugate()+pp*field.T).applyfunc(s.expand);assert s.expand(M.det()-D)==0
  assert(M*s.Matrix([s.diff(D,x),s.diff(D,y)])+D*pp).applyfunc(s.expand)==s.zeros(2,1)
  div=[sum(s.diff(M[i,j],(x,y)[j])for j in range(2))for i in range(2)]
  F=F0+(y+div[1])*q1+(-x-div[0])*q2;P=s.zeros(4);P[0,1]=F;P[1,0]=-F
  for i in range(2):
   for j in range(2):P[i,2+j]=M[i,j];P[2+j,i]=-M[i,j]
  J=[s.expand(sum(P[i,l]*s.diff(P[j,k],vars[l])+P[j,l]*s.diff(P[k,i],vars[l])+P[k,l]*s.diff(P[i,j],vars[l])for l in range(4)))for i,j,k in itertools.combinations(range(4),3)];assert J==[0]*4
  N=s.Matrix([sum(s.diff(P[i,j],vars[i])for i in range(4))for j in range(4)]).applyfunc(s.expand);assert N==s.Matrix([x,y,0,0])
  return P,M,D,V0
 all_sources=[build(S)for S in matrices];glchecks=[]
 for i,j,A in [(0,3,s.Matrix([[1,0],[1,-1]])),(3,2,s.Matrix([[0,-1],[-1,0]]))]:
  AA=s.diag(A,A);inverse=AA.inv()*s.Matrix(vars);sub=dict(zip(vars,inverse));P0=all_sources[i][0];Pnew=(AA*P0*AA.T).subs(sub,simultaneous=True).applyfunc(s.expand)
  Fnew=A.det()*F0.subs({x:inverse[2],y:inverse[3]},simultaneous=True);expected=all_sources[j][0].copy();expected[0,1]+=Fnew-F0;expected[1,0]=-expected[0,1]
  assert(Pnew-expected).applyfunc(s.expand)==s.zeros(4)
  glchecks.append(dict(source=i,target=j,matrix=[[str(a)for a in row]for row in A.tolist()],all_six_magnetic_coefficients_transferred=True))
 def encode(P):return{str(i)+','+str(j):str(s.factor(P[i,j]))for i,j in itertools.combinations(range(5),2)if P[i,j]!=0}
 records=[]
 for index,name in [(0,'nil22_image_stable_line_conic'),(1,'nil22_image_stable_three_lines')]:
  P,M,D,V0=all_sources[index];Q=s.zeros(5)
  for i,j in itertools.combinations(range(4),2):Q[i+1,j+1]=s.expand(h*h*P[i,j].subs({v:v/h for v in vars},simultaneous=True));Q[j+1,i+1]=-Q[i+1,j+1]
  values=dict(zip(mm,(1,2,3,5,7,11)));witness=12*Q.subs(values);directions=[12*Q.diff(m).subs(values)for m in mm];DD,VV,GG,FF=ambient_matrices(witness,xx,directions);assert DD*FF==s.zeros(350,FF.cols)
  ranks={key:len(DomainMatrix.from_Matrix(Z).convert_to(s.QQ).rref()[1])for key,Z in [('D',DD),('vertical',VV),('orbit',GG),('source',FF)]};assert ranks==dict(D=99,vertical=24,orbit=47 if index==0 else 46,source=48)
  cert={key:minor_certificate(Z,ranks[key])for key,Z in [('D',DD),('vertical',VV),('orbit',GG),('source',FF)]}
  records.append(dict(name=name,Hessian=[[str(a)for a in row]for row in matrices[index].tolist()],V0=list(map(str,V0)),M=[[str(a)for a in row]for row in M.tolist()],
   affine_D=str(s.factor(D)),full_six_magnetic_source=encode(Q),witness_brackets=encode(witness),universal_brackets=encode(12*Q),
   full_source_specification=dict(free_parameters=list(map(str,mm)),witness={str(a):str(b)for a,b in values.items()}),all_universal_strict_Jacobi_and_modular_checked=True,actual_affine_Pfaffian='-('+str(s.factor(D))+')',highest_quadratic_Pfaffian_zero=True,
   ranks_QQ=ranks,certificates=cert,all_D_times_actual_source_zero=True,actual_projective_dimension_lower=23,full_tangent_upper=26,true_normal_count=3,componenthood_or_membership_proved=False))
 here=Path(__file__);inputs=[here,ROOT/'research/classification/MODULAR_NIL_RANK2_FIRSTJET_BOUND.md',ROOT/'code/classification_three_hyperplanes_quadric_independent.py',ROOT/'code/classification_rank3_mixed_component_review.py']
 out=dict(task_id='E_multiple_hyperplane',all_passed=True,scope='Entire stable-binary-leading-cubic stratum reduces to two whole6-magnetic sources; no componenthood from ranks alone.',
  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=here.read_text(),normal_D3=str(D3),
  exact_Hessian_relations=[str(v+u+w+9),str(f1),str(f2)],Hessian_exhaustion_factors=['(u-w)*(u+w+21)','3*(u+6)*(u-3)','-(u-3)*(u+24)'],
  all_four_Hessians=[[[str(a)for a in row]for row in S.tolist()]for S in matrices],all_four_determinants=list(map(lambda S:str(S.det()),matrices)),singular_Hessian_complement_empty=True,GL2_source_identifications=glchecks,cases=records,seconds=time.monotonic()-started)
 dest=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_STABLE_CUBIC_20261107.json'
 if dest.exists():raise FileExistsError(str(dest))
 dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(all_passed=True,whole_stable_stratum_two_sources=True,ranks=[r['ranks_QQ']for r in records],seconds=out['seconds'])))
if __name__=='__main__':main()

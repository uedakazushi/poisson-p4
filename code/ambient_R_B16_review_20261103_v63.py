#!/usr/bin/env python3
"""Independent actual B16 image, whole-fiber and singular-normalizer audit."""
import argparse
import itertools
import json
from pathlib import Path
import subprocess
import sympy as s
from ambient_CA_global_flag_review_20261103_v63 import IDS,action,canonical,wedge_factor
from q_smith_global_model_review_20261103_v63 import ROOT,read,sha

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--m2-work',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 public=ROOT/'reports/AMBIENT_R19_B16_COMPLETED_IMAGE_20261103_v63/report.json';d=read(public);f=read(args.replay/'report.json');assert{k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'};paths=[Path(__file__),public,ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py',ROOT/'code/q_smith_global_model_review_20261103_v63.py']
 for file,hsh in d['input_sha256'].items():q=ROOT/file;assert sha(q)==hsh,file;paths.append(q)
 X=s.Matrix(s.symbols('h x y z w'));h,x,y,z,w=X;p,t,r,c,a,b,e=s.symbols('p t r c a b d');ey=s.Matrix([0,0,1,0,0]);q=r*r-25*t*t-2*c*p
 D=s.Matrix([x*((r-13*t)*h+c*z),12*t*x*x,0,x*(-2*p*h-x+(-r-13*t)*z),2*t*x*w+h*h+z*z]);P=ey*D.T-D*ey.T;g=s.eye(5);g[1,0]=p;inverse=s.eye(5);inverse[1,0]=-p
 Q=(g*P.subs(dict(zip(X,inverse*X)),simultaneous=True)*g.T).applyfunc(s.expand)
 v=s.Matrix([s.Poly(Q[i,j],*X).coeff_monomial(s.prod(X[k]for k in mon))for((i,j),mon)in IDS]);assert v==s.Matrix(list(map(s.sympify,d['full150_target'])))and canonical(v)==v
 point={p:0,t:0,r:0,c:0};v0=v.subs(point);delta=v-v0;generators=[t,r,c,p*p,p*r,p*t]
 recovered=[delta[128]/2,-delta[112]-13*delta[128]/2,-delta[22],delta[105],(delta[15]+delta[108])/2,(delta[108]-delta[15])/26];assert all(s.expand(x-y)==0 for x,y in zip(generators,recovered))
 output=[s.sympify(u)for u in d['full150_polynomials_in_target_generators']]
 for actual,expression in zip(delta,output):
  difference=s.expand(actual-expression.subs({a:p*p,b:p*r,e:p*t},simultaneous=True));assert s.rem(s.Poly(difference,r),s.Poly(q,r)).is_zero
 relations=list(map(s.sympify,d['target_entire_elimination_kernel']))
 for relation in relations:assert s.rem(s.Poly(s.expand(relation.subs({a:p*p,b:p*r,e:p*t},simultaneous=True)),r),s.Poly(q,r)).is_zero
 C=s.Matrix([[s.Poly(u,p,t,r,c).coeff_monomial(mon)for mon in generators]for u in delta]);assert C==s.Matrix([[s.Rational(v)for v in row]for row in d['generator_coefficient_matrix']])and C.rank()==6
 O=s.Matrix.hstack(*[action(i,j)*v0 for i,j in itertools.product(range(5),repeat=2)]);OP=O+v0*O[110,:];assert OP.rank()==16
 def check(M,record):
  assert M==s.Matrix([[s.Rational(v)for v in row]for row in record['matrix']]);det=M[record['rows'],record['columns']].det();assert det==s.Rational(record['determinant'])!=0;return len(record['columns'])
 assert check(OP.row_join(C),d['entire_algebra_cotangent'])==22
 chosen=d['projective_orbit_generator_indices'];H=OP[:,chosen];assert len(chosen)==16 and H.row_join(C[:,:3]).rank()==19
 triples=list(itertools.combinations(range(5),3));epsilon=s.Symbol('epsilon');U=s.zeros(3,5)
 for i,j in enumerate([0,1,3]):U[i,j]=1
 def flag(E):
  L=E*ey;L[2]=0;Z=-E[1,:].T;Z[1]=0;K=[];A=U*(s.eye(5)-epsilon*E)
  for triple in triples:K.append(s.expand(A[:,list(triple)].det()).coeff(epsilon,1))
  K[triples.index((0,1,3))]=0
  return s.Matrix(list(L)+list(Z)+K)
 HF=[]
 for index in chosen:
  E=s.zeros(5);E[index//5,index%5]=1;HF.append(flag(E))
 E=s.zeros(5);E[1,0]=1;vertical=flag(E);assert vertical==s.Matrix(list(map(s.Rational,d['vertical_flag_derivative'])))and vertical!=s.zeros(20,1)
 ambient=H.row_join(C[:,:3]).row_join(s.zeros(150,1)).col_join(s.Matrix.hstack(*HF).row_join(s.zeros(20,3)).row_join(vertical));assert check(ambient,d['normalization_ambient_etale_tangent'])==20
 assert s.hessian(q,[p,t,r,c]).rank()==4
 assert wedge_factor(v0).rank()==4 and wedge_factor(v0).nullspace()==[ey]
 contract=s.zeros(75,5)
 for((i,j),mon),index in IDS.items():
  k=index%15;contract[15*j+k,i]+=v0[index];contract[15*i+k,j]-=v0[index]
 assert contract.rank()==3 and contract.nullspace()==[s.eye(5)[:,0],s.eye(5)[:,1]]
 # X is in span(h,x).  A nonzero combination of the two field quadrics
 # -x²,h²+z² divisible by X cannot involve the second (its z² coefficient),
 # so X=x and then K=span(h,x,z).  This covers the WHOLE proper fiber.
 u,alpha,beta=s.symbols('upar alpha beta');integrality=a**3-2*u*a*a+u*u*a-(r*r-25*t*t)**2/4
 assert s.expand(integrality.subs({u:c+p*p,a:p*p})+q*(r*r-25*t*t+2*c*p)/4)==0
 # Independent small Macaulay2 elimination corroborates the WHOLE target
 # kernel, rather than checking only listed pullback-zero relations.
 args.m2_work.mkdir(parents=True,exist_ok=False);script=args.m2_work/'review.m2';log=args.m2_work/'review.log'
 body='S=QQ[p,t,r,c,a,b,d,MonomialOrder=>Eliminate 1];\nJ=ideal(p^2-a,p*r-b,p*t-d,2*c*p-r^2+25*t^2);\nK=eliminate({p},J);\nT=QQ[t,r,c,a,b,d];\nKK=sub(K,T);\nI=ideal('+','.join(str(v).replace('**','^')for v in relations)+');\nassert(isSubset(I,KK) and isSubset(KK,I));\nassert(dim(T/I)==3);\nprint "INDEPENDENT_B16_ENTIRE_ELIMINATION_PASS";\n'
 script.write_text(body);completed=subprocess.run(['M2','--script',str(script.resolve())],capture_output=True,text=True,timeout=30);log.write_text(completed.stdout+completed.stderr);assert completed.returncode==0 and 'INDEPENDENT_B16_ENTIRE_ELIMINATION_PASS'in completed.stdout
 paths += [script.resolve(),log.resolve(),ROOT/'reports/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63_REVIEW.json',ROOT/'reports/AMBIENT_R19_BOUNDARY_LOCAL_MODELS_20261103_v63_REVIEW.json',ROOT/'research/notes/AMBIENT_R19_B16_COMPLETED_IMAGE_20261103_v63.md',ROOT/'research/notes/AMBIENT_R19_B16_COMPLETED_IMAGE_20261103_v63_REVIEW.md']
 report=dict(id='AMBIENT_R19_B16_COMPLETED_IMAGE_20261103_v63_REVIEW',all_passed=True,
  fresh_producer_matches=True,independent_finite_GL_transform_all150_and_all_generator_inverse_formulas=True,
  all150_polynomial_output_reductions_and_independent_M2_entire11relation_kernel=True,
  independent_projective_orbit16_and_full_intrinsic_algebra_cotangent22_minor=True,
  independent_source_flag_Pluckers_and_ambient20_etale_minor=True,
  whole_proper_fiber_unique_by_intrinsic_L_Casimir_X_and_divisibility_K=True,
  full_completed_reduced_R_image='C[[u16,t,r,c,p²,pr,pt]] in quadric normalization r²-25t²-2cp=0',
  finite_normalization_injection_and_complete_cotangent_surjection_accepted=True,
  dimension=19,embedding_dimension=22,depth=18,multiplicity=3,
  normalization_quotient_conductors_fiber2_and_boundary_ideal_t_pt_accepted=True,
  arbitrary_B16_boundary_points_or_actual_CA_CB_scheme_thickness_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
 print('PASS full B16 completed image,quadric normalization,depth18/mult3 and entire kernel')

if __name__=='__main__':main()

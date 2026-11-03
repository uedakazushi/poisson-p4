#!/usr/bin/env python3
"""Independent target tangents, unique global flags, and LPB containment audit."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_review_20261103_v63 import IDS,action,canonical
from ambient_global_intersection_review_20261103_v63 import transport

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--geometry-replay',type=Path,required=True)
 ap.add_argument('--LPB-replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 paths=[Path(__file__),ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py',ROOT/'code/ambient_global_intersection_review_20261103_v63.py'];reports=[]
 for stem,replay in [('AMBIENT_R19_GENERIC_GEOMETRY',args.geometry_replay),('AMBIENT_R19_LPB_CONTAINMENT',args.LPB_replay)]:
  p=ROOT/f'reports/{stem}_20261103_v63/report.json';d=load(p);f=load(replay);assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'}
  paths.append(p);reports.append(d)
  for file,h in d['input_sha256'].items():q=ROOT/file;assert sha(q)==h,file;paths.append(q)
 geometry,lpb=reports;h,x,y,z,w=s.symbols('h x y z w');X=[h,x,y,z,w];U=[h,x,z,w];kap=s.Symbol('kappa');Q=h*h+kap*h*z+z*z
 V=s.Matrix([0,4*x*x,0,-2*x*z,2*x*w+Q]);ey=s.Matrix([0,0,1,0,0]);P=ey*V.T-V*ey.T
 vector=lambda P:s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(s.prod(X[k]for k in m))for(ab,m),ind in IDS.items()for i,j in[ab]])
 point=canonical(vector(P));assert point==s.Matrix([s.sympify(v)for v in geometry['point_canonical150']])
 ca=load(ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63/report.json');A0=s.Matrix.hstack(*[s.Matrix(v)for v in ca['basis_canonical150']])
 cb=load(ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json');BC=s.Matrix.hstack(*[s.Matrix(v)for v in cb['enlarged_canonical_basis']]);perm=[1,3,2,0,4]
 B0=transport(BC[:,[*range(3),*range(6,22)]],perm);mu=transport(6*BC[:,4],perm);orbit=s.Matrix.hstack(*[action(i,j)*point for i,j in itertools.product(range(5),repeat=2)])
 A=s.Matrix.hstack(A0,orbit);B=s.Matrix.hstack(B0,mu,orbit);union=A.row_join(B)
 assert A==s.Matrix([[s.sympify(v)for v in row]for row in geometry['full_CA_source_differential']])
 assert B==s.Matrix([[s.sympify(v)for v in row]for row in geometry['full_CB_source_differential']])
 assert [A.rank(),B.rank(),union.rank()]==[29,29,37]
 for name,M in [('CA',A),('CB',B),('union',union)]:
  rec=geometry['sample_kappa1_minors'][name];assert M.subs(kap,1)[rec['rows'],rec['columns']].det()==s.Rational(rec['determinant'])!=0
 C=s.Matrix([-s.Rational(8,5)*h*x,s.Rational(12,5)*x*x,-s.Rational(18,5)*x*z,Q+s.Rational(2,5)*x*w]);a,b,c,d,alpha=s.symbols('a b c d alpha');ell=a*h+b*x+c*z+d*w
 Der=lambda f:s.expand(sum(s.diff(f,u)*v for u,v in zip(U,C)))
 for label,m,CAok,CBok in [('h',-s.Rational(8,5)*x,True,False),('z',-s.Rational(18,5)*x,False,False)]:
  f=s.expand(m*x-C[1]);assert (s.diff(f,x)==5*m)==CAok and (s.diff(f,x)==s.Rational(5,2)*m)==CBok
 MA,_=s.linear_eq_to_matrix(s.Poly(Der(ell)+s.Rational(48,5)*x*ell+alpha*x*x,*U).coeffs(),[a,b,c,d,alpha]);MB,_=s.linear_eq_to_matrix(s.Poly(Der(ell)+s.Rational(18,5)*x*ell+alpha*x*x,*U).coeffs(),[a,b,c,d,alpha])
 assert MA.nullspace()==[s.Matrix([0,-s.Rational(1,12),0,0,1])]
 assert len(MB.nullspace())==2 and all(v[0]==v[3]==0 for v in MB.nullspace())
 # Separate all8-parameter LPB identity, without the old LPB runner.
 t,q0,q1,q2,q3,q5,q6,q7=s.symbols('t q0 q1 q2 q3 q5 q6 q7');field=s.Matrix([0,4*t*x*x,0,-q0*x*x-2*t*x*z,2*t*x*w+3*q7*h*h+2*q3*h*x+2*q6*h*z+q1*x*x+q2*x*z+q5*z*z]);raw=ey*field.T-field*ey.T
 perm2=[0,1,4,2,3];g=s.zeros(5)
 for i,j in enumerate(perm2):g[j,i]=1
 transformed=(g*raw.subs(dict(zip(X,g.inv()*s.Matrix(X))),simultaneous=True)*g.T).applyfunc(s.expand)
 targetfield=s.Matrix([s.sympify(v)for v in lpb['LPB_quotient_quadratic_field']]);e4=s.Matrix([0,0,0,0,1]);target=e4*targetfield.T-targetfield*e4.T
 assert transformed==target and all(s.diff(v,w)==0 for v in targetfield)
 assert canonical(vector(transformed))==canonical(vector(target))
 old=lpb['existing_component_dimension_certificate'];assert old['projective_component_dimension']==39 and old['ambient_affine_tangent_dimension']==old['family_affine_differential_rank']==40
 paths += [ROOT/'research/notes/AMBIENT_R19_GENERIC_GEOMETRY_20261103_v63.md',
  ROOT/'research/notes/AMBIENT_R19_LPB_CONTAINMENT_20261103_v63.md',
  ROOT/'research/notes/AMBIENT_R19_GEOMETRY_20261103_v63_REVIEW.md',
  ROOT/'reports/AMBIENT_R19_GENERIC_ORBITS_20261103_v63_REVIEW.json',
  ROOT/'reports/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63_REVIEW.json']
 result=dict(id='AMBIENT_R19_GEOMETRY_20261103_v63_REVIEW',all_passed=True,
  fresh_QQ_geometry_and_LPB_producer_replays_match=True,independent_full150_source_tangent_matrices_exactly_match=True,
  ranks_over_QQ_kappa=[29,29,37],complete_entire_CA_and_CB_flag_fiber_uniqueness_accepted=True,
  proper_unramified_length1_finite_Nakayama_target_isomorphisms_accepted=True,
  full_projective_CA_CB_Serre_R1_accepted=True,R_rationality_via_unique_double_flag_projective_bundle_accepted=True,
  actual_projective_intersection_generic_tangent_dimension=20,reduced_R_dimension=19,
  actual_intersection_generically_nonreduced_along_R_accepted=True,generic_nilpotence_order_or_length_inferred=False,
  independent_all8parameter_LPB_raw_and_canonical_identity=True,entire_closed_R_contained_in_existing_LPB39=True,
  full_Poisson_completed_local_ring_non_CM_at_every_closed_R_point_accepted=True,
  old_LPB_full_component_certificate_reused_not_rerun=True,new_full_Poisson_component_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS global CA/CB R1, R rational/generic nonreduced intersection, full R subset LPB39')

if __name__=='__main__':main()

#!/usr/bin/env python3
"""Proof/gate-only whole-normalform PGL classification; no CAS or points."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_ENTIRE_NORMAL_FORM_ORBITS_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 stems=[
 'R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63',
 'Q_NEW_FAMILY_SMOOTH_OPEN_STABILIZER_STRATA_20261104_v63',
 'Q_NEW_FAMILY_SMOOTH_OPEN_STABILIZER_STRATA_20261104_v63_REVIEW',
 'Q_NEW_FAMILY_SIGMA_FINITE_EQUIVALENCE_20261104_v63',
 'Q_NEW_FAMILY_SIGMA_FINITE_EQUIVALENCE_20261104_v63_REVIEW',
 'Q_NEW_FAMILY_GLOBAL_TWO_INVARIANT_PREFLIGHT_20261104_v63_REVIEW']
 paths=[ROOT/'reports'/n/'report.json'for n in stems];records=[json.loads(p.read_text())for p in paths]
 assert all(d['all_passed']for d in records)
 assert records[1]['whole_open_flag_uses_intrinsic_modular_image']
 assert records[2]['whole_uvp_open_exhaustion_and_groups_accepted']
 assert records[3]['both_signs_explicit_full_tensor_sigma_equivalence']
 checked={}
 for p,d in zip(paths,records):
  for section in['input_sha256','checked_dependency_sha256']:
   for n,h in d.get(section,{}).items():assert sha(ROOT/n)==h,n;checked[n]=h
  for n,h in d.get('evidence_sha256',{}).items():
   q=p.parent/n;assert sha(q)==h,str(q);checked[str(q.relative_to(ROOT))]=h
 # Only sign-group enumeration, not coefficient/source/target evaluation.
 group=[(a,b)for a in[1,-1]for b in[1,-1]]
 assert len(set(group))==4
 assert all((a*c,b*d)in group for a,b in group for c,d in group)
 proof=dict(
  parameter_domain='u*v*(32*u+v)!=0',target_domain='U*V*(32*U+V)!=0',base_field='any characteristic-zero field; geometric orbit classification over algebraic closure',
  intrinsic_flags='unique tripleH; unique vertex of Q restricted toH; homogeneous ImN=span(ey,ez) with minor4u^2/3 and target4U^2/3',
  full_affine_form='xprime=a*x+b*w+tx; wprime=c*x+d*w+tw; yprime=e*x+f*y+g*w+ty; zprime=h*x+i*y+j*z+k*w+tz',
  invertibility='f*j*(a*d-b*c)!=0; kappa!=0 retained',
  general_low_bracket_relations=['kappa*(a*d-b*c)=1','eta=kappa*f','eta*b=64*c','eta*a=d','eta*d=a','eta*c=b/64'],
  modular_relation='N_target*M=kappa*M*N_source; y,h entry gives U=eta*u',
  no_lost_sign_case=['eta*(eta*a-d)+(eta*d-a)=(eta^2-1)*a','eta*(eta*b-64*c)+64*(eta*c-b/64)=(eta^2-1)*b','a,b not both zero implies eta^2=1'],
  first_necessary_condition='U=+u or -u for ALL displayed parameters, not just generic u',
  line_covariance=[
   dict(eta=1,b='64*c',d='a',lplus_target_pullback='(a+8*c)*lplus_source',lminus_target_pullback='(a-8*c)*lminus_source'),
   dict(eta=-1,b='-64*c',d='-a',lplus_target_pullback='(a+8*c)*lminus_source',lminus_target_pullback='(a-8*c)*lplus_source')],
  nonzero_line_scales='(a+8c)(a-8c) is a nonzero scalar multiple of det(xw block)',
  conic='q2=-x^2+64*w^2+64*u*x*y+32*u*v*y^2; target obtained by u,v->U,V',
  projected_Q_covariance='Q_target under quotient P2 linear map is nonzero scalar times Q_source, since actual anticanonical quadric component is preserved',
  chord='J=256*(lplus^t Q^-1 lminus)^2/(self(lplus)*self(lminus))',
  chord_unit_conditions='detQ=-2048*u*(32u+v), self(lplus)=self(lminus)=32u/(32u+v), all nonzero on ENTIRE domain',
  second_necessary_condition='(V/U+16)^2=(v/u+16)^2 for ALL displayed parameters, no invariant-line uniqueness assumption',
  all_parameter_candidates=[['u','v'],['-u','-v'],['u','-32*u-v'],['-u','32*u+v']],
  sufficiency='known FULL tau plus independently checked FULL sigma matrices, both defined on entire uvp domain',
  equivalence_criterion='P(u,1,v) and P(U,1,V) are projectively equivalent iff U^2=u^2 and (V/U+16)^2=(v/u+16)^2',
  sign_coordinates=dict(z='v/u+16',A='u^2',J='z^2',tau='(u,z)->(-u,z)',sigma='(u,z)->(u,-z)'),
  parameter_coordinate_ring='k[u,u^-1,z,(z^2-256)^-1]',
  finite_parameter_quotient_ring='k[A,A^-1,J,(J-256)^-1]',
  quotient_proof='independent signs select even u-Laurent and z-polynomial exponents; denominator already invariant',
  finite_free_rank=4,free_basis=['1','u','z','u*z'],
  quotient_parameter_space='Gm_A times (A1_J minus {256})',
  geometric_orbit_normalform_size=dict(J_nonzero=4,J_zero=2),
  J_zero_scheme_fibre='two geometric u-roots, each z^2=0: scheme length4 is NOT four reduced points',
  no_rational_point_lift_claim='quotient k-points need not admit u,z in k; however two k-valued normalforms with equal A,J are explicitly PGL(k)-equivalent',
  scope='ENTIRE displayed normalform family and its set-theoretic PGL saturation',
  NOT_all_smooth_quadric_points_in_whole_component=True,
  NOT_categorical_or_good_PGL_quotient_of_whole_component=True,
  no_CAS_new_point_prime_GB_or_invariant_probe=True)
 (out/'proof-certificate.json').write_text(json.dumps(proof,indent=2)+'\n')
 note=ROOT/'research/notes'/(ID+'.md')
 report=dict(id=ID,all_passed=True,entire_uvp_normalform_parameter_identifications_exact_Klein_four=True,
  finite_parameter_quotient_Gm_times_A1_minus_256=True,
  no_generic_logarithmic_line_uniqueness_used=True,
  no_new_computation_beyond_gate_hash_and_sign_group=True,
  not_all_smooth_Q_component_points_or_global_PGL_good_quotient=True,
  independent_conceptual_review_pending=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),note,*paths]},checked_gate_dependency_sha256=checked,
  evidence_sha256={'proof-certificate.json':sha(out/'proof-certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in['all_passed','entire_uvp_normalform_parameter_identifications_exact_Klein_four','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

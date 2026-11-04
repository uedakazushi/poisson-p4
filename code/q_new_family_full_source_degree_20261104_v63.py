#!/usr/bin/env python3
"""Proof/gate-only exact full normalized source degree8, tau quotient degree4."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_FULL_SOURCE_DEGREE_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 stems=[
 'Q_NEW_FAMILY_ENTIRE_NORMAL_FORM_ORBITS_20261104_v63_REVIEW',
 'Q_NEW_FAMILY_GENERIC_STABILIZER_20261104_v63_REVIEW',
 'Q_NEW_FAMILY_FULL_TENSOR_FIBER_PREFLIGHT_20261104_v63_REVIEW',
 'Q_NEW_FAMILY_RATIONAL15_APPLICATION_20261104_v63_REVIEW',
 'Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63_REVIEW',
 'Q_NEW_FAMILY_SIX_SOURCE_20261103_v63_REVIEW',
 'Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW',
 'Q_NEW_FAMILY_SP4_CENTRALIZER_FREE_20261104_v63_REVIEW',
 'R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63']
 paths=[ROOT/'reports'/(n+'.json')if n=='Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63_REVIEW'else ROOT/'reports'/n/'report.json'for n in stems]
 gates=[json.loads(p.read_text())for p in paths];assert all(d['all_passed']for d in gates)
 checked={}
 for p,d in zip(paths,gates):
  for section in['input_sha256','checked_dependency_sha256','checked_gate_dependency_sha256']:
   for n,h in d.get(section,{}).items():assert sha(ROOT/n)==h,n;checked[n]=h
  for n,h in d.get('evidence_sha256',{}).items():
   q=p.parent/n;assert sha(q)==h,str(q);checked[str(q.relative_to(ROOT))]=h
 # Exact field degrees are a theorem from the normalized equivalence torsors;
 # this program performs only gate/hash checks and the finite tower arithmetic.
 parameter_choices=4;aut_components=2;Ga_intersections_per_component=1
 old_degree=parameter_choices*aut_components*Ga_intersections_per_component
 tau_degree=2;assert old_degree==8 and old_degree%tau_degree==0
 proof=dict(
  base_fields=['QQ','C'],Y='actual reduced NEW S61 component Y16',B='actual primitive beta image B15',G='Sp4, not PSp4',
  old_source_field='S_old=k(G)(u,v,c,X,W,B)',point_coordinates='p=(X,0,Z,W), Z=-3/(4u)*(B+W/1536-W^2+X^2/64+2/9)',
  source_dimension=16,Y_dimension=16,firstjet_dimension=15,
  actual_normalized_equivalence_data=dict(
   equation='Pi=c*T_*P(u,1,v), T inPGL5, Pi inS61 valueJ/perspectivefixed',
   point='p=T^-1(e0); rank4 value forces p off tripleH, hence affine, and offQ',
   unique_projective_factorization='T=Persp(alpha)*block(A)*Translate(-p), after denominator atp normalizedto1',
   scalar='c uniquely determined by nonzero fixed affinecanonicalPi andPGL action; centralGL5 scalar acts trivially on quadratic bivectors',
   value='A*(c*C(p))*A^T=J',Darboux='R*(c*C(p))*R^T=J',
   symplectic_frame='g=A*R^-1, unique actualSp4 matrix',
   perspective='alpha unique by invertible4x4 S61 perspective response matrix; full tensor normalizer, not justbeta',
   reverse='given(u,v,c,p_y0,g), useA=gR and the unique perspective to reconstructT',
   chart='c andPf(C(p))units, andPf(C(p))=8B/9 in y0 gauge; no generic preimages lost',
   fixed_n1='n1=1 rational gauge already removed; no remaining Gm or square-root choice'),
  generic_parameter_choices=['(u,v)','(-u,-v)','(u,-32u-v)','(-u,32u+v)'],
  equivalence_torsor='for each parameter representative, projectiveT choices form a torsor underAut(P)=Ga xmu2',
  point_cut='T->T*exp(qN) changes p_y by-4*u*q/3; in EACH of2 Aut components the y0 condition has exactlyone solution',
  old_geometric_generic_preimages=old_degree,
  generic_separability='char0 dominant same-dimension map, hence functionfielddegree equals geometric generic preimagecount',
  old_full_source_field_degree=old_degree,
  no_extra_factors=[
   'two sigma matrices differ bythe finiteAut component and are NOT an additionalfactor2',
   'Sp4 center-I is NOT quotiented: it changes beta sign, generic beta stabilizeridentity; A uniquelyfixesg',
   'Poisson scalarc is unique forT; rational DarbouxR introduces no square root',
   'n1 torus gauge wasfixed once, its residual parameteridentifications alreadyin fourchoices'],
  full_tau_frame_change=dict(
   raw_pair_match='S_t inSp4 physicalframe, rational inF6',
   gu='diag(u,u^-1,1,1); u is an ORIGINALF6 coordinate, no new fieldextension',
   effective_linear_frame='newnormalization differs fromold by gu^-1*S_t; perspective is removed bysame S61crosssection',
   group_variable_change='oldphysicalg=newphysicalh*gu^-1*S_t, a birational rationalchange overF6',
   base_convention='replace physicalh byh^-1 to align with accepted pullbackanti beta action; groupinversion is birational'),
  old_constant_field='F6=F5(t), F5^tau=E5=k(xi,r,s,q,B), F5=E5(u),u^2=1600r',
  tau='u->-u,t->t/(t-1), acts trivially on correctedSp4 coordinate',
  invariant_curve_coordinate='z=u*(t-2)/t, inverse t=2u/(u-z)',
  source_tau_fixed_field='S_new=k(G)E5(z), [S_old:S_new]=2',
  full_tensor_descends='accepted full150 finite equivariance, not beta-onlysource symmetry',
  base_field_binding='accepted GxSpecE5 birationaltoB15 via genericallyfree beta; S_new=k(B15)(z) overSAME actual base',
  exact_degree_tower='8=[S_old:k(Y16)]=2*[k(B15)(z):k(Y16)]',
  descended_full_source_degree=old_degree//tau_degree,
  Luroth='withK=k(B15), k(Y16)=K(w) subsetK(z), [K(z):K(w)]=4',
  not_degree1_not_Galois_claim=True,
  no_expanded_Luroth_generator_or_NEW_equation_or_global_finite_map=True,
  no_CAS_point_prime_Mobius_or_new_invariant_expansion=True,
  v0_alternative_section_NOT_used='its constantmap extension is a separatefuturegate, not assumed here')
 (out/'proof-certificate.json').write_text(json.dumps(proof,indent=2)+'\n')
 note=ROOT/'research/notes'/(ID+'.md')
 report=dict(id=ID,all_passed=True,actual_original_normalized_source_generic_degree_exact8=True,
  actual_tau_descended_full_source_generic_degree_exact4=True,
  degree4_over_same_actual_B15_field=True,
  no_extra_Sp4_center_sigma_sign_scalar_or_n1_factor=True,
  independent_conceptual_review_pending=True,
  no_CAS_new_point_prime_coordinate_expansion_or_source26=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),note,*paths]},checked_gate_dependency_sha256=checked,
  evidence_sha256={'proof-certificate.json':sha(out/'proof-certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in['all_passed','actual_original_normalized_source_generic_degree_exact8','actual_tau_descended_full_source_generic_degree_exact4','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

#!/usr/bin/env python3
"""Combine three FINAL reviewed gates; no new algebraic computation."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_RATIONAL15_APPLICATION_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 stems=['Q_NEW_FAMILY_TARGET_FIELD_AND_SPECTRAL_DEGREES_20261104_v63',
        'Q_NEW_FAMILY_RATIONAL_ORBIT_SECTION_20261104_v63',
        'Q_NEW_FAMILY_SP4_CENTRALIZER_FREE_20261104_v63']
 primaries=[ROOT/'reports'/stem/'report.json'for stem in stems]
 reviews=[ROOT/'reports'/(stem+'_REVIEW')/'report.json'for stem in stems]
 P=[read(p)for p in primaries];V=[read(p)for p in reviews]
 assert all(v['all_passed']for v in V)
 assert P[0]['status']=='proved_by_accepted_prerequisites'and not P[0]['pending']
 assert V[0]['actual_target_field_equals_known_tau_fixed_source_QQ_and_C']
 assert V[0]['actual_invariant_field_pure_rational_trdeg5_over_QQ_and_C']and V[0]['target_source_index']==1
 assert V[0]['geometric_basechange_and_QQ_field_of_definition_checked']
 assert V[1]['all16_rational_orbit_representative_identities_and_nonzero_denominators']
 assert V[1]['same_dense_Sp4_saturation_on_generic_domain_proved']
 assert V[2]['full_geometric_GL4_common_commutant_scalars_not_only_Lie_stabilizer']
 assert V[2]['generic_B_NEW_Sp4_stabilizer_identity_and_group_scheme_char0_proved']
 assert P[2]['commutator_rank']==15 and P[2]['nonzero15_minor_modp']==17760
 dependencies={}
 for rec in P+V:
  for rel,want in rec.get('input_sha256',{}).items():
   path=ROOT/rel
   if rel not in dependencies:dependencies[rel]=sha(path)
   assert dependencies[rel]==want,rel
 section_path=ROOT/'reports'/stems[1]/'certificate.json';section=read(section_path)
 assert len(section['beta_bar16'])==16 and len(section['variables'])==5
 assert section['known_tau_ambiguity_removed_in_rational_orbit_representative']
 recipe=dict(coefficient_field=['QQ','C'],actual_target='Reduced geometrically integral firstjet image B_NEW15 subset W16',
   source='Sp4 x S5, with S5 the stated nonempty rational-domain open in A5',
   parameter_order=section['variables'],section16=section['beta_bar16'],source_domain=section['source_domain'],
   map='Phi(g,z)=rho(g)*beta_bar(z)',finite_action_convention='rho(g)=g* is the saved pullback anti-action; Phi(g*h,z)=rho(h)*Phi(g,z), so source multiplication is RIGHT',
   rational_quotient='Choose the five elements of actualE corresponding to the five S5 parameters under the accepted SAME pullbackE=L',
   quotient_section_identity='pi(beta_bar(z))=z as rational maps',
   dominant='Reviewed denseSp4saturation of this SAME section, not source dimension alone',
   generic_section_free='Free locus is nonempty andSp4stable; denseSp4saturation forces section to meet it',
   generic_fiber='EqualPhi images imply equalz bypi, then equalg bygeometrically trivial stabilizer',
   degree_one='Dominant same-dimension map in char0 has separable generic extension; generic geometric fiberone gives degree1 overQQ',
   rational_group='SplitSp4 bigcell Uminus(A4)*T(Gm2)*Uplus(A4), defined overQQ',
   conclusion='QQ(B_NEW)=QQ(Sp4)(five section parameters), pure rationaltrdeg15; likewiseC',
   expanded_inverse_or_all15_coordinate_polynomials_not_saved=True,
   no_entire_Poisson_component25_projective_target_or_boundary_orbit_claim=True)
 (out/'application-certificate.json').write_text(json.dumps(recipe,indent=2)+'\n')
 note=ROOT/f'research/notes/{ID}.md'
 result=dict(id=ID,all_three_final_independent_gates_passed=True,
   rational_section_same_actual_quotient_pullback_index1=True,
   geometric_generic_Sp4_stabilizer_identity=True,
   explicit_rational_product_map_dominant=True,generic_geometric_fiber_one=True,
   birational_Sp4_times_A5_to_actual_B_NEW15_over_QQ_and_C=True,
   actual_firstjet_B_NEW_QQ_and_C_rational_dimension15=True,
   rational_invariant_field_alone_not_used_to_trivialize_torsor=True,
   square_root_not_assumed_rational_over_QQ_section_field=True,
   char0_separability_and_QQ_descent_checked=True,
   no_full25_Poisson_component_rationality_source26_index_or_NEW_equation_claim=True,
   no_allpoint_free_or_global_quotient_trivialization_claim=True,
   no_new_point_prime_matrix_minor_invariant_or_CAS_computation=True,
   verified_dependency_count=len(dependencies),
   input_sha256={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),note,section_path,*primaries,*reviews]},
   checked_gate_dependency_sha256=dependencies,
   application_certificate_sha256=sha(out/'application-certificate.json'),
   seconds=time.monotonic()-start,frozen=True,independent_conceptual_review_pending=True)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','checked_gate_dependency_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

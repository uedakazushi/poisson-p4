#!/usr/bin/env python3
"""Independent proof gate for complete quadratic initial space and depth19.

No producer/source arithmetic replay: only fixed kernel-basis identities,
one normalized functional, existing result arrays, and accepted peer hashes.
"""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ORIGINAL_DEPTH19_INITIAL_SPACE_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out,source_peer):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
 f={
  'source_primary':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/report.json',
  'source_values':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
  'fullT':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/actual-minimal5-by825-Schur-operator.json',
  'linear':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/actual-linear-quadratic-T-certificate.json',
  'known2':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase2-certified-quadratic-input.json',
  'duals':'reports/AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/three-quadratic-duals.json',
  'zero620':'reports/AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/zero-linear620-source-combinations.json',
  'linear_input':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json'}
 src,val,T,lin,known,dual,ker,linear_input=[read(ROOT/v)for v in f.values()]
 assert sha(ROOT/f['source_primary'])=='5233957a7ab9dac28f16d9f51085fe375ffaac1d293b82d065f0381c0062bf92'
 peer=read(source_peer);assert peer['all_passed']
 assert sha(ROOT/f['source_primary'])in json.dumps(peer),'source peer must pin the exact immutable actual-source primary'
 assert src['all_passed']
 assert T['actual_minimal_generator_presentation_shape']==[5,825]
 assert T['actual_relation_order']==['168 literal d3 images','657 remaining actual source cycles']
 assert len(T['actual_five_by825_scalar_operator_root_ids'])==5 and all(len(v)==825 for v in T['actual_five_by825_scalar_operator_root_ids'])
 assert T['known_low_coefficients']['d3_constant_and_linear']=='zero'
 cert=lin['complete_linear_rank_certificate'];assert cert['rank']==37 and lin['actual_mu_Rbar_T']==5
 assert T['actual_remaining_source_labels657']==linear_input['actual657_minimal_relation_source_labels']
 factor=cert['all_column_factorization_in_normalized_basis'];assert len(factor)==657
 piv=set(cert['selected_original_columns']);assert len(piv)==37
 nonpiv=set(range(657))-piv;seen=set();assert len(ker['all620'])==620
 for rec in ker['all620']:
  coeff={i:Q(v)for i,v in rec['combination_in657']}
  outside=set(coeff)-piv;assert len(outside)==1
  own=next(iter(outside));assert own in nonpiv and own not in seen and coeff[own]==1;seen.add(own)
  total=[Q(0)]*37
  for col,c in coeff.items():
   for row,v in factor[col]:total[row]+=c*Q(v)
  assert not any(total)
 assert seen==nonpiv
 assert dual['five_targets']==[22,23,29,30,55] and len(dual['monomials36'])==36
 lam={i:Q(v)for i,v in dual['projection3_by180_rows_sparse'][2]};assert lam=={179:Q(1)}
 assert dual['monomials36'][35]==[7,7] and 179==4*36+35
 columns=known['all464_certified_quadratic_columns_sparse'];assert len(columns)==464
 source_labels=T['actual_remaining_source_labels657'];selected_labels={source_labels[i]for i in piv}
 labels=known['source_columns464'];assert len(labels)==464
 variable_labels=[(x['source_relation'],x['variable'])for x in labels if x['kind']=='variable_times_actual_linear_relation']
 assert len(variable_labels)==37*8 and set(variable_labels)=={(i,v)for i in selected_labels for v in range(8)}
 d3_labels=[x['literal_d3_column']for x in labels if x['kind']=='certified_actual_d3_quadratic']
 assert len(d3_labels)==168 and set(d3_labels)==set(range(168))
 for col in columns:assert sum((lam.get(i,Q(0))*Q(v)for i,v in col),Q(0))==0
 values=val['620_projected_source_functionals'];assert len(values)==620 and {r['index']for r in values}==set(range(620))
 for r in values:
  assert Q(r['three_known_oldB_P_and_pivot_correction_values'][2])==0
  assert r['three_unknown_division_second_functionals'][2]==[]
  assert r['three_native_X2_functionals'][2]==[]
 # Every degree2 initial relation is R2*a0+L*a1, with a0 in ker(L).
 # Its ker(L) is exactly168 d3 constants plus the proved620 cycle basis.
 assert 168+657==825 and 657-37==620 and 5*8-37==3 and 5+3+1==9
 rankp=ROOT/'reports/AMBIENT_B18_CRITICAL_CONTACT_AND_COMPARISON_RANK4_20261104_v63/report.json';rank=read(rankp)
 assert rank['all_passed'] and rank['rank_D_K']==4 and rank['dim_k_C_ell_torsion']==4
 assert rank['mu_D_K_formula']=='length_k(T)-4' and rank['origin_depth_cone_19_iff']=='length_k(T)>8'
 gs=[ROOT/'reports'/n/'report.json'for n in [
  'AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW',
  'AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63_REVIEW',
  'AMBIENT_B18_HONEST_T_EULER_OPERATOR_20261104_v63_REVIEW',
  'AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63_REVIEW',
  'AMBIENT_B18_SPECIALIZED_COMPARISON_T_PREFLIGHT_20261104_v63_REVIEW',
  'AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63_REVIEW',
  'AMBIENT_B18_FOUR_MODULE_DEPTH_CRITERION_20261104_v63_REVIEW',
  'AMBIENT_B18_T1_T2_QUADRATIC_PHI_20261104_v63_REVIEW',
  'AMBIENT_B18_T3_ZERO_STRUCTURAL_20261104_v63_REVIEW',
  'AMBIENT_B18_W1_TAIL_BINDING_20261104_v63_REVIEW',
  'AMBIENT_CA_CRITICAL_V1_CONTACT_20261104_v63_REVIEW',
  'AMBIENT_CA_CRITICAL_V1_CONTACT_ALGEBRA_20261104_v63_REVIEW',
  'AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63_REVIEW',
  'AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63_REVIEW']]
 assert all(read(p)['all_passed']for p in gs)
 proof={'complete_actual_presentation':[5,825],'cycle_columns':657,'d3_columns':168,'linear_rank':37,
  'kernel_constant_cycle_basis620_independent_complete_by_nonpivot_triangular_form':True,
  'all620_factorization_kernel_identities_checked':True,
  'quadratic_initial_relation_exhaustion':['variable multiples of all linear relations','all168d3 quadratic parts','all620zero-linear constant cycle combinations quadratic parts'],
  'formal_coefficient_proof':'For R*a in m²F, expand a=a0+a1+O(m²); La0=0 and (Ra)2=R2a0+La1. The complete finite relation module is closed. No additional degree2 relation exists.',
  'actual_surviving_dual_sparse':[[179,'1']],'representative':'fifth minimal target times rA²',
  'known464_columns_annihilated':True,'all620_actual_third_functionals_zero_no_unknown_source':True,
  'actual_source_semantics_independently_peer_gated':True,'HF_first_three_lower_bounds':[5,3,1],
  'length_T_mod_m3_lower_bound':9,'length_T_lower_bound':9,'mu_D_K_lower_bound':5,'rank_D_K':4,
  'torsion_D_K_nonzero':True,'H0_actual_D0_nonzero':True,'original_actual_intersection_cone_depth':19,'original_actual_intersection_projective_depth':18,
  'exact_total_T_length_or_H0_length_or_torsion_exponents_claimed':False,
  'explicit_D0_socle_representative_or_nil_cube_nonzero_claimed':False,
  'no_source_arithmetic_producer_GB_point_or_jet_duplicated':True}
 ep=out/'independent-initial-space-completeness-and-depth-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',source_peer,rankp]+[ROOT/v for v in f.values()]+gs
 report={'id':STEM,'all_passed':True,'status':'INDEPENDENT_complete_quadratic_initial_space_and_original_depth19_accepted',
  'source_peer_gate':str(source_peer.relative_to(ROOT)),'source_peer_sha256':sha(source_peer),
  'all_formal_coefficient_relations_to_degree2_exhausted':True,'third_dual_nonzero_on_full_gr2_T':True,
  'length_T_mod_m3_lower_bound':9,'length_T_lower_bound':9,'mu_D_K_lower_bound':5,'rank_D_K':4,
  'DVR_torsion_K_and_H0_D0_nonzero':True,'original_cone_depth':19,'original_projective_depth':18,
  'exact_T_H0_lengths_socle_class_nil_cube_or_exhaustion_claimed':False,
  'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{ep.name:sha(ep)},
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['id','all_passed','length_T_lower_bound','original_cone_depth','original_projective_depth','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer',required=True,type=Path);a=p.parse_args();main(a.output.resolve(),a.source_peer.resolve())

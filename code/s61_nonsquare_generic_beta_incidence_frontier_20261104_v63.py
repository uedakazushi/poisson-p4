#!/usr/bin/env python3
"""Source/statement audit and small monomial countermodel; no rank/point job."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='S61_NONSQUARE_GENERIC_BETA_INCIDENCE_FRONTIER_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def monomials(n,d):
    if n==1:return [(d,)]
    return [(i,)+m for i in range(d+1) for m in monomials(n-1,d-i)]
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    files=['reports/S61_BINARY_H11_RANK_DROP_RESIDUAL_FINITE_20261102.json','reports/S61_BINARY_NUISANCE_CONTRACTION_20261102.json','reports/S61_BINARY_NUISANCE_CONTRACTION_REVIEW_20261102.json','reports/S61_BINARY_NONSQUARE_FRAME_GAP_20261103_v62.json','reports/S61_BINARY_NONSQUARE_FRAME_GAP_20261103_v62_REVIEW.json','reports/S61_NONSQUARE_FIXED_QUARTIC_CODIMENSION_TRANSFER_20261103_v63_REVIEW.json','reports/S61_BINARY_HYPERPLANE_CUBIC_CODIMENSION_20261102.json','reports/S61_NONSQUARE_RAW_H_B_FINITE_BAD_LOCUS_20261104_v63_REVIEW/report.json','reports/S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63_REVIEW/report.json','reports/S61_NONSQUARE_THETA_P80_REPLACEMENT_20261104_v63_REVIEW/report.json','reports/S61_NONSQUARE_THETA_P100_REPLACEMENT_20261104_v63_REVIEW/report.json','reports/S61_BINARY_NONSQUARE_HEAVY_CURVE_CUBIC_20261103_v62/report.json','reports/S61_BINARY_NONSQUARE_HEAVY_CURVE_QUARTIC_20261103_v62/report.json']
    rs={f:read(ROOT/f) for f in files};assert all(r['all_passed'] for r in rs.values())
    checked={}
    for f,r in rs.items():
        for p,h in r.get('input_sha256',{}).items():assert sha(ROOT/p)==h,p;checked[p]=h
    H=rs[files[0]];A=rs[files[1]];F=rs[files[3]];Fr=rs[files[4]];Q=rs[files[5]];C=rs[files[6]];B=rs[files[7]]
    assert H['unrestricted_normalized_source_dimension']==11 and H['unipotent_group_dimension']==3 and H['core_fixed_under_unipotents'] and H['unipotent_fixes_Q0_exactly']
    assert H['nuisance_indices_one_based']==[5,25,26,27,34,37,48,51] and H['nuisance_weights']==[2,4,2,4,2,2,2,2]
    assert H['relative_base_dimension']==14 and H['open_core_complement_dimension_upper']==1 and H['residual_joint_dimension_upper']==16
    assert A['equations_checked']==135 and A['all_equations_individually_weight_homogeneous'] and A['homogenization_q_degree']==2 and A['normalized_chart_t_weight']==-6
    assert F['nuisance_indices_one_based']==H['nuisance_indices_one_based'] and F['full_source_differential_shape']==[16,19] and F['beta_image_dimension']==Fr['full_beta_image_dimension']==15
    assert Q['fixed319_rank_defect_codimension_at_least2_on_U0'] and Q['same_full135_including_G6_all7_weight4_prolongations_retained'] and Q['fixed_reduced_cut_determinant_minus_C_over25_unit_on_U0']
    assert C['cubic_rank_drop_dimension_upper']==1 and B['actual_B_rankbad_finite_or_empty']
    assert rs[files[9]]['ALL_P80_full_cubic57_proved'] and rs[files[10]]['ALL_P100_full_cubic57_proved']
    # Explicit logical countermodel with the same low weighted ranks. It is
    # a monomial algebra, never an actual Poisson source or new parameter point.
    light=[tuple(int(i==j) for i in range(7)) for j in range(7)]
    q=[tuple(2*x for x in light[0]),tuple(light[0][i]+light[1][i] for i in range(7))]+[tuple(2*x for x in light[j]) for j in range(2,6)]
    products={tuple(a+b for a,b in zip(qq,l)) for qq in q for l in light};assert len(products)==41
    pure3=monomials(7,3);x7cube=(0,0,0,0,0,0,3);extra=[m for m in pure3 if m not in products and m!=x7cube][:10];assert len(extra)==10
    genuine=[m+(0,) for m in extra]+[l+(1,) for l in light[:6]];all3={m+(0,) for m in products}|set(genuine);assert len(genuine)==16 and len(all3)==57
    witness=(0,0,0,0,0,0,4,0);gens=[m+(0,) for m in q]+genuine+[(0,0,0,0,0,0,0,2)]
    assert all(not all(a<=b for a,b in zip(g,witness)) for g in gens)
    toy={'actual_Jacobi_example':False,'weights':[1,1,1,1,1,1,1,2],'quadrics':q,'genuine16':genuine,'extra_known_heavy_square_relation':gens[-1],'product_rank':41,'heavy_cubic_rank':6,'full_cubic_rank':57,'low_weight_HF':[1,7,29-6,91-57],'nonzero_weight4_witness':witness,'quartic239_surjectivity_does_not_follow':True}
    ep=out/'source-budget-conditional-theorem-and-countermodel.json';evidence={'actual_source':{'primitive_W16_indices':[4,5,12,13,14,25,26,27,34,35,36,37,48,49,50,51],'Q0':'x*w*d_y wedge d_z','fixed_coordinates':{'u12':0,'u13':0,'u4':'a','u14':'(e+1)/5','u35':1,'u36':'d','u49':'e','u50':-5},'nuisance_indices':H['nuisance_indices_one_based'],'nuisance_weights':H['nuisance_weights'],'group_dimension':10,'verified_core_fixing_free_group_fiber_dimension':3},'joint_dimension_bounds_by_core_dimension':{str(c):10+8+c-3 for c in range(4)},'one_proper_bad_source_on_core_curve_suffices':{'core_curve_times_nuisance_dimension':9,'proper_closed_bad_source_dimension_upper':8,'joint_bad_image_upper':15},'full_common_frame_gap_beta_dimension':15,'common_gap_rank_identity_not_transferred_to_other_curves':True,'finite_Bbad_requires_no_new_quartic_rank_test_for_generic_beta_exclusion':True,'old_accepted_heavy_curve_quartic_and_H11_relative_theorem_reused_no_arithmetic_rerun':True,'rawH7_open_unknown_cubic_or_quartic_defect_curves_not_excluded':True,'whole_nonsquare_generic_beta_removal_not_proved_by_current_strata':True,'all_FULLU_quartic239_not_necessary_for_conditional_generic_removal':True,'whole_affine_S61_exhaustion_not_inferred':True,'toy_only':toy}
    ep.write_text(json.dumps(evidence,indent=2)+'\n')
    notes=['research/notes/S61_BINARY_H11_RANK_DROP_RESIDUAL_FINITE_20261102.md','research/notes/S61_BINARY_NUISANCE_CONTRACTION_20261102.md','research/notes/S61_BINARY_NONSQUARE_FRAME_GAP_20261103_v62.md','research/notes/S61_BINARY_NONSQUARE_HEAVY_CURVE_CUBIC_REVIEW_20261103_v62.md','research/notes/S61_BINARY_NONSQUARE_HEAVY_CURVE_QUARTIC_20261103_v62.md','code/s61_binary_H11_rank_drop_residual_finite_20261102.py','code/s61_binary_nonsquare_frame_gap_20261103_v62.py',f'research/notes/{STEM}.md']
    inputs=[Path(__file__),*[ROOT/f for f in files+notes]]
    report={'id':STEM,'all_passed':True,'frozen':True,'verification_kind':'Read-only source/accepted-statement audit, conceptual dimension proof and monomial logical countermodel only','exact_H11_source8_and_N3_budget_bound':True,'finite_core_exception_beta_image_le15':True,'one_proper_full_relative_bad_section_per_core_curve_suffices':True,'present_information_does_not_prove_generic_removal_of_all_nonsquare':True,'no_new_rank_source_point_prime_GB_or_full_source_expansion':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},'checked_prior_source_hashes':checked,'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['all_passed','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

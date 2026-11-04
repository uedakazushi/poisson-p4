#!/usr/bin/env python3
"""Finite initial-space/depth corollary; source independent review pending.

No coefficient producer, source solve, new rank search, or high source jet.
"""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ORIGINAL_DEPTH19_COROLLARY_20261104_v63'
def dump(out,name,x):(out/name).write_text(json.dumps(x,indent=2)+'\n')
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    deps={
      'complete5':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/report.json',
      'complete5_peer':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW/report.json',
      'linear':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/report.json',
      'linear_peer':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63_REVIEW/report.json',
      'duals':'reports/AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/report.json',
      'root_adjoint':'reports/AMBIENT_B18_T_NEXT3_ROOT_ATOM_ADJOINT_20261104_v63/report.json',
      'cache_transport':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/report.json',
      'projected_m1':'reports/AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/report.json',
      'rankK4_depth_criterion':'reports/AMBIENT_B18_CRITICAL_CONTACT_AND_COMPARISON_RANK4_20261104_v63/report.json'}
    D={k:read(ROOT/p) for k,p in deps.items()}
    assert all(r['all_passed'] for r in D.values())
    assert D['complete5']['actual_complete_presentation']==[5,825]
    assert D['complete5_peer']['full_ring_111_unit_cancellations_and_source_completeness_verified']
    assert D['linear']['actual_length_T_mod_m2T']==8 and D['linear']['actual_linear_relation_rank']==37
    assert D['cache_transport']['active_division_second_pair_count']==0
    assert D['rankK4_depth_criterion']['rank_D_K']==4
    assert D['rankK4_depth_criterion']['dim_k_C_ell_torsion']==4
    assert D['rankK4_depth_criterion']['origin_depth_cone_19_iff']=='length_k(T)>8'
    for name,r in D.items():
        for p,h in r.get('input_sha256',{}).items():assert sha(ROOT/p)==h,(name,p)
        for h in r.get('input_hashes',{}).values():
            if isinstance(h,dict) and 'path' in h:assert sha(ROOT/h['path'])==h['sha256']
    paths={
      'duals':'reports/AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/three-quadratic-duals.json',
      '620':'reports/AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/zero-linear620-source-combinations.json',
      'known_quadratics':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase2-certified-quadratic-input.json',
      'linear_certificate':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/actual-linear-quadratic-T-certificate.json',
      'actual_values':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
      'source_adjoint':'reports/AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/projected-source-X2-sensitivity.json'}
    Q,C,K,L,V,S=[read(ROOT/p) for p in paths.values()]
    dual=rec(Q['projection3_by180_rows_sparse'][2]);assert dual=={179:F(1)}
    assert Q['monomials36'][35]==[7,7] and Q['five_targets'][4]==55
    assert rec(S['three_rows_sparse'][2])=={2015:F(1)} # 56th X row,r²
    assert len(K['all464_certified_quadratic_columns_sparse'])==464
    for col in K['all464_certified_quadratic_columns_sparse']:assert not rec(col).get(179,F(0))
    cert=L['complete_linear_rank_certificate'];selected=set(cert['selected_original_columns']);assert len(selected)==37
    zeros=[i for i in range(657) if i not in selected]
    assert len(C['all620'])==len(V['620_projected_source_functionals'])==len(zeros)==620
    for expected,a,b in zip(zeros,C['all620'],V['620_projected_source_functionals']):
        assert a['linear_zero_receipt']==expected
        c=rec(a['combination_in657']);assert c[expected]==1 and set(c)-{expected}<=selected
        assert all(i<=expected for i in c)
        assert b['three_unknown_division_second_functionals'][2]==[]
        assert b['three_native_X2_functionals'][2]==[]
        assert F(b['three_known_oldB_P_and_pivot_correction_values'][2])==0
    survival={'nonzero_quadratic_dual':[[179,'1']],'coordinate':'target4*rA^2, fullV source55',
      'value_on_its_quadratic_basis_vector':'1','annihilates_all464_known_initial_relations':True,
      'all620_form_a_basis_of_kernel_of_actual_linear_matrix':True,
      'annihilates_all620_ACTUAL_projected_cycle_quadratics':True,
      'native_bottom8_literal_zero_so_no_native_X2_value_needed':True,
      'all_opaque_division_second_coefficients_cancel_before_value_use':True,
      'complete_presentation_relation_split':[657,168],
      'entire_quadratic_initial_relation_space_containment':'linear relations times8vars +168 d3 quadratic columns +620 zero-linear cycle quadratic columns',
      'dim_gr2_T_lower_bound':1,'length_T_mod_m2':8,'length_T_mod_m3_lower_bound':9,
      'length_T_lower_bound':9,'mu_D_K_lower_bound':5,'rank_D_K':4,
      'DVR_torsion_K_nonzero':True,'H0_original_D0_nonzero':True,
      'original_cone_depth':19,'original_projective_depth':18,
      'actual_total_T_length':None,'actual_H0_length':None,'nilpotence_cube_status':None,
      'explicit_original_D0_socle_representative':None,
      'independent_source_transport_review_status':'PENDING','independent_initial_space_depth_review_status':'PENDING'}
    dump(out,'third-quadratic-class-survival-and-depth.json',survival)
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'status':'PRIMARY_PROOF_PENDING_INDEPENDENT_SOURCE_AND_INITIAL_SPACE_REVIEWS',
      'original_cone_depth':19,'original_projective_depth':18,'length_T_mod_m3_lower_bound':9,
      'length_T_exact':None,'H0_length_exact':None,'full_nil_cube_computed':False,
      'generic_exhaustion_fullclassification_or_fullquartic_solved':False,
      'source_review_pending':True,'pure_initial_space_depth_review_pending':True,
      'input_sha256':{p:sha(ROOT/p) for p in list(deps.values())+list(paths.values())},
      'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{'third-quadratic-class-survival-and-depth.json':sha(out/'third-quadratic-class-survival-and-depth.json')},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','status','length_T_mod_m3_lower_bound','original_cone_depth','original_projective_depth','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)

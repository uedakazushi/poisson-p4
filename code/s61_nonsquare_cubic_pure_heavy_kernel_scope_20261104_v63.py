#!/usr/bin/env python3
"""Frozen stage/projection scope audit; no coefficient/rank recomputation."""
import argparse
import hashlib
import json
import time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_CUBIC_PURE_HEAVY_KERNEL_SCOPE_20261104_v63'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def main(out):
    began=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    read=lambda p:json.loads(p.read_text())
    b=ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63'
    old=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63'
    aug=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63'
    final=ROOT/'reports/S61_NONSQUARE_CUBIC_LAST_Q9_UNIT_COMPLETION_20261104_v63'
    data=read(aug/'input.json');cp=read(old/'checkpoint.json');bc=read(b/'certificate.json');br=read(b/'report.json');fc=read(final/'certificate.json')
    assert br['generic_rank_exactly10_proved']and br['six_remainder_test_completed']
    for name in['S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63_REVIEW',
                 'S61_NONSQUARE_CUBIC_LAST_Q9_UNIT_COMPLETION_20261104_v63_REVIEW']:
        assert read(ROOT/'reports'/name/'report.json')['all_passed']
    assert data['target31_recipe']['retained31_full40_columns'][-5:]==[35,36,37,38,39]
    assert data['target31_recipe']['NO_heavy5_Schur_has_been_applied']
    assert cp['selected_minor_rows']==[8,1,2,3,0,6,7,11]
    assert cp['selected_minor_columns']==[20,10,16,12,0,23,24,26]
    assert bc['heavy5_fixed_genuine_relative_rows']==[1,2,8,10,11]
    assert '10quad and5heavy' in bc['R11_scope']
    headers={}
    for name,expected in [('completed-afterquad16x21.txt','16 21'),('completed-R11x16.txt','11 16')]:
        with(b/name).open()as f:headers[name]=f.readline().strip()
        assert headers[name]==expected
    pure_support=[j for j in range(26)if cp['matrix25x32'][11][j][0]]
    assert len(pure_support)==19
    assert len(fc['new21_target_labels'])==21
    assert sum(j<26 for j in fc['new21_target_labels'])==16
    assert sum(j>=26 for j in fc['new21_target_labels'])==5
    proof={'initial_augmented_pure_heavy_targets':[26,6],
           'old_eight_selected_pivots':list(zip(cp['selected_minor_rows'],cp['selected_minor_columns'])),
           'old_eight_classification':{'quadratic_pure':7,'genuine_heavy':1},
           'genuine_pivot':{'source_label':11,'target_label':26,'target':'H*L','not_pure':True,
                            'normalized_pure_support19':pure_support},
           'later_three_pivots_all_quadratic_pure':[[4,5],[5,21],[9,1]],
           'total_eleven_classification':{'pure':10,'heavy':1},
           'final_pure_heavy_target_counts':[16,5],
           'canonical11_headers':headers,
           'canonical11_scope':bc['R11_scope'],
           'canonical_generic_rank10_bound_applies_after_heavy5_Schur_only':True,
           'canonical_rank10_not_whole_initial_pure_projection_upper':True,
           'no_all10_pure_minors_zero_inference':True,
           'correct_kernel_criterion':{
               'field_map':'M=(P,H):k14->k16+ k5',
               'K':'ker(H), dimension14-r where r=rank(H)',
               'exact_rank_formula':'rank(M)=r+rank(P|K)',
               'full_injectivity':'rank(P|K)=14-r',
               'rank5_heavy_open':'K9; nine-by16 kernel-pure map must be injective',
               'heavy_rank_atmost4':'K dimension at least10; P16 may still inject it, so no rank13 conclusion'},
           'constant_rank_chart_recipe':{
               'choose':'r source rows I and r heavy columns C with Delta=det(H[I,C]) nonzero',
               'kernel_vectors':'Delta*e_u - sum_i (H[u,C]*adj(H[I,C]))_i*e_i, for u outside I',
               'condition':'On rank(H)=r these vectors annihilate all five heavy columns and form a kernel frame after inverting Delta',
               'actual_pure_rows':'Delta*P_u - sum_i lambda_ui*P_i; retain finite DAG, no coefficient expansion',
               'no_global_unknown_factor_unit':'Delta is only a chart denominator unless separately certified unit'},
           'logical_countermodel':{
               'matrix_pure2_heavy1':[['1','0','x'],['0','1','0']],
               'pure_projection_rank':'2 everywhere',
               'heavy_rank':'1 on D(x),0 at x=0',
               'after_generic_heavy_Schur':'one pure row remains, of rank1',
               'full_rank':'2 everywhere because the pure2 minor is1',
               'purpose':'A generic Schur remainder rank is not an upper bound on the whole fixed pure projection.'},
           'Bbad_does_not_by_itself_prove_cubic57_failure':True,
           'no_accepted_ledger_theorem_withdrawal':'The invalid implication was only proposed for audit and was never adopted.'}
    (out/'certificate.json').write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),b/'report.json',b/'certificate.json',
           b/'completed-afterquad16x21.txt',b/'completed-R11x16.txt',old/'checkpoint.json',
           aug/'input.json',final/'certificate.json',
           ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63_REVIEW/report.json',
           ROOT/'reports/S61_NONSQUARE_CUBIC_LAST_Q9_UNIT_COMPLETION_20261104_v63_REVIEW/report.json']
    r={'id':ID,'all_stage_header_source_and_pivot_scope_checks_passed':True,
       'all11_pivots_are_pure':False,'pure_pivots':10,'heavy_pivots':1,
       'whole_final_pure_projection_upper9_NOT_proved':True,
       'B_rankdrop_implies_cubic_rankdrop_NOT_proved':True,
       'correct_pointwise_kernel_H_criterion_recorded':True,
       'no_new_coefficient_expansion_rank_point_GB_or_minor_calculation':True,
       'all_accepted_full14x21_H5_and_generic_certificates_unchanged':True,
       'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
       'certificate_sha256':sha(out/'certificate.json'),
       'seconds':time.monotonic()-began,'frozen':True}
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output.resolve())

#!/usr/bin/env python3
"""Independent pure-math fixed-cut codimension audit; no coefficient computation."""
import argparse
from math import comb
import json
from pathlib import Path
import time

import ambient_B18_zero_linear_contact4_review_20261103_v63 as util

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_FIXED_QUARTIC_CODIMENSION_TRANSFER_20261103_v63_REVIEW'


def checked(path,require_passed=True):
    data=json.loads(path.read_text())
    if require_passed:assert data['all_passed'] is True
    for name,digest in data.get('input_sha256',{}).items():assert util.sha(ROOT/name)==digest,name
    for name,digest in data.get('evidence_sha256',{}).items():
        p=ROOT/name if (ROOT/name).exists() else path.parent/name
        assert util.sha(p)==digest,name
    return data


def main(output):
    start=time.monotonic();assert not output.exists()
    ap=ROOT/'reports/S61_NONSQUARE_FIXED_QUARTIC_CODIMENSION_TRANSFER_20261103_v63/report.json'
    author=checked(ap,False);assert author['all_static_premise_hash_checks_passed']
    stems=['S61_BINARY_HYPERPLANE_TRANSVERSE_CODIMENSION_20261102',
        'S61_BINARY_HYPERPLANE_TRANSVERSE_CODIMENSION_REVIEW_20261102',
        'S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102',
        'S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102',
        'S61_BINARY_HYPERPLANE_CUBIC_CODIMENSION_20261102',
        'S61_BINARY_HYPERPLANE_GENERIC_RANKS_20261102']
    paths=[ROOT/'reports'/f'{s}.json' for s in stems]
    tr,oldreview,chart,q,c,g=[checked(p) for p in paths]
    assert tr['coordinate_free_existence_of_transverse']
    assert tr['relative_completion_base_dimension']==11 and tr['relative_module_rank_on_zero_nuisance_section']==65
    assert tr['exceptional_core_dimension_upper']==1
    assert oldreview['proof_audit_accepted'] and oldreview['completed_regular_base_dimension']==11 and not oldreview['H11_in_reviewed_scope']
    assert chart['family_parameters']==[1,2,5]
    assert chart['family_parameter_minor']=='-(100*a + 18*e**2 - 9*e - 2)/25'
    assert chart['normal_weights34']==[1,1,1,1,1,2,1,1]
    assert chart['full_linear_kernel_checked'] and chart['selected_right_inverse_checked']
    assert chart['implicit_solution_exact_degree_bound']==4 and chart['full_remaining_equation_weight_bound']==6
    assert q['bad_locus_dimension_upper_on_three_core_open']==1 and q['quadratic_row_count']==6
    assert c['cubic_rank_drop_dimension_upper']==1 and c['all_core_cubic_rank_upper']==57
    assert g['all_core_cubic_rank_upper']==57
    # Independent finite counts and row-index retention, no coefficient matrices evaluated.
    weighted=[1,7,comb(8,2)+1,comb(9,3)+7,comb(10,4)+comb(8,2)+1]
    assert weighted==[1,7,29,91,239]
    low=[weighted[0],weighted[1],weighted[2]-6,weighted[3]-57]
    assert low==[1,7,23,34] and sum(low)==65
    rp=ROOT/'reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/report.json'
    rowsreport=checked(rp)
    assert rowsreport['row_span_equality_both_directions'] and rowsreport['triangular_deleted_minor']=='-a^14'
    rcpath=rp.parent/'certificate.json';rc=json.loads(rcpath.read_text())
    assert rc['original_shape']==[354,239] and rc['retained_shape']==[319,239]
    assert len(rc['genuine_quartic_original_rows_retained'])==33
    G6=[i for i,row in enumerate(rc['source_rows']) if row[0]=='cubic_normalized' and row[2]==134]
    assert len(G6)==7 and set(G6)<=set(rc['retained_row_indices'])
    assert 6*29+21*7+33==354 and 6*29+16*7+33==319
    helper=ROOT/'code/s61_flag_transverse_obstruction_20261101.py'
    text=helper.read_text();assert 'def specialize(model,beta,extra=False)' in text and "['universal_new_equation']] if extra else []" in text
    source=ROOT/'code/s61_binary_hyperplane_weighted_chart_20261102.py'
    text=source.read_text();assert 'specialize(model,beta,True)' in text and 'assert s.cancel(fdet+Ccore/25)==0' in text
    assert 'T34=(W.to_DM()*T34.to_DM()).to_Matrix()' in text and 'U34=(W.to_DM()*U34.to_DM()).to_Matrix()' in text
    notes=[ROOT/'research/notes/S61_BINARY_HYPERPLANE_TRANSVERSE_CODIMENSION_20261102.md',
        ROOT/'research/notes/S61_BINARY_HYPERPLANE_TRANSVERSE_CODIMENSION_REVIEW_20261102.md',
        ROOT/'research/notes/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.md',
        ROOT/'research/notes/S61_BINARY_HYPERPLANE_GENERIC_RANKS_20261102.md']
    inputs=[Path(__file__),Path(util.__file__),ROOT/f'research/notes/{STEM}.md',ap,*paths,rp,rcpath,helper,source,*notes]
    result=dict(id=STEM,all_passed=True,
        verification_kind='Independent conceptual proof review and static finite source/index/hash checks only; no coefficient jobs, GB, field or point computations.',
        old_H8_theorem_is_full_relative_free65_not_only_cut_length=True,
        generic_factor_length65_consumes_full_rank_and_nil_ideal_power65_zero=True,
        admissible_new_cut_finiteness_by_nil_filtration_and_free65_by_CM_proved=True,
        fixed_cut_is_actual_W_inverse_free34_coordinates_not_unbound_raw_indices=True,
        fixed_reduced_cut_determinant_minus_C_over25_unit_on_U0=True,
        exact_positive_weight_polynomial_chart_not_only_finite_jet=True,
        low_weight_exact_Hilbert_function=[1,7,23,34],low_weight_total_length=65,
        same_full135_including_G6_all7_weight4_prolongations_retained=True,
        fixed_full354_and319_weight4_image_equal_on_aFBC=True,
        fixed319_rank_defect_codimension_at_least2_on_U0=True,
        all_height1_primes_geometric_and_rational_covered=True,
        C0_outside_fixed_chart_and_only_old_existential_cut_theorem_used_there=True,
        eighty_pointwise_independent_kernel_subbundle_would_imply_global_surjectivity=True,
        eighty_directions_constructed_or_full239_everywhere_or_ambient_removal_claimed=False,
        C122_corollary_only_conditional_on_its_separate_exact_frame_rank_identity=True,
        scope='Only nuisance-zero H8 and the coherent fixed weighted chart on U0=Spec QQ[a,e,d,1/(aFBC)]; no H11 or whole S61 exhaustion.',
        input_sha256={str(p.relative_to(ROOT)):util.sha(p) for p in inputs},seconds=time.monotonic()-start)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    main(p.parse_args().output)

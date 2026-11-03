#!/usr/bin/env python3
"""Audit the completed c=2a pencil and marked-square source bridge.

The geometric implication is proved in the companion note. This producer
checks the exact prerequisites and every recursively recorded input hash.
It does not infer characteristic-zero radicals from finite-field radicals.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_marked_sigma_finite_cover_20261103_v62 import verify as conditional_verify

ROOT=Path(__file__).resolve().parents[1]


def verify():
    dependencies={};visited=set();hash_checks=[]
    def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    def audit(path):
        path=path.resolve()
        if path in visited:return
        visited.add(path)
        if path.suffix!='.json':return
        data=json.loads(path.read_text())
        if not isinstance(data,dict):return
        for key,base in [('input_sha256',ROOT),('evidence_sha256',path.parent),('output_sha256',path.parent)]:
            for name,expected in data.get(key,{}).items():
                target=base/name
                assert sha(target)==expected,(str(path),key,name)
                hash_checks.append([str(path.relative_to(ROOT)),key,name])
                if key=='input_sha256':audit(target)
    def load(name):
        path=ROOT/'reports'/name;data=json.loads(path.read_text());assert data['all_passed']
        audit(path);dependencies[str(path.relative_to(ROOT))]=sha(path);return data
    def evidence(name):
        path=ROOT/name;dependencies[name]=sha(path)

    conditional=load('S61_MARKED_SIGMA_FINITE_COVER_20261103_v62.json')
    assert conditional==conditional_verify()
    assert conditional['conditional_whole_pencil_implication_proved']
    assert conditional['remaining_conditions_not_claimed_here']==[
        'The OLD h^2+2=0 plane special-fiber discriminant-zero junction has a mixed-characteristic ambient nilpotent neighborhood.']
    mixed=load('S61_MARKED_SIGMA_OLD_DOUBLE_MIXED_CERTIFICATE_20261103_v62_FINAL/report.json')
    assert len(mixed['cases'])==1
    cut=mixed['cases'][0]
    assert cut['cut']=='delta_u4+delta_u20=0'and cut['cut_count']==17
    assert cut['ambient_Jacobian_rank']==41 and cut['normal_dimension']==3 and cut['quadratic_rank']==5
    assert cut['mixed_substitution_verified_in_all135_tangent_directions']
    assert cut['jacobian_minor_norm_prime_unit']and cut['quadratic_minor_norm_prime_unit']
    assert cut['scalar_denominators_prime_units']and cut['checked_scalar_count']==16248
    for root in ['root3083','root28920']:
        assert 'LENGTH=5'in cut['outputs'][root]and 'DIM=0'in cut['outputs'][root]
    mixed_review=load('S61_MARKED_SIGMA_OLD_DOUBLE_MIXED_CERTIFICATE_20261103_v62_REVIEW.json')
    assert mixed_review['independent_QQ_and_two_Fp_Groebner_length']==5
    assert mixed_review['exact_mixed_affine_slice_in_all62_directions']
    assert mixed_review['independent_original135_multivariate_Taylor_coefficients_match']
    family=load('S61_MARKED_SIGMA_OLD_DOUBLE_JUNCTION_FAMILY_20261103_v62.json')
    assert family['all135_identically_zero']and family['t_identically_zero']
    assert family['constant_beta']and family['coefficient_denominators_prime_units']
    lower2=load('S61_MARKED_SIGMA_OLD_DOUBLE_QUOTIENT_20261103_v62_FINAL/report.json')
    assert lower2['exact_dual_number_thickening_with_t_zero']and lower2['same_incidence_generic_length_lower']==2
    lower2_review=load('S61_MARKED_SIGMA_OLD_DOUBLE_QUOTIENT_20261103_v62_REVIEW.json')
    lower3=load('S61_NONBINARY_DOUBLE_ROOT_OPEN_20261101.json')
    assert lower3['completed_ring_rank']==3 and lower3['homogenizer_zero']
    geometry=load('S61_MARKED_SIGMA_OLD_PLANE_GEOMETRY_20261103_v62.json')
    assert all(r['full135_vanishes_on_entire_linear_H']for r in geometry['cases'])
    double=next(r for r in geometry['cases']if r['case']=='double-root')
    assert double['H_dimension']==10 and double['orbit_dimension']==7 and double['total_incidence_dimension']==17
    sign=load('S61_MARKED_SIGMA_H_SIGN_TRANSPORT_20261103_v62.json')
    assert sign['all135_homogenized_equations_preserved_up_to_sign']
    assert sign['rho_fixed']and sign['a_fixed']and sign['c_fixed']and sign['h_character']==-1
    assert sign['q25_fixed']and sign['three_orbit_cuts_and_q32_base_preserved']and sign['homogenizer_character']==1
    original=load('S61_MARKED_SIGMA_RESIDUAL_FINITE_20261102.json')
    frame=original['frame'];assert frame['ambient_orbit_cuts']==[0,19,27]
    assert frame['nonpositive_q_base']==32 and frame['orbit_normal_minor']==4
    assert all(w<=0 for w in frame['H15_base_weights'][1:])
    source=load('S61_MARKED_SIGMA_PROJECTIVE_CORE_20261102.json')
    assert source['projective_weights']==[2,2,3]
    assert source['nuisance_weights']==[2,2,4,2,4,2,2,2]
    assert source['finite_core_subset_beta_image_dimension_bound']==15 and source['zero_core_beta_image_dimension_bound']==14
    load('S61_MARKED_SIGMA_PROJECTIVE_CORE_PROOF_REVIEW_20261102.json')
    intersection=load('S61_MARKED_SIGMA_INTERSECTION_REMOVAL_20261102.json')
    assert intersection['all_general_beta_marked_b1_zero_removed']
    model=S61();one=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==1]
    base_beta=['u'+str(i+1)for i in one if i not in [3,4]];assert len(base_beta)==14
    note='research/notes/S61_MARKED_SIGMA_FULL_PENCIL_BRIDGE_20261103_v62.md'
    for path in [str(Path(__file__).relative_to(ROOT)),note,
                 'research/notes/S61_MARKED_SIGMA_OLD_DOUBLE_JUNCTION_20261103_v62.md',
                 'reports/S61_MARKED_SIGMA_OLD_DOUBLE_MIXED_CERTIFICATE_20261103_v62_FINAL/certificate.json']:
        evidence(path)
    return dict(id='S61_MARKED_SIGMA_FULL_PENCIL_BRIDGE_20261103_v62',all_passed=True,
                proof_kind='Exact recorded arithmetic prerequisites plus the geometric proof in the companion note',
                prime=32003,all_recorded_dependency_hashes_current=True,recursive_hash_checks=len(hash_checks),
                nonpositive_base_beta14=base_beta,additional_nonpositive_base='q32; the DVR uniformizer has weight zero',
                finite_orbit_cuts=[0,19,27],infinity_orbit_cuts=[0,9,27],group_restoration_minors=[4,-2],
                old_double_mixed_transverse_length_upper=5,same_incidence_generic_length_lowers=[2,3],
                old_double_raw_ambient_free_rank_over_regular17cut_base=5,old_double_ambient_t_zero_over_V=True,
                old_plane_all_geometric_special_directions_covered=True,
                old_core_positive_central_support_exhausted=True,
                actual_saturated_support_to_origin_lemma_applied=True,
                all_geometric_special_points_of_test_pencil_removed=True,
                proper_test_cover='P1_V -> P(2,2,3)_V, [s:t] -> [s^2:2s^2:t^3]',
                whole_characteristic_zero_test_pencil_removed=True,
                ample_Cartier_equation='(c-2a)^3 in O(6)',
                bad_projective_core_dimension_at_most=0,fixed_flag_bad_source_dimension_at_most=11,
                joint_bad_incidence_dimension_at_most=15,zero_core_joint_dimension_at_most=14,
                whole_marked_square_stratum_removed_for_general_beta=True,
                general_five_points_exhaustion_claimed=False,global_component_upper_bound_claimed=False,
                finite_field_radical_promoted_to_characteristic_zero=False,
                input_sha256=dependencies)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify()
    with args.output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print('PASS',report['id'],report['recursive_hash_checks'],'recursive dependency hashes')

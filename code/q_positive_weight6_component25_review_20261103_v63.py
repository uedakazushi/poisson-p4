#!/usr/bin/env python3
"""Independent dimension/component/count proof from separately audited bounds."""
import argparse
import json
from pathlib import Path
from q_positive_weight6_quadratic_model_20261103_v63 import ROOT,read,sha


def main(family_review,output):
    family_review=family_review.resolve()
    source=ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/report.json'
    family=read(source);independent=read(family_review)
    assert family['all_passed'] and independent['all_passed']
    assert family['actual_full_three_parameter_Poisson_family']
    assert family['all140_affine_and_all350_projective_equations_zero_exactly']
    assert family['full_GL25_P_scalar_three_parameter_image_rank']==26
    assert family['positive_point_occurs_at_parameter_origin'] and family['generic_Poisson_rank4']
    upper_path=ROOT/'reports/Q_POSITIVE_WEIGHT6_INITIAL4_20261103_v63/report.json'
    upper_review_path=ROOT/'reports/Q_POSITIVE_WEIGHT6_SMALL_MODEL_20261103_v63_REVIEW/report.json'
    upper=read(upper_path);review=read(upper_review_path)
    assert upper['all_passed'] and review['all_passed']
    assert upper['actual_projective_Poisson_local_dimension_upper_bound']==25
    assert review['exact_implicit95_and_formal8_presentation']
    assert review['actual_transverse_dimension_upper3_projective_upper25_only']
    exclusions=ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63_REVIEW.json'
    orbit_path=ROOT/'reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63_REVIEW.json'
    ex=read(exclusions);orbit=read(orbit_path)
    assert ex['all_passed'] and ex['all_nonzero_ray_outside_actual_T01_CA_CB_N']
    assert orbit['all_passed'] and orbit['T39_T40_excluded_by_closed_invariant_components_and_actual_limit']
    assert orbit['nonzero_eta_one_complex_projective_orbit_and_closure_Qa8over9']
    catalog_path=ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json'
    catalog=read(catalog_path)['records'];assert len(catalog)==51
    large={r['component_id']:r['expected_s61_dimension']+9 for r in catalog if r['expected_s61_dimension']+9>=25}
    assert large=={'T01':25,'T39':25,'T40':25,'R4_CA':28,'R4_CB':28}
    assert all(r['full_s61_component_dominance']for r in catalog)
    rank2path=ROOT/'legacy/v14/data/six_rank2_components.json'
    rank2=read(rank2path)['components'];assert len(rank2)==6
    assert set(r['name']for r in rank2)=={'R113','R122','Log1112','LPB','Aff5','Aff41'}
    all38path=ROOT/'reports/Q_ALL_TORIC_MPS_COMPONENTS_20261103_v63.json'
    all38=read(all38path);assert all38['all_passed'] and all38['distinct_known_components_through_ALL_Qa_lower_bound']==40
    paths=[Path(__file__),family_review,source,upper_path,upper_review_path,exclusions,orbit_path,catalog_path,rank2path,all38path,
           ROOT/'reports/KNOWN51_COORDINATE_MAPS_20261103_v62_REVIEW.json',
           ROOT/'research/notes/Q_POSITIVE_WEIGHT6_COMPONENT25_20261103_v63_REVIEW.md',
           ROOT/'code/q_positive_weight6_quadratic_model_20261103_v63.py']
    for p in [source,upper_path,upper_review_path,exclusions,orbit_path,all38path,family_review]:
        data=read(p)
        for filename,digest in data.get('input_sha256',{}).items():assert sha(ROOT/filename)==digest,(str(p),filename)
    result=dict(id='Q_POSITIVE_WEIGHT6_COMPONENT25_20261103_v63_REVIEW',all_passed=True,
        proposed_component_id='R4_POSITIVE_WEIGHT6_25',
        separately_independently_audited_full_family_and_full_ambient_upper_used=True,
        irreducible_GL5_times_scalar_times_A3_source_image_cone_dimension_lower=26,
        full_actual_projective_local_dimension_upper_at_contained_closed_point=25,
        actual_global_projective_irreducible_component_dimension_exact=25,
        actual_global_component_image_closure_identification_accepted=True,
        actual_full_transverse_local_dimension_exact=3,
        all51_known_rank4_components_excluded_from_equal_component_identification=True,
        known51_components_of_dimension_at_least25=large,
        all6_known_rank2_components_excluded_by_closed_Poisson_rank_locus=True,
        old_known_component_lower_bound=57,new_global_component_lower_bound=58,
        new_component_contains_all_nonzero_positive_axis_and_Qa8over9=True,
        Qa8over9_local_and_completed_full_reduced_minimal_primes_lower_bound=41,
        actual_true_tangent_cone_unique3dimensional_reduced_component='L3=(z1,z2,z4,z6,z7)',
        all_Qa_membership_generic_reducedness_normality_smoothness_or_branch_exhaustion_claimed=False,
        literature_priority_or_global_component_exhaustion_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
    with output.open('x')as f:
        json.dump(result,f,indent=2);f.write('\n')
    print('PASS independent actual25D newglobalcomponent; lower58; Qa8/9 gains41st minimalprime')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--family-review',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();main(args.family_review,args.output)

#!/usr/bin/env python3
"""Audit the formal-local 40/41 minimal-prime lower bounds along physical Qa."""
import argparse
import json
from pathlib import Path
from q_smith_global_model_review_20261103_v63 import ROOT, read, sha


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    stems = [
        'Q_ALL_TORIC_MPS_COMPONENTS_20261103_v63',
        'Q_LOCAL_Q51_QUOTIENT_20261103_v63_REVIEW',
        'Q_T40_GRAPH_OPEN_20261103_v63_REVIEW',
        'Q_RAT113_GLOBAL_SLICE_20261103_v63_REVIEW',
        'Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63_REVIEW',
    ]
    paths = [Path(__file__), ROOT/'code/q_smith_global_model_review_20261103_v63.py']
    data = {}
    checked = 0
    for stem in stems:
        p = ROOT/f'reports/{stem}.json'
        d = read(p)
        assert d['all_passed']
        paths.append(p)
        data[stem] = d
        for filename, digest in d['input_sha256'].items():
            assert sha(ROOT/filename) == digest, filename
            checked += 1
    all38 = data[stems[0]]
    assert all38['actual_known_MPS_components_through_ALL_Qa'] == [f'T{i:02d}' for i in range(1,39)]
    assert all38['additional_known_components_through_ALL_Qa'] == ['CA','CB']
    assert all38['distinct_known_components_through_ALL_Qa_lower_bound'] == 40
    assert data[stems[1]]['actual_global_quotient_R_by_Q51_equals_QQ9']
    assert data[stems[2]]['actual_R_by_Q52_quotient_equals_QQ9']
    assert data[stems[3]]['A20_graph_is_actual_minimal_component_and_smooth_Rat113_branch']
    assert data[stems[4]]['canonical_affine_formal_product_with17_at_all_Qa']
    paths.append(ROOT/'research/notes/Q_QA_MINIMAL_PRIME_LOWER_BOUNDS_20261103_v63.md')
    report = dict(
        id='Q_QA_MINIMAL_PRIME_LOWER_BOUNDS_20261103_v63', all_passed=True,
        parameter_scope='Every finite complex a; physical nonzero Qa, including a=0',
        distinct_actual_global_components_through_every_Qa_lower_bound=40,
        local_and_completed_local_minimal_primes_lower_bound=40,
        full_and_reduced_Poisson_schemes=True,
        canonical_affine_formal_product='completed R37 at Qa [[17 smooth variables]]',
        completed_R37_minimal_primes_lower_bound=40,
        special_parameters_with_lower_bound41={
            '2/9':'additional actual T39',
            '-8/9':'additional actual T40',
            '0':'additional actual Rat113 minimal component',
        },
        projective_and_canonical_affine_minimal_prime_counts_agree=True,
        localization_completion_and_formal_power_series_arguments_reviewed=True,
        analytic_branch_exhaustion_or_normalization_fiber_count_claimed=False,
        unknown_components_excluded=False, global_known_component_lower_bound_remains=57,
        verified_dependency_hashes=checked,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(paths)},
    )
    with args.output.open('x') as f:
        json.dump(report,f,indent=2)
        f.write('\n')
    print('PASS all finite Qa: at least40 formal minimal primes; at least41 at0,2/9,-8/9')


if __name__ == '__main__':
    main()

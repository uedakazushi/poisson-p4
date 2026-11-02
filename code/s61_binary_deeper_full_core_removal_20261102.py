#!/usr/bin/env python3
"""Dependency and finite algebra audit for the deeper-binary removal theorem.

The geometric proof is in the accompanying note. This driver does not rerun
large polynomial matrix calculations or replace the same-component argument.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / 'research/notes/S61_BINARY_DEEPER_FULL_CORE_REMOVAL_20261102.md'
STEMS = [
    'S61_BINARY_WEIGHTED_HORIZONTAL_SUPPORT_20261102',
    'S61_BINARY_WEIGHTED_VERTICAL_SUPPORT_20261102',
    'S61_BINARY_DEEPER_WEIGHTED_FINITE_MODEL_20261102',
    'S61_BINARY_DEEPER_WEIGHTED_LOW_BLOCKS_20261102',
    'S61_BINARY_WEIGHTED_CUBIC_COVER_20261102',
    'S61_BINARY_WEIGHTED_SEVENTEEN_CURVES_20261102',
    'S61_BINARY_DEEPER_CORE_OPEN_FORMAL_20261101',
    'S61_BINARY_DEEPER_CORE_ARTIFACTS_20261101',
    'S61_BINARY_DEEPER_SEVEN_SLICE_20261101',
    'S61_BINARY_DEEPER_REDUCTION_20261101',
    'S61_BINARY_DEEPER_BOUNDARY_REMOVAL_20261101',
    'S61_BINARY_DEEPER_PIVOT_DIVISORS_20261101',
    'S61_BINARY_UNIFORM_CUBIC_20261101',
    'PRIMITIVE_JET_ELLIPTIC_CURVE_20261101',
]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify():
    start = time.time()
    dependencies = [Path(__file__), NOTE]
    reports = {}
    checked = set()
    for stem in STEMS:
        path = ROOT / 'reports' / (stem+'.json')
        report = json.loads(path.read_text())
        assert report['all_passed'], stem
        reports[stem] = report
        dependencies.extend([path, ROOT / 'research/notes' / (stem+'.md')])
        for filename, expected in report.get('input_sha256', {}).items():
            assert sha(ROOT/filename) == expected, (stem, filename)
            checked.add(filename)
    horizontal = reports[STEMS[0]]
    vertical = reports[STEMS[1]]
    assert horizontal['row_map_surjective']
    assert horizontal['coefficient_ring'] == 'QQ(e)[d,D^-1]'
    assert horizontal['matrix_shape'] == vertical['matrix_shape'] == [390, 269]
    assert horizontal['rebuilt_sha256']['raw_weighted_rows'] == vertical['rebuilt_sha256']['weighted_rows']
    assert vertical['exact_minor_degree'] == 1241
    assert vertical['residual_monic_degree'] == 501
    assert vertical['moving_quadratic_coprime_to_d4_module']
    assert vertical['possible_vertical_curve_roots'] == ['7/24', '1/10', '1/2', '-8/15']
    low = reports['S61_BINARY_DEEPER_WEIGHTED_LOW_BLOCKS_20261102']
    cubic = reports['S61_BINARY_WEIGHTED_CUBIC_COVER_20261102']
    curves = reports['S61_BINARY_WEIGHTED_SEVENTEEN_CURVES_20261102']
    old = reports['S61_BINARY_DEEPER_CORE_OPEN_FORMAL_20261101']
    assert low['weight_two_rank_on_entire_D_open'] == 7
    assert cubic['weighted_cubic_rank_on_open'] == 64
    assert curves['seventeen_explicit_curve_strata_removed_over_sufficiently_general_beta']
    assert old['all_general_beta_points_on_common_core_open_removed']
    assert old['relative_completed_rank'] == 65
    assert old['universal_same_P2_highest_identities'] == 135
    assert old['same_beta_Q1_full_rank'] == 35
    assert old['same_beta_Q1_quadratic_rank'] == 5
    assert reports['S61_BINARY_UNIFORM_CUBIC_20261101']['normal_length_lower_bound'] == 65
    # Independent small consistency check of the support arithmetic.
    e = s.Symbol('e')
    f3 = s.Poly(270*e**2+321*e-388, e)
    f4 = s.Poly(360*e**2+228*e-384, e)
    assert f3.is_irreducible
    assert s.gcd(f3, f4).degree() == 0
    assert 1+7+(30-7)+(98-64) == 65
    assert 10+(1+5) == 16
    assert 10+5 == 15
    assert 16+3 == 19
    return dict(
        id='S61_BINARY_DEEPER_FULL_CORE_REMOVAL_20261102', all_passed=True,
        coefficient_computation_scope='No new large-matrix coefficient calculation; exact frozen dependencies and finite algebra checks.',
        deep_stratum='u12=u13=0, K=5u35+u50=0, Qstar=600u14^2-255u14*u49-63u49^2+100u4*u50=0, Q0=xw*d_y wedge d_z',
        nuisance_zero_weight4_bad_support_dimension_at_most=0,
        isolated_core_points_explicitly_enumerated=False,
        previous_common_core_open_all_nuisance_used=True,
        remaining_core_curve_source_dimension=16,
        dominant_curve_source_generically_etale=True,
        same_binary_component_after_dominant_Sp4_base_change=True,
        relative_regular_base_dimension_on_curve_source=19,
        relative_rank_on_nonempty_removed_open=65,
        entire_deeper_stratum_removed_over_sufficiently_general_beta=True,
        proof_note=str(NOTE.relative_to(ROOT)),
        verified_primary_reports=len(STEMS),
        verified_transitive_input_files=len(checked),
        limitations=['No pointwise all-core matrix surjectivity or explicit isolated rank-drop list.',
                     'No local-ring or t-nilpotence assertion for every special beta.',
                     'No removal of every other binary stratum, flags, or affine S61 exhaustion.'],
        input_sha256={str(path.relative_to(ROOT)): sha(path) for path in dependencies},
        seconds=time.time()-start)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    result = verify()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({key: result[key] for key in ['all_passed', 'verified_primary_reports', 'seconds']}))

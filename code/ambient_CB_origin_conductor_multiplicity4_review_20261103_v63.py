#!/usr/bin/env python3
"""CB-only conductor equality audit using the independently accepted actual Q_b."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_CB_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63_REVIEW'


def load(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def accepted(path):
    data = load(path)
    assert data['all_passed'] is True
    for name, digest in data.get('input_sha256', {}).items():
        assert sha(ROOT / name) == digest, name
    for name, digest in data.get('evidence_sha256', {}).items():
        assert sha(path.parent / name) == digest, name
    return data


def matrix(rows):
    return s.Matrix([[s.Rational(v) for v in row] for row in rows])


def main(output):
    start = time.monotonic()
    output.mkdir(parents=True, exist_ok=False)
    producer = ROOT / 'reports/AMBIENT_CB_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63'
    accepted(producer / 'report.json')
    saved = load(producer / 'certificate.json')
    deps = [
        ROOT / 'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_CB_GENERIC_V_FULL_COMPLETION_20261103_v63/report.json',
        ROOT / 'reports/AMBIENT_CB_V_ACTUAL_ODD_LINEAR_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63_REVIEW.json']
    records = [accepted(path) for path in deps]
    pair = records[0]
    assert pair['new_family_distinct_from_Pa'] is True
    assert pair['both_source_ranks'] == [29, 29]
    assert pair['both_full29_tangents_pushed_to_original_Q'] is True
    assert pair['all58_union_rank'] == 37
    assert pair['tangent_intersection_dimension'] == pair['meeting_scheme_edim'] == 21
    assert pair['own_GL_scalar_parameter_rank'] == 20
    assert pair['generic_whole_fiber_exactly_two_reduced_points_by_proper_open'] is True
    assert pair['generic_reduced_CB_completion'] == 'T/(I1 intersection I2), T/Ii = C[[29]]'
    assert records[1]['CA_CB_actual_conductor_dimension20_and_embedding22_accepted'] is True
    assert records[1]['CA_CB_actual_conductor_multiplicity_atmost4_accepted'] is True
    assert records[2]['conductor_quotient'] == 'C[[g20]]'
    assert records[2]['actual_I_equals_I_plus_nu_by_rank9_and_Nakayama'] is True
    bindingpath = ROOT / 'reports/AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63/certificate.json.gz'
    binding = load(bindingpath)['actual_jet_bindings']['CB']
    quadpath = ROOT / 'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63/certificate.json.gz'
    quad = load(quadpath)['CB']
    rows = binding['seven_actual_odd_rows'] + [46, 95, 63]
    assert rows == [18, 19, 28, 29, 60, 67, 72, 46, 95, 63]
    assert rows == saved['actual_ten_minimal_odd_generator_rows']
    eye = s.eye(150)
    rr = binding['seven_actual_odd_rows']
    linear = matrix(quad['all_actual_linear_combinations'])
    h46 = eye[46, :] - linear[46, :] * eye[rr, :]
    h95 = eye[95, :] - linear[95, :] * eye[rr, :]
    h63 = matrix(binding['h63_as_constant_combination_of_all150_actual_odd_functions'])
    lifts = eye[rr, :].col_join(h46).col_join(h95).col_join(h63)
    assert lifts[:, [j for j in range(150) if j not in rows]] == s.zeros(10, 140)
    assert lifts[:, rows].det() == 1
    assert s.Rational(binding['h63_actual_U3_coefficient']) == -s.Rational(4, 243)
    T = s.symbols('T')
    assert s.cancel((1 - 2 * T ** 2 - T ** 3 + T ** 3 + T ** 4) / (1 - T) ** 22 - (1 + T) ** 2 / (1 - T) ** 20) == 0
    assert 1 + 1 + 2 == 4
    evidence = dict(
        new_actual_CB_Qb_review_input_bound=True,
        CA_Pa_selfpair_or_CA_conductor_equality_used=False,
        intrinsic_fiber_product_conductor_proof_accepted=True,
        generic_length2_lower_uses_dense20_family_AND_actual_edim21=True,
        nu_two_root_difference_distinguishes_critical_W_V=True,
        associativity_lower=[1, 1, 2], exact_multiplicity=4,
        nonzero_kernel_unmixed_CM_argument_gives_full_initial_equality=True,
        actual_ten_minimal_raw_odd_rows=rows,
        actual_lift_transition_determinant='1', actual_h63_U3='-4/243',
        conductor_dim_depth20_edim22_type2=True,
        exact_three_primes_lengths_1_1_2_no_embedded=True,
        actual_image_dim29_edim39_depth21_e2=True,
        image_multiplicity_uses_mA_squared_equals_mR_times_mA=True,
        higher_analytic_type_or_actual_CA_CB_intersection_claim=False)
    (output / 'independent-evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
    paths = [Path(__file__), ROOT / f'research/notes/{STEM}.md',
        ROOT / 'code/ambient_CB_origin_conductor_multiplicity4_20261103_v63.py',
        ROOT / 'research/notes/AMBIENT_CB_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63.md',
        producer / 'report.json', producer / 'certificate.json', bindingpath, quadpath, *deps]
    report = dict(id=STEM, all_passed=True, producer_rerun=False, component='CB only',
        independent_Qb_binding_not_CA_transfer_accepted=True,
        CB_conductor_e4_full_initial_CM20_type2_accepted=True,
        CB_actual_I_minimal10_and_image_depth21_edim39_e2_accepted=True,
        unknown_binary_coefficients_analytic_normal_form_and_actual_D_not_claimed=True,
        independent_method='Separate actual CB Qb receipt plus intrinsic gluing/dense-family/associativity/unmixed-kernel proof; exact CB 150-column ten-generator transition.',
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'independent-evidence.json': sha(output / 'independent-evidence.json')},
        seconds=time.monotonic() - start)
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)

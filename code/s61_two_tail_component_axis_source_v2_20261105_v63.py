#!/usr/bin/env python3
"""Prepare, without executing, two projections of the saved D12 residue."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main(output):
    output.mkdir(parents=True, exist_ok=False)
    pins = {}

    def load(path):
        p = ROOT / path
        pins[path] = sha(p)
        return json.loads(p.read_text())

    candidate = load('reports/S61_TWO_TAIL_H_EXACT_QQ_ENTRY_20261105_v63/canonical-factor-candidates.json')
    product_peer = load('reports/S61_TWO_TAIL_H_EXACT_QQ_ROOT_DENSE_REVIEW_20261105_v63/report.json')
    assert product_peer['all_passed']
    residue_peer = load('reports/S61_TWO_TAIL_H_D12_ROOT_NEWTON_REVIEW_20261105_v63/report.json')
    assert residue_peer['all_passed'] and residue_peer['H_nondivides_D12_in_BASE6_independently_proved']
    assert residue_peer['D12_remainder'] == [27788, 575, 7510, 2387]
    axis = {}
    for factor in candidate['factors']:
        polynomial = [0] * (factor['degrees'][0] + 1)
        for (a, e, d), value in factor['terms']:
            if e == d == 0:
                polynomial[a] += int(value)
        axis['R7' if factor['total_degree'] == 7 else 'S4'] = polynomial
    assert axis == {'R7': [-6408, 63000, -5605000, 24000000, 0], 'S4': [160, -4000, 0]}
    axis['R7'].pop()
    axis['S4'].pop()
    # A coefficient-table check of the existing exact factor identity, not a D12 projection.
    h_axis = [0] * 5
    for i, c in enumerate(axis['R7']):
        for j, d in enumerate(axis['S4']):
            h_axis[i + j] += c * d
    assert h_axis == [-1025280, 35712000, -1148800000, 26260000000, -96000000000]
    assert all(polynomial[-1] % 32003 and polynomial[0] % 32003 for polynomial in axis.values())
    source = dict(
        prepared_source_only=True, component_values_NOT_run=True,
        same_fixed_axis='e=d=0', same_prime=32003,
        original_H_axis=h_axis, factor_axis_integer_coefficients=axis,
        saved_D12_remainder=residue_peer['D12_remainder'],
        runtime_algorithm=[
            'For each factor, invert only its nonzero constant leading coefficient modulo32003.',
            'Reduce the saved four D12 coefficients by the resulting monic factor, highest degree first.',
            'Save every subtraction, quotient and full product-plus-remainder identity.',
            'Do not treat either quotient ring as a field or remove a component.'],
        per_component_degree={'S4': 1, 'R7': 3}, total_reduction_steps_upper=4,
        hard_field_products=1000, field_temporary_bits=30, integer_ROM_bits=37,
        hard_numeric_and_metadata_owned=4096, new_RAW_cap=1048576,
        soft_seconds=45, hard_seconds=60, RAM=2000000000,
        no_new_D12_source_evaluation=True, no_new_point_prime_GCD_or_factor=True,
        source_denominator_scope='Only inherited a/B/LCD denominators; all remain units on these factor quotients because they are quotients of accepted H_axis ring.',
        positive_implication='A nonzero component projection proves that component factor does not divide D12 in the complete six-factor BASE6 ring.',
        zero_implication='Only this fixed image is inconclusive; it is not global divisibility.',
        geometric_scope='After each component is geometrically irreducible and coprime to U6, nonzero gives a nonempty dimension2 open of that selected-tail component.',
        no_actual_bad_Fitting_or_full57_conclusion=True)
    e = output / 'component-projection-plan.json'
    e.write_text(json.dumps(source, indent=2) + '\n')
    for path in [str(Path(__file__).relative_to(ROOT)),
                 'research/notes/S61_TWO_TAIL_COMPONENT_AXIS_SOURCE_20261105_v63.md']:
        pins[path] = sha(ROOT / path)
    report = dict(id=output.name, all_passed=True, frozen=True, source_only=True,
                  component_projection_values_NOT_run=True,
                  all_runtime_inputs_are_accepted_saved_tables=True,
                  input_sha256=pins, evidence_sha256={e.name: sha(e)})
    p = output / 'report.json'
    p.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(report_sha256=sha(p), component_values_run=False,
                          axis_factor_degrees=[1, 3]), indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output.resolve())

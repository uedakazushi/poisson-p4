#!/usr/bin/env python3
"""Check completed families are regular at the frozen dominance witnesses."""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def verify(output):
    folder = ROOT/'reports/MPS_SMOOTHING_COMPLETION_20261103_v62'
    top = folder/'report.json'
    index = json.loads(top.read_text())
    x = s.symbols('x0:5')
    records = []
    inputs = [Path(__file__), top]
    for record in index['records']:
        cid = record['component_id']
        completed_path = folder/f'{cid}.json'
        dominance_path = ROOT/'reports/MPS_S61_DOMINANCE_20261103_v62'/f'{cid}.json'
        raw_path = ROOT/'reports/MPS_S61_RAW_CATALOGUE_20261103_v62'/f'{cid}.json'
        completed = json.loads(completed_path.read_text())
        dominance = json.loads(dominance_path.read_text())
        raw = json.loads(raw_path.read_text())
        periods = list(map(s.Symbol, completed['period_parameters']))
        smoothings = list(map(s.Symbol, completed['smoothing_parameters']))
        period_sample = dict(zip(periods, map(s.Rational, dominance['period_sample'])))
        sample = period_sample | dict.fromkeys(smoothings, 0)
        pfaffian = s.sympify(completed['period_pfaffian_open']).subs(sample)
        assert pfaffian
        denominators = set()
        coefficient_count = 0
        for expression in completed['homogeneous_bivector'].values():
            polynomial = s.Poly(s.sympify(expression), *x)
            for coefficient in polynomial.coeffs():
                denominator = s.fraction(s.cancel(coefficient))[1]
                value = denominator.subs(sample)
                assert value and not value.free_symbols, (cid, denominator, value)
                denominators.add((str(s.factor(denominator)), str(value)))
                coefficient_count += 1
        # Check the full completed and raw families have the same value and
        # first derivative in every source direction at the exact rank point.
        for pair in set(completed['homogeneous_bivector']) | set(raw['homogeneous_bivector']):
            difference = s.cancel(s.sympify(completed['homogeneous_bivector'].get(pair, '0'))
                                  - s.sympify(raw['homogeneous_bivector'].get(pair, '0')))
            assert s.expand(difference.subs(sample)) == 0
            for variable in periods + smoothings:
                assert s.expand(s.diff(difference, variable).subs(sample)) == 0, (cid, pair, variable)
        records.append(dict(component_id=cid, period_sample=dominance['period_sample'],
                            period_pfaffian_at_sample=str(pfaffian), coefficient_count=coefficient_count,
                            denominator_values=sorted(denominators), all_coefficients_regular=True,
                            completed_value_and_all_first_derivatives_equal_raw=True))
        inputs.extend([completed_path, dominance_path, raw_path])
        print('PASS', cid, len(denominators), 'denominators', flush=True)
    result = dict(id='MPS_CORRECTION_DENOMINATORS_20261103_v62_REVIEW', all_passed=True,
                  completed_family_count=len(records), records=records,
                  consequence='All frozen raw dominance differential matrices are also differentials of their completed rational Poisson families at the same points.',
                  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})
    with output.open('x') as stream:
        json.dump(result, stream, indent=2);stream.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    verify(args.output)

#!/usr/bin/env python3
"""Direct truncated-series and split-prime audit of finite marked-square extras."""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as s

from s61_model import S61

ROOT = Path(__file__).resolve().parents[1]
PRIME = 32003


def verify(output, replay_path=None):
    source = ROOT/'reports/S61_MARKED_SIGMA_FINITE_EXTRAS_20261103_v62/report.json'
    report = json.loads(source.read_text())
    if replay_path is not None:
        replay = json.loads(replay_path.read_text())
        assert len(replay['cases']) == len(report['cases'])
        for old, new in zip(report['cases'], replay['cases']):
            for key, value in old.items():
                if key != 'seconds':
                    assert new[key] == value, (old['orbit'], key)
    model = S61()
    one = [i for i, w in enumerate(model.reduced['parameter_weights']) if w == 1]
    two = [i for i, w in enumerate(model.reduced['parameter_weights']) if w == 2]
    equations = model.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
    results = []
    for case in report['cases']:
        K = s.QQ.algebraic_field(s.sqrt(s.Rational(case['h_square'])), alias='h')
        parse = lambda v: K.from_sympy(s.sympify(v))
        degree = case['transverse_length_upper']
        series = [[K.zero]*(degree+1) for _ in range(62)]
        for i, value in zip(one, case['beta']):
            series[i][0] = parse(value)
        for i, value in zip(two, case['q']):
            series[i][0] = parse(value)
        for i, value in zip(case['variables'], case['normal_direction']):
            series[i][1] = parse(value)
        for record in case['jets']:
            for i, value in zip(case['variables'], record['implicit_lift']):
                series[i][record['degree']] = parse(value)
        def multiply(a, b):
            return [sum((a[i]*b[n-i] for i in range(n+1)), K.zero) for n in range(degree+1)]
        residues = []
        for row in equations:
            residual = [K.zero]*(degree+1)
            for mon, n, d in row:
                term = [K.convert(n)/K.convert(d)]+[K.zero]*degree
                for i in mon:
                    term = multiply(term, series[i])
                for _ in range(2-sum(i in two for i in mon)):
                    term = multiply(term, series[61])
                residual = [a+b for a, b in zip(residual, term)]
            residues.append(residual)
        assert all(not values[0] and not values[1] for values in residues)
        for record in case['jets']:
            d = record['degree']
            expected = list(map(parse, record['residuals']))
            assert [values[d] for values in residues] == expected
            if d < degree:
                assert not any(expected)
        roots = case['residue_h_roots']
        numerator, denominator = s.fraction(s.Rational(case['h_square']))
        for root in roots:
            assert (root*root*int(denominator)-int(numerator)) % PRIME == 0
        def residue(value, root):
            total = 0
            for c in value.to_list():
                assert int(c.denominator) % PRIME
                total = (total*root+int(c.numerator)*pow(int(c.denominator), -1, PRIME)) % PRIME
            return total
        # All reported point, frame, H and implicit coefficients are integral;
        # five decisive rank/obstruction witnesses are units at BOTH primes.
        integrals = list(map(parse, case['beta']+case['q']))
        for key in ['symplectic_frame', 'highest_kernel']:
            integrals.extend(parse(v) for row in case[key] for v in row)
        integrals.extend(parse(v) for v in case['normal_direction'])
        for record in case['jets']:
            integrals.extend(parse(v) for v in record['implicit_lift'])
        for value in integrals:
            for root in roots:
                residue(value, root)
        witnesses = [case[key]['determinant'] for key in ['kernel_rank_minor', 'orbit_minor', 'beta_dominance_minor']]
        witnesses.extend([case['Jacobian_determinant'], case['unit_residual']])
        witness_residues = [[residue(parse(value), root) for root in roots] for value in witnesses]
        assert all(all(row) for row in witness_residues)
        results.append(dict(orbit=case['orbit'], all135_direct_series_match=True,
                            transverse_length_upper=degree, lower_residual_degrees_zero=True,
                            all_point_frame_H_implicit_coefficients_integral=True,
                            roots=roots, five_witness_residues=witness_residues,
                            producer_replay_all_exact_fields_match=replay_path is not None))
        print('PASS', case['orbit'], 'series degree', degree, 'roots', roots, flush=True)
    for p, digest in report['input_sha256'].items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == digest
    assert report['prepared_model_input_sha256'] == 'bdaf728347aad67858dd6090d4cc48646141c3d05ae710bf20e592e88e4c49d7'
    paths = [Path(__file__), source, ROOT/'code/s61_marked_sigma_finite_extras_20261103_v62.py',
             ROOT/'research/notes/S61_MARKED_SIGMA_FINITE_EXTRAS_20261103_v62.md']
    result = dict(id='S61_MARKED_SIGMA_FINITE_EXTRAS_20261103_v62_REVIEW', all_passed=True,
                  cases=results, same_incidence_ambient_free_rank_argument_accepted=True,
                  known_T4_T7_generic_multiplicity_and_t_zero_retained_as_dependencies=True,
                  full_marked_pencil_exhaustion_asserted=False,
                  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    with output.open('x') as stream:
        json.dump(result, stream, indent=2);stream.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--replay-report', type=Path)
    args = parser.parse_args()
    verify(args.output, args.replay_report)

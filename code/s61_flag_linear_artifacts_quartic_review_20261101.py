#!/usr/bin/env python3
"""Independent saved low-matrix checks and three-core coverage proof audit."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def verify():
    start = time.time()
    primary = 'reports/S61_FLAG_LINEAR_ARTIFACTS_QUARTIC_20261101.json'
    data = json.loads((ROOT/primary).read_text())
    assert data['all_passed']
    assert all(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
               for path, digest in data['input_sha256'].items())
    h = s.Symbol('b8')
    v = s.symbols('v0:9')
    expected = [
        [v[1], h*v[2]/4+v[3], -h*v[2]/8+v[6]],
        [10*h*v[1]/11+v[2], 5*h*v[1]/11+v[3], -5*h*v[1]/22+v[6]],
        [80*h*v[1]/159+v[2], 40*h*v[1]/159+v[3], -20*h*v[1]/159+v[6]],
    ]
    for index, case in enumerate(data['cases']):
        assert case['core'] == ['-6', '-12/5', '-119/20'][index]
        q = case['lower_basis']
        A, mu, z, K = [s.Matrix(q[key]).applyfunc(s.Rational)
                        for key in ['basis', 'mu', 'socle', 'right_inverse']]
        assert A*mu == s.zeros(9, 1) and z*mu == s.ones(1, 1)
        assert A*K == s.eye(9) and K*A == s.eye(10)-mu*z
        linear = list(map(s.sympify, q['canonical_cubic_linear']))
        assert all(s.expand(a-b) == 0 for a, b in zip(linear, expected[index]))
        # The three cubic linear factors remain independent at every h.
        pc = q['cubic_pivot_radical']
        rad = [0, 1, 2, 3, 6]
        coefficients = s.Matrix([[s.diff(f, v[rad[j]]) for j in pc] for f in linear])
        assert coefficients == s.eye(3)
        substitutions = {v[int(key)]: s.sympify(coefficient)*v[target]
                         for key, (target, coefficient) in q['radical_substitutions'].items()}
        assert all(s.expand(f.subs(substitutions, simultaneous=True)) == 0 for f in linear)
        cover = case['core_cover']
        assert cover['nonzero_factors'] == ['r']
        powers = [4, 5] if index == 0 else [6, 7]
        assert [entry['gcd_factors'] for entry in cover['minor_covers']] == [[['r', k]] for k in powers]
        assert all(len(entry['minor_choices']) == 8 for entry in cover['minor_covers'])
        witness = case['zero8']
        assert witness['all_candidate_denominators_constant']
        assert [entry['rank'] for entry in witness['checks']] == [12, 13]
        assert all(s.Rational(entry['determinant']) != 0 for entry in witness['checks'])
        assert witness['source_dimension_with_Sp4'] == 6+10 == 16
        assert witness['proper_closed_failure_beta_image_bound'] == 5+10 == 15
    paths = list(data['input_sha256'])+[primary,
        'research/notes/S61_FLAG_LINEAR_ARTIFACTS_QUARTIC_20261101.md',
        str(Path(__file__).relative_to(ROOT)),
        'research/notes/S61_FLAG_LINEAR_ARTIFACTS_QUARTIC_REVIEW_20261101.md']
    return dict(id='S61_FLAG_LINEAR_ARTIFACTS_QUARTIC_REVIEW_20261101', all_passed=True,
                independent_saved_quadratic_split=True, independent_cubic_normal_forms_and_substitutions=True,
                proof_and_source_review_accepted=True, all_primary_input_hashes_verified=True,
                accepted_cores=['-6', '-12/5', '-119/20'], hilbert_function=[1, 9, 36, 97, 201], t_order=4,
                independent_raw4_or_minor_recalculation=False,
                accepted_scope='All points of the three normalized cores over sufficiently general beta.',
                limitations=['The h=0 minors are checked as saved nonzero exact values, not independently regenerated.',
                             'The h!=0 polynomial gcd cover and raw four-jets use the source-audited primary producer.',
                             'The rank-drop core -4/25, full initial equality, ambient nilpotence and saturation remain separate.'],
                input_sha256={path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},
                seconds=time.time()-start)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print('THREE LINEAR FLAG CORE FOURTH REVIEW PASS', result['seconds'])

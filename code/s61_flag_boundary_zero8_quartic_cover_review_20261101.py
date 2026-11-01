#!/usr/bin/env python3
"""Proof/source audit and saved lower-data checks for the h=0 fourth cover."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_FLAG_BOUNDARY_ZERO8_QUARTIC_COVER_20261101'


def verify():
    start = time.time()
    primary = 'reports/'+STEM+'.json'
    global_path = 'reports/S61_FLAG_BOUNDARY_GLOBAL_QUARTIC_COVER_20261101.json'
    global_review = 'reports/S61_FLAG_BOUNDARY_GLOBAL_QUARTIC_COVER_REVIEW_20261101.json'
    data = json.loads((ROOT/primary).read_text())
    previous = json.loads((ROOT/global_path).read_text())
    reviewed = json.loads((ROOT/global_review).read_text())
    assert data['all_passed'] and reviewed['all_passed']
    for record in [data, reviewed]:
        assert all(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
                   for path, digest in record['input_sha256'].items())
    # The reused quadratic and active-cubic matrices are precisely those whose
    # identities were independently checked in the separate global review.
    for name in ['quadratic', 'active_cubic']:
        a = {k: v for k, v in data['lower_basis_certificates'][name].items() if k != 'seconds'}
        b = {k: v for k, v in previous['lower_basis_certificates'][name].items() if k != 'seconds'}
        assert a == b
    r, h = s.symbols('b2 b8')
    cubic = data['lower_basis_certificates']['true_cubic']
    allowed = {r+6, r+s.Rational(12, 5), r+s.Rational(119, 20),
               r+s.Rational(4, 25), r+s.Rational(91, 400), r-2}
    factors = set()
    for row in cubic['left_inverse']:
        for expression in row:
            value = s.cancel(s.sympify(expression))
            assert not value.has(h)
            for factor, exponent in s.factor_list(s.denom(value), r)[1]:
                monic = s.Poly(factor, r).monic().as_expr()
                assert monic in allowed
                factors.add(str(monic))
    denominator = (r-2)*(r+6)**2*(r+s.Rational(4, 25))*(r+s.Rational(91, 400))
    assert s.expand(s.sympify(cubic['primitive_denominator'])-denominator) == 0
    assert data['nonzero_factors'] == ['r + 119/20', 'r + 12/5', 'r + 4/25',
                                      'r + 6', 'r + 91/400', 'r - 2']
    assert data['hilbert_function'] == [1, 9, 36, 97, 201]
    assert (data['quartic_rank'], data['augmented_rank'], data['t_order']) == (12, 13, 4)
    assert data['pure_core_rows_without_t'] == 19 and data['pure_core_columns'] == 40
    expected = [
        {'r + 119/20': 3, 'r - 2': 10, 'r + 91/400': 10, 'r + 4/25': 26, 'r + 6': 29},
        {'r + 119/20': 3, 'r - 2': 11, 'r + 91/400': 11, 'r + 4/25': 29, 'r + 6': 32},
    ]
    for cover, gcd in zip(data['minor_covers'], expected):
        assert dict(cover['gcd_factors']) == gcd
        assert len(cover['minor_choices']) == 8
    paths = list(data['input_sha256'])+[primary, global_path, global_review,
        'research/notes/'+STEM+'.md',
        'research/notes/S61_FLAG_ZERO12_LOW_JETS_20261101.md',
        'reports/S61_FLAG_ZERO12_LOW_JETS_20261101.json',
        str(Path(__file__).relative_to(ROOT)),
        'research/notes/S61_FLAG_BOUNDARY_ZERO8_QUARTIC_COVER_REVIEW_20261101.md']
    return dict(id='S61_FLAG_BOUNDARY_ZERO8_QUARTIC_COVER_REVIEW_20261101',
                all_passed=True, proof_and_source_review_accepted=True,
                exact_reused_low_matrices_equal_independently_reviewed_global_data=True,
                independent_true_cubic_left_inverse_denominator_check=True,
                cubic_left_inverse_has_no_h=True, cubic_denominator_factors=sorted(factors),
                accepted_hilbert_function=[1, 9, 36, 97, 201], t_order=4,
                accepted_scope='h=0, (r-2)*D*rho nonzero, all six nuisance values.',
                accepted_combined_scope='Every ordinary D*rho nonzero point over sufficiently general beta.',
                independent_raw4_or_minor_recalculation=False,
                limitations=['LC*CT and fourth minors are source-audited primary certificates.',
                             'No assertion on D=0, rho=0, the zero12 chart, full initial ideal or ambient t-nilpotence.'],
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
    print('ZERO8 FLAG QUARTIC COVER REVIEW PASS', result['seconds'])

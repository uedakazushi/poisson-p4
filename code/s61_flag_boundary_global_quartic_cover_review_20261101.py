#!/usr/bin/env python3
"""Independent lower-basis identities and proof/source audit of the global cover."""
import argparse
import hashlib
import json
from itertools import combinations_with_replacement
from pathlib import Path
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_FLAG_BOUNDARY_GLOBAL_QUARTIC_COVER_20261101'


def verify():
    start = time.time()
    report = 'reports/'+STEM+'.json'
    data = json.loads((ROOT/report).read_text())
    assert data['all_passed']
    assert all(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
               for path, digest in data['input_sha256'].items())
    r, h = s.symbols('b2 b8')
    field = s.QQ.frac_field(r)
    lower = data['lower_basis_certificates']
    quad = lower['quadratic']

    def matrix(rows):
        return s.Matrix(rows).applyfunc(s.sympify).to_DM().convert_to(field)

    A, mu, z, K = (matrix(quad[key]) for key in ['basis', 'mu', 'socle', 'right_inverse'])
    assert A.shape == (9, 10) and K.shape == (10, 9)
    assert (A*mu).is_zero_matrix
    assert (z*mu).to_Matrix() == s.ones(1, 1)
    assert (A*K).to_Matrix() == s.eye(9)
    assert (K*A+mu*z).to_Matrix() == s.eye(10)
    assert quad['socle_degree'] == 0

    # Rebuild the multiplication matrix from the saved nine quadratic forms.
    # This does not import or call any primary engine.
    pairs = list(combinations_with_replacement(range(4), 2))
    triples = list(combinations_with_replacement(range(4), 3))
    B = s.zeros(36, 20)
    a = A.to_Matrix()
    for row in range(9):
        for variable in range(4):
            for column, monomial in enumerate(pairs):
                target = triples.index(tuple(sorted((*monomial, variable))))
                B[4*row+variable, target] += a[row, column]
    active = lower['active_cubic']
    X = matrix(active['lift_numerator'])
    g = s.sympify(active['denominator'])
    assert (X*B.to_DM().convert_to(field)).to_Matrix().applyfunc(s.cancel) == g*s.eye(20)
    assert s.expand(g-(r+s.Rational(119, 20))**2*(r+s.Rational(91, 400))) == 0

    allowed = {r+6, r+s.Rational(12, 5), r+s.Rational(119, 20),
               r+s.Rational(4, 25), r+s.Rational(91, 400), h}
    checked_factors = set()

    def check_factors(expression):
        numerator, denominator = s.fraction(s.cancel(expression))
        for polynomial in [numerator, denominator]:
            if polynomial.is_Number:
                continue
            for factor, exponent in s.factor_list(polynomial, r, h)[1]:
                monic = s.Poly(factor, r, h).monic().as_expr()
                assert monic in allowed, monic
                checked_factors.add(str(monic))

    for content in quad['row_contents']:
        check_factors(s.sympify(content))
    check_factors(s.sympify(quad['inverse_denominator']))
    check_factors(g)
    cubic = lower['true_cubic']
    check_factors(s.sympify(cubic['primitive_denominator']))
    # The saved L is independently checked for regularity; its product with CT
    # is verified by the primary source, whose CT matrix is not in this report.
    for row in cubic['left_inverse']:
        for expression in row:
            check_factors(s.denom(s.cancel(s.sympify(expression))))

    assert data['nonzero_factors'] == ['h', 'r + 119/20', 'r + 12/5',
                                      'r + 4/25', 'r + 6', 'r + 91/400']
    assert data['hilbert_function'] == [1, 9, 36, 97, 201]
    assert data['quartic_rank'] == 12 and data['augmented_rank'] == 13
    assert data['t_order'] == 4
    assert s.binomial(8, 4)+4*s.binomial(7, 3)+s.binomial(3, 2) == 213
    assert s.binomial(12, 4)-213+12 == 294
    expected_gcds = [
        {'r + 119/20': 2, 'r + 91/400': 10, 'r + 4/25': 26, 'r + 6': 29},
        {'r + 119/20': 3, 'r + 91/400': 11, 'r + 4/25': 29, 'r + 6': 35},
    ]
    for cover, expected in zip(data['minor_covers'], expected_gcds):
        assert dict(cover['gcd_factors']) == expected
        assert len(cover['minor_choices']) == 8
    paths = list(data['input_sha256'])+[report, 'research/notes/'+STEM+'.md',
        str(Path(__file__).relative_to(ROOT)),
        'research/notes/S61_FLAG_BOUNDARY_GLOBAL_QUARTIC_COVER_REVIEW_20261101.md']
    return dict(
        id='S61_FLAG_BOUNDARY_GLOBAL_QUARTIC_COVER_REVIEW_20261101',
        all_passed=True, primary_input_hashes_verified=True,
        independent_saved_quadratic_split_identities=True,
        independent_saved_active_cubic_lift_identity=True,
        independent_lower_basis_denominator_support=True,
        proof_and_source_review_accepted=True,
        checked_lower_factors=sorted(checked_factors),
        accepted_scope='Every eight-parameter point with b8*D(b2)*(400*b2+91) nonzero.',
        hilbert_function=[1, 9, 36, 97, 201], t_order=4,
        independent_raw4_reimplementation=False,
        independent_true_cubic_product_or_fourth_minor_recalculation=False,
        limitations=[
            'Independent arithmetic is on the saved quadratic split and active cubic lift.',
            'The raw135 jets, LC*CT product and all fourth minors use the audited primary producer.',
            'No full initial ideal, ambient t-nilpotence or flag-stratum removal.'],
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
    print('GLOBAL FLAG QUARTIC COVER REVIEW PASS', result['seconds'])

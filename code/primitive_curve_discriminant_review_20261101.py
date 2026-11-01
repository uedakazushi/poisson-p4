#!/usr/bin/env python3
"""Independent first-jet and discriminant-line checks, plus formal cubic data."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import sympy as s

from weight1_structure import constraints, form_matrices

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'reports/PRIMITIVE_CURVE_DISCRIMINANT_20261101.json'
SLICE = ROOT / 'reports/PRIMITIVE_CURVE_DOMINANT_SLICE_20261101.json'


def multiply(left, right):
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            result[i + j] += x * y
    return result


def line_restriction(encoded, base, direction):
    result = [Fraction(0)] * (max(sum(mon) for mon, _ in encoded) + 1)
    for monomial, coefficient in encoded:
        term = [Fraction(coefficient)]
        for i, exponent in enumerate(monomial):
            if exponent:
                term = multiply(term, [Fraction(math.comb(exponent, j) * base[i] ** (exponent - j)
                                                * direction[i] ** j) for j in range(exponent + 1)])
        for i, value in enumerate(term):
            result[i] += value
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    record = json.loads(SOURCE.read_text())
    data = json.loads(SLICE.read_text())
    basis = s.Matrix.hstack(*s.Matrix.vstack(*constraints()).nullspace())
    e = s.eye(4)
    u, v = e[:, 0], e[:, 2]
    variations = [(e[:, 1], s.zeros(4, 1)), (e[:, 3], e[:, 1]), (s.zeros(4, 1), e[:, 3])]
    columns = []
    for col in range(16):
        tensor = form_matrices(basis[:, col])

        def trilinear(direction, left, right):
            return sum(direction[k] * (left.T * tensor[k] * right)[0] for k in range(4))

        values = [trilinear(u, u, v), trilinear(v, u, v)]
        # Differentiate the three tensor arguments directly. This construction
        # does not form or differentiate the source script's chart polynomials.
        first = [trilinear(du, u, v) + trilinear(u, du, v) + trilinear(u, u, dv)
                 for du, dv in variations]
        second = [trilinear(dv, u, v) + trilinear(v, du, v) + trilinear(v, u, dv)
                  for du, dv in variations]
        columns.append(s.Matrix(values + first + second))
    jet = s.Matrix.hstack(*columns)
    assert jet == s.Matrix(record['first_jet_matrix'])
    reverse = list(reversed(range(16)))
    pivots = [reverse[i] for i in jet[:, reverse].rref()[1]]
    assert len(pivots) == 8
    minor = jet[:, pivots]
    assert minor.det()

    base, direction = [1, 2, 3, 4, 5, 6], [2, -1, 3, 5, 7, 11]
    c4 = line_restriction(data['c4_sparse'], base, direction)
    c6 = line_restriction(data['c6_sparse'], base, direction)
    fourth_cube = multiply(multiply(c4, c4), c4)
    sixth_square = multiply(c6, c6)
    delta = [(a - b) / 1728 for a, b in zip(fourth_cube, sixth_square, strict=True)]
    assert len(delta) == 25 and delta[-1]
    t = s.Symbol('t')
    polynomial = s.Poly(sum(s.Rational(value.numerator, value.denominator) * t ** i
                            for i, value in enumerate(delta)), t)
    assert polynomial == s.Poly(s.sympify(record['delta']), t)
    denominator, integral = polynomial.clear_denoms(convert=True)
    prime = 65521
    assert int(denominator) % prime and int(integral.LC()) % prime
    reduced = s.Poly(integral.as_expr(), t, modulus=prime)
    left, right, gcd = s.gcdex(reduced, reduced.diff())
    assert gcd == s.Poly(1, t, modulus=prime)
    assert left * reduced + right * reduced.diff() == gcd
    assert s.gcd(polynomial, polynomial.diff()).degree() == 0

    # Generic double-root normal form for the normalized cubic cover.
    x, r = s.symbols('X r')
    a, b = 3 * r ** 2, -2 * r ** 3
    assert s.expand(x ** 3 - a * x - b - (x - r) ** 2 * (x + 2 * r)) == 0
    ca, cb = s.symbols('c4 c6')
    cubic = x ** 3 - 27 * ca * x - 54 * cb
    cubic_discriminant = s.discriminant(cubic, x)
    assert s.expand(cubic_discriminant - 78732 * (ca ** 3 - cb ** 2)) == 0
    assert 78732 * 1728 == 136048896

    report = {
        'all_passed': True,
        'jet_check_method': 'direct trilinear differentiation of three arguments',
        'first_jet_matrix_matches': True, 'first_jet_rank': 8,
        'independent_minor_columns': pivots, 'independent_minor_determinant': str(minor.det()),
        'line_check_method': 'Fraction coefficient arrays and binomial convolution',
        'line_degree': 24, 'line_matches_source': True,
        'second_prime': prime, 'second_prime_bezout_left': str(left.as_expr()),
        'second_prime_bezout_right': str(right.as_expr()), 'second_prime_gcd': 1,
        'rational_gcd_degree': 0,
        'cubic_discriminant_equals_136048896_Delta': True,
        'generic_discriminant_factorization': '(X-r)^2*(X+2r), r=-3b/(2a), a=27c4,b=54c6',
        'proof_review': 'Critical-zero incidence is irreducible and proper over its image; '
                        'uniform Koszul cohomology and Fisher-Sadek Theorem2.1(ii) identify '
                        'that image with V(Delta). Square-free nonconstant restriction gives multiplicity1.',
        'normalized_cover_conclusion': 'The monogenic cubic is normal by R1+S2, finite flat of degree3; '
                        'its branch divisor is reduced Delta with generic inertia a transposition.',
        'scope': 'Independent audit of the curve discriminant and the normalized finite cubic cover. '
                 'No identification of its entire boundary with the actual T01 projection.',
        'source_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in (Path(__file__), SOURCE, SLICE, ROOT / 'code/weight1_structure.py')},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps({'all_passed': True, 'jet_rank': 8, 'jet_minor': str(minor.det()),
                      'square_free_second_prime': prime, 'cubic_discriminant_scale': 136048896}))


if __name__ == '__main__':
    main()

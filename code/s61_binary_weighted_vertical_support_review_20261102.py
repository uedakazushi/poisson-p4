#!/usr/bin/env python3
"""Independent small algebra/source review; not a second large-matrix implementation."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / 'reports/S61_BINARY_WEIGHTED_VERTICAL_SUPPORT_20261102.json'
NOTE = ROOT / 'research/notes/S61_BINARY_WEIGHTED_VERTICAL_SUPPORT_REVIEW_20261102.md'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify():
    start = time.time()
    report = json.loads(MAIN.read_text())
    assert report['all_passed']
    for name, expected in report['input_sha256'].items():
        assert digest(ROOT / name) == expected, name
    e, delta = s.symbols('e delta')
    f3 = s.Poly(270*e**2+321*e-388, e)
    f4 = s.Poly(360*e**2+228*e-384, e)
    assert f3.is_irreducible
    assert s.discriminant(f3) == 522081
    assert f3.eval(2) == 1334
    assert f3.eval(-s.Rational(287, 90)) == 1334
    # Compute multiplication by e directly from e^2=(delta+388-321e)/270.
    multiplication = s.Matrix([[0, (delta+388)/270], [1, -s.Rational(321, 270)]])
    assert s.expand(multiplication.det()) == -(delta+388)/270
    assert (270*multiplication**2+321*multiplication-388*s.eye(2)).applyfunc(s.expand) == delta*s.eye(2)
    assert (delta*s.eye(2)).det() == delta**2
    assert s.factor(multiplication.subs(delta, 1334).charpoly(e).as_expr()) == (e-2)*(90*e+287)/90
    roots = [s.Rational(7, 24), s.Rational(1, 10), s.Rational(1, 2), -s.Rational(8, 15)]
    rational_factors = [s.Poly(e-root, e) for root in roots]
    prime = report['prime']
    f3mod = s.Poly(f3.as_expr(), e, modulus=prime)
    f4mod = s.Poly(f4.as_expr(), e, modulus=prime)
    assert s.gcd(f3, f4).degree() == 0
    assert s.gcd(f3mod, f4mod).degree() == 0
    for factor in rational_factors:
        integer_factor = factor.clear_denoms()[1]
        assert s.gcd(f3mod, s.Poly(integer_factor.as_expr(), e, modulus=prime)).degree() == 0
    powers = report['exact_multiplicities']
    assert powers == [70, 113, 176, 197, 92]
    assert report['infinity']['valuation'] == 366
    assert report['tropical_degree_upper'] - 366 == report['exact_minor_degree'] == 1241
    assert 1241 - sum(powers[:4]) - 2*powers[4] == report['residual_monic_degree'] == 501
    assert report['local']['quad']['valuation'] == 184
    assert report['local']['quad']['original_quadratic_multiplicity'] == 92
    for item in [report['infinity'], *report['local'].values()]:
        stages = item['stages']
        assert sum(row[2] for row in stages) == item['valuation']
        assert stages[-1][2] == 0
    h3 = s.Poly(1, e, modulus=prime)
    h4 = s.Poly(1, e, modulus=prime)
    for root, power in zip(roots, [65, 113, 138, 160]):
        root_mod = int(root.p)*pow(int(root.q), -1, prime) % prime
        term = s.Poly(e-root_mod, e, modulus=prime)**power
        h3 *= term
        h4 *= term
    h3 *= f3mod.monic()**92
    h4 *= f4mod.monic()**92
    assert h3.degree() == h4.degree() == 660
    assert s.gcd(h3, h4).degree() == 476
    assert s.gcd(f3mod, h4).degree() == 0
    quad = ROOT / 'reports/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101.json'
    known = json.loads(quad.read_text())['irreducible_minor_factors']
    for root in roots:
        assert any(s.expand(s.sympify(f).subs(e, root)) == 0 for f in known), root
    dependencies = [Path(__file__), NOTE, MAIN, quad]
    dependencies.extend(ROOT/name for name in report['input_sha256'])
    return dict(
        id='S61_BINARY_WEIGHTED_VERTICAL_SUPPORT_REVIEW_20261102', all_passed=True,
        review_scope='Independent small QQ/modular algebra, full source and proof audit; no second implementation of the 269/538 polynomial matrix peeling.',
        primary_input_hashes_checked=len(report['input_sha256']),
        regular_representation_identity_checked=True,
        quadratic_irreducibility_and_norm_factor_two_checked=True,
        rational_roots_are_known_D_divisors=True,
        d3_d4_common_modular_degree=476,
        residual_degree=501,
        finite_DVR_localization_argument_accepted=True,
        primitive_reduction_vertical_exclusion_accepted=True,
        limitations=['Large raw coefficient, kernel, and determinant identities remain dependencies of the main exact certificate.',
                     'No all-core pointwise surjectivity, nuisance extension, or ambient saturation assertion in this review.'],
        input_sha256={str(path.relative_to(ROOT)): digest(path) for path in dependencies},
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
    print(json.dumps({key: result[key] for key in ['all_passed', 'primary_input_hashes_checked', 'seconds']}))

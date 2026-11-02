#!/usr/bin/env python3
"""Exact H^2 coefficient Bezout identity in the B-free raw residuals."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
FRAME = ROOT / 'reports/S61_BINARY_HYPERPLANE_BFREE_FRAME_20261102.json'
ENGINE = ROOT / 'code/parametric_implicit_jets_20261101.cpp'
NOTE = ROOT / 'research/notes/S61_BINARY_HYPERPLANE_BFREE_HEAVY_SQUARE_20261102.md'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_heavy_square(wd, raw, denominator_input=None):
    """Verify the supplied fresh raw jet file; the optional coeff file is audited.

    By default the denominator is the certified frame's U_common_denominator.
    If denominator_input is supplied, its first sparse QQ[a,e,d] polynomial
    must equal this denominator. Return a JSON-serializable certificate.
    """
    start = time.time()
    wd = Path(wd).resolve()
    wd.mkdir(parents=True, exist_ok=True)
    raw = Path(raw).resolve()
    a, e, d = s.symbols('a e d')
    frame = json.loads(FRAME.read_text())
    assert frame['all_passed']
    for name, expected in frame['input_sha256'].items():
        assert sha(ROOT/name) == expected
    denominator = s.Poly(s.sympify(frame['U_common_denominator']), a, e, d, domain=s.QQ)
    assert denominator and denominator.degree(d) == 0
    generated = {'raw': sha(raw)}
    if denominator_input is not None:
        denominator_input = Path(denominator_input).resolve()
        with denominator_input.open() as stream:
            count = int(stream.readline())
            terms = {}
            for _ in range(count):
                z = stream.readline().split()
                assert len(z) == 4
                exponent = tuple(map(int, z[1:]))
                terms[exponent] = terms.get(exponent, s.Rational(0))+s.Rational(z[0])
            assert s.Poly.from_dict(terms, (a, e, d), domain=s.QQ) == denominator
        generated['denominator_input'] = sha(denominator_input)

    # The engine stores 45 W2 polynomials, then 135 A2 polynomials. At degree2
    # its A2 numerator is Uden*H(Ty,Ty)+J*W2, so actual A2 is this/Uden.
    target = (0, 0, 0, 0, 0, 2, 0, 0)
    coefficients = []
    with raw.open() as stream:
        assert list(map(int, stream.readline().split())) == [135, 3, 8, 3]
        for i in range(180):
            count = int(stream.readline())
            terms = {}
            for _ in range(count):
                z = stream.readline().split()
                assert len(z) == 12
                if i < 45 or tuple(map(int, z[4:])) != target:
                    continue
                exponent = tuple(map(int, z[1:4]))
                terms[exponent] = terms.get(exponent, s.Rational(0))+s.Rational(z[0])
            if i >= 45:
                coefficients.append(s.Poly.from_dict(terms, (a, e, d), domain=s.QQ))
    assert len(coefficients) == 135
    nonzero = [i for i, c in enumerate(coefficients) if c]
    assert len(nonzero) == 21
    K = s.QQ.frac_field(a, e)
    polynomials = [s.Poly(c.as_expr(), d, domain=K) for c in coefficients]
    first, second = 77, 131
    u, v, gcd = s.gcdex(polynomials[first], polynomials[second])
    assert gcd == s.Poly(1, d, domain=K)
    assert u*polynomials[first]+v*polynomials[second] == gcd
    assert all(not c.rem(gcd) for c in polynomials)

    # Audit recognizable factorizations against the fresh coefficients;
    # these formulas are not imported from the old B-inverted frame.
    F = 500*a+90*e**2+3*e-24
    B = 50*a*d+18*d*e**2-9*d*e-2*d+60*e-40
    C = 100*a+18*e**2-9*e-2
    D0 = 50*a+9*e**2+3*e-6
    P = 125*a**2*d-450*a*e-200*a-54*e**3-27*e**2+33*e+6
    Q = 15*a*d*e-10*a*d+200*a+18*e**2+6*e-12
    expected = {
        first: -3*a**2*D0*C*F*P**2/s.Integer(3125000000000),
        second: 3*a**2*D0*C*F*B*Q/s.Integer(62500000000),
    }
    for i, value in expected.items():
        assert coefficients[i] == s.Poly(value, a, e, d, domain=s.QQ)
    denexpr = denominator.as_expr()
    actual = {i: s.Poly(polynomials[i].as_expr()/denexpr, d, domain=K)
              for i in [first, second]}
    cofactors = {first: u, second: v}
    actual_cofactors = {i: s.Poly(c.as_expr()*denexpr, d, domain=K)
                       for i, c in cofactors.items()}
    assert sum((actual_cofactors[i]*actual[i] for i in actual), s.Poly(0, d, domain=K)) == gcd
    assert all(c.degree() <= 1 for c in cofactors.values())

    prime = 32003
    allowed = [s.Poly(f, a, e, domain=s.QQ).monic() for f in [a, F, D0, C, 6*e+1]]
    for c in coefficients:
        assert all(int(x.q) % prime for x in c.coeffs())
    for c in list(cofactors.values())+list(actual_cofactors.values()):
        for value in c.all_coeffs():
            num, den = s.fraction(s.cancel(value))
            assert s.degree(den, d) == 0
            for factor, _ in s.factor_list(den, a, e)[1]:
                assert s.Poly(factor, a, e, domain=s.QQ).monic() in allowed
            val = s.cancel(den.subs({a: 1, e: 2}))
            assert val.is_Rational and int(val.p) % prime and int(val.q) % prime
            assert all(int(x.q) % prime for x in s.Poly(num, a, e, domain=s.QQ).coeffs())

    paths = [Path(__file__), FRAME, ENGINE, NOTE]
    result = {
        'all_passed': True,
        'scope': 'Exact QQ(a,e)[d] Bezout=1 for the H^2 coefficients of supplied fresh '
                 'B-free ordinary-degree2 residuals. No all-core open, length65, or ambient claim.',
        'implicit_common_denominator': str(denexpr),
        'coefficient_relation': 'actual_H2_i = raw_H2_i/U_common_denominator',
        'nonzero_sources': nonzero,
        'raw_heavy_square_coefficients': [str(s.factor(c.as_expr())) for c in coefficients],
        'gcd': '1',
        'bezout_source_rows': [first, second],
        'bezout_raw': {str(i): str(s.factor(c.as_expr())) for i, c in cofactors.items()},
        'bezout_actual': {str(i): str(s.factor(c.as_expr())) for i, c in actual_cofactors.items()},
        'bezout_coefficient_degrees_in_d': [u.degree(), v.degree()],
        'raw_cofactor_common_denominator': str(s.factor(s.lcm(
            *[s.fraction(s.factor(c.as_expr()))[1] for c in cofactors.values()]))),
        'prime': prime, 'good_point': [1, 2],
        'checks': {
            'all135_coefficients_extracted': True,
            'exactly21_nonzero_coefficients': True,
            'two_raw_factorizations_from_fresh_coefficients': True,
            'raw_polynomial_Bezout_identity': True,
            'actual_polynomial_Bezout_identity': True,
            'cofactor_denominators_independent_of_d': True,
            'cofactor_denominators_supported_on_a_F_D0_C_6e1': True,
            'cofactor_denominators_good_at_specified_modular_point': True,
        },
        'input_sha256': {str(path.relative_to(ROOT)): sha(path) for path in paths},
        'generated_inputs_sha256': generated,
        'seconds': time.time()-start,
    }
    (wd/'Bfree-heavy-square-certificate.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, required=True)
    parser.add_argument('--denominator-input', type=Path)
    parser.add_argument('--workdir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run_heavy_square(args.workdir, args.raw, args.denominator_input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'all_passed': True, 'gcd': result['gcd'], 'seconds': result['seconds']}))


if __name__ == '__main__':
    main()

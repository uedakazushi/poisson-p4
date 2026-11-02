#!/usr/bin/env python3
"""Exact constant normal form of six supplied d-infinity quadratic rows."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / 'code/certificates/s61_binary_hyperplane_infinity_quadrics_20261102.json'
NOTE = ROOT / 'research/notes/S61_BINARY_HYPERPLANE_INFINITY_QUADRATIC_NORMAL_FORM_20261102.md'


def matrix_strings(matrix):
    return [[str(s.factor(v)) for v in matrix.row(i)] for i in range(matrix.rows)]


def verify():
    begin = time.time()
    data = json.loads(CERTIFICATE.read_text())
    a, e = s.symbols('a e')
    x = s.symbols('x0:8')
    light = [0, 1, 2, 3, 4, 6, 7]
    old_light = [x[i] for i in light]
    ids = [19, 60, 84, 114, 117, 133]
    assert [r['source'] for r in data['rows']] == ids
    quadrics = []
    for row in data['rows']:
        assert row['weight'] == 2
        poly = 0
        for monomial, value in row['terms']:
            assert len(monomial) == 8 and sum(monomial) == 2 and not monomial[5]
            coefficient = s.Poly(s.sympify(value), a, e, domain=s.QQ).as_expr()
            poly += coefficient * s.prod(z**n for z, n in zip(x, monomial))
        quadrics.append(s.expand(poly))

    L = 10*a*x[2] - 15*a*x[3] + (12*e+2)*x[1] + 90*x[4] - 90*x[6] + 30*x[7]
    U = (a*x[0] + 28*a*x[4] + 28*a*x[6] - 28*a*x[7]
         - 6*e*x[6] + 6*e*x[7] + 3*x[4] - 4*x[6] + 2*x[7])
    W = (-15*e*x[2] - 12*e*x[3] + 5*x[0] - 15*x[1]
         + 3*x[3] - 5*x[4] - 5*x[6] + 75*x[7])
    p = 2*x[2] + x[3]
    q = (-40*a*x[2] + 60*a*x[3] - 48*e*x[1] + 45*x[0]
         - 8*x[1] - 270*x[2] - 135*x[3] + 720*x[6] + 420*x[7])
    Z = q/27 + 5*e*p + 2*W/3
    V = 15*U + 2*a*W - 8*a*Z - (s.Rational(32, 27)*a + 3*e)*L
    C = s.expand((-4*L + 18*W - 27*Z)/18)
    P = s.Rational(15, 2)*p
    E = s.expand(3*V/(16*a))
    forms = [L, W, C, P, E]

    y = s.symbols('L W C P E r s')
    canonical = [y[0]**2, y[0]*y[1], y[0]*(y[4]-y[2]),
                 y[1]**2-y[1]*y[2]-y[0]*y[3], y[3]*y[4], y[4]*y[2]]
    generator = s.Matrix([
        [-10, 0, 0, 0, 0, 1],
        [0, 0, 1, 0, 0, 0],
        [(90*e-45)/(16*a), 0, -s.Rational(5, 8), 0, 0, -9*e/(16*a)],
        [0, s.Rational(3, 2), 0, 0, 0, 0],
        [0, -45*e/(32*a), 45/(32*a), 45/(32*a), 0, 0],
        [0, -s.Rational(75, 16), 45*(1-e)/(32*a), 0, 15/(32*a), 0],
    ])
    new_quadrics = [f.subs(dict(zip(y[:5], forms)), simultaneous=True) for f in canonical]
    for left, right in zip(new_quadrics, generator*s.Matrix(quadrics)):
        assert all(s.cancel(c) == 0 for c in s.Poly(s.expand(left-right), *old_light).coeffs())
    generator_inverse = generator.inv().applyfunc(s.factor)
    assert (generator*generator_inverse-s.eye(6)).applyfunc(s.cancel) == s.zeros(6)
    assert (generator_inverse*generator-s.eye(6)).applyfunc(s.cancel) == s.zeros(6)
    generator_det = s.factor(generator.det())
    assert generator_det == s.Rational(91125, 32768)/a**3

    # The complement is fixed, so no data-dependent pivot choice enters the theorem.
    T = s.Matrix([[s.diff(f, z) for z in old_light] for f in forms + [x[1], x[3]]])
    Ti = T.inv().applyfunc(s.factor)
    assert (T*Ti-s.eye(7)).applyfunc(s.cancel) == s.zeros(7)
    assert (Ti*T-s.eye(7)).applyfunc(s.cancel) == s.zeros(7)
    coordinate_det = s.factor(T.det())
    chart_factor = 50*a-6*e-1
    assert s.cancel(coordinate_det-s.Rational(2278125, 4)*chart_factor/a) == 0
    T8 = s.zeros(8)
    for i, output in enumerate(light):
        for j, old in enumerate(light):
            T8[output, old] = T[i, j]
    T8[5, 5] = 1
    T8i = T8.inv().applyfunc(s.factor)
    assert (T8*T8i-s.eye(8)).applyfunc(s.cancel) == s.zeros(8)
    assert (T8i*T8-s.eye(8)).applyfunc(s.cancel) == s.zeros(8)
    y8 = s.symbols('L W C P E H r s')
    assert (T8*s.Matrix(x)-s.Matrix(forms + [x[5], x[1], x[3]])).applyfunc(s.expand) == s.zeros(8, 1)

    # Every matrix is defined on the stated chart, not only in the fraction field.
    allowed = [s.Poly(a, a, e).monic(), s.Poly(chart_factor, a, e).monic()]
    for matrix in [T, Ti, T8, T8i, generator, generator_inverse]:
        for value in matrix:
            denominator = s.fraction(s.cancel(value))[1]
            for factor, _ in s.factor_list(denominator, a, e)[1]:
                assert s.Poly(factor, a, e).monic() in allowed

    paths = [Path(__file__), CERTIFICATE, NOTE]
    return {
        'all_passed': True,
        'scope': 'Exact equivalence of six supplied infinity quadrics on a*(50*a-6*e-1)!=0; '
                 'no original135 rederivation, higher-degree rank, or ambient conclusion.',
        'certificate_provenance': {k: data[k] for k in ['provenance_path', 'provenance_sha256']},
        'coefficient_field': 'QQ(a,e)',
        'localized_chart': 'QQ[a,e,1/(a*(50*a-6*e-1))]',
        'quadric_source_ids': ids,
        'old_light_order': [str(z) for z in old_light],
        'new_light_order': [str(z) for z in y],
        'new8_order': [str(z) for z in y8],
        'normal_generators': [str(f) for f in canonical],
        'active_linear_forms': [str(s.factor(f)) for f in forms],
        'complement_old_indices': [1, 3],
        'coordinate_matrix': matrix_strings(T),
        'coordinate_inverse': matrix_strings(Ti),
        'coordinate_matrix8': matrix_strings(T8),
        'coordinate_inverse8': matrix_strings(T8i),
        'coordinate_determinant': str(coordinate_det),
        'new_linear_forms_in_old8': [str(s.factor(f)) for f in T8*s.Matrix(x)],
        'old_variables_in_new8': [str(s.factor(f)) for f in T8i*s.Matrix(y8)],
        'generator_matrix': matrix_strings(generator),
        'generator_inverse': matrix_strings(generator_inverse),
        'generator_determinant': str(generator_det),
        'checks': {
            'all_six_QQ_coefficient_identities': True,
            'two_sided_generator_inverse': True,
            'two_sided_seven_variable_inverse': True,
            'two_sided_eight_variable_inverse': True,
            'all_matrix_denominators_on_stated_chart': True,
            'heavy_and_two_complement_variables_absent': True,
        },
        'input_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in paths},
        'seconds': time.time()-begin,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = verify()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'all_passed': report['all_passed'], 'seconds': report['seconds']}))


if __name__ == '__main__':
    main()

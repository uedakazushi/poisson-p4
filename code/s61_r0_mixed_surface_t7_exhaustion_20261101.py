#!/usr/bin/env python3
"""Remove the complete T7 mixed-surface orbit over sufficiently general beta."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import sympy as s
from sympy.polys.matrices import DomainMatrix as DM
from sympy.polys.rings import ring

from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_model import S61
from weight1_structure import COORDS, infinitesimal_action, sp_basis

ROOT = Path(__file__).resolve().parents[1]


def enc(matrix):
    return [[str(v) for v in row] for row in matrix.tolist()]


def run():
    started = time.time()
    old = json.loads((ROOT / 'reports/S61_R0_MIXED_SURFACE_T7_FORMAL_20261101.json').read_text())
    saved = old['universal']
    orbit = json.loads((ROOT / 'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json').read_text())['orbits']['T7']
    model = S61()
    one = [i for i, w in enumerate(model.reduced['parameter_weights']) if w == 1]
    two = [i for i, w in enumerate(model.reduced['parameter_weights']) if w == 2]
    bb = s.symbols('b:11')
    field = s.QQ.frac_field(*bb)
    poly, *b = ring(bb, s.QQ)
    H = s.Matrix(saved['beta_kernel']).applyfunc(s.Rational)
    beta_value = H * s.Matrix(bb)
    q0 = s.Matrix([s.Rational(orbit['point61'][i]) for i in two])
    fixed = next(i for i, v in enumerate(q0) if v)
    chart = [i for i in range(45) if i != fixed]
    _, equations = specialize(model, dict(zip(one, beta_value)), True)
    J = s.zeros(135, 45)
    terms = []
    for row, equation in enumerate(equations):
        ts = []
        value = 0
        for mon, coefficient in equation.items():
            if len(mon) == 2:
                i, j = mon
                value += coefficient * q0[i] * q0[j]
                if i != fixed:
                    J[row, chart.index(i)] += coefficient * q0[j]
                if j != fixed:
                    J[row, chart.index(j)] += coefficient * q0[i]
            elif len(mon) == 1:
                J[row, 44] += coefficient * q0[mon[0]]
            if fixed not in mon:
                i, j = [chart.index(k) for k in mon] + [44] * (2 - len(mon))
                ts.append((i, j, poly.from_expr(coefficient)))
        assert s.expand(value) == 0
        terms.append(ts)
    JD = J.to_DM().convert_to(field)
    assert (JD - s.Matrix(saved['jacobian']).applyfunc(s.sympify).to_DM().convert_to(field)).is_zero_matrix
    rr, pivot_rows = JD.transpose().rref(method='GJ')
    assert len(pivot_rows) == 42
    pivot_columns = saved['pivot_columns']
    assert [i for i in range(45) if i not in pivot_columns] == [19, 24, 43]
    determinant = JD.extract(pivot_rows, pivot_columns).det()
    expected_minor = -s.Rational(304054179985577710014910857, 10**14) * bb[8]**14 * bb[9]**2
    assert determinant == field.from_sympy(expected_minor)
    rr = rr.to_Matrix()
    selected_row = 103
    assert selected_row not in pivot_rows
    left = s.zeros(1, 135)
    left[0, selected_row] = 1
    for k, row in enumerate(pivot_rows):
        left[0, row] = -rr[k, selected_row]
    assert (left.to_DM().convert_to(field) * JD).is_zero_matrix

    # The common denominators are monomials. Work in QQ[b], never expand
    # products of unsimplified rational expressions.
    denominators = {1: b[8]**2*b[9]**2, 2: b[8]**4*b[9]**5, 3: b[8]**6*b[9]**8}
    left_denominator = b[8]**3

    def cleared(expression, denominator):
        numerator, den = s.fraction(s.sympify(expression))
        return (poly.from_expr(numerator) * denominator).exquo(poly.from_expr(den))

    corrections = {}
    for degree, key, column in [(1, 'tangent', 2), (2, 'quadratic_correction', 5), (3, 'cubic_correction', 9)]:
        corrections[degree] = [cleared(row[column], denominators[degree]) for row in saved[key]]
        assert corrections[degree][44] == 0
        for entry in {value for row in saved[key] for value in row}:
            den = poly.from_expr(s.fraction(s.sympify(entry))[1])
            assert len(den) == 1
            mon = next(iter(den))
            assert all(exponent == 0 for i, exponent in enumerate(mon) if i not in [8, 9])
    assert [corrections[1][i] for i in [19, 24, 43]] == [poly.zero, poly.zero, denominators[1]]
    assert all(corrections[d][i] == 0 for d in [2, 3] for i in [19, 24, 43])
    jpoly = [[poly.from_expr(value) for value in row] for row in J.tolist()]
    for row, ts in enumerate(terms):
        linear = sum((jpoly[row][i] * corrections[1][i] for i in range(45)), poly.zero)
        assert linear == 0
        quadratic = sum((co*corrections[1][i]*corrections[1][j] for i, j, co in ts), poly.zero)
        jsecond = sum((jpoly[row][i] * corrections[2][i] for i in range(45)), poly.zero)
        assert quadratic * denominators[2] + jsecond * denominators[1]**2 == 0
        cubic = sum((co*(corrections[1][i]*corrections[2][j] + corrections[2][i]*corrections[1][j]) for i, j, co in ts), poly.zero)
        jthird = sum((jpoly[row][i] * corrections[3][i] for i in range(45)), poly.zero)
        assert cubic * denominators[3] + jthird * denominators[1] * denominators[2] == 0
    print('raw135 Jacobian and pure 3-jets PASS', time.time() - started, flush=True)

    bilinear = {}
    for row, lv in enumerate(left):
        if not lv:
            continue
        ln = cleared(lv, left_denominator)
        for i, j, co in terms[row]:
            ij = tuple(sorted((i, j)))
            bilinear[ij] = bilinear.get(ij, poly.zero) + ln * co
    bilinear = {ij: co for ij, co in bilinear.items() if co}
    assert denominators[2]**2 == denominators[1] * denominators[3]
    numerator = sum((co*(corrections[2][i]*corrections[2][j] + corrections[1][i]*corrections[3][j] + corrections[3][i]*corrections[1][j]) for (i, j), co in bilinear.items()), poly.zero)
    P = 19800*bb[0]*bb[5]*bb[8] + 4480*bb[0]*bb[9]**2 - 1144*bb[10]*bb[8]*bb[9] - 1000*bb[2]*bb[8]*bb[9] - 5940*bb[6]*bb[8]**2 - 2079*bb[7]*bb[8]**2
    constant = -s.Rational(3709319, 44470277627904)
    assert numerator == poly.from_expr(constant * P**4) * b[8]**8 * b[9]**2
    quartic = constant * P**4 / (bb[8]**3 * bb[9]**8)
    assert quartic.subs(dict(zip(bb, range(2, 13)))) == s.sympify(old['point']['quartic_coefficient'])
    assert s.degree(P, bb[6]) == 1 and s.gcd(s.diff(P, bb[6]), P.subs(bb[6], 0)) == 1
    print('normal quartic fourth-power factorization PASS', time.time() - started, flush=True)

    # Reconstruct the full beta source differential using the intrinsic
    # Sp4 action, independently of any saved image-dimension assertion.
    z, variables, primitive = primitive_beta(model)
    tensor = s.Matrix([primitive[i, j].diff(z[k]) for (i, j), k in COORDS])
    basis = tensor.jacobian(variables)
    rows = list(basis.T.rref()[1])
    inverse = basis[rows, :].inv()
    action = s.Matrix.hstack(*[inverse * infinitesimal_action(tensor, A)[rows, :] for A in sp_basis()])
    action = action.xreplace(dict(zip(variables, beta_value)))
    annihilator = s.Matrix.hstack(*H.T.nullspace()).T
    normal_action = (annihilator * action).applyfunc(s.expand)
    normal_columns = [2, 4, 5, 6, 7]
    assert all(normal_action[:, i] == s.zeros(5, 1) for i in range(10) if i not in normal_columns)
    normal_minor = s.factor(normal_action[:, normal_columns].det())
    assert s.expand(normal_minor - s.Rational(27, 112000) * bb[8]**2 * P) == 0
    orbit_derivative = s.Matrix(orbit['q_orbit_differential']).applyfunc(s.Rational)
    kernel_columns = [i for i in range(10) if i not in normal_columns]
    surface_tangent = orbit_derivative[[19, 24, 43], kernel_columns]
    assert surface_tangent == s.Matrix([[0, 0, 0, 0, 1], [0, -2, 0, -4, 0], [0, 0, 0, 0, 0]])
    assert surface_tangent.rank() == 2

    exceptions = []
    for label, factor, eliminate, expected_rank in [('b8', bb[8], 8, 14), ('b9', bb[9], 9, 15), ('P', P, 6, 15)]:
        coordinate = bb[eliminate]
        value = s.solve(factor, coordinate)[0]
        others = [v for v in bb if v != coordinate]
        substitution = {coordinate: value}
        restricted_beta = beta_value.xreplace(substitution)
        source_derivative = restricted_beta.jacobian(others).row_join(action.xreplace(substitution))
        DD = source_derivative.to_DM().convert_to(field)
        reduced, pivots = DD.transpose().rref(method='GJ')
        assert len(pivots) == expected_rank
        reduced = reduced.to_Matrix()
        free = [i for i in range(16) if i not in pivots]
        kernel = s.zeros(len(free), 16)
        for i, f in enumerate(free):
            kernel[i, f] = 1
            for j, p in enumerate(pivots):
                kernel[i, p] = -reduced[j, f]
        assert (kernel.to_DM().convert_to(field) * DD).is_zero_matrix
        sample = dict(zip(others, range(2, 12)))
        sample_matrix = source_derivative.subs(sample)
        sample_columns = list(sample_matrix.rref()[1])
        sample_rows = list(sample_matrix[:, sample_columns].T.rref()[1])
        sample_minor = sample_matrix[sample_rows, sample_columns].det()
        assert len(sample_columns) == expected_rank and sample_minor != 0
        assert s.degree(factor, coordinate) == 1
        assert s.gcd(s.diff(factor, coordinate), factor.subs(coordinate, 0)) == 1
        exceptions.append(dict(label=label, factor=str(factor), elimination_variable=str(coordinate), elimination_value=str(value), geometrically_irreducible=True, beta_image_dimension=expected_rank, full_source_derivative=enc(source_derivative), left_kernel=enc(kernel), sample_values={str(k): str(v) for k, v in sample.items()}, sample_minor_rows=sample_rows, sample_minor_columns=sample_columns, sample_minor=str(sample_minor)))
    print('all three bad images proper PASS', time.time() - started, flush=True)

    paths = [Path(__file__), ROOT/'reports/S61_R0_MIXED_SURFACE_T7_FORMAL_20261101.json', ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json', ROOT/'reports/S61_R0_MIXED_SURFACE_T7_REVIEW_20261101.json', ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json', ROOT/'code/s61_model.py', ROOT/'code/s61_flag_transverse_obstruction_20261101.py', ROOT/'code/primitive_jet_elliptic_curve_20261101.py', ROOT/'code/weight1_structure.py', ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_R0_MIXED_SURFACE_T7_EXHAUSTION_20261101', all_passed=True, beta_kernel=enc(H), chart_fixed_q_index=fixed, chart_columns=chart, jacobian_rank=42, minor_rows=list(pivot_rows), minor_columns=pivot_columns, jacobian_minor=str(expected_minor), tangent_free_columns=[19,24,43], left_quartic_row=selected_row, left_kernel_row=list(map(str,left)), all_frame_denominators_supported_on=['b8','b9'], pure_first_second_third_jets_verified_all135=True, pure_jet_homogenizer_coefficients_zero=True, quartic_coefficient=str(quartic), quartic_numerator_terms=len(numerator), critical_cubic=str(P), normal_source_action=enc(normal_action), normal_source_minor_columns=normal_columns, normal_source_minor=str(normal_minor), actual_surface_tangent=enc(surface_tangent), exceptions=exceptions, all_bad_factors_exhausted=['b8','b9','P'], completed_ring_on_open='C[[u0,u1,e]]/(e^4)', homogenizer_zero_on_open=True, whole_T7_orbit_removed_over_sufficiently_general_beta=True, limitations=['The assertion is for sufficiently general beta, not every beta.', 'T5, the binary strata, other minimal tails, and full affine exhaustion are separate.'], seconds=time.time()-started, input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'all_passed': result['all_passed'], 'seconds': result['seconds']}))

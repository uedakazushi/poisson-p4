#!/usr/bin/env python3
"""A B-free weighted implicit frame, directly checked against the original 135 equations."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import sympy as s
from sympy.polys.rings import ring

from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / 'research/notes/S61_BINARY_HYPERPLANE_BFREE_FRAME_20261102.md'


def verify():
    begin = time.time()
    chart_path = ROOT / 'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json'
    light_path = ROOT / 'reports/S61_BINARY_HYPERPLANE_LINEAR_LIGHT_FRAME_20261102.json'
    old_path = ROOT / 'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json'
    chart = json.loads(chart_path.read_text())
    light_data = json.loads(light_path.read_text())
    old_data = json.loads(old_path.read_text())
    for dep in [chart, light_data, old_data]:
        assert dep['all_passed']
        for name, expected in dep.get('input_sha256', {}).items():
            assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected
    old = old_data['linear']
    R, a, e, d = ring('a,e,d', s.QQ)
    K = R.to_field()
    polynomial = lambda value: R.from_expr(s.sympify(value))
    fraction = lambda value: K.from_expr(s.sympify(value))
    F = 500*a + 90*e**2 + 3*e - 24
    B = 50*a*d + 18*d*e**2 - 9*d*e - 2*d + 60*e - 40
    C = 100*a + 18*e**2 - 9*e - 2
    D0 = 50*a + 9*e**2 + 3*e - 6
    S = 3*e - 2
    P = 125*a**2*d - 450*a*e - 200*a - 54*e**3 - 27*e**2 + 33*e + 6
    Q = 15*a*d*e - 10*a*d + 200*a + 18*e**2 + 6*e - 12
    assert 2*P + (6*e+1)*Q == 5*a*B
    assert 25*a*Q - S*P == D0*C

    T = [[polynomial(v) for v in row] for row in light_data['T_new']]
    U = [[fraction(v) for v in row] for row in chart['U34']]
    assert len(T) == len(U) == 45
    assert all(len(row) == 8 for row in T)
    assert all(len(row) == 34 for row in U)
    k = [K(row[5]) for row in T]
    assert [k[i] for i in [4, 12, 35]] == [K(P)/250, K.zero, K(Q)/50]
    lam = [K.zero]*45
    lam[4] = 250*K(S)
    lam[35] = -1250*K(a)
    lamk = sum((lam[i]*k[i] for i in [4, 35]), K.zero)
    assert lamk == -K(D0*C)
    lamU = [sum((lam[i]*U[i][j] for i in [4, 35]), K.zero) for j in range(34)]
    Unew = [[U[i][j] - k[i]*lamU[j]/lamk for j in range(34)] for i in range(45)]
    assert all(sum((lam[i]*Unew[i][j] for i in [4, 35]), K.zero) == 0
               for j in range(34))

    # Direct original equations; no partial-elimination model or raw jet cache.
    model = S61()
    beta = {int(key[1:])-1: s.sympify(value) for key, value in chart['beta'].items()}
    _, equations = specialize(model, beta, True)
    fixed = old['fixed']
    cols = old['cols']
    weights = chart['local_coordinate_weights']
    equation_weights = chart['equation_weights']
    normal_weights = chart['normal_weights34']
    J = [{} for _ in range(135)]
    quadratic_terms = 0
    for i, eq in enumerate(equations):
        for mon, value in eq.items():
            coefficient = K(polynomial(value))
            if mon == (fixed, fixed):
                assert not coefficient
            if len(mon) == 2:
                x, y = mon
                if x != fixed and y == fixed:
                    j = cols.index(x)
                    J[i][j] = J[i].get(j, K.zero) + coefficient
                if y != fixed and x == fixed:
                    j = cols.index(y)
                    J[i][j] = J[i].get(j, K.zero) + coefficient
            elif len(mon) == 1 and mon[0] == fixed:
                J[i][44] = J[i].get(44, K.zero) + coefficient
            if fixed not in mon and coefficient:
                indices = [cols.index(j) for j in mon] + [44]*(2-len(mon))
                assert sum(weights[j] for j in indices) == equation_weights[i]
                quadratic_terms += 1
        assert all(not c or weights[j] == equation_weights[i] for j, c in J[i].items())
    assert len(equations) == 135
    selected = chart['selected34_rows']
    assert len(selected) == 34
    assert all(sum((c*K(T[j][h]) for j, c in row.items()), K.zero) == 0
               for row in J for h in range(8))
    for r, i in enumerate(selected):
        for h in range(34):
            assert sum((c*U[j][h] for j, c in J[i].items()), K.zero) == int(r == h)
            assert sum((c*Unew[j][h] for j, c in J[i].items()), K.zero) == int(r == h)

    # Preserve the old light coordinate functionals and all three linear cuts.
    W = s.Matrix(old['chart_basis']).applyfunc(s.sympify)
    Wi = W.inv()
    normal_rows = [chart['free34'][j] for j in chart['normal34']]
    cut_rows = [chart['free34'][j] for j in chart['family_parameters']]
    old_pi = [[fraction(v) for v in Wi[row, :]] for row in normal_rows]
    cuts = [[fraction(v) for v in Wi[row, :]] for row in cut_rows]
    light = [0, 1, 2, 3, 4, 6, 7]
    transition = s.Matrix([
        [sum((old_pi[i][r]*K(T[r][j]) for r in range(45)), K.zero).as_expr()
         for j in light] for i in light])
    light_det = s.factor(transition.to_DM().to_field().det().as_expr())
    assert polynomial(light_det) == 2560000000*a
    light_inverse = transition.inv().applyfunc(s.cancel)
    Pi = [[K.zero]*45 for _ in range(8)]
    for i, outrow in enumerate(light):
        for j in range(45):
            Pi[outrow][j] = sum((fraction(light_inverse[i, h])*old_pi[light[h]][j]
                                for h in range(7)), K.zero)
    Pi[5] = [c/lamk for c in lam]
    for i in range(8):
        for h in range(8):
            assert sum((Pi[i][j]*K(T[j][h]) for j in range(45)), K.zero) == int(i == h)
        for h in range(34):
            assert sum((Pi[i][j]*Unew[j][h] for j in range(45)), K.zero) == 0
    assert all(sum((row[j]*K(T[j][h]) for j in range(45)), K.zero) == 0
               for row in cuts for h in range(8))
    assert all(sum((row[j]*Unew[j][h] for j in range(45)), K.zero) == 0
               for row in cuts for h in range(34))
    assert all(sum((old_pi[i][j]*(Unew[j][h]-U[j][h]) for j in range(45)), K.zero) == 0
               for i in light for h in range(34))

    # Factor-free checks of the claimed localized ring for every denominator.
    common_U = R.one
    common_Pi = R.one
    for row in Unew:
        for c in row:
            assert not (a**2*F*D0*C).div(c.denom)[1]
            assert c.denom.degree(d) == 0
            common_U = common_U.lcm(c.denom)
    for row in Pi:
        for c in row:
            assert not (a*D0*C).div(c.denom)[1]
            assert c.denom.degree(d) == 0
            common_Pi = common_Pi.lcm(c.denom)
    for i, row in enumerate(T):
        assert all(not c or weights[i] == normal_weights[j] for j, c in enumerate(row))
    for i, row in enumerate(Unew):
        assert all(not c or weights[i] == equation_weights[selected[j]]
                   for j, c in enumerate(row))
    for i, row in enumerate(Pi):
        assert all(not c or weights[j] == normal_weights[i] for j, c in enumerate(row))

    # Explicit identity on a constant basis of the 42-dimensional cut space.
    keep = [j for j in range(45) if j not in cut_rows]
    E = [[fraction(W[i, j]) for j in keep] for i in range(45)]
    assert len(keep) == 42
    assert all(sum((row[j]*E[j][h] for j in range(45)), K.zero) == 0
               for row in cuts for h in range(42))
    PiE = [[sum((Pi[i][j]*E[j][h] for j in range(45)), K.zero)
            for h in range(42)] for i in range(8)]
    JE = [[sum((c*E[j][h] for j, c in J[i].items()), K.zero)
           for h in range(42)] for i in selected]
    for i in range(45):
        for h in range(42):
            value = sum((K(T[i][j])*PiE[j][h] for j in range(8)), K.zero)
            value += sum((Unew[i][j]*JE[j][h] for j in range(34)), K.zero)
            assert value == E[i][h]

    paths = [Path(__file__), NOTE, chart_path, light_path, old_path,
             ROOT/'code/s61_model.py', ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
             ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',
             ROOT/'code/certificates/s61_binary_hyperplane_linear_light_frame_20261102.json']
    paths += [ROOT/'legacy/v02/data'/name for name in [
        'symplectic_slice_certificate.json', 'symplectic_slice_equations.json',
        'symplectic_slice_reduced_equations.json']]
    strings = lambda matrix: [[str(c.as_expr()) for c in row] for row in matrix]
    return dict(
        id='S61_BINARY_HYPERPLANE_BFREE_FRAME_20261102', all_passed=True,
        core_variables=['a', 'e', 'd'], open_factors=[str(v.as_expr()) for v in [a, F, D0, C]],
        original_equations=135, quadratic_terms_checked=quadratic_terms,
        selected_rows=selected, normal_weights=normal_weights,
        cut_rows_in_constant_chart=cut_rows, cut_space_dimension=42,
        polynomials={name: str(value.as_expr()) for name, value in
                     dict(F=F, B=B, C=C, D0=D0, S=S, P=P, Q=Q).items()},
        heavy_kernel_coordinates=[str(k[i].as_expr()) for i in [4, 12, 35]],
        lambda_nonzero_entries={str(i): str(lam[i].as_expr()) for i in [4, 35]},
        lambda_kernel=str(lamk.as_expr()), light_determinant=str(light_det),
        light_transition=[[str(c) for c in row] for row in transition.tolist()],
        T_new=strings([[K(c) for c in row] for row in T]),
        U_new=strings(Unew), Pi_new=strings(Pi), cut_matrix=strings(cuts),
        U_common_denominator=str(s.factor(common_U.as_expr())),
        Pi_common_denominator=str(s.factor(common_Pi.as_expr())),
        all_denominators_d_independent=True, original_JT_zero=True,
        original_base_point_equations_zero=True,
        original_selected_JUold_identity=True, original_selected_JUnew_identity=True,
        PiT_identity=True, PiU_zero=True, all_cuts_preserved=True,
        old_light_projection_preserved=True, cut_space_inverse_identity=True,
        all_original_weights_checked=True, all_frame_weights_checked=True,
        maximum_coordinate_weight=max(weights), maximum_equation_weight=max(equation_weights),
        exact_implicit_solution_degree_bound=4, exact_residual_degree_bound=6,
        proof_note=str(NOTE.relative_to(ROOT)),
        limitations=[
            'The new open is a*F*D0*C nonzero; D0 is an additional factor, while B is not inverted.',
            'This is a weighted coordinate and finite-reconstruction theorem, not uniform weighted rank, length65, or ambient removal.',
            'No partial10-variable work model or expanded nonlinear jet cache is used as a proof input.'],
        input_sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                      for path in paths}, seconds=time.time()-begin)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    result = verify()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({key: result[key] for key in
                      ['all_passed', 'U_common_denominator', 'Pi_common_denominator', 'seconds']}))

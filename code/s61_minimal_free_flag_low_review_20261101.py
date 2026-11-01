#!/usr/bin/env python3
"""Independent Lie-flow and raw-Jacobian/quadratic checks; replay the cubic producer."""
import argparse
import hashlib
import json
import time
from itertools import combinations_with_replacement
from pathlib import Path

import sympy as s

from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_model import S61
from weight1_structure import COORDS, infinitesimal_action, sp_basis

ROOT = Path(__file__).resolve().parents[1]


def run():
    started = time.time()
    names = ['NORMAL_FORM', 'UNIVERSAL_QUADRATIC', 'UNIVERSAL_CUBIC']
    paths = [ROOT / f'reports/S61_MINIMAL_FREE_FLAG_{name}_20261101.json' for name in names]
    normal, quadratic, cubic = [json.loads(path.read_text()) for path in paths]
    checks = 0
    for data in (normal, quadratic, cubic):
        assert data['all_passed']
        for path, expected in data['input_sha256'].items():
            assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == expected, path
            checks += 1
    model = S61()
    z, u, beta = primitive_beta(model)
    vector = s.Matrix([beta[i, j].diff(z[k]) for (i, j), k in COORDS])
    basis = vector.jacobian(u)
    coefficient_rows = list(basis.T.rref()[1])
    inverse = basis[coefficient_rows, :].inv()
    action = s.Matrix.hstack(*[inverse * infinitesimal_action(vector, A)[coefficient_rows, :]
                              for A in sp_basis()])
    old = json.loads((ROOT / 'reports/S61_FLAG_UNIFORM_OPEN_20261101.json').read_text())
    H = s.Matrix(old['flag']['beta_kernel']).applyfunc(s.Rational)
    stabilizer = s.Matrix(old['frame_change']['stabilizer']).applyfunc(s.Rational)
    b = s.symbols('b0:15')
    t = s.Symbol('t')
    base = (H * s.Matrix(b)).subs({b[0]: 0, b[1]: 0, b[2]: 0, b[3]: -b[12]})
    for j in (0, 2, 3):
        generator = (action * stabilizer[:, j]).jacobian(u)
        assert generator ** 4 == s.zeros(16)
        finite = sum((t ** d * generator ** d / s.factorial(d) for d in range(4)), s.zeros(16))
        expected = s.Matrix(normal['finite_flows'][str(j)]).applyfunc(s.sympify)
        assert (finite * base)[1:, :].applyfunc(s.expand) == expected.applyfunc(s.expand)
    for bad in normal['proper_images']:
        source = s.Matrix(bad['source_beta']).applyfunc(s.sympify)
        parameters = list(map(s.Symbol, bad['parameters']))
        differential = source.jacobian(parameters).row_join(action.subs(dict(zip(u, source))))
        assert differential == s.Matrix(bad['beta_differential']).applyfunc(s.sympify)
        kernel = s.Matrix(bad['left_kernel']).applyfunc(s.sympify)
        assert (kernel.T * differential).applyfunc(s.expand) == s.zeros(1, differential.cols)
        witness = bad['rank_witness']
        sample = dict(zip(parameters, range(1, len(parameters) + 1)))
        block = differential.subs(sample)[witness['rows_zero_based'], witness['columns_zero_based']]
        assert block.det() == s.Rational(witness['determinant']) != 0
    a = s.symbols('a0:5')
    source = s.Matrix(normal['normalized_beta16']).applyfunc(s.sympify)
    dominant = source.jacobian(a).row_join(action.subs(dict(zip(u, source)))).row_join(source)
    assert dominant.subs(dict(zip(a, range(1, 6)))).det() == -1944

    point = json.loads((ROOT / 'reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json').read_text())['strata']['free_M0']['local_data']
    one = [i for i, weight in enumerate(model.reduced['parameter_weights']) if weight == 1]
    two, equations = specialize(model, dict(zip(one, source)), True)
    q0 = s.Matrix(point['point61']).applyfunc(s.Rational)[two, :]
    fixed = two.index(point['chart_fixed_u'] - 1)
    chart = point['chart_columns']
    J = s.zeros(135, 45)
    T = s.Matrix(quadratic['tangent_basis_after_shift']).applyfunc(s.sympify)
    U = s.Matrix(quadratic['right_inverse']).applyfunc(s.sympify)
    mons = list(combinations_with_replacement(range(11), 2))
    H2 = s.zeros(135, len(mons))
    for row, equation in enumerate(equations):
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
                for column, (v, w) in enumerate(mons):
                    H2[row, column] += coefficient * (T[i, v] * T[j, w]
                        + (T[i, w] * T[j, v] if v != w else 0))
        assert s.expand(value) == 0
    assert J == s.Matrix(quadratic['jacobian']).applyfunc(s.sympify)
    assert (J * T).applyfunc(s.expand) == s.zeros(135, 11)
    assert T[-1, :] == s.zeros(1, 11)
    rows = quadratic['pivot_rows']
    assert (J[rows, :] * U).applyfunc(s.cancel) == s.eye(34)
    actual = (H2.to_DM() - J.to_DM() * U.to_DM() * H2[rows, :].to_DM()).to_Matrix()
    coefficients = s.Matrix(quadratic['quadratic_row_coefficients']).applyfunc(s.sympify)
    constant = s.Matrix(quadratic['constant_quadratic_rref']).applyfunc(s.Rational)
    assert (actual.to_DM() - coefficients.to_DM() * constant.to_DM()).is_zero_matrix
    assert constant.rank() == 12
    hessians = s.zeros(132, 11)
    for row in range(12):
        for column, (i, j) in enumerate(mons):
            hessians[11 * row + i, j] += constant[row, column]
            hessians[11 * row + j, i] += constant[row, column]
    assert hessians.rank() == 8 and hessians[:, :3] == s.zeros(132, 3)

    # This final step is explicitly a replay, not an independent cubic contraction.
    from s61_minimal_free_flag_universal_cubic_20261101 import run as cubic_run
    replay = cubic_run()
    for key in ['quadratic_product_rank', 'new_cubic_rank', 'true_initial_hilbert_function',
                'constant_cubic_polynomials', 't3_normal_form', 't_order']:
        assert replay[key] == cubic[key], key
    sources = paths + [Path(__file__), ROOT / 'code/s61_minimal_free_flag_universal_cubic_20261101.py',
                       ROOT / 'code/s61_model.py', ROOT / 'code/s61_flag_transverse_obstruction_20261101.py',
                       ROOT / 'code/primitive_jet_elliptic_curve_20261101.py', ROOT / 'code/weight1_structure.py']
    sources += [ROOT / f'research/notes/S61_MINIMAL_FREE_FLAG_{name}_20261101.md' for name in names]
    return dict(id='S61_MINIMAL_FREE_FLAG_LOW_REVIEW_20261101', all_passed=True,
                input_hashes_checked=checks, independent_finite_lie_flows=3,
                independent_proper_image_left_kernels=3, independent_dominance_minor='-1944',
                raw135_jacobian_and_quadratics_reconstructed=True, quadratic_rank=12, common_hessian_rank=8,
                cubic_producer_replayed=True, cubic_coefficients_independently_reimplemented=False,
                accepted=['Finite normal form covers all stated endpoint points over general beta.',
                          'Function-field constant quadratic and true cubic initial ideals; corrected t order3.'],
                limitations=['No independent full cubic contraction or finite fourth-flow reimplementation.',
                             'No higher initial ideal, t-nilpotence, flatness, or endpoint saturation.'],
                proof_note='research/notes/S61_MINIMAL_FREE_FLAG_LOW_REVIEW_20261101.md',
                input_sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in sources}, seconds=time.time() - started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    report = run()
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print(json.dumps(dict(all_passed=True, seconds=report['seconds'])))

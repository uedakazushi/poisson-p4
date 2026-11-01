#!/usr/bin/env python3
"""Eight core-divisor witnesses from the raw 135 equations, over F_32003.

The characteristic-zero conclusion uses Hensel lifting and the proof note.
No finite-field grid or ignored work file is an input.
"""
import argparse
import hashlib
import json
from itertools import combinations_with_replacement
from pathlib import Path
import time

import numpy as np
import sympy as s

from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_binary_deeper_uniform_quadratic_20261101 import elimination_rank, determinant_mod

ROOT = Path(__file__).resolve().parents[1]
PRIME = 32003
NUISANCE = [5, 7, 11, 13, 17]
VARIABLES = s.symbols('e d p q r h k')
MONS = {j: list(combinations_with_replacement(range(9), j)) for j in (2, 3, 4)}
INDEX = {j: {m: i for i, m in enumerate(ms)} for j, ms in MONS.items()}
MULTIPLY = {j: np.array([[INDEX[j + 1][tuple(sorted((a,) + m))] for m in MONS[j]]
                       for a in range(9)]) for j in (2, 3)}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mod(value):
    value = s.Rational(value)
    assert int(value.q) % PRIME
    return int(value.p) % PRIME * pow(int(value.q), -1, PRIME) % PRIME


def matrix_sha(matrix):
    return hashlib.sha256(json.dumps(matrix.tolist(), separators=(',', ':')).encode()).hexdigest()


def reduce_rows(matrix, cut=None, full=False):
    """Tracked row operations on the entire matrix, pivoting only before cut."""
    matrix = matrix.copy() % PRIME
    rank, pivots = 0, []
    limit = cut or matrix.shape[1]
    for column in range(limit):
        choices = np.flatnonzero(matrix[rank:, column])
        if not len(choices):
            continue
        chosen = rank + int(choices[0])
        matrix[[chosen, rank]] = matrix[[rank, chosen]]
        matrix[rank] = matrix[rank] * pow(int(matrix[rank, column]), -1, PRIME) % PRIME
        targets = (np.flatnonzero(matrix[:, column]) if full else
                   rank + 1 + np.flatnonzero(matrix[rank + 1:, column]))
        targets = targets[targets != rank]
        if len(targets):
            matrix[targets] = (matrix[targets] - matrix[targets, column, None] * matrix[rank, None, :]) % PRIME
        pivots.append(column)
        rank += 1
        if rank == len(matrix):
            break
    return matrix, rank, pivots


def inverse(matrix):
    size = len(matrix)
    augmented, rank, pivots = reduce_rows(
        np.concatenate([matrix, np.eye(size, dtype=np.int64)], axis=1), size, True)
    assert rank == size and pivots == list(range(size))
    return augmented[:, size:]


def rank_witness(matrix, expected):
    rank, columns = elimination_rank(matrix, PRIME)
    assert rank == expected
    rank2, rows = elimination_rank(matrix[:, columns].T, PRIME)
    assert rank2 == expected
    determinant = determinant_mod(matrix[np.ix_(rows, columns)], PRIME)
    assert determinant
    return dict(rank=rank, rows=rows, columns=columns, determinant_mod_prime=determinant)


def multiples(matrix, degree):
    result = np.zeros((len(matrix), 9, len(MONS[degree + 1])), dtype=np.int64)
    for a in range(9):
        result[:, a, MULTIPLY[degree][a]] = matrix
    return result.reshape(len(matrix) * 9, -1)


def prepare(linear, model):
    e, d = VARIABLES[:2]
    nuisance = dict(zip(VARIABLES[2:], NUISANCE))

    def encode(expression):
        numerator, denominator = s.cancel(s.sympify(expression).subs(nuisance)).as_numer_denom()
        return [[(monomial, mod(coefficient)) for monomial, coefficient in s.Poly(poly, e, d).terms()]
                for poly in (numerator, denominator)]

    raw = {key: [[encode(value) for value in row] for row in linear[key]] for key in ('T', 'U', 'J')}
    beta = {int(name[1:]) - 1: s.sympify(value) for name, value in linear['beta'].items()}
    _, equations = specialize(model, {i: value.subs(nuisance) for i, value in beta.items()}, True)
    assert len(equations) == 135
    raw['H'] = []
    for i, equation in enumerate(equations):
        for monomial, coefficient in equation.items():
            if linear['fixed'] in monomial:
                continue
            a, b = [linear['cols'].index(v) for v in monomial] + [44] * (2 - len(monomial))
            raw['H'].append((i, a, b, encode(coefficient if a == b else coefficient / 2)))
    return raw, beta


def evaluate(raw, core):
    e, d = core

    def value(fraction):
        parts = [sum(c * pow(e, m[0], PRIME) * pow(d, m[1], PRIME) for m, c in poly) % PRIME
                 for poly in fraction]
        assert parts[1]
        return parts[0] * pow(parts[1], -1, PRIME) % PRIME

    T, U, J = [np.array([[value(v) for v in row] for row in raw[key]], dtype=np.int64)
               for key in ('T', 'U', 'J')]
    H = np.zeros((135, 45, 45), dtype=np.int64)
    for i, a, b, coefficient in raw['H']:
        H[i, a, b] = H[i, b, a] = value(coefficient)
    return T, U, J, H


def upper65(linear, Tfull, U, J, H):
    selected = linear['pivot_rows']
    assert not np.any(J @ Tfull % PRIME)
    assert elimination_rank(Tfull, PRIME)[0] == 12
    projection = J @ U % PRIME
    assert np.array_equal(projection[selected, :], np.eye(33, dtype=np.int64))
    linear_witness = rank_witness(J, 33)
    T = Tfull[:, linear['normal']]
    assert T.shape == (45, 9) and not np.any(T[-1])
    a = np.array([i for i, j in MONS[2]])
    b = np.array([j for i, j in MONS[2]])
    scale = np.array([1 if i == j else 2 for i, j in MONS[2]])
    polar = T.T @ (H @ T % PRIME) % PRIME
    Q = polar[:, a, b] * scale[None, :] % PRIME
    W2 = -U @ Q[selected, :] % PRIME
    A2 = (Q - projection @ Q[selected, :]) % PRIME
    polar3 = 2 * T.T @ (H @ W2 % PRIME) % PRIME
    cubic = np.zeros((135, 165), dtype=np.int64)
    for a in range(9):
        cubic[:, MULTIPLY[2][a]] = (cubic[:, MULTIPLY[2][a]] + polar3[:, a, :]) % PRIME
    A3 = (cubic - projection @ cubic[selected, :]) % PRIME

    # Constant row operations change the ACTUAL residual series, including
    # all their higher tails. No claim is made that A4 itself is zero.
    augmented, rank2, pivots2 = reduce_rows(np.concatenate([A2, A3], axis=1), 45, True)
    assert rank2 == 24
    basis2 = augmented[:rank2, :45]
    tails3 = augmented[:rank2, 45:]
    genuine3 = augmented[rank2:, 45:]
    M3 = multiples(basis2, 2)
    tails4 = multiples(tails3, 3)
    tracked, rankM3, _ = reduce_rows(np.concatenate([M3, tails4], axis=1), 165)
    assert rankM3 == 129
    # A linear syzygy of the quadratic leading forms has quartic leading
    # form given by its product with the cubic tails; raw A4 contributes
    # only in degree5 after multiplication by a linear form.
    quartics = tracked[rankM3:, 165:]
    assert quartics.shape == (87, 495)
    cubic_witness = rank_witness(np.concatenate([M3, genuine3]), 131)
    fourth = np.concatenate([multiples(M3, 3), multiples(genuine3, 3), quartics])
    fourth_witness = rank_witness(fourth, 495)
    return dict(linear=linear_witness, quadratic=rank_witness(A2, 24),
                quadratic_linear_rank=129, cubic=cubic_witness, quartic=fourth_witness,
                quadratic_pivot_columns=pivots2, tracked_linear_syzygies=87,
                raw_quadratic_matrix_sha256=matrix_sha(A2), raw_cubic_matrix_sha256=matrix_sha(A3),
                actual_quartic_initial_matrix_sha256=matrix_sha(fourth),
                initial_quotient_upper_hilbert_function=[1, 9, 21, 34], transverse_colength_upper_bound=65)


def ordinary_same_component(linear, J, H):
    # The universal P2 identity is an explicit frozen dependency. Q1 has
    # coordinates q19=q24=1 in addition to the fixed coordinate q21=1.
    offset = np.zeros(45, dtype=np.int64)
    offset[linear['cols'].index(19)] = 1
    offset[linear['cols'].index(24)] = 1
    J1 = (J + 2 * (H @ offset)) % PRIME
    leading = rank_witness(J1[:, :44], 34)
    full = rank_witness(J1, 35)
    reduced, rank, pivots = reduce_rows(J1, full=True)
    free = [i for i in range(45) if i not in pivots]
    T = np.zeros((45, 10), dtype=np.int64)
    T[free, :] = np.eye(10, dtype=np.int64)
    T[pivots, :] = -reduced[:rank, free] % PRIME
    assert not np.any(J1 @ T % PRIME) and not np.any(T[-1])
    rows = full['rows']
    projection = J1[:, pivots] @ inverse(J1[np.ix_(rows, pivots)]) % PRIME
    monomials = list(combinations_with_replacement(range(10), 2))
    polar = T.T @ (H @ T % PRIME) % PRIME
    Q = np.array([[polar[i, a, b] * (1 if a == b else 2) % PRIME for a, b in monomials]
                  for i in range(135)], dtype=np.int64)
    obstruction = (Q - projection @ Q[rows, :]) % PRIME
    quadratic = rank_witness(obstruction, 5)
    hessians = []
    for row in obstruction:
        matrix = np.zeros((10, 10), dtype=np.int64)
        for value, (a, b) in zip(row, monomials):
            matrix[a, b] += value
            matrix[b, a] += value
        hessians.extend(matrix.tolist())
    hessian = rank_witness(np.array(hessians, dtype=np.int64) % PRIME, 7)
    return dict(leading=leading, full=full, quadratic=quadratic, hessian=hessian,
                smooth_actual_binary_dimension=3, same_nonempty_P2_fiber=True,
                q1_jacobian_sha256=matrix_sha(J1))


def verify():
    started = time.time()
    linear_path = ROOT / 'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json'
    factor_path = ROOT / 'reports/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101.json'
    action_path = ROOT / 'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json'
    component_path = ROOT / 'reports/S61_BINARY_DEEPER_CORE_OPEN_FORMAL_20261101.json'
    linear = json.loads(linear_path.read_text())['linear']
    factors = json.loads(factor_path.read_text())
    component = json.loads(component_path.read_text())
    assert component['all_passed'] and component['universal_same_P2_highest_identities'] == 135
    assert [(i, v) for i, v in enumerate(component['same_P2_qfamily']) if v != '0'] == [(19, 'H'), (21, '1'), (24, 'K')]
    model = S61()
    raw, beta = prepare(linear, model)
    one = [i for i, weight in enumerate(model.reduced['parameter_weights']) if weight == 1]
    bv = s.Matrix([beta[i] for i in one])
    action = s.Matrix(json.loads(action_path.read_text())['primitive_action']).applyfunc(s.sympify)
    action = action.subs({s.Symbol('u' + str(i + 1)): beta[i] for i in one}, simultaneous=True)
    source = action.row_join(bv.jacobian(VARIABLES))
    e, d = VARIABLES[:2]
    cases = []
    for index, candidate in enumerate(factors['factor_nondivisibility_witnesses']):
        factor = s.sympify(candidate['factor'])
        core = candidate['sample_core']
        values = dict(zip(VARIABLES, core + NUISANCE))
        assert mod(factor.subs(values)) == 0
        gradient = [mod(s.diff(factor, v).subs(values)) for v in (e, d)]
        assert any(gradient)
        slice_denominator = mod((4 - 33 * e).subs(values))
        assert slice_denominator
        T, U, J, H = evaluate(raw, core)
        upper = upper65(linear, T, U, J, H)
        ordinary = ordinary_same_component(linear, J, H)
        B = np.array([[mod(v) for v in row] for row in source.subs(values).tolist()], dtype=np.int64)
        domain = np.column_stack([B[:, :10], B[:, 12:], (gradient[1] * B[:, 10] - gradient[0] * B[:, 11]) % PRIME])
        dominance = rank_witness(domain, 16)
        cj = np.array([[mod(s.sympify(v).subs(values)) for v in row] for row in linear['curve_jacobian']], dtype=np.int64)
        curve = rank_witness(cj, 2)
        family_minor = mod(s.sympify(linear['family_parameter_minor']).subs(values))
        assert family_minor
        cases.append(dict(index=index, factor=candidate['factor'], total_degree=candidate['total_degree'],
                          core_residue=core, nuisance=NUISANCE, factor_gradient_mod_prime=gradient,
                          hensel_variable='e' if gradient[0] else 'd', slice_denominator_mod_prime=slice_denominator,
                          upper=upper, ordinary_same_component=ordinary,
                          source_beta=dominance, source_beta_matrix=domain.tolist(),
                          vertex_curve=curve, family_parameter_minor_mod_prime=family_minor))
        print('case', index, 'PASS', round(time.time() - started, 3), flush=True)
    assert len(cases) == 8
    note = ROOT / 'research/notes/S61_BINARY_DEEPER_CORE_ARTIFACTS_20261101.md'
    dependencies = [Path(__file__), linear_path, factor_path, action_path, component_path, note,
                    ROOT / 'code/s61_model.py', ROOT / 'code/s61_flag_transverse_obstruction_20261101.py',
                    ROOT / 'code/s61_binary_deeper_uniform_quadratic_20261101.py',
                    ROOT / 'reports/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_REVIEW_20261101.json',
                    ROOT / 'reports/S61_BINARY_UNIFORM_RANK_20261101.json',
                    ROOT / 'reports/S61_BINARY_UNIFORM_QUADRATIC_20261101.json',
                    ROOT / 'reports/S61_BINARY_UNIFORM_CUBIC_20261101.json',
                    ROOT / 'research/notes/S61_FLAG_BINARY_CONNECTION_FORMAL_20261101.md',
                    ROOT / 'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_BINARY_DEEPER_CORE_ARTIFACTS_20261101', all_passed=True, prime=PRIME,
                raw_equations=135, cases=cases, eight_QQ_irreducible_core_divisors_removed_over_general_beta=True,
                characteristic_zero_argument='Smooth Hensel lifting; integral transverse upper65; the same relative binary C3 lower65; relative free65 and etale descent; Galois conjugacy and dimension15 complement.',
                limitations=['The common-core-open quartic minor can have further, still unprocessed factors.',
                             'No finite-field grid is used as a characteristic-zero exhaustion argument.',
                             'No claim that t is nonzero or that t^2 vanishes at every special point on these divisors.'],
                proof_note=str(note.relative_to(ROOT)),
                input_sha256={str(p.relative_to(ROOT)): sha(p) for p in dependencies}, seconds=time.time() - started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    result = verify()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False)
        stream.write('\n')
    print('PASS', round(result['seconds'], 3), flush=True)

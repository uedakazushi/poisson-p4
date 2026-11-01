#!/usr/bin/env python3
"""The exact QQ transverse Artin algebra and its nonzero homogenizing class."""
import argparse
import hashlib
import json
import time
from itertools import combinations_with_replacement
from pathlib import Path

import sympy as s

from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize


def verify():
    start = time.time()
    root = Path(__file__).resolve().parents[1]
    quadratic_path = root / 'reports/S61_BINARY_QUADRATIC_OBSTRUCTION_20261101.json'
    local_path = root / 'reports/S61_BINARY_LOCAL_DIMENSION_20261101.json'
    source = json.loads(quadratic_path.read_text())
    local = json.loads(local_path.read_text())
    assert source['all_passed'] and local['all_passed']
    assert local['quartic_rank'] == 210
    m = S61()
    point = list(map(s.Rational, source['point']))
    beta = {i: point[i] for i, w in enumerate(m.reduced['parameter_weights']) if w == 1}
    _, polys = specialize(m, beta, True)
    fixed = source['chart_fixed_q']
    cols = [i for i in range(45) if i != fixed]
    J = s.Matrix(source['jacobian']).applyfunc(s.Rational)
    normal = local['remaining_normal_variables']
    T = s.Matrix(source['tangent']).applyfunc(s.Rational)[:, normal]
    rows, pivots = local['jacobian_pivot_rows'], local['jacobian_pivot_columns']
    U = s.zeros(45, 35)
    inverse = J.extract(rows, pivots).inv(method='DM')
    for i, col in enumerate(pivots):
        for j in range(35):
            U[col, j] = inverse[i, j]
    assert J.extract(rows, list(range(45))) * U == s.eye(35)
    assert J * T == s.zeros(135, 7)
    assert T[-1, :] == s.zeros(1, 7)
    JU = J * U
    mon2 = list(combinations_with_replacement(range(7), 2))
    mon3 = list(combinations_with_replacement(range(7), 3))
    index3 = {mon: i for i, mon in enumerate(mon3)}
    Q = s.zeros(135, 28)
    hessian_terms = []
    for row, f in enumerate(polys):
        terms = []
        for mon, value in f.items():
            if fixed in mon:
                continue
            i, j = [cols.index(k) for k in mon] + [44] * (2 - len(mon))
            value = s.Rational(value)
            terms.append((i, j, value))
            for n, (a, b) in enumerate(mon2):
                Q[row, n] += value * (T[i, a] * T[j, b]
                    + (T[i, b] * T[j, a] if a != b else 0))
        hessian_terms.append(terms)
    w2 = -U * Q.extract(rows, list(range(28)))
    A2 = Q - JU * Q.extract(rows, list(range(28)))
    cubic = s.zeros(135, 84)
    for row, terms in enumerate(hessian_terms):
        for i, j, value in terms:
            for a in range(7):
                for n, mon in enumerate(mon2):
                    cubic[row, index3[tuple(sorted((a,) + mon))]] += value * (
                        T[i, a] * w2[j, n] + T[j, a] * w2[i, n])
    w3 = -U * cubic.extract(rows, list(range(84)))
    A3 = cubic - JU * cubic.extract(rows, list(range(84)))
    generators = s.Matrix.hstack(A2, A3)
    multiples = []
    for row in range(135):
        for a in range(7):
            v = [s.S.Zero] * 112
            for n, mon in enumerate(mon2):
                v[28 + index3[tuple(sorted((a,) + mon))]] += A2[row, n]
            if any(v):
                multiples.append(v)
    relations = s.Matrix.vstack(generators, s.Matrix(multiples))
    _, relation_pivots = relations.to_DM().rref()
    tclass = s.Matrix.hstack(w2[-1:, :], w3[-1:, :])
    _, with_t_pivots = s.Matrix.vstack(relations, tclass).to_DM().rref()
    assert len(relation_pivots) == 55 and len(with_t_pivots) == 56
    # t has no linear part; its quadratic part is a relation, hence t is in m^3.
    assert len(A2.to_DM().rref()[1]) == 5
    assert len(s.Matrix.vstack(A2, w2[-1:, :]).to_DM().rref()[1]) == 5
    assert 1 + 7 + 28 + 84 - len(relation_pivots) == 65
    paths = [Path(__file__), quadratic_path, local_path, root / 'code/s61_model.py',
             root / 'code/s61_flag_transverse_obstruction_20261101.py',
             root / 'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_BINARY_TRANSVERSE_ALGEBRA_20261101', all_passed=True,
        field='QQ', point=source['point'], chart_fixed_q=fixed,
        linear_parameters=local['linear_parameters'], normal_variables=normal,
        jacobian_pivot_rows=rows, jacobian_pivot_columns=pivots,
        quadratic_monomials=mon2, cubic_monomials=mon3,
        residual_quadratic_matrix=[[str(x) for x in row] for row in A2.tolist()],
        residual_cubic_matrix=[[str(x) for x in row] for row in A3.tolist()],
        homogenizing_class_degree2_and3=[str(x) for x in list(tclass)],
        truncated_relation_rank=55, rank_after_adding_t=56,
        quadratic_relation_rank=5, quadratic_rank_after_adding_t2=5,
        maximal_ideal_fourth_power_zero=True, transverse_algebra_length=65,
        homogenizing_class_nonzero=True, homogenizing_class_in_cubic_socle=True,
        transverse_infinity_quotient_length=64,
        proof_note='research/notes/S61_BINARY_TRANSVERSE_ALGEBRA_20261101.md',
        limitations=['The quotient by three specified formal parameters at one rational point.',
            'No claim that the unsliced local algebra is finite free over those parameters.',
            'The homogenizing class is nonzero in the full local ring, but nilpotence there is not proved.',
            'No t-saturation or full affine-fiber exhaustion conclusion.'],
        source_sha256={str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.time() - start)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify()
    with args.output.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps({key: result[key] for key in ['all_passed', 'transverse_algebra_length',
        'homogenizing_class_nonzero', 'transverse_infinity_quotient_length', 'seconds']}))


if __name__ == '__main__':
    main()

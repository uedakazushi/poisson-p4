#!/usr/bin/env python3
"""Independent sparse-product verification of binary local cubic obstructions.

Does not import the main specialization, Hessian, or implicit-elimination code.
Rebuilds from encoded S61, verifies exact QQ quadratic rank, and independently
forms all modular cubics and the saved quartic minor.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
from itertools import combinations_with_replacement
import json
from pathlib import Path
import time
import numpy as np
import sympy as s
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61, X, J as J5

ROOT = Path(__file__).resolve().parents[1]


def verify():
    started = time.monotonic(); prime = 32003
    quadratic_path = ROOT/'reports/S61_BINARY_QUADRATIC_OBSTRUCTION_20261101.json'
    cubic_path = ROOT/'reports/S61_BINARY_LOCAL_DIMENSION_20261101.json'
    saved_q = json.loads(quadratic_path.read_text()); saved_c = json.loads(cubic_path.read_text())
    extra_path = ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    extra = json.loads(extra_path.read_text())['universal_new_equation']
    model = S61(); one = [i for i, w in enumerate(model.reduced['parameter_weights']) if w == 1]
    two = [i for i, w in enumerate(model.reduced['parameter_weights']) if w == 2]
    old_to_new = {old: new for new, old in enumerate(two)}
    x, y, z, w = X[1:]; P = X[0]**2*J5
    P[2, 3] += x*x+w*w; P[3, 2] -= x*x+w*w
    point = model.normalize(P)['parameters']
    for a, i in enumerate(one): point[i] = s.Integer(a+2)
    point[11] = point[12] = s.S.Zero
    assert [str(v) for v in point] == saved_q['point']
    polys = []
    for f in model.reduced['equations']+[extra]:
        poly = {}
        for mon, n, d in f:
            qmon = tuple(sorted(old_to_new[i] for i in mon if i in old_to_new))
            assert len(qmon) <= 2
            c = s.Rational(n, d)*s.prod(point[i] for i in mon if i not in old_to_new)
            poly[qmon] = poly.get(qmon, s.S.Zero)+c
        polys.append({mon: c for mon, c in poly.items() if c})
    q0 = [point[i] for i in two]; fixed = saved_q['chart_fixed_q']
    assert q0[fixed] == 1
    columns = [i for i in range(45) if i != fixed]
    local_index = {old: new for new, old in enumerate(columns)}
    Jac = s.zeros(135, 45); leading_jac = s.zeros(135, 45)
    # Store each genuine quadratic chart term as (coefficient, variable i, j).
    # The last chart variable, index44, is the homogenizing variable t.
    terms = []
    for r, f in enumerate(polys):
        constant = s.S.Zero; row_terms = []
        for mon, c in f.items():
            if len(mon) == 2:
                i, j = mon; constant += c*q0[i]*q0[j]
                leading_jac[r, i] += c*q0[j]; leading_jac[r, j] += c*q0[i]
            elif len(mon) == 1:
                Jac[r, 44] += c*q0[mon[0]]
            if fixed not in mon:
                indices = [local_index[i] for i in mon]+[44]*(2-len(mon))
                row_terms.append((c, *indices))
        assert constant == 0
        terms.append(row_terms)
    Jac[:, :44] = leading_jac[:, columns]
    assert Jac == s.Matrix(saved_q['jacobian']).applyfunc(s.Rational)
    assert leading_jac.rank() == 34 and Jac.rank() == 35
    T = s.Matrix(saved_q['tangent']).applyfunc(s.Rational)
    assert T.shape == (45, 10) and Jac*T == s.zeros(135, 10) and T.rank() == 10
    assert T[44, :] == s.zeros(1, 10)
    mons2 = list(combinations_with_replacement(range(10), 2)); index2 = {m:i for i,m in enumerate(mons2)}
    mons3 = list(combinations_with_replacement(range(10), 3)); index3 = {m:i for i,m in enumerate(mons3)}
    # Direct product of the two linear forms; no symmetric Hessian convention.
    Q = s.zeros(135, 55)
    for r, row in enumerate(terms):
        for c, a, b in row:
            for i in range(10):
                for j in range(10):
                    if T[a, i] and T[b, j]:
                        Q[r, index2[tuple(sorted((i, j)))]] += c*T[a, i]*T[b, j]
    cokernel = s.Matrix.hstack(*Jac.T.nullspace()).T
    qq_projected = (cokernel.to_DM()*Q.to_DM()).to_Matrix()
    assert qq_projected.to_DM().rank() == 5

    def mod(v):
        q = s.Rational(v)
        return int(q.p)*pow(int(q.q), -1, prime) % prime
    def array(M):
        return np.array([[mod(v) for v in row] for row in M.tolist()], dtype=np.int64)
    def reduced(A):
        A = A.copy() % prime; pivots = []; r = 0
        for col in range(A.shape[1]):
            possible = np.flatnonzero(A[r:, col])
            if not len(possible): continue
            pivot = r+int(possible[-1])  # Opposite pivot selection from main code.
            A[[r, pivot]] = A[[pivot, r]]
            A[r] = A[r]*pow(int(A[r, col]), -1, prime) % prime
            factors = A[:, col].copy(); factors[r] = 0
            A = (A-factors[:, None]*A[r][None, :]) % prime
            pivots.append(col); r += 1
            if r == A.shape[0]: break
        return A, pivots
    j, t, q = array(Jac), array(T), array(Q)
    R, C = saved_c['jacobian_pivot_rows'], saved_c['jacobian_pivot_columns']
    square = j[np.ix_(R, C)]
    rr, pivots = reduced(np.concatenate([square, np.eye(35, dtype=np.int64)], axis=1))
    assert pivots == list(range(35))
    inverse = rr[:, 35:]; assert np.array_equal(square@inverse % prime, np.eye(35, dtype=np.int64))
    U = np.zeros((45, 35), dtype=np.int64); U[C, :] = inverse
    JU = j@U % prime
    W = -U@q[R, :] % prime
    quadratic_residual = (q-JU@q[R, :]) % prime
    assert np.array_equal(quadratic_residual, np.array(saved_c['residual_quadratic_matrix']))
    # Polarize every original product c*z_a*z_b as c*(T_a*W_b+T_b*W_a).
    cubic = np.zeros((135, 220), dtype=np.int64)
    for r, row in enumerate(terms):
        for coefficient, a, b in row:
            coefficient = mod(coefficient)
            cross = coefficient*(t[a, :, None]*W[b, None, :]
                                  +t[b, :, None]*W[a, None, :]) % prime
            for i in range(10):
                for j2, mon in enumerate(mons2):
                    cubic[r, index3[tuple(sorted((i,)+mon))]] += cross[i, j2]
        cubic[r] %= prime
    cubic_residual = (cubic-JU@cubic[R, :]) % prime
    assert np.array_equal(cubic_residual, np.array(saved_c['residual_cubic_matrix']))
    augmented, pivots = reduced(np.concatenate([quadratic_residual, cubic_residual], axis=1))
    assert sum(v < 55 for v in pivots) == 5
    quadrics = augmented[:5, :55]; cubics = augmented[5:len(pivots), 55:]
    assert np.array_equal(quadrics, np.array(saved_c['initial_quadrics']))
    assert np.array_equal(cubics, np.array(saved_c['initial_cubics']))
    assert len(cubics) == 21
    normal = saved_c['remaining_normal_variables']; params = saved_c['linear_parameters']
    assert len(normal) == 7 and len(params) == 3 and sorted(normal+params) == list(range(10))
    mon = {d: list(combinations_with_replacement(normal, d)) for d in [1, 2, 3, 4]}
    ind = {d: {m:i for i,m in enumerate(mon[d])} for d in mon}
    nq = quadrics[:, [index2[m] for m in mon[2]]]
    nc = cubics[:, [index3[m] for m in mon[3]]]
    def products(gens, degree, output_degree):
        answer = []
        for row in gens:
            for mult in mon[output_degree-degree]:
                out = [0]*len(mon[output_degree])
                for m, c in zip(mon[degree], row, strict=True):
                    out[ind[output_degree][tuple(sorted(m+mult))]] += int(c)
                answer.append([c % prime for c in out])
        return np.array(answer, dtype=np.int64)
    mat4 = np.concatenate([products(nq, 2, 4), products(nc, 3, 4)])
    assert np.array_equal(mat4, np.array(saved_c['quartic_matrix']))
    assert len(reduced(nq)[1]) == 5
    assert len(reduced(np.concatenate([products(nq, 2, 3), nc]))[1]) == 50
    M = [[int(v) for v in row] for row in mat4[saved_c['quartic_minor_rows'], :]]
    determinant = 1
    for col in range(210):
        pivot = next(i for i in range(col, 210) if M[i][col] % prime)
        if pivot != col: M[pivot], M[col] = M[col], M[pivot]; determinant *= -1
        p = M[col][col] % prime; determinant = determinant*p % prime
        for row in range(col+1, 210):
            multiplier = M[row][col]*pow(p, -1, prime) % prime
            for a in range(col, 210): M[row][a] = (M[row][a]-multiplier*M[col][a]) % prime
    assert determinant == saved_c['quartic_minor_determinant_mod_prime'] and determinant != 0
    zz, uu, beta = primitive_beta(model)
    beta = beta.subs(dict(zip(uu, [point[i] for i in one], strict=True)))
    a, b, d = s.symbols('a b d'); v1 = s.Matrix([a, 1, 0, b]); v2 = s.Matrix([-b, 0, 1, d])
    section = s.Matrix([(v1.T*beta.subs(dict(zip(zz, v, strict=True)))*v2)[0] for v in [v1, v2]])
    origin = {a: 0, b: 0, d: 0}; assert section.subs(origin) == s.zeros(2, 1)
    curve_jac = section.jacobian([a, b, d]).subs(origin); assert curve_jac.rank() == 2
    paths = [Path(__file__), quadratic_path, cubic_path, extra_path,
             ROOT/'code/s61_model.py', ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
             model.source/'symplectic_slice_reduced_equations.json']
    return {'id': 'S61_BINARY_LOCAL_DIMENSION_REVIEW_20261101', 'all_passed': True,
        'generated_utc': datetime.now(timezone.utc).isoformat(), 'prime': prime,
        'independent_135_specialization_and_jacobian': True, 'exact_QQ_quadratic_rank': 5,
        'sparse_product_cubic_residual_matches': True, 'reduced_cubic_rows': len(cubics),
        'quartic_matrix_shape': list(mat4.shape), 'quartic_minor_determinant': determinant,
        'independent_integer_determinant_check': True, 'canonical_curve_rank': 2,
        'local_dimension': 3, 'normal_colength_upper_bound': 65,
        'lifting_reason': 'QQ quadratic rank5 and unchanged mod-p rank imply all modular quadratic cancellations lift exactly; later pivots and quartic determinant are units modulo p.',
        'scope': 'Only the 135-equation unsaturated local ring at the recorded rational binary point; no t-saturation or completed-ring determination.',
        'proof_note': 'research/notes/S61_BINARY_LOCAL_DIMENSION_REVIEW_20261101.md',
        'elapsed_seconds': round(time.monotonic()-started, 3),
        'input_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); report = verify(); args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out: json.dump(report, out, indent=2); out.write('\n')
    print(json.dumps({'all_passed': True, 'local_dimension': 3, 'seconds': report['elapsed_seconds']}))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Exact independent audit of the perspective quotient on linear bivectors.

No legacy implementation is imported.  This verifies linear algebra over QQ and
polynomial identities, not extension of every first-order jet to a global
Poisson structure or decomposition of S61.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import sympy as sp


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    y = sp.Matrix(sp.symbols('y1:5'))
    a = sp.Matrix(sp.symbols('a1:5'))
    pairs = list(itertools.combinations(range(4), 2))
    triples = list(itertools.combinations(range(4), 3))
    basis = [(i, j, k) for i, j in pairs for k in range(4)]
    coeff = sp.symbols('l0:24')
    J = sp.diag(sp.Matrix([[0, 1], [-1, 0]]),
                sp.Matrix([[0, 1], [-1, 0]]))
    L = sp.zeros(4)
    for symbol, (i, j, k) in zip(coeff, basis):
        L[i, j] += symbol * y[k]
        L[j, i] -= symbol * y[k]

    def wedge(first: sp.Matrix, second: sp.Matrix) -> sp.Matrix:
        return first * second.T - second * first.T

    def trace(poly: sp.Matrix) -> sp.Expr:
        return sp.expand(poly[0, 1] + poly[2, 3])

    def vector(poly: sp.Matrix) -> sp.Matrix:
        return sp.Matrix([sp.expand(poly[i, j]).coeff(y[k])
                          for i, j, k in basis])

    def differential(poly: sp.Matrix) -> sp.Matrix:
        return sp.Matrix([
            sp.expand(sum(J[i, h] * sp.diff(poly[j, k], y[h])
                          + J[j, h] * sp.diff(poly[k, i], y[h])
                          + J[k, h] * sp.diff(poly[i, j], y[h])
                          for h in range(4)))
            for i, j, k in triples
        ])

    f = (a.T * y)[0]
    T = -2 * f * J + wedge(y, J * a)
    D = differential(L).jacobian(coeff)
    R = sp.Matrix([trace(L).coeff(yy) for yy in y]).jacobian(coeff)
    A = vector(T).jacobian(a)
    projection = sp.eye(24) + A * R / 5
    closed_basis = sp.Matrix.hstack(*D.nullspace())
    primitive_basis = sp.Matrix.vstack(D, R).nullspace()

    # Dg_a at x=g_a^{-1}(y); no numerical point is substituted.
    h = 1 - f
    derivative = h * (sp.eye(4) - y * a.T)
    exact_constant_pushforward = derivative * J * derivative.T
    expected = h**2 * (J + wedge(y, J*a))
    generic = sp.zeros(4)
    for symbol, (i, j) in zip(sp.symbols('p0:6'), pairs):
        generic[i, j] = symbol
        generic[j, i] = -symbol
    linear_map = sp.eye(4) - y*a.T
    b = sp.Matrix(sp.symbols('b1:5'))
    radial = wedge(y, b)
    checks = {
        'closed_linear_rank_4': D.rank() == 4,
        'trace_rank_4': R.rank() == 4,
        'closed_and_primitive_rank_8': sp.Matrix.vstack(D, R).rank() == 8,
        'perspective_rank_4': A.rank() == 4,
        'perspective_is_closed': D*A == sp.zeros(4),
        'trace_of_pushforward_direction_minus_5a': R*A == -5*sp.eye(4),
        'primitive_projection_idempotent': projection**2 == projection,
        'primitive_projection_preserves_differential': D*projection == D,
        'primitive_projection_kills_trace': R*projection == sp.zeros(4, 24),
        'closed_dimension_20': closed_basis.cols == 20,
        'primitive_closed_dimension_16': len(primitive_basis) == 16,
        'projected_closed_rank_16': (projection*closed_basis).rank() == 16,
        'exact_finite_constant_pushforward': all(
            sp.expand(v) == 0 for v in exact_constant_pushforward - expected),
        'linear_trace_transformation': sp.expand(trace(L+T)-trace(L)+5*f) == 0,
        'generic_skew_pushforward_identity': all(sp.expand(v) == 0 for v in
            linear_map*generic*linear_map.T-generic-wedge(y, generic*a)),
        'radial_cubic_cancellation': all(sp.expand(v) == 0 for v in
            linear_map*radial*linear_map.T-h*radial),
        'matrix_trace_formula': sp.expand(sp.trace(J.inv()*L)/2-trace(L)) == 0,
    }
    report = {
        'date': '2026-10-30',
        'field': 'QQ',
        'scope': 'Exact linear-jet quotient and perspective identities; no component/exhaustion claim.',
        'basis': [[i+1, j+1, k+1] for i, j, k in basis],
        'basis_meaning': 'y_k partial_i wedge partial_j; all listed indices are 1-based',
        'jacobi_linear_matrix': [list(map(int, D.row(i))) for i in range(D.rows)],
        'trace_matrix': [list(map(int, R.row(i))) for i in range(R.rows)],
        'perspective_pushforward_matrix': [list(map(int, A.row(i))) for i in range(A.rows)],
        'checks': checks,
        'all_passed': all(checks.values()),
        'sympy_version': sp.__version__,
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        if args.output.exists():
            raise FileExistsError(f'Refusing to overwrite {args.output}')
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'all_passed': report['all_passed'], 'checks': len(checks),
                      'closed_dimension': closed_basis.cols,
                      'primitive_dimension': len(primitive_basis)}, indent=2))
    if not report['all_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Retain all 90 affine coefficients in the five rank-two nilpotent charts.

Only linear reduction is performed. Every residual Jacobi coefficient and
the actual Pfaffian are saved; no sampled point or unbounded ideal solve.
"""
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s
from classification_multiple_hyperplane_nilpotent import (
    PAIRS, TRIPLES, coefficients, divergence, jacobi, pf,
)


def chart(kind, constant):
    started = time.monotonic()
    xx = s.symbols('x1:5')
    x = s.Matrix(xx)
    A = s.zeros(4)
    if kind == '31':
        A[0, 1] = A[1, 2] = 1
    else:
        A[0, 1] = A[2, 3] = 1
    X = A*x
    c = s.zeros(4, 1)
    if constant is not None:
        c[constant] = 1
    Z = X+c
    mon = [s.Integer(1), *xx, *(xx[i]*xx[j] for i in range(4) for j in range(i, 4))]
    unknowns = s.symbols('u0:90')
    P, Q = s.zeros(4), s.zeros(4)
    for k, (i, j) in enumerate(PAIRS):
        P[i, j] = sum(unknowns[15*k+l]*v for l, v in enumerate(mon))
        Q[i, j] = sum(unknowns[15*k+l]*v for l, v in enumerate(mon) if l >= 5)
        P[j, i], Q[j, i] = -P[i, j], -Q[i, j]
    inv = [sum(Z[k]*s.diff(P[i, j], xx[k])-P[k, j]*s.diff(Z[i], xx[k])
               -P[i, k]*s.diff(Z[j], xx[k]) for k in range(4)) for i, j in PAIRS]
    wedge = [X[i]*Q[j, k]-X[j]*Q[i, k]+X[k]*Q[i, j] for i, j, k in TRIPLES]
    equations = coefficients(inv+wedge+list(divergence(P, xx)-Z), xx)
    mat, rhs = s.linear_eq_to_matrix(equations, unknowns)
    solset = s.linsolve((mat, rhs), unknowns)
    if not solset:
        return {'kind': kind, 'constant_index': constant, 'linear_inconsistent': True}
    sol = next(iter(solset))
    free = sorted(set().union(*(v.free_symbols for v in sol)), key=str)
    nn = s.symbols(f'n0:{len(free)}')
    rename = dict(zip(free, nn))
    sub = dict(zip(unknowns, [s.expand(v.subs(rename)) for v in sol]))
    Pn = P.subs(sub, simultaneous=True).applyfunc(s.expand)
    Qn = Q.subs(sub, simultaneous=True).applyfunc(s.expand)
    assert all(s.expand(f.subs(sub, simultaneous=True)) == 0 for f in equations)
    assert pf(Qn) == 0
    assert jacobi(Qn, xx) == [0]*4
    residual = []
    for f in coefficients(jacobi(Pn, xx), xx):
        if f == 0:
            continue
        poly = s.Poly(f, *nn).monic().as_expr()
        f = s.factor(poly)
        if f not in residual:
            residual.append(f)
    return {
        'kind': kind, 'constant_index': constant,
        'full_initial_coefficient_count': 90,
        'linear_matrix_shape': list(mat.shape),
        'linear_matrix_sha256': hashlib.sha256(str(mat).encode()).hexdigest(),
        'full_free_parameter_count': len(free),
        'free_parameters_before_renaming': list(map(str, free)),
        'full_affine_substitution': {str(k): str(v) for k, v in sub.items()},
        'full_six_brackets': [str(Pn[i, j]) for i, j in PAIRS],
        'all_distinct_jacobi_coefficients': list(map(str, residual)),
        'actual_pfaffian': str(s.factor(pf(Pn))),
        'highest_pfaffian_zero': True, 'highest_jacobi_zero': True,
        'linear_substitution_verified': True,
        'seconds': round(time.monotonic()-started, 4),
    }


if __name__ == '__main__':
    cases = [chart(k, c) for k, c in [('31', None), ('31', 2), ('31', 3), ('22', None), ('22', 1)]]
    source = Path(__file__)
    helper = source.with_name('classification_multiple_hyperplane_nilpotent.py')
    report = {
        'status': 'PASS', 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'helper_sha256': hashlib.sha256(helper.read_bytes()).hexdigest(),
        'coefficient_ring': 'QQ[n0,...,x1,x2,x3,x4], all five affine normal forms',
        'scope': 'Full linear normal forms and all residual equations, not support coverage.',
        'cases': cases,
    }
    out = Path('reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json')
    out.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', 'report': str(out), 'case_counts': [
        [c['kind'], c['constant_index'], c.get('full_free_parameter_count'),
         len(c.get('all_distinct_jacobi_coefficients', [])), c.get('seconds')]
        for c in cases]}))

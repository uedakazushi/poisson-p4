#!/usr/bin/env python3
"""Exact identities supporting the global structural reductions.

This is NOT a verifier of the 147-component classification. It checks universal
first-jet identities, deformation expansions, and a rank-two counterexample.
Python 3 and SymPy suffice; no repository producers or sampled points are used.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import platform
from pathlib import Path
import sympy as s

PAIRS = list(itertools.combinations(range(4), 2))
TRIPLES = list(itertools.combinations(range(4), 3))


def wedge(left: dict, right: dict) -> dict:
    result: dict = {}
    for a, va in left.items():
        for b, vb in right.items():
            word = a + b
            if len(set(word)) != len(word):
                continue
            sign = (-1) ** sum(word[i] > word[j]
                               for i in range(len(word)) for j in range(i + 1, len(word)))
            key = tuple(sorted(word))
            result[key] = s.expand(result.get(key, 0) + sign * va * vb)
    return result


def skew(values) -> s.Matrix:
    matrix = s.zeros(4)
    for (i, j), value in zip(PAIRS, values):
        matrix[i, j], matrix[j, i] = value, -value
    return matrix


def run() -> dict:
    counts: dict[str, int] = {}

    def check(name: str, expressions) -> None:
        items = list(expressions)
        for i, value in enumerate(items):
            if s.expand(value) != 0:
                raise AssertionError(f'{name}[{i}]: {s.expand(value)}')
        counts[name] = len(items)

    p = s.symbols('p01 p02 p03 p12 p13 p23')
    dp = {(i, j, k): s.Symbol(f'd{i}{j}_{k}') for i, j in PAIRS for k in range(4)}
    P = skew(p)

    def derivative(value, k):
        return sum(s.diff(value, a) * dp[i, j, k] for (i, j), a in zip(PAIRS, p))

    f = P[0, 1] * P[2, 3] - P[0, 2] * P[1, 3] + P[0, 3] * P[1, 2]
    Z = s.Matrix([sum(derivative(P[i, j], i) for i in range(4)) for j in range(4)])
    J = {(a, b, c): s.expand(sum(
        P[a, l] * derivative(P[b, c], l)
        + P[b, l] * derivative(P[c, a], l)
        + P[c, l] * derivative(P[a, b], l) for l in range(4))) for a, b, c in TRIPLES}
    theta = s.Matrix([(-1) ** (3 - i) * J[tuple(j for j in range(4) if j != i)]
                      for i in range(4)])
    residual = P.T * s.Matrix([derivative(f, i) for i in range(4)]) - f * Z
    check('universal_modular_jacobi', residual - P.T * theta)
    if all(s.expand(v) == 0 for v in residual + P.T * theta):
        raise AssertionError('Wrong-sign mutation was not detected')
    counts['wrong_sign_detected'] = 1

    # alpha = contraction of P with dx0 ^ dx1 ^ dx2 ^ dx3.
    alpha = {}
    for (i, j), value in zip(PAIRS, p):
        complement = tuple(k for k in range(4) if k not in (i, j))
        word = (i, j) + complement
        sign = (-1) ** sum(word[a] > word[b] for a in range(4) for b in range(a + 1, 4))
        alpha[complement] = sign * value
    dalpha = {}
    for pair, value in alpha.items():
        for k in range(4):
            term = wedge({(k,): derivative(value, k)}, {pair: 1})
            for index, coefficient in term.items():
                dalpha[index] = dalpha.get(index, 0) + coefficient
    contraction_Z = {tuple(k for k in range(4) if k != j): (-1)**j * Z[j]
                     for j in range(4)}
    check('divergence_exterior_identity', [dalpha.get(k, 0) + contraction_Z[k] for k in TRIPLES])
    g = s.symbols('g0:4')
    dg = {(i,): g[i] for i in range(4)}
    field = P.T * s.Matrix(g)
    contracted = {tuple(k for k in range(4) if k != j): (-1)**j * field[j]
                  for j in range(4)}
    product = wedge(dg, alpha)
    check('gradient_exterior_identity', [product.get(k, 0) + contracted[k] for k in TRIPLES])

    t = s.Symbol('t')
    av, bv, cv = [s.symbols(prefix + '0:6') for prefix in ('a', 'b', 'c')]
    a, b, c = [dict(zip(PAIRS, values)) for values in (av, bv, cv)]
    form_t = {k: a[k] + t*b[k] + t*t*c[k] for k in PAIRS}
    q_t = s.expand(wedge(form_t, form_t)[(0, 1, 2, 3)] / 2)
    check('pfaffian_lifting_expansion', [q_t.coeff(t, 1) - wedge(a, b)[(0, 1, 2, 3)],
        q_t.coeff(t, 2) - wedge(a, c)[(0, 1, 2, 3)] - wedge(b, b)[(0, 1, 2, 3)] / 2])
    F, phi, psi = s.symbols('F phi psi')
    grad_F, grad_phi, grad_psi = [s.symbols(prefix + '0:4') for prefix in ('F', 'phi', 'psi')]
    da, db, dc = [dict(zip(TRIPLES, s.symbols(prefix + '0:4'))) for prefix in ('da', 'db', 'dc')]
    Ft = F + t*phi + t*t*psi
    dFt = {(i,): grad_F[i] + t*grad_phi[i] + t*t*grad_psi[i] for i in range(4)}
    exterior_product = wedge(dFt, form_t)
    linear_expected = {k: F*db[k]+phi*da[k]-wedge(dict(zip([(i,) for i in range(4)], grad_F)), b).get(k,0)
        -wedge(dict(zip([(i,) for i in range(4)], grad_phi)), a).get(k,0) for k in TRIPLES}
    quadratic_expected = {k: F*dc[k]+phi*db[k]+psi*da[k]
        -wedge(dict(zip([(i,) for i in range(4)], grad_F)), c).get(k,0)
        -wedge(dict(zip([(i,) for i in range(4)], grad_phi)), b).get(k,0)
        -wedge(dict(zip([(i,) for i in range(4)], grad_psi)), a).get(k,0) for k in TRIPLES}
    expansion = {k: s.expand(Ft*(da[k]+t*db[k]+t*t*dc[k])-exterior_product.get(k,0)) for k in TRIPLES}
    check('closed_form_first_order', [expansion[k].coeff(t,1)-linear_expected[k] for k in TRIPLES])
    check('closed_form_second_order', [expansion[k].coeff(t,2)-quadratic_expected[k] for k in TRIPLES])

    # Non-Poisson rank-two bivector with vanishing modular-Pfaffian residual.
    x = s.symbols('x0:4')
    example = skew([1, x[0], 0, 0, 0, 0])
    example_f = example[0,1]*example[2,3]-example[0,2]*example[1,3]+example[0,3]*example[1,2]
    example_Z = s.Matrix([sum(s.diff(example[i,j],x[i]) for i in range(4)) for j in range(4)])
    example_J = sum(example[0,l]*s.diff(example[1,2],x[l])+example[1,l]*s.diff(example[2,0],x[l])
                    +example[2,l]*s.diff(example[0,1],x[l]) for l in range(4))
    check('rank_two_gate_counterexample', [example_f, example_J-1]
          +list(example.T*s.Matrix([s.diff(example_f,v) for v in x])-example_f*example_Z))

    return {'date':'2026-11-09','version':'v1','status':'PASS',
        'scope':'Universal polynomial identities and the necessity of the nonzero Pfaffian hypothesis only.',
        'not_verified':['147-component classification','All source recognition and noncontainment',
                        'Global generic-fibre flatness','All higher deformation obstructions'],
        'base_commit':'d970bf8dda43f08ac4d54b27c088fdf7beebbf4d',
        'python':platform.python_version(),'sympy':s.__version__,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'checks':counts,'check_count':sum(counts.values()),
        'sampling_used':False,'imports_old_producers':False}


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=run()
    text=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text,encoding='utf-8')
    print(text)


if __name__=='__main__':
    main()

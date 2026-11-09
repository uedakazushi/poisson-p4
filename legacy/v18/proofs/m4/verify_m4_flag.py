#!/usr/bin/env python3
"""Exact, standard-library verification of the exceptional m4-pencil exclusion.

This verifies finite polynomial identities and a full-column-rank minor.  It
does not certify the geometric hypotheses, completeness of the Poisson
classification, or existence of a full higher-order lift.
"""
from fractions import Fraction
from pathlib import Path
import json


def add(p, q):
    out = dict(p)
    for e, c in q.items():
        out[e] = out.get(e, 0) + c
        if not out[e]:
            del out[e]
    return out


def mul(p, q):
    out = {}
    for e, c in p.items():
        for f, d in q.items():
            g = tuple(a + b for a, b in zip(e, f))
            out[g] = out.get(g, 0) + c * d
    return {e: c for e, c in out.items() if c}


def scale(p, c):
    return {e: c * d for e, d in p.items() if c * d}


def derivative(p, i):
    out = {}
    for e, c in p.items():
        if e[i]:
            f = list(e)
            f[i] -= 1
            out[tuple(f)] = c * e[i]
    return out


def variable(n, i):
    e = [0] * n
    e[i] = 1
    return {tuple(e): 1}


def power(p, n):
    out = {(0,) * len(next(iter(p))): 1}
    for _ in range(n):
        out = mul(out, p)
    return out


def determinant(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    n = len(a)
    assert all(len(row) == n for row in a)
    answer = Fraction(1)
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != j:
            a[pivot], a[j] = a[j], a[pivot]
            answer = -answer
        p = a[j][j]
        answer *= p
        for i in range(j + 1, n):
            q = a[i][j] / p
            for k in range(j, n):
                a[i][k] -= q * a[j][k]
    return answer


def matmul(a, b):
    out = [[{} for _ in b[0]] for _ in a]
    for i in range(len(a)):
        for j in range(len(b[0])):
            for k in range(len(b)):
                out[i][j] = add(out[i][j], mul(a[i][k], b[k][j]))
    return out


def main():
    cert = json.loads(Path(__file__).with_name('m4_invariant_certificate.json').read_text())
    s, t, x, y, z = (variable(5, i) for i in range(5))
    ds = lambda f: add(mul(y, derivative(f, 2)), mul(z, derivative(f, 3)))
    dt = lambda f: add(mul(x, derivative(f, 3)), scale(mul(y, derivative(f, 4)), -1))
    d = lambda f: add(mul(s, ds(f)), mul(t, dt(f)))
    normal_mons = [(i, j, 4 - i - j) for i in range(5) for j in range(5 - i)]
    exponents = [(i, 1-i, *e) for i in (1, 0) for e in normal_mons]
    assert exponents == [tuple(e) for e in cert['input_monomials']]
    images = [d({e: 1}) for e in exponents]
    keys = sorted(set().union(*(set(p) for p in images)))
    assert keys == [tuple(e) for e in cert['output_monomials']]
    matrix = [[p.get(e, 0) for p in images] for e in keys]
    assert [len(matrix), len(matrix[0])] == cert['matrix_shape'] == [43, 30]
    minor_det = determinant([matrix[i] for i in cert['selected_rows']])
    assert minor_det == cert['minor_determinant'] == 110075314176

    q = add(scale(mul(x, z), 2), scale(power(y, 2), -1))
    r = add(scale(mul(x, z), 2), power(y, 2))
    avec = [power(z, 4), mul(power(z, 2), q), power(q, 2)]
    bvec = [power(x, 4), mul(power(x, 2), r), power(r, 2)]
    assert all(not ds(f) for f in avec)
    assert all(not dt(f) for f in bvec)
    middle = [dt(f) for f in avec] + [ds(f) for f in bvec]
    selected_mons = [mul(y, f) for f in [power(z, 3), power(x, 3),
                    mul(x, power(z, 2)), mul(power(y, 2), z),
                    mul(power(x, 2), z), mul(x, power(y, 2))]]
    selected_exps = [next(iter(f)) for f in selected_mons]
    compact = [[p.get(e, 0) for p in middle] for e in selected_exps]
    assert compact == cert['compact_matrix']
    compact_det = determinant(compact)
    assert compact_det == cert['compact_determinant'] == -147456

    # Clear denominators in G K G^{-1} = d K0.  Here H=d^2 G,
    # Ktilde=d^2 K, so H*Ktilde=d^3*K0*H, as polynomials in a,b,d.
    a, b, dd = (variable(3, i) for i in range(3))
    zero = {}
    d2, d3 = power(dd, 2), power(dd, 3)
    ad = mul(a, dd)
    h = [[d2, scale(ad, -1), add(mul(b, dd), scale(power(a, 2), 2))],
         [zero, d2, scale(ad, -1)], [zero, zero, d2]]
    kt = [[mul(a, d2), mul(b, d2), power(a, 3)],
          [d3, scale(mul(a, d2), -2), add(mul(b, d2), scale(mul(power(a, 2), dd), 3))],
          [zero, scale(d3, -1), mul(a, d2)]]
    d3k0 = [[zero, zero, zero], [d3, zero, zero], [zero, scale(d3, -1), zero]]
    assert matmul(h, kt) == matmul(d3k0, h)
    print(json.dumps({'all_passed': True, 'full_matrix_shape': [43, 30],
                      'full_rank': 30, 'full_minor_determinant': int(minor_det),
                      'compact_minor_determinant': int(compact_det),
                      'exceptional_conjugation_identity': True}, indent=2))


if __name__ == '__main__':
    main()

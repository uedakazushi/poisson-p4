#!/usr/bin/env python3
"""A dominant six-dimensional Sp4 slice and its exact Fisher invariants.

The tangent determinant certifies dominance. Fisher contractions are performed
over QQ[s0,...,s5]; numerical examples only check the implementation/conventions.
"""
from __future__ import annotations

import argparse
import hashlib
from itertools import combinations, combinations_with_replacement, permutations, product
import json
import math
from pathlib import Path
import time

import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

from elliptic_curve_invariants_20261101 import curve_invariants_from_beta
from weight1_structure import constraints, infinitesimal_action, sp_basis

ROOT = Path(__file__).resolve().parents[1]
C5_SOURCE = ROOT / 'reports/TORIC_DEGREE5_COVARIANT_20261101.json'


def evaluate(f, values):
    return sum((coefficient * math.prod(values[i] ** power for i, power in enumerate(monomial)
                if power) for monomial, coefficient in f.items()), QQ.zero)


def serialize(f):
    return [[list(monomial), str(coefficient)] for monomial, coefficient in sorted(f.items())]


def dominant_slice():
    basis = s.Matrix.hstack(*s.Matrix.vstack(*constraints()).nullspace())
    rows = list(basis.T.rref()[1])
    left = basis[rows, :].inv()
    actions = []
    for a in sp_basis():
        images = s.Matrix.hstack(*(infinitesimal_action(basis[:, i], a) for i in range(16)))
        action = left * images[rows, :]
        assert basis * action == images
        actions.append(action)
    point = s.Matrix(range(1, 17))
    orbit = s.Matrix.hstack(*(a * point for a in actions))
    augmented = orbit.row_join(point)
    rows = list(augmented.T.rref()[1])
    assert len(rows) == 11
    free = [i for i in range(16) if i not in rows]
    assert free == [11, 12, 13, 14, 15]
    slice_matrix = s.Matrix.hstack(point, *(s.eye(16)[:, i] for i in free))
    tangent = orbit.row_join(slice_matrix)
    assert tangent.det() == 8428695240
    return basis, slice_matrix, tangent


def restricted_fisher(basis, slice_matrix):
    beta = basis * slice_matrix
    r, *variables = ring(','.join(f's{i}' for i in range(6)), QQ)
    rx, *x = ring('x0,x1,x2,x3,x4', r.to_domain())
    sv = [rx.from_dict({(0,) * 5: q}) for q in variables]
    a, b, c, d, e = x
    plucker = [a, b, c, d, e, -a]
    pairs4 = list(combinations(range(4), 2))
    star = [[rx.zero for _ in range(4)] for _ in range(4)]
    for (i, j), value in zip(pairs4, [-a, -e, d, c, -b, a], strict=True):
        star[i][j], star[j][i] = value, -value
    ell = [sum((sum((QQ.from_sympy(beta[4 * h + i, k]) * sv[k] for k in range(6)), rx.zero)
                     * plucker[h] for h in range(6)), rx.zero) for i in range(4)]
    phi = [row + [ell[i]] for i, row in enumerate(star)] + [[-t for t in ell] + [rx.zero]]

    def pfaff4(indices):
        i, j, k, l = indices
        return phi[i][j] * phi[k][l] - phi[i][k] * phi[j][l] + phi[i][l] * phi[j][k]

    pfaffians = [(-1) ** i * pfaff4([j for j in range(5) if j != i]) for i in range(5)]
    assert all(not sum((phi[i][j] * pfaffians[j] for j in range(5)), rx.zero) for i in range(5))
    dp = [[p.diff(t) for p in pfaffians] for t in x]
    da = [[[phi[i][j].diff(t) for j in range(5)] for i in range(5)] for t in x]

    def scalar(left, matrix, right):
        return sum((left[i] * matrix[i][j] * right[j] for i in range(5) for j in range(5)
                    if matrix[i][j]), rx.zero)

    omega = [[rx.zero for _ in range(5)] for _ in range(5)]
    for i, j in combinations(range(5), 2):
        k, l, m = [z for z in range(5) if z not in (i, j)]
        permutation = [i, j, k, l, m]
        sign = (-1) ** sum(permutation[a] > permutation[b] for a in range(5)
                          for b in range(a + 1, 5))
        value = sign * (scalar(dp[k], da[l], dp[m]) + scalar(dp[m], da[k], dp[l])
                        + scalar(dp[l], da[m], dp[k]))
        omega[i][j], omega[j][i] = value, -value
    rr = range(5)
    keys = list(product(rr, repeat=4))
    w = {}
    for i, j, a, b in keys:
        exponent = [0] * 5
        exponent[a] += 1
        exponent[b] += 1
        w[i, j, a, b] = (2 if a == b else 1) * omega[i][j].get(tuple(exponent), r.zero)
    z = {(i, j, a, b): sum((w[i, h, u, a] * w[j, u, h, b] + w[i, h, u, b] * w[j, u, h, a]
                            for h in rr for u in rr), r.zero) for i, j, a, b in keys}
    assert all(z[i, j, a, b] == z[j, i, a, b] == z[i, j, b, a] for i, j, a, b in keys)
    pairs = list(combinations_with_replacement(rr, 2))
    multiplicity = {ij: 1 if ij[0] == ij[1] else 2 for ij in pairs}
    c4 = sum((multiplicity[i, j] * multiplicity[a, b] * z[i, j, a, b] * z[a, b, i, j]
              for i, j in pairs for a, b in pairs), r.zero) / QQ(13440)

    # T_{I,A} is symmetric in its last triple A. Summing over the distinct
    # permutations of its first triple converts the full 125-by-125 contraction
    # into a 35-by-35 contraction, with no extra multinomial factors.
    triples = list(combinations_with_replacement(rr, 3))
    distinct_permutations = {triple: set(permutations(triple)) for triple in triples}
    t = {}
    for first in triples:
        for last in triples:
            a, b, c = last
            t[first, last] = sum((z[i, j, h, a] * w[h, k, b, c] + z[i, j, h, b] * w[h, k, a, c]
                                  + z[i, j, h, c] * w[h, k, a, b]
                                  for i, j, k in distinct_permutations[first] for h in rr), r.zero)
    c6 = -sum(((1 if i == j else 2) * t[first, last] * t[last, first]
               for i, first in enumerate(triples) for j, last in enumerate(triples) if i <= j),
              r.zero) / QQ(1036800)
    assert set(map(sum, c4)) == {8} and set(map(sum, c6)) == {12}
    return r, variables, c4, c6


def restrict_c5(record, r, variables, slice_matrix):
    substitutions = [sum((QQ.from_sympy(slice_matrix[i, j]) * variables[j] for j in range(6)), r.zero)
                     for i in range(16)]
    powers = [[value ** k for k in range(6)] for value in substitutions]
    answer = []
    for entry in record['C5_sparse']:
        value = r.zero
        for monomial, coefficient in entry:
            term = r(QQ.from_sympy(s.Rational(coefficient)))
            for i, power in enumerate(monomial):
                if power:
                    term *= powers[i][power]
            value += term
        assert set(map(sum, value)) == {5}
        answer.append(value)
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    basis, slice_matrix, tangent = dominant_slice()
    r, variables, c4, c6 = restricted_fisher(basis, slice_matrix)
    source = json.loads(C5_SOURCE.read_text())
    assert s.Matrix(source['weight_basis_24_by16']) == basis
    u = restrict_c5(source, r, variables, slice_matrix)
    tests = []
    for parameters in ((1, 0, 0, 0, 0, 0), (1, 2, 3, 4, 5, 6), (2, -1, 3, -2, 5, -4)):
        beta = basis * slice_matrix * s.Matrix(parameters)
        expected4, expected6, _, _ = curve_invariants_from_beta(beta)
        assert evaluate(c4, parameters) == expected4
        assert evaluate(c6, parameters) == expected6
        b_values = slice_matrix * s.Matrix(parameters)
        for entry, restricted in zip(source['C5_sparse'], u, strict=True):
            expected = sum(s.Rational(coefficient) * math.prod(b_values[i] ** power
                            for i, power in enumerate(monomial) if power)
                           for monomial, coefficient in entry)
            assert expected == evaluate(restricted, parameters)
        tests.append({'parameters': parameters, 'c4': str(expected4), 'c6': str(expected6)})
    report = {
        'all_passed': True, 'field': 'QQ[s0,s1,s2,s3,s4,s5]',
        'variables': list(map(str, variables)),
        'weight_basis_24_by16': [[str(x) for x in row] for row in basis.tolist()],
        'slice_matrix_b': [[str(x) for x in row] for row in slice_matrix.tolist()],
        'dominance_tangent_point_s': [1, 0, 0, 0, 0, 0],
        'dominance_tangent_matrix': [[str(x) for x in row] for row in tangent.tolist()],
        'dominance_tangent_determinant': str(tangent.det()),
        'dominance_statement': 'Sp4 times the six-dimensional linear slice maps dominantly to W16',
        'primary_formula': 'https://arxiv.org/html/1406.6953v1#S2',
        'primary_formula_location': 'Fisher-Sadek Theorem2.3, canonical Pfaffian Omega construction',
        'c4_sparse': serialize(c4), 'c6_sparse': serialize(c6),
        'C5_sparse': [serialize(f) for f in u],
        'monomial_counts': {'c4': len(c4), 'c6': len(c6), 'C5_total': sum(map(len, u))},
        'independent_original_contraction_sample_checks': tests,
        'scope': 'A dominance certificate and exact restrictions of the established Fisher invariants '
                 'and C5. An invariant/covariant polynomial identity verified on this slice extends '
                 'universally. No cover identification or S61 exhaustion is asserted here.',
        'source_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (Path(__file__), C5_SOURCE, ROOT / 'code/weight1_structure.py',
                         ROOT / 'code/elliptic_curve_invariants_20261101.py')},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps({'all_passed': True, 'tangent_determinant': '8428695240',
                      'monomial_counts': report['monomial_counts'],
                      'seconds': time.time() - started, 'output': str(args.output)}, indent=2))


if __name__ == '__main__':
    main()

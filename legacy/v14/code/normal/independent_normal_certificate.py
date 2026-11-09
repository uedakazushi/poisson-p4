#!/usr/bin/env python3
"""Independent exact arithmetic certificate, using the Python standard library.

This implementation does not import the main SymPy verifier.  A modular
minor, together with an exact nonzero kernel vector, proves the characteristic
zero rank for one explicitly specified integer matrix.  No classification is
inferred from sampling over a finite field.
"""
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import json

N = 5
ZERO = (0,) * N


def clean(p):
    return {m: c for m, c in p.items() if c}


def add(*ps):
    out = {}
    for p in ps:
        for m, c in p.items():
            out[m] = out.get(m, 0) + c
    return clean(out)


def scale(p, a):
    return clean({m: c * a for m, c in p.items()})


def mul(p, q):
    out = {}
    for m, a in p.items():
        for n, b in q.items():
            e = tuple(u + v for u, v in zip(m, n))
            out[e] = out.get(e, 0) + a * b
    return clean(out)


def diff(p, i):
    out = {}
    for m, c in p.items():
        if m[i]:
            e = list(m)
            e[i] -= 1
            out[tuple(e)] = c * m[i]
    return out


def mons(n, degree):
    if n == 1:
        return [(degree,)]
    return [(i,) + tail for i in range(degree, -1, -1)
            for tail in mons(n - 1, degree - i)]


def rank_fraction(rows):
    basis = {}
    for row in rows:
        v = {j: Fraction(a) for j, a in enumerate(row) if a}
        while v:
            pivot = min(v)
            if pivot not in basis:
                a = v[pivot]
                basis[pivot] = {j: b / a for j, b in v.items()}
                break
            a = v[pivot]
            for j, b in basis[pivot].items():
                value = v.get(j, 0) - a * b
                if value:
                    v[j] = value
                else:
                    v.pop(j, None)
    return len(basis)


def modular_rank(rows, prime):
    basis = {}
    chosen_rows = []
    chosen_cols = []
    for i, row in enumerate(rows):
        v = {j: int(a) % prime for j, a in enumerate(row) if int(a) % prime}
        while v:
            pivot = min(v)
            if pivot not in basis:
                inv = pow(v[pivot], -1, prime)
                basis[pivot] = {j: b * inv % prime for j, b in v.items()}
                chosen_rows.append(i)
                chosen_cols.append(pivot)
                break
            a = v[pivot]
            for j, b in basis[pivot].items():
                value = (v.get(j, 0) - a * b) % prime
                if value:
                    v[j] = value
                else:
                    v.pop(j, None)
    return len(basis), chosen_rows, chosen_cols


def determinant_mod(matrix, prime):
    a = [[int(v) % prime for v in row] for row in matrix]
    value = 1
    for c in range(len(a)):
        pivot = next((r for r in range(c, len(a)) if a[r][c]), None)
        if pivot is None:
            return 0
        if pivot != c:
            a[c], a[pivot] = a[pivot], a[c]
            value = -value
        value = value * a[c][c] % prime
        inv = pow(a[c][c], -1, prime)
        for r in range(c + 1, len(a)):
            factor = a[r][c] * inv % prime
            if factor:
                for j in range(c, len(a)):
                    a[r][j] = (a[r][j] - factor * a[c][j]) % prime
    return value % prime


def main():
    x = []
    for i in range(N):
        e = [0] * N
        e[i] = 1
        x.append({tuple(e): 1})
    p = [[{} for _ in range(N)] for _ in range(N)]
    for i in range(N):
        j = (i + 1) % N
        p[i][j] = add(scale(mul(x[i], x[j]), -2),
                      scale(mul(x[(i + 2) % N], x[(i + 4) % N]), -10),
                      scale(mul(x[(i + 3) % N], x[(i + 3) % N]), 5))
        p[j][i] = scale(p[i][j], -1)
        j = (i + 2) % N
        p[i][j] = add(scale(mul(x[i], x[j]), -4),
                      scale(mul(x[(i + 3) % N], x[(i + 4) % N]), 10),
                      scale(mul(x[(i + 1) % N], x[(i + 1) % N]), 5))
        p[j][i] = scale(p[i][j], -1)
    pf = []
    for i in range(N):
        a, b, c, d = [j for j in range(N) if j != i]
        pf.append(scale(add(mul(p[a][b], p[c][d]),
                            scale(mul(p[a][c], p[b][d]), -1),
                            mul(p[a][d], p[b][c])), (-1) ** i))
    F = scale(add(*(mul(x[i], pf[i]) for i in range(N))), Fraction(1, 5))
    assert all(diff(F, i) == pf[i] for i in range(N))
    assert all(not add(*(diff(p[i][j], i) for i in range(N))) for j in range(N))
    for i, j, k in combinations(range(N), 3):
        assert not add(*(add(mul(p[i][l], diff(p[j][k], l)),
                             mul(p[j][l], diff(p[k][i], l)),
                             mul(p[k][l], diff(p[i][j], l))) for l in range(N)))
    assert all(not add(*(mul(p[i][j], pf[j]) for j in range(N))) for i in range(N))

    # Homogeneous Macaulay matrices certify a finite three-variable slice.
    partials = [{m[:3]: c for m, c in g.items() if m[3:] == (0, 0)} for g in pf]
    hilbert = []
    ranks = []
    matrix_shapes = []
    for degree in range(9):
        output_mons = mons(3, degree)
        column_specs = [(i, m) for i in range(5) for m in mons(3, degree - 4)] if degree >= 4 else []
        rows = []
        for e in output_mons:
            row = []
            for i, m in column_specs:
                q = tuple(a - b for a, b in zip(e, m))
                row.append(partials[i].get(q, 0) if min(q) >= 0 else 0)
            rows.append(row)
        rank = rank_fraction(rows)
        ranks.append(rank)
        hilbert.append(len(output_mons) - rank)
        matrix_shapes.append([len(output_mons), len(column_specs)])
    assert hilbert == [1, 3, 6, 10, 10, 6, 3, 1, 0]
    assert ranks[8] == 45

    # Direct 1050 x 150 matrix in a separately generated monomial convention.
    pairs = list(combinations(range(N), 2))
    quad_mons = mons(N, 2)
    sextic_mons = mons(N, 6)
    row_index = {(i, m): 210 * i + j for i in range(N) for j, m in enumerate(sextic_mons)}
    matrix = [[0] * 150 for _ in range(1050)]
    vector = []
    for pair_index, (i, j) in enumerate(pairs):
        for mon_index, m in enumerate(quad_mons):
            col = 15 * pair_index + mon_index
            vector.append(p[i][j].get(m, 0))
            for target, g, sign in [(i, pf[j], 1), (j, pf[i], -1)]:
                for n, c in g.items():
                    e = tuple(a + b for a, b in zip(m, n))
                    matrix[row_index[target, e]][col] += sign * c
    assert any(vector)
    assert all(sum(a * b for a, b in zip(row, vector)) == 0 for row in matrix)
    prime = 1000003
    rank, rids, cids = modular_rank(matrix, prime)
    assert rank == 149
    minor = [[matrix[r][c] for c in cids] for r in rids]
    determinant = determinant_mod(minor, prime)
    assert determinant != 0
    out = {
        "implementation": "Python standard library; Fraction and integer modular arithmetic",
        "exact_Pfaffian_gradient_identity": True,
        "exact_divergence_zero": True,
        "exact_Jacobi_zero": True,
        "exact_P_gradient_annihilation": True,
        "normality_certificate": {
            "slice": "x3=x4=0",
            "Macaulay_matrix_shapes_degrees_0_to_8": matrix_shapes,
            "ranks_over_Q_degrees_0_to_8": ranks,
            "Hilbert_function_degrees_0_to_8": hilbert,
            "slice_length": sum(hilbert),
            "contains_eighth_power_of_slice_maximal_ideal": True,
            "conclusion": "height J(F)>=3; Pfaffian bound gives height exactly3; F is irreducible and normal"
        },
        "annihilator_rank_certificate": {
            "shape": [1050, 150], "prime": prime,
            "rank_mod_prime": rank,
            "nonzero_minor_rows": rids,
            "nonzero_minor_columns": cids,
            "minor_determinant_mod_prime": determinant,
            "exact_nonzero_kernel_vector": vector,
            "rank_over_Q": 149,
            "meaning": "lower bound from a nonzero integer minor; upper bound from an exact kernel vector"
        },
        "quintic_terms": [[list(m), int(c)] for m, c in sorted(F.items(), reverse=True)],
        "limitations": ["one witness, not an exhaustive classification", "Buchsbaum-Eisenbud theorem is a theoretical input"]
    }
    target = Path(__file__).with_suffix('.json')
    target.write_text(json.dumps(out, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({"passed": True, "slice_length": sum(hilbert),
                      "rank_over_Q": 149, "minor_mod_prime": determinant,
                      "output": target.name}, ensure_ascii=False))


if __name__ == '__main__':
    main()

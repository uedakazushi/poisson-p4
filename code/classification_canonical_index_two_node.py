#!/usr/bin/env python3
"""One bounded whole-canonical-pair ambient restriction lift; QQ sparse arithmetic."""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import signal
import time


def add(*terms):
    out = {}
    for p in terms:
        for m, c in p.items():
            out[m] = out.get(m, 0) + c
            if not out[m]:
                del out[m]
    return out


def scale(p, c):
    return {m: v*c for m, v in p.items() if v*c}


def mul(p, q):
    out = {}
    for a, c in p.items():
        for b, d in q.items():
            m = tuple(x+y for x, y in zip(a, b))
            out[m] = out.get(m, 0) + c*d
    return {m: c for m, c in out.items() if c}


def power(p, n):
    out = {(0,)*len(next(iter(p))): 1}
    for _ in range(n):
        out = mul(out, p)
    return out


def var(n, j):
    return {tuple(int(i == j) for i in range(n)): 1}


def diff(p, j):
    out = {}
    for m, c in p.items():
        if m[j]:
            a = list(m)
            a[j] -= 1
            out[tuple(a)] = c*m[j]
    return out


def mons(n, d):
    return sorted([tuple(a.count(i) for i in range(n))
                   for a in it.combinations_with_replacement(range(n), d)], reverse=True)


def serial(p):
    return [[list(m), str(c)] for m, c in sorted(p.items(), reverse=True)]


def rem(p, f):
    """Ordinary division by ONE polynomial, graded lex, constant leading coefficient."""
    p = dict(p)
    lm = max(f, key=lambda m: (sum(m), m))
    lc = f[lm]
    out = {}
    while p:
        m = max(p, key=lambda a: (sum(a), a))
        c = p.pop(m)
        if all(a >= b for a, b in zip(m, lm)):
            shift = tuple(a-b for a, b in zip(m, lm))
            factor = Q(c, lc)
            for a, v in f.items():
                if a == lm:
                    continue
                b = tuple(x+y for x, y in zip(shift, a))
                p[b] = p.get(b, 0)-factor*v
                if not p[b]:
                    del p[b]
        else:
            out[m] = c
    return out


def raw(p, i, j):
    if i == j:
        return {}
    return p.get((i, j), {}) if i < j else scale(p.get((j, i), {}), -1)


def dehom(p):
    out = {}
    for m, c in p.items():
        out[m[:4]] = out.get(m[:4], 0)+c
    return {m: c for m, c in out.items() if c}


def axpy(a, b, c):
    for k, v in b.items():
        a[k] = a.get(k, 0)+c*v
        if not a[k]:
            del a[k]


def canonical():
    a, b, c, d, e = [var(5, j) for j in range(5)]
    q = add(mul(a, c), scale(power(b, 2), -1))
    A = add(power(c, 2), mul(b, d))
    B = add(power(a, 2), mul(b, e))
    C = add(mul(power(c, 2), e), mul(power(a, 2), d),
            scale(mul(mul(a, c), b), 2), mul(mul(b, d), e), scale(power(b, 3), -1))
    assert mul(b, C) == add(mul(A, B), scale(power(q, 2), -1))
    F = add(mul(C, q), scale(mul(c, power(A, 2)), -1),
            scale(mul(b, mul(A, B)), 2), scale(mul(a, power(B, 2)), -1))
    g1 = add(mul(b, B), scale(mul(c, A), -1))
    g2 = add(mul(b, A), scale(mul(a, B), -1))
    return F, g1, g2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    resource.setrlimit(resource.RLIMIT_AS, (768*1024**2,)*2)
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError('120s finite contract')))
    signal.alarm(120)
    start = time.monotonic()
    F5, g15, g25 = canonical()
    F, g1, g2 = map(dehom, (F5, g15, g25))
    G = mul(g1, g2)
    nu = [add(mul(g2, diff(g1, i)), scale(mul(g1, diff(g2, i)), -1)) for i in range(4)]
    target = {(i, j): rem(add(mul(diff(F, i), nu[j]), scale(mul(diff(F, j), nu[i]), -1)), F)
              for i, j in it.combinations(range(4), 2)}
    pairs = list(it.combinations(range(5), 2))
    columns = [(pair, m) for pair in pairs for m in mons(5, 2)]
    rows = {}
    labels = {}
    rhs = {}

    def put(label, col, value):
        if value:
            rows.setdefault(label, {})[col] = rows.setdefault(label, {}).get(col, 0)+value

    cache = {}
    for col, ((i, j), m5) in enumerate(columns):
        # The full divergence-zero Euler gauge: all25 coefficients, one redundancy.
        for derivative, component, sign in [(i, j, 1), (j, i, -1)]:
            if m5[derivative]:
                exponent = list(m5)
                exponent[derivative] -= 1
                put(('div', component, tuple(exponent)), col, sign*m5[derivative])
        affine = []
        if j < 4:
            affine.append(((i, j), m5[:4], 1))
        else:
            for p, q in it.combinations(range(4), 2):
                if q == i:
                    exponent = list(m5[:4]); exponent[p] += 1
                    affine.append(((p, q), tuple(exponent), 1))
                if p == i:
                    exponent = list(m5[:4]); exponent[q] += 1
                    affine.append(((p, q), tuple(exponent), -1))
        for (p, q), exponent, sign in affine:
            if exponent not in cache:
                cache[exponent] = rem(mul(G, {exponent: 1}), F)
            complement = tuple(k for k in range(4) if k not in (p, q))
            sign *= (-1)**(p+q-1)
            for m, c in cache[exponent].items():
                put(('restriction', complement, m), col, sign*c)
    for pair, polynomial in target.items():
        for m, c in polynomial.items():
            label = ('restriction', pair, m)
            rows.setdefault(label, {})
            rhs[label] = c
    keys = sorted(rows, key=repr)
    matrix = [rows[k] for k in keys]
    values = [rhs.get(k, 0) for k in keys]
    matrix_serial = [{'label': k, 'row': [[c, str(v)] for c, v in sorted(row.items()) if v],
                      'rhs': str(value)} for k, row, value in zip(keys, matrix, values)]
    matrix_hash = hashlib.sha256(json.dumps(matrix_serial, separators=(',', ':')).encode()).hexdigest()
    print(f'complete rows={len(matrix)} cols=150; building exact elimination', flush=True)
    pivots = {}
    inconsistency = None
    for index, (row, value) in enumerate(zip(matrix, values)):
        row = {c: Q(v) for c, v in row.items() if v}
        value = Q(value)
        combination = {index: Q(1)}
        while row:
            p = min(row)
            if p not in pivots:
                factor = row[p]
                row = {c: v/factor for c, v in row.items()}
                value /= factor
                combination = {k: v/factor for k, v in combination.items()}
                pivots[p] = (row, value, combination)
                break
            reduced, rv, rc = pivots[p]
            factor = row[p]
            axpy(row, reduced, -factor)
            value -= factor*rv
            axpy(combination, rc, -factor)
        else:
            if value:
                inconsistency = (value, combination)
                break
    result = {'status': 'PASS', 'method': 'One QQ sparse all150 ambient restriction lift; no GB/profile.',
              'full_ambient_column_count': 150, 'all_restriction_twoforms': 6,
              'full_row_count': len(matrix), 'rows': matrix_serial, 'matrix_sha256': matrix_hash,
              'column_order': [{'pair': pair, 'exponents': m} for pair, m in columns],
              'canonical_homogeneous_F': serial(F5), 'canonical_homogeneous_g1': serial(g15),
              'canonical_homogeneous_g2': serial(g25), 'single_divisor_leading_monomial':
                  list(max(F, key=lambda m: (sum(m), m))),
              'gauge': 'deltaP_j=sum_i partial_i P_ij=0, all25 coefficients; validEuler24 gauge',
              'source_scope': 'The entire canonicalindex2 node source is onecanonicalpair or empty by frozen whole-orbit entrance.'}
    # Exact irreducibility witness for the primitive conic in(d,e).
    conic = {}
    for m, c in F5.items():
        mon = (m[3], m[4])
        conic[mon] = conic.get(mon, 0)+c*(1**m[0])*(1**m[1])*(2**m[2])
    conic = {m: c for m, c in conic.items() if c}
    assert conic == {(2, 0): -2, (1, 1): 3, (0, 2): -1,
                     (1, 0): -13, (0, 1): 10, (0, 0): -22}
    h = [[Q(-2), Q(3, 2), Q(-13, 2)], [Q(3, 2), Q(-1), Q(5)],
         [Q(-13, 2), Q(5), Q(-22)]]
    det = sum((-1)**sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3))*
              h[0][p[0]]*h[1][p[1]]*h[2][p[2]] for p in it.permutations(range(3)))
    assert det == Q(1, 4)
    result['absolute_irreducibility_certificate'] = {
        'primitive_content_argument': 'The three quadratic coefficients in(d,e) are -b²c,b(ac+b²),-ab², gcd=b; constant atb0=-a5-c5, so fullcontent1.',
        'conic_specialization_spatial_a_b_c': [1, 1, 2],
        'conic_coefficients': [[list(m), str(c)] for m, c in sorted(conic.items())],
        'homogenized_conic_matrix': [[str(c) for c in row] for row in h], 'determinant': str(det),
        'implication': 'Generic conic nonsingular, hence irreducible over algebraicclosure C(a,b,c); primitive Gauss lemma gives absoluteirreducible homogeneousF.'}
    if inconsistency:
        value, combination = inconsistency
        check = {}
        cv = Q(0)
        for k, c in combination.items():
            axpy(check, matrix[k], c)
            cv += c*values[k]
        assert not check and cv == value and value
        result.update({'outcome': 'NO_AMBIENT_LIFT_WHOLE_CANONICAL_NODE_EMPTY',
                       'left_inconsistency_functional': [[k, str(c)] for k, c in sorted(combination.items())],
                       'functional_all150_columns_zero': True, 'functional_rhs_nonzero': str(cv),
                       'pivot_count_before_certificate': len(pivots)})
    else:
        assert len(pivots) == 150, f'unexpected lift kernel, pivots={len(pivots)}'
        solution = [Q(0)]*150
        for p in sorted(pivots, reverse=True):
            row, value, _ = pivots[p]
            solution[p] = value-sum(v*solution[c] for c, v in row.items() if c != p)
        assert all(sum(v*solution[c] for c, v in row.items()) == val for row, val in zip(matrix, values))
        P = {pair: {} for pair in pairs}
        for ((pair, m), coefficient) in zip(columns, solution):
            if coefficient:
                P[pair][m] = coefficient
        divergence = [add(*[diff(raw(P, i, j), i) for i in range(5)]) for j in range(5)]
        assert not any(divergence)
        J = {}
        for i, j, k in it.combinations(range(5), 3):
            J[(i, j, k)] = add(*[mul(raw(P, i, l), diff(raw(P, j, k), l))
                                  for l in range(5)],
                                *[mul(raw(P, j, l), diff(raw(P, k, i), l))
                                  for l in range(5)],
                                *[mul(raw(P, k, l), diff(raw(P, i, j), l))
                                  for l in range(5)])
        EJ = [add(*[scale(mul(var(5, x), J[tuple(a for a in quad if a != x)]), (-1)**j)
                    for j, x in enumerate(quad)]) for quad in it.combinations(range(5), 4)]
        assert not any(EJ)
        Pf = {}
        for i in range(5):
            j, k, l, m = [a for a in range(5) if a != i]
            pf4 = add(mul(raw(P, j, k), raw(P, l, m)),
                      scale(mul(raw(P, j, l), raw(P, k, m)), -1),
                      mul(raw(P, j, m), raw(P, k, l)))
            Pf = add(Pf, scale(mul(var(5, i), pf4), (-1)**i))
        pf_scalar = Q(0) if not Pf else Q(Pf[next(iter(F5))], F5[next(iter(F5))])
        assert add(Pf, scale(F5, -pf_scalar)) == {}
        result.update({'outcome': 'ACTUAL_RANKFOUR_ONE_ORBIT_PARENT_UNIDENTIFIED' if pf_scalar else
                       'AMBIENT_LIFT_HAS_ZERO_PF_WHOLE_CANONICAL_NODE_EMPTY',
                       'all150_solution': [str(c) for c in solution],
                       'homogeneous_bivector': {f'{i},{j}': serial(p) for (i, j), p in P.items()},
                       'all_original_linear_rows_verified': True, 'all_divergence_coefficients_zero': True,
                       'all350_projective_Jacobi_coefficients_zero': True,
                       'signed_actual_Pfaffian_126': serial(Pf), 'actual_Pfaffian_scalar_times_F': str(pf_scalar),
                       'all126_actual_Pfaffian_residuals_zero': True})
    result['seconds'] = time.monotonic()-start
    result['peak_rss_kib'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result['producer_source_bytes'] = Path(__file__).read_text()
    result['producer_source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    job = Path('research/classification/jobs/canonical_index_two_node_canonical_ambient_lift.json')
    result['job'] = {'path': str(job), 'sha256': hashlib.sha256(job.read_bytes()).hexdigest()}
    path = Path(args.output)
    with path.open('x') as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
        stream.write('\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['status', 'outcome', 'full_row_count', 'seconds', 'peak_rss_kib']}) )


if __name__ == '__main__':
    main()

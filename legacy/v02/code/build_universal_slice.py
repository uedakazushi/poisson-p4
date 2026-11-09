#!/usr/bin/env python3
"""Exact positively weighted transverse equations at x0^2*d3^d4.

Requires NumPy and the previous verification modules.  All elimination is over Q.
The output is a research certificate, not an irreducible decomposition.
"""
from collections import defaultdict, Counter
from fractions import Fraction as F
from pathlib import Path
import itertools
import json
import sys
import time

HERE = Path(__file__).resolve().parent
import numpy as np
import poisson_p4_components_check as p
import verify_rank2_components as v

ROOT = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE.parent/'data'
ROOT.mkdir(parents=True, exist_ok=True)
W = (2, 0, 0, -1, -1)
STAR = p.BIVECTOR_INDEX[((3, 4), (2, 0, 0, 0, 0))]
PSTAR = {(3, 4): {(2, 0, 0, 0, 0): 1}}


def rref(matrix, transform=False):
    a = [[F(int(x)) if not isinstance(x, F) else x for x in row] for row in matrix]
    m, n = len(a), len(a[0])
    t = [[F(i == j) for j in range(m)] for i in range(m)] if transform else None
    pivots, rank = [], 0
    for c in range(n):
        k = next((i for i in range(rank, m) if a[i][c]), None)
        if k is None:
            continue
        a[rank], a[k] = a[k], a[rank]
        if transform:
            t[rank], t[k] = t[k], t[rank]
        q = a[rank][c]
        a[rank] = [x/q for x in a[rank]]
        if transform:
            t[rank] = [x/q for x in t[rank]]
        nz = [j for j, x in enumerate(a[rank]) if x]
        for i in range(m):
            if i != rank and a[i][c]:
                q = a[i][c]
                for j in nz:
                    a[i][j] -= q*a[rank][j]
                if transform:
                    for j in range(m):
                        t[i][j] -= q*t[rank][j]
        pivots.append(c)
        rank += 1
        if rank == m:
            break
    return a, pivots, t


def trans(a):
    return list(map(list, zip(*a)))


def weight(raw_index):
    pair, alpha = p.BIVECTOR_BASIS[raw_index]
    return 6 + sum(W[i] for i in pair) - sum(a*w for a, w in zip(alpha, W))


def clean(q):
    return {m: c for m, c in q.items() if c}


def add_to(q, r, scalar=1):
    if scalar:
        for m, c in r.items():
            q[m] += scalar*c


def mul(q, r):
    out = defaultdict(F)
    for m, c in q.items():
        for n, d in r.items():
            out[tuple(sorted(m+n))] += c*d
    return clean(out)


def substitute(q, expressions):
    out = defaultdict(F)
    for m, c in q.items():
        term = {(): c}
        for k in m:
            term = mul(term, expressions[k])
        add_to(out, term)
    return clean(out)


def jacobi_equations(raw):
    # Each list entry is (x-exponents, parameter monomial, scalar).
    bv = defaultdict(list)
    for k, r in enumerate([STAR]+raw):
        pair, alpha = p.BIVECTOR_BASIS[r]
        mon = () if k == 0 else (k-1,)
        bv[pair].append((alpha, mon, 1))
        bv[pair[::-1]].append((alpha, mon, -1))
    equations = [defaultdict(F) for _ in p.FOUR_VECTOR_BASIS]
    for triple in p.TRIPLES:
        for i, j, k in (triple, (triple[1], triple[2], triple[0]),
                        (triple[2], triple[0], triple[1])):
            for h in range(5):
                for a, m, c in bv[i, h]:
                    for b, n, d in bv[j, k]:
                        if not b[h]:
                            continue
                        alpha = [u+v for u, v in zip(a, b)]
                        alpha[h] -= 1
                        for e in range(5):
                            sg = p.permutation_sign((e,)+triple)
                            if not sg:
                                continue
                            beta = alpha.copy()
                            beta[e] += 1
                            row = p.FOUR_VECTOR_INDEX[(tuple(sorted((e,)+triple)), tuple(beta))]
                            equations[row][tuple(sorted(m+n))] += sg*c*d*b[h]
    return [clean(q) for q in equations]


def encode(q):
    return [[list(m), c.numerator, c.denominator] for m, c in sorted(q.items())]


def polynomial_row_basis(equations):
    """Constant-linear Gaussian elimination in a sparse monomial basis."""
    basis = {}
    selected = []
    for i, equation in enumerate(equations):
        q = equation.copy()
        while q:
            pivot = min(q, key=lambda m: (len(m), m))
            if pivot not in basis:
                factor = q[pivot]
                q = {m: c/factor for m, c in q.items()}
                basis[pivot] = q
                selected.append(i)
                break
            factor = q[pivot]
            for m, c in basis[pivot].items():
                new = q.get(m, F(0))-factor*c
                if new:
                    q[m] = new
                else:
                    q.pop(m, None)
    return selected, list(basis.values())


def singular_polynomial(q):
    pieces = []
    for mon, coef in sorted(q.items(), key=lambda x: (len(x[0]), x[0])):
        num, den = coef.numerator, coef.denominator
        factors = []
        if abs(num) != 1 or not mon:
            factors.append(str(abs(num)))
        if den != 1:
            if not factors:
                factors.append('1')
            factors[-1] = '('+factors[-1]+'/'+str(den)+')'
        factors += ['u'+str(k+1)+(('^'+str(e)) if e != 1 else '')
                    for k, e in sorted(Counter(mon).items())]
        body = '*'.join(factors) or '1'
        pieces.append(('-' if num < 0 else ('+' if pieces else ''))+body)
    return ''.join(pieces) or '0'


def main():
    t0 = time.time()
    vertical = p.vertical_bivector_matrix()
    vr, vp, _ = rref(vertical.T.tolist())
    active = [i for i in range(150) if i not in vp]
    def quotient(vec):
        return [F(int(vec[k]))-sum(vr[a][k]*int(vec[h]) for a, h in enumerate(vp))
                for k in active]
    orbit = []
    generators = []
    star = active.index(STAR)
    for i, j in itertools.product(range(5), repeat=2):
        a = np.zeros((5, 5), dtype=np.int64)
        a[i, j] = 1
        vec = quotient(p.bivector_vector(v.linear_action_derivative(PSTAR, a)))
        vec[star] = F(0)
        orbit.append(vec)
        generators.append((i, j))
    _, op, _ = rref(orbit)
    raw = [active[j] for j in range(126) if j != star and j not in op]
    assert len(op) == 8 and len(raw) == 117
    Jfull = p.linearized_jacobi_matrix(PSTAR)
    J = Jfull[:, raw].tolist()
    _, independent_rows, _ = rref(trans(J))
    jr, jp, transform = rref([J[i] for i in independent_rows], transform=True)
    free = [i for i in range(len(raw)) if i not in jp]
    assert len(jp) == 39 and len(free) == 78
    wt = [weight(k) for k in raw]
    report = {
        'one_parameter_subgroup_weights': W,
        'star_raw_index': STAR,
        'vertical_pivots': vp,
        'orbit_pivots_raw': [active[k] for k in op],
        'orbit_weights': dict(sorted(Counter(weight(active[k]) for k in op).items())),
        'linear_Jacobi_rank_exact': len(jp),
        'affine_tangent_dimension': 87,
        'projective_tangent_dimension': 86,
        'orbit_dimension': 8,
        'transverse_tangent_dimension': 78,
        'slice_ambient_weights': dict(sorted(Counter(wt).items())),
        'eliminated_weights': dict(sorted(Counter(wt[k] for k in jp).items())),
        'tangent_slice_weights': dict(sorted(Counter(wt[k] for k in free).items())),
        'raw_slice_indices': raw,
        'raw_free_indices': [raw[k] for k in free],
        'equation_rows_used_for_elimination': independent_rows,
        'free_variable_basis': [{'pair': p.BIVECTOR_BASIS[raw[k]][0],
                                 'exponents': p.BIVECTOR_BASIS[raw[k]][1],
                                 'weight': wt[k]} for k in free],
    }
    print(json.dumps({k: report[k] for k in ('linear_Jacobi_rank_exact', 'orbit_weights',
          'tangent_slice_weights')}, indent=2), flush=True)
    equations = jacobi_equations(raw)
    # Check this independent generation against the old linearization.
    for i, q in enumerate(equations):
        assert not q.get((), 0)
        assert all(q.get((j,), 0) == J[i][j] for j in range(len(raw)))
        weights = {sum(wt[k] for k in m) for m in q}
        assert len(weights) <= 1
    pivot_equations = []
    for row in transform:
        q = defaultdict(F)
        for a, c in enumerate(row):
            add_to(q, equations[independent_rows[a]], c)
        pivot_equations.append(clean(q))
    free_index = {k: a for a, k in enumerate(free)}
    expressions = {k: {(free_index[k],): F(1)} for k in free}
    for a in sorted(range(len(jp)), key=lambda a: wt[jp[a]]):
        k = jp[a]
        q = pivot_equations[a].copy()
        assert q.pop((k,)) == 1
        # RREF eliminates all other pivot variables in the linear part.
        assert all(j in expressions for m in q for j in m), (k, wt[k])
        expressions[k] = {m: -c for m, c in substitute(q, expressions).items()}
    report['elimination_expression_count'] = len(jp)
    report['elimination_expression_max_degree'] = max(len(m) for k in jp for m in expressions[k])
    report['elimination_expression_terms'] = sum(len(expressions[k]) for k in jp)
    reduced = []
    for i, q in enumerate(equations):
        reduced.append(substitute(q, expressions))
        if i % 50 == 0:
            print('Reduced row', i, 'seconds', round(time.time()-t0, 2), flush=True)
    assert all(not reduced[i] for i in independent_rows)
    assert all(len(m) >= 2 for q in reduced for m in q)
    unique = {}
    for i, q in enumerate(reduced):
        if not q:
            continue
        first = q[min(q)]
        norm = tuple((m, c/first) for m, c in sorted(q.items()))
        unique.setdefault(norm, []).append(i)
    report['nonzero_reduced_rows'] = sum(bool(q) for q in reduced)
    report['distinct_reduced_equations_up_to_scalar'] = len(unique)
    report['reduced_equation_max_degree'] = max(len(m) for q in reduced for m in q)
    report['reduced_equation_terms'] = sum(len(q) for q in reduced)
    report['reduced_equation_degree_counts'] = dict(sorted(Counter(max(map(len, q))
             for q in reduced if q).items()))
    quadratic_parts = [{m: c for m, c in q.items() if len(m) == 2} for q in reduced]
    quadratic_selected, quadratic_basis = polynomial_row_basis(quadratic_parts)
    selected, equation_basis = polynomial_row_basis(reduced)
    report['quadratic_obstruction_span_rank'] = len(quadratic_selected)
    report['constant_linear_equation_span_rank'] = len(selected)
    report['independent_original_equation_rows'] = selected
    report['independent_equation_terms'] = sum(len(reduced[k]) for k in selected)
    report['kappa_variable_zero_based_index'] = free.index(raw.index(p.BIVECTOR_INDEX[((1,2),(2,0,0,0,0))]))
    report['kappa_variable_name'] = 'u'+str(report['kappa_variable_zero_based_index']+1)
    report['equation_weight_counts'] = dict(sorted(Counter(
        sum(wt[free[j]] for j in next(iter(reduced[k]))) for k in selected).items()))
    report['seconds'] = round(time.time()-t0, 3)
    with (ROOT/'universal_slice_certificate.json').open('w') as f:
        json.dump(report, f, indent=2)
    with (ROOT/'universal_slice_equations.json').open('w') as f:
        json.dump({'parameter_count': 78, 'parameter_weights': [wt[k] for k in free],
                   'eliminated_raw_variables': [{'raw_index': raw[k], 'polynomial': encode(expressions[k])}
                                                 for k in jp],
                   'equations': [{'original_row': k, 'polynomial': encode(reduced[k])}
                                 for k in selected],
                   'quadratic_obstruction_basis': [encode(q) for q in quadratic_basis]
                   }, f, separators=(',', ':'))
    names = ','.join('u'+str(i+1) for i in range(78))
    weights = ','.join(str(wt[k]) for k in free)
    with (ROOT/'universal_slice.sing').open('w') as f:
        f.write('// Exact transverse Poisson scheme at x0^2*d3^d4 over QQ.\n')
        f.write('// No primary decomposition or saturation is asserted to be computed.\n')
        f.write('ring R=0,('+names+'),wp('+weights+');\n')
        f.write('ideal I=\n'+',\n'.join(singular_polynomial(reduced[k]) for k in selected)+';\n')
        f.write('poly kappa='+report['kappa_variable_name']+';\n')
        f.write('// Optional expensive commands, intentionally not run by this file:\n')
        f.write('// LIB "elim.lib";\n// ideal Irank4=sat(I,ideal(kappa))[1];\n')
        f.write('// LIB "primdec.lib";\n// list minimalPrimes=minAssGTZ(I);\n')
        f.write('// A minimal prime contains kappa iff its full Poisson component has generic rank2.\n')
    print(json.dumps({k:v for k,v in report.items() if not isinstance(v,list)}, indent=2))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Independent prescribed S(2,3) conductor and rank-two sextic review.

Only the reviewer's earlier standard-library Fraction arithmetic is reused.
No producer or CAS imports; no projection search. Geometry is audited in the
accompanying Japanese note, not certified by this finite polynomial run.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools as it
import json
import math
from pathlib import Path
import resource
import time

from classification_projected_veronese_zero_review import (
    Ring, plus, times, scale, derivative, monomials, sparse_rref, det,
)

ROOT = Path(__file__).resolve().parents[1]


def substitute(p, values, out_ring):
    return plus(*(scale(times(out_ring.const(1), *(
        out_ring.power(v, k) for v, k in zip(values, e))), c)
        for e, c in p.items()))


def divide(p, q):
    """Literal lexicographic single-divisor division over QQ."""
    assert q
    quotient = {}; remainder = {}; p = dict(p)
    f = max(q); b = q[f]
    while p:
        e = max(p); a = p[e]
        if all(i >= j for i, j in zip(e, f)):
            mon = {tuple(i-j for i, j in zip(e, f)): a/b}
            quotient = plus(quotient, mon)
            p = plus(p, scale(times(mon, q), -1))
        else:
            remainder[e] = a
            del p[e]
    return quotient, remainder


def gcd_univariate(*ps):
    g = {}
    for p in ps:
        p = dict(p)
        while p:
            _, rem = divide(g, p)
            g, p = p, rem
        if g: g = scale(g, 1/g[max(g)])
    return g


def encode(p):
    return [[list(e), str(c)] for e, c in sorted(p.items(), reverse=True)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', required=True, type=Path)
    args = ap.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU, (119, 120))
    resource.setrlimit(resource.RLIMIT_AS, (1536*1024**2, 1536*1024**2))
    started = time.monotonic()
    manifest_path = ROOT/'research/classification/inputs/scroll23_conductor_packet_manifest_20261107.json'
    manifest = json.loads(manifest_path.read_text())
    for row in manifest['files']:
        assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest() == row['sha256']
    saved = json.loads((ROOT/'reports/CLASSIFICATION_SCROLL23_CONDUCTOR_WITNESS_20261107.json').read_text())

    S = Ring(['s', 't']); s, t = map(S.var, S.names)
    A = [times(s, s), times(s, t), times(t, t), {}]
    B = [S.power(s, 3), S.power(t, 3), times(s, s, t), plus(S.power(s, 3), S.power(t, 3))]
    minors = [plus(times(A[i], B[j]), scale(times(A[j], B[i]), -1))
              for i, j in it.combinations(range(4), 2)]
    assert minors == [S.parse(v) for v in saved['all_six_fiber_independence_minors']]
    dehom_minors = [{(e[1],): c for e, c in p.items()} for p in minors]
    assert gcd_univariate(*dehom_minors) == {(0,): Q(1)}
    at_infinity = [[p.get((0, d), 0) for p in row] for row, d in [(A, 2), (B, 3)]]
    assert len(sparse_rref(at_infinity)[0]) == 2
    section_matrix = [[p.get(e, Q(0)) for p in row]
                      for row, d in [(A, 2), (B, 3)] for e in monomials(2, d)]
    assert len(sparse_rref(section_matrix)[0]) == 4

    Y = Ring(['s', 't', 'u', 'v']); sy, ty, u, v = map(Y.var, Y.names)
    Z = [plus(times(u, substitute(a, [sy, ty], Y)),
              times(v, substitute(b, [sy, ty], Y))) for a, b in zip(A, B)]
    exponents = monomials(4, 5)
    columns = [times(Y.const(1), *(Y.power(z, k) for z, k in zip(Z, e))) for e in exponents]
    # Use all H0(5H) basis rows, including a row whose coefficients may all
    # vanish for this prescribed subsystem, rather than only the support union.
    row_labels = sorted(((15-i-j, j, i, 5-i)
                         for i in range(6) for j in range(16-i)), reverse=True)
    matrix = [[p.get(e, Q(0)) for p in columns] for e in row_labels]
    red, _, selected_rows = sparse_rref(matrix)
    assert len(matrix) == 81 and len(red) == 55
    free = next(j for j in range(56) if j not in red)
    kernel = [Q(0)]*56; kernel[free] = 1
    for pivot, row in red.items(): kernel[pivot] = -row.get(free, Q(0))
    denominator = math.lcm(*(q.denominator for q in kernel))
    primitive_gcd = math.gcd(*(int(q*denominator) for q in kernel))
    kernel = [q*denominator/primitive_gcd for q in kernel]
    if next(q for q in kernel if q) < 0: kernel = [-q for q in kernel]
    H = Ring(['x', 'y', 'z', 'w'])
    F = {e: q for e, q in zip(exponents, kernel) if q}
    affine_F = {e[:3]: q for e, q in F.items()}
    assert affine_F == Ring(['x', 'y', 'z']).parse(saved['exact_quintic_equation'])
    assert not substitute(F, Z, Y)
    pivots = sorted(red)
    implicit_minor = det([[matrix[i][j] for j in pivots] for i in selected_rows])
    assert implicit_minor

    R = Ring(['r', 'T']); r, T = map(R.var, R.names)
    d = plus(R.const(1), R.power(T, 3))
    z = [plus(r, R.const(1)), plus(times(r, T), R.power(T, 3)),
         plus(times(r, T, T), T), d]
    jnum = plus(times(plus(times(derivative(z[0], 1), d), scale(times(z[0], derivative(d, 1)), -1)), derivative(z[1], 0)),
                scale(times(derivative(z[0], 0), plus(times(derivative(z[1], 1), d), scale(times(z[1], derivative(d, 1)), -1))), -1))
    assert jnum == R.parse(saved['Jacobian_numerator'])
    numerator = substitute(derivative(F, 2), z, R)
    G, remainder = divide(numerator, jnum)
    assert not remainder and G == R.parse(saved['conductor_affine_polynomial'])
    assert not plus(numerator, scale(times(jnum, G), -1))
    assert max(e[0] for e in G) == 3 and max(e[1] for e in G) == 4
    assert all(j <= 6-i for i, j in G)
    global_G = {(6-i-j, j, i, 3-i): c for (i, j), c in G.items()}
    assert global_G == Y.parse(saved['conductor_global_section_3H_minus3f'])
    infinity = {e: c for e, c in global_G.items() if e[0] == 0}
    negative_section = {e: c for e, c in global_G.items() if e[3] == 0}
    assert infinity == Y.parse('-t**4*u**2*v')
    assert negative_section == Y.parse('-s*u**3*(4*s**2+6*s*t+3*t**2)')

    # Primitivity in QQ[T][r] and one full-degree separable r-specialization
    # prove squarefreeness, without multivariate factorization or gcd imports.
    content = gcd_univariate(*[{(j,): c for (ri, j), c in G.items() if ri == i} for i in range(4)])
    assert content == {(0,): Q(1)}
    G0 = {(i,): c for (i, j), c in G.items() if j == 0}
    assert max(e[0] for e in G0) == 3
    assert gcd_univariate(G0, derivative(G0, 0)) == {(0,): Q(1)}

    unknowns = [('A', i, j) for i in range(3) for j in range(5)]
    unknowns += [('B', i, j) for i in range(4) for j in range(4)]
    closed_columns = []
    for kind, i, j in unknowns:
        mon = {(i, j): Q(1)}
        closed_columns.append(plus(times(G, derivative(mon, 1)), scale(times(mon, derivative(G, 1)), -1)) if kind == 'A'
                              else plus(scale(times(G, derivative(mon, 0)), -1), times(mon, derivative(G, 0))))
    labels = sorted(set().union(*(set(p) for p in closed_columns)), reverse=True)
    closed_matrix = [[p.get(e, Q(0)) for p in closed_columns] for e in labels]
    cert = saved['closed_form_absolute_factor_certificate']
    assert unknowns == [tuple(v) for v in cert['unknown_order']]
    assert labels == [tuple(v) for v in cert['all_spatial_row_exponents']]
    assert closed_matrix == [[Q(a) for a in row] for row in cert['full_matrix']]
    closed_red, _, _ = sparse_rref(closed_matrix)
    assert len(closed_matrix) == 44 and len(closed_red) == 30
    dlog = [(derivative(G, 0) if kind == 'A' else derivative(G, 1)).get((i, j), Q(0))
            for kind, i, j in unknowns]
    assert dlog == [Q(a) for a in cert['dlog_f_vector']] and any(dlog)
    assert all(not sum(a*b for a, b in zip(row, dlog)) for row in closed_matrix)
    minor = cert['nonzero_rank_minor']
    closed_det = det([[closed_matrix[i][j] for j in minor['columns']] for i in minor['rows']])
    assert closed_det == Q(minor['determinant']) == 3486368987873280

    inputs = [row['path'] for row in manifest['files']] + [
        str(manifest_path.relative_to(ROOT)),
        'research/classification/RANK_TWO_SEXTIC_SURFACE_NONCONICAL_DIVISOR_EXCLUSION.md',
        'code/classification_projected_veronese_zero_review.py',
    ]
    result = dict(task_id='G', all_passed=True,
        review_method='Independent standard-library Fraction reconstruction of the one prescribed projection. No producer/CAS imports or extra projection search. All geometry is separately audited in the review note.',
        one_prescribed_projection_only=True,
        all_six_minors= [encode(p) for p in minors], no_finite_or_infinite_basepoint=True,
        all_four_section_rank=4,
        implicit_reconstruction=dict(matrix_shape=[81, 56], rank=55,
            coefficient_monomials=exponents, substitution_row_monomials=row_labels,
            full_matrix=[[str(a) for a in row] for row in matrix], primitive_quintic=encode(F),
            nonzero55minor=dict(rows=selected_rows, columns=pivots, determinant=str(implicit_minor)),
            exact_saved_quintic_equal=True, all_homogeneous_pullback_coefficients_zero=True),
        canonical_conductor=dict(Jacobian_numerator=encode(jnum),
            full_canonical_frame_polynomial_identity=True, scalar='1',
            affine_polynomial=encode(G), global_section=encode(global_G),
            infinity_boundary=encode(infinity), negative_section_boundary=encode(negative_section),
            whole_global_boundary_retained=True, coefficient_content='1',
            full_degree_T0_r_polynomial=encode(G0), T0_derivative_gcd='1',
            squarefree_over_C=True),
        absolute_irreducibility=dict(full_closed_form_matrix_shape=[44, 31], rank=30,
            all_matrix_entries_match=True, dlog_spans_whole_kernel=True,
            nonzero30minor_determinant=str(closed_det),
            logic='Every absolute factor contributes an independent dlog in the complete rectangular bounds. Squarefreeness and kernel1 force one absolute factor. Nonzero global boundary restrictions exclude additional components.'),
        mathematical_scope='Whole rank2 degree6 branch: smooth normalization113/122 excluded; conical014/005 excluded and023 image<=24, conditional on the independently audited relative conductor-family/residue/dimension proof in the review note. No E membership or component count.',
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        seconds=time.monotonic()-started,
        maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x') as out:
        json.dump(result, out, ensure_ascii=False, indent=2); out.write('\n')
    print(json.dumps(dict(all_passed=True, implicit_rank=55, closed_form_rank=30,
                         seconds=result['seconds'], maximum_RSS_KiB=result['maximum_RSS_KiB'])))


if __name__ == '__main__':
    main()

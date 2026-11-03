#!/usr/bin/env python3
"""One-variable physical CA/CB source4 jets and normalized even graph differences."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_NOETHER_FIBER_EVEN_QUADRATIC_20261103_v63'
Z = (F(0),)*5
ONE = (F(1),)+(F(0),)*4
PVAR = (F(0), F(1), F(0), F(0), F(0))


def add(a, b, c=1):
    return tuple(x+c*y for x, y in zip(a, b))


def scale(a, c):
    return tuple(c*x for x in a) if c else Z


def mul(a, b):
    if a == Z or b == Z:
        return Z
    out = [F(0)]*5
    for i, x in enumerate(a):
        if x:
            for j in range(5-i):
                if b[j]:
                    out[i+j] += x*b[j]
    return tuple(out)


def sumj(values):
    out = Z
    for v in values:
        out = add(out, v)
    return out


def mm(A, B):
    return [[sumj(mul(A[i][k], B[k][j]) for k in range(len(B)))
             for j in range(len(B[0]))] for i in range(len(A))]


def ma(A, B, c=1):
    return [[add(a, b, c) for a, b in zip(ar, br)] for ar, br in zip(A, B)]


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def enc(A):
    if isinstance(A, (tuple, list)):
        return [enc(v) for v in A]
    return str(A)


def dump_gz(p, value):
    with p.open('wb') as f:
        with gzip.GzipFile(fileobj=f, mode='wb', filename='', mtime=0) as g:
            g.write((json.dumps(value, separators=(',', ':'))+'\n').encode())


def main(out):
    start = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    out.mkdir(parents=True, exist_ok=False)
    reduction_path = ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json'
    first_path = ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz'
    tangent_path = ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63/certificate.json.gz'
    finite_gate = ROOT/'reports/AMBIENT_B18_ACTUAL_D_FINITE13_20261103_v63_REVIEW.json'
    r, first, tangent, finite = map(read, [reduction_path, first_path, tangent_path, finite_gate])
    assert finite['all_passed'] and finite['closed_fiber_length_allowed'] == [12, 13]
    regular = r['common28_linear_target_rows']
    fold = r['common29th_target_row']
    selected = tangent['common39_physical_coordinate_rows']
    odd = tangent['common10_odd_target_rows']
    assert len(regular) == 28 and fold == 61 and len(odd) == 10
    full_T = s.Matrix(tangent['common150_by39_normalized_tangent'])
    assert full_T[selected, :] == s.eye(39)
    h0 = [[F(q) for q in row] for row in full_T[:, 29:].tolist()]
    even = [s.Matrix(tangent['families'][name]['full150_by29_even_cotangent_embedding']) for name in ['CA', 'CB']]
    assert even[0] == even[1]
    h, x, y, z, w = s.symbols('h x y z w')
    X = s.Matrix([h, x, y, z, w])
    pairs = list(itertools.combinations(range(5), 2))
    mons = list(itertools.combinations_with_replacement(range(5), 2))
    def vec(P):
        return s.Matrix([s.Poly(P[i, j], *X).coeff_monomial(X[a]*X[b]) for i, j in pairs for a, b in mons])
    def canonical(P):
        d = s.Matrix([sum(s.diff(P[i, j], X[i]) for i in range(5)) for j in range(5)])
        return (P-(X*d.T-d*X.T)/5).applyfunc(s.expand)
    def raw(f, g, H, cb):
        P = s.zeros(5)
        P[1, 2] = f
        P[3, 4] = g if cb else g*h*h
        P[2, 3] = (-s.diff(f, x)*z/2 if cb else 0)-s.diff(H, w)
        P[2, 4] = (-s.diff(f, x)*w/2 if cb else 0)+s.diff(H, z)
        return P-P.T
    sw = s.eye(5)
    sw[0, 0] = sw[1, 1] = 0
    sw[0, 1] = sw[1, 0] = 1
    change = lambda P: (sw*P.subs({h: x, x: h}, simultaneous=True)*sw.T).applyfunc(s.expand)
    cubic_monomials = [s.prod((h, x, z, w)[i] for i in m) for m in itertools.combinations_with_replacement(range(4), 3) if any(i >= 2 for i in m)]
    ey = s.eye(5)[:, 2]
    D = s.Matrix([x*x, 0, 0, h*x, h*h+z*z])
    P0 = ey*D.T-D*ey.T
    point = list(map(F, vec(P0)))
    oi = r['orbit_GL_generator_indices']
    O = []
    for ij in oi:
        E = s.zeros(5)
        E[ij//5, ij % 5] = 1
        MX = E*X
        O.append(vec((E*P0+P0*E.T-sum((MX[k]*P0.diff(X[k]) for k in range(5)), s.zeros(5))).applyfunc(s.expand)))
    O = s.Matrix.hstack(*O)
    I = [[ONE if i == j else Z for j in range(5)] for i in range(5)]
    # Original physical divergence and strict Jacobi coefficient operations.
    div = [[0]*150 for _ in range(25)]
    cubic_spatial = list(itertools.combinations_with_replacement(range(5), 3))
    jac = [defaultdict(int) for _ in range(350)]
    pi = {pair: i for i, pair in enumerate(pairs)}
    def entry(i, j):
        return (pi[(i, j)], 1) if i < j else (pi[(j, i)], -1)
    for a, ((i, j), (m, n)) in enumerate(itertools.product(pairs, mons)):
        for k in [m, n]:
            other = n if k == m else m
            if k == i:
                div[5*j+other][a] += 1
            if k == j:
                div[5*i+other][a] -= 1
    for ti, triple in enumerate(itertools.combinations(range(5), 3)):
        for i, j, k in [triple, (triple[1], triple[2], triple[0]), (triple[2], triple[0], triple[1])]:
            for l in range(5):
                if l == i:
                    continue
                ai, sa = entry(i, l)
                bi, sb = entry(j, k)
                for am, (u, v) in enumerate(mons):
                    for bm, (a, b) in enumerate(mons):
                        count = (a == l)+(b == l)
                        if not count:
                            continue
                        other = b if a == l else a
                        row = 35*ti+cubic_spatial.index(tuple(sorted((u, v, other))))
                        ia, ib = 15*ai+am, 15*bi+bm
                        jac[row][tuple(sorted((ia, ib)))] += sa*sb*count
    jac = [{k: v for k, v in row.items() if v} for row in jac]
    records = {}
    for name in ['CA', 'CB']:
        cb = name == 'CB'
        basis = [change(canonical(raw(f, 0, 0, cb))) for f in [h*h, h*x, x*x]]
        gbasis = [change(canonical(raw(0, f, 0, cb))) for f in [h*h, h*x, x*x]] if cb else []
        basis += [change(canonical(raw(0, h*h if cb else 1, 0, cb)))]
        basis += [change(canonical(raw(0, 0, H, cb))) for H in cubic_monomials]
        C = [[F(q) for q in row] for row in s.Matrix.hstack(*map(vec, basis)).tolist()]
        CG = [[F(q) for q in row] for row in s.Matrix.hstack(*map(vec, gbasis)).tolist()] if cb else None
        case = r['families'][name]
        extras = case['source_fiber_basis_indices']
        J = O.row_join(s.Matrix(C)[:, extras])
        assert J == s.Matrix(case['linear_J'])
        Ls = J[regular, :].inv()
        assert Ls*J[regular, :] == J[regular, :]*Ls == s.eye(28)
        left = [[F(q) for q in row] for row in Ls.tolist()]
        km = [[F(q) for q in row] for row in case['source_kernel_Lie_matrix']]
        km2 = [[sum(km[i][k]*km[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
        co = list(map(F, case['source_velocity']))
        inner = [[sumj([I[i][j], scale(PVAR, km[i][j]), scale(mul(PVAR, PVAR), km2[i][j]/2)]) for j in range(5)] for i in range(5)]
        innerinv = [[sumj([I[i][j], scale(PVAR, -km[i][j]), scale(mul(PVAR, PVAR), km2[i][j]/2)]) for j in range(5)] for i in range(5)]
        assert mm(inner, innerinv) == I
        base = [F(0)]*20
        base[0] = base[7] = -1
        base[11], base[16] = F(1), F(1, 3)
        def evaluate(source):
            M = [[Z for _ in range(5)] for _ in range(5)]
            for q, ij in zip(source[:19], oi):
                M[ij//5][ij % 5] = add(M[ij//5][ij % 5], q)
            M2 = mm(M, M)
            M3 = mm(M2, M)
            M4 = mm(M3, M)
            outer = ma(I, M)
            outerinv = ma(ma(ma(ma(I, M, -1), M2), M3, -1), M4)
            G, Gi = mm(outer, inner), mm(innerinv, outerinv)
            assert mm(G, Gi) == I
            alpha = [add(scale(ONE, a), scale(PVAR, -c)) for a, c in zip(base, co)]
            for q, j in zip(source[19:], extras):
                alpha[j] = add(alpha[j], q)
            if cb:
                g = [scale(mul(alpha[3], alpha[i]), -1) for i in range(3)]
                P = [sumj([*[scale(alpha[i], C[row][i]) for i in [0, 1, 2, *range(4, 20)] if C[row][i]],
                            *[scale(g[i], CG[row][i]) for i in range(3) if CG[row][i]]]) for row in range(150)]
            else:
                g = []
                P = [sumj(scale(a, c) for a, c in zip(alpha, row) if c) for row in C]
            quad = []
            for a, b in mons:
                rr = []
                for i, j in mons:
                    q = mul(Gi[a][i], Gi[b][j])
                    if i != j:
                        q = add(q, mul(Gi[a][j], Gi[b][i]))
                    rr.append(q)
                quad.append(rr)
            changed = [[sumj(mul(P[15*k+l], quad[l][m]) for l in range(15)) for m in range(15)] for k in range(10)]
            wedge = [[add(mul(G[i][a], G[j][b]), mul(G[i][b], G[j][a]), -1) for a, b in pairs] for i, j in pairs]
            q = [sumj(mul(wedge[i][k], changed[k][m]) for k in range(10)) for i in range(10) for m in range(15)]
            return q, alpha, g
        source = [Z]*28
        for degree in range(1, 5):
            q, _, _ = evaluate(source)
            rhs = [-q[i][degree] for i in regular]
            correction = [sum(a*b for a, b in zip(row, rhs)) for row in left]
            for i, value in enumerate(correction):
                v = list(source[i])
                v[degree] += value
                source[i] = tuple(v)
            dump_gz(out/f'{name}-stage{degree}.json.gz', dict(degree=degree, source28=enc(source)))
        q, alpha, g = evaluate(source)
        assert all(q[i][0] == point[i] for i in range(150))
        assert all(v[1] == 0 for v in q)
        assert all(q[i] == scale(ONE, point[i]) for i in regular)
        assert all(sumj(scale(v, a) for a, v in zip(row, q) if a) == Z for row in div)
        assert all(sumj(scale(mul(q[a], q[b]), c) for (a, b), c in row.items()) == Z for row in jac)
        fold_jet = q[fold]
        w2, c3, d4 = fold_jet[2:]
        assert w2 == (-4 if name == 'CA' else F(3, 8))
        a1 = [v[2]/w2 for v in q]
        beta1 = [(v[3]-e*c3)/w2 for v, e in zip(q, a1)]
        a2raw = [(v[4]-e*d4-b*c3)/(w2*w2) for v, e, b in zip(q, a1, beta1)]
        a2 = [(v[4]-v[2]*d4/w2+F(3, 2)*(v[2]*c3*c3/(w2*w2)-v[3]*c3/w2))/(w2*w2) for v in q]
        assert a2 == [a-c3*b/(2*w2*w2) for a, b in zip(a2raw, beta1)]
        assert a1 == [F(q) for q in even[0][:, 28]]
        assert beta1 == [F(q) for q in s.Matrix(first['CA_CB_full_firstjet_data'][name]['actual'])[:, 28]]
        tau = (F(0), F(-1), -c3/w2, -(c3/w2)**2, F(0))
        taupow = [ONE]
        for _ in range(4):
            taupow.append(mul(taupow[-1], tau))
        compose = lambda v: sumj(scale(a, c) for a, c in zip(taupow, v))
        assert compose(fold_jet) == fold_jet
        tracehalf = scale(add(PVAR, tau), F(1, 2))
        tracezero = add(PVAR, tracehalf, -1)
        for v, e, ee, b in zip(q, a1, a2, beta1):
            average = scale(add(v, compose(v)), F(1, 2))
            assert average == sumj([scale(ONE, v[0]), scale(fold_jet, e), scale(mul(fold_jet, fold_jet), ee)])
            assert scale(add(v, compose(v), -1), F(1, 2)) == scale(mul(fold_jet, tracezero), b)
        normalized = [v-sum(a*a2[i] for a, i in zip(row, odd)) for row, v in zip(h0, a2)]
        assert all(normalized[i] == 0 for i in selected)
        records[name] = dict(source28_parameter_GL_indices=oi, source28_parameter_fixed20_indices=extras,
            source28_IFT_coefficients=enc(source), actual_fixed20_parameter_jets=enc(alpha),
            full_CB_minus_ratio_times_binary_jets=enc(g),
            all150_physical_source_coefficients_through4=enc(q), regular28_inverse=enc(left),
            fold_row61_coefficients=enc(fold_jet), trace_involution_through3=enc(tau),
            all150_trace_even_linear_z_coefficients=enc(a1), all150_trace_odd_linear_z_coefficients=enc(beta1),
            all150_raw_even_quadratic_z_coefficients=enc(a2raw),
            all150_trace_even_quadratic_z_coefficients=enc(a2),
            all150_normalized_even_graph_quadratic_z_coefficients=enc(normalized))
    delta = [F(a)-F(b) for a, b in zip(records['CA']['all150_normalized_even_graph_quadratic_z_coefficients'],
                                       records['CB']['all150_normalized_even_graph_quadratic_z_coefficients'])]
    nonbasis = [i for i in range(150) if i not in selected]
    recon = s.Matrix(r['canonical_reconstruction_from126'])
    order = selected+nonbasis
    piv = list(recon[order, :].T.rref()[1])
    free126 = [order[i] for i in piv]
    assert len(free126) == 126 and free126[:39] == selected
    graph87 = free126[39:]
    nonzero = [i for i, v in enumerate(delta) if v]
    assert bool(nonzero) == any(delta[i] for i in graph87)
    evidence = dict(regular28_physical_rows=regular, fold_physical_row=fold, common39_physical_coordinates=selected,
        common_odd10_physical_rows=odd, full150_by10_common_odd_constant_graph_matrix=enc(h0), families=records,
        all150_normalized_CA_minus_CB_even_quadratic_coefficients=enc(delta),
        all111_nonbasis_physical_rows=nonbasis, independent87_physical_graph_rows=graph87,
        independent87_difference_coefficients=enc([delta[i] for i in graph87]))
    dump_gz(out/'certificate.json.gz', evidence)
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', reduction_path, first_path, tangent_path, finite_gate]
    result = dict(id=STEM, all_passed=True, single_source_degree=4, family_count=2,
        all150_genuine_physical_outputs_and_original25_div350_Jac_through4=True,
        actual28_target_coordinates_fixed_at_P0_including_scalar_baseline=True,
        full_CB_nonlinear_minus_ratio_times_f_retained=True,
        common_actual_h0_and_normalized_even_graph_comparison=True,
        trace_even_a2_no_p5_needed_by_all150_zero_p_linear_coefficients=True,
        independent_trace_involution_composition_check_all150=True,
        full111_nonbasis_and87_independent_graph_coverage=True,
        nonzero_normalized_graph_difference_rows=nonzero,
        any_normalized_even_graph_difference=bool(nonzero),
        proposed_exact_closed_Noether_fiber_length=12 if nonzero else 13,
        application_requires_independent_source4_and_full_graph_gate=True,
        no_p5_or_higher_or_full29_parameter_even_expansion=True,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={p.name: sha(p) for p in [out/'certificate.json.gz', *[out/f'{name}-stage{i}.json.gz' for name in ['CA', 'CB'] for i in range(1, 5)]]},
        seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'nonzero_normalized_graph_difference_rows',
          'proposed_exact_closed_Noether_fiber_length', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())

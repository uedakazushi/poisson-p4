#!/usr/bin/env python3
"""One bounded bivariate physical N-source IFT and all150 cubic graph classes."""
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
STEM = 'AMBIENT_B18_N_R_DVR_CUBIC_GRAPH_20261103_v63'
MON = [(a, n-a) for n in range(4) for a in range(n, -1, -1)]
MI = {m: i for i, m in enumerate(MON)}
Z = (F(0),)*10
ONE = (F(1),)+(F(0),)*9
CONV = [(i, j, MI[(a+c, b+d)]) for i, (a, b) in enumerate(MON)
        for j, (c, d) in enumerate(MON) if a+b+c+d <= 3]


def add(a, b, c=1):
    return tuple(x+c*y for x, y in zip(a, b))


def mul(a, b):
    if a == Z or b == Z:
        return Z
    out = [F(0)]*10
    for i, j, k in CONV:
        if a[i] and b[j]:
            out[k] += a[i]*b[j]
    return tuple(out)


def scale(a, c):
    return tuple(c*x for x in a) if c else Z


def sumj(values):
    out = Z
    for a in values:
        out = add(out, a)
    return out


def mm(A, B):
    return [[sumj(mul(A[i][k], B[k][j]) for k in range(len(B)))
             for j in range(len(B[0]))] for i in range(len(A))]


def ma(A, B, c=1):
    return [[add(a, b, c) for a, b in zip(ar, br)] for ar, br in zip(A, B)]


def read(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def encode(A):
    if isinstance(A, (list, tuple)):
        return [encode(a) for a in A]
    return str(A)


def main(out):
    start = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    out.mkdir(parents=True, exist_ok=False)
    reduction_path = ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json'
    first_path = ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz'
    curve_path = ROOT/'reports/AMBIENT_B18_R_CURVE_TANGENT_20261103_v63/certificate.json.gz'
    second_path = ROOT/'reports/AMBIENT_B18_N_R_SECOND_CONTACT_20261103_v63/report.json'
    r, first, curve, second = map(read, [reduction_path, first_path, curve_path, second_path])
    assert second['all_passed'] and second['whole_N_ideal_restriction_to_R_curve_has_order_atleast3']
    h, x, y, z, w, tt = s.symbols('h x y z w t')
    X = s.Matrix([h, x, y, z, w])
    pairs = list(itertools.combinations(range(5), 2))
    mons = list(itertools.combinations_with_replacement(range(5), 2))
    def vec(P):
        return s.Matrix([s.Poly(P[i, j], *X).coeff_monomial(X[a]*X[b]) for i, j in pairs for a, b in mons])
    def canonical(P):
        d = s.Matrix([sum(s.diff(P[i, j], X[i]) for i in range(5)) for j in range(5)])
        return (P-(X*d.T-d*X.T)/5).applyfunc(s.expand)
    def raw(f, g, H):
        P = s.zeros(5)
        P[1, 2], P[3, 4] = f, g*h*h
        P[2, 3], P[2, 4] = -s.diff(H, w), s.diff(H, z)
        return P-P.T
    sw = s.eye(5)
    sw[0, 0] = sw[1, 1] = 0
    sw[0, 1] = sw[1, 0] = 1
    change = lambda P: (sw*P.subs({h: x, x: h}, simultaneous=True)*sw.T).applyfunc(s.expand)
    cubics = [s.prod((h, x, z, w)[i] for i in m) for m in itertools.combinations_with_replacement(range(4), 3) if any(i >= 2 for i in m)]
    basis = [change(canonical(raw(f, 0, 0))) for f in [h*h, h*x, x*x]]
    basis += [change(canonical(raw(0, 1, 0)))]+[change(canonical(raw(0, 0, H))) for H in cubics]
    C = [[F(q) for q in row] for row in s.Matrix.hstack(*map(vec, basis)).tolist()]
    ca = r['families']['CA']
    extras = [i for i in ca['source_fiber_basis_indices'] if i != 2]
    oi = r['orbit_GL_generator_indices']
    co = list(map(F, ca['source_velocity']))
    km = [[F(q) for q in row] for row in ca['source_kernel_Lie_matrix']]
    km2 = [[sum(km[i][k]*km[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
    I = [[ONE if i == j else Z for j in range(5)] for i in range(5)]
    p = tuple(F(int(m == (0, 1))) for m in MON)
    p2 = mul(p, p)
    inner = [[sumj([I[i][j], scale(p, km[i][j]), scale(p2, km2[i][j]/2)]) for j in range(5)] for i in range(5)]
    innerinv = [[sumj([I[i][j], scale(p, -km[i][j]), scale(p2, km2[i][j]/2)]) for j in range(5)] for i in range(5)]
    assert mm(inner, innerinv) == I
    base = [F(0)]*20
    base[0] = base[7] = -1
    base[11], base[16] = F(1), F(1, 3)
    def evaluate(source):
        M = [[Z for j in range(5)] for i in range(5)]
        for q, index in zip(source[:19], oi):
            M[index//5][index % 5] = add(M[index//5][index % 5], q)
        M2, M3 = mm(M, M), None
        M3 = mm(M2, M)
        outer, outerinv = ma(I, M), ma(ma(ma(I, M, -1), M2), M3, -1)
        G, Gi = mm(outer, inner), mm(innerinv, outerinv)
        assert mm(G, Gi) == I
        alpha = [add(scale(ONE, q), scale(p, -c)) for q, c in zip(base, co)]
        for q, j in zip(source[19:], extras):
            alpha[j] = add(alpha[j], q)
        assert alpha[2] == Z
        P = [sumj(scale(a, c) for a, c in zip(alpha, row) if c) for row in C]
        # Pullback of each quadratic monomial by the actual inverse matrix.
        quad = []
        for a, b in mons:
            row = []
            for i, j in mons:
                q = mul(Gi[a][i], Gi[b][j])
                if i != j:
                    q = add(q, mul(Gi[a][j], Gi[b][i]))
                row.append(q)
            quad.append(row)
        changed = [[sumj(mul(P[15*k+l], quad[l][m]) for l in range(15)) for m in range(15)] for k in range(10)]
        wedge = [[add(mul(G[i][a], G[j][b]), mul(G[i][b], G[j][a]), -1) for a, b in pairs] for i, j in pairs]
        return [sumj(mul(wedge[i][k], changed[k][m]) for k in range(10)) for i in range(10) for m in range(15)]
    J = s.Matrix(ca['linear_J'])[:, [i for i in range(28) if i != 20]]
    regular = first['N_own_linear27_rows']
    fold = first['N_own_branch_row']
    left = [[F(q) for q in row] for row in J[regular, :].inv().tolist()]
    target = []
    for q in curve['physical_curve150']:
        poly = s.Poly(q, tt)
        target.append(tuple(F(poly.nth(a)) if b == 0 else F(0) for a, b in MON))
    source = [Z for _ in range(27)]
    stage = []
    for degree in range(1, 4):
        q = evaluate(source)
        for m, exps in enumerate(MON):
            if sum(exps) != degree:
                continue
            rhs = [target[i][m]-q[i][m] for i in regular]
            correction = [sum(a*b for a, b in zip(row, rhs)) for row in left]
            for i, a in enumerate(correction):
                v = list(source[i])
                v[m] += a
                source[i] = tuple(v)
        stage.append(dict(degree=degree, source27=encode(source), seconds=time.monotonic()-start))
        with gzip.open(out/f'stage{degree}.json.gz', 'wt') as f:
            json.dump(stage[-1], f)
    q = evaluate(source)
    assert all(q[i] == target[i] for i in regular)
    previous = read(second_path.parent/'certificate.json.gz')
    sparse = lambda v: {str(i): str(a) for i, a in enumerate(v) if a}
    assert sparse([target[i][MI[(2, 0)]]-q[i][MI[(2, 0)]] for i in range(150)]) == previous['source27_projected_delta']
    assert sparse([v[MI[(1, 1)]] for v in q]) == previous['source27_projected_mixed']
    assert sparse([v[MI[(0, 2)]] for v in q]) == previous['source27_projected_kernel_square']
    f = add(q[fold], target[fold], -1)
    assert f[MI[(0, 2)]] == -3 and f[MI[(1, 1)]] == 48
    assert all(f[i] == 0 for i, mon in enumerate(MON) if sum(mon) <= 2 and mon not in [(0, 2), (1, 1)])
    a30, a21, a12, a03 = [f[MI[m]] for m in [(3, 0), (2, 1), (1, 2), (0, 3)]]
    A, B = a30/3, (a21+16*a12+256*a03)/3
    wp = add(add(p2, scale(mul(p, tuple(F(int(m == (1, 0))) for m in MON)), -16)),
             tuple(-A if m == (3, 0) else -B if m == (2, 1) else F(0) for m in MON))
    unit = tuple(-3 if m == (0, 0) else a12+16*a03 if m == (1, 0) else a03 if m == (0, 1) else F(0) for m in MON)
    assert mul(unit, wp) == f
    def reduce(v):
        aa, bb = [F(0)]*4, [F(0)]*3
        for value, (i, j) in zip(v, MON):
            if not value:
                continue
            if j == 0:
                aa[i] += value
            elif j == 1:
                bb[i] += value
            elif (i, j) == (0, 2):
                bb[1] += 16*value
                aa[3] += A*value
                bb[2] += B*value
            elif (i, j) == (1, 2):
                bb[2] += 16*value
            elif (i, j) == (0, 3):
                bb[2] += 256*value
            else:
                raise AssertionError((i, j))
        return aa, bb
    reduced = [reduce(v) for v in q]
    aa, bb = [v[0] for v in reduced], [v[1] for v in reduced]
    assert reduce(f) == ([F(0)]*4, [F(0)]*3)
    for i in regular+[fold]:
        assert aa[i] == [target[i][MI[(j, 0)]] for j in range(4)] and bb[i] == [F(0)]*3
    assert all(b[0] == 0 for b in bb)
    nrows = curve['N37_coordinate_rows']
    odd = nrows[28:]
    assert odd == [18, 19, 28, 29, 60, 61, 67, 46, 72]
    TN = [[F(a) for a in row] for row in curve['full150_by37_actual_N_normalized_tangent']]
    h0 = [row[28:] for row in TN]
    ztarget = [[target[i][MI[(j, 0)]]-aa[i][j] for j in range(4)] for i in range(150)]
    assert all(z[0] == z[1] == 0 for z in ztarget)
    l = [bb[i][1] for i in odd]
    nonzero = next(i for i, a in enumerate(l) if a)
    mean_raw = ztarget[odd[nonzero]][2]/l[nonzero]
    assert mean_raw == F(300, 23)
    assert all(bb[i][1] == 0 and ztarget[i][2] == 0 for i in odd[-2:])
    assert all(ztarget[i][2] == mean_raw*bb[i][1] for i in range(150))
    assert all(bb[j][1] == sum(a*bb[i][1] for a, i in zip(h0[j], odd)) for j in range(150))
    cubic = []
    bres = []
    for j in range(150):
        r2 = bb[j][2]-sum(a*bb[i][2] for a, i in zip(h0[j], odd))
        d3 = ztarget[j][3]-sum(a*ztarget[i][3] for a, i in zip(h0[j], odd))
        bres.append(r2)
        cubic.append(d3-mean_raw*r2)
    # Independent trace-zero rewriting of the same formula.
    # wp=p²-(16T+B T²)p-A T³, so trace(p)/2=8T+B T²/2.
    trace_aa = [[a[0], a[1], a[2]+8*b[1], a[3]+8*b[2]+B*b[1]/2] for a, b in zip(aa, bb)]
    trace_z = [[target[i][MI[(j, 0)]]-trace_aa[i][j] for j in range(4)] for i in range(150)]
    mean_trace = mean_raw-8
    assert mean_trace == F(116, 23)
    trace_cubic = [trace_z[j][3]-sum(a*trace_z[i][3] for a, i in zip(h0[j], odd))-mean_trace*bres[j] for j in range(150)]
    assert cubic == trace_cubic
    assert all(cubic[i] == 0 for i in nrows)
    nonbasis = [i for i in range(150) if i not in nrows]
    discrepancy = mean_trace**2-64
    assert discrepancy == -F(20400, 529)
    evidence = dict(
        bivariate_monomial_order=MON, source27_parameter_GL_indices=oi,
        source27_parameter_fixed_N_fiber_indices=extras,
        source27_IFT_coefficients=encode(source), regular27_physical_rows=regular,
        regular27_inverse=encode(left), fold_physical_row=fold,
        full150_physical_source_coefficients_afterIFT=encode(q),
        full150_actual_R_target_coefficients=encode(target),
        raw_fold_equation_coefficients=encode(f), Weierstrass_monic_coefficients=encode(wp),
        Weierstrass_unit_coefficients=encode(unit),
        full150_remainders_even_a=encode(aa), full150_remainders_odd_b=encode(bb),
        actual_N37_physical_coordinate_rows=nrows, actual_odd9_physical_rows=odd,
        actual_full150_by9_minimal_odd_class_constants=encode(h0),
        raw_odd_mean_second_coefficient=str(mean_raw), trace_zero_odd_mean_second_coefficient=str(mean_trace),
        full150_trace_zero_even_a=encode(trace_aa),
        full150_b_minus_h0_bminimal_second_coefficients=encode(bres),
        full150_linear_graph_cubic_restrictions=encode(cubic),
        all113_nonbasis_physical_rows=nonbasis,
        all113_nonbasis_cubic_restrictions=encode([cubic[i] for i in nonbasis]),
        trace_zero_discriminant_second_coefficient='64', square_discrepancy=str(discrepancy),
    )
    with (out/'certificate.json.gz').open('wb') as f0:
        with gzip.GzipFile(fileobj=f0, mode='wb', filename='', mtime=0) as g:
            g.write((json.dumps(evidence, separators=(',', ':'))+'\n').encode())
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', reduction_path, first_path,
             curve_path, curve_path.parent/'report.json', second_path, second_path.parent/'certificate.json.gz']
    result = dict(
        id=STEM, all_passed=True, source='genuine actual N19 fixed space with GL19 and kernel1',
        total_bivariate_degree=3, physical_coordinate_count=150,
        exact_IFT_regular_target_count=27, exact_fold_Weierstrass_unit_identity_through3=True,
        full150_even_odd_remainder_coefficients_saved=True,
        actual_minimal_N_coordinate_count=37, actual_minimal_odd_generator_count=9,
        constant_division_coefficients_from_accepted_actual_N37=True,
        raw_mean_second_coefficient=str(mean_raw), trace_zero_mean_second_coefficient=str(mean_trace),
        raw_and_trace_zero_graph_cubic_formulas_equal_all150=True,
        all_nonbasis_physical_coordinate_count=113,
        nonzero_cubic_graph_rows=[i for i, q in enumerate(cubic) if q],
        all150_linear_graph_cubic_restrictions_zero=all(q == 0 for q in cubic),
        exact_k_h_application_requires_separate_full_N81_h4_gate=True,
        no_higher_degree_or_multibase_expansion=True,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={p.name: sha(p) for p in [out/'certificate.json.gz',
                         *[out/f'stage{i}.json.gz' for i in range(1, 4)]]},
        seconds=time.monotonic()-start,
    )
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'nonzero_cubic_graph_rows',
          'all150_linear_graph_cubic_restrictions_zero', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())

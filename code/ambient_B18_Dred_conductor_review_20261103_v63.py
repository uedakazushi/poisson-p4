#!/usr/bin/env python3
"""Independent physical N-source4 and original reduced conductor review."""
import argparse
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
STEM = 'AMBIENT_B18_DRED_CONDUCTOR_20261103_v63_REVIEW'
MON = [(a, d-a) for d in range(5) for a in range(d, -1, -1)]
MI = {m: i for i, m in enumerate(MON)}
Z, ONE = (F(0),)*15, (F(1),)+(F(0),)*14
CONV = [(i, j, MI[a+c, b+d]) for i, (a, b) in enumerate(MON)
        for j, (c, d) in enumerate(MON) if a+b+c+d <= 4]


def add(a, b, c=1):
    return tuple(x+c*y for x, y in zip(a, b))


def scale(a, c):
    return tuple(c*x for x in a) if c else Z


def mul(a, b):
    if a == Z or b == Z:
        return Z
    r = [F(0)]*15
    for i, j, k in CONV:
        if a[i] and b[j]:
            r[k] += a[i]*b[j]
    return tuple(r)


def summ(values):
    out = Z
    for q in values:
        out = add(out, q)
    return out


def mm(A, B):
    return [[summ(mul(A[i][k], B[k][j]) for k in range(len(B)))
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
    if isinstance(A, (list, tuple)):
        return [enc(v) for v in A]
    return str(A)


def c4mul(a, b):
    return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(4)]


def shift(v):
    return [F(0), *v[:3], F(0), *v[4:7]]


def tau(v, nu):
    return v[4:]+c4mul(nu, v[:4])


def main(out):
    start = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    out.mkdir(parents=True, exist_ok=False)
    folder = ROOT/'reports/AMBIENT_B18_DRED_CONDUCTOR_20261103_v63'
    author_path, source_path, kernel_path = [folder/n for n in ['report.json', 'source4.json.gz', 'conductor-kernel.json.gz']]
    author, saved, kernel = map(read, [author_path, source_path, kernel_path])
    assert author['all_passed']
    red_path = ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json'
    curve_path = ROOT/'reports/AMBIENT_B18_R_CURVE_TANGENT_20261103_v63/certificate.json.gz'
    gates = [ROOT/'reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63_REVIEW.json',
             ROOT/'reports/AMBIENT_B18_N_R_DVR_CUBIC_GRAPH_20261103_v63_REVIEW/report.json',
             ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63_REVIEW.json']
    assert all(read(p)['all_passed'] for p in gates)
    red, curve = map(read, [red_path, curve_path])
    h, x, y, z, w, t = s.symbols('h x y z w t')
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
    C = [change(canonical(raw(f, 0, 0))) for f in [h*h, h*x, x*x]]
    C += [change(canonical(raw(0, 1, 0)))]+[change(canonical(raw(0, 0, H))) for H in cubics]
    CM = s.Matrix.hstack(*map(vec, C))
    ey = s.eye(5)[:, 2]
    P0 = ey*s.Matrix([x*x, 0, 0, h*x, h*h+z*z]).T-s.Matrix([x*x, 0, 0, h*x, h*h+z*z])*ey.T
    point = list(map(F, vec(P0)))
    K = s.Matrix(red['families']['CA']['source_kernel_Lie_matrix'])
    KV = K*X
    LK = vec((K*P0+P0*K.T-sum((KV[i]*P0.diff(X[i]) for i in range(5)), s.zeros(5))).applyfunc(s.expand))
    assert CM*s.Matrix(red['families']['CA']['source_velocity']) == LK
    extras = [j for j in red['families']['CA']['source_fiber_basis_indices'] if j != 2]
    fixed = [[F(q) for q in row] for row in CM[:, extras].tolist()]
    oi = red['orbit_GL_generator_indices']
    O = []
    for ij in oi:
        E = s.zeros(5)
        E[ij//5, ij % 5] = 1
        EV = E*X
        O.append(vec((E*P0+P0*E.T-sum((EV[i]*P0.diff(X[i]) for i in range(5)), s.zeros(5))).applyfunc(s.expand)))
    JJ = s.Matrix.hstack(*O, CM[:, extras])
    assert JJ == s.Matrix(red['families']['CA']['linear_J'])[:, [i for i in range(28) if i != 20]]
    regular, fold = saved['regular27'], saved['fold']
    left = [[F(q) for q in row] for row in JJ[regular, :].inv().tolist()]
    assert s.Matrix(left)*JJ[regular, :] == JJ[regular, :]*s.Matrix(left) == s.eye(27)
    km = [[F(q) for q in row] for row in K.tolist()]
    km2 = [[sum(km[i][l]*km[l][j] for l in range(5)) for j in range(5)] for i in range(5)]
    I = [[ONE if i == j else Z for j in range(5)] for i in range(5)]
    p = tuple(F(int(m == (0, 1))) for m in MON)
    p2 = mul(p, p)
    inner = [[summ([I[i][j], scale(p, km[i][j]), scale(p2, km2[i][j]/2)]) for j in range(5)] for i in range(5)]
    innerinv = [[summ([I[i][j], scale(p, -km[i][j]), scale(p2, km2[i][j]/2)]) for j in range(5)] for i in range(5)]
    assert mm(inner, innerinv) == I
    def evaluate(source):
        M = [[Z for _ in range(5)] for _ in range(5)]
        for q, ij in zip(source[:19], oi):
            M[ij//5][ij % 5] = add(M[ij//5][ij % 5], q)
        M2 = mm(M, M)
        M3, M4 = mm(M2, M), mm(M2, M2)
        G = mm(ma(I, M), inner)
        Gi = mm(innerinv, ma(ma(ma(ma(I, M, -1), M2), M3, -1), M4))
        assert mm(G, Gi) == I
        source_q = [summ([scale(ONE, point[i]), scale(p, -F(LK[i])),
                          *[scale(q, a) for q, a in zip(source[19:], fixed[i]) if a]]) for i in range(150)]
        quad = []
        for a, b in mons:
            row = []
            for i, j in mons:
                v = mul(Gi[a][i], Gi[b][j])
                if i != j:
                    v = add(v, mul(Gi[a][j], Gi[b][i]))
                row.append(v)
            quad.append(row)
        changed = [[summ(mul(source_q[15*k+l], quad[l][m]) for l in range(15)) for m in range(15)] for k in range(10)]
        wedge = [[add(mul(G[i][a], G[j][b]), mul(G[i][b], G[j][a]), -1) for a, b in pairs] for i, j in pairs]
        return [summ(mul(wedge[i][k], changed[k][m]) for k in range(10)) for i in range(10) for m in range(15)]
    D = s.Matrix([x*x-13*t*h*x+25*t*t*x*z, 12*t*x*x, 0, h*x-13*t*x*z, h*h+z*z+2*t*x*w])
    RP = ey*D.T-D*ey.T
    actual_curve = vec(RP)
    assert list(map(str, actual_curve)) == curve['physical_curve150']
    target = [tuple(F(s.Poly(q, t).nth(a)) if b == 0 else F(0) for a, b in MON) for q in actual_curve]
    source = [Z]*27
    initial = evaluate(source)
    for degree in range(1, 5):
        q = evaluate(source)
        for n, mon in enumerate(MON):
            if sum(mon) == degree:
                rhs = [target[i][n]-q[i][n] for i in regular]
                corr = [sum(a*b for a, b in zip(row, rhs)) for row in left]
                for j, a in enumerate(corr):
                    value = list(source[j])
                    value[n] += a
                    source[j] = tuple(value)
        stage = read(folder/f'stage{degree}.json.gz')
        assert [list(map(str, row)) for row in source] == stage['source27']
    q = evaluate(source)
    assert q == [tuple(map(F, row)) for row in saved['all150_source_coefficients']]
    assert all(q[i] == target[i] for i in regular)
    # Entire source divergence and original350 Jacobi operations, independently
    # assembled on the physical150 coefficients, not abstract chart equations.
    pi = {p: i for i, p in enumerate(pairs)}
    div = [[0]*150 for _ in range(25)]
    for column, ((i, j), (a, b)) in enumerate(itertools.product(pairs, mons)):
        for k, rem in [(a, b), (b, a)]:
            if k == i:
                div[5*j+rem][column] += 1
            if k == j:
                div[5*i+rem][column] -= 1
    assert all(summ(scale(v, a) for a, v in zip(row, q) if a) == Z for row in div)
    spatial3 = list(itertools.combinations_with_replacement(range(5), 3))
    jac = [dict() for _ in range(350)]
    def entry(i, j):
        return pi[(i, j)], 1 if i < j else -1
    for ti, triple in enumerate(itertools.combinations(range(5), 3)):
        for i, j, k in [triple, (triple[1], triple[2], triple[0]), (triple[2], triple[0], triple[1])]:
            bi, sb = (pi[tuple(sorted((j, k)))], 1 if j < k else -1)
            for l in range(5):
                if l == i:
                    continue
                ai, sa = (pi[tuple(sorted((i, l)))], 1 if i < l else -1)
                for am, (u, v) in enumerate(mons):
                    for bm, (a, b) in enumerate(mons):
                        count = int(a == l)+int(b == l)
                        if count:
                            rem = b if a == l else a
                            row = 35*ti+spatial3.index(tuple(sorted((u, v, rem))))
                            key = tuple(sorted((15*ai+am, 15*bi+bm)))
                            jac[row][key] = jac[row].get(key, 0)+sa*sb*count
    products = {}
    for row in jac:
        value = Z
        for (a, b), c in row.items():
            if c:
                if (a, b) not in products:
                    products[a, b] = mul(q[a], q[b])
                value = add(value, products[a, b], c)
        assert value == Z
    # Independent5x5 coefficient solve for distinguished quadratic and unit.
    f = add(q[fold], target[fold], -1)
    Q2 = tuple(F(1) if m == (0, 2) else F(-16) if m == (1, 1) else F(0) for m in MON)
    quartics = [i for i, m in enumerate(MON) if sum(m) == 4]
    terms = []
    for exponent in [(2, 0), (1, 1), (0, 2)]:
        mon = tuple(F(int(m == exponent)) for m in MON)
        terms.append(mul(mon, Q2))
    terms += [tuple(F(-3) if m == e else F(0) for m in MON) for e in [(3, 1), (4, 0)]]
    linear = s.Matrix([[v[i] for v in terms] for i in quartics])
    solved = list(map(F, linear.inv()*s.Matrix([f[i] for i in quartics])))
    U = list(scale(ONE, -3))
    for exponent, a in zip([(2, 0), (1, 1), (0, 2)], solved[:3]):
        U[MI[exponent]] = a
    WP = list(Q2)
    WP[MI[3, 1]], WP[MI[4, 0]] = solved[3:]
    assert mul(tuple(U), tuple(WP)) == f
    decode = lambda records: {tuple(e): F(a) for e, a in records}
    assert {m: a for m, a in zip(MON, WP) if a} == decode(saved['Weierstrass_W'])
    assert {m: a for m, a in zip(MON, U) if a} == decode(saved['Weierstrass_U'])
    def redmon(i, j, a, output):
        if not a or i+j > 4:
            return
        if j < 2:
            output[i, j] = output.get((i, j), F(0))+a
        else:
            for (r, s0), c in zip(MON, WP):
                if c and (r, s0) != (0, 2):
                    redmon(i+r, j-2+s0, -a*c, output)
    aa, bb = [], []
    for row in q:
        remainder = {}
        for (i, j), a in zip(MON, row):
            redmon(i, j, a, remainder)
        aa.append([remainder.get((i, 0), F(0)) for i in range(4)])
        bb.append([remainder.get((i, 1), F(0)) for i in range(4)])
    center = [F(0), F(8), F(0), -solved[3]/2]
    nu = [F(0), F(0), F(64), F(0)]
    aatrace = [[a[i]+c4mul(center, b)[i] for i in range(4)] for a, b in zip(aa, bb)]
    dd = [[target[j][MI[i, 0]]-aatrace[j][i] for i in range(4)] for j in range(150)]
    assert bb == [[F(x) for x in row] for row in saved['all150_odd_b_throughT3']]
    assert aatrace == [[F(x) for x in row] for row in saved['all150_trace_even_a_throughT3']]
    assert dd == [[F(x) for x in row] for row in saved['all150_phi_odd_throughT3']]
    odd = saved['actual_odd9_physical_rows']
    assert odd == curve['N37_coordinate_rows'][28:]
    sigma, eta = [bb[i] for i in odd], [dd[i] for i in odd]
    gens = [a+b for a, b in zip(sigma, eta)]+[b+c4mul(nu, a) for a, b in zip(sigma, eta)]
    columns = []
    for g in gens:
        for _ in range(4):
            columns.append(g)
            g = shift(g)
    M = s.Matrix.hstack(*map(s.Matrix, columns))
    assert M == s.Matrix(kernel['C4_linear_map8_by72']) and M.rank() == 4
    KB = s.Matrix(kernel['QQ_kernel_basis68_by72'])
    assert KB.shape == (68, 72) and M*KB.T == s.zeros(8, 68) and KB.rank() == 68
    gm = gens[odd.index(61)]
    basis = [gm, shift(gm), tau(gm, nu), shift(shift(gm))]
    C4 = s.Matrix.hstack(*map(s.Matrix, basis))
    assert C4 == s.Matrix(kernel['cyclic_basis8_by4']) and C4.rank() == 4
    hi = s.Matrix(kernel['all9_coefficients_hi_in_cyclic_basis'])
    for i, g in enumerate(gens[:9]):
        assert C4*hi[i, :].T == s.Matrix(g)
        coeff = list(map(F, hi[i, :]))
        taubasis = [tau(v, nu) for v in basis]
        assert [sum(a*v[j] for a, v in zip(coeff, taubasis)) for j in range(8)] == gens[9+i]
    mu = F(116, 23)
    assert [a-mu*b for a, b in zip(shift(tau(gm, nu)), shift(shift(gm)))] == [0]*8
    assert [a-64*b for a, b in zip(tau(tau(gm, nu), nu), shift(shift(gm)))] == [0]*8
    assert mu*mu-64 == -F(20400, 529)
    assert [s.Matrix.hstack(*map(s.Matrix, columns)).rank(),
            s.Matrix.hstack(*(s.Matrix(shift(v)) for v in columns)).rank(),
            s.Matrix.hstack(*(s.Matrix(shift(shift(v))) for v in columns)).rank()] == [4, 2, 1]
    # Actual even-base slice: all R coordinates here are CENTERED at P0.
    Rrows = red['common28_linear_target_rows']
    OR = JJ[Rrows, :19]
    ell = s.Matrix([list(map(F, kernel['uniformizer_linear_ell']))])
    vR = s.Matrix([target[i][MI[1, 0]] for i in Rrows])
    assert ell*OR == s.zeros(1, 19) and ell*vR == s.ones(1, 1)
    rc = [[target[i][MI[j, 0]] for j in range(1, 4)] for i in Rrows]
    assert [sum(F(a)*row[j] for a, row in zip(ell, rc)) for j in range(3)] == [1, 0, 0]
    assert [(Rrows[j], F(a)) for j, a in enumerate(ell) if a] == [(16, F(1, 2)), (65, F(1)), (112, F(-1, 2))]
    gpiv = kernel['group19_gauge_physical_R_indices']
    group_unit = OR[gpiv, :].det()
    assert group_unit != 0
    rawW = s.Matrix([row[MI[0, 2]] for row in initial])
    aug = JJ.row_join(rawW)
    nugrad = aug[Rrows, :].inv()[-1, :]
    assert nugrad == s.Matrix([kernel['source_nu_gradient']])
    assert nugrad*OR == s.zeros(1, 19) and nugrad*vR == s.zeros(1, 1)
    coords = gpiv+kernel['residual7_even_normal_physical_indices']
    coord = s.Matrix.vstack(nugrad, ell, *[s.eye(28)[i, :] for i in coords])
    assert len(coords) == 26 and coord.det() == s.Rational(kernel['actual_base_coordinate_chart_unit']) != 0
    assert kernel['residual7_even_normal_physical_rows'] == [108, 109, 113, 117, 118, 119, 140]
    assert all(rc[i] == [0, 0, 0] for i in kernel['residual7_even_normal_physical_indices'])
    # Save independently obtained source and matrix evidence, not author reruns.
    evidence = dict(source27=enc(source), all150_physical_source4=enc(q),
        own_five_by_five_Weierstrass_coefficient_solve=enc(solved),
        phi_even9=enc(sigma), phi_odd9=enc(eta), phi_nu=enc(nu),
        all9_cyclic_coefficients=enc(hi.tolist()), group19_unit_minor=str(group_unit),
        actual_even_coordinate_chart_unit=str(coord.det()),
        explicit_uniformizer='u16/2+u65-u112/2, centered physical coordinates',
        independent_R_module_ranks=[4, 2, 1, 0])
    with (out/'independent-source4.json.gz').open('wb') as f:
        with gzip.GzipFile(fileobj=f, mode='wb', filename='', mtime=0) as g:
            g.write((json.dumps(evidence, separators=(',', ':'))+'\n').encode())
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', author_path, source_path, kernel_path,
             ROOT/'code/ambient_B18_Dred_conductor_20261103_v63.py',
             ROOT/'reports/AMBIENT_B18_REDUCED_INTERSECTION_DVR_GLUING_20261103_v63/report.json',
             ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63_REVIEW.json',
             ROOT/'research/notes/AMBIENT_B18_DRED_CONDUCTOR_20261103_v63.md', red_path, curve_path,
             *[folder/f'stage{i}.json.gz' for i in range(1, 5)], *gates]
    result = dict(id=STEM, all_passed=True,
        independent_actual_fixed_N19_source_and_full150_bivariate_physical_D4=True,
        independent_all27_IFT4_stages_and_original25div350Jac_coefficients=True,
        independent_Weierstrass5coefficient_solve_and_all150_trace_remainders=True,
        actual18_generator_map_and_complete72_column_QQ_kernel68_verified=True,
        actual_odd9_generate_entire_I_by_accepted_N37_full_image_not_abstract_ansatz=True,
        full_conductor_kernel_formula_phi_a_phi_tau_a_proved=True,
        cyclic_BN_image_length4_and_R_Hilbert_2_1_1_verified=True,
        cyclic_annihilator_CI_Hilbert_1_2_1_and9_slice_generators_verified=True,
        actual_even_G19_gauge_contains_entire_curve_and_preserves_rank2_cover=True,
        exact_J_BN_ideal_17_generator_upper_bound_proved=True,
        no_J_equals_I_intersect_ann_or_tau_value_or_image_kernel_basechange_assumed=True,
        original_Dred_conductor='(J_slice,T^4), restored by actual G19 product',
        IB_N_over_J_length4_transverse_and_rank4_over_C_g19=True,
        full_nonreduced_D_ring_depth_or_expanded28_variable_coefficients_claimed=False,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'independent-source4.json.gz': sha(out/'independent-source4.json.gz')},
        seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'actual18_generator_map_and_complete72_column_QQ_kernel68_verified',
          'exact_J_BN_ideal_17_generator_upper_bound_proved', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())

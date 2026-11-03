#!/usr/bin/env python3
"""Independent cached-frame QQ(b) physical coefficient-array audit, no98s worker."""
import argparse
import gzip
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix
import ambient_CB_self_pair_actual_D_twojet_review_20261103_v63 as physical
from ambient_CA_origin_self_flag_4jet_review_20261103_v63 import reconstruct_CA

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_CB1_CA_RELATIVE_CI4_20261103_v63_REVIEW'
K, b, X = physical.K, physical.b, physical.X
h, x, y, z, w = X


def load(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)


def sha(path):
    H = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            H.update(block)
    return H.hexdigest()


def accepted(path):
    d = load(path)
    assert d['all_passed'] is True
    for name, digest in d.get('input_sha256', {}).items():
        assert sha(ROOT / name) == digest, name
    for name, digest in d.get('evidence_sha256', {}).items():
        assert sha(path.parent / name) == digest, name
    return d


def parse_matrix(rows):
    return [[physical.coerce(v) for v in row] for row in rows]


def mv(A, v):
    return [sum((c * v[j] for j, c in enumerate(row) if c), K.zero) for row in A]


def transpose(A):
    return list(map(list, zip(*A)))


def subvec(a, b0):
    return [x0 - y0 for x0, y0 in zip(a, b0)]


def determinant(A):
    n = len(A)
    value = K.zero
    for permutation in it.permutations(range(n)):
        inversions = sum(permutation[i] > permutation[j] for i in range(n) for j in range(i + 1, n))
        term = K.one if inversions % 2 == 0 else -K.one
        for i, j in enumerate(permutation):
            term *= A[i][j]
        value += term
    return value


def main(output):
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    signal.alarm(90)
    output.mkdir(parents=True, exist_ok=False)
    folder = ROOT / 'reports/AMBIENT_CB1_CA_RELATIVE_CI4_20261103_v63'
    pointfolder = ROOT / 'reports/AMBIENT_CB1_CA_TRANSVERSE_CI4_20261103_v63'
    accepted(folder / 'report.json')
    accepted(pointfolder / 'report.json')
    data, quotient, point = load(folder / 'generic-quadratic.json.gz'), load(folder / 'quotient-map.json'), load(pointfolder / 'certificate.json')
    gates = [ROOT / 'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_CB_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63_REVIEW/report.json']
    for path in gates:
        accepted(path)
    rawpath = ROOT / 'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63/certificate.json.gz'
    globalpath = ROOT / 'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json'
    raw, globalcb = load(rawpath), load(globalpath)
    def canon(P):
        E = s.Matrix(X)
        div = s.Matrix([sum(s.diff(P[i, j], X[i]) for i in range(5)) for j in range(5)])
        return (P - (E * div.T - div * E.T) / 5).applyfunc(s.expand)
    def rawcb(F, G, H):
        P = s.zeros(5)
        P[1, 2], P[3, 4] = F, G
        P[2, 3] = -s.diff(F, x) * z / 2 - s.diff(H, w)
        P[2, 4] = -s.diff(F, x) * w / 2 + s.diff(H, z)
        return P - P.T
    cubics = [s.prod((h, x, z, w)[j] for j in m) for m in it.combinations_with_replacement(range(4), 3) if any(j >= 2 for j in m)]
    binary = [h * h, h * x, x * x]
    basis = [physical.coefficients(canon(rawcb(F, 0, 0))) for F in binary]
    basis += [physical.coefficients(canon(rawcb(0, G, 0))) for G in binary]
    basis += [physical.coefficients(canon(rawcb(0, 0, H))) for H in cubics]
    assert basis == [[physical.coerce(c) for c in col] for col in globalcb['enlarged_canonical_basis']]
    fixed = [*basis[:3], *basis[6:]]
    flags = [ij for ij in it.product(range(5), repeat=2) if list(ij) not in globalcb['parabolic_GL_positions']]
    ops = [physical.tensor_factor_operator(i, j) for i, j in it.product(range(5), repeat=2)]
    def lie(i, j, q):
        wedge, coordinate = ops[5 * i + j]
        result = [K.zero] * 150
        for op, ip, coefficient in wedge:
            for mon in range(15):
                result[15 * op + mon] += coefficient * q[15 * ip + mon]
        for old, new, coefficient in coordinate:
            for pair in range(10):
                result[15 * pair + new] += coefficient * q[15 * pair + old]
        return result
    combine = lambda columns, values: [sum((column[i] * c for column, c in zip(columns, values) if c and column[i]), K.zero) for i in range(150)]
    datum = raw['actual_two_CB_sources'][0]
    G, Gi = parse_matrix(datum['G']), parse_matrix(datum['G_inverse'])
    back = physical.constant_map(Gi, G)
    F0 = s.sympify(datum['f'], locals=dict(zip(('b', 'h', 'x', 'y', 'z', 'w'), (b, *X))))
    H0 = s.sympify(datum['H'], locals=dict(zip(('b', 'h', 'x', 'y', 'z', 'w'), (b, *X))))
    values = [physical.coerce(s.Poly(F0, h, x).coeff_monomial(m)) for m in binary]
    values += [physical.coerce(s.Poly(H0, h, x, z, w).coeff_monomial(m)) for m in cubics]
    localpoint = combine(fixed, values)
    ratio = combine(basis[3:6], values[:3])
    point150 = [physical.coerce(v) for v in raw['full_polynomial_Q150']]
    assert back(localpoint) == point150
    JCBcolumns = [back(v) for v in [*fixed, ratio, *[lie(i, j, localpoint) for i, j in flags]]]
    capath = ROOT / 'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
    ca = load(capath)
    C = reconstruct_CA()
    assert C == s.Matrix([[s.Rational(c) for c in row] for row in ca['fixed_fiber_linear_basis150_by20']])
    CAcols = [[physical.coerce(C[i, j]) for i in range(150)] for j in range(20)]
    caflags = list(map(tuple, ca['flag_positions']))
    JCAcolumns = CAcols + [lie(i, j, point150) for i, j in caflags]
    A = transpose(JCBcolumns + [[-c for c in col] for col in JCAcolumns])
    normals = transpose(parse_matrix(data['normal58_by2']))
    frame = parse_matrix(data['frame58_by58'])
    pc, pr = data['implicit36_source_columns'], data['implicit36_target_rows']
    assert len(pc) == len(pr) == 36
    assert [row[56:] for row in frame] == transpose(normals)
    assert [row[20:56] for row in frame] == [[K.one if i == j else K.zero for j in pc] for i in range(58)]
    for direction in transpose(frame)[0:20] + normals:
        assert not any(mv(A, direction))
    inverse = parse_matrix(data['implicit36_inverse'])
    selected = [[A[i][j] for j in pc] for i in pr]
    for col in range(36):
        assert mv(selected, [inverse[i][col] for i in range(36)]) == [K.one if i == col else K.zero for i in range(36)]
    print('independent physical sources/frame/inverse bound', time.monotonic() - start, flush=True)

    def source_jet(direction, isCB):
        positions = flags if isCB else caflags
        Ag = [[K.zero for _ in range(5)] for _ in range(5)]
        for coefficient, (i, j) in zip(direction[20:], positions):
            Ag[i][j] += coefficient
        A2 = physical.matmul(Ag, Ag)
        exp = [[[K.one if i == j else K.zero, Ag[i][j], A2[i][j] / 2] for j in range(5)] for i in range(5)]
        invexp = [[[K.one if i == j else K.zero, -Ag[i][j], A2[i][j] / 2] for j in range(5)] for i in range(5)]
        if isCB:
            q1 = physical.add(combine(fixed, direction[:19]), physical.scale(ratio, direction[19]))
            q2 = physical.scale(combine(basis[3:6], direction[:3]), direction[19])
            result = physical.series_push2(exp, invexp, transpose([localpoint, q1, q2]))
            return [back(v) for v in result]
        return physical.series_push2(exp, invexp, transpose([point150, combine(CAcols, direction[:20]), [K.zero] * 150]))
    def second(direction):
        cb, ca0 = source_jet(direction[:29], True), source_jet(direction[29:], False)
        assert cb[0] == ca0[0] == point150 and cb[1] == ca0[1]
        return subvec(cb[2], ca0[2])
    assert normals[0][19] == normals[1][19] == K.zero
    q00, q11 = second(normals[0]), second(normals[1])
    q01 = subvec(subvec(second([a + c for a, c in zip(*normals)]), q00), q11)
    rawcols = [q00, q01, q11]
    assert transpose(rawcols) == parse_matrix(data['raw150_by3'])
    corrections = [mv(inverse, [v[i] for i in pr]) for v in rawcols]
    corrected = [[rawcols[j][i] - sum((A[i][col] * corrections[j][k] for k, col in enumerate(pc) if A[i][col] and corrections[j][k]), K.zero) for j in range(3)] for i in range(150)]
    assert corrected == parse_matrix(data['corrected150_by3'])
    selected_rows = [63, 68]
    rows = [corrected[i] for i in selected_rows]
    pivot = next((j, k) for j, k in it.combinations(range(3), 2) if rows[0][j] * rows[1][k] - rows[0][k] * rows[1][j])
    j, k = pivot
    det = rows[0][j] * rows[1][k] - rows[0][k] * rows[1][j]
    for v in corrected:
        c0 = (v[j] * rows[1][k] - v[k] * rows[1][j]) / det
        c1 = (rows[0][j] * v[k] - rows[0][k] * v[j]) / det
        assert v == [c0 * rows[0][i] + c1 * rows[1][i] for i in range(3)]
    cubic = [[*rows[0], K.zero], [K.zero, *rows[0]], [*rows[1], K.zero], [K.zero, *rows[1]]]
    cubicdet = determinant(cubic)
    assert cubicdet == physical.coerce(1024 * b ** 4 / (25 * (b * b - 3) ** 4))
    print('all150 QQ(b) raw/corrected two-normal quadratics independently PASS', time.monotonic() - start, flush=True)

    # Geometrically bind the CB meeting direction modulo the20 actual base directions.
    d2path = ROOT / 'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63/certificate.json'
    ucb = [physical.coerce(v) for v in load(d2path)['u1']]
    base_columns = [mv(transpose(JCBcolumns), direction[:29]) for direction in transpose(frame)[:20]]
    normal_physical = [mv(transpose(JCBcolumns), direction[:29]) for direction in normals]
    c0, c1 = map(physical.coerce, quotient['normal2'])
    assert c0 == physical.coerce(-b ** 3 + b) and c1 == physical.coerce(-2 * b * b)
    residue = [v - c0 * q - c1 * t0 for v, q, t0 in zip(mv(transpose(JCBcolumns), ucb), *normal_physical)]
    base_matrix = physical.matrix(base_columns)
    chosen = list(DomainMatrix.from_Matrix(base_matrix.subs(b, 1).T).convert_to(s.QQ).rref()[1])
    assert len(chosen) == 20
    baseinv = DomainMatrix.from_Matrix(base_matrix[chosen, :]).convert_to(K).inv().to_Matrix()
    coeff = mv(parse_matrix(baseinv.tolist()), [residue[i] for i in chosen])
    assert mv(transpose(base_columns), coeff) == residue
    eta2 = [c1 * c1, -2 * c0 * c1, c0 * c0]
    alpha, beta = physical.coerce(9 * b * b * (b * b - 3) ** 2 / 8), physical.coerce(15 * b * b * (b * b - 3) ** 2 / 4)
    assert eta2 == [alpha * rows[0][i] + beta * rows[1][i] for i in range(3)]
    assert [physical.coerce(b * b), physical.coerce(-2 * b), K.one] == [physical.coerce((b * b - 3) ** 2 / 24) * c for c in rows[0]]
    kb = physical.coerce(b)
    assert determinant([[c1, -c0], [kb, -K.one]]) == physical.coerce(-b * b * (b * b - 3))
    assert kb * c0 - c1 == physical.coerce(-b * b * (b * b - 3))
    # All cached arrays specialize to the independently bound point calculation.
    for key in ['frame58_by58', 'normal58_by2', 'implicit36_inverse', 'raw150_by3', 'corrected150_by3']:
        a = s.Matrix([[physical.expr(c).subs(b, 1) for c in row] for row in parse_matrix(data[key])])
        assert a == s.Matrix([[s.Rational(c) for c in row] for row in point[key]])
    assert DomainMatrix.from_Matrix(s.Matrix(data['frame58_by58']).subs(b, 1)).convert_to(s.QQ).det() == -2
    point_rows = [[physical.expr(c).subs(b, 1) for c in row] for row in rows]
    assert point_rows == [[6, -12, 6], [-s.Rational(23, 15), s.Rational(18, 5), -s.Rational(9, 5)]]
    assert physical.expr(cubicdet).subs(b, 1) == s.Rational(64, 25)
    evidence = dict(all22_CB_basis_and_actual_CA20_rebuilt=True,
        all150_actual_sources_and_full58_frame_and_36_inverse_bound=True,
        independent_coefficient_array_raw_and_corrected150_by3_QQb=True,
        completed98second_worker_not_repeated=True,
        full150_quadratic_rank=2, independent_rows=selected_rows,
        degree3_unit=str(physical.expr(cubicdet)), point_degree3_unit='64/25',
        actual_CB_meeting_direction_mod20_base_bound=True,
        eta_square_and_second_isotropic_square_QQb_identities=True,
        closed_ideal_equality_and_relative_R_free4_by_Nakayama_accepted=True,
        actual_quotient_diagram='R[u,v]/(u²,v²) -> R[v]/v², u->0',
        quotient_kernel_square_zero_R_free_rank=2,
        full_actual_D_Tor_nil_primary_or_label_swap_claim=False)
    (output / 'independent-evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
    paths = [Path(__file__), ROOT / f'research/notes/{STEM}.md',
        ROOT / 'code/ambient_CB_self_pair_actual_D_twojet_review_20261103_v63.py',
        ROOT / 'code/ambient_CA_origin_self_flag_4jet_review_20261103_v63.py',
        folder / 'report.json', folder / 'generic-quadratic.json.gz', folder / 'quotient-map.json',
        pointfolder / 'report.json', pointfolder / 'certificate.json', rawpath, globalpath, capath, d2path, *gates]
    report = dict(id=STEM, all_passed=True, producer_or98s_worker_rerun=False,
        independent_all150_QQb_physical_two_normal_jets_and_rank2=True,
        closed_CI4_and_full_relative_free4_accepted=True,
        actual_double_meeting_quotient_diagram_bound_and_accepted=True,
        full_actual_D_or_Tor_nil_primary_claim=False,
        independent_method='Accepted coefficient-array evaluator on new cached frame; literal inverse-coordinate source jets; complete150 QQ(b) projection/rank2; actual quotient tangent modulo20 base; relative Nakayama proof.',
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'independent-evidence.json': sha(output / 'independent-evidence.json')},
        seconds=time.monotonic() - start)
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    signal.alarm(0)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output)

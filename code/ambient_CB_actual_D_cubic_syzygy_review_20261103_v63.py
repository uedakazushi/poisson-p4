#!/usr/bin/env python3
"""Own sparse Lie D3, cached inverse constraints, and all150 cubic-socle audit."""
import argparse
from collections import Counter
import gzip
import itertools as it
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s
import ambient_CB1_CA_relative_CI4_review_20261103_v63 as rec

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_CB_ACTUAL_D_CUBIC_SYZYGY_20261103_v63_REVIEW'
K, b, physical = rec.K, rec.b, rec.physical
load, sha, accepted, matrix, transpose, mv = rec.load, rec.sha, rec.accepted, rec.parse_matrix, rec.transpose, rec.mv
PAIRS, MONS = list(it.combinations(range(5), 2)), list(it.combinations_with_replacement(range(5), 2))
PI, MI = {p: i for i, p in enumerate(PAIRS)}, {p: i for i, p in enumerate(MONS)}
ZERO = [K.zero] * 150


def va(*vectors):
    return [sum((v[i] for v in vectors), K.zero) for i in range(150)]


def vs(vector, coefficient):
    return [coefficient * x for x in vector]


def sparse_operator(i, j):
    """Differentiate exterior indices and inverse-coordinate monomials directly."""
    rows = [Counter() for _ in range(150)]
    for p, (k, l) in enumerate(PAIRS):
        for m, (u, v) in enumerate(MONS):
            source = 15 * p + m
            if k == j and i != l:
                rows[15 * PI[tuple(sorted((i, l)))] + m][source] += 1 if i < l else -1
            if l == j and k != i:
                rows[15 * PI[tuple(sorted((k, i)))] + m][source] += 1 if k < i else -1
            for old, other in ((u, v), (v, u)):
                if old == i:
                    rows[15 * p + MI[tuple(sorted((j, other)))]] [source] -= 1
    return [[(j0, c) for j0, c in row.items() if c] for row in rows]


def apply(op, vector):
    return [sum((c * vector[j] for j, c in row if vector[j]), K.zero) for row in op]


class Source:
    def __init__(self, kind, raw, globalcb, ca, C):
        self.kind = kind
        if kind == 'CA':
            self.basis = [[physical.coerce(C[i, j]) for i in range(150)] for j in range(20)]
            self.base = list(map(physical.coerce, ca['center_fixed20']))
            self.base[1] -= physical.coerce(2 * b)
            self.base[9] -= physical.coerce(b)
            self.flags = list(map(tuple, ca['flag_positions']))
            self.back = lambda x: x[:]
        else:
            elementary = [[physical.coerce(v) for v in row] for row in globalcb['enlarged_canonical_basis']]
            self.basis, self.g = [*elementary[:3], *elementary[6:]], elementary[3:6]
            item = raw['actual_two_CB_sources'][int(kind[-1]) - 1]
            X = physical.X
            h, x, y, z, w = X
            local = dict(zip(('b', 'h', 'x', 'y', 'z', 'w'), (b, *X)))
            f, H = s.sympify(item['f'], locals=local), s.sympify(item['H'], locals=local)
            cubics = [s.prod((h, x, z, w)[i] for i in mon) for mon in it.combinations_with_replacement(range(4), 3) if any(i >= 2 for i in mon)]
            self.base = [physical.coerce(s.Poly(f, h, x).coeff_monomial(m)) for m in (h * h, h * x, x * x)]
            self.base += [physical.coerce(s.Poly(H, h, x, z, w).coeff_monomial(m)) for m in cubics]
            self.flags = [p for p in it.product(range(5), repeat=2) if list(p) not in globalcb['parabolic_GL_positions']]
            self.back = physical.constant_map(matrix(item['G_inverse']), matrix(item['G']))
        self.ops = [sparse_operator(*p) for p in self.flags]
        self.point = self.back(self.combine(self.base))
        assert self.point == list(map(physical.coerce, raw['full_polynomial_Q150']))
        columns = self.basis[:]
        if kind != 'CA':
            columns.append(self.combine(self.base[:3], self.g))
        rawpoint = self.combine(self.base)
        columns += [apply(op, rawpoint) for op in self.ops]
        self.J = transpose([self.back(v) for v in columns])

    def combine(self, coefficients, basis=None):
        columns = self.basis if basis is None else basis
        return [sum((c * v[i] for c, v in zip(coefficients, columns) if c and v[i]), K.zero) for i in range(150)]

    def L(self, flags, vector):
        result = ZERO[:]
        for coefficient, op in zip(flags, self.ops):
            if coefficient:
                result = va(result, vs(apply(op, vector), coefficient))
        return result

    def fibre(self, jets, degree):
        count = 20 if self.kind == 'CA' else 19
        coefficients = [row[degree] + (self.base[i] if degree == 0 else K.zero) for i, row in enumerate(jets[:count])]
        result = self.combine(coefficients)
        if self.kind != 'CA':
            g = [sum((jets[19][k] * (jets[j][degree - k] + (self.base[j] if degree == k else K.zero)) for k in range(degree + 1)), K.zero) for j in range(3)]
            result = va(result, self.combine(g, self.g))
        return result

    def jet(self, jets, order):
        """Closed D3 exp(L1*t+L2*t²+L3*t³), not author's evaluator."""
        F = [self.fibre(jets, d) for d in range(order + 1)]
        A = [[row[d] for row in jets[20:]] for d in range(1, order + 1)]
        L1F0 = self.L(A[0], F[0])
        out = [F[0], va(F[1], L1F0)]
        L1F1, L2F0 = self.L(A[0], F[1]), self.L(A[1], F[0])
        L11F0 = self.L(A[0], L1F0)
        out.append(va(F[2], L1F1, L2F0, vs(L11F0, K.one / 2)))
        if order == 3:
            out.append(va(F[3], self.L(A[0], F[2]), self.L(A[1], F[1]), vs(self.L(A[0], L1F1), K.one / 2), self.L(A[2], F[0]), vs(va(self.L(A[0], L2F0), self.L(A[1], L1F0)), K.one / 2), vs(self.L(A[0], L11F0), K.one / 6)))
        return [self.back(v) for v in out]

    def mixed2(self, u, v):
        Fu, Fv = self.combine(u[:19] if self.kind != 'CA' else u[:20]), self.combine(v[:19] if self.kind != 'CA' else v[:20])
        extra = ZERO[:]
        if self.kind != 'CA':
            Fu = va(Fu, self.combine([u[19] * c for c in self.base[:3]], self.g))
            Fv = va(Fv, self.combine([v[19] * c for c in self.base[:3]], self.g))
            extra = self.combine([u[19] * v[i] + v[19] * u[i] for i in range(3)], self.g)
        F0 = self.combine(self.base)
        result = va(extra, self.L(u[20:], Fv), self.L(v[20:], Fu), vs(va(self.L(u[20:], self.L(v[20:], F0)), self.L(v[20:], self.L(u[20:], F0))), K.one / 2))
        return self.back(result)


def main(out):
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    signal.alarm(120)
    out.mkdir(parents=True, exist_ok=False)
    folder = ROOT / 'reports/AMBIENT_CB_ACTUAL_D_CUBIC_SYZYGY_20261103_v63'
    accepted(folder / 'report.json')
    D, setup = load(folder / 'generic-cubic.json.gz'), load(folder / 'generic-setup.json.gz')
    deps = [ROOT / 'reports/AMBIENT_CB1_CA_RELATIVE_CI4_20261103_v63_REVIEW/report.json', ROOT / 'reports/AMBIENT_CB_ACTUAL_D_QUADRATIC_GATE_20261103_v63/report.json', ROOT / 'reports/AMBIENT_CB_SELF_PAIR_COMMON_BASE_ETALE_20261103_v63/report.json']
    for path in deps:
        accepted(path)
    rawpath = ROOT / 'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63/certificate.json.gz'
    globalpath = ROOT / 'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json'
    capath = ROOT / 'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
    raw, globalcb, ca = map(load, (rawpath, globalpath, capath))
    C = rec.reconstruct_CA()
    sources = [Source(kind, raw, globalcb, ca, C) for kind in ('CB1', 'CB2', 'CA')]
    J = [source.J for source in sources]
    assert J == [matrix(q) for q in setup['source_Jacobians150_by29']]
    R, psel, j0, minimum, ri = [setup[key] for key in ('base29_rows', 'base20_indices', 'f_row', 'minimal8_linear_b_rows', 'rho_base_index')]
    cc = list(map(physical.coerce, setup['constant_odd_row_normalization']))
    inv = {mode: matrix(setup[f'inverse{mode}']) for mode in (58, 87)}
    def norm(v):
        return [q if i == j0 else q - cc[i] * v[j0] for i, q in enumerate(v)]
    def constraints(mode, values):
        if mode == 58:
            A, B = values
            an, bn = norm(A), norm(B)
            return [A[i] - B[i] for i in R] + [B[R[i]] for i in psel] + [(an[i] - bn[i]) / 2 for i in minimum] + [B[R[ri]]]
        A, B, C0 = values
        an, bn, cn = map(norm, values)
        return [A[i] - C0[i] for i in R] + [B[i] - C0[i] for i in R] + [C0[R[i]] for i in psel] + [an[i] - cn[i] for i in minimum[1:]] + [cn[j0] - bn[j0], C0[R[ri]]]
    for mode in (58, 87):
        columns = []
        for i in range(mode):
            values = [ZERO[:] for _ in range(mode // 29)]
            values[i // 29] = [row[i % 29] for row in J[i // 29]]
            columns.append(constraints(mode, values))
        N = transpose(columns)
        assert N == matrix(setup[f'Jacobian{mode}'])
        for i in range(mode):
            assert mv(N, [row[i] for row in inv[mode]]) == [K.one if j == i else K.zero for j in range(mode)]
    print('own sparse physical Jacobians and BOTH58/87 constraints/inverses PASS', time.monotonic() - start, flush=True)
    archived = []
    def curve(mode, slope):
        direction = [K.zero] * mode
        direction[49 if mode == 58 else 85] = K.one
        direction[57 if mode == 58 else 86] = slope
        jets = [[K.zero] * 4 for _ in range(mode)]
        for row, value in zip(jets, mv(inv[mode], direction)):
            row[1] = value
        for degree in (2, 3):
            values = [sources[k].jet(jets[29 * k:29 * (k + 1)], degree) for k in range(mode // 29)]
            correction = mv(inv[mode], [-q for q in constraints(mode, [v[degree] for v in values])])
            for row, q in zip(jets, correction):
                row[degree] += q
            for k, value in enumerate(values):
                value[degree] = va(value[degree], mv(J[k], correction[29 * k:29 * (k + 1)]))
        for degree in (1, 2, 3):
            assert constraints(mode, [v[degree] for v in values]) == (direction if degree == 1 else [K.zero] * mode)
        archived.append(dict(mode=mode, slope=str(slope), jets=jets, physical_values=values))
        return values
    slopes = list(map(physical.coerce, (0, 1, -1, 2)))
    samples = {mode: [curve(mode, slope) for slope in slopes] for mode in (58, 87)}
    print('eight own D3 slope curves/all source constraints PASS', time.monotonic() - start, flush=True)
    def homogeneous(values, degree):
        f0, f1, fm, f2 = values
        a0 = f0
        if degree == 0:
            coefficients = [a0]
        elif degree == 1:
            coefficients = [a0, f1 - f0]
        elif degree == 2:
            coefficients = [a0, (f1 - fm) / 2, (f1 + fm) / 2 - f0]
        else:
            a2, odd = (f1 + fm) / 2 - f0, (f1 - fm) / 2
            a3 = (f2 - f0 - 4 * a2 - 2 * odd) / 6
            coefficients = [a0, odd - a3, a2, a3]
        for slope, value in zip(slopes, values):
            assert sum((q * (slope ** j if j else K.one) for j, q in enumerate(coefficients)), K.zero) == value
        return {(degree - j, j): q for j, q in enumerate(coefficients) if q}
    def add(*polynomials):
        result = {}
        for p in polynomials:
            for word, q in p.items():
                result[word] = result.get(word, K.zero) + q
        return {word: q for word, q in result.items() if q}
    def scale(p, q):
        return {word: q * value for word, value in p.items() if q * value}
    def multiply(p, q, cap=3):
        result = {}
        for (i, j), v in p.items():
            for (k, l), u in q.items():
                if i + j + k + l <= cap:
                    word = (i + k, j + l)
                    result[word] = result.get(word, K.zero) + v * u
        return {word: q for word, q in result.items() if q}
    def poly(data):
        return {(i, j): physical.coerce(q) for i, j, q in data if physical.coerce(q)}
    one, rho = {(0, 0): K.one}, {(0, 1): K.one}
    def compose(p, x0, cap):
        result = {}
        for (i, j), q in p.items():
            term = one
            for _ in range(i):
                term = multiply(term, x0, cap)
            for _ in range(j):
                term = multiply(term, rho, cap)
            result = add(result, scale(term, q))
        return {word: q for word, q in result.items() if sum(word) <= cap}
    def coefficient(mode, sample, degree, kind):
        A, B = [norm(value[degree]) for value in sample[:2]]
        if kind == 'b':
            return [(q - r0) / 2 for q, r0 in zip(A, B)]
        C0 = norm(sample[2][degree])
        return [q - r0 for q, r0 in zip(C0, B)] if kind == 'c' else [q - r0 for q, r0 in zip(A, C0)]
    bplane = [add(*(homogeneous([coefficient(58, v, d, 'b')[j] for v in samples[58]], d) for d in (1, 2, 3))) for j in range(150)]
    bc = [add(*(homogeneous([coefficient(87, v, d, 'b')[j] for v in samples[87]], d) for d in (1, 2, 3))) for j in range(150)]
    cp = [add(*(homogeneous([coefficient(87, v, d + 1, 'c')[j] for v in samples[87]], d) for d in (0, 1, 2))) for j in range(150)]
    zp = [add(*(homogeneous([coefficient(87, v, d, 'z')[j] for v in samples[87]], d) for d in (1, 2, 3))) for j in range(150)]
    for values, key in ((bplane, 'conductor_plane_b150'), (bc, 'CA_plane_b150'), (cp, 'CA_plane_c150')):
        assert values == list(map(poly, D[key])), key
    # Fourteen mixed Hessians directly, no42 extra nonlinear curve callbacks.
    directions = []
    beta = [[K.zero] * 8 for _ in range(150)]
    for i in range(8):
        e = [K.zero] * 58
        e[49 + i] = K.one
        u = mv(inv[58], e)
        directions.append(u)
        values = [mv(J[k], u[k * 29:(k + 1) * 29]) for k in range(2)]
        bb = [(q - r0) / 2 for q, r0 in zip(*map(norm, values))]
        for j, q in enumerate(bb):
            beta[j][i] = q
    er = [K.zero] * 58
    er[57] = K.one
    urho = mv(inv[58], er)
    mixed = []
    for i in range(1, 8):
        pair = []
        for v in (directions[0], urho):
            values = [sources[k].mixed2(directions[i][k * 29:(k + 1) * 29], v[k * 29:(k + 1) * 29]) for k in range(2)]
            correction = mv(inv[58], [-q for q in constraints(58, values)])
            values = [va(q, mv(J[k], correction[k * 29:(k + 1) * 29])) for k, q in enumerate(values)]
            assert not any(constraints(58, values))
            pair.append([(q - r0) / 2 for q, r0 in zip(*map(norm, values))])
        mixed.append(pair)
    j8 = D['quadratic_generator_row']
    assert j0 == 18 and j8 == 63 and len(set(R)) == 29 and set(R).isdisjoint(minimum)
    assert bplane[j0] == {(1, 0): K.one} and all(not bplane[j] for j in minimum[1:])
    b8 = add(bplane[j8], {(1, 0): -beta[j8][0]})
    a8, a83 = b8[(0, 2)], b8.get((0, 3), K.zero)
    H8, H0, Hi = [], [], []
    for j in range(150):
        h0 = bplane[j].get((0, 2), K.zero) / a8
        h1 = (bplane[j].get((0, 3), K.zero) - h0 * a83) / a8
        H8.append({word: q for word, q in {(0, 0): h0, (0, 1): h1}.items() if q})
        numerator = add(bplane[j], scale(multiply(H8[j], b8), -K.one))
        assert all(i > 0 for i, k in numerator)
        H0.append({(i - 1, k): q for (i, k), q in numerator.items()})
        Hi.append([{word: q for word, q in {(0, 0): beta[j][i], (1, 0): mixed[i - 1][0][j] - h0 * mixed[i - 1][0][j8], (0, 1): mixed[i - 1][1][j] - h0 * mixed[i - 1][1][j8]}.items() if q} for i in range(1, 8)])
    assert H0 == list(map(poly, D['division_h0'])) and H8 == list(map(poly, D['division_h8']))
    assert Hi == [[poly(p) for p in row] for row in D['division_h_linear']]
    assert cp[j0] == one and all(not p.get((0, 0), K.zero) for j, p in enumerate(cp) if j != j0)
    c8 = add(cp[j8], *(scale(cp[j], -beta[j8][i]) for i, j in enumerate(minimum)))
    b8ca = add(bc[j8], *(scale(bc[j], -beta[j8][i]) for i, j in enumerate(minimum)))
    assert not c8.get((0, 0), K.zero)
    defects = []
    for j in range(150):
        h0 = compose(H0[j], bc[j0], 2)
        his = [compose(p, bc[j0], 1) for p in Hi[j]]
        predicted = add(multiply(h0, bc[j0]), multiply(H8[j], b8ca), *(multiply(p, bc[k]) for p, k in zip(his, minimum[1:])))
        assert not add(bc[j], scale(predicted, -K.one)), j
        defect = add(cp[j], scale(h0, -K.one), scale(multiply(H8[j], c8, 2), -K.one), *(scale(multiply(p, cp[k], 2), -K.one) for p, k in zip(his, minimum[1:])))
        assert all(sum(word) == 2 for word in defect)
        defects.append(defect)
    assert defects == list(map(poly, D['all150_syzygy_defect_quadratics']))
    rows = D['actual_Z4_quadratic_rows']
    assert rows == [18, 63]
    Z = [[zp[j].get(word, K.zero) for word in ((2, 0), (1, 1), (0, 2))] for j in rows]
    jj, kk = next((i, j) for i, j in it.combinations(range(3), 2) if Z[0][i] * Z[1][j] - Z[0][j] * Z[1][i])
    det = Z[0][jj] * Z[1][kk] - Z[0][kk] * Z[1][jj]
    for p in defects:
        v = [p.get(word, K.zero) for word in ((2, 0), (1, 1), (0, 2))]
        c0, c1 = (v[jj] * Z[1][kk] - v[kk] * Z[1][jj]) / det, (Z[0][jj] * v[kk] - Z[0][kk] * v[jj]) / det
        assert v == [c0 * Z[0][i] + c1 * Z[1][i] for i in range(3)]
    assert all(not bc[j] and not cp[j] for j in R)
    print('all150 b/c/H cofactors/defects and121-coordinate coverage PASS', time.monotonic() - start, flush=True)
    serial = lambda obj: str(obj) if not isinstance(obj, (dict, list)) else {k: serial(v) for k, v in obj.items()} if isinstance(obj, dict) else [serial(v) for v in obj]
    with (out / 'independent-source-curves.json.gz').open('wb') as f:
        with gzip.GzipFile(fileobj=f, mode='wb', filename='', mtime=0) as g:
            g.write((json.dumps(serial(archived), separators=(',', ':')) + '\n').encode())
    evidence = dict(own_index_sparse_Lie_operator_not_author_helper=True, own_all150_D3_mu_times_f_fibre_convolution=True,
        both_actual58_and87_constraint_matrices_and_inverses_verified=True, essential_source_curves_reconstructed=8,
        mixed_Hessian_projections_direct=14, old82_callbacks_or_worker_not_repeated=True,
        all150_conductor_b_and_CA_b_c_coefficients_exactly_match=True, all150_full_b_syzygies_through3=True,
        all121_nonbase_odd_coordinates_included=True, c_quadratic_generator_constant_zero=True,
        linear_h_degree2_quadratic_h_degree1_sufficiency_accepted=True,
        all150_defect_quadratics_in_actual_Z4_ideal_over_QQb=True, cubic_socle_all150_identically_zero=True,
        b_symbolic_no_parameter_sampling=True, full_D_application_not_claimed_in_this_gate=True)
    (out / 'independent-evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
    paths = [Path(__file__), ROOT / f'research/notes/{STEM}.md', ROOT / 'code/ambient_CB1_CA_relative_CI4_review_20261103_v63.py',
        folder / 'report.json', folder / 'generic-cubic.json.gz', folder / 'generic-setup.json.gz', rawpath, globalpath, capath, *deps]
    report = dict(id=STEM, all_passed=True, producer_or_original_worker_rerun=False,
        independent_original_three_source_all150_QQb_D3_binding=True, all58_87_constraints_and_full121_coverage_verified=True,
        all150_regseq_syzygy_degree_bounds_and_cubic_socle_zero_accepted=True,
        full_D_theorem_requires_separate_application_review=True,
        independent_method='Own coefficient-index Lie150; explicit noncommutative D3; cached inverse linear correction; eight exact slope curves plus14 bilinear Hessians; all150 regseq identities.',
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={name: sha(out / name) for name in ('independent-evidence.json', 'independent-source-curves.json.gz')}, seconds=time.monotonic() - start)
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    signal.alarm(0)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output)

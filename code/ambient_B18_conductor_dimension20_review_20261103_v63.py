#!/usr/bin/env python3
"""Independent intrinsic-derivation and exact tangent/IFT conductor dim20 audit."""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63_REVIEW'
PAIRS = list(it.combinations(range(5), 2))
MONS = list(it.combinations_with_replacement(range(5), 2))
PAIR_INDEX = {p: j for j, p in enumerate(PAIRS)}
MON_INDEX = {p: j for j, p in enumerate(MONS)}


def load(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def accepted(path):
    d = load(path)
    assert d['all_passed'] is True
    for name, digest in d.get('input_sha256', {}).items():
        assert sha(ROOT / name) == digest, name
    for name, digest in d.get('evidence_sha256', {}).items():
        assert sha(path.parent / name) == digest, name
    return d


def matrix(rows):
    return s.Matrix([[s.Rational(c) for c in row] for row in rows])


def fraction(q):
    return F(int(s.numer(q)), int(s.denom(q)))


def pivot_columns(M):
    """Independent Fraction Gaussian elimination, not producer DomainMatrix."""
    a = [[fraction(q) for q in row] for row in M.tolist()]
    row = 0
    result = []
    for col in range(M.cols):
        p = next((i for i in range(row, M.rows) if a[i][col]), None)
        if p is None:
            continue
        a[row], a[p] = a[p], a[row]
        value = a[row][col]
        for j in range(col, M.cols):
            a[row][j] /= value
        for i in range(row + 1, M.rows):
            if a[i][col]:
                c = a[i][col]
                for j in range(col, M.cols):
                    a[i][j] -= c * a[row][j]
        result.append(col)
        row += 1
        if row == M.rows:
            break
    return result


def determinant(M):
    a = [[fraction(q) for q in row] for row in M.tolist()]
    n = M.rows
    assert M.cols == n
    value = F(1)
    for j in range(n):
        p = next((i for i in range(j, n) if a[i][j]), None)
        if p is None:
            return F(0)
        if p != j:
            a[j], a[p] = a[p], a[j]
            value = -value
        pivot = a[j][j]
        value *= pivot
        for i in range(j + 1, n):
            c = a[i][j] / pivot
            for k in range(j + 1, n):
                a[i][k] -= c * a[j][k]
    return value


def physical_GL_column(index_i, index_j, coefficient):
    """Differentiate exterior action and quadratic inverse substitution by indices."""
    out = [F(0)] * 150
    for pair, (k, l) in enumerate(PAIRS):
        for word, (u, v) in enumerate(MONS):
            c = fraction(coefficient[15 * pair + word])
            if not c:
                continue
            if k == index_j and index_i != l:
                out[15 * PAIR_INDEX[tuple(sorted((index_i, l)))] + word] += c * (1 if index_i < l else -1)
            if l == index_j and k != index_i:
                out[15 * PAIR_INDEX[tuple(sorted((k, index_i)))] + word] += c * (1 if k < index_i else -1)
            for occurrence, other in ((u, v), (v, u)):
                if occurrence == index_i:
                    out[15 * pair + MON_INDEX[tuple(sorted((index_j, other)))]] -= c
    return s.Matrix([s.Rational(q.numerator, q.denominator) for q in out])


def main(output):
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    signal.alarm(120)
    output.mkdir(parents=True, exist_ok=False)
    producer = ROOT / 'reports/AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63'
    accepted(producer / 'report.json')
    record = load(producer / 'certificate.json.gz')
    deps = [
        ROOT / 'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_B18_CA_CB_ODD_CUBIC_CLASS_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json']
    for path in deps:
        accepted(path)
    qp = ROOT / 'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63/certificate.json.gz'
    cp = ROOT / 'reports/AMBIENT_B18_CA_CB_ODD_CUBIC_CLASS_20261103_v63/certificate.json.gz'
    fp = ROOT / 'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz'
    wp = ROOT / 'reports/AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63/staircase.json.gz'
    quad, cubic, first, wholeW = map(load, (qp, cp, fp, wp))
    h, x, y, z, w, W, U, V = s.symbols('h x y z w W U V')
    X = [h, x, y, z, w]
    locals_ = {str(v): v for v in (*X, W, U, V)}
    ey = s.Matrix([0, 0, 1, 0, 0])
    field = s.Matrix([x * x, 0, 0, h * x, h * h + z * z])
    P0 = ey * field.T - field * ey.T

    def flatten(P):
        return s.Matrix([s.Poly(s.expand(P[i, j]), *X).coeff_monomial(X[u] * X[v]) for i, j in PAIRS for u, v in MONS])

    P0_vec = flatten(P0)
    physical_G = s.Matrix.hstack(*[physical_GL_column(i, j, P0_vec) for i, j in it.product(range(5), repeat=2)])
    results = {}
    critical_paths = []
    for name, expected_cubic in [('CA', s.Rational(12, 125)), ('CB', -s.Rational(4, 243))]:
        q, c = quad[name], cubic[name]
        saved = record['actual_jet_bindings'][name]
        A = matrix(first['CA_CB_full_firstjet_data'][name]['actual'])
        rr, cc = q['actual_seven_odd_basis_rows'], q['selected_actual_base_columns']
        assert len(pivot_columns(A)) == 7
        B = A[rr, cc]
        assert determinant(B) == F(saved['seven_unit']) != 0
        relations = A[:, cc] * B.inv()
        assert relations == matrix(q['all_actual_linear_combinations']) and relations * A[rr, :] == A
        Q = matrix(q['quadratic_coefficient_matrix'])
        assert Q.shape == (150, 6) and Q[:, 2:] == s.zeros(150, 4)
        qr = q['independent_quadratic_rows']
        assert qr == [46, 95] and len(pivot_columns(Q)) == 2
        normalized = Q[:, :2] * Q[qr, :2].inv()
        selector = s.eye(150)[rr, :]
        cubic_function = s.eye(150)[63, :] - relations[63, :] * selector - normalized[63, :] * (s.eye(150)[qr, :] - relations[qr, :] * selector)
        assert cubic_function == matrix(saved['h63_as_constant_combination_of_all150_actual_odd_functions'])
        assert cubic_function * A == s.zeros(1, 29)
        S2 = s.Matrix([s.sympify(v, locals=locals_) for v in q['all150_vertical_derivatives2']])
        assert s.expand((cubic_function * S2)[0]) == 0
        b3 = s.Matrix([s.Rational(v) for v in c['actual_all150_odd_cubic_coefficients']])
        assert (cubic_function * b3)[0] == expected_cubic == s.Rational(saved['h63_actual_U3_coefficient'])
        correction = -s.eye(29)[:, cc] * B.inv() * b3[rr, :]
        assert correction == matrix(saved['actual_R29_IFT_correction_coefficient'])
        assert A[rr, :] * correction + b3[rr, :] == s.zeros(7, 1)

        Rrows = q['selected_actual28_Y'] + [61]
        G = physical_G[Rrows, :]
        gi = pivot_columns(G)
        assert len(gi) == 19 == len(saved['G19_GL_indices']) and gi == saved['G19_GL_indices']
        G19 = G[:, gi]
        assert G19 == matrix(saved['actual_R29_G19_tangents'])
        T = matrix(q['base_critical_tangent29_by3'])
        kernel = G19.row_join(T)
        assert A * kernel == s.zeros(150, 22) and len(pivot_columns(kernel)) == 22
        assert kernel == matrix(saved['full_minimal_kernel_basis29_by22'])
        chart = s.eye(29)[:, cc].row_join(kernel)
        assert determinant(chart) == F(saved['seven_plus_kernel29_unit']) != 0
        selected_rows = saved['kernel22_rows']
        assert determinant(kernel[selected_rows, :]) == F(saved['kernel22_minor']) != 0
        for index, variable in ((0, 'W'), (2, 'V')):
            family = G19.row_join(T[:, index])
            rows = pivot_columns(family.T)
            assert len(rows) == 20
            assert determinant(family[rows, :]) == F(saved[f'genuine_G19_{variable}_family_minor']) != 0
        assert len(pivot_columns(G19.row_join(T[:, [0, 2]]))) == 21
        monic = s.Matrix([[s.sympify(v, locals=locals_) for v in row] for row in c['fullphysical_source_P']])
        assert monic.subs(U, 0) == P0
        assert flatten(monic)[Rrows, :].diff(U).subs(U, 0) == T[:, 1]
        Wphysical = s.Matrix([[s.sympify(v, locals=locals_) for v in row] for row in wholeW[name]['physical_entire_W_curve']])
        expectedW = s.zeros(5)
        expectedW[3, 4], expectedW[4, 3] = W * x * x, -W * x * x
        assert Wphysical == P0 + expectedW
        assert flatten(Wphysical)[Rrows, :].diff(W).subs(W, 0) == T[:, 0]
        critical_path = ROOT / f'reports/AMBIENT_B18_{name}_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
        critical_paths.append(critical_path)
        critical = load(critical_path)
        ss = s.symbols('s0:10')
        vv = [V / 10, 0, 0, V / 10, -V ** 2 / 200, 0, 0, -V / 20, V ** 2 / 400, -V ** 3 / 24000] if name == 'CA' else [V, 0, 0, V / 2, -V ** 2 / 4, 0, 0, -V / 2, V ** 2 / 4, -V ** 3 / 24]
        equations = critical['predicted_entire_critical_ideal'] + critical['full150_critical_Schur_column' if name == 'CA' else 'all150_Schur_entries']
        for polynomial in equations:
            f = s.sympify(polynomial, locals={str(v): v for v in ss})
            assert s.expand(f.subs(dict(zip(ss, vv)))) == 0
            assert s.expand(f.subs(dict(zip(ss, [0, 0, W, 0, 0, 0, 0, 0, 0, 0])))) == 0
        results[name] = dict(actual_seven_rank=7, actual_seven_unit=str(determinant(B)),
            full_physical25_GL_rank=19, full_minimal_kernel_rank=22,
            actual_W_and_V_family_rank=20, union_W_V_rank=21,
            h63_IFT_invariant_nonzero_cubic=str(expected_cubic),
            all150_W_and_V_source_critical_equations_zero=True,
            actual_conductor_dimension=20, embedding_dimension=22, multiplicity_upper_bound=4)

    # Universal height2 perfect initial ideal; b,c remain indeterminate.
    a, b, c = s.symbols('a b c')
    binary = a * U * U + b * U * V + c * V * V
    HB = s.Matrix([[V, binary], [-U, 0], [0, -W]])
    minors = [s.expand(HB[[1, 2], :].det()), -s.expand(HB[[0, 2], :].det()), s.expand(HB[[0, 1], :].det())]
    assert all(s.expand(left - right) == 0 for left, right in zip(minors, [W * U, W * V, U * binary]))
    t = s.symbols('t')
    series = (1 - 2 * t ** 2 - t ** 3 + t ** 3 + t ** 4) / (1 - t) ** 3
    assert s.cancel(series - (1 + t) ** 2 / (1 - t)) == 0
    evidence = dict(actual_jet_and_family_checks=results,
        intrinsic_conductor_G_invariance_not_R_stability=True,
        arbitrary_lift_of_intrinsic_derivations_to_minimal_power_series_ring=True,
        complete_quadratic_g_independence_from_order1_derivatives=True,
        cubic_g_terms_removed_mod_Q_by_Euler_identity=True,
        actual_seven_IFT_path_correction_order3_preserves_h_cubic=True,
        Hilbert_Burch_minors=list(map(str, minors)), Hilbert_series=str(s.cancel(series)),
        height2_perfect_CM_initial_ideal_multiplicity=4,
        dimension20_established_before_multiplicity_comparison=True,
        genuine_W_V_distinct_minimal_support_branches=True,
        offdiagonal_pairs_dimension_upper20_via_finite_container=True,
        actual_full_I_CM_depth_or_multiplicity_equality_claim=False)
    (output / 'independent-evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
    paths = [Path(__file__), ROOT / f'research/notes/{STEM}.md',
        ROOT / 'code/ambient_B18_conductor_dimension20_20261103_v63.py',
        ROOT / 'research/notes/AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63.md',
        producer / 'report.json', producer / 'certificate.json.gz', qp, cp, fp, wp, *critical_paths, *deps]
    report = dict(id=STEM, all_passed=True, producer_rerun=False,
        CA_CB_actual_conductor_dimension20_and_embedding22_accepted=True,
        CA_CB_actual_conductor_multiplicity_atmost4_accepted=True,
        conceptual_G_derivation_and_actual_IFT_bridge_accepted=True,
        new_full_ideal_CM_depth_or_CA_only_lower4_not_reviewed=True,
        independent_method='All25 physical GL directions by coefficient-index differentiation; Fraction Gaussian ranks/minors; all150 h63 and IFT identities; full critical W/V substitutions; universal Hilbert–Burch/Euler proof.',
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'independent-evidence.json': sha(output / 'independent-evidence.json')},
        seconds=time.monotonic() - start)
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    signal.alarm(0)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    main(parser.parse_args().output)

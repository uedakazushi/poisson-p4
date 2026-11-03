#!/usr/bin/env python3
"""Independent old-R even/odd physical cotangent embedding and B18 tangent audit."""
import argparse
import gzip
import hashlib
import itertools as it
import json
from pathlib import Path
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix
from ambient_CA_origin_self_flag_4jet_review_20261103_v63 import Ring, direct_push, fraction
from ambient_B18_conductor_dimension20_review_20261103_v63 import physical_GL_column

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63_REVIEW'
PAIRS = list(it.combinations(range(5), 2))
MONS = list(it.combinations_with_replacement(range(5), 2))


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
        artifact = ROOT / name if (ROOT / name).exists() else path.parent / name
        assert sha(artifact) == digest, name
    return d


def rational_matrix(rows):
    return s.Matrix([[s.Rational(v) for v in row] for row in rows])


def main(out):
    start = time.monotonic()
    out.mkdir(parents=True, exist_ok=False)
    folder = ROOT / 'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63'
    accepted(folder / 'report.json')
    saved = load(folder / 'certificate.json.gz')
    gates = [ROOT / 'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_B18_CA_CB_ODD_CUBIC_CLASS_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_CA_CB_MINIMAL_LENGTHS_AND_R_CM_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63_REVIEW.json']
    records = [accepted(path) for path in gates]
    assert records[4]['actual_CA_CB_separate_full_formal_models_accepted']
    assert records[5]['Hilbert_Samuel_multiplicity'] == 2
    assert records[6]['generic_N_length1_by_common_reduced_Z_base_and_flat_functionfield_A10_intersection_A9']
    paths = [ROOT / 'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json',
        ROOT / 'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz',
        ROOT / 'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63/certificate.json.gz',
        ROOT / 'reports/AMBIENT_B18_CA_CB_ODD_CUBIC_CLASS_20261103_v63/certificate.json.gz']
    reduction, first, quad, cubic = map(load, paths)
    h, x, y, z, w, p = s.symbols('h x y z w p')
    X = s.Matrix([h, x, y, z, w])
    def flatten(P):
        return s.Matrix([s.Poly(s.expand(P[i, j]), *X).coeff_monomial(X[a] * X[b]) for i, j in PAIRS for a, b in MONS])
    def canon(P):
        d = s.Matrix([sum(s.diff(P[i, j], X[i]) for i in range(5)) for j in range(5)])
        return (P - (X * d.T - d * X.T) / 5).applyfunc(s.expand)
    def raw(f, g, H, cb):
        P = s.zeros(5)
        P[1, 2], P[3, 4] = f, g if cb else g * h * h
        P[2, 3] = (-s.diff(f, x) * z / 2 if cb else 0) - s.diff(H, w)
        P[2, 4] = (-s.diff(f, x) * w / 2 if cb else 0) + s.diff(H, z)
        return P - P.T
    swap = s.zeros(5)
    for i, j in enumerate((1, 0, 2, 3, 4)):
        swap[i, j] = 1
    change = lambda P: (swap * P.subs({h: x, x: h}, simultaneous=True) * swap.T).applyfunc(s.expand)
    cubics = [s.prod((h, x, z, w)[i] for i in indices) for indices in it.combinations_with_replacement(range(4), 3) if any(i >= 2 for i in indices)]
    ey = s.eye(5)[:, 2]
    field = s.Matrix([x * x, 0, 0, h * x, h * h + z * z])
    P0 = ey * field.T - field * ey.T
    point = flatten(P0)
    # Independent index-level derivative; no producer matrix Lie function.
    O = s.Matrix.hstack(*[physical_GL_column(i, j, point) for i, j in it.product(range(5), repeat=2)])
    Rrows = reduction['common28_linear_target_rows'] + [reduction['common29th_target_row']]
    assert Rrows == saved['common29_even_target_rows']
    planes, results = [], {}
    for name, expected_minor, expected_cubic in [('CA', s.Rational(9216, 25), s.Rational(12, 125)), ('CB', -s.Rational(576, 25), -s.Rational(4, 243))]:
        cb = name == 'CB'
        data = reduction['families'][name]
        record = saved['families'][name]
        basis = [change(canon(raw(f, 0, s.Integer(0), cb))) for f in (h * h, h * x, x * x)]
        basis += [change(canon(raw(0, h * h if cb else 1, s.Integer(0), cb)))]
        basis += [change(canon(raw(0, 0, H, cb))) for H in cubics]
        C = s.Matrix.hstack(*map(flatten, basis))
        J = O[:, reduction['orbit_GL_generator_indices']].row_join(C[:, data['source_fiber_basis_indices']])
        assert J == rational_matrix(data['linear_J']) == rational_matrix(record['full150_by28_physical_even_linear_J'])
        M = rational_matrix(data['source_kernel_Lie_matrix'])
        assert M ** 3 == s.zeros(5)
        center = s.zeros(20, 1)
        center[0], center[7], center[11], center[16] = -1, -1, 1, s.Rational(1, 3)
        alpha = center - p * s.Matrix([s.Rational(v) for v in data['source_velocity']])
        f = sum(alpha[i] * q for i, q in enumerate((h * h, h * x, x * x)))
        H = sum(q * m for q, m in zip(alpha[4:, 0], cubics))
        source = flatten(change(canon(raw(f, -alpha[3] * f if cb else alpha[3], H, cb))))
        source_series = [sum((Ring(fraction(s.Poly(q, p).nth(d))).shift(d) for d in range(3)), Ring()) for q in source]
        # Direct exterior action and inverse-coordinate substitution over our
        # own truncated coefficient ring, not producer's SymPy finite push.
        G = [[Ring(i == j) + Ring(fraction(M[i, j])).shift(1) + Ring(fraction((M * M)[i, j] / 2)).shift(2) for j in range(5)] for i in range(5)]
        pushed = direct_push(G, source_series)
        qcoefficient = lambda v, n: s.Rational(v.get((n, 0, 0, 0, 0), 0))
        assert s.Matrix([qcoefficient(v, 0) for v in pushed]) == point
        assert not any(qcoefficient(v, 1) for v in pushed)
        W = s.Matrix([qcoefficient(v, 2) for v in pushed])
        assert W == s.Matrix([s.Rational(v) for v in data['raw_second_coefficient']])
        assert W == s.Matrix([s.Rational(v) for v in record['physical_pure_kernel_second_W']])
        even = J.row_join(W)
        block = even[Rrows, :]
        assert DomainMatrix.from_Matrix(block).convert_to(s.QQ).det() == expected_minor
        E = even * block.inv()
        assert E[Rrows, :] == s.eye(29)
        assert E == rational_matrix(record['full150_by29_even_cotangent_embedding'])
        A = rational_matrix(first['CA_CB_full_firstjet_data'][name]['actual'])
        q = quad[name]
        rr, cc = q['actual_seven_odd_basis_rows'], q['selected_actual_base_columns']
        L = A[:, cc] * A[rr, cc].inv()
        assert L * A[rr, :] == A
        Qraw = rational_matrix(q['quadratic_coefficient_matrix'])[:, :2]
        qr = q['independent_quadratic_rows']
        Q = Qraw * Qraw[qr, :].inv()
        values = s.Matrix([s.Rational(v) for v in cubic[name]['actual_all150_odd_cubic_coefficients']])
        remainder = values - L * values[rr, :]
        remainder -= Q * remainder[qr, :]
        assert [i for i, v in enumerate(remainder) if v] == [63] and remainder[63] == expected_cubic
        B = L.row_join(Q).row_join(remainder / expected_cubic)
        assert B == rational_matrix(record['full150_by10_actual_odd_cotangent_embedding'])
        assert B[Rrows, :] == s.zeros(29, 10)
        T = E.row_join(B)
        assert T == rational_matrix(record['full150_by39_actual_tangent'])
        odd = rr + qr + [63]
        assert odd == saved['common10_odd_target_rows'] and len(set(Rrows + odd)) == 39
        # Independent block normalization avoids producer's39x39 inversion.
        odd_inverse = B[odd, :].inv()
        assert B[odd, :].det() == s.Rational(record['selected39_minor']) == 1
        Bn = B * odd_inverse
        normalized = (E - Bn * E[odd, :]).row_join(Bn)
        assert normalized[Rrows + odd, :] == s.eye(39)
        assert normalized == rational_matrix(saved['common150_by39_normalized_tangent'])
        planes.append(normalized)
        results[name] = dict(even_minor=str(expected_minor), pure_kernel_linear_zero=True,
            direct_physical_even_W_all150=True, normalized_all150_by39_entries_match=True,
            old_R_odd_basis_dimension=10)
    assert planes[0] == planes[1]
    assert point[106] == 1
    assert planes[0] * point[Rrows + odd, :] == point
    evidence = dict(families=results,
        true_kernel_square_even_shift_uses_R_order2_t_order1_I_in_m=True,
        old_R_mu10_from_intrinsic_edim39_and_exact_maximal_square=True,
        actual_I_mod_mI_basis_complete_not_arbitrary_Taylor_span=True,
        actual_CA_CB_same39plane_in_physical150=True, actual_D_tangent_dimension=39,
        actual_D_dim28_e2_from_top_N_only=True, depth_upper20_from_minimal_R20=True,
        projective_dim27_edim38_e2_depth_upper19=True,
        Dred_tangent_nil_linear_visibility_fullring_or_exact_depth_claim=False)
    (out / 'independent-evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
    inputs = [Path(__file__), ROOT / f'research/notes/{STEM}.md',
        ROOT / 'code/ambient_CA_origin_self_flag_4jet_review_20261103_v63.py', ROOT / 'code/ambient_B18_conductor_dimension20_review_20261103_v63.py',
        ROOT / 'code/ambient_B18_actual_intersection_tangent_20261103_v63.py', ROOT / 'research/notes/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63.md',
        folder / 'report.json', folder / 'certificate.json.gz', *paths, *gates]
    report = dict(id=STEM, all_passed=True, producer_rerun=False,
        independent_all150_by39_even_odd_embeddings_and_normalized_equality=True,
        genuine_even_kernel_square_and_old_R_mu10_application_accepted=True,
        actual_D_cone_dim28_edim39_e2_depth_atmost20_accepted=True,
        Dred_nilvisibility_fullring_exact_depth_not_inferred=True,
        independent_method='Own index GL columns; truncated-ring literal physical pushforward of pure kernel curve; block39 normalization; accepted actual I/mI jet classes.',
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in inputs},
        evidence_sha256={'independent-evidence.json': sha(out / 'independent-evidence.json')}, seconds=time.monotonic() - start)
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output)

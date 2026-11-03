#!/usr/bin/env python3
"""A literal R curve through the B18 origin and its actual N tangent."""
import argparse
import gzip
import hashlib
import itertools
import json
import signal
import time
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_R_CURVE_TANGENT_20261103_v63'


def read(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def main(out):
    start = time.monotonic()
    signal.alarm(60)
    out.mkdir(parents=True, exist_ok=False)
    reduction_path = ROOT / 'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json'
    first_path = ROOT / 'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz'
    quad_path = ROOT / 'reports/AMBIENT_B18_N_ODD_QUADRATIC_JETS_20261103_v63/certificate.json.gz'
    full_path = ROOT / 'reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63/report.json'
    spectral_path = ROOT / 'reports/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63/report.json'
    reviews = [
        ROOT / 'reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_B18_N_ODD_QUADRATIC_JETS_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_R19_BOUNDARY_LOCAL_MODELS_20261103_v63_REVIEW.json',
    ]
    assert all(read(p)['all_passed'] for p in reviews)
    reduction, first, quad, full, spectral = map(read, [reduction_path, first_path, quad_path, full_path, spectral_path])
    assert full['cone_dimension'] == 28 and full['cone_embedding_dimension'] == 37
    assert full['domain_and_analytically_irreducible']
    h, x, y, z, w, t = s.symbols('h x y z w t')
    X = s.Matrix([h, x, y, z, w])
    pairs = list(itertools.combinations(range(5), 2))
    mons = list(itertools.combinations_with_replacement(range(5), 2))
    ey = s.eye(5)[:, 2]
    # Exact source specialization r=0, b=d=1, e=0, c=25t²,
    # G=h²+z².  The source quadric r²+cd-25t² is identically zero.
    D = s.Matrix([x*x - 13*t*h*x + 25*t*t*x*z, 12*t*x*x, 0,
                  h*x - 13*t*x*z, h*h + z*z + 2*t*x*w])
    P = ey * D.T - D * ey.T
    v = s.Matrix([s.Poly(P[i, j], *X).coeff_monomial(X[u] * X[v])
                  for i, j in pairs for u, v in mons])
    p0 = v.subs(t, 0)
    dt = v.diff(t).subs(t, 0)
    assert [(i, q) for i, q in enumerate(dt) if q] == [(16, 13), (65, -12), (112, -13), (128, 2)]
    div = s.Matrix([sum(s.diff(P[i, j], X[i]) for i in range(5)) for j in range(5)])
    assert div == s.zeros(5, 1)
    strict = []
    for i, j, k in itertools.combinations(range(5), 3):
        q = s.expand(sum(P[i, l]*s.diff(P[j, k], X[l]) + P[j, l]*s.diff(P[k, i], X[l])
                         + P[k, l]*s.diff(P[i, j], X[l]) for l in range(5)))
        strict.append(q)
        assert q == 0
    # Bind the curve to the full stored spectral basis, rather than using a
    # sample of a normal form or merely a Poisson tangent condition.
    basis = s.Matrix.hstack(*(s.Matrix(c) for c in spectral['spectral_basis_raw150']))
    assert spectral['spectral_parameters'] == ['t', 'r', 'b', 'c', 'd', 'e', 'g0', 'g1', 'g2', 'g3', 'g4', 'g5']
    source = s.Matrix([t, 0, 1, 25*t*t, 1, 0, 1, 0, 0, 0, 0, 1])
    assert (basis * source - v).applyfunc(s.expand) == s.zeros(150, 1)
    assert s.expand(source[1]**2 + source[3]*source[4] - 25*source[0]**2) == 0
    Jca = s.Matrix(reduction['families']['CA']['linear_J'])
    O = Jca[:, :19]
    # The exact N restriction is k1=0, i.e. delete source column20.
    Jn = Jca[:, [i for i in range(28) if i != 20]]
    W = s.Matrix(first['N_full_firstjet_data']['W'])
    Y = reduction['common28_linear_target_rows']
    E0 = Jn.row_join(W)
    even_minor = s.factor(E0[Y, :].det())
    assert even_minor == s.Rational(2304, 5)
    E = E0 * E0[Y, :].inv()
    L = s.Matrix(quad['all_linear_row_combinations'])
    V, Wv = s.symbols('V W')
    Q0 = s.Matrix([[s.Poly(q, V, Wv).coeff_monomial(V*Wv),
                    s.Poly(q, V, Wv).coeff_monomial(V*V)]
                   for q in quad['all_zero_firstjet_generator_quadratics']])
    qr = quad['two_independent_zero_firstjet_rows']
    assert qr == [46, 72]
    Q = Q0 * Q0[qr, :].inv()
    odd7 = quad['seven_independent_actual_odd_rows']
    rows = Y + odd7 + qr
    T0 = E.row_join(L).row_join(Q)
    assert len(set(rows)) == 37 and T0[rows, :].det() == 1
    TN = T0 * T0[rows, :].inv()
    assert TN[rows, :] == s.eye(37)
    assert TN * dt[rows, :] == dt
    # A stronger first-order inclusion also holds: this R direction lies in
    # the rank27 source differential, but not in the rank19 group orbit.
    assert O.rank() == 19 and O.row_join(dt).rank() == 20
    assert Jn.rank() == 27 and Jn.row_join(dt).rank() == 27
    TR = O.row_join(dt)
    cols = list(TR.rref()[1])
    rr = list(TR[:, cols].T.rref()[1])
    minor20 = s.factor(TR[rr, cols].det())
    assert minor20 != 0
    assert TN.row_join(TR).rank() == 37
    # Full R source: the quadric is smooth (dq/dc=1); the accepted entire
    # boundary-fiber theorem gives one geometric flag at a=b=1.  Twenty
    # independent image tangents then make that whole fiber length1.
    evidence = dict(
        physical_curve150=list(map(str, v)), spectral_parameter_curve=list(map(str, source)),
        physical_R_curve_derivative150=list(map(str, dt)),
        N_even28_minor=str(even_minor), N37_coordinate_rows=rows,
        full150_by37_actual_N_normalized_tangent=[[str(q) for q in row] for row in TN.tolist()],
        R_curve_coordinates_in_N37=list(map(str, dt[rows, :])),
        orbit_plus_R_curve20_minor_rows=rr, orbit_plus_R_curve20_minor_columns=cols,
        orbit_plus_R_curve20_minor=str(minor20),
    )
    with (out / 'certificate.json.gz').open('wb') as f:
        with gzip.GzipFile(fileobj=f, mode='wb', filename='', mtime=0) as g:
            g.write((json.dumps(evidence, separators=(',', ':')) + '\n').encode())
    paths = [Path(__file__), ROOT / f'research/notes/{STEM}.md', reduction_path, first_path,
             quad_path, full_path, spectral_path, *reviews]
    result = dict(
        id=STEM, all_passed=True,
        exact_polynomial_curve_in_actual_proper_R_source=True,
        all25_divergence_and350_strict_Jacobi_coefficients_identically_zero=True,
        R_curve_physical_tangent_sparse={'16': '13', '65': '-12', '112': '-13', '128': '2'},
        physical_orbit_plus_curve_rank=20, R_cone_smooth_at_P0=True,
        R_G19_transverse_branch_is_regular_DVR=True,
        N_completed_branch_is_domain_by_accepted_full_completion=True,
        actual_N_tangent_dimension=37, actual_R_tangent_dimension=20,
        actual_R_tangent_contained_in_actual_N_tangent=True,
        actual_R_curve_tangent_even_contained_in_rank27_N_source_differential=True,
        actual_N_plus_R_tangent_span_dimension=37,
        N37_curve_coordinate_expression='13*E_Y16 - 12*E_Y65 - 13*E_Y112',
        exact_reduced_union_tangent_or_nilradical_visibility_from_this_inclusion_claimed=False,
        R_N_scheme_contact_length_or_full_origin_D_type_claimed=False,
        scope='Physical curve and full component tangents; sum of branch tangents need not equal tangent of their reduced union.',
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'certificate.json.gz': sha(out / 'certificate.json.gz')},
        seconds=time.monotonic() - start,
    )
    (out / 'report.json').write_text(json.dumps(result, indent=2) + '\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'physical_orbit_plus_curve_rank',
        'actual_N_tangent_dimension', 'actual_R_tangent_contained_in_actual_N_tangent', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output.resolve())

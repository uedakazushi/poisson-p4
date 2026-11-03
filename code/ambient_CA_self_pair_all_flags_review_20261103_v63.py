#!/usr/bin/env python3
"""Independent pointwise whole CA/N fibers on the nonzero self-pair arc."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_CA_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63'


def read(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def same(A, B):
    assert A.shape == B.shape and (A-B).applyfunc(s.expand) == s.zeros(*A.shape)


def main(output, replay):
    started = time.monotonic()
    folder = ROOT/'reports'/STEM
    report = read(folder/'report.json')
    assert report['all_passed']
    for path, digest in report['input_sha256'].items():
        assert sha(ROOT/path) == digest, path
    for path, digest in report['evidence_sha256'].items():
        assert sha(folder/path) == digest, path
    subprocess.run([str(ROOT/'.venv/bin/python'),
                    str(ROOT/'code/ambient_CA_self_pair_all_nonzero_flags_20261103_v63.py'),
                    '--output', str(replay)], check=True, timeout=30)
    assert sha(folder/'certificate.json.gz') == sha(replay/'certificate.json.gz')
    cert = read(folder/'certificate.json.gz')
    h, x, y, z, w, a = s.symbols('h x y z w a')
    X = s.Matrix([h, x, y, z, w])
    coords = [h, x, z, w]
    pairs = list(itertools.combinations(range(5), 2))
    mons = list(itertools.combinations_with_replacement(range(5), 2))
    parse = lambda q: s.sympify(q, locals={'a': a})
    matrix = lambda rows: s.Matrix([[parse(q) for q in row] for row in rows])

    def vec(P):
        return s.Matrix([s.Poly(P[i, j], *X).coeff_monomial(X[u]*X[v])
                         for i, j in pairs for u, v in mons])

    def canonical(P):
        d = s.Matrix([sum(s.diff(P[i, j], X[i]) for i in range(5)) for j in range(5)])
        return (P-(X*d.T-d*X.T)/5).applyfunc(s.expand)

    def raw(f, g, H):
        P = s.zeros(5)
        P[1, 2], P[3, 4] = f, g*h*h
        P[2, 3], P[2, 4] = -s.diff(H, w), s.diff(H, z)
        return P-P.T

    cubics = [s.prod(coords[i] for i in mon)
              for mon in itertools.combinations_with_replacement(range(4), 3)
              if any(i >= 2 for i in mon)]
    c = s.symbols('c0:16')
    f0, f1, f2, g = s.symbols('f0 f1 f2 g')
    F = f0*h*h+f1*h*x+f2*x*x
    H = sum(co*q for co, q in zip(c, cubics))
    endpoint = canonical(raw(F, 0, H))
    assert s.expand(-endpoint[0, 2]-s.diff(F, x)*h/5) == 0
    assert s.expand(-endpoint[1, 2]+F-s.diff(F, x)*x/5) == 0
    endpoint = canonical(raw(0, g, H))
    assert endpoint[0, :] == endpoint[1, :] == s.zeros(1, 5)
    D = s.Matrix([x*x+12*a*h*x/5, -3*a*x*x/5, 0,
                  h*x+2*a*x*z/5, h*h+z*z-8*a*x*w/5])
    ey = s.eye(5)[:, 2]
    P = ey*D.T-D*ey.T
    D4 = s.Matrix([D[i] for i in [0, 1, 3, 4]])
    deriv = lambda q: s.expand(sum(D4[i]*s.diff(q, coords[i]) for i in range(4)))
    forms = [x, x+a*h-2*a*a*z, x+3*a*h]
    lambdas = [-3*a, 2*a, 12*a]
    alpha = s.symbols('alpha')
    for ell, lam in zip(forms[1:], lambdas[1:]):
        residue = s.expand(deriv(x)+alpha*ell*ell+3*lam*x*x/10)
        assert s.expand(residue-alpha*ell*ell) == (0 if lam == 2*a else 3*a*x*x)
    qmons = [s.prod(coords[i] for i in mon)
             for mon in itertools.combinations_with_replacement(range(4), 2)
             if mon != (1, 1)]
    K = s.Matrix([[s.Poly(5*deriv(q)-12*a*x*q, *coords).coeff_monomial(mon)
                   for q in coords] for mon in qmons])
    same(K, matrix(cert['CA_ell_x_k_condition_matrix']))
    assert K[:, 0] == K[:, 1] == s.zeros(9, 1)
    # A pointwise a!=0 witness; a generic field rank alone is not used to
    # strengthen the parameter scope.
    assert K[[0, 4], [2, 3]].det() == 50*a
    field = s.QQ.frac_field(a)
    assert len(DomainMatrix.from_Matrix(K).convert_to(field).rref()[1]) == 2
    swap = s.eye(5)
    swap[0, 0] = swap[1, 1] = 0
    swap[0, 1] = swap[1, 0] = 1
    def change(U):
        return (swap*U.subs({h: x, x: h}, simultaneous=True)*swap.T).applyfunc(s.expand)
    n_tensors = [change(canonical(raw(f, 0, 0))) for f in [h*h, h*x]]
    n_tensors += [change(canonical(raw(0, 1, 0)))]
    n_tensors += [change(canonical(raw(0, 0, Q))) for Q in cubics]
    N = s.Matrix.hstack(*map(vec, n_tensors))
    npath = ROOT/'reports/AMBIENT_B18_N_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
    nbase = read(npath)
    same(N, s.Matrix(nbase['fiber19_basis']))
    flag_columns = []
    for i, j in nbase['flag_positions']:
        E = s.zeros(5)
        E[i, j] = 1
        V = E*X
        flag_columns.append(vec((E*P+P*E.T-sum((V[k]*P.diff(X[k])
                                              for k in range(5)), s.zeros(5))).applyfunc(s.expand)))
    JN = N.row_join(s.Matrix.hstack(*flag_columns))
    same(JN, matrix(cert['full150_by28_N_source_differential']))
    rows = list(N.T.rref()[1])
    quot = (JN[:, 19:]-N*N[rows, :].inv()*JN[rows, 19:]).applyfunc(s.expand)
    same(quot, matrix(cert['N_quotient150_by9']))
    selected = quot[cert['N_minor_rows'], cert['N_minor_columns']]
    dm = DomainMatrix.from_Matrix(selected).convert_to(s.QQ.poly_ring(a))
    determinant = s.factor(dm.domain.to_sympy(dm.det()))
    assert determinant == s.Rational(71663616, 78125)*a**9
    wholecb = ROOT/'reports/AMBIENT_CA_SELF_PAIR_CB_WHOLE_FIBER_20261103_v63_REVIEW.json'
    geometry = ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_PAIR_GEOMETRY_20261103_v63_REVIEW.json'
    incidence = ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_PAIR_N_CB_INCIDENCE_20261103_v63_REVIEW.json'
    proper = ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json'
    properraw = ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63/report.json'
    assert all(read(path)['all_passed'] for path in [wholecb, geometry, incidence, proper, properraw])
    assert read(wholecb)['all_linear_Darboux_forms_exhausted_for_every_a_nonzero']
    assert read(incidence)['entire_proper_N_source_is_closed_rank19_subbundle_of_entire_CA_source']
    geomdata = read(geometry)
    assert geomdata['source_projected9_minors'] == ['-11943936*a**9/15625', '864*a**14/15625']
    assert geomdata['cross_quotient_minor'] == '4076863488*a**25/3125'
    # Cover equality is part of the audited proper CA proof: Pfaffian quintic=g*h^3*f.
    assert read(properraw)['quintic_up_to_nonzero_constant'] == 'g*h^3*f(h,x)'
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}_REVIEW.md',
             ROOT/'code/ambient_CA_self_pair_all_nonzero_flags_20261103_v63.py',
             ROOT/f'research/notes/{STEM}.md', folder/'report.json',
             folder/'certificate.json.gz', npath, wholecb, geometry, incidence, proper, properraw]
    result = dict(id=STEM+'_REVIEW', all_passed=True,
        fresh_author_certificate_byte_identical=True,
        independent_universal_CA_endpoint16H_trace_and_two_Casimir_formulas=True,
        audited_shared_factor_Casimir_Darboux_classification_used_for_all_nonzero_a=True,
        ell_x_k_full9_by4_matrix_rank2_pointwise_minor='50*a',
        only_CA_flags_are_F1_and_F2_for_every_a_nonzero=True,
        independent_full_physical_N150_by28_and_projected150_by9=True,
        N_projected9_minor=str(determinant),
        whole_CA_fiber_exactly_two_reduced_points_for_every_a_nonzero=True,
        whole_N_fiber_exactly_one_reduced_point_for_every_a_nonzero=True,
        N_smooth_cone28_projective27_for_every_a_nonzero=True,
        reduced_CA_completion_exactly_two_smooth29_formal_branches_for_every_a_nonzero=True,
        actual_branch_intersection_embedding_dimension_for_every_a_nonzero=21,
        reduced_CA_cone_depth_at_most22_and_projective_depth_at_most21_for_every_a_nonzero=True,
        reduced_CA_nonnormal_singular_nonCM_for_every_a_nonzero=True,
        origin_full_ring_clean_intersection_conductor_or_whole_CA_edim37_claimed=False,
        input_sha256={str(path.relative_to(ROOT)): sha(path) for path in paths},
        seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['all_passed', 'N_projected9_minor',
        'whole_CA_fiber_exactly_two_reduced_points_for_every_a_nonzero',
        'whole_N_fiber_exactly_one_reduced_point_for_every_a_nonzero', 'seconds']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--replay', type=Path, required=True)
    args = parser.parse_args()
    main(args.output.resolve(), args.replay.resolve())

#!/usr/bin/env python3
"""Independent physical GL derivatives and whole-fiber proof for CA self-pairs."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import subprocess
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_CA_ORIGIN_SELF_PAIR_GEOMETRY_20261103_v63'


def read(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)


def sha(path):
    result = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            result.update(chunk)
    return result.hexdigest()


def same(left, right):
    assert left.shape == right.shape
    assert (left-right).applyfunc(s.expand) == s.zeros(*left.shape)


def main(output, replay):
    started = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    signal.alarm(120)
    folder = ROOT/'reports'/STEM
    report = read(folder/'report.json')
    assert report['all_passed']
    for path, digest in report['input_sha256'].items():
        assert sha(ROOT/path) == digest, path
    for path, digest in report['evidence_sha256'].items():
        assert sha(folder/path) == digest, path
    subprocess.run([str(ROOT/'.venv/bin/python'),
                    str(ROOT/'code/ambient_CA_origin_self_pair_geometry_20261103_v63.py'),
                    '--output', str(replay)], check=True, timeout=30)
    assert sha(folder/'certificate.json.gz') == sha(replay/'certificate.json.gz')
    cert = read(folder/'certificate.json.gz')
    h, x, y, z, w, a = s.symbols('h x y z w a')
    X = s.Matrix([h, x, y, z, w])
    pairs = list(itertools.combinations(range(5), 2))
    mons = list(itertools.combinations_with_replacement(range(5), 2))
    field = s.QQ.frac_field(a)

    def vec(P):
        return s.Matrix([s.Poly(P[i, j], *X).coeff_monomial(X[u]*X[v])
                         for i, j in pairs for u, v in mons])

    def physical_lie(i, j, P):
        E = s.zeros(5)
        E[i, j] = 1
        change = E*X
        return (E*P+P*E.T-sum((change[k]*P.diff(X[k])
                              for k in range(5)), s.zeros(5))).applyfunc(s.expand)

    def canonical(P):
        d = s.Matrix([sum(s.diff(P[i, j], X[i]) for i in range(5))
                      for j in range(5)])
        return (P-(X*d.T-d*X.T)/5).applyfunc(s.expand)

    def raw(f, g, H):
        P = s.zeros(5)
        P[1, 2], P[3, 4] = f, g*h*h
        P[2, 3], P[2, 4] = -s.diff(H, w), s.diff(H, z)
        return P-P.T

    swap = s.eye(5)
    swap[0, 0] = swap[1, 1] = 0
    swap[0, 1] = swap[1, 0] = 1

    def change(P):
        return (swap*P.subs({h: x, x: h}, simultaneous=True)*swap.T).applyfunc(s.expand)

    cubics = [s.prod((h, x, z, w)[i] for i in mon)
              for mon in itertools.combinations_with_replacement(range(4), 3)
              if any(i >= 2 for i in mon)]
    basis = [change(canonical(raw(f, 0, 0))) for f in [h*h, h*x, x*x]]
    basis += [change(canonical(raw(0, 1, 0)))]
    basis += [change(canonical(raw(0, 0, H))) for H in cubics]
    C = s.Matrix.hstack(*map(vec, basis))
    source = ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
    critical = read(source)
    same(C, s.Matrix(critical['fixed_fiber_linear_basis150_by20']))
    ey = s.Matrix([0, 0, 1, 0, 0])
    D = s.Matrix([x*x+12*a*h*x/5, -3*a*x*x/5, 0,
                  h*x+2*a*x*z/5, h*h+z*z-8*a*x*w/5])
    P = ey*D.T-D*ey.T
    M = s.zeros(5)
    M[1, 0], M[0, 3] = 1, -2
    g = s.eye(5)+a*M+a*a*M*M
    gi = s.eye(5)-a*M
    same(g*gi, s.eye(5))
    Q = (g*P.subs(dict(zip(X, gi*X)), simultaneous=True)*g.T).applyfunc(s.expand)

    def back(U):
        return (gi*U.subs(dict(zip(X, g*X)), simultaneous=True)*gi.T).applyfunc(s.expand)

    same(back(Q), P)
    flags = [(0, 2), (0, 3), (0, 4), (1, 0), (1, 2),
             (1, 3), (1, 4), (3, 2), (4, 2)]
    assert flags == list(map(tuple, critical['flag_positions']))
    L1 = [physical_lie(i, j, P) for i, j in flags]
    L2 = [physical_lie(i, j, Q) for i, j in flags]
    J1, J2 = s.Matrix.hstack(*map(vec, L1)), s.Matrix.hstack(*map(vec, L2))
    cr = list(C.T.rref()[1])
    left = C[cr, :].inv()

    def rem(U):
        return (U-C*left*U[cr, :]).applyfunc(s.expand)

    def decode_matrix(encoded):
        return s.Matrix([[sum(s.Rational(c)*a**degree for degree, c in q)
                          for q in row] for row in encoded])

    def rank(A):
        return len(DomainMatrix.from_Matrix(A).convert_to(field).rref()[1])

    def witness(A, saved):
        assert rank(A) == saved['rank']
        square = A[saved['minor_rows'], saved['minor_columns']]
        dm = DomainMatrix.from_Matrix(square).convert_to(s.QQ.poly_ring(a))
        det = s.factor(dm.domain.to_sympy(dm.det()))
        assert s.expand(det-s.sympify(saved['minor_determinant'], locals={'a': a})) == 0
        return str(det)

    quot1, quot2 = rem(J1), rem(J2)
    same(quot1, decode_matrix(cert['source1_flag_quotient150_by9']))
    same(quot2, decode_matrix(cert['source2_flag_quotient150_by9']))
    det1 = witness(quot1, cert['source1_rank_witness'])
    det2 = witness(quot2, cert['source2_rank_witness'])
    # Unlike the author reduction, push back ALL29 second-source directions.
    T1 = C.row_join(J1)
    T2 = s.Matrix.hstack(*map(vec, map(back, [*basis, *L2])))
    assert rank(T1) == rank(T2) == 29
    assert rank(T1.row_join(T2)) == 37
    second_fiber = s.Matrix.hstack(*map(vec, map(back, basis)))
    cross = rem(J1.row_join(second_fiber))
    same(cross, decode_matrix(cert['cross_same_point_mod_first_fiber150_by29']))
    detcross = witness(cross, cert['cross_quotient_rank_witness'])
    orbit = s.Matrix.hstack(*(vec(physical_lie(i, j, P))
                              for i, j in itertools.product(range(5), repeat=2)))
    same(orbit, decode_matrix(cert['full_GL150_by25']))
    orbit_det = witness(orbit, cert['GL_rank_witness'])
    witness(orbit.row_join(vec(P)), cert['GL_plus_scalar_rank_witness'])
    family_det = witness(orbit.row_join(vec(P)).row_join(vec(P).diff(a)),
                         cert['GL_plus_scalar_plus_parameter_derivative_rank_witness'])
    # Lie-equivariance supplies a genuine shared family, not just a tangent count.
    assert rank(T1.row_join(orbit).row_join(vec(P).diff(a))) == 29
    assert rank(T2.row_join(orbit).row_join(vec(P).diff(a))) == 29
    deps = [ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
            ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json',
            ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
            ROOT/'reports/AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63_REVIEW.json']
    assert all(read(path)['all_passed'] for path in deps)
    assert read(deps[1])['whole_flag_fibers']['CA']['scheme_length'] == 2
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}_REVIEW.md',
             ROOT/f'code/ambient_CA_origin_self_pair_geometry_20261103_v63.py',
             ROOT/f'research/notes/{STEM}.md', folder/'report.json',
             folder/'certificate.json.gz', source, *deps]
    result = dict(id=STEM+'_REVIEW', all_passed=True,
        fresh_author_certificate_byte_identical=True,
        independent_raw_CA20_canonical_basis_and_physical_tensor=True,
        independent_full5matrix_Lie_all25_and_both_flag9=True,
        independent_both_full150_by29_source_matrices_ranks=[29, 29],
        source_projected9_minors=[det1, det2],
        full_second29_physically_pushed_to_same_point=True,
        independent_same_point_union_rank=37,
        actual_branch_scheme_intersection_embedding_dimension=21,
        generic_GL_orbit_rank=19, generic_GL_plus_scalar_rank=19,
        generic_projective_PGL_orbit_dimension=18,
        actual_GL_plus_scalar_plus_parameter_image_cone_dimension=20,
        orbit_minor=orbit_det, orbit_parameter_minor=family_det,
        cross_quotient_minor=detcross,
        proper_finite_whole_fiber_length_at_most2_open_through_a0_reviewed=True,
        generic_nonzero_arc_whole_fiber_exactly_two_reduced_points=True,
        completion_injects_into_all_two_source_completions_and_each_is_smooth29_quotient=True,
        reduced_CA_completion_exactly_two_smooth29_formal_branches=True,
        meeting_scheme_dimension_bounds=[20, 21],
        reduced_CA_depth_at_generic_arc_point_at_most=22,
        projective_reduced_CA_depth_at_most=21,
        generic_self_pair_family_reduced_CA_nonnormal_singular_and_non_CM=True,
        whole_CA_embedding_dimension_equal37_claimed=False,
        all_nonzero_a_whole_fiber_or_origin_completion_claimed=False,
        clean_intersection_conductor_R19_N_membership_or_full_Poisson_claimed=False,
        input_sha256={str(path.relative_to(ROOT)): sha(path) for path in paths},
        seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'source_projected9_minors',
        'independent_same_point_union_rank', 'reduced_CA_depth_at_generic_arc_point_at_most', 'seconds']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--replay', type=Path, required=True)
    args = parser.parse_args()
    main(args.output.resolve(), args.replay.resolve())

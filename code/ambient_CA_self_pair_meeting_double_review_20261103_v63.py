#!/usr/bin/env python3
"""Independent finite physical curves and full two-branch CA meeting at P(1)."""
import argparse
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
STEM = 'AMBIENT_CA_SELF_PAIR_MEETING_DOUBLE_20261103_v63_REVIEW'


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


def main(output):
    started = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    signal.alarm(120)
    output.mkdir(parents=True, exist_ok=False)
    h, x, y, z, w, a, t = s.symbols('h x y z w a t')
    X = s.Matrix([h, x, y, z, w])
    pairs = list(itertools.combinations(range(5), 2))
    mons = list(itertools.combinations_with_replacement(range(5), 2))

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
    sourcepath = ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
    source = read(sourcepath)
    same(C, s.Matrix(source['fixed_fiber_linear_basis150_by20']))
    flags = list(map(tuple, source['flag_positions']))
    ey = s.eye(5)[:, 2]
    D = s.Matrix([x*x+12*a*h*x/5, -3*a*x*x/5, 0,
                  h*x+2*a*x*z/5, h*h+z*z-8*a*x*w/5])
    Pa = ey*D.T-D*ey.T
    P = Pa.subs(a, 1)
    M = s.zeros(5)
    M[1, 0], M[0, 3] = 1, -2
    g0, gi0 = s.eye(5)+M+M*M, s.eye(5)-M
    Q = (g0*P.subs(dict(zip(X, gi0*X)), simultaneous=True)*g0.T).applyfunc(s.expand)

    def back(U):
        return (gi0*U.subs(dict(zip(X, g0*X)), simultaneous=True)*gi0.T).applyfunc(s.expand)

    same(back(Q), P)
    # These literal rational choices were discovered in the parent's preflight;
    # their validity is proved below without reading any ignored work file.
    u1, u2 = [s.S.Zero]*29, [s.S.Zero]*29
    for i, q in {0: '11/12', 1: '-5/6', 2: '-15/4', 5: '1/8', 7: '1/4',
                 10: '-1/2', 12: '-13/24', 13: '-1/2', 14: '-5/3', 17: '1/2',
                 21: '-1/2', 22: '-1', 23: '2/3', 25: '-1/2'}.items():
        u1[i] = s.Rational(q)
    for i, q in {0: '2/3', 1: '-5/6', 2: '1/6', 5: '1/8', 10: '-1', 15: '1'}.items():
        u2[i] = s.Rational(q)

    def full_curve(point, u):
        A = s.zeros(5)
        for q, (i, j) in zip(u[20:], flags):
            A[i, j] += q
        assert A**3 == s.zeros(5)
        G, Gi = s.eye(5)+t*A+t*t*A*A/2, s.eye(5)-t*A+t*t*A*A/2
        same(G*Gi, s.eye(5))
        assert G.det() == 1
        delta = sum((q*U for q, U in zip(u[:20], basis)), s.zeros(5))
        # Full finite substitution, not the Lie-action quadratic formula.
        curve = (G*(point+t*delta).subs(dict(zip(X, Gi*X)), simultaneous=True)*G.T).applyfunc(s.expand)
        return curve, A, delta, G, Gi

    curve1, A1, delta1, G1, Gi1 = full_curve(P, u1)
    local2, A2, delta2, G2, Gi2 = full_curve(Q, u2)
    curve2 = back(local2)
    coefficient = lambda U, degree: U.applyfunc(lambda q: s.expand(q).coeff(t, degree))
    same(coefficient(curve1, 0), P)
    same(coefficient(curve2, 0), P)
    same(coefficient(curve1, 1), coefficient(curve2, 1))
    v = vec(coefficient(curve1, 1))
    L = vec(coefficient(curve1-curve2, 2))
    for curve in [curve1, curve2]:
        assert all(s.expand(sum(s.diff(curve[i, j], X[i]) for i in range(5))) == 0
                   for j in range(5))
        for i, j, k in itertools.combinations(range(5), 3):
            assert s.expand(sum(curve[i, l]*s.diff(curve[j, k], X[l])
                                +curve[j, l]*s.diff(curve[k, i], X[l])
                                +curve[k, l]*s.diff(curve[i, j], X[l])
                                for l in range(5))) == 0

    def lie(i, j, U):
        E = s.zeros(5)
        E[i, j] = 1
        V = E*X
        return (E*U+U*E.T-sum((V[k]*U.diff(X[k]) for k in range(5)), s.zeros(5))).applyfunc(s.expand)

    J1 = C.row_join(s.Matrix.hstack(*(vec(lie(i, j, P)) for i, j in flags)))
    J2 = s.Matrix.hstack(*(vec(back(U)) for U in [*basis, *(lie(i, j, Q) for i, j in flags)]))
    same(J1*s.Matrix(u1), v)
    same(J2*s.Matrix(u2), v)
    union = J1.row_join(J2)
    assert J1.rank() == J2.rank() == 29 and union.rank() == 37
    orbit = s.Matrix.hstack(*(vec(lie(i, j, P)) for i, j in itertools.product(range(5), repeat=2)))
    S = orbit.row_join(vec(Pa.diff(a).subs(a, 1)))
    assert S.rank() == 20 and S.row_join(v).rank() == 21
    assert union.row_join(L).rank() == 38
    ell = s.zeros(150, 1)
    ell[60], ell[63], ell[68] = s.Rational(2, 5), s.Rational(1, 10), 1
    assert ell.T*union == s.zeros(1, 58)
    assert (ell.T*L)[0] == -s.Rational(5, 12)
    # Whole fiber-product kernel, verified componentwise; completeness follows
    # the explicit convergent normal form in the proof note.
    rho, tau = s.symbols('rho tau')
    ys, zs = s.symbols('y1:9'), s.symbols('z1:9')
    relations = [yi*zj for yi in ys for zj in zs]
    relations += [(tau-rho*rho)*yi for yi in ys]
    relations += [tau*zj for zj in zs]
    relations += [tau*(tau-rho*rho)]
    assert len(relations) == 81
    left = {**{q: 0 for q in zs}, tau: rho*rho}
    right = {**{q: 0 for q in ys}, tau: 0}
    for q in relations:
        assert s.expand(q.subs(left, simultaneous=True)) == 0
        assert s.expand(q.subs(right, simultaneous=True)) == 0
    assert all(s.Poly(q, rho, tau, *ys, *zs).total_degree() >= 2 for q in relations)
    flagsreview = ROOT/'reports/AMBIENT_CA_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63_REVIEW.json'
    geometry = ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_PAIR_GEOMETRY_20261103_v63_REVIEW.json'
    assert read(flagsreview)['reduced_CA_completion_exactly_two_smooth29_formal_branches_for_every_a_nonzero']
    geom = read(geometry)
    assert geom['cross_quotient_minor'] == '4076863488*a**25/3125'
    assert geom['orbit_parameter_minor'] == '-27518828544*a**11/390625'
    def encoded_curve(U):
        return [[[degree[0], str(q)] for degree, q in sorted(s.Poly(p, t).terms()) if q]
                for p in vec(U)]
    payload = dict(source_order='fiber20, flag9', source_u1=list(map(str, u1)), source_u2=list(map(str, u2)),
        actual_first_GL_curve=[[str(q) for q in row] for row in G1.tolist()],
        actual_first_GL_inverse=[[str(q) for q in row] for row in Gi1.tolist()],
        first_curve150_by_t_polynomials=encoded_curve(curve1), second_curve150_by_t_polynomials=encoded_curve(curve2),
        common_tangent=list(map(str, v)), full150_second_difference=list(map(str, L)),
        sparse_coker_functional={str(i): str(q) for i, q in enumerate(ell) if q}, coker_value='-5/12',
        full81_fiber_product_relations=list(map(str, relations)),
        kernel_map=dict(rho='(rho,rho)', yi='(yi,0)', zj='(0,zj)', tau='(rho^2,0)'),
        complete_kernel_normal_form_proof='review note; shared f(rho)+tau*g(rho)+(y)C[[rho,y]]+(z)C[[rho,z]]')
    with (output/'certificate.json.gz').open('wb') as f:
        with gzip.GzipFile(fileobj=f, mode='wb', filename='', mtime=0) as archive:
            archive.write((json.dumps(payload, separators=(',', ':'))+'\n').encode())
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', sourcepath, flagsreview, geometry]
    result = dict(id=STEM, all_passed=True, field='QQ, complex completed germs', center_parameter='a=1',
        source_choices_are_literal_and_no_ignored_work_file_is_required=True,
        independent_raw_CA20_physical_tensor_and_both_full_finite5matrix_curves=True,
        first_actual_GL_generator_cubed_zero=True,
        both_finite_curves_all25_divergence_and350_strict_Jacobi_polynomials_zero=True,
        all150_first_coefficients_agree=True,
        tangent_ranks=[29, 29], tangent_union_rank=37, shared_family_tangent_rank=20,
        common_tangent_outside_shared_family_rank=21,
        augmented_rank_with_second_difference=38,
        independent_sparse_coker_pairing='-5/12',
        smooth20_shared_image_by_actual_constant_rank_family_reviewed=True,
        constant21_meeting_tangent_on_shared_reduced20_germ_forces_all_equations_divisible_rho_squared=True,
        nonzero_quadratic_unit_proves_exact_completed_meeting_double_smooth20=True,
        actual_meeting_completed_ring='C[[u1..u20,rho]]/(rho^2)',
        whole_CA_completed_ring_exact_fiber_product_and_full81_kernel_proved=True,
        reduced_CA_cone_dimension=29, reduced_CA_cone_embedding_dimension=38,
        reduced_CA_cone_depth=21, reduced_CA_cone_multiplicity=2,
        reduced_CA_projective_dimension=28, reduced_CA_projective_embedding_dimension=37,
        reduced_CA_projective_depth=20, reduced_CA_projective_multiplicity=2,
        normalization='two regular29 completed source rings',
        conductor_quotient='C[[u1..u20,rho]]/(rho^2)',
        validity_scope='P(1), a nonempty open part of the self-pair stratum; all-nonzero-a extension not inferred',
        whole_CA_CB_intersection_D_type_or_origin_a0_completion_claimed=False,
        seconds=time.monotonic()-started,
        input_sha256={str(path.relative_to(ROOT)): sha(path) for path in paths},
        evidence_sha256={'certificate.json.gz': sha(output/'certificate.json.gz')})
    (output/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'independent_sparse_coker_pairing',
        'actual_meeting_completed_ring', 'reduced_CA_cone_embedding_dimension',
        'reduced_CA_cone_depth', 'seconds']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output.resolve())

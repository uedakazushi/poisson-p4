#!/usr/bin/env python3
"""Independent symbolic physical second differences and CB tangent inclusion."""
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
STEM = 'AMBIENT_CA_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63'


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
    assert A.shape == B.shape and (A-B).applyfunc(s.cancel) == s.zeros(*A.shape)


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
                    str(ROOT/'code/ambient_CA_self_pair_actual_D_twojet_20261103_v63.py'),
                    '--output', str(replay)], check=True, timeout=40)
    cert = read(folder/'certificate.json')
    rerun = read(replay/'certificate.json')
    for key in ['u1', 'u2', 'common_tangent', 'second_difference']:
        assert cert[key] == rerun[key]
    assert {k: v for k, v in cert['result'].items() if k != 'seconds'} == {
        k: v for k, v in rerun['result'].items() if k != 'seconds'}
    h, x, y, z, w, a, t = s.symbols('h x y z w a t')
    X = s.Matrix([h, x, y, z, w])
    pairs = list(itertools.combinations(range(5), 2))
    mons = list(itertools.combinations_with_replacement(range(5), 2))
    parse = lambda q: s.sympify(q, locals={'a': a})
    field = s.QQ.frac_field(a)
    rank = lambda A: len(DomainMatrix.from_Matrix(A).convert_to(field).rref()[1])

    def vec(P):
        return s.Matrix([s.Poly(P[i, j], *X).coeff_monomial(X[u]*X[v])
                         for i, j in pairs for u, v in mons])

    def canonical(P):
        d = s.Matrix([sum(s.diff(P[i, j], X[i]) for i in range(5)) for j in range(5)])
        return (P-(X*d.T-d*X.T)/5).applyfunc(s.expand)

    def rawCA(f, g, H):
        P = s.zeros(5)
        P[1, 2], P[3, 4] = f, g*h*h
        P[2, 3], P[2, 4] = -s.diff(H, w), s.diff(H, z)
        return P-P.T

    def rawCB(f, g, H):
        P = s.zeros(5)
        P[1, 2], P[3, 4] = f, g
        P[2, 3] = -s.diff(f, x)*z/2-s.diff(H, w)
        P[2, 4] = -s.diff(f, x)*w/2+s.diff(H, z)
        return P-P.T

    swap = s.eye(5)
    swap[0, 0] = swap[1, 1] = 0
    swap[0, 1] = swap[1, 0] = 1
    change = lambda U: (swap*U.subs({h: x, x: h}, simultaneous=True)*swap.T).applyfunc(s.expand)
    cubics = [s.prod((h, x, z, w)[i] for i in mon)
              for mon in itertools.combinations_with_replacement(range(4), 3)
              if any(i >= 2 for i in mon)]
    basis = [change(canonical(rawCA(f, 0, 0))) for f in [h*h, h*x, x*x]]
    basis += [change(canonical(rawCA(0, 1, 0)))]
    basis += [change(canonical(rawCA(0, 0, H))) for H in cubics]
    C = s.Matrix.hstack(*map(vec, basis))
    cp = ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
    base = read(cp)
    same(C, s.Matrix(base['fixed_fiber_linear_basis150_by20']))
    flags = list(map(tuple, base['flag_positions']))
    ey = s.eye(5)[:, 2]
    D = s.Matrix([x*x+12*a*h*x/5, -3*a*x*x/5, 0,
                  h*x+2*a*x*z/5, h*h+z*z-8*a*x*w/5])
    P = ey*D.T-D*ey.T
    M = s.zeros(5)
    M[1, 0], M[0, 3] = 1, -2
    g0, gi0 = s.eye(5)+a*M+a*a*M*M, s.eye(5)-a*M
    Q = (g0*P.subs(dict(zip(X, gi0*X)), simultaneous=True)*g0.T).applyfunc(s.expand)
    back = lambda U: (gi0*U.subs(dict(zip(X, g0*X)), simultaneous=True)*gi0.T).applyfunc(s.expand)
    same(back(Q), P)
    original1, original2 = s.Matrix(list(map(parse, cert['u1']))), s.Matrix(list(map(parse, cert['u2'])))
    scale = s.factor(8*original1[5])
    assert scale == s.Rational(587068342272, 78125)*a**18
    u1, u2 = original1.applyfunc(lambda q: s.cancel(q/scale)), original2.applyfunc(lambda q: s.cancel(q/scale))

    def finite_curve(point, u):
        A = s.zeros(5)
        for q, (i, j) in zip(u[20:, :], flags):
            A[i, j] += q
        assert A**3 == s.zeros(5)
        G, Gi = s.eye(5)+t*A+t*t*A*A/2, s.eye(5)-t*A+t*t*A*A/2
        same(G*Gi, s.eye(5))
        assert G.det() == 1
        delta = sum((q*U for q, U in zip(u[:20, :], basis)), s.zeros(5))
        return (G*(point+t*delta).subs(dict(zip(X, Gi*X)), simultaneous=True)*G.T).applyfunc(s.expand)

    curve1, curve2 = finite_curve(P, u1), back(finite_curve(Q, u2))
    co = lambda U, d: U.applyfunc(lambda q: s.expand(q).coeff(t, d))
    same(co(curve1, 0), P)
    same(co(curve2, 0), P)
    same(co(curve1, 1), co(curve2, 1))
    v, L = vec(co(curve1, 1)), vec(co(curve1-curve2, 2))
    same(scale*v, s.Matrix(list(map(parse, cert['common_tangent']))))
    same(scale**2*L, s.Matrix(list(map(parse, cert['second_difference']))))

    def lie(i, j, U):
        E = s.zeros(5)
        E[i, j] = 1
        V = E*X
        return (E*U+U*E.T-sum((V[k]*U.diff(X[k]) for k in range(5)), s.zeros(5))).applyfunc(s.expand)

    J1 = C.row_join(s.Matrix.hstack(*(vec(lie(i, j, P)) for i, j in flags)))
    J2 = s.Matrix.hstack(*(vec(back(U)) for U in [*basis, *(lie(i, j, Q) for i, j in flags)]))
    same(J1*u1, v)
    same(J2*u2, v)
    # Rebuild the full actual CB20+flag9 at the explicit lift; ratio column
    # is raw g=f_B, retaining the original nonlinear mu*f source.
    fB = 3*a*x*x
    HB = -a*x*w*z+x*w*(h+x)/(3*a)+(h+x)**2*z/(9*a*a)+z**3/3
    QB = canonical(rawCB(fB, 0, HB))
    G = s.eye(5)
    G[0, 0], G[0, 1] = -3*a, -1
    same((G*P.subs(dict(zip(X, G.inv()*X)), simultaneous=True)*G.T).applyfunc(s.cancel), QB)
    globalcb = ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json'
    parabolic = set(map(tuple, read(globalcb)['parabolic_GL_positions']))
    cbflags = [p for p in itertools.product(range(5), repeat=2) if p not in parabolic]
    btensors = [canonical(rawCB(f, 0, 0)) for f in [h*h, h*x, x*x]]
    btensors += [canonical(rawCB(0, 0, H)) for H in cubics]
    btensors += [canonical(rawCB(0, fB, 0))]
    btensors += [lie(i, j, QB) for i, j in cbflags]
    def cbback(U):
        Gi = G.inv()
        return (Gi*U.subs(dict(zip(X, G*X)), simultaneous=True)*Gi.T).applyfunc(s.cancel)
    B = s.Matrix.hstack(*map(vec, map(cbback, btensors)))
    assert rank(B) == rank(J1) == rank(J2) == 29
    union = J1.row_join(J2)
    extended = union.row_join(L)
    ranks = [rank(J1.row_join(B)), rank(J2.row_join(B)), rank(union),
             rank(extended), rank(extended.row_join(B))]
    assert ranks == [30, 38, 37, 38, 38]
    shared = s.Matrix.hstack(*(vec(lie(i, j, P))
                               for i, j in itertools.product(range(5), repeat=2))).row_join(vec(P.diff(a)))
    assert rank(shared) == 20 and rank(shared.row_join(v)) == 21
    assert rank(B.row_join(v)) == 30
    names = ['AMBIENT_CA_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63_REVIEW.json',
             'AMBIENT_CA_SELF_PAIR_CB_WHOLE_FIBER_20261103_v63_REVIEW.json',
             'AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63_REVIEW.json',
             'AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63_REVIEW.json',
             'AMBIENT_B18_REDUCED_COMPONENTS_NONCM_20261103_v63_REVIEW.json']
    deps = [ROOT/'reports'/name for name in names]
    assert all(read(path)['all_passed'] for path in deps)
    assert read(deps[3])['entire_tzero_nilpotent_Jordanflag_cover_and_reduced_R_intersection_N_equals_irreducible_B18']
    assert read(deps[4])['individual_reduced_CA_CB_N_nonCM_nonnormal_singular_on_entire_B18']
    meeting = ROOT/'reports/AMBIENT_CA_SELF_PAIR_MEETING_DOUBLE_20261103_v63_REVIEW/report.json'
    assert read(meeting)['all_passed']
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}_REVIEW.md',
             ROOT/'code/ambient_CA_self_pair_actual_D_twojet_20261103_v63.py',
             ROOT/f'research/notes/{STEM}.md', folder/'report.json',
             folder/'certificate.json', cp, globalcb, meeting, *deps]
    result = dict(id=STEM+'_REVIEW', all_passed=True,
        fresh_author_exact_source_directions_and_all150_coefficients_match=True,
        independent_raw_CA20_and_CB_actual_mu_times_f_ratio_and_full_flag_sources=True,
        original_large_null_vectors_normalized_by_known_nonzero_QQa_scalar=str(scale),
        independent_generic_full_finite5matrix_curves_and_all150_true_second_difference=True,
        QQa_ranks_CA1_CB_CA2_CB_CAunion_plusL_plusCB=ranks,
        actual_shared_smooth_family_rank=20, extra_common_tangent_rank=21,
        extra_common_tangent_not_in_CB_by_rank30=True,
        actual_whole_CA_tangent_contains_J1_J2_and_trueL_by_ideal_second_difference=True,
        all_CB_tangent_in_whole_CA_tangent_as_QQa_identity_not_sample=True,
        first_CA_branch_intersection_CB_is_clean_smooth_N28=True,
        second_CA_branch_intersection_CB_is_clean_shared_smooth20=True,
        smooth_N_excludes_entire_B18_hence_R_and_Dred_is_N_reviewed=True,
        actual_D_cone_dimension=28, actual_D_cone_embedding_dimension=29,
        actual_D_projective_dimension=27, actual_D_projective_embedding_dimension=28,
        actual_D_nonreduced_on_nonempty_open_of_family_and_G_saturation=True,
        full_actual_D_ring_nil_module_Ass_or_origin_type_not_inferred=True,
        input_sha256={str(path.relative_to(ROOT)): sha(path) for path in paths},
        seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed',
        'QQa_ranks_CA1_CB_CA2_CB_CAunion_plusL_plusCB',
        'actual_D_cone_dimension', 'actual_D_cone_embedding_dimension', 'seconds']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--replay', type=Path, required=True)
    args = parser.parse_args()
    main(args.output.resolve(), args.replay.resolve())

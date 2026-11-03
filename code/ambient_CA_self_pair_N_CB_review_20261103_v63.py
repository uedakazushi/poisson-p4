#!/usr/bin/env python3
"""Physical N/CB incidence and the whole generic N fiber of a CA self-pair."""
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
STEM = 'AMBIENT_CA_ORIGIN_SELF_PAIR_N_CB_INCIDENCE_20261103_v63'


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
    for name, digest in report['input_sha256'].items():
        assert sha(ROOT/name) == digest, name
    for name, digest in report['evidence_sha256'].items():
        assert sha(folder/name) == digest, name
    subprocess.run([str(ROOT/'.venv/bin/python'),
                    str(ROOT/'code/ambient_CA_origin_self_pair_N_CB_incidence_20261103_v63.py'),
                    '--output', str(replay)], check=True, timeout=30)
    assert sha(folder/'certificate.json.gz') == sha(replay/'certificate.json.gz')
    data = read(folder/'certificate.json.gz')
    h, x, y, z, w, a = s.symbols('h x y z w a')
    X = s.Matrix([h, x, y, z, w])
    pairs = list(itertools.combinations(range(5), 2))
    mons = list(itertools.combinations_with_replacement(range(5), 2))

    def vec(P):
        return s.Matrix([s.Poly(P[i, j], *X).coeff_monomial(X[u]*X[v])
                         for i, j in pairs for u, v in mons])

    def canonical(P):
        div = s.Matrix([sum(s.diff(P[i, j], X[i]) for i in range(5))
                        for j in range(5)])
        return (P-(X*div.T-div*X.T)/5).applyfunc(s.expand)

    def raw(f, g, H):
        P = s.zeros(5)
        P[1, 2], P[3, 4] = f, g*h*h
        P[2, 3], P[2, 4] = -s.diff(H, w), s.diff(H, z)
        return P-P.T

    def lie(i, j, P):
        E = s.zeros(5)
        E[i, j] = 1
        V = E*X
        return (E*P+P*E.T-sum((V[k]*P.diff(X[k]) for k in range(5)),
                              s.zeros(5))).applyfunc(s.expand)

    swap = s.eye(5)
    swap[0, 0] = swap[1, 1] = 0
    swap[0, 1] = swap[1, 0] = 1

    def change(P):
        return (swap*P.subs({h: x, x: h}, simultaneous=True)*swap.T).applyfunc(s.expand)

    cubics = [s.prod((h, x, z, w)[i] for i in mon)
              for mon in itertools.combinations_with_replacement(range(4), 3)
              if any(i >= 2 for i in mon)]
    ca_tensors = [change(canonical(raw(f, 0, 0))) for f in [h*h, h*x, x*x]]
    ca_tensors += [change(canonical(raw(0, 1, 0)))]
    ca_tensors += [change(canonical(raw(0, 0, H))) for H in cubics]
    C = s.Matrix.hstack(*map(vec, ca_tensors))
    indices = [0, 1, *range(3, 20)]
    n_tensors = [ca_tensors[i] for i in indices]
    N = s.Matrix.hstack(*map(vec, n_tensors))
    np = ROOT/'reports/AMBIENT_B18_N_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
    cp = ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
    same(C, s.Matrix(read(cp)['fixed_fiber_linear_basis150_by20']))
    same(N, s.Matrix(read(np)['fiber19_basis']))
    flags = list(map(tuple, read(cp)['flag_positions']))
    nr = list(N.T.rref()[1])
    nleft = N[nr, :].inv()
    # The same flag parabolic preserves the19-dimensional subspace globally.
    # Together with its connectedness this verifies the associated closed subbundle.
    parabolic = [ij for ij in itertools.product(range(5), repeat=2) if ij not in flags]
    assert len(parabolic) == 16
    for i, j in parabolic:
        B = s.Matrix.hstack(*(vec(lie(i, j, P)) for P in n_tensors))
        same(B, N*nleft*B[nr, :])
    ey = s.Matrix([0, 0, 1, 0, 0])
    D = s.Matrix([x*x+12*a*h*x/5, -3*a*x*x/5, 0,
                  h*x+2*a*x*z/5, h*h+z*z-8*a*x*w/5])
    P = ey*D.T-D*ey.T
    ncoords = nleft*vec(P)[nr, :]
    same(N*ncoords, vec(P))
    same(ncoords, s.Matrix([s.sympify(q, locals={'a': a})
                            for q in data['first_source_N19_coordinates']]))
    M = s.zeros(5)
    M[1, 0], M[0, 3] = 1, -2
    g = s.eye(5)+a*M+a*a*M*M
    gi = s.eye(5)-a*M
    Q = (g*P.subs(dict(zip(X, gi*X)), simultaneous=True)*g.T).applyfunc(s.expand)
    cr = list(C.T.rref()[1])
    left = C[cr, :].inv()
    first = left*vec(P)[cr, :]
    second = left*vec(Q)[cr, :]
    same(C*second, vec(Q))
    assert first[2] == 0 and s.expand(second[2]+a*a) == 0
    JN = N.row_join(s.Matrix.hstack(*(vec(lie(i, j, P)) for i, j in flags)))
    field = s.QQ.frac_field(a)
    pc = list(DomainMatrix.from_Matrix(JN).convert_to(field).rref()[1])
    pr = list(DomainMatrix.from_Matrix(JN[:, pc].T).convert_to(field).rref()[1])
    assert len(pc) == len(pr) == 28
    dm = DomainMatrix.from_Matrix(JN[pr, pc]).convert_to(s.QQ.poly_ring(a))
    determinant = s.factor(dm.domain.to_sympy(dm.det()))
    assert determinant != 0
    # Actual CB formula, independently reconstructed in the physical coordinates.
    G = s.eye(5)
    G[0, 0], G[0, 1] = -3*a, -1
    assert G.det() == -3*a
    Gi = G.inv()
    changed = (G*P.subs(dict(zip(X, Gi*X)), simultaneous=True)*G.T).applyfunc(s.cancel)
    fB = 3*a*x*x
    H = -a*x*w*z+x*w*(h+x)/(3*a)+(h+x)**2*z/(9*a*a)+z**3/3
    RB = s.zeros(5)
    RB[1, 2] = fB
    RB[2, 3] = -s.diff(fB, x)*z/2-s.diff(H, w)
    RB[2, 4] = -s.diff(fB, x)*w/2+s.diff(H, z)
    CB = canonical(RB-RB.T)
    same(changed, CB)
    same(vec(changed), s.Matrix([s.sympify(q, locals={'a': a})
                                 for q in data['all150_changed_P']]))
    for tensor in [P, changed]:
        assert all(s.cancel(sum(s.diff(tensor[i, j], X[i]) for i in range(5))) == 0
                   for j in range(5))
        for i, j, k in itertools.combinations(range(5), 3):
            assert s.cancel(sum(tensor[i, l]*s.diff(tensor[j, k], X[l])
                                +tensor[j, l]*s.diff(tensor[k, i], X[l])
                                +tensor[k, l]*s.diff(tensor[i, j], X[l])
                                for l in range(5))) == 0
    geometry = ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_PAIR_GEOMETRY_20261103_v63_REVIEW.json'
    globaln = ROOT/'reports/AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63_REVIEW.json'
    assert read(geometry)['all_passed']
    assert read(globaln)['all_passed'] and read(globaln)['proper_full_N19_subspace_invariance_verified']
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}_REVIEW.md',
             ROOT/f'code/ambient_CA_origin_self_pair_N_CB_incidence_20261103_v63.py',
             ROOT/f'research/notes/{STEM}.md', folder/'report.json',
             folder/'certificate.json.gz', np, cp, geometry, globaln,
             ROOT/'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json']
    result = dict(id=STEM+'_REVIEW', all_passed=True,
        fresh_author_certificate_byte_identical=True,
        own_physical_raw_CA_and_N19_equals_CA20_deleted_column2=True,
        all16_connected_flag_parabolic_Lie_actions_preserve_N19=True,
        entire_proper_N_source_is_closed_rank19_subbundle_of_entire_CA_source=True,
        actual_N_source_membership_for_all_a=True,
        first_CA_flag_N_coefficient2='0', second_CA_flag_N_coefficient2='-a**2',
        explicit_CB_membership_on_Da_all150_finite_coordinate_identity=True,
        actual_CB_change_determinant='-3*a',
        original_and_changed_all25_divergence_and350_Jacobi_coefficients_zero=True,
        independent_N_source150_by28_rank=28,
        N_source_rank_minor=dict(rows=pr, columns=pc, determinant=str(determinant)),
        generic_nonzero_arc_whole_N_fiber_exactly_one_reduced_point=True,
        finite_N_source_equals_reduced_target_by_residue_length1_and_Nakayama=True,
        generic_nonzero_arc_reduced_N_is_smooth_cone28_projective27=True,
        all_nonzero_a_whole_fiber_or_CB_smoothness_claimed=False,
        R_membership_full_origin_CA_I_conductor_or_full_intersection_type_claimed=False,
        input_sha256={str(path.relative_to(ROOT)): sha(path) for path in paths},
        seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed',
        'independent_N_source150_by28_rank',
        'generic_nonzero_arc_reduced_N_is_smooth_cone28_projective27', 'seconds']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--replay', type=Path, required=True)
    args = parser.parse_args()
    main(args.output.resolve(), args.replay.resolve())

#!/usr/bin/env python3
"""Actual all150 second jets for R versus the N normalizer, including its kernel."""
import argparse
import gzip
import hashlib
import itertools
import json
import resource
import signal
import time
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_N_R_SECOND_CONTACT_20261103_v63'


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def main(out):
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    signal.alarm(60)
    out.mkdir(parents=True, exist_ok=False)
    reduction_path = ROOT / 'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json'
    first_path = ROOT / 'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz'
    curve_path = ROOT / 'reports/AMBIENT_B18_R_CURVE_TANGENT_20261103_v63/certificate.json.gz'
    reviews = [
        ROOT / 'reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
    ]
    assert all(read(p)['all_passed'] for p in reviews)
    reduction, first, curve = map(read, [reduction_path, first_path, curve_path])
    h, x, y, z, w, t, c = s.symbols('h x y z w t c')
    X = s.Matrix([h, x, y, z, w])
    pairs = list(itertools.combinations(range(5), 2))
    mons = list(itertools.combinations_with_replacement(range(5), 2))

    def vec(P):
        return s.Matrix([s.Poly(P[i, j], *X).coeff_monomial(X[a]*X[b])
                         for i, j in pairs for a, b in mons])

    def canonical(P):
        d = s.Matrix([sum(s.diff(P[i, j], X[i]) for i in range(5)) for j in range(5)])
        return (P - (X*d.T-d*X.T)/5).applyfunc(s.expand)

    def raw(f, g, H):
        P = s.zeros(5)
        P[1, 2] = f
        P[3, 4] = g*h*h
        P[2, 3] = -s.diff(H, w)
        P[2, 4] = s.diff(H, z)
        return P-P.T

    swap = s.eye(5)
    swap[0, 0] = swap[1, 1] = 0
    swap[0, 1] = swap[1, 0] = 1

    def change(P):
        return (swap*P.subs({h: x, x: h}, simultaneous=True)*swap.T).applyfunc(s.expand)

    cubics = [s.prod((h, x, z, w)[i] for i in m)
              for m in itertools.combinations_with_replacement(range(4), 3) if any(i >= 2 for i in m)]
    basis = [change(canonical(raw(f, 0, 0))) for f in [h*h, h*x, x*x]]
    basis += [change(canonical(raw(0, 1, 0)))]
    basis += [change(canonical(raw(0, 0, H))) for H in cubics]
    C = s.Matrix.hstack(*map(vec, basis))
    ey = s.eye(5)[:, 2]
    D = s.Matrix([x*x, 0, 0, h*x, h*h+z*z])
    P0 = ey*D.T-D*ey.T
    point = vec(P0)
    ca = reduction['families']['CA']
    J = s.Matrix(ca['linear_J'])[:, [i for i in range(28) if i != 20]]
    rr = first['N_own_linear27_rows']
    left = J[rr, :].inv()
    v = s.Matrix(curve['physical_R_curve_derivative150'])
    velocity = left*v[rr, :]
    assert J*velocity == v
    outerM = s.zeros(5)
    for q, ij in zip(velocity[:19], reduction['orbit_GL_generator_indices']):
        outerM[ij//5, ij % 5] += q
    extras = [i for i in ca['source_fiber_basis_indices'] if i != 2]
    assert len(extras) == 8 and J[:, 19:] == C[:, extras]
    fixed_velocity = s.zeros(20, 1)
    for q, j in zip(velocity[19:], extras):
        fixed_velocity[j] = q
    assert fixed_velocity[2] == 0
    kernelM = s.Matrix(ca['source_kernel_Lie_matrix'])
    co = s.Matrix(ca['source_velocity'])
    assert kernelM**3 == s.zeros(5) and co[2] == 0
    base = s.zeros(20, 1)
    base[0] = base[7] = -1
    base[11] = 1
    base[16] = s.Rational(1, 3)
    assert C*base == point
    alpha = base+t*fixed_velocity-c*t*co
    assert alpha[2] == 0  # Entire fixed source is in N's actual19-space.
    source = sum((q*P for q, P in zip(alpha, basis)), s.zeros(5))
    # Genuine normalizer arc exp(t*outerM) exp(c*t*kernelM) acting on source.
    # Only its coefficients through2 are computed; no third jet is executed.
    outer = s.eye(5)+t*outerM+t*t*outerM**2/2
    outer_inverse = s.eye(5)-t*outerM+t*t*outerM**2/2
    kernel = s.eye(5)+c*t*kernelM+c*c*t*t*kernelM**2/2
    kernel_inverse = s.eye(5)-c*t*kernelM+c*c*t*t*kernelM**2/2
    G, Gi = outer*kernel, kernel_inverse*outer_inverse
    truncated = lambda q: sum(s.Poly(s.expand(q), t).nth(i)*t**i for i in range(3))
    assert (G*Gi-s.eye(5)).applyfunc(truncated) == s.zeros(5)
    physical = (G*source.subs(dict(zip(X, Gi*X)), simultaneous=True)*G.T).applyfunc(truncated)
    q = vec(physical)
    N0 = q.applyfunc(lambda v: s.Poly(v, t).nth(0))
    N1 = q.applyfunc(lambda v: s.Poly(v, t).nth(1))
    N2 = q.applyfunc(lambda v: s.Poly(v, t).nth(2))
    assert N0 == point and (N1-v).applyfunc(s.expand) == s.zeros(150, 1)
    N20 = N2.subs(c, 0)
    mixed = N2.diff(c).subs(c, 0)
    W = N2.applyfunc(lambda q: s.Poly(q, c).nth(2))
    assert (N2-N20-c*mixed-c*c*W).applyfunc(s.expand) == s.zeros(150, 1)
    assert W == s.Matrix(first['N_full_firstjet_data']['W'])
    assert mixed == s.Matrix(first['N_full_firstjet_data']['mixed'])*velocity
    # Separate infinitesimal-action calculation is checked against the direct
    # inverse-coordinate physical pushforward, for every one of150 coefficients.
    def lie(M, P):
        MX = M*X
        return (M*P+P*M.T-sum((MX[k]*P.diff(X[k]) for k in range(5)), s.zeros(5))).applyfunc(s.expand)
    fixed = sum((a*P for a, P in zip(fixed_velocity, basis)), s.zeros(5))
    assert N20 == vec(lie(outerM, lie(outerM, P0))/2+lie(outerM, fixed))
    R2 = s.Matrix([s.Poly(a, t).nth(2) for a in curve['physical_curve150']])
    delta = R2-N20
    TN = s.Matrix(curve['full150_by37_actual_N_normalized_tangent'])
    nrows = curve['N37_coordinate_rows']
    assert delta-TN*delta[nrows, :] == s.zeros(150, 1)
    project = lambda V: (V-J*left*V[rr, :]).applyfunc(s.expand)
    d, m, ww = map(project, [delta, mixed, W])
    sparse = lambda V: {str(i): str(q) for i, q in enumerate(V) if q != 0}
    assert sparse(d) == {'61': '600', '124': '600'}
    assert sparse(ww) == {'27': '-3', '61': '-6', '124': '-6'}
    assert sparse(m) == {'27': '48', '61': '142', '124': '142'}
    mixed_square_minor = s.Matrix.hstack(ww, m)[[27, 61], :].det()
    assert mixed_square_minor == -138
    obstruction = d-c*m-c*c*ww
    p, q0 = 3*c*c-48*c, 6*c*c-142*c+600
    assert obstruction[27] == p and obstruction[61] == q0 and obstruction[124] == q0
    assert all(obstruction[i] == 0 for i in range(150) if i not in [27, 61, 124])
    bezout_u, bezout_v, gcd = s.gcdex(p, q0, c)
    assert gcd == 1 and s.expand(bezout_u*p+bezout_v*q0) == 1
    square_coefficient, mixed_coefficient = s.Rational(4800, 23), s.Rational(300, 23)
    assert d == square_coefficient*ww+mixed_coefficient*m
    discrepancy = mixed_coefficient**2-square_coefficient
    assert discrepancy == -s.Rational(20400, 529)
    # Source coefficients, full matrices and curvature cofactor certificates
    # are saved; the source-residual rows are not actual N conormals.
    evidence = dict(
        source_velocity27=list(map(str, velocity)), fixed_N20_velocity=list(map(str, fixed_velocity)),
        source_outer_GL_matrix=[[str(q) for q in row] for row in outerM.tolist()],
        source_kernel_GL_matrix=[[str(q) for q in row] for row in kernelM.tolist()],
        source_kernel_fixed20_compensation=list(map(str, co)),
        N2_kernel_zero150=list(map(str, N20)), N2_mixed150=list(map(str, mixed)),
        N2_kernel_square150=list(map(str, W)), R2_physical150=list(map(str, R2)),
        actual_curvature_difference150=list(map(str, delta)),
        source27_unit_rows=rr, source27_left_inverse=[[str(q) for q in row] for row in left.tolist()],
        source27_candidate_acceleration_kernel_zero=list(map(str, left*delta[rr, :])),
        source27_projected_delta=sparse(d), source27_projected_mixed=sparse(m),
        source27_projected_kernel_square=sparse(ww),
        full150_projected_second_lift_obstruction=list(map(str, obstruction)),
        second_lift_obstruction_polynomials=[str(p), str(q0)],
        exact_Bezout_identity_coefficients=[str(bezout_u), str(bezout_v)],
        delta_W_mixed_decomposition=[str(square_coefficient), str(mixed_coefficient)],
        scalar_square_compatibility_discrepancy=str(discrepancy),
        projected_kernel_square_mixed_rank2_minor=str(mixed_square_minor),
    )
    with (out/'certificate.json.gz').open('wb') as f:
        with gzip.GzipFile(fileobj=f, mode='wb', filename='', mtime=0) as g:
            g.write((json.dumps(evidence, separators=(',', ':'))+'\n').encode())
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', reduction_path, first_path,
             curve_path, curve_path.parent/'report.json', *reviews]
    result = dict(
        id=STEM, all_passed=True,
        independent_fixed_N19_basis_and_full150_physical_second_push=True,
        entire_source_keeps_deleted_binary_C2_zero=True,
        source_is_genuine_formal_N_arc_not_truncated_polynomial_claim=True,
        all150_mixed_kernel_and_kernel_square_coefficients_bound=True,
        R2_minus_N2_in_actual_TN37=True, R2_minus_N2_not_in_source_JN27=True,
        whole_N_ideal_restriction_to_R_curve_has_order_atleast3=True,
        source_kernel_second_lift_polynomials=[str(p), str(q0)],
        source_kernel_second_lift_ideal_is_unit_over_QQ=True,
        no_second_order_lift_in_the_actual_etale_N_normalizer_chart=True,
        discriminant_square_discrepancy=str(discrepancy),
        mixed_class_nonzero_mod_source_J_and_kernel_square=True,
        projected_kernel_square_mixed_rank2_minor=str(mixed_square_minor),
        exact_k_h_reduced_embedding_dimension_or_order4_relation_claimed=False,
        no_third_jet_executed=True,
        scope='Image curvature yields k>=3; source second-lift failure is not a defining-ideal obstruction of order2.',
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'certificate.json.gz': sha(out/'certificate.json.gz')},
        seconds=time.monotonic()-start,
    )
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'R2_minus_N2_in_actual_TN37',
          'R2_minus_N2_not_in_source_JN27', 'source_kernel_second_lift_ideal_is_unit_over_QQ', 'seconds']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    main(parser.parse_args().output.resolve())

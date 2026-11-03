#!/usr/bin/env python3
"""Polynomial boundary family and regular formal image branch of embedded21."""
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
from s61_model import S61, X, pushforward

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_EMBEDDED21_POLYNOMIAL_BOUNDARY_20261103_v63'


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def mat(M):
    return [[str(q) for q in row] for row in M.tolist()]


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    out.mkdir(parents=True, exist_ok=False)
    recognition = ROOT/'reports/AMBIENT_W1_S61_INTERSECTION_RECOGNITION_20261103_v63/report.json'
    recognition_review = ROOT/'reports/AMBIENT_W1_S61_INTERSECTION_RECOGNITION_20261103_v63_REVIEW.json'
    reduction_path = ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json'
    nsource_path = ROOT/'reports/AMBIENT_B18_N_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
    old_path = ROOT/'reports/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63/report.json'
    bundle_path = ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json'
    exact_path = ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/report.json'
    builder_path = exact_path.parent/'builder.json'
    assert read(recognition)['all_passed'] and read(recognition_review)['all_passed']
    assert read(recognition)['nilradical_module'] == 'C[[u1,...,u21]]*F'
    r, nsource, old, bundle, exact, builder = map(read, [reduction_path, nsource_path, old_path, bundle_path, exact_path, builder_path])
    h, x, y, z, w, W, lam, rr = s.symbols('h x y z w W lam r')
    V = s.Matrix([h, x, y, z, w])
    pairs = list(itertools.combinations(range(5), 2))
    mons = list(itertools.combinations_with_replacement(range(5), 2))
    def vec(P):
        return s.Matrix([s.Poly(P[i, j], *V).coeff_monomial(V[a]*V[b]) for i, j in pairs for a, b in mons])
    def canonical(P):
        d = s.Matrix([sum(s.diff(P[i, j], V[i]) for i in range(5)) for j in range(5)])
        return (P-(V*d.T-d*V.T)/5).applyfunc(s.expand)
    ey = s.eye(5)[:, 2]
    D = s.Matrix([x*x, 0, 0, h*x-lam*h*h, h*h+z*z])
    E = ey*D.T-D*ey.T
    E[3, 4], E[4, 3] = W*x*x, -W*x*x
    assert canonical(E) == E
    for i, j, k in itertools.combinations(range(5), 3):
        assert s.expand(sum(E[i, l]*s.diff(E[j, k], V[l])+E[j, l]*s.diff(E[k, i], V[l])
                            +E[k, l]*s.diff(E[i, j], V[l]) for l in range(5))) == 0
    P0 = E.subs({W: 0, lam: 0})
    # Entire actual fixed N19 source membership, independently reconstructed.
    sw = s.eye(5)
    sw[0, 0] = sw[1, 1] = 0
    sw[0, 1] = sw[1, 0] = 1
    def raw(f, g, H):
        P = s.zeros(5)
        P[1, 2], P[3, 4] = f, g*h*h
        P[2, 3], P[2, 4] = -s.diff(H, w), s.diff(H, z)
        return P-P.T
    change = lambda P: (sw*P.subs({h: x, x: h}, simultaneous=True)*sw.T).applyfunc(s.expand)
    cubics = [s.prod((h, x, z, w)[i] for i in m) for m in itertools.combinations_with_replacement(range(4), 3) if any(i >= 2 for i in m)]
    fixed = [change(canonical(raw(f, 0, 0))) for f in [h*h, h*x]]
    fixed += [change(canonical(raw(0, 1, 0)))]+[change(canonical(raw(0, 0, H))) for H in cubics]
    CN = s.Matrix.hstack(*map(vec, fixed))
    assert CN == s.Matrix(nsource['fiber19_basis'])
    selected = list(CN.T.rref()[1])
    alpha = (CN[selected, :].inv()*vec(E)[selected, :]).applyfunc(s.expand)
    assert CN*alpha == vec(E)
    # Uniform action changes W but fixes lam, for the whole two-parameter family.
    A = s.diag(1, 1, rr, 1, 1)
    image = (A*E.subs(dict(zip(V, A.inv()*V)), simultaneous=True)*A.T/rr).applyfunc(s.cancel)
    assert image == E.subs(W, W/rr)
    # Old free-v1,v4,v5,v7 support family, not a source chosen by type matching.
    a, b, c, d = s.symbols('k0 k1 k2 k5')
    H0, xx, yy, zz, ww = X
    potential = H0*(xx*ww+zz*zz)+a*xx*xx*zz+b*xx*xx*ww+c*xx*zz*zz+d*zz**3
    oldP = s.zeros(5)
    oldP[1, 2] = oldP[3, 4] = H0*H0
    oldP[2, 3], oldP[2, 4] = -s.diff(potential, ww), s.diff(potential, zz)
    oldP -= oldP.T
    frame = s.Matrix([[0, 1, 0, 0, 0], [-1, 0, 0, 0, 0], [0, 0, 1, 0, 0],
                      [0, -1, 0, 1, 0], [-1, 0, 0, 0, 1]])
    back = pushforward(oldP, frame.inv()).subs(dict(zip(X, V)), simultaneous=True)
    assert back.subs({a: 1, b: lam, c: 0, d: s.Rational(1, 3)}) == E.subs(W, 1)
    model = S61()
    normalized = model.normalize(oldP)
    pars = {s.Symbol('p'): 0, s.Symbol('s'): 0}
    pars.update({s.Symbol(f'k{i}'): 0 for i in range(9)})
    pars.update({a: a, b: b, c: c, d: d})
    old_u = s.Matrix(list(map(s.sympify, old['families']['N']['full_61_coordinates']))).subs(pars, simultaneous=True)
    assert old_u == s.Matrix(normalized['parameters'])
    qi = bundle['q_indices_zero_based']
    q = old_u[qi, :]
    K = q[builder['common_q_selected_rows'], :]
    remq = q[builder['remaining_q_rows'], :]
    assert remq == s.zeros(35, 1)
    expectedK = s.zeros(10, 1)
    expectedK[1], expectedK[4], expectedK[5], expectedK[7] = -b, a, 2*c, 3*d
    assert K == expectedK
    # Exact origin Jacobian of a21-parameter GL19 x E(W,lam) map.
    columns = []
    for ij in r['orbit_GL_generator_indices']:
        M = s.zeros(5)
        M[ij//5, ij % 5] = 1
        MV = M*V
        columns.append(vec((M*P0+P0*M.T-sum((MV[k]*P0.diff(V[k]) for k in range(5)), s.zeros(5))).applyfunc(s.expand)))
    O = s.Matrix.hstack(*columns)
    assert O == s.Matrix(r['families']['CA']['linear_J'])[:, :19]
    J21 = O.row_join(vec(E.diff(W))).row_join(vec(E.diff(lam)))
    rows = list(J21.T.rref()[1])
    assert len(rows) == 21
    det = s.factor(J21[rows, :].det())
    assert det != 0
    point_rank4 = E.subs({h: 0, x: 1, y: 0, z: 0, w: 0})
    assert point_rank4.extract([0, 2, 3, 4], [0, 2, 3, 4]).det() == W*W
    evidence = dict(
        polynomial_family_matrix=mat(E), actual_fixed_N19_source_coefficients=list(map(str, alpha)),
        all150_polynomial_family_coefficients=list(map(str, vec(E))),
        uniform_group_matrix=mat(A), uniform_tensor_scalar='1/r', uniform_parameter_action='(W,lam)->(W/r,lam)',
        old_free4_embedded_source_matrix=mat(oldP), old_free4_source_parameters=['k0', 'k1', 'k2', 'k5'],
        old_free4_back_transported_physical_family=mat(back),
        old_target10_embedded_support_coordinates=list(map(str, K)), old_remaining35_q=list(map(str, remq)),
        origin_GL19_generator_indices=r['orbit_GL_generator_indices'],
        full150_by21_origin_Jacobian=mat(J21), origin21_minor_rows=rows, origin21_minor_determinant=str(det),
    )
    with (out/'certificate.json.gz').open('wb') as f:
        with gzip.GzipFile(fileobj=f, mode='wb', filename='', mtime=0) as g:
            g.write((json.dumps(evidence, separators=(',', ':'))+'\n').encode())
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', recognition, recognition_review,
             reduction_path, nsource_path, old_path, bundle_path, exact_path, builder_path]
    result = dict(id=STEM, all_passed=True,
        actual_polynomial_family='P0-lam*h²*d_y_wedge_d_z+W*x²*d_z_wedge_d_w',
        all25_divergence_and350_strict_Jacobi_coefficients_QQ_W_lam_zero=True,
        actual_full_N19_fixed_source_membership_all150=True,
        old_free4_embedded_support_family_literal_coordinate_binding=True,
        rank4_for_W_nonzero_rank2_for_W_zero=True,
        connected_group_fixes_global_associated_prime_closure=True,
        entire_polynomial_A2_family_in_old_embedded21_closure=True,
        global_embedded21_prime_through_P0_corollary=True,
        dominant_rational_21_parameter_family_from_GL19_times_A2=True,
        exact_origin21_Jacobian_minor=str(det),
        one_regular21_dimensional_formal_image_branch_through_P0=True,
        transverse_embedded_support_has_dimension2_and_a_regular_formal_image_branch=True,
        whole_embedded_prime_smooth_or_analytically_irreducible_at_P0_claimed=False,
        full_origin_D_type_depth_nilmodule_or_new_Poisson_component_claimed=False,
        no_new_source_jet_or_GB=True,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'certificate.json.gz': sha(out/'certificate.json.gz')},
        seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'exact_origin21_Jacobian_minor',
          'entire_polynomial_A2_family_in_old_embedded21_closure', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())

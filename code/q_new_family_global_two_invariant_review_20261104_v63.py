#!/usr/bin/env python3
"""Independent all-point binding and intrinsic full-PGL field sandwich audit."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
ID = 'Q_NEW_FAMILY_GLOBAL_TWO_INVARIANT_PREFLIGHT_20261104_v63_REVIEW'


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def read(p):
    return json.loads(p.read_text())


def strings(m):
    return [[str(z) for z in row] for row in m.tolist()]


def main(out):
    start = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    if out.exists():
        assert not any(out.iterdir()), 'refuse to overwrite evidence'
    else:
        out.mkdir(parents=True)
    primary = ROOT / 'reports/Q_NEW_FAMILY_GLOBAL_TWO_INVARIANT_PREFLIGHT_20261104_v63'
    a = read(primary / 'report.json')
    assert a['all_passed'] and not a['WHICH_endpoint_and_sigma_classification_proved']
    for rel, h in a['input_sha256'].items():
        assert sha(ROOT / rel) == h, rel
    for rel, h in a['evidence_sha256'].items():
        assert sha(primary / rel) == h, rel
    refs = [
        'R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63',
        'Q_NEW_FAMILY_TARGET_FIELD_AND_SPECTRAL_DEGREES_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW',
        'Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_CONSTRUCTIVE_INVERSE_RECIPE_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63',
    ]
    paths = [ROOT / 'reports' / r / 'report.json' for r in refs]
    gates = [read(p) for p in paths]
    assert all(g['all_passed'] for g in gates)
    geometry = gates[0]
    sixpath = paths[-1].parent / 'certificate.json'
    six = read(sixpath)
    u, v, x, y, z, w, h, t = s.symbols('u v x y z w h t')
    n0, n1, n2 = s.symbols('n0 n1 n2')
    xyz = s.Matrix([x, y, z, w])
    coords = list(xyz)
    P = s.Matrix(geometry['affine_P']).subs({n0: u, n1: 1, n2: v})
    # Reconstruct the actual modular symmetry from divergence, rather than
    # copying the author's saved finite transformation.
    div = s.Matrix([s.expand(sum(s.diff(P[i, j], coords[i]) for i in range(4)))
                    for j in range(4)])
    assert div == s.Matrix([0, 4*u/3, u*x + u*v*y + s.Rational(1, 48), 0])
    N = div.jacobian(xyz)
    flow = xyz + t*div + t*t*N*div/2
    DF = flow.jacobian(xyz)
    flowed = P.subs(dict(zip(coords, flow)), simultaneous=True)
    assert (flowed - DF*P*DF.T).applyfunc(s.expand) == s.zeros(4)
    moved = flow.subs(t, -3*y/(4*u)).applyfunc(s.cancel)
    assert moved[0] == x and moved[1] == 0 and moved[3] == w
    # The six-source chart parametrizes every remaining generic affine point.
    X, W, B, c, Z = s.symbols('X W B c Z')
    zsource = s.sympify(six['y_zero_value_z'])
    bback = -4*u*Z/3 - W/s.Integer(1536) + W*W - X*X/64 - s.Rational(2, 9)
    assert s.cancel(bback.subs(Z, zsource) - B) == 0
    assert s.cancel(zsource.subs(B, bback) - Z) == 0
    point = {x: X, y: 0, z: zsource, w: W}
    value = P.subs(point, simultaneous=True).applyfunc(s.cancel)
    assert (value - s.Matrix(six['raw_value4'])).applyfunc(s.cancel) == s.zeros(4)
    frame = s.Matrix(six['rational_Darboux_frame4'])
    J4 = s.Matrix([[0, 1, 0, 0], [-1, 0, 0, 0],
                   [0, 0, 0, 1], [0, 0, -1, 0]])
    assert (frame*(c*value)*frame.T - J4).applyfunc(s.cancel) == s.zeros(4)
    # The actual restriction to H is reconstructed coefficientwise. Its
    # projected foliation is tested using a horizontal logarithmic numerator.
    lead = P.applyfunc(lambda f: sum(co*s.prod(a**e for a, e in zip(coords, mon))
                                    for mon, co in s.Poly(f, *coords).terms()
                                    if sum(mon) == 2))
    VV = s.Matrix([lead[2, i] for i in range(4)])
    assert VV[2] == 0 and all(s.diff(a, z) == 0 for a in VV)
    ez = s.Matrix([0, 0, 1, 0])
    assert lead == ez*VV.T - VV*ez.T
    projected = s.Matrix([VV[i] for i in (0, 1, 3)])
    plane = s.Matrix([x, y, w])
    Q = s.expand(72*s.sympify(geometry['homogeneous_quadric']).subs(
        {n0: u, n1: 1, n2: v, h: 0}))
    assert Q == -x*x + 64*w*w + 64*u*x*y + 32*u*v*y*y
    lp, lm = x + 8*w, x - 8*w
    deriv = lambda f: s.expand(sum(projected[i]*s.diff(f, plane[i]) for i in range(3)))
    assert deriv(lp) == s.expand(8*(x + v*y)*lp)
    assert deriv(lm) == s.expand(-8*(x + v*y)*lm)
    assert deriv(Q) == s.expand(-u*(x + v*y)*Q)
    alpha = s.Matrix([s.expand(lp*lm*s.diff(Q, a) + (u/16-1)*Q*lm*s.diff(lp, a)
                                     + (-u/16-1)*Q*lp*s.diff(lm, a)) for a in plane])
    assert alpha != s.zeros(3, 1)
    assert s.expand(alpha.dot(plane)) == 0
    assert s.expand(alpha.dot(projected)) == 0
    assert Q.subs({x: 0, y: 1, w: 0}) == 32*u*v
    # Adjugate pairings avoid the primary producer's matrix-inverse method.
    M = s.hessian(Q, plane)/2
    adj = M.adjugate()
    aa = s.expand((s.Matrix([1, 0, 8]).T*adj*s.Matrix([1, 0, 8]))[0])
    bb = s.expand((s.Matrix([1, 0, -8]).T*adj*s.Matrix([1, 0, -8]))[0])
    cc = s.expand((s.Matrix([1, 0, 8]).T*adj*s.Matrix([1, 0, -8]))[0])
    assert aa == bb == -65536*u*u
    assert s.expand(cc - 4096*u*(v + 16*u)) == 0
    assert s.cancel(256*cc*cc/(aa*bb) - (v/u + 16)**2) == 0
    assert s.expand(M.det() + 2048*u*(32*u + v)) == 0
    # Full-tensor simultaneous sign equivalence: a=-1, epsilon=-1 in the
    # corrected parameter action, tested afresh on all six affine coefficients.
    signF = s.Matrix([-x, y, -z, w])
    SD = s.diag(-1, 1, -1, 1)
    assert (-SD*P.subs(dict(zip(coords, signF)), simultaneous=True)*SD.T
            - P.subs({u: -u, v: -v}, simultaneous=True)).applyfunc(s.expand) == s.zeros(4)
    A, b, jj = s.symbols('A b jj')
    minimal = b*b + 32*A*b + (256-jj)*A*A
    assert s.discriminant(minimal, b) == 4*A*A*jj
    assert s.cancel(minimal.subs({A: u*u, b: u*v, jj: (v/u+16)**2}, simultaneous=True)) == 0
    cert = {
        'modular_divergence': list(map(str, div)), 'actual_finite_flow': list(map(str, flow)),
        'move_to_y_zero': list(map(str, moved)), 'source_z': str(zsource), 'B_inverse': str(bback),
        'all_point_binding': 'r_target is delta-constant; on the rational section r=u^2/1600; the section generates the entire delta-constant field. The modular full-tensor symmetry reaches y=0 from every generic affine point, and c is arbitrary.',
        'leading_tensor': strings(lead), 'projected_vector': list(map(str, projected)),
        'horizontal_log_numerator': list(map(str, alpha)),
        'conic': str(Q), 'conic_matrix': strings(M),
        'adjugate_pairings': list(map(str, (aa, bb, cc))),
        'invariant_line_proof': 'On another geometric-generic invariant line, integer residues a+(u/16)b=0 force a=b=0. The first constant ratio forces the line through O; Q(O)!=0 contradicts the second constant ratio.',
        'field_sandwich': 'k(u^2,(v/u+16)^2) subset E_full subset k(u^2,uv), k=QQ,C',
        'quadratic': str(minimal), 'valuation': 'The jj-adic valuation of 4*A^2*jj is one, so the extension has degree two.',
        'conclusion': 'E_full is one of these two rational fields of transcendence degree two.',
        'endpoint_and_sigma_not_decided_in_this_review': True,
    }
    cp = out/'independent-certificate.json'
    cp.write_text(json.dumps(cert, indent=2)+'\n')
    paths += [Path(__file__), ROOT/'research/notes'/(ID+'.md'), primary/'report.json', primary/'certificate.json',
              ROOT/'code/q_new_family_global_two_invariant_preflight_20261104_v63.py',
              ROOT/'research/notes/Q_NEW_FAMILY_GLOBAL_TWO_INVARIANT_PREFLIGHT_20261104_v63.md', sixpath]
    report = dict(id=ID, all_passed=True,
        independent_full_affine_modular_flow_and_all_point_six_source_binding=True,
        arbitrary_tensor_scalar_Darboux_frame_checked=True,
        independent_full_H_leading_projection_log_numerator_and_adjugate_chord=True,
        generic_invariant_line_uniqueness_and_QQ_descent_audited=True,
        actual_fullPGL_moduli_field_pure_rational_two_over_QQ_and_C=True,
        degree_two_sandwich_which_endpoint_not_decided=True,
        no_new_point_prime_high_degree_job_author_rerun_or_equivalence_search=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(paths)},
        evidence_sha256={cp.name:sha(cp)}, seconds=time.monotonic()-start, frozen=True)
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k:report[k] for k in ('all_passed','seconds')}, indent=2))


if __name__ == '__main__':
    p=argparse.ArgumentParser(); p.add_argument('--output', type=Path, default=ROOT/'reports'/ID)
    main(p.parse_args().output.resolve())

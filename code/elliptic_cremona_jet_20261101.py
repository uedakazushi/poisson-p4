#!/usr/bin/env python3
"""Universal quadratic-Cremona identification of the two elliptic first-jet covers.

Exact polynomial identities over QQ[lambda,x0,...,x4] and QQ[lambda,p1,...,p4]
prove the Poisson map and primitive first-jet transformation. Together with the
certified cubic inverse and equal parameter redundancies, this identifies the
two generic component degrees. The companion note explains the degree-one
conclusion using the separate q(5,2) reconstruction proof.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
from pathlib import Path
import time

import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

from elliptic_six_variable_fiber import affine_family, normalized_jet, transform_jet, J
from elliptic_weight1_dominance import symplectic_frame

ROOT = Path(__file__).resolve().parents[1]


def tr(a):
    return list(map(list, zip(*a, strict=True)))


def mul(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def plus(a, b):
    return [[x+y for x, y in zip(ar, br, strict=True)]
            for ar, br in zip(a, b, strict=True)]


def scale(c, a):
    return [[c*x for x in row] for row in a]


def determinant(a):
    if len(a) == 1:
        return a[0][0]
    return sum((-1)**j*a[0][j]*determinant(
        [[a[i][k] for k in range(len(a)) if k != j] for i in range(1, len(a))])
               for j in range(len(a)))


def adjugate(a):
    n = len(a)
    return [[(-1)**(i+j)*determinant(
        [[a[r][c] for c in range(n) if c != i] for r in range(n) if r != j])
             for j in range(n)] for i in range(n)]


def family_bar(kind, x, lam):
    """lambda^3 times the exact homogeneous q(5,k) bivector."""
    if kind == 1:
        a1, a2, a3 = (-3*lam**5+1)/5, -2*lam**2, lam
        b1, b2, b3 = (-lam**5-3)/5, 2*lam**3, lam**4
    else:
        a1, a2, a3 = (2*lam**5+1)/5, lam**4, -lam**2
        b1, b2, b3 = (-lam**5+2)/5, -lam, lam**3
    out = [[lam.ring.zero for _ in range(5)] for _ in range(5)]
    for i in range(5):
        j = (i+1) % 5
        out[i][j] = a1*x[i]*x[j]+a2*x[(i+2) % 5]*x[(i+4) % 5]+a3*x[(i+3) % 5]**2
        out[j][i] = -out[i][j]
        j = (i+2) % 5
        out[i][j] = b1*x[i]*x[j]+b2*x[(i+3) % 5]*x[(i+4) % 5]+b3*x[(i+1) % 5]**2
        out[j][i] = -out[i][j]
    return out


def exact_polynomial_checks():
    started = time.monotonic()
    base, lam, *x = ring('lam,x0,x1,x2,x3,x4', QQ)
    q = [lam*x[i]**2+lam**2*x[(i-1) % 5]*x[(i+1) % 5]
         - x[(i-2) % 5]*x[(i+2) % 5] for i in range(5)]
    jac = [[f.diff(coordinate) for coordinate in x] for f in q]
    b1 = family_bar(1, x, lam)
    b2_at_q = family_bar(2, q, lam)
    difference = plus(mul(mul(jac, b1), tr(jac)), b2_at_q)
    for i, j, k in combinations(range(5), 3):
        # Euler-wedge zero means equality of induced projective bivectors.
        assert q[i]*difference[j][k]-q[j]*difference[i][k]+q[k]*difference[i][j] == 0
    affine = lambda f: f.evaluate(x[0], 1)
    z = [affine(coordinate) for coordinate in x[1:]]
    affine_ring = z[0].ring
    zero = affine_ring.zero
    k = [[affine(b1[i+1][j+1]-x[i+1]*b1[0][j+1]+x[j+1]*b1[0][i+1])
          for j in range(4)] for i in range(4)]
    pf = k[0][1]*k[2][3]-k[0][2]*k[1][3]+k[0][3]*k[1][2]
    n = [[zero, k[2][3], -k[1][3], k[1][2]],
         [-k[2][3], zero, k[0][3], -k[0][2]],
         [k[1][3], -k[0][3], zero, k[0][1]],
         [-k[1][2], k[0][2], -k[0][1], zero]]
    kn = mul(k, n)
    assert all(kn[i][j] == (-pf if i == j else 0) for i in range(4) for j in range(4))
    # C_a / pf^2 is the intrinsic primitive tensor for lambda^3 P1.
    c = [[[-pf*n[b][d].diff(z[a])+QQ(3, 5)*pf.diff(z[a])*n[b][d]
            - QQ(1, 5)*n[b][a]*pf.diff(z[d])+QQ(1, 5)*pf.diff(z[b])*n[d][a]
           for d in range(4)] for b in range(4)] for a in range(4)]
    j = [[affine(f) for f in row] for row in jac]
    adj = adjugate(j)
    det = determinant(j)
    assert det == 2*pf
    product = mul(j, adj)
    assert all(product[i][h] == (det if i == h else 0) for i in range(5) for h in range(5))
    ap = [[-z[i]]+[affine_ring.one if i == h else zero for h in range(4)] for i in range(4)]
    hessians = [[[affine(q[i].diff(x[a+1]).diff(x[b+1])) for b in range(4)]
                 for i in range(5)] for a in range(4)]
    g = [scale(QQ(1, 2), mul(mul(ap, adj), h)) for h in hessians]
    traces = [sum(m[i][i] for i in range(4)) for m in g]
    for a in range(4):
        raw = scale(-1, plus(mul(tr(g[a]), n), mul(n, g[a])))
        primitive = [[raw[b][d]+QQ(2, 5)*traces[a]*n[b][d]
                      + QQ(1, 5)*n[b][a]*traces[d]-QQ(1, 5)*traces[b]*n[d][a]
                      for d in range(4)] for b in range(4)]
        assert all(primitive[b][d] == QQ(3, 2)*c[a][b][d]
                   for b in range(4) for d in range(4))
    return {'field': 'QQ', 'projective_Poisson_identity_components': 10,
            'homogeneous_family_rescaling': 'Bbar_k=lambda^3 B_k; Qbar=lambda Q',
            'Jacobian_determinant_identity': 'det Jac(Qbar)=2 Pf(Pbar_1)',
            'inverse_form_adjugate_entries': 16,
            'Jacobian_adjugate_entries': 25,
            'universal_primitive_tensor_entries': 64,
            'raw_primitive_tensor_identity': 'f^* beta_2 = beta_1/2',
            'all_polynomial_identities_zero': True,
            'elapsed_seconds': time.monotonic()-started}


def normalized_frame_checks():
    lam, z, p1, _ = affine_family(1)
    _, _, p2, _ = affine_family(2)
    xx = [s.Integer(1), *z]
    quadrics = [xx[i]**2+lam*xx[(i-1) % 5]*xx[(i+1) % 5]
                - xx[(i-2) % 5]*xx[(i+2) % 5]/lam for i in range(5)]
    f = s.Matrix([quadrics[i]/quadrics[0] for i in range(1, 5)])
    df = f.jacobian(z)
    cases = []
    for parameter, point, scalar in [(2, [1, 2, 3, 4], 1), (3, [2, 3, 4, 5], 5)]:
        scalar = s.Integer(scalar)
        subst = dict(zip(z, point, strict=True)) | {lam: s.Integer(parameter)}
        image, derivative = f.subs(subst), df.subs(subst)
        target_subst = dict(zip(z, image, strict=True)) | {lam: s.Integer(parameter)}
        k1 = scalar*p1.subs(subst)
        k2 = -4*scalar*p2.subs(target_subst, simultaneous=True)
        f1, f2 = symplectic_frame(k1), symplectic_frame(k2)
        matrix = -2*f2*derivative*f1.inv()
        assert matrix*J*matrix.T == J
        beta1 = normalized_jet(1, parameter, point, scalar)
        beta2 = normalized_jet(2, parameter, list(image), -4*scalar)
        assert transform_jet(beta1, matrix) == beta2
        cases.append({'lambda': parameter, 'point': point, 'scalar': str(scalar),
                      'image_point': [str(v) for v in image], 'target_scalar': str(-4*scalar),
                      'frame_relation': 'S0=-2 F2 Df F1^-1; F2 normalizes -4c P2',
                      'S0': [[str(v) for v in row] for row in matrix.tolist()],
                      'symplectic_and_jet_identities': True})
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = {'id': 'ELLIPTIC_CREMONA_JET_20261101',
              'generated_utc': datetime.now(timezone.utc).isoformat(),
              'proof_note': 'research/notes/ELLIPTIC_GENERIC_DEGREES_20261101.md',
              'all_passed': True, 'universal_checks': exact_polynomial_checks(),
              'frame_tests': normalized_frame_checks(),
              'parameter_transfer': '(lambda,p,c,S) -> (lambda,f(p),-4c,S*S0^-1)',
              'dependencies_for_degree_conclusion': [
                  'research/notes/ELLIPTIC_Q52_RECONSTRUCTION_REVIEW_20261101.md',
                  'reports/HESSE_QUADRO_CUBIC_INVERSE_20261101.json',
                  'research/notes/ELLIPTIC_PARAMETER_REDUNDANCY_20261101.md'],
              'scope': 'The universal identities identify the two rational first-jet parameter covers. '
                       'Degree one of T40 and the equal3000-fold redundancies then give degree one of T39. '
                       'No fullS61 exhaustion or finiteness assertion.',
              'input_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                               for path in [Path(__file__).resolve(),
                                            ROOT/'code/elliptic_weight1_dominance.py',
                                            ROOT/'code/elliptic_six_variable_fiber.py']}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print('PASS: universal Poisson and primitive-jet Cremona identities; rational covers identified.')


if __name__ == '__main__':
    main()

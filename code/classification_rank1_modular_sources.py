#!/usr/bin/env python3
"""Bounded exact identities for complete rank-one modular source reductions."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import sympy as S

ROOT = Path(__file__).resolve().parents[1]


def zero(expr):
    return S.cancel(S.expand(expr)) == 0


def matrix(entries, n=4):
    p = S.zeros(n)
    for (i, j), value in zip(itertools.combinations(range(n), 2), entries):
        p[i, j], p[j, i] = value, -value
    return p


def jacobi(p, coords):
    return [S.expand(sum(p[i, a]*S.diff(p[j, k], coords[a])
                         +p[j, a]*S.diff(p[k, i], coords[a])
                         +p[k, a]*S.diff(p[i, j], coords[a])
                         for a in range(len(coords))))
            for i, j, k in itertools.combinations(range(len(coords)), 3)]


def divergence(p, coords):
    return [S.expand(sum(S.diff(p[i, j], coords[i])
                         for i in range(len(coords)))) for j in range(len(coords))]


def pf(p):
    return S.expand(p[0, 1]*p[2, 3]-p[0, 2]*p[1, 3]+p[0, 3]*p[1, 2])


def polynomial(coords, degree, prefix, exact=False):
    mons = [S.prod(coords[i] for i in m)
            for d in ([degree] if exact else range(degree+1))
            for m in itertools.combinations_with_replacement(range(len(coords)), d)]
    coefficients = S.symbols(f'{prefix}0:{len(mons)}')
    return sum(c*m for c, m in zip(coefficients, mons)), coefficients


def main(output):
    started = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (1536*1024**2, 1536*1024**2))
    signal.alarm(60)
    x, y, z, w, h = S.symbols('x y z w h')
    A, d0, f0, d = S.symbols('A d0 f0 d')
    k, ks = polynomial((x, w), 2, 'k')
    zx = x*(A*z+d0+f0*x+d*w/3)
    zw = w*(-2*A*z-2*d0+(-1-2*f0)*x-2*d*w/3)
    zz = A*z*z/2+z*(d0+d*w)+k
    p = matrix([-zx, x, 0, zz+y, zw, 2*w])
    assert all(zero(e) for e in jacobi(p, (x, y, z, w)))
    assert divergence(p, (x, y, z, w)) == [0, x, 0, 0]
    assert zero(pf(p)-x*x*w)
    ten_source = dict(parameters=list(map(str, (A, d0, f0, d, *ks))),
                      coefficients=[str(p[i, j]) for i, j in itertools.combinations(range(4), 2)],
                      universal_jacobi_zero=True, modular='x*d/dy', pfaffian='x**2*w')

    # Universal affine source with constant modular field, retaining both residuals.
    t, r, s = S.symbols('t r s')
    u, _ = polynomial((s, w), 1, 'u')
    v, _ = polynomial((s, w), 1, 'v')
    g, _ = polynomial((s, w), 1, 'g')
    F, _ = polynomial((s, w), 2, 'F')
    H, _ = polynomial((s, w), 3, 'H')
    divu = S.diff(u, s)+S.diff(v, w)
    fun = r*g+H
    pp = matrix([-F-(1+divu)*r, u, v,
                 S.diff(fun, w), -S.diff(fun, s), g])
    residual0 = S.expand(u*S.diff(g, s)+v*S.diff(g, w)-divu*g)
    residual1 = S.expand(u*S.diff(H, s)+v*S.diff(H, w)+H-F*g)
    residual = r*residual0+residual1
    expected = [S.diff(residual, w), -S.diff(residual, s),
                S.diff(residual, r), 0]
    assert all(zero(a-b) for a, b in zip(jacobi(pp, (t, r, s, w)), expected))
    assert divergence(pp, (t, r, s, w)) == [1, 0, 0, 0]
    assert zero(pf(pp) - (residual-fun))
    transformed = []
    for branch, old in [('chain_top', (t, r+t*t/2, s, w)),
                        ('trivial_block', (s, r+s*t, t, w))]:
        inverse = ({t:x, r:y-x*x/2, s:z, w:w} if branch == 'chain_top'
                   else {t:z, r:y-x*z, s:x, w:w})
        derivative = S.Matrix(old).jacobian((t, r, s, w))
        original = (derivative*pp*derivative.T).subs(inverse, simultaneous=True).applyfunc(S.expand)
        assert all(S.Poly(original[i, j], x, y, z, w).total_degree() <= 2
                   for i, j in itertools.combinations(range(4), 2))
        target = [1, x, 0, 0] if branch == 'chain_top' else [0, x, 1, 0]
        assert divergence(original, (x, y, z, w)) == target
        transformed.append(dict(branch=branch, all_six_coefficients=[str(original[i, j])
                               for i, j in itertools.combinations(range(4), 2)],
                               degree_at_most_two=True, exact_modular=list(map(str, target))))

    # Explicit projective CB recognition on the dense a*b != 0 chart.
    a, b, HH, XX = S.symbols('a b H0 X0', nonzero=True)
    potential, _ = polynomial((h, x, z, w), 3, 'c', exact=True)
    G, L = a*x+b*h, 2*a*x-b*h
    raw = S.zeros(5)
    raw[1, 2] = G*L/(6*a*a)
    raw[3, 4] = h*G
    raw[2, 3] = -L*z/(12*a)-S.diff(potential, w)
    raw[2, 4] = -L*w/(12*a)+S.diff(potential, z)
    raw -= raw.T
    oldcoords = S.Matrix([h, x, y, z, w])
    newcoords = S.Matrix([L, G, y, z, w])
    derivative = newcoords.jacobian(oldcoords)
    change = {h:(2*XX-HH)/(3*b), x:(XX+HH)/(3*a)}
    moved = (derivative*raw*derivative.T).subs(change, simultaneous=True)
    ee = S.Matrix([HH, XX, y, z, w])
    vv = S.Matrix([0, 0, XX/(3*a), 0, 0])
    moved = (moved-ee*vv.T+vv*ee.T).applyfunc(S.expand)
    newf = XX*(HH-2*XX)/(6*a)
    newpotential = potential.subs(change, simultaneous=True)
    expected_cb = S.zeros(5)
    expected_cb[1, 2] = newf
    expected_cb[3, 4] = -2*a*newf/b
    expected_cb[2, 3] = -S.diff(newf, XX)*z/2-S.diff(newpotential, w)
    expected_cb[2, 4] = -S.diff(newf, XX)*w/2+S.diff(newpotential, z)
    expected_cb -= expected_cb.T
    assert all(zero(v) for v in moved-expected_cb)
    report = dict(all_passed=True, task_id='E_multiple_hyperplane',
                  exact_domain='QQ with all displayed parameters; a and b inverted only for CB coordinate change',
                  ten_parameter_source=ten_source,
                  constant_modular_reduction=dict(jacobi_equals_jacobian_of_two_residuals=True,
                       residual0=str(residual0), residual1=str(residual1),
                       both_full_inverse_sources=transformed),
                  cb_recognition=dict(all_twenty_cubic_coefficients_retained=True,
                       all_ten_bivector_coefficients_identical=True,
                       f=str(newf), ratio='-2*a/b',
                       boundary_b_zero='retained by closure of the original irreducible source; not evaluated in rational coordinates'),
                  global_rank1_coverage_claimed=False, seconds=time.monotonic()-started)
    report['input_sha256']={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                           for path in [Path(__file__).resolve(), ROOT/'sources/v06.zip']}
    output.write_text(json.dumps(report, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k:v for k, v in report.items()
                      if k not in ['ten_parameter_source', 'constant_modular_reduction', 'input_sha256']}, indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)

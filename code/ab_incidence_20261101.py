#!/usr/bin/env python3
"""Intrinsic AB incidence: the fixed-L core is the minimal nilpotent sl3 cone."""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
import sympy as s
from toric_zero_jet_dimension_20261101 import encode, witness

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    x, y, z, w = s.symbols('x y z w')
    q = s.Matrix([x, z])
    Y = [x, y, z, w]
    ss, hh, vv = s.symbols('s0:3'), s.symbols('h0:5'), s.symbols('v0:3')
    S = ss[0]*x*x+ss[1]*x*z+ss[2]*z*z
    H = sum(hh[i]*x**(4-i)*z**i for i in range(5))
    J2 = s.Matrix([[0, 1], [-1, 0]])
    J = s.diag(J2, J2)
    F = q*S/4+J2*s.Matrix([s.diff(H, x), s.diff(H, z)])
    mixed = F.jacobian(q)
    vertical = vv[0]*x*x+vv[1]*x*z+vv[2]*z*z
    Q = s.zeros(4)
    for i, ii in enumerate((0, 2)):
        for j, jj in enumerate((1, 3)):
            Q[ii, jj] = mixed[i, j]
            Q[jj, ii] = -mixed[i, j]
    Q[1, 3], Q[3, 1] = vertical, -vertical
    pairs = list(combinations(range(4), 2))

    def jacobi(P, constant=None):
        out = []
        for i, j, k in combinations(range(4), 3):
            out.append(s.expand(sum((constant if constant is not None else P)[a, l]*s.diff(P[b, c], Y[l])
                                    for a, b, c in ((i, j, k), (j, k, i), (k, i, j)) for l in range(4))))
        return out

    # All compatible tensors constant along the fixed L=<dy,dw> form this 11-space.
    variables = s.symbols('r0:18')
    raw = s.zeros(4)
    for i, (a, b) in enumerate(pairs):
        raw[a, b] = sum(variables[3*i+j]*m for j, m in enumerate((x*x, x*z, z*z)))
        raw[b, a] = -raw[a, b]
    conditions = [c for poly in jacobi(raw, J) for c in s.Poly(poly, x, z).coeffs()]
    constraint = s.Matrix(conditions).jacobian(variables)
    assert constraint.rank() == 7
    assert jacobi(Q, J) == [0]*4
    coefficients = s.Matrix([s.Poly(s.expand(Q[i, j]), x, z).coeff_monomial(m)
                             for i, j in pairs for m in (x*x, x*z, z*z)])
    assert coefficients.jacobian([*ss, *hh, *vv]).rank() == 11

    def Nof(s_poly, h_poly):
        sp, hp = s.Poly(s.expand(s_poly), x, z), s.Poly(s.expand(h_poly), x, z)
        a = [sp.coeff_monomial(m) for m in (x*x, x*z, z*z)]
        h = [hp.coeff_monomial(x**(4-i)*z**i) for i in range(5)]
        LH = s.Matrix([[h[2]/6, -h[1]/4, h[0]], [h[3]/2, -h[2]/3, h[1]/2],
                       [h[4], -h[3]/4, h[2]/6]])
        Cs = s.Matrix([[2*a[1], -2*a[0], 0], [4*a[2], 0, -4*a[0]],
                       [0, 2*a[2], -2*a[1]]])
        return LH+Cs/96

    N = Nof(S, H)
    assert s.trace(N) == 0
    assert N.reshape(9, 1).jacobian([*ss, *hh]).rank() == 8
    covariance = []
    for i, j in ((0, 0), (0, 1), (1, 0), (1, 1)):
        B = s.zeros(2)
        B[i, j] = 1
        ds = s.Matrix([s.diff(S, x), s.diff(S, z)]).dot(B*q)
        dh = s.Matrix([s.diff(H, x), s.diff(H, z)]).dot(B*q)-s.trace(B)*H
        rho = s.Matrix.hstack(*[s.Matrix([s.Poly(s.Matrix([s.diff(m, x), s.diff(m, z)]).dot(B*q), x, z).coeff_monomial(n)
                                         for n in (x*x, x*z, z*z)]) for m in (x*x, x*z, z*z)])
        derivative = Nof(ds, dh)
        assert (derivative-s.trace(B)*N-rho*N+N*rho).applyfunc(s.expand) == s.zeros(3)
        covariance.append(dict(gl2_generator=encode(B), quadratic_representation=encode(rho)))
    base_sub = dict(zip([*ss, *hh], [1, 0, 1, 0, s.Rational(1, 12), 0, -s.Rational(1, 12), 0]))
    base = N.subs(base_sub)
    assert base.rank() == 1 and base*base == s.zeros(3)
    orbit = s.Matrix.hstack(*[(s.trace(s.Matrix(c['gl2_generator']))*base
                              +s.Matrix(c['quadratic_representation'])*base
                              -base*s.Matrix(c['quadratic_representation'])).reshape(9, 1) for c in covariance])
    orbit_witness = witness(orbit)
    assert orbit_witness['rank'] == 4
    minors = [s.expand(N[list(r), list(c)].det()) for r in combinations(range(3), 2) for c in combinations(range(3), 2)]
    gb = s.groebner(minors, *ss, *hh, order='grevlex')
    jac = jacobi(Q)
    jac_coeff = [s.expand(c) for poly in jac for c in s.Poly(poly, *Y).coeffs()]
    assert all(gb.reduce(c)[1] == 0 for c in jac_coeff)

    # All K-dimension-three tensors in a fixed L chart have rank-one core.
    a, c, d = s.symbols('a c d')
    one_S = -a*x*x
    one_H = c*x**4/12-a*x**3*z/12
    one_N = Nof(one_S, one_H)
    assert all(s.expand(one_N[list(r), list(k)].det()) == 0
               for r in combinations(range(3), 2) for k in combinations(range(3), 2))
    one_F = q*one_S/4+J2*s.Matrix([s.diff(one_H, x), s.diff(one_H, z)])
    assert (one_F-s.Matrix([-a*x**3/3, -c*x**3/3])).applyfunc(s.expand) == s.zeros(2, 1)
    # A L-preserving change of complement leaves the mixed block unchanged.
    c0, c1, c2 = s.symbols('c0 c1 c2')
    C = s.Matrix([[c0, c1], [c1, c2]])
    shear = s.eye(4)
    for i, ii in enumerate((1, 3)):
        for j, jj in enumerate((0, 2)):
            shear[ii, jj] = C[i, j]
    assert shear*J*shear.T == J
    changed = shear*Q*shear.T
    assert changed[[0, 2], [1, 3]] == mixed
    paths = [Path(__file__), ROOT/'research/notes/TORIC_ZERO_JET_FIBER_20261101.md',
             ROOT/'research/notes/TORIC_ZERO_JET_EXHAUSTION_20261101.md']
    report = dict(id='AB_INCIDENCE_20261101', all_passed=True, field='QQ',
                  fixed_L_compatibility_rank=7, fixed_L_compatible_dimension=11,
                  cubic_vector=encode(F), mixed_block=encode(mixed), sl3_linear_identification=encode(N),
                  gl2_covariance_generators=covariance, base_rank_one_matrix=encode(base),
                  base_GL2_orbit_witness=orbit_witness, rank_one_minors=list(map(str, minors)),
                  rank_one_groebner_basis=[str(p.as_expr()) for p in gb.polys],
                  Jacobi_coefficients_reduce_to_zero=True, K_three_core=encode(one_N),
                  proper_model='Over LG(2,4), the fiber is the normal four-dimensional minimal nilpotent sl3 cone times affine3.',
                  geometric_fibers={'translation_kernel_dimension_2': 'one reduced point',
                                    'translation_kernel_dimension_3': 'P1, reduced', 'Q_zero': 'LG(2,4), reduced'},
                  singular_locus='Two irreducible six-dimensional loci B0 (constant Lagrangian plane factor) and B1 (one-variable tensor); intersection irreducible dimension five.',
                  limitations=['The proper incidence is not the normalization map: it has positive-dimensional fibers.',
                               'No global normality of AB or scheme structure of the actual T01 zero-jet fiber is asserted.'],
                  proof_note='research/notes/AB_INCIDENCE_20261101.md',
                  source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps(dict(all_passed=True, fixed_L_core='minimal_nilpotent_sl3_cone', fibers=['point', 'P1', 'LG'])))


if __name__ == '__main__':
    main()

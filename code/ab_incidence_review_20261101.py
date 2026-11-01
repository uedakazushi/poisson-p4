#!/usr/bin/env python3
"""Finite GL2 covariance and independent cubic-coordinate audit of AB incidence."""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def verify():
    x, z, y, w = s.symbols('x z y w')
    coords = [x, y, z, w]
    q = s.Matrix([x, z])
    a = s.symbols('a0:8')
    cubic = [x**3, x*x*z, x*z*z, z**3]
    F = s.Matrix([sum(a[i]*cubic[i] for i in range(4)),
                  sum(a[i+4]*cubic[i] for i in range(4))])
    S = s.diff(F[0], x)+s.diff(F[1], z)
    D = F-q*S/4
    H = s.expand((z*D[0]-x*D[1])/4)
    assert s.expand(s.diff(H, z)-D[0]) == 0 and s.expand(-s.diff(H, x)-D[1]) == 0
    def Nof(S, H):
        sp, hp = s.Poly(s.expand(S), x, z), s.Poly(s.expand(H), x, z)
        ss = [sp.coeff_monomial(v) for v in [x*x, x*z, z*z]]
        h = [hp.coeff_monomial(x**(4-i)*z**i) for i in range(5)]
        return s.Matrix([[h[2]/6+ss[1]/48, -h[1]/4-ss[0]/48, h[0]],
                         [h[3]/2+ss[2]/24, -h[2]/3, h[1]/2-ss[0]/24],
                         [h[4], -h[3]/4+ss[2]/48, h[2]/6-ss[1]/48]])
    N = Nof(S, H)
    assert s.trace(N) == 0 and N.reshape(9, 1).jacobian(a).rank() == 8
    tau, r, v = s.symbols('tau r v', nonzero=True)
    transformations = [s.Matrix([[1, tau], [0, 1]]), s.Matrix([[1, 0], [tau, 1]]), s.diag(r, v)]
    quadratic = [x*x, x*z, z*z]
    for T in transformations:
        sub = dict(zip(q, T*q))
        rho = s.Matrix.hstack(*[s.Matrix([s.Poly(m.subs(sub, simultaneous=True), x, z).coeff_monomial(n)
                                          for n in quadratic]) for m in quadratic])
        changed = Nof(S.subs(sub, simultaneous=True), H.subs(sub, simultaneous=True)/T.det())
        assert (changed-T.det()*rho*N*rho.inv()).applyfunc(s.cancel) == s.zeros(3)
    J = s.diag(s.Matrix([[0, 1], [-1, 0]]), s.Matrix([[0, 1], [-1, 0]]))
    g0, g1, g2 = s.symbols('g0 g1 g2')
    Q = s.zeros(4)
    mixed = F.jacobian(q)
    for i, qi in enumerate([0, 2]):
        for j, pj in enumerate([1, 3]):
            Q[qi, pj] = mixed[i, j]; Q[pj, qi] = -mixed[i, j]
    Q[1, 3] = g0*x*x+g1*x*z+g2*z*z; Q[3, 1] = -Q[1, 3]
    def jac(P, first):
        return [s.expand(sum(first[i, l]*s.diff(P[j, k], coords[l])
                           +first[j, l]*s.diff(P[k, i], coords[l])
                           +first[k, l]*s.diff(P[i, j], coords[l]) for l in range(4)))
                for i, j, k in combinations(range(4), 3)]
    assert jac(Q, J) == [0]*4
    rawvars = s.symbols('r0:18'); raw = s.zeros(4)
    for n, (i, j) in enumerate(combinations(range(4), 2)):
        raw[i, j] = sum(rawvars[3*n+k]*quadratic[k] for k in range(3)); raw[j, i] = -raw[i, j]
    rows = [coef for p in jac(raw, J) for coef in s.Poly(p, x, z).coeffs()]
    assert s.Matrix(rows).jacobian(rawvars).rank() == 7
    coeffs = [s.Poly(Q[i, j], x, z).coeff_monomial(mon) for i, j in combinations(range(4), 2) for mon in quadratic]
    assert s.Matrix(coeffs).jacobian([*a, g0, g1, g2]).rank() == 11
    minors = [s.expand(N[list(i), list(j)].det()) for i in combinations(range(3), 2) for j in combinations(range(3), 2)]
    gb = s.groebner(minors, *a, order='lex')
    jac_coeffs = [c for p in jac(Q, Q) for c in s.Poly(p, *coords).coeffs()]
    assert all(gb.reduce(c)[1] == 0 for c in jac_coeffs)
    base = dict.fromkeys(a, 0); base[a[0]] = base[a[7]] = s.Rational(1, 3)
    Nb = N.subs(base)
    assert Nb.rank() == 1 and Nb*Nb == s.zeros(3)
    orbit = []
    for T in transformations[:2]+[s.diag(r, 1), s.diag(1, r)]:
        sub = dict(zip(q, T*q)); changed = Nof(S.subs(sub, simultaneous=True), H.subs(sub, simultaneous=True)/T.det()).subs(base)
        variable = tau if changed.has(tau) else r
        orbit.append(changed.diff(variable).subs({tau: 0, r: 1}).reshape(9, 1))
    assert s.Matrix.hstack(*orbit).rank() == 4
    aa, cc = s.symbols('aa cc')
    one = {**dict.fromkeys(a, 0), a[0]: -aa/3, a[4]: -cc/3}
    assert all(s.expand(f.subs(one)) == 0 for f in minors)
    source = ROOT/'reports/AB_INCIDENCE_20261101.json'
    paths = [Path(__file__), source, ROOT/'research/notes/AB_INCIDENCE_20261101.md']
    return dict(id='AB_INCIDENCE_REVIEW_20261101', all_passed=True,
        independent_coordinates='eight cubic-vector coefficients; H recovered by Euler contraction',
        finite_GL2_generators_verified=['upper unipotent', 'lower unipotent', 'all diagonal matrices'],
        linear_rank=8, full_fixed_L_compatibility_rank=7, fixed_L_parameter_rank=11,
        generic_GL2_orbit_rank=4, all_Jacobi_coefficients_in_rank_one_ideal=True,
        alternative_lex_groebner_basis_size=len(gb.polys), all_K3_rank_one_identities=True,
        proof_note='research/notes/AB_INCIDENCE_REVIEW_20261101.md',
        input_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); result = verify()
    with args.output.open('x') as out:
        json.dump(result, out, indent=2); out.write('\n')
    print(json.dumps(dict(all_passed=True)))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Strict stabilizers and AB degenerations in both minimal-trace families."""
import argparse
import hashlib
import json
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_model import X, quotient_coefficients
from toric_zero_jet_dimension_20261101 import encode, witness

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    x0, x, y, z, w = X
    Y = s.Matrix([x, y, z, w])
    h = s.symbols('h0:9')
    G = sum(a*b for a, b in zip(h, [x*x*z, x*x*w, x*z*z, x*z*w, x*w*w,
                                    z**3, z*z*w, z*w*w, w**3]))
    pairs = list(combinations(range(5), 2))
    monomials = [X[i]*X[j] for i, j in combinations_with_replacement(range(5), 2)]

    def vector(matrix):
        return s.Matrix([s.Poly(s.expand(matrix[i, j]), *X).coeff_monomial(mon)
                         for i, j in pairs for mon in monomials])

    def global_action(P):
        cols = []
        for i in range(5):
            for j in range(5):
                a = s.zeros(5)
                a[i, j] = 1
                cols.append(s.Matrix(quotient_coefficients(a*P+P*a.T-X[j]*P.diff(X[i]))))
        return s.Matrix.hstack(*cols)

    def elementary(i, j):
        A = s.zeros(5)
        A[i, j] = 1
        return A

    records = {}
    a, b, k, lam = s.symbols('a b k lam')
    shear = s.Matrix([[1, 0, 0, 0], [0, 1, -b, a], [-a, 0, 1, 0], [-b, 0, 0, 1]])
    J = s.diag(s.Matrix([[0, 1], [-1, 0]]), s.Matrix([[0, 1], [-1, 0]]))
    assert shear*J*shear.T == J
    for kind, eps in (('M0', 1), ('M1', -1)):
        P = s.zeros(5)
        P[1, 2] = x0*x0+x*x
        P[3, 4] = x0*x0+(x*x if kind == 'M1' else 0)
        P[2, 3] = -s.diff(G, w)-(x*z if kind == 'M1' else 0)
        P[2, 4] = s.diff(G, z)-(x*w if kind == 'M1' else 0)
        P -= P.T
        # Three distinct divisor hyperplanes force rows x0,x to vanish modulo Euler.
        actions = []
        for i in range(2, 5):
            for j in range(5):
                A = elementary(i, j)
                actions.append(A*P+P*A.T-X[j]*P.diff(X[i]))
        full = s.Matrix.hstack(*[vector(A) for A in actions])
        constant = s.Matrix.hstack(*[s.Matrix([A[i, j].coeff(x0, 2) for i, j in pairs])
                                    for A in actions])
        K = s.Matrix.hstack(*constant.nullspace())
        assert K.cols == 9
        reduced = full*K
        indices = [(p, m) for p in pairs for m in monomials]
        linear_rows = [i for i, (_, m) in enumerate(indices) if s.degree(m, x0) == 1]
        quadratic_rows = [i for i, (_, m) in enumerate(indices) if s.degree(m, x0) == 0]
        linear = reduced[linear_rows, :]
        stack = s.Matrix.vstack(linear.subs(dict.fromkeys(h, 0)), *[linear.diff(v) for v in h])
        no_translation = s.Matrix.hstack(*stack.nullspace())
        assert no_translation.cols == 7
        quad = reduced[quadratic_rows, :]*no_translation
        stack = s.Matrix.vstack(quad.subs(dict.fromkeys(h, 0)), *[quad.diff(v) for v in h])
        M = quad[:, list(stack.rref()[1])]
        M = M[[i for i in range(M.rows) if any(M[i, :])], :]
        assert M.shape == (12, 5)
        distinct = M.subs(dict(zip(h[5:], [lam, 0, 0, lam])))
        assert s.factor(distinct[[3, 4, 5, 7], [0, 1, 3, 4]].det()) == 972*lam**4
        double = M.subs(dict(zip(h[5:], [0, 1, 0, 0])))
        assert double[[1, 3, 4, 9], [0, 1, 3, 4]].det() == 24
        triple = M.subs(dict(zip(h[5:], [1, 0, 0, 0])))
        rr, cc = [3, 7, 9], [1, 0, 4]
        assert triple[rr, cc].det() == -162
        residual = (triple[:, [2, 3]]-triple[:, cc]*triple[rr, cc].inv()*triple[rr, [2, 3]]).applyfunc(s.factor)
        assert residual[6, 0] == -h[3]-eps
        assert residual[0, 1] == -h[3]*(h[3]-eps)/3
        assert residual[0, 1].subs(h[3], -eps) == -s.Rational(2, 3)
        rank_one = M.subs(dict(zip(h[2:], [1, 0, 0, 0, 0, 0, 0])))
        assert rank_one[[0, 1, 6, 7], [0, 1, 2, 4]].det() == 8
        semisimple = M.subs(dict(zip(h, [0, 0, 0, k, 0, 0, 0, 0, 0])))
        assert s.factor(semisimple[[0, 2, 6, 7], [0, 1, 2, 3]].det()) == 4*k*k*(k-1)*(k+1)
        resonance_minor = None
        if kind == 'M0':
            resonance = M.subs(dict(zip(h, [0, k, 0, 1, 0, 0, 0, 0, 0])))
            resonance_minor = s.factor(resonance[[0, 2, 6, 7], [1, 2, 3, 4]].det())
            assert resonance_minor == -8*k
        base = P[1:5, 1:5].subs(x0, 1).subs(dict.fromkeys(h[5:], 0))
        changed = (shear*base.subs(dict(zip(Y, shear.inv()*Y)), simultaneous=True)*shear.T).applyfunc(s.expand)
        new_h0 = s.Poly(changed[1, 3]+(x*w if kind == 'M1' else 0), *Y).coeff_monomial(x*x)
        new_h1 = s.Poly(-changed[1, 2]-(x*z if kind == 'M1' else 0), *Y).coeff_monomial(x*x)
        assert s.expand(new_h0-2*a*h[2]-b*(h[3]+eps)-h[0]) == 0
        assert s.expand(new_h1-a*(h[3]-eps)-2*b*h[4]-h[1]) == 0
        assert (changed-base.subs({h[0]: new_h0, h[1]: new_h1}, simultaneous=True)).applyfunc(s.expand) == s.zeros(4)
        point = dict.fromkeys(h, 0)
        if kind == 'M0':
            point[h[3]] = 1
            basis = [elementary(2, 0), elementary(2, 1),
                     elementary(3, 1)-elementary(2, 4), s.diag(0, 0, 0, 1, -1)/2]
            expected_rank = 20
        else:
            basis = [elementary(2, 0), elementary(2, 1),
                     s.diag(0, 0, 0, 1, -1)/2, elementary(3, 4), elementary(4, 3)]
            expected_rank = 19
        action = global_action(P.subs(point))
        rank = witness(action)
        assert rank['rank'] == expected_rank
        for A in basis+[s.eye(5)]:
            assert action*A.reshape(25, 1) == s.zeros(action.rows, 1)
        B = s.Matrix.hstack(*[A.reshape(25, 1) for A in basis])
        assert B.rank() == len(basis)
        brackets = []
        for i in range(len(basis)):
            for j in range(i):
                bracket = basis[i]*basis[j]-basis[j]*basis[i]
                coeff = B.gauss_jordan_solve(bracket.reshape(25, 1))[0]
                if any(coeff):
                    brackets.append(dict(i=i, j=j, coefficients=list(map(str, coeff))))
        assert brackets
        if kind == 'M1':
            assert all(basis[i]*basis[j]-basis[j]*basis[i] == s.zeros(5) for i in range(2) for j in range(5))
            assert basis[3]*basis[4]-basis[4]*basis[3] == 2*basis[2]
        records[kind] = dict(strict_action=encode(M), distinct_root_minor=str(972*lam**4),
                             double_root_minor=24, triple_pivot=-162, triple_residual=encode(residual),
                             rank_one_H2_minor=8, semisimple_H2_minor=str(4*k*k*(k-1)*(k+1)),
                             resonant_M0_minor=str(resonance_minor),
                             H0_shear_coefficients=list(map(str, [new_h0, new_h1])),
                             exceptional_strict_action_rank=rank,
                             exceptional_strict_basis=[encode(A) for A in basis], brackets=brackets)

    # Position coordinates (x,w), momenta (y,-z): their cotangent lift is symplectic.
    t, alpha, f0, f1, f2 = s.symbols('t alpha f0 f1 f2')
    R, Z = s.symbols('R Z')
    F = f0*x*x+f1*x*w+f2*w*w
    deformations = {}
    for kind, A, expected in (
            ('M0', s.Matrix([[s.I, 0], [alpha*t, t]]), s.Matrix([[x*x, 0], [-alpha*x*x, 0]])),
            ('M1', s.Matrix([[s.I, 0], [s.I-alpha*t/2, t]]),
             s.Matrix([[x*x, 0], [2*x*w-alpha*x*x, x*x]]))):
        q = A*s.Matrix([x, w])
        mixed = (A.inv()*s.diag(-q[0]**2, -q[1]**2)*A).applyfunc(s.expand)
        assert mixed.subs(t, 0) == expected
        iq = A.inv()*s.Matrix([R, Z])
        g = F.subs({x: iq[0], w: iq[1]}, simultaneous=True)/A.det()
        assert s.cancel(A.det()*g.subs({R: q[0], Z: q[1]}, simultaneous=True)-F) == 0
        deformations[kind] = dict(position_change=encode(A), mixed_block=encode(mixed),
                                  limiting_mixed_block=encode(expected), vertical_limit=str(F),
                                  source_vertical=str(g))
    paths = [Path(__file__), ROOT/'code/s61_model.py', ROOT/'research/notes/TORIC_ZERO_JET_FIBER_20261101.md',
             ROOT/'research/notes/TORIC_STRICT_STABILIZER_20261101.md']
    report = dict(id='TORIC_ZERO_JET_MINIMAL_20261101', all_passed=True, field='QQ and QQ(i)',
                  minimal_families=records, universal_symplectic_shear=encode(shear), AB_deformations=deformations,
                  conclusion='The reduced actual T01 zero-first-jet fiber intersects each historical minimal-trace family M0,M1 only inside closed AB.',
                  limitations=['Reduced-set statement, not scheme intersection.', 'The third minimal-trace component and remaining nil22/zero-modular loci are not exhausted here.'],
                  proof_note='research/notes/TORIC_ZERO_JET_MINIMAL_20261101.md',
                  source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps(dict(all_passed=True, M0_M1_intersections_contained_in_AB=True)))


if __name__ == '__main__':
    main()

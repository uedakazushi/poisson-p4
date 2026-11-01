#!/usr/bin/env python3
"""Exact strict-action and degeneration certificates for F0 intersect R0 subset AB."""
import argparse
import hashlib
import json
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_model import X, quotient_coefficients
from weight1_structure import sp_basis
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
    monomials = [x*x*z, x*x*w, x*z*z, x*z*w, x*w*w,
                 z**3, z*z*w, z*w*w, w**3]
    G = sum(a*b for a, b in zip(h, monomials))
    Q = s.zeros(4)
    Q[1, 2], Q[1, 3] = -s.diff(G, w), s.diff(G, z)
    Q -= Q.T
    quadratic = [Y[i]*Y[j] for i, j in combinations_with_replacement(range(4), 2)]

    def vector(matrix):
        return s.Matrix([s.Poly(s.expand(matrix[i, j]), *Y).coeff_monomial(m)
                         for i, j in combinations(range(4), 2) for m in quadratic])

    embedding = vector(Q).jacobian(h)
    rows = list(embedding.T.rref()[1])
    inverse = embedding[rows, :].inv()
    symplectic = sp_basis()
    condition = s.Matrix.hstack(*(a[:, 1][[0, 2, 3], :] for a in symplectic))
    parabolic = [sum((v[i]*symplectic[i] for i in range(10)), s.zeros(4))
                 for v in condition.nullspace()]
    assert len(parabolic) == 7
    cols = []
    for a in parabolic:
        d = vector(a*Q+Q*a.T-sum(((a*Y)[k]*Q.diff(Y[k]) for k in range(4)), s.zeros(4)))
        cols.append(inverse*d[rows, :])
        assert embedding*cols[-1] == d
    raw = s.Matrix.hstack(*cols)
    independent = list(s.Matrix.vstack(*[raw.diff(a) for a in h]).rref()[1])
    M = raw[:, independent]
    assert M.shape == (9, 6)
    minors = []
    for name, tail, rr, cc, expected in (
            ('three_distinct_binary_roots', [1, 0, 0, 1], [2, 4, 5, 6, 7, 8], list(range(6)), 486),
            ('double_binary_root', [0, 1, 0, 0], [2, 3, 5, 6, 7], [0, 1, 2, 3, 5], 4)):
        minor = M.subs(dict(zip(h[5:], tail)))[rr, cc]
        assert s.expand(minor.det()) == expected
        minors.append(dict(case=name, rows=rr, columns=cc, constant_determinant=expected))
    T = M.subs(dict(zip(h[5:], [1, 0, 0, 0])))
    rr, cc = [2, 5, 6], [0, 2, 5]
    assert T[rr, cc].det() == 9
    rest = [1, 3, 4]
    N = (T[:, rest]-T[:, cc]*T[rr, cc].inv()*T[rr, rest]).applyfunc(s.factor)
    # Rank M<=4 implies rank N<=1. These successive minors force h4=h3=h1=0.
    obstruction = [s.factor(N[[1, 4], [0, 2]].det()),
                   s.factor(N[[1, 3], [1, 2]].subs(h[4], 0).det()),
                   s.factor(N[[0, 1], [1, 2]].subs({h[4]: 0, h[3]: 0}).det())]
    assert obstruction == [16*h[4]**2, -2*h[3]**3, 10*h[1]**2]
    zero_H_rank2 = M.subs(dict(zip(h[2:], [0, 1, 0, 0, 0, 0, 0])))
    assert s.expand(zero_H_rank2[[0, 1, 2, 3, 4], [0, 1, 2, 3, 5]].det()) == 2

    # Full projective strict stabilizer of the entire binary-double-root orbit.
    double = s.zeros(5)
    double[1, 2] = double[3, 4] = x0*x0
    double[2, 3], double[2, 4] = -z*z, 2*z*w
    double -= double.T
    action_cols = []
    for i in range(5):
        for j in range(5):
            a = s.zeros(5)
            a[i, j] = 1
            action_cols.append(s.Matrix(quotient_coefficients(a*double+double*a.T-X[j]*double.diff(X[i]))))
    action = s.Matrix.hstack(*action_cols)
    rank = witness(action)
    assert rank['rank'] == 20
    basis = []
    for i, j in ((1, 0), (2, 0), (2, 1)):
        a = s.zeros(5)
        a[i, j] = 1
        basis.append(a)
    basis.append(s.diag(0, -1, 1, 1, -1)/2)
    for a in basis+[s.eye(5)]:
        assert action*a.reshape(25, 1) == s.zeros(action.rows, 1)
    assert s.Matrix.hstack(*[a.reshape(25, 1) for a in basis+[s.eye(5)]]).rank() == 5
    bracket = basis[2]*basis[0]-basis[0]*basis[2]
    assert bracket == basis[1] and bracket != s.zeros(5)

    # Every extra-translation family with d_x=0 is a limit of the AB normal form.
    t, d, f0, f1, f2 = s.symbols('t d f0 f1 f2', nonzero=True)
    A = s.Matrix([[t, 0], [2*t, d*t**3]])
    q = s.Matrix([x, z])
    source_q = A*q
    mixed = (A.inv()*s.diag(-source_q[0]**2, -source_q[1]**2)*A).applyfunc(s.expand)
    expected = s.Matrix([[-t*t*x*x, 0],
                         [-6*x*x/d-8*t*t*x*z-2*d*t**4*z*z, -(2*t*x+d*t**3*z)**2]])
    assert (mixed-expected).applyfunc(s.expand) == s.zeros(2)
    F = f0*x*x+f1*x*z+f2*z*z
    R, Z = s.symbols('R Z')
    inverse_q = A.inv()*s.Matrix([R, Z])
    source_vertical = F.subs({x: inverse_q[0], z: inverse_q[1]}, simultaneous=True)/A.det()
    vertical = s.cancel(A.det()*source_vertical.subs({R: source_q[0], Z: source_q[1]}, simultaneous=True))
    assert vertical == F
    # The cotangent lift preserves the standard Poisson matrix exactly.
    S = s.zeros(4)
    for i, ii in enumerate((0, 2)):
        for j, jj in enumerate((0, 2)):
            S[ii, jj] = A[i, j]
    AiT = A.inv().T
    for i, ii in enumerate((1, 3)):
        for j, jj in enumerate((1, 3)):
            S[ii, jj] = AiT[i, j]
    J = s.diag(s.Matrix([[0, 1], [-1, 0]]), s.Matrix([[0, 1], [-1, 0]]))
    assert (S*J*S.T-J).applyfunc(s.cancel) == s.zeros(4)
    limit = mixed.subs(t, 0)
    assert limit == s.Matrix([[0, 0], [-6*x*x/d, 0]])

    report = dict(id='TORIC_ZERO_JET_R0_20261101', all_passed=True, field='QQ',
                  strict_R0_action=encode(M), binary_root_constant_minors=minors,
                  triple_root_schur=dict(pivot_rows=rr, pivot_columns=cc, pivot_determinant=9,
                                        residual=encode(N), successive_rank_one_obstructions=list(map(str, obstruction))),
                  zero_binary_cubic_rank_two_quadratic_constant_minor=2,
                  binary_double_strict_action_rank=rank,
                  binary_double_strict_stabilizer_basis=[encode(a) for a in basis],
                  nonzero_bracket=encode(bracket),
                  AB_deformation=dict(position_change=encode(A), symplectic_cotangent_lift=encode(S),
                                      mixed_block=encode(mixed), source_vertical=str(source_vertical),
                                      limiting_Q_yz=str(6*x*x/d), limiting_Q_yw=str(F)),
                  conclusion='The reduced actual T01 zero-first-jet fiber intersects the historical R0 locus only inside the closed AB component.',
                  limitations=['Set-theoretic reduced-locus statement, not a scheme intersection computation.',
                               'Does not classify the minimal-trace or remaining nil22 strata.',
                               'Uses the existing closed AB inclusion and the four-dimensional abelian strict-stabilizer necessity; no wedge-image limit assertion.'],
                  proof_note='research/notes/TORIC_ZERO_JET_R0_20261101.md',
                  source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in (Path(__file__), ROOT/'code/s61_model.py', ROOT/'code/weight1_structure.py',
                                           ROOT/'research/notes/TORIC_ZERO_JET_FIBER_20261101.md',
                                           ROOT/'research/notes/TORIC_STRICT_STABILIZER_20261101.md')})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps(dict(all_passed=True, R0_intersection_contained_in_AB=True)))


if __name__ == '__main__':
    main()

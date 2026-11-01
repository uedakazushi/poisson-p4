#!/usr/bin/env python3
"""Exceptional P1 normal bundle and complete coordinate differentials at x^2 dxy."""
import argparse
import hashlib
import json
from itertools import combinations, combinations_with_replacement
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
    x, y, z, w, t = s.symbols('x y z w t')
    Y, q = s.Matrix([x, y, z, w]), s.Matrix([x, z])
    ss, hh = s.symbols('s0:3'), s.symbols('h0:5')
    pars = [*ss, *hh]
    S = ss[0]*x*x+ss[1]*x*z+ss[2]*z*z
    H = sum(hh[i]*x**(4-i)*z**i for i in range(5))
    J2 = s.Matrix([[0, 1], [-1, 0]])

    def Nof(S, H):
        a = [s.Poly(s.expand(S), x, z).coeff_monomial(m) for m in (x*x, x*z, z*z)]
        h = [s.Poly(s.expand(H), x, z).coeff_monomial(x**(4-i)*z**i) for i in range(5)]
        return s.Matrix([[h[2]/6, -h[1]/4, h[0]], [h[3]/2, -h[2]/3, h[1]/2],
                         [h[4], -h[3]/4, h[2]/6]])+s.Matrix([[2*a[1], -2*a[0], 0],
                         [4*a[2], 0, -4*a[0]], [0, 2*a[2], -2*a[1]]])/96

    N = Nof(S, H)
    E = N.reshape(9, 1).jacobian(pars)
    inverse = (E.T*E).inv()*E.T
    assert inverse*E == s.eye(8)
    F = q*S/4+J2*s.Matrix([s.diff(H, x), s.diff(H, z)])
    mixed = F.jacobian(q)

    def Qof(values, vertical=0):
        out = s.zeros(4)
        A = mixed.subs(dict(zip(pars, values)), simultaneous=True)
        for i, ii in enumerate((0, 2)):
            for j, jj in enumerate((1, 3)):
                out[ii, jj], out[jj, ii] = A[i, j], -A[i, j]
        out[1, 3], out[3, 1] = vertical, -vertical
        return out

    baseN = s.zeros(3)
    baseN[0, 1] = -s.Rational(1, 24)
    baseQ = Qof(inverse*baseN.reshape(9, 1))
    assert baseQ[0, 1] == x*x and sum(bool(baseQ[i, j]) for i, j in combinations(range(4), 2)) == 1
    a, b, d, e = s.symbols('a b d e')
    chartN = s.Matrix([1, a, b])*s.Matrix([[-a*d-b*e, d, e]])
    at = {a: 0, b: 0, d: -s.Rational(1, 24), e: 0}
    core = [Qof(inverse*chartN.diff(v).subs(at).reshape(9, 1)) for v in (a, b, d, e)]
    core += [Qof([0]*8, g) for g in (x*x, x*z, z*z)]

    def upper_shear(B):
        A = s.eye(4)
        for i, ii in enumerate((0, 2)):
            for j, jj in enumerate((1, 3)):
                A[ii, jj] = B[i, j]
        return A

    b0, b1 = s.symbols('b0 b1')
    B = s.Matrix([[b0, b1], [b1, t]])
    Bprime = s.Matrix([[b0-b1*b1/t, b1/t], [b1/t, -1/t]])
    rotation = s.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]])
    h = ((rotation*upper_shear(Bprime)).inv()*upper_shear(B)).applyfunc(s.cancel)
    assert h == s.Matrix([[1, 0, -b1/t, 0], [0, 1, 0, 0], [0, 0, 1/t, 0], [0, b1, 1, t]])
    h0 = h.subs({b0: 0, b1: 0})
    K = h.diff(b1).subs({b0: 0, b1: 0})*h0.inv()

    def push(A, Q):
        return (A*Q.subs(dict(zip(Y, A.inv()*Y)), simultaneous=True)*A.T).applyfunc(s.expand)

    def core_coordinates(Q):
        A = Q[[0, 2], [1, 3]]
        F = A*q/3
        S = s.diff(F[0], x)+s.diff(F[1], z)
        R = F-q*S/4
        H = -s.Poly(R[1], x, z).coeff_monomial(x**3)*x**4/4
        for i in range(1, 5):
            H += s.Poly(R[0], x, z).coeff_monomial(x**(4-i)*z**(i-1))*x**(4-i)*z**i/i
        N = Nof(S, H)
        return s.Matrix([-24*N[1, 1], -24*N[2, 1], N[0, 1], N[0, 2],
                         *[s.Poly(Q[1, 3], x, z).coeff_monomial(m) for m in (x*x, x*z, z*z)]])

    columns = []
    for k in range(9):
        out = s.zeros(9, 1)
        if k == 0:
            out[0] = 1
        elif k == 1:
            out[1] = 1/t
            derivative = K*baseQ+baseQ*K.T-sum(((K*Y)[i]*baseQ.diff(Y[i]) for i in range(4)), s.zeros(4))
            out[2:, 0] = core_coordinates(derivative)
        else:
            out[2:, 0] = core_coordinates(push(h0, core[k-2]))
        columns.append(out)
    transition = s.Matrix.hstack(*columns).applyfunc(s.factor)
    expected = s.diag(1, 1/t, t, t*t, 1, 1/t, t, t*t, t**3)
    expected[2, 1], expected[6, 5] = 2, 12
    assert transition == expected and s.factor(transition.det()) == t**7
    gauges = []
    for k in (2, 12):
        block = s.Matrix([[1/t, 0], [k, t]])
        on_zero = s.Matrix([[1, t], [0, -k]])
        on_infinity = s.Matrix([[1/t, 1], [k, 0]])
        assert block*on_zero == on_infinity
        assert on_zero.det() == on_infinity.det() == -k
        gauges.append(dict(k=k, zero_chart=encode(on_zero), infinity_chart=encode(on_infinity)))

    pairs = list(combinations(range(4), 2))
    monomials = [Y[i]*Y[j] for i, j in combinations_with_replacement(range(4), 2)]

    def vector(Q):
        return s.Matrix([s.Poly(s.expand(Q[i, j]), *Y).coeff_monomial(m) for i, j in pairs for m in monomials])

    SB = upper_shear(s.Matrix([[0, 0], [0, t]]))
    differentials = []
    for i, j in ((0, 0), (0, 1)):
        A = s.zeros(4)
        A[(0, 2)[i], (1, 3)[j]] = A[(0, 2)[j], (1, 3)[i]] = 1
        differentials.append(vector(A*baseQ+baseQ*A.T-sum(((A*Y)[k]*baseQ.diff(Y[k]) for k in range(4)), s.zeros(4))))
    differentials += [vector(push(SB, C)) for C in core]
    D = s.Matrix.hstack(*differentials)
    span_columns, labels = [], []
    for j in range(D.cols):
        for degree in range(5):
            col = D[:, j].applyfunc(lambda entry: s.Poly(entry, t).coeff_monomial(t**degree))
            if any(col):
                span_columns.append(col)
                labels.append([j, degree])
    span = s.Matrix.hstack(*span_columns)
    rank = witness(span)
    assert rank['rank'] == 16
    source = [Path(__file__), ROOT/'research/notes/AB_INCIDENCE_20261101.md']
    report = dict(id='AB_ONE_VARIABLE_LOCAL_20261101', all_passed=True, field='QQ(t)',
                  base_Q='x^2 partial_x wedge partial_y', rank_one_core_chart=encode(chartN),
                  Lagrangian_chart_transition=encode(Bprime), coordinate_comparison=encode(h),
                  normal_coordinate_order=['b0', 'b1', 'a', 'b', 'd+1/24', 'e', 'g0', 'g1', 'g2'],
                  exceptional_normal_transition=encode(transition), trivializing_two_by_two_gauges=gauges,
                  normal_bundle='O^6 plus O(-2)^2 plus O(-3) on P1',
                  full_coordinate_differential=encode(D), coefficient_span_labels=labels,
                  coordinate_differential_span_witness=rank, conormal_global_sections_dimension=16,
                  conclusion='AB is normal along the nonzero-trace one-variable orbit; local dimension10, embedding dimension16, multiplicity7; tangent cone A6 times the affine rational normal scroll S(2,2,3) cone.',
                  limitations=['Uses the completion-generation lemma proved in the note, not the false reduced-fibers normality shortcut.',
                               'No claim that the completed local ring is isomorphic to its tangent cone.',
                               'Global AB normality and the zero-trace intersection remain unresolved.'],
                  proof_note='research/notes/AB_ONE_VARIABLE_LOCAL_20261101.md',
                  source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps(dict(all_passed=True, coordinate_rank=16, conormal_h0=16)))


if __name__ == '__main__':
    main()

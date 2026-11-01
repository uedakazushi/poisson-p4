#!/usr/bin/env python3
"""Small automorphism-action matrices proving the sharp toric zero-jet dimension.

The dimension argument uses binary-cubic strata, not a Groebner calculation.
"""
import argparse
import hashlib
import json
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_model import X, quotient_coefficients
from weight1_structure import sp_basis

ROOT = Path(__file__).resolve().parents[1]


def encode(matrix):
    return [[str(x) for x in row] for row in matrix.tolist()]


def witness(matrix):
    cols = list(matrix.rref()[1])
    rows = list(matrix[:, cols].T.rref()[1])
    minor = matrix.extract(rows, cols)
    return dict(rank=len(cols), rows=rows, columns=cols,
                determinant=str(minor.det()), matrix=encode(minor))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    h = s.symbols('h0:9')
    x0, x, y, z, w = X
    spatial = s.Matrix([x, y, z, w])
    monomials = [x*x*z, x*x*w, x*z*z, x*z*w, x*w*w,
                 z**3, z*z*w, z*w*w, w**3]
    G = sum(a*b for a, b in zip(h, monomials))
    Q = s.zeros(4)
    Q[1, 2], Q[1, 3] = -s.diff(G, w), s.diff(G, z)
    Q -= Q.T
    quadratic_monomials = [spatial[i]*spatial[j]
                           for i, j in combinations_with_replacement(range(4), 2)]
    def vector4(matrix):
        return s.Matrix([s.Poly(s.expand(matrix[i, j]), *spatial).coeff_monomial(mon)
                         for i, j in combinations(range(4), 2) for mon in quadratic_monomials])
    embedding = vector4(Q).jacobian(h)
    embedding_rows = list(embedding.T.rref()[1])
    inverse = embedding[embedding_rows, :].inv()
    symplectic = sp_basis()
    factor_constraint = s.Matrix.hstack(*(a[:, 1][[0, 2, 3], :] for a in symplectic))
    parabolic = [sum((v[i]*symplectic[i] for i in range(10)), s.zeros(4))
                 for v in factor_constraint.nullspace()]
    assert len(parabolic) == 7
    columns = []
    for a in parabolic:
        derivative = vector4(a*Q+Q*a.T-sum(((a*spatial)[k]*Q.diff(spatial[k])
                                           for k in range(4)), s.zeros(4)))
        column = inverse*derivative[embedding_rows, :]
        assert embedding*column == derivative
        columns.append(column)
    columns.append(s.Matrix(h))
    raw = s.Matrix.hstack(*columns)
    independent = list(s.Matrix.vstack(*[raw.diff(a) for a in h]).rref()[1])
    r0 = raw[:, independent]
    assert r0.shape == (9, 7)
    point = dict.fromkeys(h, 0)
    point[h[0]], point[h[6]] = 1, 1  # x^2 z+z^2 w: double binary root
    r0_witness = witness(r0.subs(point))
    assert r0_witness['rank'] == 6

    pairs5 = list(combinations(range(5), 2))
    mon5 = [X[i]*X[j] for i, j in combinations_with_replacement(range(5), 2)]
    def vector5(matrix):
        return s.Matrix([s.Poly(s.expand(matrix[i, j]), *X).coeff_monomial(mon)
                         for i, j in pairs5 for mon in mon5])
    minimal = {}
    for kind in ('M0', 'M1'):
        P = s.zeros(5)
        P[1, 2] = x0*x0+x*x
        P[3, 4] = x0*x0+(x*x if kind == 'M1' else 0)
        P[2, 3] = -s.diff(G, w)-(x*z if kind == 'M1' else 0)
        P[2, 4] = s.diff(G, z)-(x*w if kind == 'M1' else 0)
        P -= P.T
        actions = []
        # Preserving each of the three divisor hyperplanes fixes x0 and x
        # after subtracting a scalar matrix. Only the last three rows remain.
        for i in range(2, 5):
            for j in range(5):
                a = s.zeros(5)
                a[i, j] = 1
                actions.append(a*P+P*a.T-X[j]*P.diff(X[i]))
        actions.append(P)
        all_actions = s.Matrix.hstack(*[vector5(a) for a in actions])
        constant = s.Matrix.hstack(*[s.Matrix([a[i, j].coeff(x0, 2) for i, j in pairs5])
                                    for a in actions])
        kernel = s.Matrix.hstack(*constant.nullspace())
        assert kernel.cols == 10
        reduced = all_actions*kernel
        linear_rows = [i for i, (pair, mon) in enumerate((p, m) for p in pairs5 for m in mon5)
                       if s.degree(mon, x0) == 1]
        quadratic_rows = [i for i, (pair, mon) in enumerate((p, m) for p in pairs5 for m in mon5)
                          if s.degree(mon, x0) == 0]
        linear = reduced[linear_rows, :]
        linear_stack = s.Matrix.vstack(linear.subs(dict.fromkeys(h, 0)), *[linear.diff(a) for a in h])
        no_translation = s.Matrix.hstack(*linear_stack.nullspace())
        assert no_translation.cols == 8
        quad = reduced[quadratic_rows, :]*no_translation
        stack = s.Matrix.vstack(quad.subs(dict.fromkeys(h, 0)), *[quad.diff(a) for a in h])
        quad = quad[:, list(stack.rref()[1])]
        quad = quad[[i for i in range(quad.rows) if any(quad[i, :])], :]
        assert quad.cols == 6
        triple_point = dict.fromkeys(h, 0)
        triple_point[h[3]], triple_point[h[5]] = 2, 1  # 2xzw+z^3
        rank = witness(quad.subs(triple_point))
        assert rank['rank'] == 5
        minimal[kind] = dict(matrix=encode(quad), triple_root_rank_witness=rank,
                             constant_solution_dimension=10, no_translation_dimension=8)

    # The universal lower rank estimates use GL2 orbits of binary cubics and
    # the two independent first derivatives (translation-by-x directions).
    binary_checks = []
    for name, cubic, orbit_rank, derivative_rank in (
            ('three_distinct', z**3+w**3, 4, 2), ('double', z*z*w, 3, 2),
            ('triple', z**3, 2, 1)):
        def coeff(poly, degree):
            return s.Matrix([s.Poly(poly, z, w).coeff_monomial(z**i*w**(degree-i))
                             for i in range(degree+1)])
        action = s.Matrix.hstack(*(coeff(var*s.diff(cubic, deriv), 3)
                                  for var in (z, w) for deriv in (z, w)))
        derivatives = s.Matrix.hstack(*(coeff(s.diff(cubic, deriv), 2) for deriv in (z, w)))
        assert action.rank() == orbit_rank and derivatives.rank() == derivative_rank
        binary_checks.append(dict(type=name, orbit_rank=orbit_rank,
                                  derivative_rank=derivative_rank))
    # A generic affine R0 example is actually outside T01, by stabilizer rank.
    sample = dict(zip(h, range(1, 10)))
    P = s.zeros(5)
    P[1, 2] = P[3, 4] = x0*x0
    P[2, 3], P[2, 4] = -s.diff(G, w).subs(sample), s.diff(G, z).subs(sample)
    P -= P.T
    global_columns = []
    for i in range(5):
        for j in range(5):
            a = s.zeros(5)
            a[i, j] = 1
            global_columns.append(s.Matrix(quotient_coefficients(a*P+P*a.T-X[j]*P.diff(X[i]))))
    global_columns.append(s.Matrix(quotient_coefficients(P)))
    global_rank = witness(s.Matrix.hstack(*global_columns))
    assert global_rank['rank'] == 23
    paths = [Path(__file__), ROOT/'code/s61_model.py', ROOT/'code/weight1_structure.py',
             ROOT/'reports/TORIC_ZERO_JET_FIBER_20261101.json']
    report = dict(id='TORIC_ZERO_JET_DIMENSION_20261101', all_passed=True,
                  field='QQ[h0,...,h8]', potential=str(G),
                  r0_action_matrix=encode(r0), r0_double_root_rank_witness=r0_witness,
                  minimal_action_matrices=minimal, binary_cubic_checks=binary_checks,
                  generic_R0_projective_orbit_rank_witness=global_rank,
                  generic_R0_projective_orbit_dimension=22, generic_R0_stabilizer_dimension=2,
                  necessary_T01_stabilizer_dimension_at_least=4,
                  sharp_reduced_zero_jet_fiber_dimension=10,
                  proof_note='research/notes/TORIC_ZERO_JET_DIMENSION_20261101.md',
                  limitations=['No exact irreducible decomposition of the zero-first-jet fiber.',
                               'Dimension bounds use the old exhaustive compatible-quadratic classification and the geometric argument in the note.',
                               'No Groebner output or finite-field rank is used for the sharp dimension proof.',
                               'No full S61 generic exhaustion or new ambient component count.'],
                  source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps({'all_passed': True, 'zero_jet_fiber_dimension': 10}))


if __name__ == '__main__':
    main()

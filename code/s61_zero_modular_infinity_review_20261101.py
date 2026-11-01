#!/usr/bin/env python3
"""Independent matrix review of C0/L0 containment in R0 at generic beta."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
from pathlib import Path
import time
import sympy as s
from s61_model import S61, X, J as J5
from weight1_structure import J

ROOT = Path(__file__).resolve().parents[1]


def matrix(encoded):
    return s.Matrix([[s.sympify(v) for v in row] for row in encoded])


def verify():
    started = time.monotonic(); model = S61()
    one = [i for i, w in enumerate(model.reduced['parameter_weights']) if w == 1]
    two = [i for i, w in enumerate(model.reduced['parameter_weights']) if w == 2]
    one_index = {v: i for i, v in enumerate(one)}
    def beta_matrix(point):
        B = s.zeros(134, 16)
        for row, equation in enumerate(model.reduced['equations']):
            for mon, n, d in equation:
                qmon = [i for i in mon if i in two]
                if len(qmon) != 2: continue
                bmon = [i for i in mon if i in one]
                assert len(bmon) <= 1
                if bmon:
                    B[row, one_index[bmon[0]]] += s.Rational(n, d)*s.prod(point[i] for i in qmon)
        return B.applyfunc(s.factor)
    z = s.Matrix(X[1:]); cpath = ROOT/'reports/S61_C0_INFINITY_20261101.json'
    cdata = json.loads(cpath.read_text()); cchecks = []
    # Hand-constructed basis sp4 = J * symmetric matrices.
    symplectic_basis = []
    for i in range(4):
        for j in range(i, 4):
            H = s.zeros(4); H[i, j] = H[j, i] = 1
            symplectic_basis.append(J*H)
    for case in cdata['cases']:
        A, B = matrix(case['A']), matrix(case['B'])
        assert A*B == B*A and A.T*J+J*A == s.zeros(4) and B.T*J+J*B == s.zeros(4)
        av, bv = s.Matrix([0, *list(A*z)]), s.Matrix([0, *list(B*z)])
        Q = av*bv.T-bv*av.T
        point = model.normalize(X[0]**2*J5+Q)['parameters']
        assert [str(v) for v in point] == case['S61_coordinates']
        M = beta_matrix(point); assert M == matrix(case['beta_leading_matrix']) and M.rank() == 12
        v = s.Matrix(s.symbols('v0:4')); expressions = []
        for i, j, k in combinations(range(4), 3):
            expressions += s.Poly(s.expand(v[i]*Q[j+1, k+1]-v[j]*Q[i+1, k+1]+v[k]*Q[i+1, j+1]), *z).coeffs()
        assert s.Matrix(expressions).jacobian(v).rank() == 4
        witness = case['beta_rank_minor']; determinant = M.extract(witness['rows_zero_based'], witness['columns_zero_based']).det()
        assert determinant == s.Rational(witness['determinant']) and determinant != 0
        cchecks.append({'case': case['label'], 'beta_rank': 12, 'minor': str(determinant)})
    # The two further regular Jordan types have a rank-one element in their
    # two-dimensional centralizer. Verify convenient representatives directly.
    S = s.diag(1, -1, 0, 0); N1 = s.zeros(4); N1[3, 2] = 1
    N4 = s.zeros(4); N4[2, 0] = 1; N4[3, 2] = 1; N4[1, 3] = -1
    for A, rank_one in [(S+N1, N1), (N4, N4**3)]:
        assert A.T*J+J*A == s.zeros(4) and rank_one.rank() == 1 and A*rank_one == rank_one*A
        comm = s.Matrix.hstack(*[(A*C-C*A).reshape(16, 1) for C in symplectic_basis])
        assert comm.rank() == 8

    lpath = ROOT/'reports/S61_L0_INFINITY_20261101.json'
    ldata = json.loads(lpath.read_text()); q = s.symbols('q0:45'); p = [s.S.Zero]*61
    for i, v in zip(two, q, strict=True): p[i] = v
    bivector = model.reconstruct(p).subs(X[0], 1); Q = bivector[1:, 1:]-J
    assert Q == matrix(ldata['compatible_quadratic_matrix'])
    t = s.Symbol('t'); f = s.symbols('f0:3'); lchecks = []
    for case in ldata['cases']:
        U = s.Matrix(list(map(s.sympify, case['U'])))
        equations = s.Poly(Q[0, 1]+Q[2, 3], *z).coeffs()
        for i, j, k in combinations(range(4), 3):
            equations += s.Poly(s.expand(U[i]*Q[j, k]-U[j]*Q[i, k]+U[k]*Q[i, j]), *z).coeffs()
        M = s.Matrix(equations).jacobian(q); assert M == matrix(case['constraints'])
        N = matrix(case['kernel']); assert N.shape == (45, 4)
        assert (M*N).applyfunc(s.cancel) == s.zeros(M.rows, 4)
        assert N[case['free_q_columns'], :] == s.eye(4)
        determinant = s.factor(M.extract(case['constraint_minor_rows'], case['constraint_minor_columns']).det())
        assert determinant == s.sympify(case['constraint_minor_determinant'])
        if case['label'] == 'off_diagonal':
            assert determinant == t**13
            for v in N:
                denominator = s.Poly(s.denom(s.cancel(v)), t)
                assert len(denominator.terms()) == 1  # No hidden nonzero-t poles.
        else: assert determinant != 0 and not determinant.free_symbols
        c3 = s.Symbol('c3'); qp = Q.subs(dict(zip(q, N*s.Matrix([*f, c3]), strict=True)), simultaneous=True).applyfunc(s.factor)
        assert (qp-matrix(case['parametric_Q'])).applyfunc(s.cancel) == s.zeros(4)
        assert all(qp[i, j].subs(c3, 0) == 0 for i, j in combinations(range(4), 2) if (i, j) != (1, 3))
        point = [s.S.Zero]*61
        for i, v in zip(two, N*s.Matrix([*f, 1]), strict=True): point[i] = v
        B = beta_matrix(point); assert (B-matrix(case['beta_leading_matrix'])).applyfunc(s.cancel) == s.zeros(134, 16)
        if case['label'] == 'off_diagonal':
            used = set(); roots = set(); rank_sum = 0
            for block in case['triangular_blocks']:
                rr, cc, rank = block['rows'], block['columns'], block['generic_rank']
                assert not (used & set(cc))
                assert all(B[i, j] == 0 for i in rr for j in range(16) if j not in used | set(cc))
                D = B.extract(rr, cc); assert not any(v.has(*f) for v in D)
                gcd = s.S.Zero
                for rs in combinations(range(D.rows), rank):
                    for cs in combinations(range(D.cols), rank):
                        numerator, denominator = s.fraction(s.cancel(D.extract(rs, cs).det()))
                        assert len(s.Poly(denominator, t).terms()) == 1
                        gcd = s.gcd(gcd, numerator)
                while s.Poly(gcd, t).eval(0) == 0: gcd = s.cancel(gcd/t)
                gcd = s.Poly(gcd, t).monic().as_expr()
                assert gcd == s.sympify(block['critical_polynomial_away_from_t_zero'])
                current = set(s.solve(gcd, t)); assert not (roots & current)
                for a in current: assert D.subs(t, a).rank() >= rank-1
                roots |= current; used |= set(cc); rank_sum += rank
            assert rank_sum == 12 and len(roots) == 10
            lchecks.append({'case': case['label'], 'uniform_rank_lower_bound': 11,
                'distinct_exception_values': sorted(map(str, roots)), 'no_hidden_denominator_zero': True})
        else:
            value = s.factor(B.extract(case['beta_minor_rows'], case['beta_minor_columns']).det())
            assert value == s.Rational(4886335125, 64)
            lchecks.append({'case': case['label'], 'uniform_rank_lower_bound': 12, 'minor': str(value)})
    paths = [Path(__file__), cpath, lpath, ROOT/'code/s61_model.py', ROOT/'code/weight1_structure.py',
             model.source/'symplectic_slice_reduced_equations.json']
    return {'id': 'S61_ZERO_MODULAR_INFINITY_REVIEW_20261101', 'all_passed': True,
        'generated_utc': datetime.now(timezone.utc).isoformat(), 'C0_checks': cchecks, 'L0_checks': lchecks,
        'other_regular_jordan_types_have_rank_one_centralizer_element': True,
        'incidence_dimension_upper_bounds': {'C0_nonR0': 14, 'L0_nonR0': 14},
        'proof_note': 'research/notes/S61_ZERO_MODULAR_INFINITY_REVIEW_20261101.md',
        'scope': 'Generic-beta reduced zero-modular unsaturated infinity is contained in R0, using the historical R0/C0/L0 exhaustion; no nonzero-modular or affine-fiber exhaustion.',
        'elapsed_seconds': round(time.monotonic()-started, 3),
        'input_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); report = verify(); args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out: json.dump(report, out, indent=2); out.write('\n')
    print(json.dumps({'all_passed': True, 'seconds': report['elapsed_seconds']}))


if __name__ == '__main__': main()

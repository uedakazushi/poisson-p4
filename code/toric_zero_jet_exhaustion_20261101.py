#!/usr/bin/env python3
"""Uniform strict-stabilizer exclusion of zero-modular C0/L0 outside R0."""
import argparse
import hashlib
import json
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_model import X
from toric_zero_jet_dimension_20261101 import encode, witness

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    x, y, z, w = X[1:]
    Y = s.Matrix([x, y, z, w])
    t, c3 = s.symbols('t c3')
    f = s.symbols('f0:3')
    pairs = list(combinations(range(4), 2))
    mons2 = [Y[i]*Y[j] for i, j in combinations_with_replacement(range(4), 2)]

    def vector(Q, degree=2):
        mons = mons2 if degree == 2 else list(Y)
        return s.Matrix([s.Poly(s.expand(Q[i, j]), *Y).coeff_monomial(m)
                         for i, j in pairs for m in mons])

    # Order (q,p)=(x,z;y,w). sp4 consists of (a,b;c,-a^t), b,c symmetric.
    order = [0, 2, 1, 3]
    basis = []
    for typ in ('upper_symmetric', 'gl2', 'lower_symmetric'):
        ij = list(combinations_with_replacement(range(2), 2)) if typ != 'gl2' else [(0, 0), (0, 1), (1, 0), (1, 1)]
        for i, j in ij:
            A = s.zeros(4)
            if typ == 'gl2':
                A[i, j], A[j+2, i+2] = 1, -1
            else:
                row, col = (0, 2) if typ == 'upper_symmetric' else (2, 0)
                A[row+i, col+j] = A[row+j, col+i] = 1
            B = s.zeros(4)
            for r in range(4):
                for c in range(4):
                    B[order[r], order[c]] = A[r, c]
            basis.append(B)
    J = s.diag(s.Matrix([[0, 1], [-1, 0]]), s.Matrix([[0, 1], [-1, 0]]))
    assert all(A*J+J*A.T == s.zeros(4) for A in basis)
    assert s.Matrix.hstack(*[A.reshape(16, 1) for A in basis]).rank() == 10

    def action(Q):
        return s.Matrix.hstack(*[vector(A*Q+Q*A.T-sum(((A*Y)[k]*Q.diff(Y[k])
                                                     for k in range(4)), s.zeros(4))) for A in basis])

    source_path = ROOT/'reports/S61_L0_INFINITY_20261101.json'
    source = json.loads(source_path.read_text())
    records = []
    for case in source['cases']:
        Q = s.Matrix(case['parametric_Q']).subs(c3, 1)
        assert s.expand(Q[0, 1]+Q[2, 3]) == 0
        assert s.expand(Q[0, 1]*Q[2, 3]-Q[0, 2]*Q[1, 3]+Q[0, 3]*Q[1, 2]) == 0
        M = action(Q)
        blocks = []
        for start, end in ((0, 3), (3, 7), (7, 10)):
            rows = [i for i in range(M.rows)
                    if all(M[i, j] == 0 for j in range(end, 10))
                    and all(not M[i, j].has(*f) for j in range(start, end))
                    and any(M[i, j] for j in range(start, end))]
            cols = list(range(start, end))
            D = M[rows, cols]
            blocks.append(dict(rows=rows, columns=cols, matrix=encode(D)))
        B, A, C = [s.Matrix(block['matrix']) for block in blocks]
        # The three nested row blocks are lower triangular in the ordered column groups.
        assert B.cols == 3 and A.cols == 4 and C.cols == 3
        if case['label'] == 'off_diagonal':
            assert B[[1, 2, 5], :].det() == -2
            assert s.factor(A[[0, 1, 2], [0, 1, 2]].det()) == (t+1)**2/t
            assert A[1, 0] == 1
            expected_C = s.diag(1-3*t, -2*(t-1)**2/t, (t-3)/t)
            assert (C-expected_C).applyfunc(s.factor) == s.zeros(3)
            # At most one diagonal entry vanishes on t!=0; all three are nonzero at -1.
            numerators = [1-3*t, (t-1)**2, t-3]
            assert all(s.degree(s.gcd(a, b), t) == 0 for a, b in combinations(numerators, 2))
            assert C.subs(t, -1).det() != 0
            orbit_bounds = dict(t_not_minus_one=8, t_minus_one=7)
            # Translation coefficients: dx,dz first, then dy,dw.
            T = s.Matrix.hstack(*[vector(Q.diff(v), 1) for v in Y])
            # Mixed coefficients force dx=dz=0; only the last two entries can drop.
            assert Q[0, 3] == -t*x*x and Q[1, 2] == -z*z/t
            assert s.expand(Q[1, 3].diff(y)-(2*t-1)*x) == 0
            assert s.expand(Q[1, 3].diff(w)-(2/t-1)*z) == 0
            assert s.expand((2*t-1)+2*(2-t)) == 3
            assert T.subs(t, -1).rank() == 4
            translation_bounds = dict(all_t_nonzero_kernel_at_most=1,
                                      t_minus_one_kernel=0, possible_drop_parameters=['1/2', '2'])
        else:
            b_witness, a_witness, c_witness = witness(B), witness(A), witness(C)
            assert (b_witness['rank'], a_witness['rank'], c_witness['rank']) == (3, 3, 3)
            assert C.det() == 128
            T = s.Matrix.hstack(*[vector(Q.diff(v), 1) for v in Y])
            sample = T.subs(dict.fromkeys(f, 0))
            cc = list(sample.rref()[1])
            rr = list(sample[:, cc].T.rref()[1])
            determinant = s.expand(T[rr, cc].det())
            assert len(cc) == 4 and determinant != 0 and not determinant.free_symbols
            orbit_bounds = dict(all_parameters=9)
            translation_bounds = dict(all_parameters_kernel=0, rows=rr, columns=cc,
                                      constant_minor=str(determinant))
            blocks[0]['rank_witness'] = b_witness
            blocks[1]['rank_witness'] = a_witness
            blocks[2]['rank_witness'] = c_witness
        records.append(dict(label=case['label'], Q=encode(Q), sp4_action=encode(M),
                            triangular_blocks=blocks, orbit_rank_lower_bounds=orbit_bounds,
                            translation_action=encode(T), translation_kernel_bounds=translation_bounds,
                            strict_stabilizer_upper_bound=3 if case['label'] == 'off_diagonal' else 1))
    H1, H2 = s.Matrix([x, -y, 0, 0]), s.Matrix([0, 0, z, -w])
    Qc = H1*H2.T-H2*H1.T
    cartan_orbit = witness(action(Qc))
    cartan_translation = witness(s.Matrix.hstack(*[vector(Qc.diff(v), 1) for v in Y]))
    assert cartan_orbit['rank'] == 8 and cartan_translation['rank'] == 4
    H, U = s.Matrix([-x, y, z, -w]), s.Matrix([0, z, 0, x])
    Qj = H*U.T-U*H.T
    assert Qj == s.Matrix(records[0]['Q']).subs({t: 1, **dict.fromkeys(f, 0)})
    paths = [Path(__file__), source_path, ROOT/'research/notes/S61_L0_INFINITY_20261101.md',
             ROOT/'research/notes/S61_C0_INFINITY_20261101.md',
             ROOT/'research/notes/S61_ZERO_MODULAR_INFINITY_REVIEW_20261101.md',
             ROOT/'research/notes/TORIC_ZERO_JET_NONZERO_TRACE_20261101.md',
             ROOT/'research/notes/TORIC_ZERO_JET_R0_20261101.md',
             ROOT/'research/notes/TORIC_ZERO_JET_FIBER_20261101.md']
    report = dict(id='TORIC_ZERO_JET_EXHAUSTION_20261101', all_passed=True, field='QQ(t,f0,f1,f2), all t!=0',
                  sp4_basis=[encode(A) for A in basis], L0_uniform_exclusions=records,
                  C0_Cartan_orbit_witness=cartan_orbit, C0_Cartan_translation_witness=cartan_translation,
                  C0_Cartan_strict_stabilizer_dimension=2, C0_Jordan_representative_is_L0_t_one=True,
                  conclusion='The reduced actual T01 zero-first-jet fiber equals the closed historical AB component, irreducible of dimension ten.',
                  limitations=['Reduced equality, not a scheme-theoretic fiber or radical ideal presentation.',
                               'Uses the historical exhaustive zero-modular decomposition R0 union C0 union L0 and the independently reviewed proper-incidence boundary parameterizations.',
                               'No full S61 generic exhaustion follows.'],
                  proof_note='research/notes/TORIC_ZERO_JET_EXHAUSTION_20261101.md',
                  source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps(dict(all_passed=True, reduced_zero_jet_fiber='closed_AB', dimension=10)))


if __name__ == '__main__':
    main()

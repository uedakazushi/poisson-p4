#!/usr/bin/env python3
"""Independent full action-matrix and combined-minor audit of F0=closed AB."""
import argparse
import hashlib
import json
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import sympy as s
from weight1_structure import J

ROOT = Path(__file__).resolve().parents[1]


def verify():
    src = ROOT/'reports/TORIC_ZERO_JET_EXHAUSTION_20261101.json'
    data = json.loads(src.read_text())
    x, y, z, w, t = s.symbols('x y z w t')
    f0, f1, f2 = s.symbols('f0 f1 f2')
    coords = (x, y, z, w)
    Y = s.Matrix(coords)
    pairs = list(combinations(range(4), 2))
    mons = [coords[i]*coords[j] for i, j in combinations_with_replacement(range(4), 2)]
    basis = [s.Matrix(A) for A in data['sp4_basis']]
    assert len(basis) == 10 and s.Matrix.hstack(*[A.reshape(16, 1) for A in basis]).rank() == 10
    assert all(A*J+J*A.T == s.zeros(4) for A in basis)
    F = f0*x*x+f1*x*z+f2*z*z
    def tensor(values):
        Q = s.zeros(4)
        for (i, j), value in zip(pairs, values):
            Q[i, j] = value; Q[j, i] = -value
        return Q
    # The two written normal forms, without importing the author's kernel routine.
    forms = [tensor([-x*z, 0, -t*x*x, -z*z/t, (2*t-1)*x*y+(2/t-1)*z*w+F, x*z]),
             tensor([-x*(x+z), 0, x*x, (x+z)**2, -3*x*y-4*x*w-3*z*w+F, x*(x+z)])]
    def vec(Q, degree=2):
        monomials = mons if degree == 2 else coords
        return s.Matrix([s.Poly(s.expand(Q[i, j]), *coords).coeff_monomial(mon)
                         for i, j in pairs for mon in monomials])
    def action(Q):
        columns = []
        for A in basis:
            field = A*Y
            values = []
            for i, j in pairs:
                values.append(s.expand(sum(A[i, k]*Q[k, j]+A[j, k]*Q[i, k]
                    -field[k]*s.diff(Q[i, j], coords[k]) for k in range(4))))
            columns.append(vec(tensor(values)))
        return s.Matrix.hstack(*columns)
    results = []
    for Q, record in zip(forms, data['L0_uniform_exclusions']):
        assert (Q-s.Matrix(record['Q']).subs(dict(zip(s.symbols('x1:5'), coords)))).applyfunc(s.cancel) == s.zeros(4)
        M = action(Q)
        assert (M-s.Matrix(record['sp4_action'])).applyfunc(s.cancel) == s.zeros(60, 10)
        blocks = record['triangular_blocks']
        for block in blocks:
            assert M[block['rows'], block['columns']] == s.Matrix(block['matrix'])
            assert all(M[i, j] == 0 for i in block['rows'] for j in range(max(block['columns'])+1, 10))
        if record['label'] == 'off_diagonal':
            rb = [blocks[0]['rows'][i] for i in [1, 2, 5]]
            ra = [blocks[1]['rows'][i] for i in [0, 1, 2]]
            combined = []
            for indices in combinations(range(3), 2):
                rc = [blocks[2]['rows'][i] for i in indices]
                columns = [0, 1, 2, 3, 4, 5]+[7+i for i in indices]
                determinant = s.factor(M[rb+ra+rc, columns].det())
                assert determinant != 0 and not determinant.has(f0, f1, f2)
                combined.append(determinant)
            gcd = s.factor(s.gcd_list([s.fraction(d)[0] for d in combined]))
            assert set(s.solve(gcd, t)) == {-1}
            # At the only remaining common zero, select the rank-one a block
            # and all three c rows, giving a single full seven-dimensional minor.
            rows = rb+[blocks[1]['rows'][1]]+blocks[2]['rows']
            columns = [0, 1, 2, 3, 7, 8, 9]
            exceptional = s.factor(M[rows, columns].subs(t, -1).det())
            assert exceptional != 0 and not exceptional.free_symbols
            # Derivatives of Qxw and Qyz kill vx,vz before the remaining two.
            assert Q[0, 3] == -t*x*x and Q[1, 2] == -z*z/t
            assert s.diff(Q[1, 3], y) == (2*t-1)*x
            assert s.diff(Q[1, 3], w) == (2/t-1)*z
            assert s.gcd(2*t-1, 2-t) == 1
            results.append(dict(label=record['label'], combined_rank8_minors=list(map(str, combined)),
                common_numerator_gcd=str(gcd), exceptional_rank7_minor=str(exceptional),
                all_t_nonzero_strict_stabilizer_dimension_upper_bound=3))
        else:
            rows = []; columns = []
            for block in blocks:
                witness = block['rank_witness']
                assert witness['rank'] == 3
                rows += [block['rows'][i] for i in witness['rows']]
                columns += [block['columns'][i] for i in witness['columns']]
            determinant = s.factor(M[rows, columns].det())
            assert determinant != 0 and not determinant.free_symbols
            translation = s.Matrix.hstack(*[vec(Q.diff(c), 1) for c in coords])
            witness = record['translation_kernel_bounds']
            trdet = s.factor(translation[witness['rows'], witness['columns']].det())
            assert trdet != 0 and not trdet.free_symbols
            results.append(dict(label=record['label'], combined_rank9_minor=str(determinant),
                translation_rank4_minor=str(trdet), strict_stabilizer_dimension_upper_bound=1))
    H1, H2 = s.Matrix([x, -y, 0, 0]), s.Matrix([0, 0, z, -w])
    cartan = H1*H2.T-H2*H1.T
    assert action(cartan).rank() == 8
    assert s.Matrix.hstack(*[vec(cartan.diff(c), 1) for c in coords]).rank() == 4
    H, U = s.Matrix([-x, y, z, -w]), s.Matrix([0, z, 0, x])
    assert H*U.T-U*H.T == forms[0].subs({t: 1, f0: 0, f1: 0, f2: 0})
    paths = [Path(__file__), src, ROOT/'code/weight1_structure.py',
        ROOT/'research/notes/S61_ZERO_MODULAR_INFINITY_REVIEW_20261101.md',
        ROOT/'research/notes/TORIC_ZERO_JET_R0_REVIEW_20261101.md',
        ROOT/'research/notes/TORIC_ZERO_JET_MINIMAL_REVIEW_20261101.md',
        ROOT/'research/notes/TORIC_ZERO_JET_NONZERO_TRACE_REVIEW_20261101.md']
    return dict(id='TORIC_ZERO_JET_EXHAUSTION_REVIEW_20261101', all_passed=True,
        independent_combined_minors=results, C0_Cartan_strict_stabilizer_dimension=2,
        accepted_conclusion='The actual reduced T01 beta=0 fiber equals closed AB and is irreducible of dimension10.',
        scope='Reduced sets, conditional on the explicit historical zero-modular exhaustion; no scheme or full S61 conclusion.',
        proof_note='research/notes/TORIC_ZERO_JET_EXHAUSTION_REVIEW_20261101.md',
        input_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); result = verify()
    with args.output.open('x') as out:
        json.dump(result, out, indent=2); out.write('\n')
    print(json.dumps(dict(all_passed=True)))


if __name__ == '__main__':
    main()

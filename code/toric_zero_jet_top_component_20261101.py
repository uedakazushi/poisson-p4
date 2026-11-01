#!/usr/bin/env python3
"""Certificates isolating AB as the unique maximal-dimensional zero-jet component."""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
import sympy as s
from s61_model import S61, X, quotient_coefficients
from toric_zero_jet_dimension_20261101 import witness, encode

ROOT = Path(__file__).resolve().parents[1]


def schur(matrix, rows, cols):
    remaining = [j for j in range(matrix.cols) if j not in cols]
    pivot = matrix.extract(rows, cols)
    result = (matrix[:, remaining]-matrix[:, cols]*pivot.inv()*matrix[rows, remaining]).applyfunc(s.factor)
    return pivot.det(), result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    source = ROOT/'reports/TORIC_ZERO_JET_DIMENSION_20261101.json'
    previous = json.loads(source.read_text())
    h = s.symbols('h0:9')
    r0 = s.Matrix(previous['r0_action_matrix'])
    double = r0.subs(dict(zip(h[5:], [0, 1, 0, 0])))
    determinant, residual = schur(double, [2, 3, 5, 6, 7], [0, 1, 2, 3, 5])
    assert determinant == 4
    residual = residual[[i for i in range(residual.rows) if any(residual[i, :])], :]
    expected = s.Matrix([[2*(h[0]-h[2]*h[3]), -2*(h[0]-h[2]*h[3])],
                         [4*h[1]-h[3]**2, -(4*h[1]-4*h[2]*h[4]-h[3]**2)/2],
                         [4*h[4], -h[4]]])
    assert (residual-expected).applyfunc(s.expand) == s.zeros(3, 2)
    triple_point = dict(zip(h, [1, 2, 3, 4, 5, 1, 0, 0, 0]))
    triple_witness = witness(r0.subs(triple_point))
    assert triple_witness['rank'] == 7
    minimal = {}
    for kind in ('M0', 'M1'):
        matrix = s.Matrix(previous['minimal_action_matrices'][kind]['matrix'])
        matrix = matrix.subs(dict(zip(h[5:], [1, 0, 0, 0])))
        determinant, reduced = schur(matrix, [3, 7, 9], [1, 0, 4])
        assert determinant != 0 and not determinant.free_symbols
        reduced = reduced[[i for i in range(reduced.rows) if any(reduced[i, :])], :]
        assert reduced.shape == (5, 3)
        # Two rows alone force h4=0 when the residual rank is at most one.
        assert s.expand(reduced.extract([1, 2], [1, 2]).det()+4*h[4]**2/3) == 0
        two_rows = reduced.subs(h[4], 0)[[0, 3], :]
        minimal[kind] = dict(pivot_determinant=str(determinant),
                             residual=encode(reduced), two_rows_after_h4_zero=encode(two_rows))

    # One nonfactorable divisor witness in each non-AB nil22 saturation.
    historical_path = ROOT/'legacy/v07/data/nil22_classification_summary.json'
    historical = json.loads(historical_path.read_text())
    xx, yy, zz, ww = s.symbols('x y z w')
    tests = []
    for name, entry in historical.items():
        if not isinstance(entry, dict) or 'Q' not in entry:
            continue
        q = {}
        for key, value in entry['Q'].items():
            pair = tuple(int(v.strip()) for v in key.strip('()').split(','))
            q[pair] = s.sympify(value)
        pf = s.factor(q[0, 1]*q[2, 3]-q[0, 2]*q[1, 3]+q[0, 3]*q[1, 2], extension=True)
        record = dict(stratum=name, pfaffian=str(pf))
        if pf.free_symbols.intersection({yy, ww}):
            factorable = False
            record['obstruction'] = 'Pfaffian depends on variables outside the rank-two trace support.'
        elif pf == 0:
            factorable = False
            record['obstruction'] = 'Pfaffian zero would force the trace to have rank at most one.'
        else:
            coefficient, factors = s.factor_list(pf, xx, zz, extension=s.sqrt(3))
            if any(exponent % 2 for _, exponent in factors):
                factorable = False
                record['obstruction'] = 'Binary quartic has a factor with odd multiplicity.'
                record['factorization'] = [[str(f), int(exponent)] for f, exponent in factors]
            else:
                square = s.sqrt(coefficient)*s.prod(f**(exponent//2) for f, exponent in factors)
                determinants = []
                for sign in (1, -1):
                    quadratic = s.Poly(s.expand(xx*zz+sign*2*square), xx, zz)
                    determinants.append(s.simplify(quadratic.coeff_monomial(xx**2)*quadratic.coeff_monomial(zz**2)
                                                   -quadratic.coeff_monomial(xx*zz)**2/4))
                factorable = all(value == 0 for value in determinants)
                record['square_root_quartic'] = str(square)
                record['two_square_determinants'] = list(map(str, determinants))
        assert factorable == (name == 'AB')
        record['passes_necessary_complete_factorization'] = factorable
        tests.append(record)

    # The zero-modular L0 component has a point with stabilizer only three.
    x0, x, y, z, w = X
    U = s.Matrix([0, z, 0, 3*x])
    V = s.Matrix([-3*x, 5*y, z, w])
    Q = V*U.T-U*V.T
    P = s.zeros(5)
    P[1:5, 1:5] = Q
    for i, j in ((1, 2), (3, 4)):
        P[i, j] += x0*x0
        P[j, i] -= x0*x0
    model = S61()
    model.verify_point(model.normalize(P))
    columns = []
    for i in range(5):
        for j in range(5):
            a = s.zeros(5)
            a[i, j] = 1
            columns.append(s.Matrix(quotient_coefficients(a*P+P*a.T-X[j]*P.diff(X[i]))))
    columns.append(s.Matrix(quotient_coefficients(P)))
    L0_witness = witness(s.Matrix.hstack(*columns))
    assert L0_witness['rank'] == 22
    paths = [Path(__file__), source, historical_path, ROOT/'code/s61_model.py',
             ROOT/'legacy/v07/minimal_thickness.tex', ROOT/'legacy/v07/proofs/zero_reduced_components_theorem.txt']
    report = dict(id='TORIC_ZERO_JET_TOP_COMPONENT_20261101', all_passed=True,
                  field='QQ[h0,...,h4] with exact algebraic factor tests',
                  R0_double_pivot_determinant='4', R0_double_residual=encode(residual),
                  R0_double_fixed_H_dimension=2, R0_triple_rank_witness=triple_witness,
                  minimal_triple_reductions=minimal, minimal_triple_fixed_H_dimension=2,
                  nil22_factorization_tests=tests, L0_orbit_witness=L0_witness,
                  L0_projective_stabilizer_dimension=3,
                  conclusion='AB is the unique ten-dimensional irreducible component of the actual reduced T01 zero-first-jet fiber.',
                  limitations=['Additional irreducible components of dimension at most nine remain unexcluded.',
                               'No exact ideal or full reduced-fiber exhaustion is asserted.',
                               'The historical exhaustive Q45 classification remains a dependency.'],
                  proof_note='research/notes/TORIC_ZERO_JET_TOP_COMPONENT_20261101.md',
                  source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps({'all_passed': True, 'unique_top_dimensional_component': 'AB'}))


if __name__ == '__main__':
    main()

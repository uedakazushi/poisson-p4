#!/usr/bin/env python3
"""Independent raw134 multiplicity-count Taylor2 and generic T01-fiber review."""
import argparse
import hashlib
import json
import time
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations_with_replacement as cwr
from math import comb
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def matrix(rows):
    return s.Matrix(rows).applyfunc(s.Rational)


def verify():
    start = time.time()
    report_path = 'reports/TORIC_ZERO_FIBER_GENERIC_REDUCED_20261102.json'
    main = json.loads((ROOT / report_path).read_text())
    assert main['all_passed']
    for path, digest in main['input_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest, path
    old = json.loads((ROOT / 'reports/TORIC_ZERO_JET_FIBER_20261101.json').read_text())
    sub = {s.Symbol(a): s.Rational(b) for a, b in old['parameter_sample'].items()}
    point = [Fraction(str(s.sympify(a).subs(sub))) for a in old['S61_parameters']]
    assert list(map(str, point)) == main['point']
    encoded = json.loads((ROOT / 'legacy/v02/data/symplectic_slice_reduced_equations.json').read_text())
    one = [i for i, w in enumerate(encoded['parameter_weights']) if w == 1]
    two = [i for i, w in enumerate(encoded['parameter_weights']) if w == 2]
    assert all(point[i] == 0 for i in one)

    # Taylor coefficients from exponent multiplicities/binomial coefficients;
    # no occurrence-index differentiation loop or main producer is used.
    def coefficient(counts, removed):
        value = Fraction(1)
        for i, power in counts.items():
            r = removed.get(i, 0)
            if r > power or (not point[i] and power > r):
                return Fraction(0)
            value *= comb(power, r) * point[i] ** (power - r)
        return value

    jac = [[Fraction(0) for _ in range(61)] for _ in range(134)]
    quadratic = []
    for row, equation in enumerate(encoded['equations']):
        h = defaultdict(Fraction)
        constant = Fraction(0)
        for mon, num, den in equation:
            if sum(not point[i] for i in mon) > 2:
                continue
            counts = Counter(mon)
            scalar = Fraction(num, den)
            constant += scalar * coefficient(counts, {})
            for i in counts:
                jac[row][i] += scalar * coefficient(counts, {i: 1})
            for i, j in cwr(sorted(counts), 2):
                h[i, j] += scalar * coefficient(counts, Counter((i, j)))
        assert constant == 0
        quadratic.append({ij: value for ij, value in h.items() if value})
    J = s.Matrix(jac)
    assert J.rank() == 43 and J[:, one].rank() == 10 and J[:, two].rank() == 33
    saved_minor = J.extract(main['Jacobian_minor_rows'], main['Jacobian_minor_columns']).det()
    assert saved_minor == s.Rational(main['Jacobian_minor_determinant'])
    U = matrix(main['tangent_basis'])
    assert U.rank() == 18 and J * U == s.zeros(134, 18)
    mons = list(cwr(range(18), 2))
    index = {mon: i for i, mon in enumerate(mons)}
    support = [[(j, U[i, j]) for j in range(18) if U[i, j]] for i in range(61)]
    A = s.zeros(134, len(mons))
    for row, h in enumerate(quadratic):
        for (i, j), value in h.items():
            for a, ca in support[i]:
                for b, cb in support[j]:
                    A[row, index[tuple(sorted((a, b)))]] += s.Rational(value) * ca * cb

    # Reverse both pivot preferences to use a different implicit-elimination
    # presentation, then compare whole quadratic row spaces.
    columns = [60 - i for i in J[:, ::-1].rref()[1]]
    rows = [133 - i for i in J[::-1, columns].T.rref()[1]]
    block = J.extract(rows, columns)
    assert len(rows) == len(columns) == 43 and block.det() != 0
    projection = J[:, columns] * block.inv()
    assert projection * J[rows, :] == J
    residual = A - projection * A[rows, :]
    stored = matrix(main['quadratic_basis'])
    assert residual.rank() == stored.rank() == residual.col_join(stored).rank() == 7

    C = matrix(main['coordinate_frame'])
    assert C.det() != 0
    x = s.Matrix(s.symbols('x0:6'))
    v = C[:, :6] * x
    restricted = []
    for row in stored.tolist():
        restricted.append(s.expand(sum(c * v[i] * v[j] for c, (i, j) in zip(row, mons))))
    assert restricted == [s.sympify(f) for f in main['active_quadratics']]
    all_hessians = []
    fullvars = s.symbols('v0:18')
    for row in stored.tolist():
        hh = s.zeros(18)
        for value, (i, j) in zip(row, mons):
            hh[i, j] += value
            hh[j, i] += value
        assert hh * C[:, 6:] == s.zeros(18, 12)
        all_hessians.append(hh)
    assert s.Matrix.vstack(*all_hessians).rank() == 6

    product = [s.sympify(f) for f in main['product_coordinates']]
    assert s.Matrix(product).jacobian(x).det() == 8
    u, v0, a, b, c, d = product
    generators = [u*u, u*v0, v0*v0, u*a, u*b, v0*c, v0*d]
    change = matrix(main['quadratic_row_change'])
    assert change.det() != 0
    assert all(s.expand(f) == 0 for f in change * s.Matrix(generators) - s.Matrix(restricted))
    radical = C.inv()[[2, 5], :]
    assert radical == matrix(main['radical_forms'])
    V = U * s.Matrix.hstack(*radical.nullspace())
    assert V.rank() == 16 and V == matrix(main['T01_tangent_basis'])
    assert V[one, :].rank() == 6
    beta_minor = V.extract(main['beta_minor_rows'], main['beta_minor_columns']).det()
    assert beta_minor == s.Rational(main['beta_minor_determinant']) == -s.Rational(3, 4)

    # Independent monomial counting for the active product ideal.
    product_mons = [(0, 0), (0, 1), (1, 1), (0, 2), (0, 3), (1, 4), (1, 5)]
    active_hf = []
    for degree in range(9):
        standard = 0
        for mon in cwr(range(6), degree):
            count = Counter(mon)
            if not any(all(count[i] >= e for i, e in Counter(g).items()) for g in product_mons):
                standard += 1
        expected = comb(degree + 3, 3) + (2 * degree if degree else 0)
        assert standard == expected
        active_hf.append(standard)
    z = s.Symbol('z')
    numerator = s.expand(1 + 2*z*(1-z)**2)
    assert numerator == 1 + 2*z - 4*z**2 + 2*z**3 and numerator.subs(z, 1) == 1

    paths = sorted(set(main['input_sha256']) | {
        report_path, str(Path(__file__).relative_to(ROOT)),
        'research/notes/TORIC_ZERO_FIBER_GENERIC_REDUCED_REVIEW_20261102.md'})
    return dict(id='TORIC_ZERO_FIBER_GENERIC_REDUCED_REVIEW_20261102', all_passed=True,
                independent_raw134_Taylor2=True, Jacobian_rank=43,
                beta_column_rank=10, q_column_rank=33, quadratic_row_space_rank=7,
                independently_selected_pivot_rows=rows, independently_selected_pivot_columns=columns,
                common_Hessian_kernel_dimension=12, active_product_coordinate_determinant=8,
                product_ideal_Hilbert_values=active_hf, lower_initial_dimension=16,
                lower_initial_multiplicity=1, T01_beta_differential_rank=6,
                beta_minor_determinant=str(beta_minor), T01_regular_dimension=16,
                zero_fiber_regular_dimension=10, zero_fiber_generic_length=1,
                multiplicity_and_fiber_proof_accepted=True,
                limitations=['Known reduced irreducible T01 and its membership at this AB point are frozen dependencies.',
                             'No identification of the full surrounding S61 local ring or full initial ideal.',
                             'No assertion that the T01 zero fiber is reduced at every AB point.'],
                input_sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
                seconds=time.time() - start)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print('TORIC ZERO FIBER GENERIC REDUCED REVIEW PASS', result['seconds'])

#!/usr/bin/env python3
"""Construct the unique degree-five equivariant kernel covariant of toric L.

All representation, covariance and kernel computations are over QQ. The two
residue samples fix a useful scale; they do not prove a universal E[2] identity.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from itertools import combinations_with_replacement, permutations, product
import json
import math
from pathlib import Path
import time

import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from toric_algebra_fiber import TRIPLES
from weight1_structure import COORDS, J, constraints, infinitesimal_action, sp_basis

ROOT = Path(__file__).resolve().parents[1]
LINEAR_SOURCE = ROOT / 'reports/TORIC_UNIVERSAL_THREE_PLANE_20261101.json'


def symmetric_tensor_action(a):
    """Action on symmetric tensor entries, without multinomial factors."""
    n = s.zeros(20)
    for row, triple in enumerate(TRIPLES):
        for position, i in enumerate(triple):
            for j in range(4):
                if a[j, i]:
                    changed = list(triple)
                    changed[position] = j
                    n[row, TRIPLES.index(tuple(sorted(changed)))] += a[j, i]
    return n


def sparse_action(matrix):
    return [[(j, QQ.from_sympy(matrix[i, j])) for j in range(matrix.cols)
             if matrix[i, j]] for i in range(matrix.rows)]


def derivative(f, action):
    answer = f.ring.zero
    for monomial, coefficient in f.items():
        for i, multiplicity in enumerate(monomial):
            if multiplicity:
                for j, value in action[i]:
                    exponent = list(monomial)
                    exponent[i] -= 1
                    exponent[j] += 1
                    key = tuple(exponent)
                    answer[key] = answer.get(key, QQ.zero) + coefficient * multiplicity * value
    answer.strip_zero()
    return answer


def evaluate(f, values):
    return sum((coefficient * math.prod(values[i] ** power
                for i, power in enumerate(monomial) if power)
                for monomial, coefficient in f.items()), QQ.zero)


def construct():
    basis = s.Matrix.hstack(*s.Matrix.vstack(*constraints()).nullspace())
    pivot_rows = list(basis.T.rref()[1])
    left = basis[pivot_rows, :].inv()

    def input_action(a):
        images = s.Matrix.hstack(*(infinitesimal_action(basis[:, i], a) for i in range(16)))
        answer = left * images[pivot_rows, :]
        assert basis * answer == images
        return answer

    standard_weights = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    weights = []
    for col in range(16):
        found = {tuple(standard_weights[i][a] + standard_weights[j][a]
                       + standard_weights[k][a] for a in range(2))
                 for value, ((i, j), k) in zip(basis[:, col], COORDS, strict=True) if value}
        assert len(found) == 1
        weights.append(found.pop())

    raising1 = s.zeros(4)
    raising1[0, 2], raising1[3, 1] = 1, -1
    raising2 = s.zeros(4)
    raising2[2, 3] = 1
    raising = [input_action(a) for a in (raising1, raising2)]
    assert all(a.T * J + J * a == s.zeros(4) for a in (raising1, raising2))

    monomials = {}
    character = Counter()
    for indices in combinations_with_replacement(range(16), 5):
        exponent = tuple(indices.count(i) for i in range(16))
        weight = tuple(sum(weights[i][a] for i in indices) for a in range(2))
        monomials.setdefault(weight, []).append(exponent)
        character[weight] += 1
    highest_weight, rho = (3, 0), (2, 1)
    multiplicity = 0
    for permutation in permutations((0, 1)):
        for signs in product((-1, 1), repeat=2):
            transformed_rho = tuple(signs[i] * rho[permutation[i]] for i in range(2))
            sign = (1 if permutation == (0, 1) else -1) * signs[0] * signs[1]
            multiplicity += sign * character[tuple(highest_weight[i] + rho[i]
                                                    - transformed_rho[i] for i in range(2))]
    source_monomials = monomials[highest_weight]
    blocks = []
    for action, weight in zip(raising, ((4, -1), (3, 2)), strict=True):
        targets = {monomial: i for i, monomial in enumerate(monomials[weight])}
        matrix = s.zeros(len(targets), len(source_monomials))
        for col, monomial in enumerate(source_monomials):
            for i, exponent in enumerate(monomial):
                if exponent:
                    for j in range(16):
                        if action[i, j]:
                            changed = list(monomial)
                            changed[i] -= 1
                            changed[j] += 1
                            matrix[targets[tuple(changed)], col] += exponent * action[i, j]
        blocks.append(matrix)
    highest_kernel = s.Matrix.vstack(*blocks).nullspace()
    assert len(highest_kernel) == multiplicity == 6

    r, *b = ring(','.join(f'b{i}' for i in range(16)), QQ)
    lowering1 = raising1.T
    lowering2 = raising2.T
    input_lowering = [sparse_action(input_action(a)) for a in (lowering1, lowering2)]
    output_lowering = [symmetric_tensor_action(a) for a in (lowering1, lowering2)]
    rows, words, queue = [s.eye(20)[0, :]], [()], [0]
    while queue and len(rows) < 20:
        index = queue.pop(0)
        for k, action in enumerate(output_lowering):
            row = rows[index] * action
            if s.Matrix.vstack(*rows, row).rank() > len(rows):
                rows.append(row)
                words.append(words[index] + (k,))
                queue.append(len(rows) - 1)
    assert len(rows) == 20
    change = s.Matrix.vstack(*rows).inv()
    covariants = []
    for vector in highest_kernel:
        forms = {(): r.from_dict({monomial: QQ.from_sympy(coefficient)
                                  for monomial, coefficient in zip(source_monomials, vector,
                                                                  strict=True) if coefficient})}
        for word in words[1:]:
            forms[word] = derivative(forms[word[:-1]], input_lowering[word[-1]])
        values = [forms[word] for word in words]
        covariants.append([sum((QQ.from_sympy(change[i, j]) * values[j]
                               for j in range(20) if change[i, j]), r.zero) for i in range(20)])

    # Verify all six maps, rather than relying only on highest-weight theory.
    for a in sp_basis():
        input_op = sparse_action(input_action(a))
        output_op = symmetric_tensor_action(a)
        for covariant in covariants:
            for i in range(20):
                assert derivative(covariant[i], input_op) == sum((QQ.from_sympy(output_op[i, j])
                    * covariant[j] for j in range(20) if output_op[i, j]), r.zero)

    z, u, beta = primitive_beta(S61())
    flat = s.Matrix([s.diff(beta[i, j], z[k]) for (i, j), k in COORDS])
    b_in_u = left * flat[pivot_rows, 0].jacobian(u)
    assert basis * b_in_u == flat.jacobian(u)
    u_in_b = b_in_u.inv()
    substitution = {symbol: sum(u_in_b[i, j] * s.Symbol(f'b{j}') for j in range(16))
                    for i, symbol in enumerate(u)}
    linear_record = json.loads(LINEAR_SOURCE.read_text())
    assert linear_record['beta_variables'] == list(map(str, u))
    linear = s.Matrix(linear_record['all51_homogeneous_linear_consequences'])
    linear_b = [[r.from_expr(s.expand(linear[i, j].subs(substitution, simultaneous=True)))
                 for j in range(20)] for i in range(51)]
    test_b = list(range(1, 17))
    evaluated_l = s.Matrix([[s.Rational(evaluate(f, test_b)) for f in row] for row in linear_b])
    evaluated_maps = s.Matrix.hstack(*(s.Matrix([s.Rational(evaluate(f, test_b)) for f in covariant])
                                      for covariant in covariants))
    constraints_on_maps = evaluated_l * evaluated_maps
    kernel = constraints_on_maps.nullspace()
    assert len(kernel) == 1
    denominator = s.ilcm(*[v.q for v in kernel[0]])
    integers = [int(v * denominator) for v in kernel[0]]
    divisor = math.gcd(*integers)
    integers = [v // divisor for v in integers]
    assert integers == [-46, -11, -48, -90, 30, 40]
    # This scale agrees with Tr(c X) in two independent E[2] cubic examples.
    coefficients = [-9 * value for value in integers]
    c5 = [sum((QQ(a) * covariant[i] for a, covariant in zip(coefficients, covariants,
                                                           strict=True) if a), r.zero)
          for i in range(20)]
    assert any(c5)
    assert all(set(map(sum, f)) == {5} for f in c5)
    assert all(not sum((a * b for a, b in zip(row, c5, strict=True)), r.zero) for row in linear_b)

    # This first test vector was independently obtained as a trace in an étale
    # cubic algebra, not used to determine the one-dimensional kernel direction.
    residue_u = list(range(1, 17))
    expected_residue = [-70399062, 31236138, 38352672, -74420748, -12597516,
        52856826, -26759490, 13422786, 44744802, 18882354, -12246912,
        20286792, -11313840, -13744650, 3120486, -20196570, 237146922,
        -99516438, 13566846, -75751074]
    assert [evaluate(f, b_in_u * s.Matrix(residue_u)) for f in c5] == expected_residue
    pivot_columns = list(constraints_on_maps.rref()[1])
    pivot_rows_constraints = list(constraints_on_maps[:, pivot_columns].T.rref()[1])
    minor = constraints_on_maps.extract(pivot_rows_constraints, pivot_columns)
    assert minor.shape == (5, 5) and minor.det()
    report = {
        'all_passed': True, 'field': 'QQ',
        'weight_basis_24_by16': [[str(x) for x in row] for row in basis.tolist()],
        'weight_coordinates_b': list(map(str, b)), 'weights': weights,
        'beta_variables': list(map(str, u)),
        'b_in_u': [[str(x) for x in row] for row in b_in_u.tolist()],
        'c_variables': linear_record['c_variables'], 'symmetric_tensor_triples': TRIPLES,
        'tensor_convention': 'c_ijk are symmetric tensor entries, not cubic monomial coefficients',
        'highest_weight': highest_weight, 'highest_weight_monomial_count': len(source_monomials),
        'raising_matrix_shape': [sum(block.rows for block in blocks), len(source_monomials)],
        'Hom_Sym5W_to_Sym3V_dimension': multiplicity,
        'highest_weight_monomials': source_monomials,
        'highest_weight_kernel_sparse': [[[i, str(x)] for i, x in enumerate(v) if x]
                                         for v in highest_kernel],
        'lowering_words': words,
        'all_six_covariants_all_ten_Lie_generators_verified': True,
        'infinitesimal_equivariance_polynomial_identities_checked': 1200,
        'kernel_test_b': test_b,
        'constraint_on_six_maps_rank': 5,
        'constraint_minor_rows': pivot_rows_constraints,
        'constraint_minor_columns': pivot_columns,
        'constraint_minor_determinant': str(minor.det()),
        'coefficients_in_six_covariants': coefficients,
        'equivariant_degree5_kernel_dimension': 1,
        'all51_universal_linear_kernel_identities': True,
        'C5_total_monomials': sum(map(len, c5)),
        'C5_sparse': [[[list(monomial), str(coefficient)]
                       for monomial, coefficient in sorted(f.items())] for f in c5],
        'residue_scale_test_beta': residue_u, 'residue_scale_test_C5': expected_residue,
        'scope': 'Unique degree5 Sp4-equivariant polynomial map in ker L, up to scalar. '
                 'This does not establish a universal identification of the toric cover with '
                 'nonzero two-torsion, nor full S61 fiber exhaustion.',
    }
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    report = construct()
    report['source_sha256'] = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (Path(__file__), LINEAR_SOURCE, ROOT / 'code/weight1_structure.py',
                     ROOT / 'code/primitive_jet_elliptic_curve_20261101.py',
                     ROOT / 'code/toric_algebra_fiber.py')}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps({'all_passed': True, 'Hom_dimension': 6,
                      'equivariant_kernel_dimension': 1,
                      'C5_monomials': report['C5_total_monomials'],
                      'seconds': time.time() - started, 'output': str(args.output)}, indent=2))


if __name__ == '__main__':
    main()

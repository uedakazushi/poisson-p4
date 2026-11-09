#!/usr/bin/env python3
"""Certify six generic-rank-two components of the FULL Poisson scheme of P^4.

Run with Python 3 and NumPy:
    python3 verify_rank2_components.py > ../data/rank2_certificates.json

All local imports are sibling files. No network access or files are required
apart from the supplied Python modules, and the only output is JSON on stdout.

For each example, J(delta)=E wedge [P,delta] is a 350-by-150 integer matrix.
V comprises the 25 generators E wedge A of the 24-dimensional vertical kernel.
A comprises actual derivatives of an explicitly parametrized irreducible family.
The exact identities Jacobi(P)=0, JV=JA=0 are checked. Nonzero minors modulo
1000003 certify rank lower bounds over Q. Consequently
    dim(full ambient tangent) <= 126-rank(J mod p),
    dim(family image) >= rank([V A] mod p)-24.
Equality proves smoothness of the chosen point in the full affine Poisson
scheme and proves that the family closure is an irreducible component.

The six affine/projective dimensions are:
 R113: 37/36; R122: 29/28; Log1112: 29/28;
 LPB: 40/39; Aff5: 23/22; Aff41: 23/22.
The two affine-Lie families use coordinate changes, scale, and relative block
weights in H. Their projective Jordan types distinguish their generic foliations.

The full-scheme assertions are certificates supplied here; they should not be
attributed to statements concerning only the rank-two Poisson locus.
References: https://arxiv.org/html/2207.12880v5 and
https://arxiv.org/pdf/2408.10057 (Lie-subalgebra components and foliations).
"""

from collections import defaultdict
import itertools
import json

import numpy as np

import poisson_p4_components_check as p
from lpb_family import build_lpb_case


def polynomial_field(matrix):
    """Convert a 5-by-5 matrix into the coefficients of its linear field."""
    return [
        {tuple(int(k == j) for k in range(5)): int(matrix[i, j])
         for j in range(5) if matrix[i, j]}
        for i in range(5)
    ]


def wedge_fields(first, second):
    """Wedge polynomial vector fields; also works with rational coefficients."""
    result = {}
    for i, j in p.PAIRS:
        coefficient = defaultdict(int)
        p.add_polynomial(coefficient,
                         p.polynomial_product(first[i], second[j]))
        p.add_polynomial(coefficient,
                         p.polynomial_product(first[j], second[i]), -1)
        result[(i, j)] = dict(coefficient)
    return result


def linear_action_derivative(bivector, matrix):
    """Lie derivative L_(matrix*x) of a quadratic bivector."""
    field = polynomial_field(matrix)
    result = {pair: defaultdict(int) for pair in p.PAIRS}
    for i, j in p.PAIRS:
        for k in range(5):
            p.add_polynomial(
                result[(i, j)],
                p.polynomial_product(
                    field[k], p.polynomial_derivative(bivector.get((i, j), {}), k)
                ),
            )
            p.add_polynomial(result[(i, j)], p.bivector_entry(bivector, k, j),
                             -int(matrix[i, k]))
            p.add_polynomial(result[(i, j)], p.bivector_entry(bivector, i, k),
                             -int(matrix[j, k]))
    return {pair: dict(polynomial) for pair, polynomial in result.items()}


def exterior_derivative(form):
    degree = len(next(iter(form)))
    result = {indices: defaultdict(int)
              for indices in itertools.combinations(range(5), degree + 1)}
    for indices, polynomial in form.items():
        for direction in range(5):
            extended = (direction,) + indices
            sign = p.permutation_sign(extended)
            if sign:
                p.add_polynomial(result[tuple(sorted(extended))],
                                 p.polynomial_derivative(polynomial, direction), sign)
    return result


def log_bivector(polynomials, parameters):
    """Log(1,1,1,2) using P dual to d(omega).

    omega=(product F_i)*sum B_ij dlog(F_i) wedge dlog(F_j),
    where B*(1,1,1,2)^T=0. For four logarithmic coordinates the horizontal
    coefficient matrix has rank at most two, so the induced bivector has rank 2.
    A constant normalization by 1/5 is unnecessary for the certificate.
    """
    a, b, c = parameters
    coefficients = {
        (0, 1): 4 * a, (0, 2): 4 * b, (1, 2): 4 * c,
        (0, 3): -2 * (a + b), (1, 3): 2 * (a - c),
        (2, 3): 2 * (b + c),
    }
    omega = {indices: defaultdict(int) for indices in p.PAIRS}
    for (i, j), scalar in coefficients.items():
        if not scalar:
            continue
        product = {(0,) * 5: 1}
        for k in range(4):
            if k not in (i, j):
                product = p.polynomial_product(product, polynomials[k])
        for u, v in p.PAIRS:
            term = defaultdict(int)
            p.add_polynomial(term, p.polynomial_product(
                p.polynomial_derivative(polynomials[i], u),
                p.polynomial_derivative(polynomials[j], v)))
            p.add_polynomial(term, p.polynomial_product(
                p.polynomial_derivative(polynomials[i], v),
                p.polynomial_derivative(polynomials[j], u)), -1)
            p.add_polynomial(omega[(u, v)],
                             p.polynomial_product(product, term), scalar)
    three_form = exterior_derivative(omega)
    bivector = {}
    for pair in p.PAIRS:
        complement = tuple(i for i in range(5) if i not in pair)
        sign = p.permutation_sign(pair + complement)
        bivector[pair] = {monomial: sign * coefficient
                         for monomial, coefficient in three_form[complement].items()
                         if coefficient}
    return bivector


def build_jacobian_cases():
    cubic = p.diagonal_polynomial(3, [1] * 5)
    cubic[(1, 0, 1, 1, 0)] = 1
    cubic[(0, 1, 0, 1, 1)] = 1
    first_quadric = p.diagonal_polynomial(2, [1] * 5)
    second_quadric = p.diagonal_polynomial(2, [1, 2, 3, 4, 5])
    specifications = [
        ('R113', [cubic, p.coordinate_polynomial(0), p.coordinate_polynomial(1)], 89, 61),
        ('R122', [first_quadric, second_quadric, p.coordinate_polynomial(0)], 97, 53),
    ]
    for name, polynomials, jacobi_rank, family_rank in specifications:
        yield {
            'name': name,
            'bivector': p.jacobian_bivector(*polynomials),
            'parameter_derivative': p.parameter_derivative_matrix(polynomials),
            'expected_jacobi_rank': jacobi_rank,
            'expected_vertical_plus_parameter_rank': family_rank,
            'affine_component_dimension': family_rank - 24,
            'projective_component_dimension': family_rank - 25,
        }


def build_log_case():
    quadric = p.diagonal_polynomial(2, [1, 2, 3, 5, 7])
    quadric[(1, 0, 0, 1, 0)] = 1
    quadric[(0, 1, 0, 0, 1)] = 1
    polynomials = [p.coordinate_polynomial(i) for i in range(3)] + [quadric]
    parameters = [1, 7, 17]
    columns = []
    # The construction is multilinear in the four polynomials and in B.
    for k, polynomial in enumerate(polynomials):
        degree = sum(next(iter(polynomial)))
        for monomial in p.homogeneous_monomials(degree):
            varied = polynomials.copy()
            varied[k] = {monomial: 1}
            columns.append(p.bivector_vector(log_bivector(varied, parameters)))
    for k in range(3):
        columns.append(p.bivector_vector(log_bivector(
            polynomials, [int(i == k) for i in range(3)])))
    return {
        'name': 'Log1112',
        'bivector': log_bivector(polynomials, parameters),
        'parameter_derivative': np.column_stack(columns),
        'expected_jacobi_rank': 97,
        'expected_vertical_plus_parameter_rank': 53,
        'affine_component_dimension': 29,
        'projective_component_dimension': 28,
    }


def build_affine_case(blocks):
    nilpotent = np.zeros((5, 5), dtype=np.int64)
    grading = np.zeros((5, 5), dtype=np.int64)
    offset = 0
    for block_number, size in enumerate(blocks):
        for i in range(size):
            grading[offset + i, offset + i] = i + 7 * block_number
            if i + 1 < size:
                nilpotent[offset + i + 1, offset + i] = 1
        offset += size
    bivector = wedge_fields(polynomial_field(nilpotent), polynomial_field(grading))
    columns = [p.bivector_vector(bivector)]
    for i, j in itertools.product(range(5), repeat=2):
        matrix = np.zeros((5, 5), dtype=np.int64)
        matrix[i, j] = 1
        columns.append(p.bivector_vector(linear_action_derivative(bivector, matrix)))
    offset = 0
    for size in blocks:
        block_scalar = np.zeros((5, 5), dtype=np.int64)
        for i in range(offset, offset + size):
            block_scalar[i, i] = 1
        columns.append(p.bivector_vector(wedge_fields(
            polynomial_field(nilpotent), polynomial_field(block_scalar))))
        offset += size
    return {
        'name': 'Aff' + ''.join(map(str, blocks)),
        'jordan_blocks': list(blocks),
        'grading_diagonal': grading.diagonal().tolist(),
        'bivector': bivector,
        'parameter_derivative': np.column_stack(columns),
        'expected_jacobi_rank': 103,
        'expected_vertical_plus_parameter_rank': 47,
        'affine_component_dimension': 23,
        'projective_component_dimension': 22,
    }


def assert_integer_product_zero(first, second):
    """Bound every integer matrix product entry before using int64 arithmetic."""
    bound = (max(abs(int(v)) for v in first.flat)
             * max(abs(int(v)) for v in second.flat) * first.shape[1])
    assert bound < np.iinfo(np.int64).max
    assert not np.any(first @ second)


def verify_case(case, vertical):
    bivector = case['bivector']
    derivatives = case['parameter_derivative']
    assert all(not value for value in p.affine_jacobi_coefficients(bivector).values())
    jacobi = p.linearized_jacobi_matrix(bivector)
    assert_integer_product_zero(jacobi, vertical)
    assert_integer_product_zero(jacobi, derivatives)
    jacobi_certificate = p.modular_rank_certificate(jacobi)
    family_certificate = p.modular_rank_certificate(np.column_stack([vertical, derivatives]))
    assert jacobi_certificate['rank_mod_prime'] == case['expected_jacobi_rank']
    assert family_certificate['rank_mod_prime'] == case['expected_vertical_plus_parameter_rank']
    tangent_bound = 126 - jacobi_certificate['rank_mod_prime']
    family_bound = family_certificate['rank_mod_prime'] - 24
    assert tangent_bound == family_bound == case['affine_component_dimension']
    metadata = {key: value for key, value in case.items()
                if key not in ('bivector', 'parameter_derivative')}
    metadata.update({
        'affine_jacobi_exact_zero': True,
        'J_vertical_exact_zero': True,
        'J_parameter_derivative_exact_zero': True,
        'ambient_affine_tangent_dimension': tangent_bound,
        'family_affine_differential_rank': family_bound,
        'jacobi_matrix_certificate': jacobi_certificate,
        'vertical_and_family_matrix_certificate': family_certificate,
    })
    return metadata


def main():
    vertical = p.vertical_bivector_matrix()
    vertical_certificate = p.modular_rank_certificate(vertical)
    assert vertical_certificate['rank_mod_prime'] == 24
    cases = list(build_jacobian_cases()) + [
        build_log_case(), build_lpb_case(p),
        build_affine_case((5,)), build_affine_case((4, 1)),
    ]
    report = {
        'field': 'C; integer matrices certified via nonzero minors modulo a prime',
        'prime': p.PRIME,
        'quadratic_affine_bivector_dimension': 150,
        'projective_global_bivector_dimension': 126,
        'vertical_certificate': vertical_certificate,
        'components': [verify_case(case, vertical) for case in cases],
    }
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()

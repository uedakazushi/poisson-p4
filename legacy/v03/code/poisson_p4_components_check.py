#!/usr/bin/env python3
"""Exact certificates for two irreducible components of Pois(P^4_C).

Requirements and usage
----------------------
Python 3 and NumPy only. Run:

    python3 poisson_p4_components_check.py

The script prints one JSON report and writes no files. All polynomial operations
use integers; Gaussian elimination uses the prime field F_1000003. No floating
point calculation, random example, or external computer algebra system is used.

Construction and examples
-------------------------
On C^5 let Omega = dx0 ^ ... ^ dx4 and E = sum xi*d/dxi. Define a quadratic
bivector P by contraction with Omega:

    i_P Omega = dF ^ dG ^ dH,       deg(F) + deg(G) + deg(H) = 5.

These Jacobian brackets are Poisson, have rank at most two, and descend to P^4.
The two points checked here are:

  * Cubic: F = sum xi^3 + x0*x2*x3 + x1*x3*x4, G = x0, H = x1.
  * Quadrics: F = sum xi^2, G = sum (i+1)*xi^2, H = x0.

The cubic is smooth over C: at a nonzero common zero of its derivatives, take
a coordinate with maximum absolute value M. Its derivative would equate a term
of absolute value 3*M^2 with a sum of absolute value at most 2*M^2. Both quadrics
are smooth because all their diagonal coefficients are nonzero.

What is certified
-----------------
There are 150 homogeneous quadratic bivectors on C^5. Their quotient by the
24-dimensional space of vertical bivectors E ^ A, with A a linear vector field,
is H^0(P^4, wedge^2 T), of dimension 126. In this quotient the linearized Poisson
equation at P is exactly

    J_P(delta) = E ^ [P, delta] = 0.

We express J_P as a 350-by-150 integer matrix. Let V be the matrix of the 25
generators E ^ (xj*d/dxi), whose span has dimension 24. Let A be the matrix of
derivatives of the Jacobian construction with respect to every coefficient of
F, G, and H. The script verifies the full affine Jacobi identity and the integer
identities J_P*V = J_P*A = 0.

The ranks of J_P modulo 1000003 are 89 and 97. Nonzero minors therefore give
upper bounds 37 and 29 for the characteristic-zero ambient tangent dimensions.
The ranks of [V A] modulo the same prime are 61 and 53. Their nonzero minors give
lower bounds 37 and 29 for the characteristic-zero differentials of the parameter
maps after quotienting by V. The parameter spaces are smooth and irreducible,
so their image closures have at least these dimensions. Equality with the tangent
upper bounds proves that the chosen points are smooth points of the FULL affine
Poisson scheme and that the two image closures are irreducible components.

Thus the affine component dimensions are 37 and 29; their projectivizations have
dimensions 36 and 28. This does not identify the Poisson space of one fixed
hypersurface with a component: the hypersurface data vary in these constructions.

The minor determinants in the report are determinants of submatrices of the
original integer matrices, reduced modulo the prime. Thus the computation is a
characteristic-zero rank certificate, rather than a numerical rank estimate.

References and scope
--------------------
* Brent Pym, Poisson Structures and Lie Algebroids in Complex Geometry,
  Sections 7.3--7.4 (especially Example 7.4.5):
  https://www.math.mcgill.ca/bpym/files/bpym-thesis.pdf
* Mauricio Correa and Alan Muniz, Holomorphic foliations of degree two and
  arbitrary dimension, v5 (2026), Example 2.2 and Section 5:
  https://arxiv.org/html/2207.12880v5

The latter reference discusses rational foliation components and the rank-two
Poisson locus. The full-Poisson-scheme assertion above is justified here by the
separate tangent-space certificates, not attributed to that reference.
"""

from collections import defaultdict
import itertools
import json

import numpy as np


NUM_VARIABLES = 5
PRIME = 1_000_003
PAIRS = list(itertools.combinations(range(NUM_VARIABLES), 2))
TRIPLES = list(itertools.combinations(range(NUM_VARIABLES), 3))
QUADRUPLES = list(itertools.combinations(range(NUM_VARIABLES), 4))


def homogeneous_monomials(degree):
    """Exponent vectors, ordered by combinations with replacement of variables."""
    result = []
    for factors in itertools.combinations_with_replacement(
        range(NUM_VARIABLES), degree
    ):
        exponents = [0] * NUM_VARIABLES
        for variable in factors:
            exponents[variable] += 1
        result.append(tuple(exponents))
    return result


DEGREE_TWO = homogeneous_monomials(2)
DEGREE_FOUR = homogeneous_monomials(4)
BIVECTOR_BASIS = list(itertools.product(PAIRS, DEGREE_TWO))
BIVECTOR_INDEX = {basis: index for index, basis in enumerate(BIVECTOR_BASIS)}
FOUR_VECTOR_BASIS = list(itertools.product(QUADRUPLES, DEGREE_FOUR))
FOUR_VECTOR_INDEX = {
    basis: index for index, basis in enumerate(FOUR_VECTOR_BASIS)
}


def polynomial_derivative(polynomial, variable):
    result = {}
    for exponents, coefficient in polynomial.items():
        exponent = exponents[variable]
        if exponent:
            differentiated = list(exponents)
            differentiated[variable] -= 1
            result[tuple(differentiated)] = coefficient * exponent
    return result


def polynomial_product(first, second):
    result = defaultdict(int)
    for first_exponents, first_coefficient in first.items():
        for second_exponents, second_coefficient in second.items():
            exponents = tuple(
                a + b for a, b in zip(first_exponents, second_exponents)
            )
            result[exponents] += first_coefficient * second_coefficient
    return {exponents: coefficient for exponents, coefficient in result.items()
            if coefficient}


def add_polynomial(destination, source, scalar=1):
    for exponents, coefficient in source.items():
        destination[exponents] += scalar * coefficient


def permutation_sign(indices):
    if len(set(indices)) != len(indices):
        return 0
    inversions = sum(
        indices[i] > indices[j]
        for i in range(len(indices))
        for j in range(i + 1, len(indices))
    )
    return -1 if inversions % 2 else 1


def bivector_entry(bivector, first, second):
    if first < second:
        return bivector.get((first, second), {})
    if first > second:
        return {exponents: -coefficient for exponents, coefficient
                in bivector.get((second, first), {}).items()}
    return {}


def jacobian_bivector(first, second, third):
    gradients = [
        [polynomial_derivative(polynomial, variable)
         for variable in range(NUM_VARIABLES)]
        for polynomial in (first, second, third)
    ]
    bivector = {}
    for pair in PAIRS:
        remaining = [variable for variable in range(NUM_VARIABLES)
                     if variable not in pair]
        coefficient = defaultdict(int)
        for permutation in itertools.permutations(remaining):
            term = polynomial_product(
                polynomial_product(gradients[0][permutation[0]],
                                   gradients[1][permutation[1]]),
                gradients[2][permutation[2]],
            )
            add_polynomial(coefficient, term, permutation_sign(pair + permutation))
        bivector[pair] = dict(coefficient)
    return bivector


def bivector_vector(bivector):
    result = np.zeros(len(BIVECTOR_BASIS), dtype=np.int64)
    for pair, polynomial in bivector.items():
        for exponents, coefficient in polynomial.items():
            if coefficient:
                result[BIVECTOR_INDEX[(pair, exponents)]] += coefficient
    return result


def affine_jacobi_coefficients(bivector):
    """Coefficients of {xi,{xj,xk}} + cyclic terms, over the integers."""
    result = {}
    for triple in TRIPLES:
        coefficient = defaultdict(int)
        cyclic_orders = (triple, (triple[1], triple[2], triple[0]),
                         (triple[2], triple[0], triple[1]))
        for first, second, third in cyclic_orders:
            for variable in range(NUM_VARIABLES):
                add_polynomial(coefficient, polynomial_product(
                    bivector_entry(bivector, first, variable),
                    polynomial_derivative(
                        bivector_entry(bivector, second, third), variable),
                ))
        result[triple] = {exponents: value for exponents, value
                          in coefficient.items() if value}
    return result


def linearized_jacobi_matrix(bivector):
    """Matrix of delta -> E wedge [P, delta] in the indicated monomial bases."""
    matrix = np.zeros((len(FOUR_VECTOR_BASIS), len(BIVECTOR_BASIS)),
                      dtype=np.int64)
    for column, (pair, exponents) in enumerate(BIVECTOR_BASIS):
        variation = {pair: {exponents: 1}}
        for triple in TRIPLES:
            coefficient = defaultdict(int)
            cyclic_orders = (triple, (triple[1], triple[2], triple[0]),
                             (triple[2], triple[0], triple[1]))
            for first, second, third in cyclic_orders:
                for variable in range(NUM_VARIABLES):
                    add_polynomial(coefficient, polynomial_product(
                        bivector_entry(bivector, first, variable),
                        polynomial_derivative(
                            bivector_entry(variation, second, third), variable),
                    ))
                    add_polynomial(coefficient, polynomial_product(
                        bivector_entry(variation, first, variable),
                        polynomial_derivative(
                            bivector_entry(bivector, second, third), variable),
                    ))
            for variable in range(NUM_VARIABLES):
                exterior_indices = (variable,) + triple
                sign = permutation_sign(exterior_indices)
                if not sign:
                    continue
                quadruple = tuple(sorted(exterior_indices))
                for monomial, value in coefficient.items():
                    multiplied = list(monomial)
                    multiplied[variable] += 1
                    row = FOUR_VECTOR_INDEX[(quadruple, tuple(multiplied))]
                    matrix[row, column] += sign * value
    return matrix


def vertical_bivector_matrix():
    """Columns are the 25 bivectors E wedge (xj*d/dxi)."""
    columns = []
    for direction, coefficient_variable in itertools.product(range(5), repeat=2):
        bivector = {pair: defaultdict(int) for pair in PAIRS}
        for euler_direction in range(NUM_VARIABLES):
            if euler_direction == direction:
                continue
            exponents = [0] * NUM_VARIABLES
            exponents[euler_direction] += 1
            exponents[coefficient_variable] += 1
            pair = tuple(sorted((euler_direction, direction)))
            sign = 1 if euler_direction < direction else -1
            bivector[pair][tuple(exponents)] += sign
        columns.append(bivector_vector(bivector))
    return np.column_stack(columns)


def parameter_derivative_matrix(polynomials):
    """Derivative with respect to every coefficient of each defining polynomial."""
    columns = []
    for polynomial_index, polynomial in enumerate(polynomials):
        degree = sum(next(iter(polynomial)))
        for exponents in homogeneous_monomials(degree):
            variation = list(polynomials)
            variation[polynomial_index] = {exponents: 1}
            columns.append(bivector_vector(jacobian_bivector(*variation)))
    return np.column_stack(columns)


def determinant_mod_prime(matrix, prime):
    """Exact modular determinant; entries stay below prime squared < 2^63."""
    matrix = matrix.copy() % prime
    determinant = 1
    for column in range(matrix.shape[0]):
        candidates = np.flatnonzero(matrix[column:, column])
        if not len(candidates):
            return 0
        pivot_row = int(candidates[0]) + column
        if pivot_row != column:
            matrix[[column, pivot_row]] = matrix[[pivot_row, column]]
            determinant = -determinant
        pivot = int(matrix[column, column])
        determinant = determinant * pivot % prime
        if column + 1 < matrix.shape[0]:
            factors = matrix[column + 1:, column] * pow(pivot, -1, prime) % prime
            matrix[column + 1:] = (
                matrix[column + 1:] - factors[:, None] * matrix[column][None, :]
            ) % prime
    return int(determinant % prime)


def modular_rank_certificate(matrix, prime=PRIME):
    """Rank and a nonzero minor from the original matrix, both modulo prime."""
    reduced = matrix.copy() % prime
    original_rows = list(range(matrix.shape[0]))
    minor_rows, minor_columns = [], []
    rank = 0
    for column in range(matrix.shape[1]):
        candidates = np.flatnonzero(reduced[rank:, column])
        if not len(candidates):
            continue
        pivot_row = int(candidates[0]) + rank
        if pivot_row != rank:
            reduced[[rank, pivot_row]] = reduced[[pivot_row, rank]]
            original_rows[rank], original_rows[pivot_row] = (
                original_rows[pivot_row], original_rows[rank]
            )
        minor_rows.append(original_rows[rank])
        minor_columns.append(column)
        pivot = int(reduced[rank, column])
        reduced[rank] = reduced[rank] * pow(pivot, -1, prime) % prime
        lower_rows = np.flatnonzero(reduced[rank + 1:, column]) + rank + 1
        reduced[lower_rows] = (
            reduced[lower_rows]
            - reduced[lower_rows, column, None] * reduced[rank][None, :]
        ) % prime
        rank += 1
        if rank == min(matrix.shape):
            break
    minor = matrix[np.ix_(minor_rows, minor_columns)]
    minor_determinant = determinant_mod_prime(minor, prime)
    assert minor_determinant != 0, "The claimed rank minor must be nonzero."
    return {
        "rank_mod_prime": rank,
        "nonzero_minor_size": rank,
        "minor_determinant_mod_prime": minor_determinant,
    }


def coordinate_polynomial(variable):
    exponents = [0] * NUM_VARIABLES
    exponents[variable] = 1
    return {tuple(exponents): 1}


def diagonal_polynomial(degree, coefficients):
    return {
        tuple(degree if position == variable else 0
              for position in range(NUM_VARIABLES)): coefficient
        for variable, coefficient in enumerate(coefficients)
    }


def verify_example(name, polynomials, vertical_matrix, expected_dimension):
    bivector = jacobian_bivector(*polynomials)
    assert all(not coefficient
               for coefficient in affine_jacobi_coefficients(bivector).values()), (
        "The full affine Jacobi identity failed over the integers."
    )
    jacobi_matrix = linearized_jacobi_matrix(bivector)
    parameter_matrix = parameter_derivative_matrix(polynomials)
    assert not np.any(jacobi_matrix @ vertical_matrix), "J*V must vanish exactly."
    assert not np.any(jacobi_matrix @ parameter_matrix), "J*A must vanish exactly."

    jacobi_certificate = modular_rank_certificate(jacobi_matrix)
    family_certificate = modular_rank_certificate(
        np.column_stack((vertical_matrix, parameter_matrix))
    )
    tangent_upper_bound = 126 - jacobi_certificate["rank_mod_prime"]
    family_lower_bound = family_certificate["rank_mod_prime"] - 24
    assert tangent_upper_bound == family_lower_bound == expected_dimension
    return {
        "example": name,
        "exact_integer_checks": {
            "full_affine_Jacobi_zero": True,
            "J_times_vertical_matrix_zero": True,
            "J_times_parameter_derivatives_zero": True,
        },
        "linearized_Jacobi_matrix_shape": list(jacobi_matrix.shape),
        "linearized_Jacobi_certificate": jacobi_certificate,
        "vertical_plus_parameter_certificate": family_certificate,
        "characteristic_zero_tangent_upper_bound": tangent_upper_bound,
        "characteristic_zero_family_dimension_lower_bound": family_lower_bound,
        "affine_component_dimension": expected_dimension,
        "projectivized_component_dimension": expected_dimension - 1,
        "conclusion": "Smooth point of the full affine Poisson scheme; "
                      "the Jacobian family closure is an irreducible component.",
    }


def main():
    vertical_matrix = vertical_bivector_matrix()
    vertical_certificate = modular_rank_certificate(vertical_matrix)
    assert vertical_certificate["rank_mod_prime"] == 24

    cubic = diagonal_polynomial(3, [1] * NUM_VARIABLES)
    cubic[(1, 0, 1, 1, 0)] = 1
    cubic[(0, 1, 0, 1, 1)] = 1
    first_quadric = diagonal_polynomial(2, [1] * NUM_VARIABLES)
    second_quadric = diagonal_polynomial(2, [1, 2, 3, 4, 5])
    report = {
        "prime": PRIME,
        "quadratic_bivector_dimension_on_C5": 150,
        "vertical_certificate": vertical_certificate,
        "global_bivector_dimension_on_P4": 126,
        "examples": [
            verify_example(
                "(3,1,1): sum xi^3 + x0*x2*x3 + x1*x3*x4; x0; x1",
                [cubic, coordinate_polynomial(0), coordinate_polynomial(1)],
                vertical_matrix, 37,
            ),
            verify_example(
                "(2,2,1): sum xi^2; sum (i+1)*xi^2; x0",
                [first_quadric, second_quadric, coordinate_polynomial(0)],
                vertical_matrix, 29,
            ),
        ],
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

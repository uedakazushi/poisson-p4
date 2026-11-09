"""LPB rank-two component: case builder for the exact P^4 Poisson checker.

This module does not execute a computation or read/write files on import.
Call ``build_lpb_case(backend)`` with the module poisson_p4_components_check as
backend. It returns the quadratic bivector, all parameter derivatives needed
for the tangent sandwich, and the expected dimensions.

Geometric family
----------------
Choose a projection center p in P^4 and a quadratic vector field on the
four-dimensional quotient vector space, modulo radial fields l*E_4. The vector
field parameter space has dimension 4*10-4=36, while the center has dimension 4.
Its associated bivector is e_4 wedge v in a chart centered at [e_4]. Thus the
family has affine dimension 40 and projective dimension 39 once the differential
is seen to have rank 40.

The example is the Jouanolou vector field
v=x1^2*d/dx0+x2^2*d/dx1+x3^2*d/dx2+x0^2*d/dx3.
It is independent of x4, so e_4 wedge v is Poisson. The 40 coefficient
variations of v and four center variations give a concrete parameter derivative.
For center variation j, use g_t(x)=x+t*x4*e_j. The induced derivative is
    e_j wedge v - e_4 wedge (x4 * partial_j(v)).

Certified once at prime 1000003 using the backend:
J rank 86, nonzero minor determinant 607350;
[V A] rank 64, nonzero minor determinant 1000002;
V rank 24. Hence tangent upper bound 40 equals family lower bound 40.
The full affine Jacobi identity and J*V=J*A=0 vanish over the integers.
An independent implementation gives J rank 86 also modulo 101.
"""

from collections import defaultdict
import numpy as np


def build_lpb_case(backend):
    """Return a case dict; no assumptions about a persistent filesystem path."""
    p = backend

    def wedge_with_e4(vector_field):
        bivector = {pair: {} for pair in p.PAIRS}
        for i in range(4):
            bivector[(i, 4)] = {
                exponent: -coefficient
                for exponent, coefficient in vector_field[i].items()
            }
        return bivector

    def add_wedge(bivector, direction, vector_field, factor=1):
        for i, polynomial in enumerate(vector_field):
            if i == direction:
                continue
            pair = tuple(sorted((direction, i)))
            sign = factor * (1 if direction < i else -1)
            for exponent, coefficient in polynomial.items():
                bivector[pair][exponent] += sign * coefficient

    vector_field = []
    for i in range(4):
        exponent = [0] * 5
        exponent[(i + 1) % 4] = 2
        vector_field.append({tuple(exponent): 1})
    bivector = wedge_with_e4(vector_field)

    columns = []
    for i in range(4):
        for exponent in p.DEGREE_TWO:
            if exponent[4]:
                continue
            variation = [{} for _ in range(4)]
            variation[i] = {exponent: 1}
            columns.append(p.bivector_vector(wedge_with_e4(variation)))

    for direction in range(4):
        variation = {pair: defaultdict(int) for pair in p.PAIRS}
        add_wedge(variation, direction, vector_field)
        differentiated = [
            p.polynomial_product(
                p.coordinate_polynomial(4),
                p.polynomial_derivative(polynomial, direction),
            )
            for polynomial in vector_field
        ]
        add_wedge(variation, 4, differentiated, -1)
        columns.append(p.bivector_vector(variation))

    return {
        'name': 'LPB',
        'bivector': bivector,
        'parameter_derivative': np.column_stack(columns),
        'expected_jacobi_rank': 86,
        'expected_vertical_plus_parameter_rank': 64,
        'affine_component_dimension': 40,
        'projective_component_dimension': 39,
        'coefficient_parameter_count': 40,
        'center_parameter_count': 4,
    }

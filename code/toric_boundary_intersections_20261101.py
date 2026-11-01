#!/usr/bin/env python3
"""A universal T01/T40 boundary intersection and an eight-dimensional missing orbit."""
import argparse
import hashlib
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path

import sympy as s

from s61_model import S61, X
from weight1_structure import J, sp_basis

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    l0, l1, l2 = s.symbols('l0 l1 l2')
    ell = s.Matrix([l0, l1, l2, 1])
    e = s.eye(4)[:, 3]
    trace_form = 4 * ell * ell.T
    multiplication = [ell[i] * s.eye(4) + s.eye(4)[:, i] * ell.T - 5 * ell[i] * e * ell.T
                      for i in range(4)]
    regular = [s.Matrix([[0]]).row_join(trace_form[i, :]).col_join(
               s.eye(4)[:, i].row_join(multiplication[i])) for i in range(4)]
    fields = s.Matrix.hstack(*(matrix * s.Matrix(X) for matrix in regular))
    homogeneous = (fields * J * fields.T).applyfunc(s.expand)
    model = S61()
    normalized = model.normalize(homogeneous)
    parameters = normalized['parameters']
    one = [i for i, weight in enumerate(model.reduced['parameter_weights']) if weight == 1]
    two = [i for i, weight in enumerate(model.reduced['parameter_weights']) if weight == 2]
    u = s.symbols('u1:62')
    substitution = {u[i]: parameters[i] for i in one}
    q52_path = ROOT / 'reports/Q52_UNIVERSAL_SECTION_20261101.json'
    q51_path = ROOT / 'reports/Q51_UNIVERSAL_SECTION_20261101.json'
    q52 = json.loads(q52_path.read_text())
    assert all(s.expand(s.sympify(value).subs(substitution) - parameters[i]) == 0
               for value, i in zip(q52['quadratic_section_q'], two, strict=True))
    q51 = json.loads(q51_path.read_text())
    difference = [s.expand(s.sympify(value).subs(substitution) - parameters[i])
                  for value, i in zip(q51['quadratic_section_q'], two, strict=True)]
    sample = {l0: 1, l1: 2, l2: 3}
    nonzero = next(i for i, value in enumerate(difference) if value.subs(sample))

    # Independent orbit computation on the entire quadratic affine bivector.
    y = s.Matrix(s.symbols('y0:4'))
    quadratic = s.zeros(4)
    quadratic[0, 1], quadratic[1, 0] = -y[0] ** 2, y[0] ** 2
    quadratic[2, 3], quadratic[3, 2] = -y[2] ** 2, y[2] ** 2
    pairs = list(combinations(range(4), 2))
    monomials = [y[i] * y[j] for i, j in combinations_with_replacement(range(4), 2)]
    columns = []
    for action in sp_basis():
        derivative = action * quadratic + quadratic * action.T - sum(
            ((action * y)[k] * quadratic.diff(y[k]) for k in range(4)), s.zeros(4))
        columns.append(s.Matrix([s.Poly(s.expand(derivative[i, j]), *y).coeff_monomial(monomial)
                                 for i, j in pairs for monomial in monomials]))
    tangent = s.Matrix.hstack(*columns)
    pivot_columns = list(tangent.rref()[1])
    pivot_rows = list(tangent[:, pivot_columns].T.rref()[1])
    assert len(pivot_rows) == len(pivot_columns) == 8
    minor = tangent.extract(pivot_rows, pivot_columns)
    assert minor.det()
    missing_path = ROOT / 'reports/TORIC_MISSING_BOUNDARY_20261101.json'
    missing = json.loads(missing_path.read_text())
    boundary = [s.sympify(value).subs({s.Symbol('a'): 1, s.Symbol('b'): 1})
                for value in missing['boundary_S61_coordinates']]
    normalized_quadratic = s.zeros(5)
    for i, j in pairs:
        normalized_quadratic[i + 1, j + 1] = s.expand(J[i, j] * X[0] ** 2
            + quadratic[i, j].subs(dict(zip(y, X[1:5], strict=True)), simultaneous=True))
        normalized_quadratic[j + 1, i + 1] = -normalized_quadratic[i + 1, j + 1]
    assert model.normalize(normalized_quadratic)['parameters'] == boundary
    report = {
        'all_passed': True, 'field': 'QQ[l0,l1,l2] and exact QQ orbit ranks',
        'A1_A4sq_parameterization': {'e': list(map(str, e)), 'ell': list(map(str, ell))},
        'A1_A4sq_S61_coordinates': list(map(str, parameters)),
        'all45_universal_Q52_graph_identities': True,
        'A1_A4sq_image_closure_is_six_dimensional_subvariety_of_T01_intersect_T40': True,
        'not_contained_in_T39_sample': {'ell': [1, 2, 3, 1], 'q_index_zero_based': nonzero,
                                      'Q51_minus_actual_q': str(difference[nonzero].subs(sample))},
        'missing_zero_jet_quadratic': [[str(value) for value in row] for row in quadratic.tolist()],
        'orbit_tangent_matrix': [[str(value) for value in row] for row in tangent.tolist()],
        'orbit_rank': 8, 'orbit_minor_rows': pivot_rows, 'orbit_minor_columns': pivot_columns,
        'orbit_minor_determinant': str(minor.det()),
        'actual_T01_zero_jet_fiber_dimension_lower_bound': 8,
        'entire_eight_dimensional_orbit_is_outside_the_finite_algebra_model_image': True,
        'proof_note': 'research/notes/TORIC_BOUNDARY_INTERSECTIONS_20261101.md',
        'limitations': ['No upper bound or component decomposition of the T01 zero-jet fiber.',
                        'No assertion that the six-dimensional subvariety is an intersection component.',
                        'No full S61 exhaustion or full T01 boundary classification.'],
        'source_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in (Path(__file__), q52_path, q51_path, missing_path,
                                      ROOT / 'code/s61_model.py', ROOT / 'code/weight1_structure.py',
                                      ROOT / 'reports/TORIC_ALGEBRA_BOUNDARY_FIBERS_20261101.json')},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps({'all_passed': True, 'universal_product_boundary_in_T40': True,
                      'missing_zero_jet_orbit_dimension': 8, 'minor': str(minor.det())}))


if __name__ == '__main__':
    main()

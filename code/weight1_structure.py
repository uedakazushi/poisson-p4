#!/usr/bin/env python3
"""Exact linear algebra for the primitive weight-one representation of S61.

The source is independent of the archived S61 elimination. All ranks and minors
are computed over QQ. The generic stabilizer calculation concerns the entire
16-dimensional first-jet representation, not the image of the Poisson scheme.
Run: .venv/bin/python code/weight1_structure.py --output reports/weight1.json
"""
from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path

import sympy as s

PAIRS = list(combinations(range(4), 2))
TRIPLES = list(combinations(range(4), 3))
COORDS = [(pair, k) for pair in PAIRS for k in range(4)]
J = s.Matrix([[0, 1, 0, 0], [-1, 0, 0, 0],
              [0, 0, 0, 1], [0, 0, -1, 0]])


def sign(values):
    if len(set(values)) != len(values):
        return 0
    return (-1) ** sum(values[i] > values[j]
                       for i in range(len(values))
                       for j in range(i + 1, len(values)))


def constraints():
    """d beta = 0 and Lambda beta = 0, on linear two-forms."""
    d = s.zeros(4, 24)
    trace = s.zeros(4, 24)
    for col, ((i, j), k) in enumerate(COORDS):
        if k not in (i, j):
            d[TRIPLES.index(tuple(sorted((k, i, j)))), col] = sign((k, i, j))
        if (i, j) in ((0, 1), (2, 3)):
            trace[k, col] = 1
    return d, trace


def form_matrices(vector):
    matrices = [s.zeros(4) for _ in range(4)]
    for value, ((i, j), k) in zip(vector, COORDS):
        matrices[k][i, j] = value
        matrices[k][j, i] = -value
    return matrices


def flatten(matrices):
    return s.Matrix([matrices[k][i, j] for (i, j), k in COORDS])


def coefficient_map(vector):
    """Five by four map U -> Lambda^2_0 U*: omit dependent (3,4)."""
    matrices = form_matrices(vector)
    return s.Matrix([[matrices[k][i, j] for k in range(4)]
                     for i, j in PAIRS if (i, j) != (2, 3)])


def sp_basis():
    basis = []
    for i in range(4):
        for j in range(i, 4):
            symmetric = s.zeros(4)
            symmetric[i, j] = symmetric[j, i] = 1
            a = J * symmetric
            assert a.T * J + J * a == s.zeros(4)
            basis.append(a)
    return basis


def infinitesimal_action(vector, a):
    """Lie derivative of beta under y -> exp(t A)y (pullback)."""
    b = form_matrices(vector)
    result = []
    for k in range(4):
        value = a.T * b[k] + b[k] * a
        for i in range(4):
            value += a[i, k] * b[i]
        result.append(value)
    return flatten(result)


def minor_witness(matrix):
    columns = list(matrix.rref()[1])
    rows = list(matrix[:, columns].T.rref()[1])
    minor = matrix.extract(rows, columns)
    determinant = minor.det()
    assert determinant
    return {"rank": len(columns), "rows_zero_based": rows,
            "columns_zero_based": columns, "determinant": str(determinant),
            "matrix": [[str(x) for x in row] for row in minor.tolist()]}


def verify():
    d, trace = constraints()
    constraint = d.col_join(trace)
    kernel = s.Matrix.hstack(*constraint.nullspace())
    assert (d.rank(), trace.rank(), constraint.rank(), kernel.shape[1]) == (4, 4, 8, 16)
    basis = sp_basis()
    for a in basis:
        for col in range(16):
            assert constraint * infinitesimal_action(kernel[:, col], a) == s.zeros(8, 1)

    # A fixed integer witness, with no numerical or finite-field rank inference.
    parameters = s.Matrix(list(range(1, 17)))
    generic = kernel * parameters
    action = s.Matrix.hstack(*(infinitesimal_action(generic, a) for a in basis))
    assert action.rank() == 10
    assert action.row_join(generic).rank() == 11
    assert coefficient_map(generic).rank() == 4

    # Highest-weight/rank-one representative beta = y1 dy1 wedge dy3.
    rank_one = s.zeros(24, 1)
    rank_one[COORDS.index(((0, 2), 0))] = 1
    assert constraint * rank_one == s.zeros(8, 1)
    assert coefficient_map(rank_one).rank() == 1
    rank_one_action = s.Matrix.hstack(*(infinitesimal_action(rank_one, a) for a in basis))
    assert rank_one_action.rank() == 5

    # Raise indices. The representative extends with no higher affine terms.
    y = s.symbols("y1:5")
    b = sum((y[k] * b for k, b in enumerate(form_matrices(rank_one))), s.zeros(4))
    linear = -J * b * J
    poisson = J + linear
    jacobi = []
    for i, j, k in TRIPLES:
        expr = sum(poisson[i, h] * s.diff(poisson[j, k], y[h])
                   + poisson[j, h] * s.diff(poisson[k, i], y[h])
                   + poisson[k, h] * s.diff(poisson[i, j], y[h])
                   for h in range(4))
        jacobi.append(s.expand(expr))
    assert jacobi == [0] * 4
    assert s.expand(poisson[0, 1]*poisson[2, 3]
                    - poisson[0, 2]*poisson[1, 3]
                    + poisson[0, 3]*poisson[1, 2]) == 1

    return {
        "all_passed": True,
        "field": "QQ",
        "coordinate_order": [{"pair_one_based": [i+1, j+1], "linear_variable_one_based": k+1}
                             for (i, j), k in COORDS],
        "closure_rank": 4, "combined_closure_primitive_rank": 8,
        "primitive_closed_dimension": 16,
        "kernel_basis_columns": [[str(x) for x in row] for row in kernel.tolist()],
        "sp4_lie_basis": [[[str(x) for x in row] for row in a.tolist()] for a in basis],
        "all_ten_generators_preserve_all_sixteen_basis_vectors": True,
        "generic_witness_parameters": list(range(1, 17)),
        "generic_witness_coefficients": [str(x) for x in generic],
        "generic_coefficient_rank_witness": minor_witness(coefficient_map(generic)),
        "generic_sp4_orbit_rank_witness": minor_witness(action),
        "generic_sp4_times_scalars_orbit_rank_witness": minor_witness(action.row_join(generic)),
        "rank_one_coefficients": [str(x) for x in rank_one],
        "rank_one_sp4_orbit_rank_witness": minor_witness(rank_one_action),
        "rank_one_poisson_extension": "J + y1*d2 wedge d4",
        "rank_one_jacobi_coefficients": [str(x) for x in jacobi],
        "rank_one_affine_pfaffian": "1",
        "scope": "First-jet representation and rank-one stratum, not S61 exhaustion. "
                 "The generic rank witness is not claimed to extend to a Poisson structure.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = verify()
    payload = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        # A new checkpoint should never silently replace old evidence.
        with args.output.open("x") as out:
            out.write(payload)
    print(json.dumps({key: report[key] for key in (
        "all_passed", "field", "primitive_closed_dimension", "rank_one_poisson_extension", "scope")}, indent=2))
    if args.output:
        print(args.output)


if __name__ == "__main__":
    main()

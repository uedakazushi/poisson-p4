#!/usr/bin/env python3
"""Specialize the full S61 to a nonzero primitive rank-one first jet.

Retains all 45 weight-two coordinates. Constant row operations expose linear
consequences; no Groebner decomposition is attempted. Inputs stay immutable.
"""
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import argparse

import sympy as s
from s61_model import S61, X, J
from weight1_structure import minor_witness


def clean(poly):
    return {mon: c for mon, c in poly.items() if c}


def add(poly, other, scalar):
    for mon, c in other.items():
        value = poly.get(mon, Fraction(0)) + scalar*c
        if value:
            poly[mon] = value
        else:
            poly.pop(mon, None)


def row_reduce(polynomials):
    basis = {}
    for original in polynomials:
        row = dict(original)
        while row:
            pivot = min(row, key=lambda mon: (-len(mon), mon))
            if pivot not in basis:
                leading = row[pivot]
                basis[pivot] = {mon: c/leading for mon, c in row.items()}
                break
            add(row, basis[pivot], -row[pivot])
    return basis


def encode(poly):
    return [[list(mon), c.numerator, c.denominator]
            for mon, c in sorted(poly.items())]


def run():
    model = S61()
    bivector = X[0]**2*J
    bivector[2, 4] += X[0]*X[1]
    bivector[4, 2] -= X[0]*X[1]
    normalized = model.normalize(bivector)
    model.verify_point(normalized)
    weights = model.reduced["parameter_weights"]
    weight_one = [i for i, w in enumerate(weights) if w == 1]
    weight_two = [i for i, w in enumerate(weights) if w == 2]
    variable_index = {old: new for new, old in enumerate(weight_two)}
    equations = []
    for equation in model.reduced["equations"]:
        result = {}
        for monomial, numerator, denominator in equation:
            coefficient = Fraction(numerator, denominator)
            free = []
            for i in monomial:
                if i in variable_index:
                    free.append(variable_index[i])
                else:
                    coefficient *= Fraction(normalized["parameters"][i])
            free = tuple(free)
            result[free] = result.get(free, Fraction(0)) + coefficient
        equations.append(clean(result))
    basis = row_reduce(equations)
    linear = [row for pivot, row in basis.items() if len(pivot) <= 1]
    assert all(max(map(len, row), default=0) <= 1 for row in linear)
    linear_basis = row_reduce(linear)
    assert () not in linear_basis, "Empty fiber would contradict the exact Poisson point"
    assert all(normalized["parameters"][i] == 0 for i in weight_two)
    jacobian = s.zeros(134, 45)
    for row, poly in enumerate(equations):
        assert not poly.get((), 0)
        for mon, coefficient in poly.items():
            if len(mon) == 1:
                jacobian[row, mon[0]] = s.Rational(coefficient.numerator, coefficient.denominator)
    assert jacobian.rank() == 7
    # Keep the full specialized equations and all exact linear consequences.
    # No nonlinear consequence or radical claim is inferred.
    source_paths = [model.source/name for name in (
        "symplectic_slice_certificate.json", "symplectic_slice_equations.json",
        "symplectic_slice_reduced_equations.json")]
    return {
        "all_passed": True, "field": "QQ", "parameter_count": 45,
        "base_bivector": "x0^2*(d1^d2+d3^d4)+x0*x1*d2^d4",
        "normalized_weight_one": {f"u{i+1}": str(normalized["parameters"][i]) for i in weight_one},
        "known_point_weight_two": [str(normalized["parameters"][i]) for i in weight_two],
        "perspective_parameters": [str(x) for x in normalized["perspective"]],
        "variable_names": [f"u{i+1}" for i in weight_two],
        "equations": [encode(poly) for poly in equations],
        "degree_distribution_including_zero": dict(sorted(Counter(
            max(map(len, poly), default=-1) for poly in equations).items())),
        "constant_linear_span_rank": len(basis),
        "linear_consequence_rank": len(linear_basis),
        "linear_consequences": [encode(row) for row in linear_basis.values()],
        "jacobian_at_known_point_rank_witness": minor_witness(jacobian),
        "fiber_zariski_tangent_dimension_at_known_point": 38,
        "constant_row_basis": [encode(row) for row in basis.values()],
        "sources": [{"path": str(path.relative_to(Path(__file__).resolve().parents[1])),
                     "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths],
        "scope": "Exact fiber in the archived S61 gauge at a nonzero primitive rank-one first jet. "
                 "No radical, component, dimension, or exhaustion assertion.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    report = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as out:
        json.dump(report, out, indent=2)
        out.write("\n")
    print(json.dumps({k: report[k] for k in (
        "all_passed", "normalized_weight_one", "perspective_parameters",
        "degree_distribution_including_zero", "constant_linear_span_rank", "linear_consequence_rank")}, indent=2))


if __name__ == "__main__":
    main()

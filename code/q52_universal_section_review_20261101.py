#!/usr/bin/env python3
"""Independent unrestricted-coefficient review of the universal q52 section.

The original solver starts from the 45 reduced weight-two coordinates. This
review instead solves for all 60 affine quadratic coefficients using the graph
condition and the degree-one Jacobi identity. The saved section is read only
after that independent solve. No implementation from the original solver is
imported. Arithmetic is over QQ, with all 16 first-jet coordinates symbolic.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path
import time

import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

from s61_model import S61, X

ROOT = Path(__file__).resolve().parents[1]
SOURCE_REPORT = ROOT / "reports/Q52_UNIVERSAL_SECTION_20261101.json"
J = s.Matrix([[0, 1, 0, 0], [-1, 0, 0, 0],
              [0, 0, 0, 1], [0, 0, -1, 0]])
PAIRS = list(combinations(range(4), 2))


def star(p):
    """Pfaffian adjugate: star(p) = -Pf(p)*p^{-1}."""
    return s.Matrix([[0, p[2, 3], -p[1, 3], p[1, 2]],
                     [-p[2, 3], 0, p[0, 3], -p[0, 2]],
                     [p[1, 3], -p[0, 3], 0, p[0, 1]],
                     [-p[1, 2], p[0, 2], -p[0, 1], 0]])


def affine(b, y):
    sub = dict(zip(X, [1, *y], strict=True))
    return s.Matrix(4, 4, lambda i, j: s.expand(
        (b[i + 1, j + 1] - X[i + 1] * b[0, j + 1]
         + X[j + 1] * b[0, i + 1]).subs(sub)))


def part(f, y, degree):
    return sum(c * s.prod(v**e for v, e in zip(y, mon, strict=True))
               for mon, c in s.Poly(f, *y).terms() if sum(mon) == degree)


def matrix_strings(a):
    return [[str(v) for v in row] for row in a.tolist()]


def verify():
    started = time.monotonic()
    model = S61()
    u = s.symbols("u1:62")
    weights = model.reduced["parameter_weights"]
    one = [i for i, w in enumerate(weights) if w == 1]
    two = [i for i, w in enumerate(weights) if w == 2]
    y = s.symbols("y0:4")
    base = model.reconstruct([u[i] if i in one else 0 for i in range(61)])
    linear = affine(base, y).applyfunc(lambda f: part(f, y, 1))
    beta = -J * linear * J
    assert all(s.expand(s.diff(beta[i, j], y[k])
                        + s.diff(beta[j, k], y[i])
                        + s.diff(beta[k, i], y[j])) == 0
               for i, j, k in combinations(range(4), 3))
    jet_matrix = s.Matrix([s.diff(beta[i, j], v) for i, j in PAIRS for v in y])
    assert jet_matrix.jacobian([u[i] for i in one]).rank() == 16

    # Sixty unrestricted quadratic coefficients, independent of S61's
    # previously eliminated lower Jacobi equations.
    monomials = [y[a] * y[b]
                 for a, b in combinations_with_replacement(range(4), 2)]
    unknowns = s.symbols("z0:60")
    quadratic = s.zeros(4)
    for h, (i, j) in enumerate(PAIRS):
        quadratic[i, j] = sum(unknowns[10 * h + k] * v
                              for k, v in enumerate(monomials))
        quadratic[j, i] = -quadratic[i, j]

    r = s.symbols("r0:5")
    z = s.Matrix([[0, r[0], r[1], r[2]],
                  [-r[0], 0, r[3], r[4]],
                  [-r[1], -r[3], 0, -r[0]],
                  [-r[2], -r[4], r[0], 0]])
    ell = s.Matrix([sum(s.diff(beta[i, j], y[k]) * z[i, j]
                        for i, j in PAIRS) for k in range(4)])
    w = -J * ell
    yy = s.Matrix(y)
    linear_graph = yy * w.T - w * yy.T
    graph_expression = s.expand(sum(
        star(quadratic)[i, j] * z[i, j]
        + star(linear)[i, j] * linear_graph[i, j] for i, j in PAIRS))
    graph_poly = s.Poly(graph_expression, *y, *r)
    graph_equations = graph_poly.coeffs()
    graph_a, graph_b = s.linear_eq_to_matrix(graph_equations, unknowns)
    assert graph_a.shape == (50, 60) and not graph_a.free_symbols
    assert graph_a.rank() == 50

    # Degree-one part of [J+L+Q+C,J+L+Q+C]; the cubic term cannot enter.
    jacobi_equations = []
    for i, j, k in combinations(range(4), 3):
        expression = sum(
            J[i, a] * s.diff(quadratic[j, k], y[a])
            + J[j, a] * s.diff(quadratic[k, i], y[a])
            + J[k, a] * s.diff(quadratic[i, j], y[a])
            + linear[i, a] * s.diff(linear[j, k], y[a])
            + linear[j, a] * s.diff(linear[k, i], y[a])
            + linear[k, a] * s.diff(linear[i, j], y[a])
            for a in range(4))
        jacobi_equations.extend(s.Poly(s.expand(expression), *y).coeffs())
    a, b = s.linear_eq_to_matrix(graph_equations + jacobi_equations, unknowns)
    assert a.shape == (66, 60) and not a.free_symbols
    rows = list(a.T.rref()[1])
    assert len(rows) == 60
    minor = a[rows, :]
    solution = (minor.inv() * b[rows, 0]).applyfunc(s.expand)
    assert (a * solution - b).applyfunc(s.expand) == s.zeros(66, 1)
    assert all(s.Poly(f, *[u[i] for i in one]).total_degree() == 2
               for f in solution)

    # Original 45-coordinate linear matrix is independently recovered by
    # restricting the graph map to the homogeneous lower-Jacobi kernel.
    kernel_affine = affine(model.reconstruct(
        [u[i] if i in two else 0 for i in range(61)]), y)
    kernel_coefficients = s.Matrix([
        s.Poly(kernel_affine[i, j], *y).coeff_monomial(mon)
        for i, j in PAIRS for mon in monomials])
    embedding = kernel_coefficients.jacobian([u[i] for i in two])
    restricted = graph_a * embedding
    restricted_rows = list(restricted.T.rref()[1])
    assert embedding.shape == (60, 45) and embedding.rank() == 45
    assert restricted.shape == (50, 45) and len(restricted_rows) == 45

    # Only now read and compare the original section.
    original = json.loads(SOURCE_REPORT.read_text())
    assert original["weight1_u_indices_zero_based"] == one
    assert original["weight2_u_indices_zero_based"] == two
    qvalues = [s.sympify(f) for f in original["quadratic_section_q"]]
    values = [None] * 61
    for i in one:
        values[i] = u[i]
    for i, f in zip(two, qvalues, strict=True):
        values[i] = f
    full = affine(model.reconstruct(values), y)
    independent_q = quadratic.subs(dict(zip(unknowns, solution, strict=True)))
    assert all(s.expand(part(full[i, j], y, 2) - independent_q[i, j]) == 0
               for i, j in PAIRS)
    full_degree = max(s.Poly(full[i, j], *y).total_degree() for i, j in PAIRS)
    assert full_degree == 3
    annihilation = s.expand(sum(star(full)[i, j]
                               * (z + linear_graph)[i, j] for i, j in PAIRS))
    assert annihilation == 0

    # Direct full affine Jacobi, independently of all 134 reduced equations.
    polys, *gens = ring(",".join(map(str, [u[i] for i in one] + list(y))), QQ)
    p = {(i, j): polys.from_expr(full[i, j])
         for i in range(4) for j in range(4)}
    ygens = gens[-4:]
    for i, j, k in combinations(range(4), 3):
        jacobi = sum((p[i, a] * p[j, k].diff(ygens[a])
                      + p[j, a] * p[k, i].diff(ygens[a])
                      + p[k, a] * p[i, j].diff(ygens[a])
                      for a in range(4)), polys.zero)
        assert not jacobi
    polyvalues = [polys.from_expr(f) for f in values]
    for equation in model.reduced["equations"]:
        result = polys.zero
        for mon, numerator, denominator in equation:
            term = polys(QQ(numerator, denominator))
            for i in mon:
                term *= polyvalues[i]
            result += term
        assert not result

    paths = [Path(__file__), ROOT / "code/s61_model.py", SOURCE_REPORT]
    paths += [model.source / filename for filename in (
        "symplectic_slice_certificate.json", "symplectic_slice_equations.json",
        "symplectic_slice_reduced_equations.json")]
    return {
        "id": "Q52_UNIVERSAL_SECTION_REVIEW_20261101",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "all_passed": True,
        "field": "QQ; all 16 first-jet variables symbolic",
        "independence": "60 unrestricted quadratic coefficients, then graph and Jacobi",
        "closed_first_jet_rank": 16,
        "graph_matrix_shape": list(graph_a.shape), "graph_matrix_rank": 50,
        "graph_matrix": matrix_strings(graph_a),
        "graph_rhs": list(map(str, graph_b)),
        "combined_matrix_shape": list(a.shape), "combined_matrix_rank": 60,
        "combined_matrix": matrix_strings(a), "combined_rhs": list(map(str, b)),
        "combined_pivot_rows": rows, "combined_pivot_determinant": str(minor.det()),
        "independent_affine_quadratic_coefficients": list(map(str, solution)),
        "recovered_50_by_45_matrix": matrix_strings(restricted),
        "recovered_50_by_45_rank": 45,
        "recovered_pivot_rows": restricted_rows,
        "recovered_pivot_determinant": str(restricted[restricted_rows, :].det()),
        "all_60_quadratic_coefficients_match": True,
        "full_affine_bivector_degree": full_degree,
        "full_graph_annihilation": True, "direct_affine_jacobi": True,
        "all_134_reduced_equations": True,
        "proof_note": "research/notes/Q52_UNIVERSAL_SECTION_REVIEW_20261101.md",
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "input_sha256": {str(path.relative_to(ROOT)):
                         hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = verify()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
    print(json.dumps({"all_passed": report["all_passed"],
                      "elapsed_seconds": report["elapsed_seconds"],
                      "output": str(args.output)}))


if __name__ == "__main__":
    main()

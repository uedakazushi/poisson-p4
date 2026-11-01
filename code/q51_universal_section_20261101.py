#!/usr/bin/env python3
"""Universal quadratic q(5,1) section and the three-covariant section ansatz.

This proves polynomial identities for all 16 first-jet variables. The finite
Groebner calculation only determines constant coefficients in a covariant
ansatz; it is not an elimination of the full S61 fiber.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
from itertools import combinations, combinations_with_replacement, permutations, product
import json
from pathlib import Path
import time

import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

from s61_model import S61, X
from weight1_structure import COORDS, constraints

ROOT = Path(__file__).resolve().parents[1]
J = s.Matrix([[0, 1, 0, 0], [-1, 0, 0, 0],
              [0, 0, 0, 1], [0, 0, -1, 0]])
PAIRS = list(combinations(range(4), 2))
Q52_REPORT = ROOT / "reports/Q52_UNIVERSAL_SECTION_20261101.json"
WITNESS_REPORT = ROOT / "reports/S61_KNOWN_FIBER_AUDIT_20261031.json"


def affine(b, y):
    substitution = dict(zip(X, [1, *y], strict=True))
    return s.Matrix(4, 4, lambda i, j: s.expand(
        (b[i + 1, j + 1] - X[i + 1] * b[0, j + 1]
         + X[j + 1] * b[0, i + 1]).subs(substitution)))


def exterior_derivative(alpha, y):
    return s.Matrix(4, 4, lambda i, j:
                    s.expand(s.diff(alpha[j], y[i]) - s.diff(alpha[i], y[j])))


def evaluate_polynomial(encoded, values, polynomial_ring):
    answer = polynomial_ring.zero
    for monomial, numerator, denominator in encoded:
        term = polynomial_ring(QQ(numerator, denominator))
        for i in monomial:
            term *= values[i]
        answer += term
    return answer


def character_certificate():
    """Exact characters using the closed-form Koszul sequence and Weyl alternant."""
    standard = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    kernel = s.Matrix.vstack(*constraints()).nullspace()
    w_weights = []
    for column in kernel:
        weights = set()
        for coefficient, ((i, j), k) in zip(column, COORDS, strict=True):
            if coefficient:
                weights.add(tuple(standard[i][h] + standard[j][h] + standard[k][h]
                                  for h in range(2)))
        assert len(weights) == 1
        w_weights.append(weights.pop())
    assert len(w_weights) == 16

    def add(v, w):
        return v[0] + w[0], v[1] + w[1]

    def multiply(left, right):
        result = Counter()
        for v, a in left.items():
            for w, b in right.items():
                result[add(v, w)] += a * b
        return result

    sym_w = Counter(add(v, w) for i, v in enumerate(w_weights) for w in w_weights[i:])
    sym_v = Counter(add(v, w) for i, v in enumerate(standard) for w in standard[i:])
    wedge2 = Counter(add(standard[i], standard[j]) for i, j in combinations(range(4), 2))
    wedge3 = Counter(tuple(sum(standard[i][h] for i in indices) for h in range(2))
                     for indices in combinations(range(4), 3))
    # 0 -> K45 -> Sym²(V*)⊗Λ²(V*) -> V*⊗Λ³(V*) -> Λ⁴(V*) -> 0.
    k_weights = multiply(sym_v, wedge2)
    for weight, multiplicity in multiply(Counter(standard), wedge3).items():
        k_weights[weight] -= multiplicity
    k_weights[0, 0] += 1
    assert min(k_weights.values()) >= 0 and sum(k_weights.values()) == 45
    rho = (2, 1)

    def decomposition(character):
        result = []
        for first in range(9):
            for second in range(first + 1):
                value = 0
                for perm in permutations((0, 1)):
                    for signs in product((-1, 1), repeat=2):
                        wr = tuple(signs[i] * rho[perm[i]] for i in range(2))
                        sign = (1 if perm == (0, 1) else -1) * signs[0] * signs[1]
                        value += sign * character[first + rho[0] - wr[0],
                                                  second + rho[1] - wr[1]]
                if value:
                    assert value > 0
                    a, b = first - second, second
                    dimension = (a + 1) * (b + 1) * (a + b + 2) * (a + 2*b + 3) // 6
                    result.append({"Dynkin_highest_weight": [a, b],
                                   "multiplicity": value, "dimension": dimension})
        return result

    sym_decomposition = decomposition(sym_w)
    k_decomposition = decomposition(k_weights)
    assert sym_decomposition == [
        {"Dynkin_highest_weight": [2, 0], "multiplicity": 2, "dimension": 10},
        {"Dynkin_highest_weight": [2, 1], "multiplicity": 1, "dimension": 35},
        {"Dynkin_highest_weight": [2, 2], "multiplicity": 1, "dimension": 81}]
    assert k_decomposition == [
        {"Dynkin_highest_weight": [2, 0], "multiplicity": 1, "dimension": 10},
        {"Dynkin_highest_weight": [2, 1], "multiplicity": 1, "dimension": 35}]
    return {"W_weights_with_multiplicity": [list(v) for v in w_weights],
            "Sym2W_decomposition": sym_decomposition, "K45_decomposition": k_decomposition,
            "equivariant_quadratic_map_dimension": 3}


def verify():
    started = time.monotonic()
    model = S61()
    u = s.symbols("u1:62")
    one = [i for i, w in enumerate(model.reduced["parameter_weights"]) if w == 1]
    two = [i for i, w in enumerate(model.reduced["parameter_weights"]) if w == 2]
    y = s.Matrix(s.symbols("y0:4"))
    base = affine(model.reconstruct([u[i] if i in one else 0 for i in range(61)]), y)
    zero = dict.fromkeys(y, 0)
    linear = sum((base.diff(v).subs(zero) * v for v in y), s.zeros(4))
    alpha = s.Matrix([s.diff(linear[0, 1] + linear[2, 3], v)/5 for v in y])
    vector = J * alpha
    primitive = linear - 2*alpha.dot(y)*J + y*vector.T - vector*y.T
    beta = (-J * primitive * J).applyfunc(s.expand)
    assert s.expand(beta[0, 1] + beta[2, 3]) == 0

    # Three closed quadratic two-forms, homogeneous of degree two in beta.
    xfield = -J * beta * y
    dx = xfield.jacobian(y)
    eta1 = (dx.T * J * dx).applyfunc(s.expand)
    f = s.expand(beta[0, 1]*beta[2, 3] - beta[0, 2]*beta[1, 3]
                 + beta[0, 3]*beta[1, 2])
    eta2 = exterior_derivative(-f * J * y, y)
    nz = [(i, k, J[i, k]) for i in range(4) for k in range(4) if J[i, k]]
    second_moment = s.Matrix(4, 4, lambda a, b: sum(
        s.diff(beta[a, i], y[j])*s.diff(beta[b, k], y[l])*v*w
        for i, k, v in nz for j, l, w in nz))
    g = s.expand((y.T * second_moment * y)[0])
    eta3 = exterior_derivative(-g * J * y, y)
    for eta in (eta1, eta2, eta3):
        assert all(s.expand(s.diff(eta[i, j], y[k]) + s.diff(eta[j, k], y[i])
                            + s.diff(eta[k, i], y[j])) == 0
                   for i, j, k in combinations(range(4), 3))

    monomials = [y[i]*y[j] for i, j in combinations_with_replacement(range(4), 2)]
    quadratic_kernel = affine(model.reconstruct(
        [u[i] if i in two else 0 for i in range(61)]), y)
    embedding = s.Matrix([s.Poly(quadratic_kernel[i, j], *y).coeff_monomial(mon)
                          for i, j in PAIRS for mon in monomials]).jacobian([u[i] for i in two])
    rows = list(embedding.T.rref()[1])
    left = embedding[rows, :].inv()
    covariants = []
    for eta in (eta1, eta2, eta3):
        quadratic = -J * eta * J
        coefficients = s.Matrix([s.Poly(quadratic[i, j], *y).coeff_monomial(mon)
                                 for i, j in PAIRS for mon in monomials])
        covariant = (left * coefficients[rows, 0]).applyfunc(s.expand)
        assert (embedding * covariant - coefficients).applyfunc(s.expand) == s.zeros(60, 1)
        covariants.append(covariant)

    q52 = s.Matrix(list(map(s.sympify, json.loads(Q52_REPORT.read_text())["quadratic_section_q"])))
    sample = dict(zip([u[i] for i in one], range(1, 17), strict=True))
    independence = s.Matrix.hstack(*covariants).subs(sample)
    independence_rows = list(independence.T.rref()[1])
    assert len(independence_rows) == 3
    a, b, c = s.symbols("a b c")
    trial = (q52 + a*covariants[0] + b*covariants[1] + c*covariants[2]).subs(sample)
    ansatz_ring, *_ = ring("a,b,c", QQ)
    sample_values = [None] * 61
    for i in one:
        sample_values[i] = ansatz_ring(sample[u[i]])
    for i, expression in zip(two, trial, strict=True):
        sample_values[i] = ansatz_ring.from_expr(s.expand(expression))
    residuals = [evaluate_polynomial(eq, sample_values, ansatz_ring)
                 for eq in model.reduced["equations"]]
    groebner = s.groebner([f.as_expr() for f in residuals if f], a, b, c, order="lex")
    expected = [a - 20*c/3, b - 8*c, c**2 + 3*c/32]
    assert [poly.as_expr() for poly in groebner.polys] == expected

    constants = [-s.Rational(5, 8), -s.Rational(3, 4), -s.Rational(3, 32)]
    q51 = (q52 + sum((v*cov for v, cov in zip(constants, covariants, strict=True)),
                     s.zeros(45, 1))).applyfunc(s.expand)
    polynomial_ring, *variables = ring(",".join(str(u[i]) for i in one), QQ)
    values = [None] * 61
    for i, variable in zip(one, variables, strict=True):
        values[i] = variable
    for i, expression in zip(two, q51, strict=True):
        values[i] = polynomial_ring.from_expr(expression)
    for equation in model.reduced["equations"]:
        assert not evaluate_polynomial(equation, values, polynomial_ring)
    expression_values = [value.as_expr() for value in values]
    full = affine(model.reconstruct(expression_values), y)
    full_ring, *full_variables = ring(",".join(map(str, [u[i] for i in one] + list(y))), QQ)
    p = {(i, j): full_ring.from_expr(full[i, j]) for i in range(4) for j in range(4)}
    yvariables = full_variables[-4:]
    for i, j, k in combinations(range(4), 3):
        assert not sum((p[i, h]*p[j, k].diff(yvariables[h])
                        + p[j, h]*p[k, i].diff(yvariables[h])
                        + p[k, h]*p[i, j].diff(yvariables[h])
                        for h in range(4)), full_ring.zero)

    # Identify the component at an already known T39 point, checking its full
    # S61 fixed-jet Jacobian again. This uses no generic count inference.
    known = json.loads(WITNESS_REPORT.read_text())
    witness = next(point for point in known["points"] if point["component_id"] == "T39")
    point = list(map(s.Rational, witness["S61_coordinates_u1_through_u61"]))
    substitution = {u[i]: point[i] for i in one}
    assert all(expression.subs(substitution) == point[i]
               for i, expression in zip(two, q51, strict=True))
    assert any(expression.subs(substitution) != point[i]
               for i, expression in zip(two, q52, strict=True))
    positions = {old: new for new, old in enumerate(two)}
    jacobian = s.zeros(134, 45)
    for row, equation in enumerate(model.reduced["equations"]):
        for monomial, numerator, denominator in equation:
            for position, index in enumerate(monomial):
                if index in positions:
                    jacobian[row, positions[index]] += s.Rational(numerator, denominator) * s.prod(
                        point[v] for j, v in enumerate(monomial) if j != position)
    jacobian_rows = list(jacobian.T.rref()[1])
    assert len(jacobian_rows) == 45
    character = character_certificate()
    paths = [Path(__file__), ROOT / "code/s61_model.py", ROOT / "code/weight1_structure.py",
             Q52_REPORT, WITNESS_REPORT]
    paths += [model.source / name for name in (
        "symplectic_slice_certificate.json", "symplectic_slice_equations.json",
        "symplectic_slice_reduced_equations.json")]
    return {
        "id": "Q51_UNIVERSAL_SECTION_20261101", "all_passed": True,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "field": "QQ; all16 first-jet variables symbolic",
        "weight1_u_indices_zero_based": one, "weight2_u_indices_zero_based": two,
        "covariants_in_q_coordinates": [list(map(str, cov)) for cov in covariants],
        "covariant_independence_sample_first_jet": list(range(1, 17)),
        "covariant_independence_rows": independence_rows,
        "covariant_independence_minor": str(independence[independence_rows, :].det()),
        "character_certificate": character,
        "sample_ansatz_groebner": list(map(str, expected)),
        "q51_constants": list(map(str, constants)), "quadratic_section_q": list(map(str, q51)),
        "section_monomial_count": sum(len(value) for value in values if value not in variables),
        "universal_134_equations": True, "direct_full_affine_Jacobi": True,
        "full_affine_bivector_degree": max(s.Poly(full[i, j], *y).total_degree() for i, j in PAIRS),
        "T39_known_witness_match": True, "q52_distinct_at_witness": True,
        "full_S61_fixed_jet_Jacobian_rank_at_T39_witness": 45,
        "Jacobian_pivot_rows": jacobian_rows,
        "Jacobian_pivot_determinant": str(jacobian[jacobian_rows, :].det()),
        "proof_note": "research/notes/Q51_UNIVERSAL_SECTION_20261101.md",
        "scope": "T39 reduced S61 closure is A16 graph; the two graded Sp4-equivariant polynomial sections are q51 and q52. No full-fiber exhaustion claim.",
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "input_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in paths},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")
    print(json.dumps({"all_passed": True, "elapsed_seconds": result["elapsed_seconds"],
                      "output": str(args.output)}))


if __name__ == "__main__":
    main()

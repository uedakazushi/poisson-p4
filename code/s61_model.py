#!/usr/bin/env python3
"""Exact coordinate adapter for the immutable v02 S61 model.

Coefficient arithmetic is over QQ or a symbolic rational-function field.
This module constructs no prime decomposition. See code/S61_README.md.
"""
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
import json
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "legacy/v02/data"
X = s.symbols("x0:5")
PAIRS = list(combinations(range(5), 2))
MONOMIALS = [s.prod(X[i] for i in pair)
             for pair in combinations_with_replacement(range(5), 2)]
BASIS = list(product(PAIRS, MONOMIALS))
J = s.zeros(5)
J[1, 2] = J[3, 4] = 1
J[2, 1] = J[4, 3] = -1


def canonical(value):
    return s.cancel(value)


def coefficients(bivector):
    """150 homogeneous quadratic coefficients, v02's exact ordering."""
    if bivector.shape != (5, 5):
        raise ValueError("Expected a 5 by 5 bivector matrix")
    if any(canonical(bivector[i, j] + bivector[j, i]) != 0
           for i in range(5) for j in range(i, 5)):
        raise ValueError("Expected a skew-symmetric bivector matrix")
    polynomials = {pair: s.Poly(bivector[pair[0], pair[1]], *X)
                   for pair in PAIRS}
    if any(sum(powers) != 2 for polynomial in polynomials.values()
           for powers, coefficient in polynomial.terms() if coefficient != 0):
        raise ValueError("Expected homogeneous quadratic coefficients in x0,...,x4")
    return [polynomials[pair].coeff_monomial(mon) for pair, mon in BASIS]


def bivector_from_coefficients(values):
    out = s.zeros(5)
    for value, (pair, mon) in zip(values, BASIS, strict=True):
        out[pair[0], pair[1]] += value * mon
    for i, j in PAIRS:
        out[j, i] = -out[i, j]
    return out


@lru_cache(maxsize=1)
def vertical_rref():
    """Independent construction of E wedge End(V), without old code imports."""
    rows = []
    euler = s.Matrix(X)
    for direction, variable in product(range(5), repeat=2):
        field = s.zeros(5, 1)
        field[direction] = X[variable]
        rows.append(coefficients(euler * field.T - field * euler.T))
    reduced, pivots = s.Matrix(rows).rref()
    assert len(pivots) == 24
    return reduced[:24, :], pivots


def quotient_coefficients(bivector):
    values = coefficients(bivector)
    reduced, pivots = vertical_rref()
    return [canonical(values[k] - sum(reduced[i, k] * values[p]
                                     for i, p in enumerate(pivots)))
            for k in range(150)]


def pushforward(bivector, matrix):
    """M_*P(x)=M P(M^{-1}x) M^t; M must be invertible."""
    inverse = matrix.inv()
    substitution = dict(zip(X, inverse * s.Matrix(X), strict=True))
    out = matrix * bivector.subs(substitution, simultaneous=True) * matrix.T
    return out.applyfunc(s.expand)


def jacobi_projective(bivector):
    """Sparse coefficients of E wedge Jac(P), independent of legacy code.

    Vanishing is equivalent to vanishing of all 350 ambient coefficients.
    """
    result = []
    for quad in combinations(range(5), 4):
        expression = 0
        for pos, e in enumerate(quad):
            i, j, k = quad[:pos] + quad[pos + 1:]
            jac = sum(bivector[i, h] * s.diff(bivector[j, k], X[h])
                      + bivector[j, h] * s.diff(bivector[k, i], X[h])
                      + bivector[k, h] * s.diff(bivector[i, j], X[h])
                      for h in range(5))
            expression += (-1)**pos * X[e] * jac
        result.extend(canonical(c) for c in s.Poly(s.expand(expression), *X).coeffs())
    return result


class S61:
    def __init__(self, source=SOURCE):
        self.source = Path(source)
        self.certificate = json.loads((self.source / "symplectic_slice_certificate.json").read_text())
        self.data = json.loads((self.source / "symplectic_slice_equations.json").read_text())
        self.reduced = json.loads((self.source / "symplectic_slice_reduced_equations.json").read_text())
        assert list(vertical_rref()[1]) == self.certificate["vertical_pivots"]

    @staticmethod
    def evaluate(encoded, parameters):
        monomials = {(): s.S.One}
        total = 0
        for mon, numerator, denominator in encoded:
            key = tuple(mon)
            if key not in monomials:
                monomials[key] = s.prod(parameters[i] for i in mon)
            total += s.Rational(numerator, denominator) * monomials[key]
        return canonical(total)

    def reconstruct(self, parameters):
        if len(parameters) != 61:
            raise ValueError("S61 requires exactly 61 parameters")
        values = coefficients(X[0]**2 * J)
        for index, value in zip(self.certificate["raw_free_indices"], parameters, strict=True):
            values[index] = value
        for entry in self.data["eliminated_raw_variables"]:
            values[entry["raw_index"]] = self.evaluate(entry["polynomial"], parameters)
        return bivector_from_coefficients(values)

    @lru_cache(maxsize=1)
    def perspective_matrix(self):
        a = s.symbols("a1:5")
        matrix = s.eye(5)
        for j in range(1, 5):
            matrix[0, j] = a[j - 1]
        values = quotient_coefficients(pushforward(X[0]**2 * J, matrix))
        selected = [values[k] for k in self.certificate["perspective_orbit_pivots_raw"]]
        response = s.Matrix(selected).jacobian(a)
        assert not response.free_symbols and response.det() != 0
        return response

    def normalize(self, bivector):
        """Remove perspective coordinates when P([e0])=J; return all data.

        This fixes the evaluation and projective scale in advance. It does not
        choose a rank-four point or normalize a general value automatically.
        No Poisson hypothesis is needed for the perspective operation itself.
        """
        initial = quotient_coefficients(bivector)
        base = coefficients(X[0]**2 * J)
        for index in self.certificate["constant_raw_indices"]:
            if canonical(initial[index] - base[index]) != 0:
                raise ValueError("First fix the value at [e0] to d12+d34")
        selected = s.Matrix([initial[k] for k in self.certificate["perspective_orbit_pivots_raw"]])
        parameters = -self.perspective_matrix().inv() * selected
        matrix = s.eye(5)
        for j in range(1, 5):
            matrix[0, j] = canonical(parameters[j - 1])
        normalized = quotient_coefficients(pushforward(bivector, matrix))
        assert all(normalized[k] == 0 for k in self.certificate["perspective_orbit_pivots_raw"])
        return {
            "parameters": [normalized[k] for k in self.certificate["raw_free_indices"]],
            "perspective": list(parameters),
            "matrix": matrix,
            "coefficients": normalized,
        }

    def verify_point(self, normalized, check_equations=True):
        parameters = normalized["parameters"]
        reconstruction = coefficients(self.reconstruct(parameters))
        assert all(canonical(a - b) == 0 for a, b in
                   zip(reconstruction, normalized["coefficients"], strict=True))
        if check_equations:
            self.verify_equations(parameters)
        assert all(v == 0 for v in jacobi_projective(bivector_from_coefficients(reconstruction)))
        return {"reconstruction_coefficients": 150,
                "reduced_equations": 134 if check_equations else None,
                "projective_jacobi": True}

    def verify_equations(self, parameters):
        """Exact QQ-polynomial verification after a common denominator clears.

        The common denominator is nonzero in the coefficient field. We do not
        infer anything at its zero locus. Sparse polynomial rings avoid costly
        repeated simplification of a large rational-expression sum.
        """
        symbols = sorted(set().union(*(s.sympify(v).free_symbols for v in parameters)), key=str)
        if not symbols:
            assert all(self.evaluate(eq, parameters) == 0 for eq in self.reduced["equations"])
            return
        denominator = s.S.One
        for value in parameters:
            denominator = s.lcm(denominator, s.denom(canonical(value)))
        ring = s.polys.rings.ring(symbols, s.QQ)[0]
        denominator_poly = ring.from_expr(denominator)
        numerators = [ring.from_expr(canonical(value*denominator)) for value in parameters]

        @lru_cache(maxsize=4096)
        def monomial(indices):
            if not indices:
                return ring.one
            if any(not numerators[i] for i in indices):
                return ring.zero
            return numerators[indices[0]]*monomial(indices[1:])

        powers = [denominator_poly**k for k in range(6)]
        for equation_index, equation in enumerate(self.reduced["equations"]):
            degree = max(len(m) for m, _, _ in equation)
            total = ring.zero
            for indices, num, den in equation:
                term = monomial(tuple(indices))
                if term:
                    total += ring.domain.convert(s.Rational(num, den))*term*powers[degree-len(indices)]
            assert not total, f"Nonzero reduced equation {equation_index}"


def toric_family(a, b, c, d, e, f):
    """Full six-parameter toric family on a*(a*f-b*e+c*d) != 0.

    Move (1:1:1:1:1) to [e0], then apply explicit symplectic row operations.
    The returned bivector has value J. All six original parameters vary.
    """
    a, b, c, d, e, f = map(s.sympify, (a, b, c, d, e, f))
    delta = a*f - b*e + c*d
    if a == 0 or delta == 0:
        raise ValueError("This chart requires a != 0 and Pfaffian != 0")
    out = s.zeros(5)
    for (i, j), value in zip(combinations(range(1, 5), 2), (a, b, c, d, e, f), strict=True):
        out[i, j] = value*(X[0] + X[i])*(X[0] + X[j])
        out[j, i] = -out[i, j]
    frame = s.eye(5)
    frame[1:5, 1:5] = s.Matrix([
        [1, 0, 0, 0], [0, 1/a, 0, 0],
        [d/a, -b/a, 1, 0], [e/delta, -c/delta, 0, a/delta],
    ])
    return pushforward(out, frame), frame, delta

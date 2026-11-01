#!/usr/bin/env python3
"""Exact intersection of the two universal elliptic sections.

The input Groebner candidate is NOT trusted: its degree-two span, every
nontrivial Buchberger pair and its monomial Hilbert series are verified here.
No Macaulay2 runtime is required to verify this certificate.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
from itertools import combinations
import json
from pathlib import Path
import time

import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

from s61_model import S61
from q51_universal_section_20261101 import J, affine

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "reports/ELLIPTIC_SECTION_INTERSECTION_GB_INPUT_20261101.json"
Q51 = ROOT / "reports/Q51_UNIVERSAL_SECTION_20261101.json"
Q52 = ROOT / "reports/Q52_UNIVERSAL_SECTION_20261101.json"


def echelon(polynomials):
    pivots = {}
    for polynomial in polynomials:
        row = dict(polynomial)
        while row:
            pivot = max(row)
            coefficient = row[pivot]
            if pivot not in pivots:
                pivots[pivot] = {key: value/coefficient for key, value in row.items()}
                break
            for key, value in pivots[pivot].items():
                updated = row.get(key, QQ.zero) - coefficient*value
                if updated:
                    row[key] = updated
                elif key in row:
                    del row[key]
    return pivots


def remainder(polynomial, pivots):
    row = dict(polynomial)
    while row:
        pivot = max(row)
        if pivot not in pivots:
            return row
        coefficient = row[pivot]
        for key, value in pivots[pivot].items():
            updated = row.get(key, QQ.zero) - coefficient*value
            if updated:
                row[key] = updated
            elif key in row:
                del row[key]
    return row


def minimal_monomials(generators):
    answer = []
    for monomial in sorted(set(generators), key=lambda v: (sum(v), v)):
        if not any(all(a <= b for a, b in zip(old, monomial, strict=True)) for old in answer):
            answer.append(monomial)
    return tuple(answer)


def polynomial_add(left, right):
    result = [0] * max(len(left), len(right))
    for i, value in enumerate(left):
        result[i] += value
    for i, value in enumerate(right):
        result[i] += value
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


@lru_cache(maxsize=None)
def monomial_hilbert_numerator(generators, variable_count):
    """Return numerator with denominator (1-t)^variable_count.

    Uses 0 -> R/(I:x)(-1) -> R/I -> R/(I,x) -> 0. All arithmetic
    is integer arithmetic on coefficient tuples, not sampled Hilbert values.
    """
    if not generators:
        return (1,)
    if sum(generators[0]) == 0:
        return (0,)
    variable = max(range(variable_count),
                   key=lambda i: sum(monomial[i] > 0 for monomial in generators))
    eliminated = minimal_monomials(
        monomial[:variable] + monomial[variable + 1:]
        for monomial in generators if monomial[variable] == 0)
    colon = minimal_monomials(
        tuple(value - (index == variable and value > 0)
              for index, value in enumerate(monomial)) for monomial in generators)
    left = monomial_hilbert_numerator(eliminated, variable_count - 1)
    right = monomial_hilbert_numerator(colon, variable_count)
    return polynomial_add(polynomial_add(left, (0,) + tuple(-v for v in left)),
                          (0,) + right)


def verify():
    started = time.monotonic()
    model = S61()
    u = s.symbols("u1:62")
    one = [i for i, weight in enumerate(model.reduced["parameter_weights"]) if weight == 1]
    variables = [u[i] for i in one]
    polynomial_ring, *generators = ring(",".join(map(str, variables)), QQ, order="grevlex")
    y = s.Matrix(s.symbols("y0:4"))
    base = affine(model.reconstruct([u[i] if i in one else 0 for i in range(61)]), y)
    zero = dict.fromkeys(y, 0)
    linear = sum((base.diff(v).subs(zero)*v for v in y), s.zeros(4))
    alpha = s.Matrix([s.diff(linear[0, 1] + linear[2, 3], v)/5 for v in y])
    vector = J*alpha
    beta = (-J*(linear - 2*alpha.dot(y)*J + y*vector.T - vector*y.T)*J).applyfunc(s.expand)
    pairs = list(combinations(range(4), 2))[:-1]
    coefficient_matrix = s.Matrix([[s.diff(beta[i, j], v) for v in y] for i, j in pairs])
    rank_one_generators = [polynomial_ring.from_expr(
        s.expand(coefficient_matrix.extract(rows, columns).det()))
        for rows in combinations(range(5), 2) for columns in combinations(range(4), 2)]
    q51 = json.loads(Q51.read_text())
    q52 = json.loads(Q52.read_text())
    difference_generators = [polynomial_ring.from_expr(s.sympify(a) - s.sympify(b))
                             for a, b in zip(q51["quadratic_section_q"],
                                             q52["quadratic_section_q"], strict=True)]
    i2 = echelon(difference_generators)
    j2 = echelon(rank_one_generators)
    assert len(i2) == 45 and len(j2) == 55
    assert all(not remainder(f, j2) for f in difference_generators)
    i3 = echelon(x*f for x in generators for f in difference_generators)
    assert len(i3) == 560
    assert all(not remainder(x*f, i3) for x in generators for f in rank_one_generators)

    candidate = json.loads(INPUT.read_text())
    assert candidate["variables"] == list(map(str, variables))
    groebner = [polynomial_ring.from_expr(s.sympify(f.replace("^", "**")))
                for f in candidate["groebner"]]
    assert len(groebner) == 55 and all(set(map(sum, f)) == {2} for f in groebner)
    g2 = echelon(groebner)
    assert len(g2) == 55
    assert all(not remainder(f, j2) for f in groebner)
    assert all(not remainder(f, g2) for f in rank_one_generators)
    groebner = [f.monic() for f in groebner]
    checked_pairs = 0
    coprime_pairs = 0
    for f, g in combinations(groebner, 2):
        lm_f, lm_g = f.LM, g.LM
        if all(min(a, b) == 0 for a, b in zip(lm_f, lm_g, strict=True)):
            coprime_pairs += 1
            continue
        lcm = tuple(max(a, b) for a, b in zip(lm_f, lm_g, strict=True))
        spoly = (f.mul_monom(tuple(a-b for a, b in zip(lcm, lm_f, strict=True)))
                 - g.mul_monom(tuple(a-b for a, b in zip(lcm, lm_g, strict=True))))
        assert not spoly.rem(groebner)
        checked_pairs += 1

    leading = minimal_monomials(f.LM for f in groebner)
    numerator = monomial_hilbert_numerator(leading, 16)
    t = s.Symbol("t")
    expected = s.Poly((1+11*t+11*t**2+t**3)*(1-t)**11, t)
    assert list(numerator) == [expected.nth(i) for i in range(expected.degree()+1)]
    # The nilradical quotient is killed by m; its degree2 dimension is 55-45.
    intersection_numerator = s.expand(expected.as_expr() + 10*t**2*(1-t)**16)
    paths = [Path(__file__), INPUT, Q51, Q52, ROOT / "code/s61_model.py",
             ROOT / "code/q51_universal_section_20261101.py"]
    paths += [model.source / name for name in (
        "symplectic_slice_certificate.json", "symplectic_slice_equations.json",
        "symplectic_slice_reduced_equations.json")]
    return {
        "id": "ELLIPTIC_SECTION_INTERSECTION_20261101", "all_passed": True,
        "generated_utc": datetime.now(timezone.utc).isoformat(), "field": "QQ",
        "variables": list(map(str, variables)),
        "primitive_coefficient_matrix": [[str(v) for v in row] for row in coefficient_matrix.tolist()],
        "intersection_generators": [str(f.as_expr()) for f in difference_generators],
        "rank_one_minor_generators": [str(f.as_expr()) for f in rank_one_generators],
        "degree2_intersection_rank": 45, "degree2_rank_one_rank": 55,
        "degree3_intersection_rank": 560,
        "I_subset_J": True, "maximal_ideal_times_J_subset_I": True,
        "nilradical_quotient_dimension": 10, "nilradical_degree": 2,
        "groebner_basis_span_matches_J2": True,
        "buchberger_pairs_checked": checked_pairs, "coprime_pairs_skipped": coprime_pairs,
        "hilbert_numerator_J_denominator_1_minus_t_power16": list(numerator),
        "hilbert_numerator_I_denominator_1_minus_t_power16": str(intersection_numerator),
        "hilbert_series_J": "(1+11*t+11*t^2+t^3)/(1-t)^5",
        "hilbert_series_I": "(1+11*t+11*t^2+t^3)/(1-t)^5 + 10*t^2",
        "hilbert_recursion_states": monomial_hilbert_numerator.cache_info().currsize,
        "dimension": 5, "degree": 24,
        "proof_note": "research/notes/ELLIPTIC_SECTION_INTERSECTION_20261101.md",
        "scope": "Scheme-theoretic intersection of the two reduced elliptic graph components; no classification of their surrounding S61 local rings.",
        "elapsed_seconds": round(time.monotonic()-started, 3),
        "input_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in paths},
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
    print(json.dumps({"all_passed": True, "elapsed_seconds": report["elapsed_seconds"],
                      "output": str(args.output)}))


if __name__ == "__main__":
    main()

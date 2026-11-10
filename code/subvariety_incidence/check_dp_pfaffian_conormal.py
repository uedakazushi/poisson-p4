#!/usr/bin/env python3
"""Exact Hilbert-function test for the anticanonical-Pfaffian incidence of (2,2) del Pezzo surfaces.

The geometric inclusion Pf(pi) in I_S I_C is proved in the accompanying TeX;
this script independently checks the dimension statements on a smooth rational
diagonal pencil over QQ. It is not a test for existence of rank-four extensions.
"""
import argparse
import itertools as it
import json
from pathlib import Path

import sympy as s
from sympy.polys.matrices import DomainMatrix

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, default=Path("work/pfaffian_conormal_criterion.json"))
parser.add_argument("--check", type=Path)
args = parser.parse_args()
assert s.__version__ == "1.14.0"

x = s.symbols("x0:5")
q0 = sum(z * z for z in x)
q1 = sum(i * z * z for i, z in enumerate(x))
h = sum(x)

def monomials(d):
    return [s.prod(x[j] for j in indices) for indices in it.combinations_with_replacement(range(5), d)]

def exact_rank(polynomials, degree):
    exponents = [s.Poly(m, *x).monoms()[0] for m in monomials(degree)]
    columns = []
    for p in polynomials:
        terms = dict(s.Poly(p, *x).terms())
        columns.append(s.Matrix([terms.get(exp, 0) for exp in exponents]))
    return DomainMatrix.from_Matrix(s.Matrix.hstack(*columns)).convert_to(s.QQ).rank()

I5 = [q * m for q in (q0, q1) for m in monomials(3)]
I2_5 = [q * m for q in (q0*q0, q0*q1, q1*q1) for m in monomials(1)]
hI4 = [h * q * m for q in (q0, q1) for m in monomials(2)]
rI = exact_rank(I5, 5)
rI2 = exact_rank(I2_5, 5)
rHI4 = exact_rank(hI4, 5)
rProd = exact_rank(I2_5 + hI4, 5)
rIC2 = exact_rank([q0, q1] + [h*xx for xx in x], 2)
assert (rI, rI2, rHI4, rProd, rIC2) == (65, 15, 29, 41, 7)

report = {
    "date": "2026-11-11", "version": "v17", "status": "PASS",
    "sympy": s.__version__, "field": "QQ", "ambient": "P4",
    "S": "smooth complete intersection (2,2)",
    "C": "anticanonical divisor (2,2,1)",
    "C_degree": 4, "C_arithmetic_genus": 1,
    "example_q0": str(q0), "example_q1": str(q1), "example_h": str(h),
    "quadrics_through_C_dimension": rIC2,
    "smooth_elliptic_quartic_in_hyperplane_parameter_dimension": 20,
    "smooth_del_pezzo_S_with_fixed_C_affine_fibre_dimension": 10,
    "incidence_pairs_S_and_h_dimension": 30,
    "degree5_total_polynomial_dimension": len(monomials(5)),
    "degree5_ideal_S_dimension": rI,
    "degree5_ideal_S_squared_dimension": rI2,
    "degree5_h_times_ideal_S_degree4_dimension": rHI4,
    "degree5_ideal_S_times_ideal_C_dimension": rProd,
    "degree5_extra_linear_conditions_on_quintics_containing_S": rI-rProd,
    "universal_proof": "For S lci and C=V(h)|S an effective Cartier divisor, I_S cap I_C^2=I_S I_C. If pi preserves S and pi|S corresponds to h, pi vanishes as a bivector on C; hence Pf(pi) in I_C^2 and Pf(pi) in I_S. Hilbert degree counts hold for every smooth S and h not vanishing identically on S.",
    "scope": "Necessary Pfaffian condition for rank-four ambient extensions, not sufficient for existence; no new rank-four example or universal rigidity proven."
}
if args.check is not None:
    assert report == json.loads(args.check.read_text(encoding="utf-8")), "certificate mismatch"
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print("PASS: (I_S)_5=65, (I_S I_C)_5=41; 24 independent linear conditions")

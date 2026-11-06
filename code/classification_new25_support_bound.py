#!/usr/bin/env python3
"""A small support-dimension certificate for the existing NEW25 proof.

Only eight already-certified ordinary initial forms are needed.  This does not
compute the full initial ideal or assert reducedness of the local Poisson ring.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import resource
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
SOURCE = "reports/Q_POSITIVE_WEIGHT6_INITIAL4_20261103_v63/certificate.json.gz"


def bareiss(matrix):
    rows = [list(map(int, row)) for row in matrix]
    previous, sign, n = 1, 1, len(rows)
    for k in range(n-1):
        pivot = next(i for i in range(k, n) if rows[i][k])
        if pivot != k:
            rows[pivot], rows[k] = rows[k], rows[pivot]
            sign = -sign
        value = rows[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                numerator = rows[i][j]*value-rows[i][k]*rows[k][j]
                assert numerator % previous == 0
                rows[i][j] = numerator//previous
            rows[i][k] = 0
        previous = value
    return sign*rows[-1][-1]


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU, (50, 51))
    resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))
    begun = time.monotonic()
    source = ROOT/SOURCE
    data = json.loads(gzip.decompress(source.read_bytes()))
    z = s.symbols("z0:8")
    forms = [s.sympify(v) for key in ("quadratic_generators", "cubic_new_generators", "quartic_new_generators")
             for v in data[key]]
    assert len(forms) == 21
    a, b, c = z[1], z[5], z[7]
    # Fixed exact cofactors, independently obtained from degree-three products.
    cofactors = [
        [40*a/21+720*b/7-268*c/21,
         488*a/21-1440*b/7+232*c/3,
         -152*a/7+480*b/7-536*c/7],
        [-11*a/42+20*b/7+5*c/21,
         19*a/21-40*b/7+4*c/3,
         -8*a/21+40*b/21-26*c/21],
    ]
    for target, row in zip((a**3, c**3), cofactors):
        assert s.expand(sum(v*forms[k] for v, k in zip(row, (2, 3, 4)))-target) == 0
    q0 = s.factor(forms[0].subs({a:0, c:0}))
    assert q0 == -64*b*(z[4]+2*z[6])
    branch_a = {a:0, c:0, b:0}
    f, g = [s.Poly(forms[k].subs(branch_a), z[3], z[4], z[6]).clear_denoms()[1].primitive()[1].as_expr()
            for k in (16, 18)]
    # A common coefficient-only divisor of f must divide its leading coefficient
    # -3531600*(z4+2*z6); its constant coefficient does not have that factor.
    lc = s.Poly(f, z[3]).LC()
    constant = f.subs(z[3], 0)
    assert s.expand(lc+3531600*(z[4]+2*z[6])) == 0
    assert s.expand(constant.subs(z[4], -2*z[6])+194040*z[6]**4) == 0
    # A nonzero specialized Sylvester determinant excludes a common divisor
    # of positive z3 degree in QQ(z4,z6)[z3].  No Groebner basis is used.
    ff, gg = [s.Poly(v.subs({z[4]:1, z[6]:0}), z[3]) for v in (f, g)]
    assert ff.degree() == gg.degree() == 3
    sylvester = []
    for coeffs in (ff.all_coeffs(), gg.all_coeffs()):
        for shift in range(3):
            sylvester.append([0]*shift+list(map(int, coeffs))+[0]*(2-shift))
    determinant = bareiss(sylvester)
    assert determinant == -27062984035539000000*2558850183801
    assert s.Matrix(sylvester).det() == determinant
    branch_b = {a:0, c:0, z[4]:-2*z[6]}
    cubic = s.factor(forms[8].subs(branch_b, simultaneous=True))
    assert s.expand(cubic+s.Rational(8, 3)*b*z[6]**2) == 0
    quartic = s.factor(forms[17].subs({a:0,c:0,z[4]:0,z[6]:0}))
    assert quartic == 1433272320*b*z[2]**3
    report = dict(
        date_JST="2026-11-07", task_id="A", all_passed=True,
        selected_initial_indices_zero_based=[0,2,3,4,8,16,17,18],
        selected_initial_forms={str(k):str(forms[k]) for k in (0,2,3,4,8,16,17,18)},
        power_targets=[str(a**3),str(c**3)],
        cofactor_generator_indices=[2,3,4],
        power_cofactors=[[str(v) for v in row] for row in cofactors],
        quadratic_restriction=str(q0),
        branch_a=dict(variables=[str(v) for v in (z[3],z[4],z[6])],
                      primitive_quartics=[str(f),str(g)],
                      specialized_sylvester_matrix=sylvester,
                      determinant=str(determinant),
                      leading_coefficient=str(lc),
                      constant_on_leading_factor="-194040*z6**4",
                      coprime_over_C=True, dimension_upper=3),
        branch_b=dict(cubic_restriction=str(cubic),
                      quartic_at_z4_z6_zero=str(quartic),
                      additional_linear_support_ideal="(z1,z2,z4,z6,z7)",
                      dimension_upper=3),
        transverse_support_dimension_upper=3,
        projective_smooth_orbit_factor_dependency=22,
        actual_projective_local_dimension_upper=25,
        dependencies=["Existing independent 350-equation formal8 presentation and initial-element inclusion",
                      "Existing 22-dimensional projective smooth orbit factor"],
        not_claimed=["Full initial ideal generation or Hilbert function", "Reducedness, smoothness or branch exhaustion", "Global Poisson exhaustion"],
        input_sha256={SOURCE:hashlib.sha256(source.read_bytes()).hexdigest(),
                      "code/classification_new25_support_bound.py":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        seconds=time.monotonic()-begun)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(dict(all_passed=True, transverse_dimension_upper=3,
                         sylvester_determinant=str(determinant),seconds=report["seconds"])))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args().output)

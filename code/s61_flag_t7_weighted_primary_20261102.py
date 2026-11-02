"""Exact primary decomposition of the fixed-beta T7 flag weighted curve."""
import argparse
import hashlib
import json
import os
import subprocess
import tempfile
import time
from pathlib import Path

import sympy as sym

from s61_flag_linear_artifact_weighted_engine_20261102 import build
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_model import X, J
from weight1_structure import sp_basis

ROOT = Path(__file__).resolve().parents[1]
STEM = "S61_FLAG_T7_WEIGHTED_PRIMARY_20261102"
NOTE = f"research/notes/{STEM}.md"
CORE = -sym.Rational(91, 400)
WEIGHTS = (1, 1, 2, 2, 1, 2, 2)
POINTS_REPORT = "reports/S61_FLAG_REMAINING_CORE_WEIGHTED_POINTS_20261102.json"
BOUNDARY_REPORT = "reports/S61_FLAG_T7_BOUNDARY_20261101.json"
FRAME_REPORT = "reports/S61_FLAG_UNIFORM_OPEN_20261101.json"
ROW_ENGINE = "code/s61_flag_ordinary_weighted_rows_20261102.cpp"
M2 = "/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary"


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def matrix(values):
    return sym.Matrix(values).applyfunc(sym.sympify)


def actual_reduced_line(model, frame, tangent, inverse):
    previous = json.loads((ROOT / FRAME_REPORT).read_text())
    stabilizer = matrix(previous["frame_change"]["stabilizer"])
    derivative = matrix(previous["frame_change"]["beta_derivative"])
    parameters = sym.symbols("b0:15")
    values = sym.zeros(15, 1)
    values[0] = values[12] = 1
    values[2] = CORE
    spatial_change = sym.eye(4)
    finite_word = []
    for column in [4, 2]:
        if column == 4:
            scale = (sym.Rational(7, 20) * values[12] - values[3]) / values[0]
        else:
            scale = -4 * (values[13] + 5 * values[8]) / (3 * values[12])
        generator = sum((stabilizer[i, column] * g for i, g in enumerate(sp_basis())), sym.zeros(4))
        assert generator**2 == sym.zeros(4)
        spatial_change *= sym.eye(4) + scale * generator
        action = derivative[:, column].jacobian(parameters)
        assert action**4 == sym.zeros(15)
        values = sum((scale**k * action**k / sym.factorial(k) for k in range(4)), sym.zeros(15)) * values
        finite_word.append([column, str(scale)])
    assert finite_word == [[4, "7/20"], [2, "0"]]
    assert values[1] == values[2] == 0 and 20 * values[3] == 7 * values[12] and values[13] == -5 * values[8]
    epsilon = sym.Symbol("e")
    x, y, z, w = X[1:]
    spatial = sym.Matrix([x, y, z, w])
    bivector = sym.zeros(4)
    bivector[1, 2] = 3 * epsilon * w**2
    bivector[1, 3] = -x**2
    bivector -= bivector.T
    pulled = spatial_change * bivector.subs(dict(zip(spatial, spatial_change.inv() * spatial)), simultaneous=True) * spatial_change.T
    ambient = X[0]**2 * J
    ambient[1:, 1:] += pulled
    one = [i for i, weight in enumerate(model.reduced["parameter_weights"]) if weight == 1]
    two = [i for i, weight in enumerate(model.reduced["parameter_weights"]) if weight == 2]
    q = matrix(model.normalize(ambient)["parameters"])[two, :]
    assert q[25] == -1 and q[19] == q[26] == 0
    beta = dict.fromkeys(one, sym.S.Zero)
    beta.update({4: 1, 12: CORE, 48: 1})
    _, equations = specialize(model, beta, True)
    for equation in equations:
        value = sum(coefficient * sym.prod(q[j] for j in monomial) for monomial, coefficient in equation.items() if len(monomial) == 2)
        assert sym.expand(value) == 0
    free = []
    for i in range(7):
        row = next(j for j in range(45) if list(tangent.row(j)) == list(sym.eye(7).row(i)) and not any(inverse.row(j)))
        assert row < 44
        free.append(frame["chart"][row])
    normal = [q[i] for i in free]
    assert normal == [0, 0, -sym.Rational(147, 200) * epsilon, -sym.Rational(21, 20) * epsilon, 0, 0, 0]
    return dict(word=finite_word, normal_q_indices=free, normal_coordinates=list(map(str, normal)), raw_q_support={str(i): str(v) for i, v in enumerate(q) if v}, all135_at_t_zero=True)


def residuals(temp, frame, environment):
    denominator = temp / "den"
    output = temp / "rows"
    executable = temp / "rows-native"
    denominator.write_text("\n".join("1 " + frame[key] for key in ["ud", "td"]) + "\n")
    include = ROOT / "work/flint-headers"
    subprocess.run(["c++", "-O3", "-std=c++17", "-I" + str(include / "flint/3.5.0/include"), "-I" + str(include / "mpfr/4.2.2/include"), "-I/opt/homebrew/include", str(ROOT / ROW_ENGINE), str(ROOT / "work/m2-runtime/lib/libflint.24.0.dylib"), "-L/opt/homebrew/lib", "-lgmp", "-o", str(executable)], check=True)
    subprocess.run([str(executable), str(temp / "raw/raw.txt"), str(denominator), str(output)], env=environment, check=True, capture_output=True, timeout=60)

    def read_constant(stream):
        terms = stream.readline().split()
        assert len(terms) - 1 == int(terms[0]) and len(terms) <= 2
        return terms[1] if len(terms) == 2 else "0"

    expressions = []
    with output.open() as stream:
        assert int(stream.readline()) == 136
        for i in range(136):
            index, weight, count = map(int, stream.readline().split())
            assert index == i and sym.Rational(read_constant(stream))
            content = sym.Rational(read_constant(stream))
            assert not count or content
            terms = []
            for _ in range(count):
                exponents = tuple(map(int, stream.readline().split()))
                assert len(exponents) == 7 and sum(e * w for e, w in zip(exponents, WEIGHTS)) == weight
                coefficient = read_constant(stream)
                monomial = "".join(f"*x{j}^{e}" for j, e in enumerate(exponents) if e)
                terms.append(f"({coefficient}){monomial}")
            expressions.append("+".join(terms) or "0")
    return expressions


# The graph and square-zero coordinates give the primary structure directly.
M2_CHECKS = """
I=trim I;
h=x0+2*x1+3*x4;
g1=1080612*x3+732300*x0*x4-5925684*x1*x4-9697992*x4^2-57725*x5+1043668*x6;
g2=1350765*x2+1098060*x0*x4-20592192*x1*x4-37235472*x4^2-611885*x5+4146964*x6;
N=ideal(x0+4*x4,x1+2*x4,x5-2*x6,2*x4^2-x6);
C=ideal(g1,g2)+N^2;
Pc=ideal(x1+2*x4,x0+4*x4,x5-2*x6,2*x4^2-x6,2*x3+x6,x2+2*x6);
E=ideal(45*x1-284*x4,18*x0-49*x4,x6,x5,x4^2,10*x2-7*x3);
Pe=ideal(x4,x1,x0,x6,x5,10*x2-7*x3);
assert(C==saturate(I,ideal h));assert(E==saturate(I,ideal(2*x3+x6)));
assert(radical C==Pc);assert(radical E==Pe);
assert(C+E==E+ideal(x3));
U=intersect(C,E);K=U/I;
assert(not isSubset(C,E) and not isSubset(E,C));
assert(tval%C==0 and tval%E==0 and tval%I!=0);
assert(dim K==0 and numColumns basis K==106);
assert(apply(toList(0..7),i->hilbertFunction(i,K))=={0,0,4,16,34,43,9,0});
assert(dim(R/(C+E))==0 and numColumns basis(R/(C+E))==2);
print("T7_PRIMARY_PASS");
"""


def verify():
    start = time.monotonic()
    points = json.loads((ROOT / POINTS_REPORT).read_text())
    boundary = json.loads((ROOT / BOUNDARY_REPORT).read_text())
    for prior in [points, boundary]:
        assert prior["all_passed"]
        for path, digest in prior["input_sha256"].items():
            assert sha(path) == digest, path
    with tempfile.TemporaryDirectory(prefix="flag-t7-primary-") as directory:
        temp = Path(directory)
        frame, model, _, _, tangent, inverse = build(temp / "raw", CORE)
        assert frame == points["points"][2]["frame"]
        line = actual_reduced_line(model, frame, tangent, inverse)
        environment = os.environ | {"DYLD_LIBRARY_PATH": str(ROOT / "work/m2-runtime/lib") + ":/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib"}
        expressions = residuals(temp, frame, environment)
        script = temp / "check.m2"
        script.write_text("R=QQ[x0,x1,x2,x3,x4,x5,x6,Degrees=>{1,1,2,2,1,2,2}];\nI=ideal(" + ",".join(expressions[:135]) + ");\ntval=" + expressions[135] + ";\n" + M2_CHECKS)
        run = subprocess.run([M2, "--script", str(script)], env=environment, capture_output=True, text=True, timeout=100)
        assert run.returncode == 0 and "T7_PRIMARY_PASS" in run.stdout, (run.stdout, run.stderr)
    inputs = [str(Path(__file__).relative_to(ROOT)), NOTE, POINTS_REPORT, BOUNDARY_REPORT, FRAME_REPORT, ROW_ENGINE, *points["input_sha256"], *boundary["input_sha256"]]
    return dict(id=STEM, all_passed=True, field="QQ", core=str(CORE), normal_weights=WEIGHTS, actual_T7_curve=line,
                flag_primary="ideal(g1,g2)+(x0+4x4,x1+2x4,x5-2x6,2x4^2-x6)^2",
                T7_primary="(45x1-284x4,18x0-49x4,x6,x5,x4^2,10x2-7x3)", generic_lengths=[5,2],
                branch_parameters_weights=[1,2], T7_squarezero_weight=1, scheme_intersection_Hilbert=[1,1],
                vertex_primary="I+(x0,x1,x2,x3,x4,x5,x6)^7", primary_intersection_from_weight_bound=True, irredundancy_from_distinct_minimal_primes_and_nonzero_K=True,
                embedded_difference_Hilbert=[0,0,4,16,34,43,9], embedded_difference_length=106,
                t_in_both_isolated_primary_ideals=True, t_nonzero_in_S=True,
                scope="Specified fixed beta, zero nuisance, weighted two-coordinate cut. Relative extension is separate.",
                proof_note=NOTE, input_sha256={p: sha(p) for p in dict.fromkeys(inputs)}, seconds=time.monotonic()-start)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    result_path = parser.parse_args().output
    result = verify()
    with result_path.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(f"PASS {result['seconds']:.3f}s, {len(result['input_sha256'])} hashes")

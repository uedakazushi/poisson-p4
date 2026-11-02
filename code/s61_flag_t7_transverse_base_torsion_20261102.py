"""The actual DVR transverse ring along a canonical T7 nuisance direction."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import sympy as sym

from s61_flag_transverse_obstruction_20261101 import specialize
from s61_model import S61

ROOT = Path(__file__).resolve().parents[1]
STEM = "S61_FLAG_T7_TRANSVERSE_BASE_TORSION_20261102"
NOTE = f"research/notes/{STEM}.md"
PRIMARY = "reports/S61_FLAG_T7_WEIGHTED_PRIMARY_20261102.json"
ORDINARY = "reports/S61_FLAG_ORDINARY_WEIGHTED_FINITE_KERNEL_20261102.json"
C = -sym.Rational(5121, 30976)


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def verify():
    start = time.monotonic()
    primary = json.loads((ROOT / PRIMARY).read_text())
    for path, digest in primary["input_sha256"].items():
        assert sha(path) == digest, path
    ordinary = json.loads((ROOT / ORDINARY).read_text())
    weights = dict(zip(ordinary["frame"]["chart"], ordinary["raw_coordinate_weights"]))
    weights.update({25: 0, 45: 5})
    model = S61()
    one = [i for i, w in enumerate(model.reduced["parameter_weights"]) if w == 1]
    tau, epsilon = sym.symbols("tau e")
    beta = dict.fromkeys(one, sym.S.Zero)
    beta.update({4: 1, 12: -sym.Rational(91, 400), 48: 1, 34: -tau / 5, 49: tau})
    _, equations = specialize(model, beta, True)
    point = sym.zeros(45, 1)
    for key, value in primary["actual_T7_curve"]["raw_q_support"].items():
        point[int(key)] = sym.sympify(value)
    for equation in equations:
        at_point = sum(c * sym.prod(point[j] for j in mon) for mon, c in equation.items() if len(mon) == 2)
        assert sym.expand(at_point) == 0
        characters = set()
        for mon, c in equation.items():
            for (degree,), coefficient in sym.Poly(c, tau).terms():
                if coefficient:
                    characters.add(sum(weights[i] for i in mon) + (2-len(mon))*weights[45] - degree)
        assert len(characters) <= 1
    assert weights[24] == 2 and weights[35] == 1
    point = point.subs(epsilon, 1)
    frozen = {19, 24, 25, 26}
    variables = [j for j in range(45) if j not in frozen] + [45]
    location = {j: i for i, j in enumerate(variables)}
    t_index = location[45]
    eta_index = location[35]
    jacobian = sym.zeros(135, 42)
    quadratics = []
    for row, equation in enumerate(equations):
        h = {}
        for mon, coefficient in equation.items():
            j, k = (list(mon) + [45, 45])[:2]
            for left, right in [(j, k), (k, j)]:
                if left in location and right < 45:
                    jacobian[row, location[left]] += coefficient * point[right]
            if j in location and k in location:
                h[location[j], location[k]] = coefficient
        quadratics.append(h)
    jacobian = jacobian.applyfunc(sym.cancel)
    central = jacobian.subs(tau, 0)
    tangent_vectors = central.nullspace()
    assert len(tangent_vectors) == 1
    tangent = tangent_vectors[0]
    assert tangent[eta_index] == 1 and tangent[t_index] == 0
    pivots = list(central.to_DM().rref()[1])
    assert eta_index not in pivots
    selected = list(central[:, pivots].T.to_DM().rref()[1])
    inverse = central[selected, pivots].inv()
    right = sym.zeros(42, 41)
    for i, j in enumerate(pivots):
        right[j, :] = inverse[i, :]
    assert central[selected, :] * right == sym.eye(41)
    left = sym.Matrix.hstack(*central.T.nullspace()).T
    assert left * central == sym.zeros(94, 42)

    depth = 4
    linear = [jacobian.applyfunc(lambda v: sym.Poly(v, tau).nth(k)) for k in range(depth+1)]
    hseries = []
    for h in quadratics:
        hseries.append([(i, j, [sym.Poly(c, tau).nth(k) for k in range(depth+1)]) for (i, j), c in h.items()])

    def coefficient(order, path):
        out = sym.zeros(135, 1)
        for k in range(1, min(len(path), order+1)):
            out += linear[order-k] * path[k]
        for row, terms in enumerate(hseries):
            for i, j, cs in terms:
                for k in range(1, min(len(path), order)):
                    if path[k][i]:
                        for l in range(1, min(len(path), order-k+1)):
                            if path[l][j] and cs[order-k-l]:
                                out[row] += cs[order-k-l] * path[k][i] * path[l][j]
        return out

    pure_eta = sym.Matrix([sum(cs[0]*tangent[i]*tangent[j] for i,j,cs in terms) for terms in hseries])
    mixed = linear[1] * tangent
    forms = left * pure_eta.row_join(mixed)
    pivot = next(i for i in range(forms.rows) if forms[i, 0])
    assert forms.rank() == 1 and -forms[pivot, 1]/forms[pivot, 0] == C
    assert forms[73, 0] == -sym.Rational(30976, 1644875) and forms[73, 1] == -sym.Rational(5121, 1644875)
    path = [sym.zeros(42, 1), C * tangent]
    choices, obstructions = {}, {}
    for order in range(2, depth+1):
        current = coefficient(order, path)
        if order >= 3:
            changed = path.copy()
            changed[-1] = changed[-1] + tangent
            variation = coefficient(order, changed) - current
            normal = left * current
            direction = left * variation
            assert direction[pivot] != 0
            free_coefficient = -normal[pivot]/direction[pivot]
            choices[order-1] = str(free_coefficient)
            path[-1] += free_coefficient * tangent
            current += free_coefficient * variation
        residual = left * current
        obstructions[order] = {str(j): str(v) for j,v in enumerate(residual) if v}
        if order < depth:
            assert residual == sym.zeros(94, 1)
            path.append(-right * current[selected, :])
        else:
            assert obstructions[4] == {"73": "291826677388437/68777861447680000"}
    assert choices[2] == "0" and all(v[t_index] == 0 for v in path)
    inputs = [str(Path(__file__).relative_to(ROOT)), NOTE, PRIMARY, ORDINARY, *primary["input_sha256"]]
    return dict(id=STEM, all_passed=True, field="QQ", normalized_core="-91/400", path="b8=-tau/5,b13=tau, other nuisance=0",
                transverse_cuts="q19=q26=0,q24=3 at e=1,q25=-1", central_Jacobian_rank=len(pivots),
                parameter="eta=q35", central_tangent={str(variables[i]):str(v) for i,v in enumerate(tangent) if v},
                pivot_columns=pivots, selected_rows=selected, scalar_implicit_pivot=pivot,
                second_root_slope=str(C), scalar_path_coefficients=choices, projected_residuals=obstructions,
                homogenizer_along_second_test=[str(v[t_index]) for v in path],
                DVR_algebra="QQ[[tau]][eta]/(eta^2-c*tau*eta,tau^3*eta),c=-5121/30976",
                free_rank=1, torsion_length=3, special_length=2, homogenizer_zero=True,
                torus_weights=dict(e=2,tau=-1,eta=1,homogenizer=5), generic_e_transferred_by_quadratic_field_extension=True,
                canonical_two_cut_A_flat_T7_double_lift_excluded=True, full_core_relative_removal_claimed=False,
                proof_note=NOTE,input_sha256={p:sha(p) for p in dict.fromkeys(inputs)},seconds=time.monotonic()-start)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    output = parser.parse_args().output
    result = verify()
    with output.open("x") as stream:
        json.dump(result,stream,indent=2)
        stream.write("\n")
    print(f"PASS {result['seconds']:.3f}s, {len(result['input_sha256'])} hashes")

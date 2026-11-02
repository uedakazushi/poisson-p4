"""The nine-dimensional base and the height-one curvilinear chart."""
import json
from pathlib import Path

import sympy as sym

from s61_flag_transverse_obstruction_20261101 import specialize
from s61_model import S61, X

ROOT = Path(__file__).resolve().parents[1]
H11_REPORT = "reports/S61_FLAG_T7_BOUNDARY_20261101.json"
FLAG_REPORT = "reports/S61_FLAG_UNIFORM_OPEN_20261101.json"
PRIMARY = "reports/S61_FLAG_T7_WEIGHTED_PRIMARY_20261102.json"
FORMAL = "reports/S61_R0_MIXED_SURFACE_T7_FORMAL_20261101.json"


def matrix(data):
    return sym.Matrix(data).applyfunc(sym.sympify)


def base_data():
    m = S61()
    h11 = matrix(json.loads((ROOT / H11_REPORT).read_text())["fixed_beta_kernel"])
    incidence = json.loads((ROOT / FLAG_REPORT).read_text())
    b = sym.symbols("b0:15")
    flag = matrix(incidence["flag"]["beta_kernel"])
    differential = matrix(incidence["frame_change"]["beta_derivative"])[:, 4].jacobian(b)
    v = sym.symbols("v0:9")
    h11_parameters = sym.Matrix([1, *v[:7], 1, *v[7:]])
    scale = -sym.Rational(7, 20)
    flow = sum((scale**j / sym.factorial(j) * differential**j for j in range(4)), sym.zeros(15))
    beta = flag * flow * (h11 * h11_parameters)[1:, :]
    point = sym.zeros(46, 1)
    curve = json.loads((ROOT / PRIMARY).read_text())["actual_T7_curve"]["raw_q_support"]
    for i, value in curve.items():
        point[int(i)] = sym.sympify(value).subs(sym.Symbol("e"), 1)
    one = [i for i, w in enumerate(m.reduced["parameter_weights"]) if w == 1]
    two = [i for i, w in enumerate(m.reduced["parameter_weights"]) if w == 2]
    return m, v, beta, one, two, point


def verify_height_one(data):
    m, v, beta, one, two, point = data
    original = json.loads((ROOT / H11_REPORT).read_text())
    formal = json.loads((ROOT / FORMAL).read_text())["universal"]
    spatial_point = matrix(original["projective_q"])
    e = sym.Symbol("e")
    epsilon_tangent = spatial_point.diff(e)
    spatial_point = spatial_point.subs(e, 1)
    old_fixed = next(i for i, a in enumerate(spatial_point) if a)
    old_chart = [i for i in range(45) if i != old_fixed]
    pure = sym.zeros(45, 1)
    last_column = matrix(formal["tangent"])[:, 2]
    assert last_column[-1] == 0
    for i, j in enumerate(old_chart):
        pure[j] = last_column[i]
    old_b = sym.symbols("b0:11")
    old_sub = dict(zip(old_b, [1, *v[:7], 1, *v[7:]]))
    pure = (pure + pure[25] * spatial_point).applyfunc(lambda x: sym.cancel(x.subs(old_sub)))
    shear = sym.zeros(45, 1)
    shear[19] = 1
    change = sym.eye(5)
    change[4, 3] = -sym.Rational(7, 20)
    transform = sym.zeros(45, 45)
    needed = [i for i in range(45) if any(col[i] for col in [spatial_point, pure, epsilon_tangent, shear])]
    for i in needed:
        z = sym.zeros(61, 1)
        z[two[i]] = 1
        q = m.reconstruct(list(z))
        q = change * q.subs(dict(zip(X, change.inv() * sym.Matrix(X))), simultaneous=True) * change.T
        transform[:, i] = sym.Matrix(m.normalize(q)["parameters"])[two, :]
    center = transform * spatial_point
    assert center == point[:45, :]
    def normalized(col):
        moved = transform * col
        return moved + moved[25] * center
    tangents = normalized(epsilon_tangent).row_join(normalized(shear))
    transverse = normalized(pure)
    cut_minor = tangents[[24, 19], :].det()
    assert cut_minor == 3
    transverse -= tangents * tangents[[24, 19], :].inv() * transverse[[24, 19], :]
    coefficient = sym.factor(transverse[26] / transverse[35])
    denominator = 312*v[8]*v[7] + 312*v[1]*v[7] - 5400*v[4] + 1620*v[5] + 567*v[6] - 1216*v[7]**2
    assert sym.cancel(coefficient + 324*v[7] / (5*denominator)) == 0
    assert denominator.subs(v[7], 0) != 0
    # At this rational point on v7=0 the 42 eliminated columns are independent.
    sample = dict(zip(v, [1, 2, 3, 4, 5, 6, 7, 0, 10]))
    _, equations = specialize(m, dict(zip(one, beta.subs(sample))), True)
    columns = [i for i in range(46) if i not in [19, 24, 25, 35]]
    assert len(columns) == 42
    pos = {q:i for i,q in enumerate(columns)}
    jac = sym.zeros(135, 42)
    for i, row in enumerate(equations):
        for mon, coeff in row.items():
            j, k = (list(mon) + [45, 45])[:2]
            for q, u in [(j, k), (k, j)]:
                if q in pos and point[u]:
                    jac[i, pos[q]] += coeff * point[u]
    pivots = list(jac.T.to_DM().rref()[1])
    assert len(pivots) == 42
    det = jac[pivots, :].det(method="domain-ge")
    assert det
    return dict(height_one_prime="v7=a9", physical_cuts=["q19=0","q24=3"], regular_free_coordinate="q35",
                point_on_divisor=list(sample.values()), pivot_rows=pivots, pivot_columns=columns,
                pivot_determinant=str(det), reduced_cut_minor=str(cut_minor),
                q26_linear_in_q35=str(coefficient), valuation=1, nil_ideal_presentation_DVR_length=3)

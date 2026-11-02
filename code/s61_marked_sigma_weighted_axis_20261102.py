"""Exact positive-weight residue algebra on the marked-square axis."""
import json
import os
from pathlib import Path
import subprocess
import tempfile

import sympy as sym
from sympy.polys.rings import ring

from s61_model import BASIS, S61, X
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_flag_t7_weighted_primary_20261102 import M2

ROOT = Path(__file__).resolve().parents[1]
FLAG = "reports/S61_FLAG_UNIFORM_OPEN_20261101.json"
ACTION = "reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json"


def build():
    model = S61()
    weights = model.reduced["parameter_weights"]
    one = [i for i, w in enumerate(weights) if w == 1]
    two = [i for i, w in enumerate(weights) if w == 2]
    flag = json.loads((ROOT / FLAG).read_text())["flag"]
    base = sym.Matrix(flag["beta_kernel"]).applyfunc(sym.Rational)
    assert base == sym.eye(16)[:, 1:]
    power = [0, 1, -1, 1, -1]
    combined = []
    for degree, i in zip(weights, model.certificate["raw_free_indices"]):
        pair, monomial = BASIS[i]
        exponent = sym.Poly(monomial, *X).monoms()[0]
        raw = sum(power[j] for j in pair) - sum(a*b for a, b in zip(power, exponent))
        combined.append(degree - raw)
    extra = json.loads((ROOT / "reports/S61_UNIVERSAL_WEIGHT6_20261031.json").read_text())
    for row in model.reduced["equations"] + [extra["universal_new_equation"]]:
        assert len({sum(combined[i] for i in mon) for mon, _, _ in row}) == 1
    fixed = 25
    chart = [i for i in range(45) if i != fixed]
    assert combined[two[fixed]] == 6 and all(w % 2 == 0 for w in combined)
    qweight = [(6 - combined[two[j]]) // 2 for j in chart] + [3]
    assert [chart[i] for i in range(44) if qweight[i] == 0] == [27, 32]
    base_weight = [-combined[i] // 2 for i in one[1:]]
    assert base_weight[0] == 1 and max(base_weight[1:]) == 0

    rho = sym.Symbol("rho")
    center = sym.zeros(15, 1)
    center[1] = center[5] = 1
    center[6], center[14] = 3, 2
    assert all(base_weight[i] == 0 for i, v in enumerate(center) if v)
    beta = base * (center + rho*sym.eye(15)[:, 0])
    _, source = specialize(model, dict(zip(one, beta)), True)
    full_linear = sym.Matrix([
        [-row.get(tuple(sorted((fixed, j))), 0) for j in chart] + [-row.get((fixed,), 0)]
        for row in source])
    linear = full_linear.subs(rho, 0)
    cuts = [chart.index(j) for j in [0, 19, 27, 32]]
    columns = [i for i in range(45) if i not in cuts]
    reduced, pivots = linear[:, columns].to_DM().to_field().rref()
    reduced = reduced.to_Matrix()
    free = [j for j in range(len(columns)) if j not in pivots]
    tangent = sym.zeros(45, len(free))
    for k, j in enumerate(free):
        tangent[columns[j], k] = 1
    for i, j in enumerate(pivots):
        for k, z in enumerate(free):
            tangent[columns[j], k] = -reduced[i, z]
    selected = list(linear[:, [columns[j] for j in pivots]].T.to_DM().rref()[1])
    left = linear[selected, [columns[j] for j in pivots]].inv(method="DM")
    inverse = sym.zeros(45, len(selected))
    for k, j in enumerate(pivots):
        inverse[columns[j], :] = left[k, :]
    assert len(selected) == 28 and len(free) == 13
    assert linear * tangent == sym.zeros(135, 13) and linear[selected, :] * inverse == sym.eye(28)
    names = [f"v{i}" for i in range(13)] + ["rho"]
    nw = [qweight[columns[j]] for j in free] + [1]
    assert min(nw) == 1
    algebra, *x = ring(",".join(names), sym.QQ)

    def scalar_poly(value):
        return sum((sym.QQ.convert(c)*x[-1]**k[0] for k, c in sym.Poly(value, rho).terms()), algebra.zero)

    nonlinear = []
    for row in source:
        terms = []
        for mon, c in row.items():
            ids = [chart.index(j) for j in mon if j != fixed] + [44]*(2 - len(mon))
            c *= (-1)**mon.count(fixed)
            if not ids:
                assert c == 0
                continue
            if len(ids) == 1:
                c = sym.expand(c - c.subs(rho, 0))
            if c:
                terms.append((ids, scalar_poly(c)))
        nonlinear.append(terms)
    y = [sum((sym.QQ.convert(tangent[i,j])*x[j] for j in range(13) if tangent[i,j]), algebra.zero)
         for i in range(45)]

    def higher(i):
        result = algebra.zero
        for ids, c in nonlinear[i]:
            if len(ids) == 1:
                result += c*y[ids[0]]
            else:
                result += c*(y[ids[0]]*y[ids[1]])
        return result

    for d in range(1, max(qweight) + 1):
        for k in columns:
            if qweight[k] == d:
                y[k] -= sum((sym.QQ.convert(inverse[k,j])*higher(selected[j])
                             for j in range(28) if inverse[k,j]), algebra.zero)
    for i, value in enumerate(y):
        assert all(sum(a*b for a, b in zip(mon, nw)) == qweight[i] for mon in value)
    residuals = [sum((sym.QQ.convert(linear[i,j])*y[j] for j in columns if linear[i,j]), algebra.zero)
                 + higher(i) for i in range(135)]
    assert not any(residuals[i] for i in selected)
    assert all(len({sum(a*b for a, b in zip(mon, nw)) for mon in value}) <= 1 for value in residuals)
    assert y[-1] == x[12] and all(not y[j] for j in cuts)

    act = sym.Matrix(json.loads((ROOT / ACTION).read_text())["primitive_action"]).applyfunc(sym.sympify)
    b = base*center
    act = act.subs({sym.Symbol(f"u{i+1}"): v for i, v in zip(one, b)})
    orbit = sym.Matrix(flag["q_orbit_differential"]).applyfunc(sym.Rational)
    four = orbit[[chart.index(i) for i in [0, 19, 27]], :].col_join(act[0,:])
    assert four[:,[4,5,6,9]].det() == 4
    h13 = base[:,[i for i in range(15) if i not in [0,3,12]]].row_join(base[:,12]+sym.Rational(4,5)*base[:,3])
    assert h13.row_join(act).to_DM().rank() == 16
    return names, nw, residuals, dict(center_beta16=list(map(str,b)), q_weights=qweight,
                H15_base_weights=base_weight, positive_base="u5", nonpositive_q_base=32,
                ambient_orbit_cuts=[0,19,27], selected_rows=selected,
                selected_columns=[columns[j] for j in pivots], free_chart_columns=[columns[j] for j in free],
                normal_weights_with_rho=nw, transverse_Jacobian_rank=28,
                exact_residual_terms=sum(map(len,residuals)), homogenizer=names[12],
                orbit_normal_minor_columns=[4,5,6,9], orbit_normal_minor=4, H13_source_rank=16)


def artin_certificate(names, weights, residuals):
    def encode(polynomial):
        return "+".join(f"({c})" + "".join("*"+n+(f"^{d}" if d>1 else "")
                         for n,d in zip(names,mon) if d) for mon,c in polynomial.items()) or "0"
    source = "R=QQ["+ ",".join(names) +",Degrees=>{"+ ",".join(map(str,weights)) +"}];\n"
    source += "I=ideal("+ ",".join(encode(f) for f in residuals if f) +");\n"
    source += """
    P=ideal apply(0..12,i->R_i);
    A=saturate(I,ideal(rho)); U=I+ideal(rho);
    assert(dim(R/I)==1 and dim(R/U)==0 and numColumns basis(R/U)==70);
    assert(radical A==P and A==saturate(A,ideal(rho)));
    assert(numColumns basis(R/(A+ideal(rho)))==4 and v12%A==0);
    assert(toString(hilbertSeries(R/I,Reduce=>true))=="(1+7*T+23*T^2+35*T^3-32*T^4-30*T^5)/((1-T))");
    assert(v12%I!=0 and v12^2%I==0);
    print("MARKED_SIGMA_AXIS_PASS");
    """
    with tempfile.TemporaryDirectory(prefix="marked-sigma-axis-") as directory:
        path=Path(directory)/"exact.m2"
        path.write_text(source)
        env=os.environ|{"DYLD_LIBRARY_PATH":str(ROOT/"work/m2-runtime/lib")+":/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib"}
        run=subprocess.run([M2,"--script",str(path)],capture_output=True,text=True,env=env,timeout=100)
        assert run.returncode==0 and "MARKED_SIGMA_AXIS_PASS" in run.stdout,(run.stdout,run.stderr)
    return dict(central_dimension=1, rho_zero_length=70, rho_saturated_radical="(v0,...,v12)",
                rho_saturated_rank=4, t_in_saturation=True, central_t_nonzero=True, central_t_square_zero=True,
                central_Hilbert_numerator_over_one_minus_z=[1,7,23,35,-32,-30])

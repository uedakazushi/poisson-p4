"""Nine-parameter relative T7 transverse algebra and the height-one comparison."""
import argparse
import hashlib
import itertools
import json
import time
from functools import lru_cache
from pathlib import Path

import sympy as sym
from sympy import QQ
from sympy.polys.rings import ring
from s61_model import BASIS, X

from s61_flag_t7_relative_torsion_geometry_20261102 import ROOT, PRIMARY, FORMAL, base_data, verify_height_one

STEM = "S61_FLAG_T7_RELATIVE_TORSION_FULL_BASE_20261102"
NOTE = f"research/notes/{STEM}.md"
DVR = "reports/S61_FLAG_T7_TRANSVERSE_BASE_TORSION_20261102.json"
EXTRA = "reports/S61_UNIVERSAL_WEIGHT6_20261031.json"
C = -sym.Rational(5121, 30976)


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def verify():
    start = time.monotonic()
    data = base_data()
    m, variables, beta, one, two, point = data
    height_one = verify_height_one(data)
    prior = json.loads((ROOT / DVR).read_text())
    assert prior["all_passed"]
    for p, digest in prior["input_sha256"].items():
        assert sha(p) == digest
    polynomial_ring, *v = ring(",".join(str(x) for x in variables), QQ)
    scalar = lambda x: QQ.from_sympy(sym.Rational(x))
    zero = polynomial_ring.zero
    constant = lambda x: polynomial_ring.ground_new(scalar(x))
    parameter_forms = {i:polynomial_ring.from_expr(beta[j]) for j,i in enumerate(one)}
    quad_index = {i:j for j,i in enumerate(two)}
    chart = [i for i in range(46) if i not in [19, 24, 25, 26]]
    chart_index = {j:i for i,j in enumerate(chart)}
    equations = m.reduced["equations"] + [json.loads((ROOT / EXTRA).read_text())["universal_new_equation"]]

    @lru_cache(None)
    def coefficient_parts(monomial):
        product = polynomial_ring.one
        for j in monomial:
            product *= parameter_forms[j]
        homogeneous = [{} for _ in range(5)]
        for a, coeff in product.items():
            if sum(a) <= 4:
                homogeneous[sum(a)][a] = coeff
        return [polynomial_ring.from_dict(g) for g in homogeneous]

    spatial_weights = (0, 1, -1, 0, 0)
    combined = []
    for i, original in enumerate(m.certificate["raw_free_indices"]):
        pair, mon = BASIS[original]
        powers = sym.Poly(mon, *X).monoms()[0]
        raw_weight = sum(spatial_weights[j] for j in pair) - sum(x*y for x,y in zip(spatial_weights,powers))
        combined.append(m.reduced["parameter_weights"][i]-raw_weight)
    base_weights = []
    at_zero = dict.fromkeys(variables, 0)
    for j, value in enumerate(one):
        if beta[j].subs(at_zero):
            assert combined[value] == 0
    for variable in variables:
        observed = {-combined[i] for j,i in enumerate(one) if beta[j].diff(variable)}
        assert len(observed) == 1
        base_weights.append(observed.pop())
    assert max(base_weights) <= 0 and base_weights[7] == -1
    for equation in equations:
        characters = set()
        for mon, numerator, _ in equation:
            if numerator:
                spatial_degree = sum(i in quad_index for i in mon)
                characters.add(5*(2-spatial_degree) + sum((5-combined[i] if i in quad_index else -combined[i]) for i in mon))
        assert len(characters) == 1
    pairs = []
    central = sym.zeros(135, 42)
    for row, equation in enumerate(equations):
        terms = {}
        for mon, numerator, denominator in equation:
            bmon = tuple(i for i in mon if i in parameter_forms)
            qmon = [quad_index[i] for i in mon if i in quad_index]
            j, k = (qmon + [45, 45])[:2]
            coeff = scalar(sym.Rational(numerator, denominator))
            parts = coefficient_parts(bmon)
            target = terms.setdefault((j, k), [zero for _ in range(5)])
            for degree, term in enumerate(parts):
                target[degree] += coeff * term
            base_constant = QQ.to_sympy(coeff * parts[0].get((0,)*9, QQ.zero))
            for left, right in [(j,k),(k,j)]:
                if left in chart_index and point[right]:
                    central[row,chart_index[left]] += base_constant * point[right]
        pairs.append([(j,k,h) for (j,k),h in terms.items() if any(h)])
    kernel = central.nullspace()
    assert len(kernel) == 1
    tangent = kernel[0]
    assert tangent[chart_index[35]] == 1
    cols = list(central.to_DM().rref()[1])
    rows = list(central[:,cols].T.to_DM().rref()[1])
    inverse = central[rows,cols].inv()
    right = sym.zeros(42, 41)
    for i,j in enumerate(cols):
        right[j,:] = inverse[i,:]
    left = sym.Matrix.hstack(*central.T.nullspace()).T
    assert central[rows,:] * right == sym.eye(41) and left * central == sym.zeros(94,42)

    def linear(mat, polys):
        return [sum((scalar(mat[i,j])*polys[j] for j in range(mat.cols) if mat[i,j] and polys[j]),zero) for i in range(mat.rows)]

    def pack(polys):
        return [polys[chart_index[i]] if i in chart_index else zero for i in range(46)]

    def coefficient(degree, path):
        packed = [[constant(x) for x in point]] + [pack(q) for q in path]
        output = []
        for terms in pairs:
            val = zero
            for i,j,parts in terms:
                for k in range(min(degree+1,len(packed))):
                    if not packed[k][i]:
                        continue
                    for l in range(min(degree-k+1,len(packed))):
                        if packed[l][j] and parts[degree-k-l]:
                            val += parts[degree-k-l]*packed[k][i]*packed[l][j]
            output.append(val)
        return output

    # The fixed projective Q is a section across the entire H11 source.
    # The accepted exact H11 identity supplies all degrees.  Also test every
    # coefficient of the section through the degrees used below.
    for n in range(5):
        assert not any(coefficient(n, []))
    # Unit linear obstruction after dividing all equations by eta.
    pure = []
    normal = pack([constant(x) for x in tangent])
    for terms in pairs:
        z = sum((parts[0]*normal[i]*normal[j] for i,j,parts in terms),zero)
        pure.append(QQ.to_sympy(z.get((0,)*9,QQ.zero)))
    col = left * sym.Matrix(pure)
    pivot = next(i for i,z in enumerate(col) if z)
    multiplier = C*col
    path = [[scalar(C*x)*v[7] for x in tangent]]
    free_coefficients = {}
    final = None
    for order in range(2,5):
        raw = coefficient(order,path)
        projected = linear(left,raw)
        if order >= 3:
            # Obtain and check the response to a variation in the newly free
            # coefficient, retaining the complete nine-variable polynomial.
            trial = [list(z) for z in path]
            test = v[0]**(order-1)
            trial[-1] = [z+scalar(x)*test for z,x in zip(trial[-1],tangent)]
            changed = linear(left,coefficient(order,trial))
            assert all(b-a == scalar(g)*v[7]*test for a,b,g in zip(projected,changed,multiplier))
            choice = (-projected[pivot]).exquo(scalar(multiplier[pivot])*v[7])
            free_coefficients[str(order-1)] = str(choice)
            path[-1] = [z+scalar(x)*choice for z,x in zip(path[-1],tangent)]
            raw = coefficient(order,path)
            projected = linear(left,raw)
        if order < 4:
            assert not any(projected)
        else:
            # eta has first part C*v7.  Uniqueness of the pivot root shows
            # these lowest divided residuals equal the cubic initial forms
            # up to this fixed nonzero scalar.
            final = [z.exquo(v[7]) for z in projected]
        path.append([-z for z in linear(right,[raw[i] for i in rows])])
    support = [z for z in final if z]
    assert support and all(z.exquo(v[7]**3).is_ground for z in support)
    assert all(not q[chart_index[45]] for q in path)
    paths = [str(Path(__file__).relative_to(ROOT)), NOTE,
             "code/s61_flag_t7_relative_torsion_geometry_20261102.py", PRIMARY, FORMAL, DVR, EXTRA,
             "reports/S61_FLAG_T7_FLAT4_COVER_20261101.json", "research/notes/S61_FLAG_T7_FLAT4_20261101.md",
             *prior["input_sha256"]]
    return dict(id=STEM,all_passed=True,base_variables=[str(x) for x in variables],
                H11_parameters=["1",*map(str,variables[:7]),"1",*map(str,variables[7:])],
                normalized_beta=list(map(str,beta)), central_Q_support=primary_support(point),
                cuts=["q19=q26=0","q24=3","q25=-1"], central_rank=len(cols), scalar_pivot=pivot,
                second_root_first_term="-5121*v7/30976", higher_root_coefficients=free_coefficients,
                bottom_ideal_degree1_rank=0,bottom_ideal_degree2_rank=0,bottom_ideal_cubic_rank=1,
                true_cubic_initial_generator="v7^3",cubic_witness=[str(z) for z in support],
                height_one_comparison=height_one,base_weights=base_weights,exact_base_ideal="(v7^3)",second_root_is_v7_times_unit=True,
                full_algebra="A plus A/(v7^3) eta, eta^2=f(v) eta, f=v7 times a unit",
                homogenizer_zero=True,special_E_generic_point_only=True,flag_endpoint_ambient_removal=False,
                proof_note=NOTE,input_sha256={p:sha(p) for p in dict.fromkeys(paths)},seconds=time.monotonic()-start)


def primary_support(q):
    return {str(i):str(v) for i,v in enumerate(q) if v}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    destination = parser.parse_args().output
    report = verify()
    with destination.open("x") as stream:
        json.dump(report,stream,indent=2)
        stream.write("\n")
    print(f"PASS {report['seconds']:.3f}s, {len(report['input_sha256'])} hashes")

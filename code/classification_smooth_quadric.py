#!/usr/bin/env python3
"""Full invariant quadratic lifts for a smooth quadric normal Hamiltonian.

This finite entrance calculation does not assert component membership.
All ten homogeneous brackets and all 350 Jacobi coefficient slots are kept.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s


def run(ambient=False):
    started = time.monotonic()
    x = s.symbols("x0:5")
    weights = (0, 1, -1, 2, -2)
    q = x[0] ** 2 + x[1] * x[2] + x[3] * x[4]
    pairs = list(itertools.combinations(range(5), 2))
    monomials = [(x[i] * x[j], weights[i] + weights[j])
                 for i in range(5) for j in range(i, 5)]
    unknowns = []
    P = s.zeros(5)
    for i, j in pairs:
        for mon, weight in monomials:
            if weight == weights[i] + weights[j]:
                c = s.Symbol(f"p{len(unknowns)}")
                unknowns.append(c)
                P[i, j] += c * mon
        P[j, i] = -P[i, j]

    linear_polynomials = [sum(s.diff(P[i, j], x[i]) for i in range(5))
                          for j in range(5)]
    # Convention: X(f) = {q,f}/q, hence dq_i P^{ij}=q X^j.
    linear_polynomials += [sum(s.diff(q, x[i]) * P[i, j] for i in range(5))
                          - q * weights[j] * x[j] for j in range(5)]
    equations = []
    for poly in linear_polynomials:
        equations += s.Poly(s.expand(poly), *x).coeffs()
    M, b = s.linear_eq_to_matrix(equations, unknowns)
    solution = next(iter(s.linsolve((M, b), unknowns)))
    substitution = dict(zip(unknowns, solution))
    P = P.applyfunc(lambda f: s.expand(f.subs(substitution, simultaneous=True)))
    parameters = sorted(set().union(*(f.free_symbols for f in P)) - set(x), key=str)
    jacobi = []
    for i, j, k in itertools.combinations(range(5), 3):
        poly = s.expand(sum(P[i, r] * s.diff(P[j, k], x[r])
                            + P[j, r] * s.diff(P[k, i], x[r])
                            + P[k, r] * s.diff(P[i, j], x[r])
                            for r in range(5)))
        for coefficient in s.Poly(poly, *x).coeffs():
            if coefficient != 0:
                jacobi.append(s.factor(coefficient))
    jacobi = sorted(set(jacobi), key=str)
    pfaffians = []
    for missing in range(5):
        a, b, c, d = [i for i in range(5) if i != missing]
        pfaffians.append(s.factor((-1) ** missing * (
            P[a, b] * P[c, d] - P[a, c] * P[b, d] + P[a, d] * P[b, c])))
    a, b, u, v, t = (s.Symbol(name) for name in ("p12", "p14", "p15", "p16", "p18"))
    conditions = [a*b, a*(u+1), b*(v+t+3), t*(u+v-1),
                  8*(u+v)*t + 30*(u+v) + 7*t + 20]
    coefficient_monomials = sorted(set().union(*(
        set(s.Poly(f, *parameters).monoms()) for f in jacobi + conditions)))
    def coefficient_matrix(polys):
        rows = [s.Poly(f, *parameters).as_dict() for f in polys]
        return s.Matrix([[row.get(mon, 0) for mon in coefficient_monomials] for row in rows])
    assert coefficient_matrix(jacobi).rank() == 5
    assert coefficient_matrix(conditions).rank() == 5
    assert coefficient_matrix(jacobi + conditions).rank() == 5
    z = s.Symbol("z")
    branches = [
        ("I_central", {a:0,b:0,u:z,v:-s.Rational(2,3)-z,t:0}),
        ("I_a", {a:z,b:0,u:-1,v:s.Rational(1,3),t:0}),
        ("I_b", {a:0,b:z,u:s.Rational(7,3),v:-3,t:0}),
        ("II_central", {a:0,b:0,u:z,v:1-z,t:-s.Rational(10,3)}),
        ("II_a", {a:z,b:0,u:-1,v:2,t:-s.Rational(10,3)}),
        ("II_b", {a:0,b:z,u:s.Rational(2,3),v:s.Rational(1,3),t:-s.Rational(10,3)}),
    ]
    branch_reports = []
    for name, sub in branches:
        assert all(s.expand(f.subs(sub, simultaneous=True)) == 0 for f in jacobi)
        pf = [s.expand(f.subs(sub, simultaneous=True)) for f in pfaffians]
        D = s.factor(sum(x[i]*pf[i] for i in range(5))/5)
        assert all(s.expand(s.diff(D,x[i])-pf[i]) == 0 for i in range(5))
        cubic = s.cancel(D/q)
        assert s.denom(cubic).free_symbols.isdisjoint(x)
        degeneration_weights = ((0,0,3,1,1) if name == "I_b" else (0,3,0,1,1)
                                if name.startswith("I_") else (0,0,2,2,1)
                                if name == "II_b" else (0,2,0,1,2))
        limit = s.zeros(5)
        exponents = []
        branch_P = P.subs(sub, simultaneous=True)
        for i,j in pairs:
            for mon,coefficient in s.Poly(branch_P[i,j],*x).terms():
                power = degeneration_weights[i]+degeneration_weights[j]-sum(
                    mon[k]*degeneration_weights[k] for k in range(5))
                assert power >= 0
                exponents.append(power)
                if power == 0:
                    limit[i,j] += coefficient*s.prod(x[k]**mon[k] for k in range(5))
            limit[j,i] = -limit[i,j]
            assert s.cancel(limit[i,j]/(x[i]*x[j])).free_symbols.isdisjoint(x)
        limit_pf = []
        for missing in range(5):
            i,j,k,ell = [r for r in range(5) if r != missing]
            limit_pf.append((-1)**missing*(limit[i,j]*limit[k,ell]-limit[i,k]*limit[j,ell]+limit[i,ell]*limit[j,k]))
        limit_D = s.factor(sum(x[i]*limit_pf[i] for i in range(5))/5)
        assert s.expand(limit_D-s.Rational(10,3)*s.prod(x)) == 0
        branch_reports.append({"name":name, "substitution":{str(k):str(v) for k,v in sub.items()},
                               "actual_divisor":str(D), "cubic_factor":str(s.factor(cubic)),
                               "full_jacobi_zero":True, "all_pfaffian_minors_equal_gradient":True,
                               "diagonal_degeneration_weights":list(degeneration_weights),
                               "nonnegative_parameter_powers":sorted(set(exponents)),
                               "limit_is_log_canonical":True,"limit_actual_divisor":str(limit_D)})
    ambient_reports = []
    if ambient:
        from sympy.polys.matrices import DomainMatrix
        from classification_rank3_mixed_component_review import (
            COLUMNS, jacobi as full_jacobi, sparse, vector, minor_certificate)
        for name, sub in branches:
            if name.endswith("central"):
                continue
            point_sub = {k:s.sympify(v).subs(z,1) for k,v in sub.items()}
            point = P.subs(point_sub, simultaneous=True)
            scale = s.ilcm(*[s.denom(c) for f in point for c in s.Poly(f,*x).coeffs()])
            point = point.applyfunc(lambda f:s.expand(scale*f))
            assert all(s.denom(c) == 1 for f in point for c in s.Poly(f,*x).coeffs())
            ps = sparse(point, x)
            assert full_jacobi(ps) == s.zeros(350,1)
            columns = []
            for ij, exponents in COLUMNS:
                variation = {key:dict(value) for key,value in ps.items()}
                variation.setdefault(ij,{})[exponents] = variation.get(ij,{}).get(exponents,0)+1
                columns.append(full_jacobi(variation))
            D = s.Matrix.hstack(*columns)
            Xcoord = s.Matrix(x)
            vertical, orbit = [], []
            for i,j in itertools.product(range(5),repeat=2):
                field = s.zeros(5,1)
                field[i] = x[j]
                vertical.append(vector(sparse(Xcoord*field.T-field*Xcoord.T,x)))
                E = s.zeros(5)
                E[i,j] = 1
                orbit.append(vector(sparse(E*point+point*E.T-x[j]*point.diff(x[i]),x)))
            F = s.Matrix.hstack(*vertical,*orbit,vector(ps))
            assert D*F == s.zeros(350,F.cols)
            ranks = {label:len(DomainMatrix.from_Matrix(mat).convert_to(s.QQ).rref()[1])
                     for label,mat in (("D",D),("orbit_vertical_scale",F))}
            ambient_reports.append({"source":name,"point_parameter":1,"representative_scale":str(scale),"ranks_QQ":ranks,
                "orbit_dimension":ranks["orbit_vertical_scale"]-25,
                "projective_tangent_dimension":125-ranks["D"],
                "minors":{label:minor_certificate(mat,ranks[label])
                          for label,mat in (("D",D),("orbit_vertical_scale",F))}})
    source = Path(__file__).read_text()
    return {
        "scope": "complete affine linear invariant-lift source; all residual Jacobi retained",
        "case": "smooth q, normal Hamiltonian diag(0,1,-1,2,-2)",
        "coefficient_field": "QQ",
        "q": str(q), "weights": list(weights),
        "invariant_quadratic_unknowns": len(unknowns),
        "linear_matrix_shape": list(M.shape), "linear_rank": M.rank(),
        "parameters": list(map(str, parameters)),
        "brackets": {f"{i}{j}": str(P[i, j]) for i, j in pairs},
        "residual_jacobi_coefficients": list(map(str, jacobi)),
        "equivalent_five_conditions": list(map(str, conditions)),
        "condition_span_rank":5,
        "six_complete_sources":branch_reports,
        "four_nonzero_orbit_ambient_tests":ambient_reports,
        "signed_pfaffian_minors": list(map(str, pfaffians)),
        "elapsed_seconds": time.monotonic() - started,
        "producer_sha256": hashlib.sha256(source.encode()).hexdigest(),
        "producer_source": source,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--ambient", action="store_true")
    args = parser.parse_args()
    result = run(ambient=args.ambient)
    Path(args.output).write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in {"producer_source", "brackets", "signed_pfaffian_minors",
                                   "residual_jacobi_coefficients", "four_nonzero_orbit_ambient_tests"}},
                     ensure_ascii=False, indent=2))
    if args.ambient:
        print(json.dumps([{k:v for k,v in report.items() if k != "minors"}
                          for report in result["four_nonzero_orbit_ambient_tests"]], indent=2))

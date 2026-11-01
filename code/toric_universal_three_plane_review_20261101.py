#!/usr/bin/env python3
"""Independent coefficient/perspective audit of the universal known-five plane.

Rebuilds multiplication and associativity from tensor indices. The original
derive(), primitive_beta(), multiplication helpers, and normalize() are not
imported. The finite-neighborhood argument is reviewed in the attached note.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path
import time

import sympy as s
from s61_model import S61, X, BASIS, vertical_rref

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "reports/TORIC_UNIVERSAL_THREE_PLANE_20261101.json"
FIVE = ROOT / "reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json"
J = s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])


def matrix(data):
    return s.Matrix([[s.sympify(v) for v in row] for row in data])


def verify():
    started = time.monotonic()
    model = S61()
    certificate = json.loads(MAIN.read_text())
    five = json.loads(FIVE.read_text())
    one = certificate["weight1_u_indices_zero_based"]
    two = certificate["weight2_u_indices_zero_based"]
    u = list(map(s.Symbol, certificate["beta_variables"]))
    c = list(map(s.Symbol, certificate["c_variables"]))
    all_u = [s.S.Zero] * 61
    for i, value in zip(one, u, strict=True):
        all_u[i] = value
    y = s.Matrix(s.symbols("y0:4"))
    b = model.reconstruct(all_u)
    substitution = dict(zip(X, [1,*y], strict=True))
    p = s.Matrix(4,4,lambda i,j:s.expand(
        (b[i+1,j+1]-X[i+1]*b[0,j+1]+X[j+1]*b[0,i+1]).subs(substitution)))
    origin = dict.fromkeys(y, 0)
    old_linear = sum((p.diff(v).subs(origin)*v for v in y), s.zeros(4))
    alpha = s.Matrix([s.diff(old_linear[0,1]+old_linear[2,3],v)/5 for v in y])
    ja = J*alpha
    linear = old_linear-2*alpha.dot(y)*J+y*ja.T-ja*y.T
    beta = (-J*linear*J).applyfunc(s.expand)
    triples = list(combinations_with_replacement(range(4),3))
    cubic = dict(zip(triples,c,strict=True))
    # Explicit independent tensor construction; no toric model helper.
    multiplication = []
    for i in range(4):
        multiplication.append(s.Matrix(4,4,lambda a,j:sum(
            J[a,k]*((s.diff(beta[j,k],y[i])+s.diff(beta[i,k],y[j]))/3
                    +cubic[tuple(sorted((i,j,k)))]) for k in range(4))))
    assert all(s.trace(m)==0 for m in multiplication)
    assert all(multiplication[i][:,j]==multiplication[j][:,i]
               for i in range(4) for j in range(4))
    form = s.Matrix(4,4,lambda i,j:s.expand(s.trace(multiplication[i]*multiplication[j])/3))
    equations = []
    for i,j in combinations(range(4),2):
        commutator = multiplication[i]*multiplication[j]-multiplication[j]*multiplication[i]
        commutator += s.eye(4)[:,i]*form[j,:]-s.eye(4)[:,j]*form[i,:]
        equations.extend(s.expand(v) for v in commutator)
    cm = list(combinations_with_replacement(range(20),2))
    h = s.Matrix([[s.Poly(f,*c).coeff_monomial(c[a]*c[b]) for a,b in cm] for f in equations])
    assert not h.free_symbols
    hrows = certificate["H_pivot_rows"]
    hcols = certificate["H_pivot_columns"]
    assert h[hrows,hcols].det() != 0

    # Independently validate all 51 null relations; they span the full nullspace.
    relations = s.zeros(51,96)
    for row, entries in enumerate(certificate["all51_constant_combinations_of96"]):
        for col,value in entries:
            relations[row,col] = s.Rational(value)
    assert relations.rank()==51 and relations*h==s.zeros(51,210)
    l = matrix(certificate["all51_homogeneous_linear_consequences"])
    assert (relations*s.Matrix(equations)-l*s.Matrix(c)).applyfunc(s.expand)==s.zeros(51,1)

    pairs = list(combinations(range(4),2))
    ym = list(combinations_with_replacement(range(4),2))
    my = sum((y[i]*multiplication[i] for i in range(4)),s.zeros(4))
    # Product-rule coefficient of t² in (I+t my-t² y y^t form)J(... )^t.
    frame2 = -y*(y.T*form)
    p2 = (my*J*my.T+frame2*J+J*frame2.T).applyfunc(s.expand)
    coefficients = s.Matrix([s.Poly(p2[i,j],*y).coeff_monomial(y[a]*y[b])
                             for i,j in pairs for a,b in ym])
    symbol = s.Matrix([[s.Poly(f,*c).coeff_monomial(c[a]*c[b]) for a,b in cm]
                      for f in coefficients])
    reduction = matrix(certificate["quadratic_frame_symbol_reduction_matrix"])
    assert reduction*h[hrows,:]==symbol and symbol.rank()==45
    reduced = (coefficients-reduction*s.Matrix([equations[i] for i in hrows])).applyfunc(s.expand)
    assert all(s.Poly(f,*c).total_degree()<=1 for f in reduced)

    # Independent inverse perspective on a two-jet. The map from primitive
    # coordinates z to old coordinates is y=z/(1-alpha(z)). Its inverse is
    # z=y/(1+alpha(y)), so the degree-two correction is the expression below.
    correction = (alpha.dot(y)*linear+y*alpha.T*linear+linear*alpha*y.T
                  +alpha.dot(y)**2*J
                  +2*alpha.dot(y)*(y*alpha.T*J+J*alpha*y.T)).applyfunc(s.expand)
    old_coefficients = reduced+s.Matrix([
        s.Poly(correction[i,j],*y).coeff_monomial(y[a]*y[b]) for i,j in pairs for a,b in ym])
    # Constant quotient by E wedge End(V), applied to homogeneous quadratic
    # spatial coefficients. No call to S61.normalize() is made.
    reduced_vertical,pivots = vertical_rref()
    position = {pivot:row for row,pivot in enumerate(pivots)}
    quotient = s.zeros(45,60)
    free_raw = [model.certificate["raw_free_indices"][i] for i in two]
    for column,((i,j),(a,b)) in enumerate((pair,mon) for pair in pairs for mon in ym):
        raw = BASIS.index(((i+1,j+1),X[a+1]*X[b+1]))
        for row,target in enumerate(free_raw):
            quotient[row,column] = int(target==raw)
            if raw in position:
                quotient[row,column] -= reduced_vertical[position[raw],target]
    q = (quotient*old_coefficients).applyfunc(s.expand)
    center = s.Matrix(list(map(s.sympify,certificate["q_center"])))
    n = matrix(certificate["q_linear_matrix_N"])
    assert (q-center-n*s.Matrix(c)).applyfunc(s.expand)==s.zeros(45,1)
    section_paths = [ROOT/f"reports/{name}_UNIVERSAL_SECTION_20261101.json" for name in ("Q51","Q52")]
    q51,q52 = [s.Matrix(list(map(s.sympify,json.loads(path.read_text())["quadratic_section_q"])))
               for path in section_paths]
    assert (3*center-2*q51-q52).applyfunc(s.expand)==s.zeros(45,1)

    toric = s.Matrix(list(map(s.sympify,five["toric_S61_coordinates"])))
    sample = dict(zip(u,[toric[i] for i in one],strict=True))
    recipe = certificate["generic_plane_recipe"]
    l0,n0,center0 = l.subs(sample),n.subs(sample),center.subs(sample)
    kernel = matrix(recipe["sample_kernel_T"])
    assert l0.rank()==18 and l0*kernel==s.zeros(51,2) and kernel.rank()==2
    directions = (n0*kernel).row_join((q52-center).subs(sample))
    assert directions.rank()==3
    # Match the older geometric plane by one explicit affine coordinate change.
    tau,rho,sigma = s.symbols("tau rho sigma")
    embedding = s.Matrix(list(map(s.sympify,five["geometric_plane_embedding_S61_coordinates"])))
    target = s.Matrix([embedding[i] for i in two])-center0
    rows = list(directions.T.rref()[1])
    coordinates = directions[rows,:].inv()*target[rows,0]
    assert (directions*coordinates-target).applyfunc(s.expand)==s.zeros(45,1)
    linear_change = coordinates.jacobian((tau,rho,sigma))
    assert linear_change.det()!=0

    # Directly redo all134 restrictions in the old plane and verify the saved
    # CONSTANT row combinations into the five simple quadrics.
    from sympy.polys.rings import ring
    from sympy.polys.domains import QQ
    pr,*_ = ring("tau,rho,sigma",QQ)
    values = [pr.from_expr(f) for f in embedding]
    simple = [pr.from_expr(s.sympify(f)) for f in five["five_simple_quadrics"]]
    row_combinations = matrix(five["all134_geometric_quadric_span_coefficients"])
    restricted = []
    for equation in model.reduced["equations"]:
        value = pr.zero
        for mon,num,den in equation:
            term = pr(QQ(num,den))
            for i in mon:
                term *= values[i]
            value += term
        restricted.append(value)
    assert row_combinations.rank()==5
    assert all(restricted[i]==sum((pr.from_expr(row_combinations[i,j])*simple[j]
                                  for j in range(5)),pr.zero) for i in range(134))
    cubic_monomials = list(combinations_with_replacement(range(3),3))
    variables = (tau,rho,sigma)
    tops = [sum(co*s.prod(v**n for v,n in zip(variables,mon,strict=True))
                for mon,co in s.Poly(f.as_expr(),*variables).terms() if sum(mon)==2) for f in simple]
    cubic_matrix = s.Matrix([[s.Poly(v*f,*variables).coeff_monomial(s.prod(variables[j] for j in mon))
                             for mon in cubic_monomials] for f in tops for v in variables])
    cubic_rows = list(cubic_matrix.T.rref()[1])
    assert len(cubic_rows)==10
    paths = [Path(__file__),MAIN,FIVE,ROOT/"code/s61_model.py",*section_paths]
    paths += [model.source/name for name in ("symplectic_slice_certificate.json",
               "symplectic_slice_equations.json","symplectic_slice_reduced_equations.json")]
    return {"id":"TORIC_UNIVERSAL_THREE_PLANE_REVIEW_20261101","all_passed":True,
            "generated_utc":datetime.now(timezone.utc).isoformat(),"field":"QQ[16 jet variables,20 cubic variables]",
            "independently_rebuilt_associativity":True,"leading_rank":45,"all51_null_relations_verified":True,
            "full36variable_quadratic_symbol_reduction":True,"independent_inverse_perspective_q_reconstruction":True,
            "universal_center_identity":True,"sample_L_rank":18,"sample_plane_direction_rank":3,
            "sample_old_plane_affine_coordinate_change":list(map(str,coordinates)),
            "sample_old_plane_coordinate_change_determinant":str(linear_change.det()),
            "all134_restrictions_rebuilt_and_constant_span_verified":True,
            "leading_cubic_rank":10,"leading_cubic_minor":str(cubic_matrix[cubic_rows,:].det()),
            "proof_note":"research/notes/TORIC_UNIVERSAL_THREE_PLANE_REVIEW_20261101.md",
            "scope":"Generic reduced length5 only for the intersection with the specified rational three-plane; no containment of the full S61 fiber.",
            "elapsed_seconds":round(time.monotonic()-started,3),
            "input_sha256":{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    report=verify()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(report,stream,indent=2);stream.write("\n")
    print(json.dumps({"all_passed":True,"elapsed_seconds":report["elapsed_seconds"],"output":str(args.output)}))


if __name__=="__main__":
    main()

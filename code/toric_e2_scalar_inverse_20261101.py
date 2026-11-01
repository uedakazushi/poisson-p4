#!/usr/bin/env python3
"""Universal scalar inverse and trace identities for the toric Jacobian-2 cover.

Verifies both contractions on a dominant six-dimensional slice, independently
checks slice dominance and uniqueness of the rational degree-nine companion,
and proves the remaining cubic algebra identities with formal c4,c6.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
from pathlib import Path
import time

import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

from weight1_structure import J, form_matrices, infinitesimal_action, sp_basis
from toric_degree3_scalar_20261101 import TRIPLES, multiplicity

ROOT = Path(__file__).resolve().parents[1]
D_PATH = ROOT/"reports/TORIC_DEGREE3_SCALAR_20261101.json"
U_PATH = ROOT/"reports/TORIC_DEGREE5_COVARIANT_20261101.json"
SLICE_PATH = ROOT/"reports/PRIMITIVE_CURVE_DOMINANT_SLICE_20261101.json"
V_PATH = ROOT/"reports/TORIC_JACOBIAN_TWO_TORSION_20261101.json"


def verify():
    started = time.monotonic()
    d_report = json.loads(D_PATH.read_text())
    u_report = json.loads(U_PATH.read_text())
    slice_report = json.loads(SLICE_PATH.read_text())
    v_report = json.loads(V_PATH.read_text())
    basis = s.Matrix(d_report["weight_basis_24_by16"])
    assert basis==s.Matrix(u_report["weight_basis_24_by16"])
    assert basis==s.Matrix(slice_report["weight_basis_24_by16"])
    slice_matrix = s.Matrix(slice_report["slice_matrix_b"])
    assert slice_matrix==s.Matrix(v_report["slice_matrix_b"])
    point = slice_matrix*s.Matrix([1,0,0,0,0,0])
    basis_rows = list(basis.T.rref()[1])
    inverse = basis[basis_rows,:].inv()
    tangent_columns = [inverse*infinitesimal_action(basis*point,a)[basis_rows,0]
                       for a in sp_basis()]
    tangent = s.Matrix.hstack(*tangent_columns,slice_matrix)
    tangent_det = tangent.det()
    assert tangent_det != 0

    polynomial_ring,*z = ring("s0,s1,s2,s3,s4,s5",QQ)
    images = [sum((polynomial_ring(QQ.from_sympy(slice_matrix[i,j]))*z[j]
                   for j in range(6) if slice_matrix[i,j]),polynomial_ring.zero) for i in range(16)]
    powers = [[value**degree for degree in range(6)] for value in images]
    cache = {}

    def restrict(polynomials):
        result = []
        for encoded in polynomials:
            value = polynomial_ring.zero
            for monomial,coefficient in encoded:
                key = tuple(monomial)
                if key not in cache:
                    monomial_value = polynomial_ring.one
                    for i,power in enumerate(monomial):
                        if power:
                            monomial_value *= powers[i][power]
                    cache[key] = monomial_value
                value += polynomial_ring(QQ.from_sympy(s.Rational(coefficient)))*cache[key]
            result.append(value)
        return result

    def decode(encoded):
        return polynomial_ring.from_dict({tuple(mon):QQ.from_sympy(s.Rational(co)) for mon,co in encoded})

    d = restrict(d_report["D3_sparse"])
    u = restrict(u_report["C5_sparse"])
    assert u==list(map(decode,slice_report["C5_sparse"]))
    v = list(map(decode,v_report["C9_sparse"]))
    c4 = decode(slice_report["c4_sparse"])
    c6 = decode(slice_report["c6_sparse"])
    mate = [1,0,3,2]
    signs = [1,-1,1,-1]

    def pair(left,right):
        result = polynomial_ring.zero
        for triple,value in zip(TRIPLES,left,strict=True):
            index = TRIPLES.index(tuple(sorted(mate[i] for i in triple)))
            sign = signs[triple[0]]*signs[triple[1]]*signs[triple[2]]
            result += multiplicity(triple)*sign*value*right[index]
        return result

    assert pair(d,u)==54*c4
    assert pair(d,v)==162*c6

    # Independent generic-injectivity witness for the linear associativity B.
    beta = form_matrices(basis*point)
    base = [s.Matrix(4,4,lambda a,j:sum(
        J[a,k]*(beta[i][j,k]+beta[j][i,k])/3 for k in range(4))) for i in range(4)]
    cubic_directions = []
    for triple in TRIPLES:
        cubic_directions.append([s.Matrix(4,4,lambda a,j:sum(
            J[a,k]*int(tuple(sorted((i,j,k)))==triple) for k in range(4))) for i in range(4)])
    b_columns = []
    for direction in cubic_directions:
        dq = s.Matrix(4,4,lambda i,j:s.trace(direction[i]*base[j]+base[i]*direction[j])/3)
        column = []
        for i,j in combinations(range(4),2):
            residual = (direction[i]*base[j]+base[i]*direction[j]
                        -direction[j]*base[i]-base[j]*direction[i]
                        +s.eye(4)[:,i]*dq[j,:]-s.eye(4)[:,j]*dq[i,:])
            column.extend(list(residual))
        b_columns.append(s.Matrix(column))
    b0 = s.Matrix.hstack(*b_columns)
    b_rows = list(b0.T.rref()[1])
    assert len(b_rows)==20
    b_minor = b0[b_rows,:].det()
    assert b_minor != 0

    # Formal identities require no specialized beta or sample root.
    x,a,b = s.symbols("X a b")
    f = x**3-a*x-b
    derivative = 3*x*x-a
    discriminant = 4*a**3-27*b*b
    inverse_numerator = 6*a*x*x-9*b*x-4*a*a
    assert s.expand(derivative*inverse_numerator-(18*a*x-27*b)*f-discriminant)==0
    assert s.expand(x*derivative-(2*a*x+3*b)-3*f)==0
    multiplication_x = s.Matrix([[0,0,b],[1,0,a],[0,1,0]])
    inverse_matrix = (6*a*multiplication_x**2-9*b*multiplication_x-4*a*a*s.eye(3))/discriminant
    assert (inverse_matrix*(3*multiplication_x**2-a*s.eye(3))-s.eye(3)).applyfunc(s.cancel)==s.zeros(3)
    assert s.cancel(s.trace(inverse_matrix))==0
    assert s.cancel(s.trace(multiplication_x*inverse_matrix))==0
    paths = [Path(__file__),D_PATH,U_PATH,SLICE_PATH,V_PATH,
             ROOT/"code/toric_degree3_scalar_20261101.py",ROOT/"code/weight1_structure.py"]
    return {"id":"TORIC_E2_SCALAR_INVERSE_20261101","all_passed":True,
            "generated_utc":datetime.now(timezone.utc).isoformat(),
            "dominance_tangent_matrix":[[str(v) for v in row] for row in tangent.tolist()],
            "dominance_tangent_determinant":str(tangent_det),
            "D3_on_slice_sparse":[[[list(mon),str(co)] for mon,co in sorted(f.items())] for f in d],
            "independent_U_restriction_checked":True,"D_pair_U_equals_54c4":True,"D_pair_V_equals_162c6":True,
            "linear_associativity_B_generic_rank":20,"B_rank_witness_rows":b_rows,"B_rank_witness_determinant":str(b_minor),
            "formal_cubic_inverse_Bezout":True,"formal_scalar_inverse":True,
            "formal_trace_inverse_derivative":0,"formal_trace_X_over_derivative":0,
            "universal_toric_barycenter":"Tr(c)=0; Tr(q_toric)=3*q_center=2*Q51+Q52",
            "proof_note":"research/notes/TORIC_E2_SCALAR_INVERSE_20261101.md",
            "dependency":"The companion report proves all96 associativity identities and existence/covariance of the unique rational V; this script independently verifies its scalar contractions and generic uniqueness witness.",
            "scope":"Inverse of the generic toric/Jacobian nonzero2torsion identification and its trace corollary; no exhaustion of the full S61 fiber.",
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

#!/usr/bin/env python3
"""Degree-three covariant scalar candidate for the toric E[2] identification.

The exported functions allow an arbitrary Poisson matrix. The certificate
proves covariance and the exact saved sample identities; a universal cubic
identity needs the separate dense-source calculation, not these samples.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from itertools import combinations_with_replacement, product
import json
import math
from pathlib import Path

import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

from weight1_structure import J, COORDS, form_matrices, infinitesimal_action, sp_basis

ROOT = Path(__file__).resolve().parents[1]
TRIPLES = list(combinations_with_replacement(range(4), 3))


def multiplicity(triple):
    return math.factorial(3)//math.prod(math.factorial(triple.count(i)) for i in range(4))


def scalar_D3(beta_matrices, poisson_matrix=J):
    """Symmetric tensor entries of D=(9/2) beta(y)(y,Pi M y)-36 beta(y)(y,Pi N y).

    beta_matrices[a][i,j] is beta_ij,a. Tensor entries, not raw cubic
    coefficients, are returned in combinations_with_replacement order.
    """
    t = beta_matrices
    pi = poisson_matrix
    nz = [(i,k,pi[i,k]) for i in range(4) for k in range(4) if pi[i,k] != 0]
    m = s.Matrix(4,4,lambda a,b:sum(t[a][i,j]*t[b][k,l]*v*w
                                 for i,k,v in nz for j,l,w in nz))
    n = s.Matrix(4,4,lambda a,b:sum(t[j][a,i]*t[l][b,k]*v*w
                                 for i,k,v in nz for j,l,w in nz))
    y = s.Matrix(s.symbols("d_y0:4"))
    beta = sum((y[i]*t[i] for i in range(4)),s.zeros(4))
    operator = pi*(s.Rational(9,2)*m-36*n)
    polynomial = s.Poly(s.expand((y.T*beta*operator*y)[0]),*y)
    return s.Matrix([polynomial.coeff_monomial(s.prod(y[i] for i in triple))/multiplicity(triple)
                     for triple in TRIPLES])


def symmetric_cubic_pair(left, right, poisson_matrix=J):
    """Full contraction D_ijk C_abc Pi_ia Pi_jb Pi_kc."""
    pi = poisson_matrix
    right_values = {triple:right[i] for i,triple in enumerate(TRIPLES)}
    total = 0
    for triple,coefficient in zip(TRIPLES,left,strict=True):
        i,j,k = triple
        total += multiplicity(triple)*coefficient*sum(
            pi[i,a]*pi[j,b]*pi[k,c]*right_values[tuple(sorted((a,b,c)))]
            for a,b,c in product(range(4),repeat=3)
            if pi[i,a] != 0 and pi[j,b] != 0 and pi[k,c] != 0)
    return s.expand(total)


def verify():
    u_path = ROOT/"reports/TORIC_DEGREE5_COVARIANT_20261101.json"
    u_report = json.loads(u_path.read_text())
    basis = s.Matrix(u_report["weight_basis_24_by16"])
    variables = s.symbols("b0:16")
    value = basis*s.Matrix(variables)
    tensor = form_matrices(value)
    d = scalar_D3(tensor)
    assert all(s.Poly(f,*variables).total_degree()==3 for f in d if f)
    basis_rows = list(basis.T.rref()[1])
    inverse = basis[basis_rows,:].inv()
    polynomial_ring,*b = ring(",".join(map(str,variables)),QQ)
    dp = [polynomial_ring.from_expr(f) for f in d]
    checks = 0
    for operator in sp_basis():
        action = s.Matrix.hstack(*(infinitesimal_action(basis[:,i],operator) for i in range(16)))
        source = inverse*action[basis_rows,:]
        assert basis*source==action
        directions = [sum((polynomial_ring(QQ.from_sympy(source[i,j]))*b[j]
                           for j in range(16) if source[i,j]),polynomial_ring.zero) for i in range(16)]
        for row,triple in enumerate(TRIPLES):
            derivative = sum((dp[row].diff(b[i])*directions[i] for i in range(16)),polynomial_ring.zero)
            target = polynomial_ring.zero
            for index,i in enumerate(triple):
                for j in range(4):
                    if operator[j,i]:
                        changed = list(triple)
                        changed[index] = j
                        target += polynomial_ring(QQ.from_sympy(operator[j,i]))*dp[TRIPLES.index(tuple(sorted(changed)))]
            assert derivative==target
            checks += 1
    sample_path = ROOT/"reports/TORIC_TWO_TORSION_SAMPLES_20261101.json"
    samples = json.loads(sample_path.read_text())["samples"]
    results = []
    x = s.Symbol("t")
    for sample in samples:
        actual_d = scalar_D3(form_matrices(list(map(s.Rational,sample["primitive_beta"]))))
        coordinates = list(map(s.sympify,sample["c_in_t"]))
        cubic = s.sympify(sample["toric_cubic"])
        scalar = s.rem(symmetric_cubic_pair(actual_d,coordinates),cubic,x)
        expected = s.sympify(sample["E2_X_in_t"])
        assert s.expand(scalar-expected)==0
        residue_u = list(map(s.Rational,sample["C5_residue"]))
        residue_v = list(map(s.Rational,sample["C9_residue"]))
        assert symmetric_cubic_pair(actual_d,residue_u)==54*s.Rational(sample["canonical_c4"])
        assert symmetric_cubic_pair(actual_d,residue_v)==162*s.Rational(sample["canonical_c6"])
        results.append({"beta_u_values":sample["beta_u_values"],"scalar_X_in_t":str(scalar),
                        "D_pair_U_equals_54c4":True,"D_pair_V_equals_162c6":True})
    paths = [Path(__file__),u_path,sample_path,ROOT/"code/weight1_structure.py"]
    return {"id":"TORIC_DEGREE3_SCALAR_20261101","all_passed":True,
            "generated_utc":datetime.now(timezone.utc).isoformat(),
            "b_variables":list(map(str,variables)),"b_in_u":u_report["b_in_u"],
            "weight_basis_24_by16":u_report["weight_basis_24_by16"],
            "D3_tensor":list(map(str,d)),
            "D3_sparse":[[[list(mon),str(co)] for mon,co in sorted(f.items())] for f in dp],
            "total_monomials":sum(map(len,dp)),"universal_Lie_covariance_checks":checks,
            "scalar_formula":"X = <D3(beta), c> using full symmetric tensor contraction with Pi^tensor3",
            "exact_sample_identifications":results,
            "scope":"Universal cubic covariance and exact two-sample normalization; the universal E[2] polynomial identity is a separate task.",
            "input_sha256":{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    report=verify()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(report,stream,indent=2);stream.write("\n")
    print(json.dumps({"all_passed":True,"universal_covariance_checks":report["universal_Lie_covariance_checks"],"output":str(args.output)}))


if __name__=="__main__":
    main()

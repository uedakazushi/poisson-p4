#!/usr/bin/env python3
"""Exact finite identities for the E beta=0 support crosswalk, not Q45 exhaustion."""
from pathlib import Path
from itertools import combinations, combinations_with_replacement
import argparse
import hashlib
import json
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [
    "legacy/v07/quadratic_classification.tex",
    "legacy/v03/p4_poisson_completeness_v3.tex",
    "legacy/v05/triple_hyperplane.tex",
    "legacy/v06/ambient_families.tex",
    "legacy/v06/nilpotent22.tex",
    "legacy/v06/zero_log_component.tex",
    "legacy/v07/nil22_v7.tex",
    "legacy/v07/minimal_thickness.tex",
    "legacy/v07/dplus_fullp4.tex",
    "legacy/v08/sections/cplus_component.tex",
    "legacy/v08/sections/cminus_component.tex",
    "legacy/v08/sections/ab_toric.tex",
    "legacy/v08/sections/u1_t16.tex",
    "legacy/v08/data/reproduced/u1_t16_certificate.json",
    "legacy/v08/data/reproduced/u1_t16_independent_audit.json",
    "legacy/v08/data/reproduced/ab_toric_certificate.json",
    "legacy/v08/data/reproduced/cplus_component_certificate.json",
    "legacy/v08/data/reproduced/independent_cminus_audit.json",
    "code/additional_rank4_families_20261103_v62.py",
    "research/notes/MPS_CATALOGUE_20261103_v62_REVIEW.md",
    "legacy/v23/proofs/literature/SS_geometry_and_MPS_identification.txt",
]
h, x, y, z, w = s.symbols("h x y z w")
V = [x, y, z, w]
PAIR = list(combinations(range(4), 2))
MON = [V[i] * V[j] for i, j in combinations_with_replacement(range(4), 2)]
J = s.zeros(4)
J[0, 1] = J[2, 3] = 1
J[1, 0] = J[3, 2] = -1


def four(values):
    out = s.zeros(4)
    for (i, j), value in zip(PAIR, values, strict=True):
        out[i, j], out[j, i] = value, -value
    return out


def pf(P):
    return s.expand(P[0, 1] * P[2, 3] - P[0, 2] * P[1, 3] + P[0, 3] * P[1, 2])


def encode(P):
    return {f"{i + 1},{j + 1}": str(s.expand(P[i, j])) for i, j in PAIR}


def coeff(P):
    polys = [s.Poly(s.expand(P[i, j]), *V) for i, j in PAIR]
    return s.Matrix([p.coeff_monomial(m) for p in polys for m in MON])


def orbit_minor(Q):
    """Ten sp4 columns plus actual Q-scaling, no linearized Jacobi upper bound."""
    cols, labels = [], []
    for i, j in combinations_with_replacement(range(4), 2):
        S = s.zeros(4)
        S[i, j] = S[j, i] = 1
        A = J * S
        assert A * J + J * A.T == s.zeros(4)
        Av = A * s.Matrix(V)
        dq = A * Q + Q * A.T
        dq -= sum((Av[k] * Q.diff(V[k]) for k in range(4)), s.zeros(4))
        cols.append(coeff(dq))
        labels.append(f"J*S{i}{j}")
    cols.append(coeff(Q))
    labels.append("Q_scalar")
    matrix = s.Matrix.hstack(*cols)
    domain_matrix = DomainMatrix.from_Matrix(matrix, extension=True).to_field()
    _, col_ids = domain_matrix.rref()
    assert len(col_ids) == 10
    selected = domain_matrix.extract(range(60), list(col_ids))
    _, row_ids = selected.transpose().rref()
    small = domain_matrix.extract(list(row_ids), list(col_ids))
    det = small.domain.to_sympy(small.det())
    assert det != 0
    return dict(group="Sp4 x Gm", rows=list(row_ids), columns=list(col_ids),
                determinant=str(det), labels=labels,
                original60by11_matrix=[[str(v) for v in row] for row in matrix.tolist()],
                exact_coefficient_field=str(domain_matrix.domain),
                upper_dimension_dependency="this Q45 component has dimension10 by Q45_22; its group orbit is contained in it")


def main(output):
    begun = time.monotonic()
    A, B, D, E = s.symbols("A B D E")
    def N(a, b, c, d, e, p=0, r=0):
        return four([a*x*y+b*z*w, e*x*z, (c-e)*x*w+p*y*y,
                     (-b-e)*y*z+r*x*x, (b-c+e)*y*w, c*x*y+d*z*w])
    nonnil = {
        "G1": N(A, 0, 0, D, E), "G2": N(A, B, 0, B, E),
        "T1_Q45": N(1, 0, 0, 3, 1, 1),
        "T2_Q45": N(s.Rational(1,2), 0, s.Rational(1,2), 3, 1, 1),
        "T3_Q45": N(1, s.Rational(3,2), 0, s.Rational(3,2), -s.Rational(1,2), 1),
        "T4_Q45": N(-1, 0, 2, 3, 1, 1, 3),
        "E1": four([x*y, x*z, -2*x*y-x*w, -y*z, y*y+y*w, z*w]),
        "E2": four([x*y/2, x*z, -2*x*y-x*w/2, -y*z, y*y+y*w/2, x*y/2+z*w]),
        "E3": four([x*y, x*z/2, -x*w/2+y*z, -y*z/2, y*w/2, z*w]),
        "E4": four([x*y/2, x*z/2, y*z, -y*z/2, 0, x*y/2+z*w]),
    }
    records = []
    expected_kappa = [A*D, A*B, 3, s.Rational(3,2), s.Rational(3,2), -3, 1, s.Rational(1,2), 1, s.Rational(1,2)]
    for (name, Q), kappa in zip(nonnil.items(), expected_kappa, strict=True):
        weights = ({"G2": [2, 0, 0, 1, 1], "T3_Q45": [3, 0, 2, 3, 0]}
                   .get(name, [3, 2, 2, 1, 0]))
        P = h*h*J+Q
        leading = s.zeros(4)
        terms = []
        for i,j in PAIR:
            for mon,c in s.Poly(s.expand(P[i,j]),h,*V).terms():
                if c == 0:
                    continue
                weight = sum(m*d for m,d in zip(mon,weights))-weights[i+1]-weights[j+1]
                assert weight >= 0, (name, (i,j), mon, c, weight)
                if weight == 0:
                    assert mon[i+1] == mon[j+1] == 1 and sum(mon) == 2
                    leading[i,j] += c*V[i]*V[j]
                    leading[j,i] = -leading[i,j]
                terms.append(dict(pair=[i+1,j+1], monomial=list(mon),coefficient=str(c),weight=weight))
        assert s.expand(pf(leading)-kappa*x*y*z*w)==0
        assert kappa != 0
        records.append(dict(Q45_component=name, polynomial_family="weighted pullback of h^2J+Q, valid for epsilon!=0; nonnegative exponents extend to0",
                            pullback_weights=weights, all_terms=terms, limit=encode(leading),
                            limit_Pfaffian=str(kappa*x*y*z*w),
                            projective_limit_degeneracy_divisor=str(kappa*h*x*y*z*w),
                            generic_parameter_open=str(kappa)+" != 0",
                            known_target="union of the40 MPS components via a component containing the orbit closure and its SNC limit"))

    # Whole9-parameter minimal nilpotent source, with exact potential formulas.
    t = s.symbols("t1:10")
    minimal = []
    for c, target in [(s.Rational(0),"R4_CA"),(s.Rational(1,4),"R4_CB")]:
        q23 = t[0]*x*x+(-2*c-t[5])*x*z-t[7]*z*z/2-2*t[8]*z*w+t[1]*x*w+t[2]*w*w
        q24 = t[3]*x*x+t[4]*x*z+t[5]*x*w+t[6]*z*z+t[7]*z*w+t[8]*w*w
        H = (t[3]*x*x*z+t[4]*x*z*z/2+(t[5]+c)*x*z*w+t[6]*z**3/3
             +t[7]*z*z*w/2+t[8]*z*w*w-t[0]*x*x*w-t[1]*x*w*w/2-t[2]*w**3/3)
        source = h*h*J+four([(s.Rational(1,2)-c)*x*x,0,0,q23,q24,c*x*x])
        f = h*h+(s.Rational(1,2)-c)*x*x
        extra = 0 if c==0 else -s.diff(f,x)/2
        target_P = four([f,0,0,extra*z-s.diff(H,w),extra*w+s.diff(H,z),h*h if c==0 else f])
        assert all(s.expand(v)==0 for v in source-target_P)
        minimal.append(dict(Q45_component="M0" if c==0 else "M1", known_target=target,
                            all9_parameters_retained=True,c=str(c),potential=str(H),f=str(f),all6_coefficient_differences_zero=True))

    # Scalar/Sp4 dense orbit proof for two representatives whose full image
    # would not follow from a single point without the dimension argument.
    q_open = four([-x*x+x*z-x*w,-x*x,0,3*x*y-4*x*z+z*z-3*z*w,
                   x*x+x*y-2*x*z+2*x*w+z*z/2-2*z*w+w*w,x*x+x*w])
    q_a0 = four([3*x*z/4,0,z*z+3*x*x/4,-z*z/4,
                 x*x+2*x*z+3*z*z-3*x*y/4+z*w/4,x*z/4])
    q_m = four([3*x*x-6*x*z,0,-2*x*x,3*x*z-18*z*z,
                6*x*x-2*x*y-3*x*w+30*z*w,6*x*z])
    q_cplus = four([x*x+8*x*z-z*z,0,5*x*x-4*x*z-z*z,
                    x*x+4*x*z-5*z*z,x*x+2*(z-x)*(y+w),-x*x+8*x*z+z*z])
    rr = s.sqrt(3)
    cc, dd, uu, vv = rr/18, s.Rational(1,2)+rr/18, s.Rational(1,6), -s.Rational(1,36)
    q_dpm = four([(1-dd)*x*z-vv*x*x-uu*z*z,0,
                 -3*cc*z*z-4*uu*x*z+(cc-3*dd+2)*x*x/2,
                 -cc*x*x/12-4*vv*x*z-(cc+dd)*z*z/2,
                 x*x+2*x*z+3*z*z+(2*dd-cc-1)*x*y+cc*z*w+2*uu*y*z+2*vv*x*w,
                 dd*x*z+vv*x*x+uu*z*z])
    dense = []
    for name,target,Q in [("nil22_open","R4_H_QUARTIC_ESSENTIAL5",q_open),
                          ("A0","R4_A0",q_a0),("M","R4_MINIMAL_M",q_m),
                          ("Cplus","R4_CPLUS_EXT",q_cplus),("Dpm","R4_DPM",q_dpm)]:
        dense.append(dict(Q45_component=name,known_target=target,representative=encode(Q),orbit_certificate=orbit_minor(Q)))

    body=dict(task_id="E_zero",scope="new exact identities for the beta0 crosswalk; Q45_22 exhaustion and known componenthood remain explicit input theorems",
              all_passed=True,nonnil_ten=records,whole_minimal_two=minimal,dense_orbits=dense,
              no_new_global_exhaustion_or_component_count=True,
              input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES},
              code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.monotonic()-begun)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open("x") as f:json.dump(body,f,ensure_ascii=False,indent=2);f.write("\n")
    print(json.dumps(dict(all_passed=True,nonnil=10,whole_minimal=2,dense_orbits=[(r['Q45_component'],r['orbit_certificate']['determinant'])for r in dense],seconds=body['seconds'])))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args();main(args.output)

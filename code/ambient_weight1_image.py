#!/usr/bin/env python3
"""Exact first-jet image certificates for the existing CA and CB components.

Their complex set-theoretic primitive first-jet images equal the same
nine-dimensional closed irreducible subset. No scheme-image claim is made.
The dimension upper bound uses the line-incidence proof
in research/notes/AMBIENT_WEIGHT1_IMAGE_20261030.md, not sampled ranks alone.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import sympy as s

from weight1_structure import (COORDS,J,PAIRS,coefficient_map,constraints,
                               form_matrices,infinitesimal_action,
                               minor_witness,sp_basis)

ROOT=Path(__file__).resolve().parents[1]


def flatten_form(matrix,y):
    return s.Matrix([s.diff(matrix[i,j],y[k]) for (i,j),k in COORDS])


def primitive(linear,y):
    tau=linear[0,1]+linear[2,3]
    alpha=s.Matrix([s.diff(tau,z)/5 for z in y])
    v=J*alpha
    return s.expand(linear-2*alpha.dot(y)*J+y*v.T-v*y.T)


def common_plane():
    y=s.Matrix(s.symbols("x y z w"));x,_,z,w=y
    t,A,B,C,D,F=s.symbols("t A B C D F")
    params=[t,A,B,C,D,F]
    h=A*x*z+B*x*w+C*z*z+D*z*w+F*w*w
    common=s.zeros(4)
    common[0,1]=t*x;common[2,3]=-t*x
    common[1,2]=-s.diff(h,w)+t*z/2
    common[1,3]=s.diff(h,z)+t*w/2
    for i,j in PAIRS:common[j,i]=-common[i,j]
    for kind,p in [("A",5*t/2),("B",-5*t)]:
        raw=s.zeros(4)
        raw[0,1]=p*x
        raw[2,3]=0 if kind=="A" else p*x
        raw[1,2]=-s.diff(h,w)-(0 if kind=="A" else p*z/2)
        raw[1,3]=s.diff(h,z)-(0 if kind=="A" else p*w/2)
        for i,j in PAIRS:raw[j,i]=-raw[i,j]
        assert primitive(raw,y)==common
    vector=flatten_form(-J*common*J,y)
    assert s.Matrix.vstack(*constraints())*vector==s.zeros(8,1)
    return params,vector


def intrinsic_fiber_constraints():
    """K_beta(e2)=0 and i_e2 beta in span(x dx), plus closed/primitive."""
    unknown=s.Matrix(s.symbols("b0:24"))
    matrices=form_matrices(unknown)
    equations=list(s.Matrix.vstack(*constraints())*unknown)
    equations += [matrices[1][i,j] for i,j in PAIRS]
    equations += [matrices[k][1,j] for k in range(4) for j in range(4)
                  if (k,j)!=(0,0)]
    matrix,_=s.linear_eq_to_matrix(equations,list(unknown))
    return matrix


def general_constant_frame_check():
    a,d,e,g=s.symbols("a d e g",nonzero=True)
    lam=s.Matrix([[0,a,0,0],[-a,0,d,e],[0,-d,0,g],[0,-e,-g,0]])
    frame=s.Matrix([[1,0,0,0],[0,1/a,0,0],
                    [d/a,0,1,0],[e/(a*g),0,0,1/g]])
    assert (frame*lam*frame.T-J).applyfunc(s.cancel)==s.zeros(4)
    y=s.Matrix(s.symbols("x y z w"));x,_,z,w=y
    f1,g1,ux,uz,uw,vx,vz=s.symbols("f1 g1 ux uz uw vx vz")
    raw=s.zeros(4)
    raw[0,1]=f1*x;raw[2,3]=g1*x
    raw[1,2]=ux*x+uz*z+uw*w
    # Degree-zero Jacobi: a*g1+g*(uz+vw)=0.
    raw[1,3]=vx*x+vz*z+(-a*g1/g-uz)*w
    for i,j in PAIRS:raw[j,i]=-raw[i,j]
    shifted=raw.subs(dict(zip(y,frame.inv()*y)),simultaneous=True)
    linear=(frame*shifted*frame.T).applyfunc(s.cancel)
    assert s.cancel(linear[0,1]-f1*x/a)==0
    assert s.cancel(linear[2,3]-g1*x/g)==0
    normalized=primitive(linear,y).applyfunc(s.cancel)
    vector=flatten_form(-J*normalized*J,y)
    intrinsic=intrinsic_fiber_constraints()
    assert (intrinsic*vector).applyfunc(s.cancel)==s.zeros(intrinsic.rows,1)
    return {"constant_matrix_order":"a=P12(0), d=P23(0), e=P24(0), g=P34(0)",
            "open_condition":"a*g != 0",
            "frame":[[str(x) for x in row] for row in frame.tolist()],
            "primitive_slope":"(2*f1/a-3*g1/g)/5"}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"reports/AMBIENT_WEIGHT1_IMAGE_20261030.json")
    args=parser.parse_args()
    params,vector=common_plane()
    basis=vector.jacobian(params)
    intrinsic=intrinsic_fiber_constraints()
    assert intrinsic.rank()==18 and basis.rank()==6
    assert intrinsic*basis==s.zeros(intrinsic.rows,6)

    lie=sp_basis()
    # Lie algebra of the line parabolic preserving C e2 has dimension seven.
    line_action=s.Matrix.hstack(*(s.Matrix([a[i,1] for i in (0,2,3)]) for a in lie))
    stabilizer_parameters=line_action.nullspace()
    assert len(stabilizer_parameters)==7
    for coeff in stabilizer_parameters:
        a=sum((c*b for c,b in zip(coeff,lie)),s.zeros(4))
        for col in range(6):
            assert intrinsic*infinitesimal_action(basis[:,col],a)==s.zeros(intrinsic.rows,1)

    sample=[2,3,5,7,11,13]
    value=vector.subs(dict(zip(params,sample)))
    action=s.Matrix.hstack(*(infinitesimal_action(value,a) for a in lie))
    differential=basis.row_join(action)
    assert differential.rank()==9
    assert action.rank()==8
    assert coefficient_map(value).rank()==3
    assert coefficient_map(value).nullspace()==[s.Matrix([0,1,0,0])]
    simple_parameters=[1,0,0,s.Rational(1,2),s.Rational(1,2),s.Rational(1,2)]
    simple=vector.subs(dict(zip(params,simple_parameters)))
    simple_minor=coefficient_map(simple).extract([0,1,2],[0,2,3])
    assert simple_minor.det()==1
    frame_report=general_constant_frame_check()
    report={
        "id":"AMBIENT_WEIGHT1_IMAGE_20261030",
        "timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "all_passed":True,"field":"QQ",
        "parameter_order":[str(x) for x in params],
        "common_primitive_linear_bivector":"t*x*(d12-d34)+d2 wedge (X_h2+(t/2)*(z*d3+w*d4))",
        "potential":"h2=A*x*z+B*x*w+C*z^2+D*z*w+F*w^2",
        "original_slopes":{"CA":"p=5*t/2","CB":"p=-5*t"},
        "common_primitive_form_coefficients":[str(x) for x in vector],
        "intrinsic_fiber_constraint_rank":18,
        "intrinsic_fiber_dimension":6,
        "line_parabolic_dimension":7,
        "line_parabolic_preserves_fiber":True,
        "incidence_base_dimension":3,
        "geometric_image_dimension_upper_bound":9,
        "sample_parameters":sample,
        "sample_combined_differential_witness":minor_witness(differential),
        "sample_sp4_orbit_rank":8,
        "sample_coefficient_rank_witness":minor_witness(coefficient_map(value)),
        "sample_coefficient_kernel":[0,1,0,0],
        "simple_rank_three_witness":{
            "parameters":[str(x) for x in simple_parameters],
            "form":"x dx wedge dy+(z+w)dx wedge dz+w dx wedge dw-x dz wedge dw",
            "coefficient_map_rows_zero_based":[0,1,2],
            "coefficient_map_columns_zero_based":[0,2,3],
            "minor":[[str(x) for x in row] for row in simple_minor.tolist()],
            "determinant":"1"},
        "arbitrary_nonzero_constant_frame":frame_report,
        "claim":"The complex set-theoretic primitive first-jet images of CA and CB equal the same closed irreducible reduced subset Z of dimension 9; Z is contained in coefficient rank <= 3, hence proper in V16.",
        "image_scope":"Complex points and reduced closed subsets only; no scheme-theoretic image or fiber structure is asserted.",
        "dimension_dependency":"Using bundled component dimension28 and S61 dimension difference9, the associated19dim slice components have generic fiber dimension10 over Z.",
        "limitations":["Does not identify or classify the generic10dim fibers.",
                       "Does not determine CA intersection CB; equal firstjets do not imply equal full structures.",
                       "No new irreducible component is counted."],
        "proof_note":"research/notes/AMBIENT_WEIGHT1_IMAGE_20261030.md",
        "source_sha256":{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in [Path(__file__),ROOT/"code/weight1_structure.py",ROOT/"legacy/v06/ambient_families.tex"]},
        "reproduction_command":"PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_weight1_image.py",
    }
    with args.output.open("x") as stream:
        json.dump(report,stream,ensure_ascii=False,indent=2);stream.write("\n")
    print(json.dumps({"all_passed":True,"image_dimension":9,
                      "sample_minor":report["sample_combined_differential_witness"]["determinant"],
                      "report":str(args.output)}))


if __name__=="__main__":main()

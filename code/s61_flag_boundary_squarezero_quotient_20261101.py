#!/usr/bin/env python3
"""A closed flat square-zero length-five quotient at every flag with u5 nonzero."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
from s61_flag_transverse_obstruction_20261101 import specialize

ROOT=Path(__file__).resolve().parents[1]


def verify():
    start=time.time();model=S61();u=s.symbols('u1:62')
    candidate='code/data/s61_flag_squarezero_kernel_20261101.json'
    source='reports/S61_FLAG_UNIFORM_OPEN_20261101.json'
    K=s.Matrix(json.loads((ROOT/candidate).read_text())['kernel']).applyfunc(s.sympify)
    assert K.shape==(45,8)
    one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    fixed=two.index(37);assert fixed==25
    beta={i:u[i] for i in one};beta[3]=0
    _,equations=specialize(model,beta,True)
    J=s.Matrix([[f.get(tuple(sorted((fixed,j))),0)*(2 if j==fixed else 1)
                 for j in range(45)] for f in equations])
    assert all(not f.get((fixed,fixed),0) for f in equations)
    assert (J*K).applyfunc(s.cancel)==s.zeros(135,8)
    for entry in K:
        denominator=s.denom(s.cancel(entry))
        assert denominator.free_symbols<={u[4]}
        assert len(s.Poly(denominator,u[4]).terms())==1
    kernel_rows=[0,19,20,21,25,26,35,38]
    KM=K[kernel_rows,:].to_DM().to_field()
    kernel_minor=s.factor(KM.domain.to_sympy(KM.det()))
    assert kernel_minor==s.Rational(80,21)*u[4]**4
    radial=s.zeros(45,1);radial[fixed]=1
    assert K[:,0]==radial
    z,bs,B=primitive_beta(model)
    vec=s.Matrix([B[i,j].diff(z[k]) for (i,j),k in COORDS])
    basis=vec.jacobian(bs);rows=list(basis.T.rref()[1]);iv=basis[rows,:].inv()
    action=s.Matrix.hstack(*[iv*infinitesimal_action(vec,A)[rows,:] for A in sp_basis()])
    action=action.subs(dict(zip(bs,[beta[i] for i in one])))
    normal=action[one.index(3),:]
    expected=s.Matrix([[0,0,0,0,-2*u[25],-u[4],-u[13]+2*u[48],0,0,-u[11]]])
    assert normal==expected
    pivot=5;allowed=s.zeros(10,9);columns=[j for j in range(10) if j!=pivot]
    for col,j in enumerate(columns):
        allowed[j,col]=1;allowed[pivot,col]=-normal[j]/normal[pivot]
    assert normal*allowed==s.zeros(1,9)
    old=json.loads((ROOT/source).read_text())['flag']
    O0=s.Matrix(old['q_orbit_differential']).applyfunc(s.Rational)
    O=s.zeros(45,10);chart=[j for j in range(45) if j!=fixed]
    for i,j in enumerate(chart):O[j,:]=O0[i,:]
    OC=(O*allowed).applyfunc(s.cancel)
    selected=list(OC.to_DM().to_field().rref()[1]);assert len(selected)==3
    tangent=OC[:,selected]
    field=s.QQ.frac_field(*sorted(K.row_join(tangent).free_symbols,key=str))
    coeff=(KM.convert_to(field).inv()*tangent[kernel_rows,:].to_DM().convert_to(field)).to_Matrix().applyfunc(s.cancel)
    target=s.zeros(8,3);target[1,0]=-2;target[4,1]=1;target[2,2]=2/u[4]
    assert coeff==target and (K*coeff-tangent).applyfunc(s.cancel)==s.zeros(45,3)
    radial8=s.zeros(8,1);radial8[0]=1
    normal_columns=[3,5,6,7]
    complement=radial8.row_join(coeff).row_join(s.eye(8)[:,normal_columns])
    complement_minor=s.factor(complement.det());assert complement_minor==-4/u[4]
    N=K[:,normal_columns]
    assert N[fixed,:]==s.zeros(1,4)
    assert (J*N).applyfunc(s.cancel)==s.zeros(135,4)
    # Every leading equation is quadratic in q: the verified value and
    # derivative imply its exact vanishing in the square-zero algebra.
    enc=lambda a:[[str(x) for x in row] for row in a.tolist()]
    paths=[str(Path(__file__).relative_to(ROOT)),candidate,source,
        'code/s61_model.py','code/primitive_jet_elliptic_curve_20261101.py',
        'code/weight1_structure.py','code/s61_flag_transverse_obstruction_20261101.py',
        'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',
        'reports/S61_FLAG_ZERO12_RANK_FIVE_20261101.json',
        'reports/S61_FLAG_TRUE_CORE_RANK_FIVE_20261101.json',
        'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',
        'research/notes/S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101.md']
    return dict(id='S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101',all_passed=True,
        open='u4=0,u5!=0 at the canonical projective flag',
        all_135_leading_values_and_1080_derivative_coefficients_zero=True,
        kernel=enc(K),kernel_minor_rows=kernel_rows,kernel_minor=str(kernel_minor),
        denominators_only_powers_of_u5=True,incidence_normal=list(map(str,normal)),
        compensated_sp_directions=enc(allowed),actual_C_tangent_selected_columns=selected,
        actual_C_tangent=enc(tangent),actual_C_kernel_coefficients=enc(coeff),
        radial_C_normal_complement_minor=str(complement_minor),
        squarezero_normal_columns=normal_columns,squarezero_normal_directions=enc(N),
        fixed_beta_closed_flat_quotient='k[[z0,z1,z2]][epsilon0,epsilon1,epsilon2,epsilon3]/(epsilon)^2',
        homogenizer_in_quotient_zero=True,
        finite_rank_five_sample_consequence='The kernel over the same regular three-parameter base is finite torsion and has at most44 generators.',
        limitations=['This constructs an actual closed quotient, not an isomorphism with the ambient completion.',
                     'No ambient homogenizer nilpotence, exact kernel presentation or whole flag saturation removal follows.'],
        proof_note='research/notes/S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101.md',
        input_sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},
        seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=verify()
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print('FLAG CLOSED SQUAREZERO QUOTIENT PASS',result['seconds'])

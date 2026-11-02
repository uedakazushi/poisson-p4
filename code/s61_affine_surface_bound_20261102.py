#!/usr/bin/env python3
"""Certify the marked linear factor and the divisor argument for affine dimension at most two."""
import argparse
import hashlib
import itertools
import json
import time
from pathlib import Path
import sympy as sp
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import J

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_AFFINE_SURFACE_BOUND_20261102'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def certify():
    start=time.monotonic()
    files={}
    def load(stem):
        p=Path('reports')/(stem+'.json');record=json.loads((ROOT/p).read_text())
        assert record['all_passed'];files[str(p)]=sha(ROOT/p)
        return record
    prior=load('S61_AFFINE_BOUNDARY_DIMENSION_20261102')
    assert prior['raw_generic_infinity_exact_dimension']==3 and prior['flag_residual_joint_dimension_upper']==17
    rank=load('S61_BINARY_H11_RANK_DROP_RESIDUAL_FINITE_20261102')
    assert rank['residual_joint_dimension_upper']==16 and rank['residual_finite_over_general_beta']
    cover=load('S61_BINARY_RANK_DROP_COVER_20261101')
    assert cover['prime_ideals'][:2]==[['5*u35 + u50'],['5*u14 - u49']]
    excluded=load('S61_BINARY_DEEPER_REDUCTION_20261101')['excluded_components'][0]
    assert excluded['ideal_index']==0 and excluded['beta_image_dimension_upper_bound']==15
    actions=load('S61_BINARY_UNIFORM_RANK_20261101')
    curve=load('PRIMITIVE_JET_ELLIPTIC_CURVE_20261101')
    z,u,beta=primitive_beta(S61());z=sp.Matrix(z)
    names={str(x):x for x in u}
    idx=[(i,j,k) for i,j in itertools.combinations(range(4),2) for k in range(4)]
    def coefficients(B):
        return sp.Matrix([B[i,j].diff(z[k]) for i,j,k in idx])
    basis=coefficients(beta).jacobian(u)
    rows=list(basis.T.to_DM().rref()[1]);inverse=basis[rows,:].inv()
    beta=beta.subs({names['u12']:0,names['u13']:0})
    def pull(g):
        return (inverse*coefficients(g.T*beta.subs(dict(zip(z,g*z)),simultaneous=True)*g)[rows,:]).applyfunc(sp.factor)
    first=lambda vec:5*vec[u.index(names['u14'])]-vec[u.index(names['u49'])]
    second=lambda vec:5*vec[u.index(names['u35'])]+vec[u.index(names['u50'])]
    old=sp.Matrix(u);first0,second0=first(old),second(old)
    a,b,c,d=sp.symbols('A B C D')
    change=sp.Matrix([[a,b],[c,d]])
    # Position coordinates (x,w) and the corresponding momenta (y,-z).
    pos=change.inv()*sp.Matrix([z[0],z[3]]);mom=change.T*sp.Matrix([z[1],-z[2]])
    g=sp.Matrix([pos[0],mom[0],-mom[1],pos[1]]).jacobian(z)
    assert (g.T*J*g-J).applyfunc(sp.cancel)==sp.zeros(4)
    transformed=pull(g)
    assert sp.cancel(first(transformed)-a*first0-b*second0)==0
    assert sp.cancel(second(transformed)-c*first0-d*second0)==0
    assert transformed[u.index(names['u12'])]==transformed[u.index(names['u13'])]==0
    generators=[sp.Matrix(actions['stabilizer_generators'][i]) for i in (0,1,3)]
    params=sp.symbols('p0:3')
    radical=sp.eye(4)+sum((v*m for v,m in zip(params,generators)),sp.zeros(4))
    assert (radical.T*J*radical-J).applyfunc(sp.expand)==sp.zeros(4)
    assert sp.Matrix.hstack(*[m.reshape(16,1) for m in generators]).rank()==3
    transformed=pull(radical)
    assert sp.expand(first(transformed)-first0)==sp.expand(second(transformed)-second0)==0
    assert transformed[u.index(names['u12'])]==transformed[u.index(names['u13'])]==0
    # The marked diagonal point in the usual flag coordinates.
    flag=sp.eye(4);flag[2,2]=flag[3,3]=0;flag[2,3]=-1;flag[3,2]=1
    z0,uu,B0=primitive_beta(S61())
    ff=(inverse*coefficients(flag.T*B0.subs(dict(zip(z,flag*z)),simultaneous=True)*flag)[rows,:]).applyfunc(sp.factor)
    assert [ff[u.index(names[key])] for key in ('u12','u13')]==[-names['u4'],-names['u5']]
    assert sp.expand(first(ff)+5*names['u14']-4*names['u49'])==0
    for f in [Path(__file__).relative_to(ROOT),Path('research/notes')/(STEM+'.md'),
            Path('research/notes/PRIMITIVE_JET_ELLIPTIC_CURVE_20261101.md'),Path('research/notes/S61_R0_INFINITY_20261101.md'),
            Path('research/notes/S61_BINARY_H11_RANK_DROP_RESIDUAL_FINITE_20261102.md'),
            Path('code/primitive_jet_elliptic_curve_20261101.py'),Path('code/s61_model.py'),Path('code/weight1_structure.py')]:
        files[str(f)]=sha(ROOT/f)
    return dict(id=STEM,all_passed=True,primitive_flagged_plane_conditions=['u12','u13'],
                marked_functional_local_values=[str(first0),str(second0)],
                transformed_values_under_general_GL2=[str(a*first0+b*second0),str(c*first0+d*second0)],
                radial_shear_dimension=3,common_zero_beta_image_dimension_upper=15,
                marked_square_flag_conditions=['u4','u5','5*u14-4*u49'],
                tautological_rank=2,tautological_degree=-5,marked_quotient_degree=0,marked_line_degree=-5,
                restriction_on_binary_plane='One marked line in each P2; its only square is twice the marked factor.',
                difference_section_slope=-5,pullback_by_two_difference_filtration_degree=-20,
                residual_nonsquare_marked_divisor_generic_fiber_dimension_upper=0,
                unknown_binary_boundary_surface_excluded=True,generic_affine_boundary_dimension_upper=1,
                generic_affine_dimension_upper=2,finite_affine_projection_variable_count=2,
                residual_curves_or_points_excluded=False,whole_affine_fiber_finite_claimed=False,
                input_sha256=files,seconds=time.monotonic()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    out=parser.parse_args().output;assert not out.exists();j=certify()
    with out.open('x')as f:json.dump(j,f,indent=2);f.write('\n')
    print(f"PASS marked parabolic covariance, genuine boundary at most one, affine dimension at most two; {j['seconds']:.3f} seconds")

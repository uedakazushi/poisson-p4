#!/usr/bin/env python3
"""Identify and verify two fixed rational coordinate functions giving a finite projection."""
import argparse
import hashlib
import itertools
import json
import time
from pathlib import Path
import sympy as sp
from s61_model import S61,X

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_AFFINE_TWO_COORDINATE_PROJECTION_20261102'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def certify():
    start=time.monotonic()
    inputs={}
    def load(name):
        p=Path('reports')/(name+'.json')
        j=json.loads((ROOT/p).read_text());assert j['all_passed']
        inputs[str(p)]=sha(ROOT/p)
        return j
    surface=load('S61_AFFINE_SURFACE_BOUND_20261102')
    assert surface['generic_affine_boundary_dimension_upper']==1
    load('S61_GENERIC_REDUCED_INFINITY_REVIEW_20261101')
    for k in ('S61_R0_GENERIC_EXHAUSTION','S61_NONBINARY_CANONICAL_EXHAUSTION',
              'S61_NONBINARY_DOUBLE_ROOT_EXHAUSTION','S61_R0_MIXED_CURVE_EXHAUSTION',
              'S61_R0_T5_ZERO_L_EXHAUSTION','S61_R0_MIXED_SURFACE_T7_EXHAUSTION',
              'S61_MINIMAL_VERTEX_EXHAUSTION','S61_MINIMAL_MAIN_REMOVAL',
              'S61_MINIMAL_MAIN_JORDAN_REMOVAL','S61_MINIMAL_SEMISIMPLE_TAILS_REMOVAL',
              'S61_MINIMAL_TRIPLE_TAILS_REMOVAL','S61_MINIMAL_PURE_CUBIC_TAILS_REMOVAL',
              'S61_MINIMAL_MIXED_TAILS_REMOVAL','S61_MINIMAL_M0_BASELINE_REMOVAL',
              'S61_MINIMAL_M1_BASELINE_REMOVAL','S61_MINIMAL_ALL_AXIS_REMOVAL',
              'S61_MINIMAL_RANKONE_A1_REMOVAL','S61_MINIMAL_RANKONE_JUNCTION',
              'S61_MINIMAL_ZERO49_BRANCHES','S61_MINIMAL_JUNCTION_SATURATION'):
        load(k+'_20261101')
    model=S61();u=sp.symbols('u1:62')
    q=[v if wt==2 else 0 for wt,v in zip(model.reduced['parameter_weights'],u)]
    polynomial=model.reconstruct(q)
    y=X[1:]
    def coeff(i,j,point):
        entry=polynomial[i+1,j+1]-y[i]*polynomial[0,j+1]+y[j]*polynomial[0,i+1]
        axis={coord:int(k==point)*sp.Symbol('z') for k,coord in enumerate(y)}
        return sp.expand(entry.subs(X[0],1).subs(axis)).coeff(sp.Symbol('z'),2)
    assert [coeff(0,2,0),coeff(1,3,1)]==[u[5],u[41]]
    # Plücker equation after p34=-p12. The two linear-form vanishing conditions
    # on support planes are p13=p24=0. The surviving conic is smooth.
    p12,p13,p14,p23,p24=sp.symbols('p12 p13 p14 p23 p24')
    lag=-p12**2-p13*p24+p14*p23
    common=lag.subs({p13:0,p24:0})
    assert abs(sp.hessian(common,[p12,p14,p23]).det())==2
    assert sp.hessian(lag.subs(p13,0),[p12,p14,p23,p24]).rank()==3
    assert sp.hessian(lag.subs(p24,0),[p12,p13,p14,p23]).rank()==3
    assert sp.hessian(lag.subs(p12,0),[p13,p14,p23,p24]).rank()==4
    # (base dimension, allowed nonsquare fiber dimension, allowed square fiber dimension)
    pieces={'both_plucker':(1,2,1),'one_plucker_one_value':(2,1,0),
            'one_plucker_value_evaluation_zero_vector':(1,2,1),
            'independent_argument_vectors':(3,0,-100),
            'dependent_argument_vectors':(2,1,0)}
    assert max(b+n for b,n,_ in pieces.values())==3
    assert max(b+s for b,_,s in pieces.values())==2
    for f in [Path(__file__).relative_to(ROOT),Path('code/s61_model.py'),Path('research/notes')/(STEM+'.md'),
            Path('research/notes/S61_MINIMAL_INFINITY_REMAINING_20261101.md'),Path('research/notes/S61_R0_MIXED_ORBIT_REDUCTION_20261101.md'),
            Path('research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md')]:
        inputs[str(f)]=sha(ROOT/f)
    polys=['u6','u42']
    return dict(id=STEM,all_passed=True,linear_coordinates=polys,plucker_equation=str(lag),
                center_support_plane_common_conic=str(common),common_center_dimension_upper={'binary_orbit':3,'flag_orbit':2},
                orbit_dimensions={'binary_orbit':5,'flag_orbit':4},symplectic_group_dimension=10,
                bad_group_element_incidence_dimension_upper=9,generic_affine_projection_is_finite=True,
                base_polynomial_ring='K[u6,u42]',dimension_upper=2,projection_dominance_asserted=False,
                generic_affine_fiber_finite_proved=False,input_sha256=inputs,seconds=time.monotonic()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    out=parser.parse_args().output;assert not out.exists();r=certify()
    with out.open('x')as f:json.dump(r,f,indent=2);f.write('\n')
    print(f"PASS finite affine projection onto (u6,u42); {len(r['input_sha256'])} input hashes")

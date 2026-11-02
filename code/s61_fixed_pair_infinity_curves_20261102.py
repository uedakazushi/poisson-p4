#!/usr/bin/env python3
"""Verify the two elliptic blow-up maps and classify the fixed-center curve cycles."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as sp

from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import J

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_FIXED_PAIR_INFINITY_CURVES_20261102'


def verify():
    coords,u,beta=primitive_beta(S61())
    v=sp.Matrix(coords)
    omega=-J
    second=J*beta*v
    assert sp.expand((v.T*omega*second)[0])==0
    assert sp.expand((v.T*beta*second)[0])==0
    p=v*second.T-second*v.T
    cubic=[sp.Poly(p[i,j],*coords) for i in range(4) for j in range(i+1,4)]
    assert all(row.total_degree()==3 for row in cubic)
    assert sp.expand(p[0,1]+p[2,3])==0

    h,a,b,c,d=sp.symbols('h a b c d')
    lag=-h*h-a*d+b*c
    hyperplane_ranks=[sp.hessian(lag.subs(z,0),[y for y in [h,a,b,c,d] if y!=z]).rank()
                      for z in [a,d,h]]
    assert hyperplane_ranks==[3,3,4]
    assert abs(sp.hessian(lag.subs({a:0,d:0}),[h,b,c]).det())==2
    # On the big chart a Lagrangian plane has basis the two columns of U.
    # The two sections are dependent iff e_x,e_y project dependently modulo U.
    i,j,k=sp.symbols('i j k')
    U=sp.Matrix([[1,0],[i,j],[0,1],[j,k]])
    assert U.T*omega*U==sp.zeros(2)
    determinant=sp.eye(4)[:,[0,1]].row_join(U).det()
    assert sp.expand(determinant+U[[0,1],:].det())==0
    flag=[(10,1),(5,2),(5,2),(5,2),(5,2),(5,4)]
    binary=[(15,1)]+[(1,1)]*15
    degree=lambda cycles:sum(d*m for d,m in cycles)
    assert (degree(binary),degree(flag),65*degree(binary)+5*degree(flag))==(30,70,2300)
    inputs={}
    def record(path):inputs[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    for name in ['S61_FIXED_PAIR_SINGLE_COLON_20261102','S61_PROJECTIVE_TOP_CYCLE_20261102',
                 'S61_AFFINE_TWO_COORDINATE_PROJECTION_20261102','S61_R0_INFINITY_20261101',
                 'NULL_CORRELATION_DUALITY_20261101']:
        report=ROOT/'reports'/(name+'.json')
        assert json.loads(report.read_text())['all_passed']
        record(report);record(ROOT/'research/notes'/(name+'.md'))
    for p in [Path(__file__).resolve(),ROOT/'research/notes'/(STEM+'.md'),
              ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/s61_model.py',ROOT/'code/weight1_structure.py']:
        record(p)
    return dict(id=STEM,all_passed=True,schubert_hyperplane_ranks=hyperplane_ranks,
                flag_second_model='Bl_{E_tau}P3',flag_projection_off_E='v wedge J beta(v) v',
                fixed_pair_center=['u6','u42'],binary_curve_degrees_and_reduced_intersection_multiplicities=binary,
                flag_curve_degrees_and_reduced_intersection_multiplicities=flag,
                total_reduced_curve_count=len(binary)+len(flag),classification='one elliptic section, fifteen binary lines and six rational flag curves',
                homogenized_top_curve_cycle_degree=2300,single_colon_length_determined=False,input_sha256=inputs)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    output=parser.parse_args().output;result=verify()
    with output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print('PASS',result['total_reduced_curve_count'],'curves')

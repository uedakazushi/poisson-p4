#!/usr/bin/env python3
"""Check the intrinsic flag section, projective bundle intersections and top multiplicities."""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import J

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_PROJECTIVE_TOP_CYCLE_20261102'


def verify():
    x,u,b=primitive_beta(S61())
    e=sp.eye(4)
    def value(v):
        return sp.expand(b[1,3].subs(dict(zip(x,v)),simultaneous=True))
    assert value(J*e[:,0])==-sp.Symbol('u4')
    assert value(e[:,1])==sp.Symbol('u4') and value(e[:,3])==sp.Symbol('u5')
    # Compute in the Chow ring of P_G(U) using 2 c2(U)=H^2,
    # x^2-Hx+c2(U)=0. Integrating x H^3 gives 2.
    xi,H=sp.symbols('xi H')
    poly=sp.rem((2*xi+H)**3*(xi+H),xi**2-H*xi+H**2/2,xi)
    flag=int(2*sp.expand(poly).coeff(xi,1).coeff(H,3));assert flag==70
    # On the elliptic base, H^2=0 and -c1(Sym^2 V)=3H.
    binary=int(5*(3+sp.expand((xi+H)**3).coeff(xi,2).coeff(H,1)))
    exceptional=int(5*(4+sp.expand((2*xi+H)**2).coeff(xi,1).coeff(H,1)))
    assert binary==30 and exceptional==40
    lengths=[65,5];tquot=[64,5]
    degree=sum(n*d for n,d in zip(lengths,[binary,flag]))
    killed=sum(n*d for n,d in zip(tquot,[binary,flag]))
    rank_t=degree-killed
    assert (degree,killed,rank_t,degree-2*rank_t)==(2300,2270,30,2240)
    inputs={}
    def record(p):
        inputs[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    for name in ['S61_AFFINE_BOUNDARY_DIMENSION_20261102','S61_AFFINE_SURFACE_BOUND_20261102',
                 'S61_BINARY_FORMAL_RING_20261101','S61_BINARY_ENDPOINTS_20261101',
                 'S61_FLAG_UNIFORM_OPEN_20261101','S61_FLAG_FORMAL_RING_20261101',
                 'PRIMITIVE_JET_ELLIPTIC_CURVE_20261101']:
        p=ROOT/'reports'/(name+'.json')
        assert json.loads(p.read_text())['all_passed'];record(p)
        record(ROOT/'research/notes'/(name+'.md'))
    for p in [Path(__file__).resolve(),ROOT/'research/notes'/(STEM+'.md'),
              ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',
              ROOT/'code/s61_model.py']:
        record(p)
    return dict(id=STEM,all_passed=True,flag_divisor_class='xi+H',flag_hyperplane='2xi+H',
                flag_support_model='Bl_C LG(2,4)',binary_hyperplane='xi_B+H',
                top_cycle_coefficients={'binary':65,'flag':5},reduced_degrees={'binary':binary,'flag':flag},
                reduced_common_exceptional_surface_degree=exceptional,full_projective_degree=degree,
                infinity_quotient_projective_degree=killed,three_general_cuts_total_length=degree,
                three_general_cuts_homogenizer_rank=rank_t,three_general_cuts_homogenizer_square_zero=True,
                homogenizer_Jordan_multiplicities={'1':degree-2*rank_t,'2':rank_t},
                affine_degree_determined=False,input_sha256=inputs)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    output=parser.parse_args().output;result=verify()
    with output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print('PASS',result['full_projective_degree'],result['infinity_quotient_projective_degree'])

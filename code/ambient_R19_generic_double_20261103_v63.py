#!/usr/bin/env python3
"""Exact second-order obstruction gives generic primary length2 along R19.

All calculations use actual CA/CB source families and full canonical150
coordinates, including their genuine finite coordinate transformations.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,raw,canonical,vector
from ambient_CB_global_flag_model_20261103_v63 import raw as cbraw
from ambient_rank_two_flag_intersection_table_20261103_v63 import lie,permutation

ROOT=Path(__file__).resolve().parents[1]

def act(values,M):
    return sum((M[i,j]*lie(values,i,j) for i,j in itertools.product(range(5),repeat=2) if M[i,j]),s.zeros(150,1))

def push(P,g):return (g*P.subs(dict(zip(X,g.inv()*s.Matrix(X))),simultaneous=True)*g.T).applyfunc(s.expand)

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic();h,x,y,z,w=X;kappa,eps=s.symbols('kappa eps')
    geometry_path=ROOT/'reports/AMBIENT_R19_GENERIC_GEOMETRY_20261103_v63/report.json';geo=json.loads(geometry_path.read_text())
    A=s.Matrix([list(map(s.sympify,row)) for row in geo['full_CA_source_differential']]);B=s.Matrix([list(map(s.sympify,row)) for row in geo['full_CB_source_differential']]);W=A.row_join(B)
    point=s.Matrix(list(map(s.sympify,geo['point_canonical150'])))
    uA=s.zeros(45,1);uA[17]=3;uA[28]=2
    uB=s.zeros(45,1);uB[2]=-21;uB[13]=1
    v=A*uA;assert v==B*uB
    MA=s.Matrix(5,5,list(uA[20:,0]));MB=s.Matrix(5,5,list(uB[20:,0]));dm=uB[19]
    assert MA[1,3]==2 and MA*MA==s.zeros(5) and MB==s.zeros(5) and dm==0
    ca=json.loads((ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63/report.json').read_text())
    cb=json.loads((ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json').read_text())
    CA=s.Matrix.hstack(*(s.Matrix(list(map(s.Rational,c))) for c in ca['basis_canonical150']))
    CB0=s.Matrix.hstack(*(s.Matrix(list(map(s.Rational,c))) for i,c in enumerate(cb['enlarged_canonical_basis']) if i<3 or i>=6))
    perm=tuple(geo['CB_source_permutation']);CB=s.Matrix.hstack(*(permutation(CB0[:,i],perm) for i in range(19)))
    dfA=CA*uA[:20,0];dfB=CB*uB[:19,0]+B[:,19]*dm
    binary=[h*h,h*x,x*x];df=sum(uB[i]*binary[i] for i in range(3))
    mixed=dm*permutation(vector(canonical(cbraw(0,df,0))),perm)
    order2A=act(dfA,MA)+(act(act(point,MA),MA)-act(point,MA*MA))/2
    order2B=act(dfB,MB)+(act(act(point,MB),MB)-act(point,MB*MB))/2+mixed
    difference=order2A-order2B;expected=s.zeros(150,1);expected[72]=-18;assert difference==expected
    covector=s.zeros(150,1);covector[63]=-1;covector[72]=kappa
    assert covector.T*W==s.zeros(1,W.cols)
    obstruction=s.expand((covector.T*difference)[0]);assert obstruction==-18*kappa
    # The chosen first-order direction is not tangent to reduced R.
    orbit=s.Matrix.hstack(*(lie(point,i,j) for i,j in itertools.product(range(5),repeat=2)))
    rdirs=[vector(canonical(raw(-4*x*x,0,2*x*z*w)))]+[vector(canonical(raw(0,0,H))) for H in [w*x*x,x*x*z,x*z*z/2,2*h*x*z,z**3/3,h*z*z,3*h*h*z]]
    R=s.Matrix.hstack(*rdirs).row_join(orbit);assert R.rank()==20 and R.row_join(v).rank()==21
    # Direct exact finite curves independently confirm all150 second-order coefficients.
    cubics=list(map(s.sympify,cb['cubic_monomials']));assert cubics[10]==x*z*w
    source_coeff=list(map(s.sympify,geo['CB_source19_coefficients']));f0=sum(source_coeff[i]*binary[i] for i in range(3));H0=sum(source_coeff[i+3]*m for i,m in enumerate(cubics))
    baseCA=raw(-4*x*x,0,h*h*z+kappa*h*z*z/2+z**3/3+2*x*z*w)
    curveCA=canonical(push(baseCA+eps*raw(0,0,3*z*z*w),s.eye(5)+eps*MA))
    gp=s.zeros(5)
    for i,j in enumerate(perm):gp[j,i]=1
    # g=mu*f is multiplied before extracting eps², even though dm=0 here.
    sourceCB=cbraw(f0+eps*df,(eps*dm)*(f0+eps*df),H0+eps*cubics[10])
    curveCB=canonical(push(push(sourceCB,gp),s.eye(5)+eps*MB))
    assert vector(curveCA-curveCB)==eps*eps*difference
    assert vector(curveCA.diff(eps).subs(eps,0))==v
    assert W.rank()==37 and W.row_join(difference).rank()==38
    paths=[Path(__file__),geometry_path,ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',
        ROOT/'code/ambient_CB_global_flag_model_20261103_v63.py',ROOT/'code/ambient_rank_two_flag_intersection_table_20261103_v63.py',
        ROOT/'research/notes/AMBIENT_R19_GENERIC_GEOMETRY_20261103_v63.md']
    body=dict(id='AMBIENT_R19_GENERIC_DOUBLE_20261103_v63',all_passed=True,field='QQ(kappa)',
        first_order_CA_source=list(map(str,uA)),first_order_CB_source=list(map(str,uB)),
        common_first_order_canonical150=list(map(str,v)),CA_group_velocity=MA.tolist(),CB_group_velocity=MB.tolist(),
        CB_mu_velocity=str(dm),CB_mixed_mu_f_term_computed=True,
        actual_second_order_difference_canonical150=list(map(str,difference)),
        tangent_sum_annihilator_canonical150=list(map(str,covector)),obstruction_value=str(obstruction),
        exact_finite_curve_difference='eps²*(-18*z²*d_x wedge d_y)',
        no_omitted_higher_eps_terms=True,tangent_sum_rank=37,obstruction_augmented_rank=38,
        reduced_R_cone_tangent_rank=20,extra_direction_augmented_R_rank=21,
        generic_primary_algebra='K(R)[eps]/(eps²)',generic_primary_length=2,
        scope='Actual genericprimarystructure alongR19, not a globalprimarydecomposition; kappa=0 is not decided bythisobstruction.',
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(body,indent=2,default=str)+'\n');print('R19 EXACT GENERIC DOUBLE PASS',body['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)

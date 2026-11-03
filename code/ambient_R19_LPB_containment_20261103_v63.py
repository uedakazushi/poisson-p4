#!/usr/bin/env python3
"""Exact full8-parameter R source embedding into the existing LPB39 family."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,canonical,vector

ROOT=Path(__file__).resolve().parents[1]

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic();h,x,y,z,w=X
    t,q0,q1,q2,q3,q5,q6,q7=s.symbols('t q0 q1 q2 q3 q5 q6 q7')
    field=s.Matrix([0,4*t*x*x,0,-q0*x*x-2*t*x*z,2*t*x*w+3*q7*h*h+2*q3*h*x+2*q6*h*z+q1*x*x+q2*x*z+q5*z*z])
    ey=s.Matrix([0,0,1,0,0]);P=ey*field.T-field*ey.T
    perm=[0,1,4,2,3];g=s.zeros(5)
    for i,j in enumerate(perm):g[j,i]=1
    transformed=(g*P.subs(dict(zip(X,g.inv()*s.Matrix(X))),simultaneous=True)*g.T).applyfunc(s.expand)
    quotient=s.Matrix([0,4*t*x*x,-q0*x*x-2*t*x*y,2*t*x*z+3*q7*h*h+2*q3*h*x+2*q6*h*y+q1*x*x+q2*x*y+q5*y*y,0])
    e4=s.Matrix([0,0,0,0,1]);LPB=e4*quotient.T-quotient*e4.T
    assert transformed==LPB and all(s.diff(v,w)==0 for v in quotient)
    assert canonical(transformed-LPB)==s.zeros(5)
    old_path=ROOT/'legacy/v14/data/six_rank2_components.json';old=json.loads(old_path.read_text())
    records=old if isinstance(old,list) else old.get('components',old.get('cases',[]))
    if not isinstance(records,list):records=[]
    lpb=next(v for v in records if v.get('name')=='LPB')
    assert lpb['projective_component_dimension']==39 and lpb['affine_component_dimension']==40
    assert lpb['ambient_affine_tangent_dimension']==lpb['family_affine_differential_rank']==40
    paths=[Path(__file__),ROOT/'legacy/v14/code/rank2_lower_bounds/lpb_family.py',old_path,
        ROOT/'legacy/v14/sections/rank2_global_eight.tex',ROOT/'research/notes/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63.md']
    body=dict(id='AMBIENT_R19_LPB_CONTAINMENT_20261103_v63',all_passed=True,field='QQ(all8parameters)',
        coordinate_permutation=perm,LPB_center='e4',LPB_quotient_quadratic_field=list(map(str,quotient)),
        full8parameter_raw_and_canonical_identity=True,existing_LPB_component_projective_dimension=39,
        existing_component_dimension_certificate=lpb,
        entire_R_containment='R19 subset LPB39, by exactsourceembedding and closedPGL-invariance',
        full_Poisson_completed_local_consequence='NON Cohen-Macaulay at EVERYclosedpointofR, because actualminimalcomponentsCA28/CB28 andLPB39 have unequallocaldimensions',
        component_count_change=0,
        validation_scope='OldLPBcomponentcertificate reused and hashed; itslargeJacobian calculation was notrerun. The newallparameterLPBsourceidentity isexactQQ.',
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(body,indent=2)+'\n');print('ENTIRE R19 IN LPB39 PASS',body['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)

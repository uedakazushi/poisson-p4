#!/usr/bin/env python3
"""The eleven already proved non-MPS rank-four families in the common X order.

This is a coordinate catalogue, not a new component or exhaustion theorem.
Exact symbolic Jacobi and fresh image differential certificates are executable.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s

from s61_model import X,J,PAIRS,coefficients,jacobi_projective,vertical_rref
from mps_s61_dominance_20261103_v62 import rank_minor,PRIME
ROOT=Path(__file__).resolve().parents[1]

def four(values):
    P=s.zeros(5)
    for (i,j),value in zip(itertools.combinations(range(1,5),2),values):P[i,j]=value;P[j,i]=-value
    return P

def families():
    h,x,y,z,w=X;result=[]
    def add(name,P,parameters,dimension,path,anchor,field_relation=None):
        result.append(dict(component_id=name,P=P,parameters=parameters,projective_dimension=dimension,
                           source=path,source_anchor=anchor,field_relation=field_relation))
    P=four([-x*x+2*h*y,2*x*y+3*h*z,5*h*w+y*y+2*x*z,h*w+5*y*y+4*x*z,2*y*z+4*x*w,3*z*z-8*y*w])
    add('R4_DOUBLE_SMOOTH_QUADRIC',P,[],22,'legacy/v03/p4_poisson_completeness_v3.tex','sec:new, prop:tangent, thm:quadric')
    P=h*h*J+four([3*x*w,-3*x*x,3*x*z,3*x*y-z*w,-3*y*z+w*w,-x*w+z*z])
    add('R4_TRIPLE_HYPERPLANE_QUADRIC',P,[],22,'legacy/v05/triple_hyperplane.tex','sec:new')
    a=s.symbols('a0:3');b=s.symbols('b0:16');f=a[0]*h*h+a[1]*h*x+a[2]*x*x
    mons=[s.prod((h,x,z,w)[i] for i in ids) for ids in itertools.combinations_with_replacement(range(4),3)
          if any(i>=2 for i in ids)];assert len(mons)==16
    potential=sum(c*m for c,m in zip(b,mons))
    for name,g,extra in [('R4_CA',h*h,0),('R4_CB',f,-s.diff(f,x)/2)]:
        P=four([f,0,0,extra*z-s.diff(potential,w),extra*w+s.diff(potential,z),g])
        add(name,P,[*a,*b],28,'legacy/v06/ambient_families.tex','eq:PA, eq:PB, thm:AB')
    P=h*h*J+four([-x*x+x*z-x*w,-x*x,0,3*x*y-4*x*z+z*z-3*z*w,
                   x*x+x*y-2*x*z+2*x*w+z*z/2-2*z*w+w*w,x*x+x*w])
    add('R4_H_QUARTIC_ESSENTIAL5',P,[],23,'legacy/v06/nilpotent22.tex','subsec:quarticcomponent')
    a,b=s.symbols('a b');P=h*h*J+four([-a*b*x*z,0,-b*b*x*x,-a*a*z*z,b*(2*b-a)*x*y-a*(b-2*a)*z*w,a*b*x*z])
    add('R4_ZERO_LOG_5H',P,[a,b],23,'legacy/v06/zero_log_component.tex','eq:zeroQ, thm:zero_log')
    P=h*h*J+four([3*x*x-6*x*z,0,-2*x*x,3*x*z-18*z*z,6*x*x-2*x*y-3*x*w+30*z*w,6*x*z])
    add('R4_MINIMAL_M',P,[],23,'legacy/v07/minimal_thickness.tex','thm:minimal-full')
    P=h*h*J+four([3*x*z/4,0,z*z+3*x*x/4,-z*z/4,x*x+2*x*z+3*z*z-3*x*y/4+z*w/4,x*z/4])
    add('R4_A0',P,[],23,'legacy/v07/nil22_v7.tex','A0 full P4 subsection')
    r=s.Symbol('r');u=s.Rational(1,6);v=-s.Rational(1,36);c=r/18;d=s.Rational(1,2)+c;a=-3*c;b=-c/12
    P=h*h*J+four([(1-d)*x*z-v*x*x-u*z*z,0,a*z*z-4*u*x*z+(c-3*d+2)*x*x/2,
        b*x*x-4*v*x*z-(c+d)*z*z/2,x*x+2*x*z+3*z*z+(2*d-c-1)*x*y+c*z*w+2*u*y*z+2*v*x*w,d*x*z+v*x*x+u*z*z])
    add('R4_DPM',P,[],23,'legacy/v07/dplus_fullp4.tex','thm:dplus-full; common normal form in nil22_v7.tex',r*r-3)
    c,beta=s.symbols('c beta');P0=h*h*J+four([x*x+8*x*z-z*z,0,5*x*x-4*x*z-z*z,
        x*x+4*x*z-5*z*z,x*x+2*(z-x)*(y+w),-x*x+8*x*z+z*z])
    H=h+4*x-4*z;A=3*h*x-5*h*z+7*x*x-12*x*z+5*z*z;B=5*h*x-3*h*z+5*x*x-12*x*z+7*z*z
    U=s.Matrix([2*h*H,-A,0,B,0]);V=s.Matrix([0,0,-1,0,1])
    P=P0+(c+1)*(U*V.T-V*U.T)/5
    P[2,4]+=(c+1)*(2*h+18*(x-z))*(y+w)/5+beta*h*x;P[4,2]=-P[2,4]
    add('R4_CPLUS_EXT',P,[c,beta],24,'legacy/v08/sections/cplus_component.tex','eq:cplus-universal-family, thm:cplus-family')
    beta,gamma=s.symbols('beta gamma');P=four([h*h+x*x,0,0,-x*z,-x*w-z*z-beta*x*z-gamma*h*z,h*h+x*x-2*z*z])
    add('R4_CMINUS_MAG',P,[beta,gamma],23,'legacy/v08/sections/cminus_component.tex','eq:cminusmag-family, thm:cminusmag')
    return result

def verify(record):
    start=time.monotonic();P=record['P'];relation=record['field_relation'];r=s.Symbol('r')
    residuals=jacobi_projective(P)
    if relation is not None:residuals=[s.rem(s.together(v).as_numer_denom()[0],relation,r) for v in residuals]
    assert all(v==0 for v in residuals)
    # Rational/algebraic source points. The same reduction mod32003 is an
    # honest residue-field specialization of QQ or QQ[r]/(r^2-3).
    point={v:s.Integer(i+2) for i,v in enumerate(record['parameters'])}
    embedding=None
    if relation is not None:embedding=int(s.sqrt_mod(3,PRIME));point[r]=s.Integer(embedding)
    sample=P.subs(point);columns=[];labels=[]
    for i in range(5):
        for j in range(5):
            g=s.zeros(5);g[i,j]=1
            derivative=g*sample+sample*g.T-sum(((g*s.Matrix(X))[k]*sample.diff(X[k]) for k in range(5)),s.zeros(5))
            columns.append(s.Matrix(coefficients(derivative)));labels.append(f'gl_{i}_{j}')
    columns.append(s.Matrix(coefficients(sample)));labels.append('overall_scalar')
    for variable in record['parameters']:
        columns.append(s.Matrix(coefficients(P.diff(variable).subs(point))));labels.append(str(variable))
    matrix=s.Matrix.hstack(*columns);vertical,pivots=vertical_rref();free=[i for i in range(150) if i not in pivots]
    quotient=matrix[free,:]-vertical[:,free].T*matrix[list(pivots),:];minor=rank_minor(quotient)
    assert minor['rank']==record['projective_dimension']+1,(record['component_id'],minor['rank'])
    return dict(component_id=record['component_id'],generic_projective_jacobi=True,
        coefficient_field='QQ[r]/(r^2-3)' if relation is not None else 'QQ',field_relation=str(relation) if relation is not None else None,
        source_parameters=list(map(str,record['parameters'])),source_sample={str(k):str(v) for k,v in point.items()},
        residue_field_embedding_r_mod32003=embedding,projective_image_dimension=record['projective_dimension'],
        s61_dimension=record['projective_dimension']-9,lower_bound_minor=minor,
        differential_column_order=labels,differential_matrix_at_residue_lift_QQ=[[str(v) for v in row] for row in quotient.tolist()],
        component_upper_bound_and_identification_dependency=dict(path=record['source'],anchor=record['source_anchor']),
        homogeneous_bivector={f'{i},{j}':str(s.expand(P[i,j])) for i,j in PAIRS},seconds=time.monotonic()-start)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False);records=[]
    for record in families():
        body=verify(record);path=args.output/(record['component_id']+'.json');path.write_text(json.dumps(body,indent=2)+'\n')
        records.append(dict(component_id=record['component_id'],path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        print(record['component_id'],'GENERIC JACOBI + IMAGE LOWER PASS',body['seconds'],flush=True)
    sources={Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/mps_s61_dominance_20261103_v62.py',*[ROOT/r['source'] for r in families()],ROOT/'legacy/v07/nil22_v7.tex'}
    (args.output/'report.json').write_text(json.dumps(dict(id='ADDITIONAL_RANK4_FAMILIES_20261103_v62',all_passed=True,records=records,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources)},
        new_components_claimed=0,global_exhaustion_claimed=False),indent=2)+'\n')

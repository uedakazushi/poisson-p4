#!/usr/bin/env python3
"""A global27dim common CA/CB family containing M and the distinguished Q."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from additional_rank4_families_20261103_v62 import four,ROOT
from s61_model import X,PAIRS,S61,coefficients,pushforward,jacobi_projective,vertical_rref
from mps_s61_dominance_20261103_v62 import rank_minor

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic();h,x,y,z,w=X;p=s.Symbol('p');a=s.symbols('a0:16')
    monomials=[s.prod((h,x,z,w)[i] for i in ids) for ids in itertools.combinations_with_replacement(range(4),3)
        if any(i>=2 for i in ids)];potential=sum(c*m for c,m in zip(a,monomials))
    f=h*(h+p*x);PA=four([f,0,0,-s.diff(potential,w),s.diff(potential,z),h*h])
    assert all(v==0 for v in jacobi_projective(PA))
    M=s.eye(5);M[0,1]=p;transformed=pushforward(PA,M);oldh=h-p*x
    fB=oldh*oldh;potentialB=s.expand(potential.subs(h,oldh));radial=-s.diff(fB,x)/2
    PB=four([fB,0,0,radial*z-s.diff(potentialB,w),radial*w+s.diff(potentialB,z),fB])
    vector=s.Matrix([0,0,p*oldh,0,0]);euler=s.Matrix(X);vertical=euler*vector.T-vector*euler.T
    assert (transformed-vertical-PB).applyfunc(s.expand)==s.zeros(5)
    assert all(v==0 for v in jacobi_projective(PB))
    # Fresh actual differential, including every projective coordinate direction.
    parameters=[p,*a];sample={p:s.Integer(2)}|{v:s.Integer(i+2) for i,v in enumerate(a)};P=PA.subs(sample)
    columns=[];labels=[]
    for i in range(5):
        for j in range(5):
            g=s.zeros(5);g[i,j]=1
            derivative=g*P+P*g.T-sum(((g*s.Matrix(X))[k]*P.diff(X[k]) for k in range(5)),s.zeros(5))
            columns.append(s.Matrix(coefficients(derivative)));labels.append(f'gl_{i}_{j}')
    columns.append(s.Matrix(coefficients(P)));labels.append('overall_scalar')
    for variable in parameters:columns.append(s.Matrix(coefficients(PA.diff(variable).subs(sample))));labels.append(str(variable))
    matrix=s.Matrix.hstack(*columns);rref,pivots=vertical_rref();free=[i for i in range(150) if i not in pivots]
    quotient=matrix[free,:]-rref[:,free].T*matrix[list(pivots),:];minor=rank_minor(quotient);assert minor['rank']==28
    # Exact coordinate transport of the v24 distinguished endpoint.
    t,u,v,ww,zz=X;Q=s.zeros(5)
    for (i,j),value in {(1,4):t*t,(2,3):-s.Rational(2,9)*t*t,(2,4):s.Rational(4,3)*t*u,
        (3,4):s.Rational(2,3)*t*v+u*u}.items():Q[i,j]=value;Q[j,i]=-value
    transport=s.zeros(5);transport[0,0]=transport[1,1]=transport[2,4]=transport[3,2]=1;transport[4,3]=-s.Rational(9,2)
    target_potential=s.Rational(4,3)*h*x*w+s.Rational(3,2)*h*z*z+s.Rational(9,2)*x*x*z
    point={p:0}|{c:s.Poly(target_potential,*X).coeff_monomial(mon) for c,mon in zip(a,monomials)}
    assert (pushforward(Q,transport)-PA.subs(point)).applyfunc(s.expand)==s.zeros(5)
    # At p=0 the whole sixteen-parameter magnetic family is retained.
    assert PA.subs(p,0)==four([h*h,0,0,-s.diff(potential,w),s.diff(potential,z),h*h])
    payload=dict(id='AMBIENT_INTERSECTION_DIMENSION27_20261103_v62',all_passed=True,
        generic_projective_jacobi=True,universal_CB_coordinate_identity=True,
        parameter_order=list(map(str,parameters)),cubic_monomials=list(map(str,monomials)),
        homogeneous_family={f'{i},{j}':str(PA[i,j]) for i,j in PAIRS},
        perspective_matrix=[[str(v) for v in row] for row in M.tolist()],
        CB_binary_quadratic=str(fB),CB_potential=str(potentialB),Euler_vertical_vector=list(map(str,vector)),
        differential_sample={str(k):str(v) for k,v in sample.items()},differential_column_order=labels,
        quotient_raw_indices=free,differential_matrix_QQ=[[str(v) for v in row] for row in quotient.tolist()],
        affine_nonzero_minor=minor,projective_image_dimension_lower_bound=27,
        projective_image_dimension_upper_bound=27,
        upper_dependency='CA and CB are distinct irreducible28dim components; their intersection has dimension at most27.',
        distinguished_Q_transport=[[str(v) for v in row] for row in transport.tolist()],
        distinguished_Q_parameter_values={str(k):str(v) for k,v in point.items()},Q_exactly_in_family=True,
        entire_magnetic_M_contained=True,
        conclusions=['The irreducible PGL-saturated closure N is a27dim component of CA intersection CB.',
            'The local dimension at the distinguished Q of CA intersection CB is exactly27.',
            'The previously known irreducible26dim magnetic M is a proper closed subset of N, so is not an intersection component.'],
        limitations=['No uniqueness of the27dim intersection component through Q.','No completed local ring or scheme-theoretic intersection ideal.','No global Poisson component count upper.'],
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in [Path(__file__),
            ROOT/'code/additional_rank4_families_20261103_v62.py',ROOT/'code/s61_model.py',ROOT/'code/mps_s61_dominance_20261103_v62.py',
            ROOT/'legacy/v06/ambient_families.tex',ROOT/'legacy/v24/sections/incidence.tex',ROOT/'legacy/v24/sections/setup.tex']},
        seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(payload,indent=2)+'\n');print('L1 EXACT DIMENSION27 PASS',minor,payload['seconds'],flush=True)

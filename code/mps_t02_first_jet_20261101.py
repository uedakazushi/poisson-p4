#!/usr/bin/env python3
"""Exact T02 single-edge family and a rank15 primitive first-jet witness.

The upper bound15 uses the existing MPS component identification/dimension.
No equation of its first-jet image hypersurface or full fiber is constructed.
"""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import sympy as s

from s61_model import X,jacobi_projective
from elliptic_weight1_dominance import projective_affine,symplectic_frame,primitive_vector
from weight1_structure import J,infinitesimal_action,sp_basis,minor_witness,constraints

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    parameters=s.symbols('a d e f b');a,d,e,f,b=parameters
    bivector=s.zeros(5);bivector[0,1]=b*X[2]**2;bivector[1,0]=-bivector[0,1]
    toric_values=[a,2*d+a,2*e+a,d,e,f]
    for (i,j),v in zip(combinations(range(1,5),2),toric_values):
        bivector[i,j]=v*X[i]*X[j];bivector[j,i]=-bivector[i,j]
    assert all(v==0 for v in jacobi_projective(bivector))
    y=s.symbols('y0:4');affine=projective_affine(bivector,X,y)
    parameter_values=[2,3,5,7,-2];point_values=[2,1,1,1]
    substitution=dict(zip(parameters,parameter_values))|dict(zip(y,point_values))
    constant=affine.subs(substitution);frame=symplectic_frame(constant);inverse=frame.inv()
    derivatives=[affine.diff(v).subs(substitution) for v in y]
    linear=[frame*sum((inverse[h,j]*derivatives[h] for h in range(4)),s.zeros(4))*frame.T
            for j in range(4)]
    value=primitive_vector(linear);columns=[]
    for variable in [*parameters,*y]:
        change_constant=affine.diff(variable).subs(substitution)
        change_derivatives=[affine.diff(v).diff(variable).subs(substitution) for v in y]
        h=frame*change_constant*frame.T*J/2
        assert h*J+J*h.T+frame*change_constant*frame.T==s.zeros(4)
        change_linear=[h*linear[j]+linear[j]*h.T
            -sum((h[k,j]*linear[k] for k in range(4)),s.zeros(4))
            +frame*sum((inverse[k,j]*change_derivatives[k] for k in range(4)),s.zeros(4))*frame.T
            for j in range(4)]
        columns.append(primitive_vector(change_linear))
    differential=s.Matrix.hstack(*columns,*[infinitesimal_action(value,g) for g in sp_basis()])
    assert s.Matrix.vstack(*constraints())*differential==s.zeros(8,19)
    witness=minor_witness(differential)
    assert witness['rank']==15 and s.Rational(witness['determinant'])==s.Rational(1843456,378125)

    toric=s.zeros(4)
    for (i,j),v in zip(combinations(range(4),2),toric_values):
        toric[i,j]=v;toric[j,i]=-v
    period=-toric.subs(substitution).inv();biresidue=s.zeros(5);biresidue[1:5,1:5]=period
    for i in range(1,5):
        biresidue[i,0]=-sum(biresidue[i,j] for j in range(1,5));biresidue[0,i]=-biresidue[i,0]
    edges=[]
    for i,j in combinations(range(5),2):
        if not biresidue[i,j]:continue
        opposite=[k for k in range(5) if k not in (i,j)]
        orders=[s.cancel((biresidue[j,k]+biresidue[k,i])/biresidue[i,j]) for k in opposite]
        if all(v.is_Integer and v>=0 for v in orders):
            edges.append({'edge':[i,j],'opposite_vertices':opposite,'orders':list(map(int,orders))})
    assert edges==[{'edge':[0,1],'opposite_vertices':[2,3,4],'orders':[2,0,0]}]
    report={
        'id':'MPS_T02_FIRST_JET_20261101','all_passed':True,'component_id':'T02',
        'family_field':'QQ(a,d,e,f,b)','parameter_order':list(map(str,parameters)),
        'homogeneous_bivector':[[str(v) for v in row] for row in bivector.tolist()],
        'universal_projective_Jacobi_checked':True,
        'sample_parameters':parameter_values,'sample_affine_point':point_values,
        'sample_biresidue':[[str(v) for v in row] for row in biresidue.tolist()],
        'sample_smoothable_edges':edges,'sample_primitive_beta':list(map(str,value)),
        'differential_column_order':list(map(str,[*parameters,*y]))+[f'sp4_{i}' for i in range(10)],
        'horizontal_frame_differential':[[str(v) for v in row] for row in differential.tolist()],
        'nonzero_rank15_minor':witness,
        'image_dimension_lower_bound':15,'image_dimension_upper_bound':15,
        'upper_bound_dependency':'Existing MPS T02 projective component dimension24, hence S61 dimension15; component identification uses the single smoothable edge with orders(2,0,0).',
        'scope':'The reduced closed primitive first-jet image is an irreducible hypersurface in W16; its defining polynomial and generic fiber degree are not determined.',
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/elliptic_weight1_dominance.py',
                      ROOT/'code/weight1_structure.py',ROOT/'reports/MPS40_CROSSWALK_20261031.json')}
    }
    with args.output.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps({'all_passed':True,'component':'T02','image_dimension':15,
                      'minor':witness['determinant']}))


if __name__=='__main__':main()

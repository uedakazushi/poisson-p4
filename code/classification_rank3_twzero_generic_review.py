#!/usr/bin/env python3
"""One full-ambient test of the entire irreducible T(w)=0 source."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_twzero_source import projective_chart
from classification_rank3_mixed_component_review import (
    COLUMNS, jacobi, sparse, vector, minor_certificate)

ROOT=Path(__file__).resolve().parents[1]


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic()
    source,xx,params=projective_chart()
    values=dict(zip(params,(1,1,1,2,2,2,2,1,1,1,1)))
    P=source.subs(values);ps=sparse(P,xx)
    assert all(not q for q in jacobi(ps,raw=True).values())
    directions=[source.diff(p).subs(values) for p in params]
    scales=[]
    for i,R in enumerate(directions):
        scale=s.ilcm(1,*[s.denom(c) for ij in itertools.combinations(range(5),2)
                          for c in s.Poly(R[ij],*xx).coeffs()])
        directions[i]=R*scale;scales.append(str(scale))
    assert all(c.q==1 for A in [P,*directions]
               for ij in itertools.combinations(range(5),2)
               for c in s.Poly(A[ij],*xx).coeffs())
    columns=[]
    for ij,e in COLUMNS:
        variation={key:dict(value)for key,value in ps.items()}
        variation.setdefault(ij,{})[e]=variation.get(ij,{}).get(e,0)+1
        columns.append(jacobi(variation))
    D=s.Matrix.hstack(*columns)
    X=s.Matrix(xx);vcols=[];ocols=[]
    for a,b in itertools.product(range(5),repeat=2):
        field=s.zeros(5,1);field[a]=xx[b]
        vcols.append(vector(sparse(X*field.T-field*X.T,xx)))
        E=s.zeros(5);E[a,b]=1
        ocols.append(vector(sparse(E*P+P*E.T-xx[b]*P.diff(xx[a]),xx)))
    V=s.Matrix.hstack(*vcols);M=V.row_join(s.Matrix.hstack(*ocols)).row_join(vector(ps))
    F=M.row_join(s.Matrix.hstack(*(vector(sparse(R,xx))for R in directions)))
    assert D*F==s.zeros(350,F.cols)
    mats={'full350_D':D,'vertical':V,'orbit_vertical_scale':M,'actual_full_source':F}
    ranks={name:len(DomainMatrix.from_Matrix(A).convert_to(s.QQ).rref()[1])for name,A in mats.items()}
    certs={name:minor_certificate(A,ranks[name])for name,A in mats.items()}
    lower=ranks['actual_full_source']-25;upper=125-ranks['full350_D']
    report=dict(task_id='E_quadric3',all_passed=True,
        source_parameters=list(map(str,params)),point_substitution={str(k):str(v)for k,v in values.items()},
        point_brackets={str(ij):str(P[ij])for ij in itertools.combinations(range(5),2)if P[ij]!=0},
        all11_actual_parameter_directions=[{str(ij):str(R[ij])for ij in itertools.combinations(range(5),2)if R[ij]!=0}for R in directions],
        parameter_direction_integer_scales=scales,
        raw_Jacobi_all_ten_zero=True,full350_times_all_actual_directions_zero=True,
        ranks_QQ=ranks,certificates=certs,projective_family_dimension_lower=lower,
        projective_tangent_dimension_upper=upper,
        smooth_full_component_proved=(lower==upper),
        global_exhaustion_not_claimed=True,
        input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()for name in
            ['code/classification_rank3_twzero_generic_review.py',
             'code/classification_rank3_twzero_source.py',
             'code/classification_rank3_mixed_component_review.py']},
        seconds=time.monotonic()-begun)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(dict(ranks_QQ=ranks,lower=lower,upper=upper,
                         component_proved=lower==upper,seconds=report['seconds'])))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output)

#!/usr/bin/env python3
"""Exact characteristic-zero image-dimension certificates for MPS T01--T38.

Nonzero minors are recorded modulo32003 from rational matrices with unit
denominators. They prove rank lower bounds over QQ, never a rank upper bound.
The matching upper dimension and component label use the existing MPS local
deformation theorem and the checked smoothing-diagram crosswalk.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s

from mps_s61_catalogue_20261103_v62 import raw_family,DIAGRAMS,ROOT
from s61_model import X,coefficients,vertical_rref

PRIME=32003


def rank_minor(matrix):
    original=[[int(v.p)*pow(int(v.q),-1,PRIME)%PRIME for v in row] for row in matrix.tolist()]
    data=[row[:] for row in original];labels=list(range(matrix.rows));pivots=[];rows=[];r=0
    for c in range(matrix.cols):
        pivot=next((i for i in range(r,matrix.rows) if data[i][c]),None)
        if pivot is None:continue
        data[r],data[pivot]=data[pivot],data[r];labels[r],labels[pivot]=labels[pivot],labels[r]
        rows.append(labels[r]);pivots.append(c)
        inverse=pow(data[r][c],-1,PRIME)
        data[r]=[v*inverse%PRIME for v in data[r]]
        for i in range(r+1,matrix.rows):
            if data[i][c]:
                factor=data[i][c];data[i]=[(a-factor*b)%PRIME for a,b in zip(data[i],data[r])]
        r+=1
    sub=[[original[i][j] for j in pivots] for i in rows];det=1
    for c in range(r):
        p=next(i for i in range(c,r) if sub[i][c])
        if p!=c:sub[p],sub[c]=sub[c],sub[p];det=-det
        det=det*sub[c][c]%PRIME;inv=pow(sub[c][c],-1,PRIME)
        for i in range(c+1,r):
            if sub[i][c]:
                factor=sub[i][c]*inv%PRIME
                sub[i]=[(a-factor*b)%PRIME for a,b in zip(sub[i],sub[c])]
    assert det%PRIME
    return dict(rank=r,rows=rows,columns=pivots,determinant_mod_32003=det%PRIME)


def certificate(record):
    started=time.monotonic()
    periods,smoothings,P,pf,A=raw_family(record)
    Bexample=s.Matrix(record['biresidue_example'])[1:5,1:5]
    v=s.Matrix([Bexample[i,j] for i,j in itertools.combinations(range(4),2)])
    base=s.Matrix([[s.Rational(row[i]) for row in record['period_basis']] for i in range(6)])
    values=base.gauss_jordan_solve(v)[0]
    point=dict(zip(periods,values))|dict.fromkeys(smoothings,0)
    assert pf.subs(point)!=0
    toric=P.subs(point)
    raw=[];labels=[]
    for i in range(5):
        for j in range(5):
            g=s.zeros(5);g[i,j]=1
            derivative=g*toric+toric*g.T-sum(((g*s.Matrix(X))[k]*toric.diff(X[k])
                                               for k in range(5)),s.zeros(5))
            raw.append(s.Matrix(coefficients(derivative)));labels.append(f'pgl_{i}_{j}')
    for variable in [*periods,*smoothings]:
        raw.append(s.Matrix(coefficients(P.diff(variable).subs(point))));labels.append(str(variable))
    matrix=s.Matrix.hstack(*raw)
    vertical,pivots=vertical_rref();free=[i for i in range(150) if i not in pivots]
    quotient=matrix[free,:]-vertical[:,free].T*matrix[list(pivots),:]
    witness=rank_minor(quotient)
    expected=20+len(periods)+len(smoothings)
    assert witness['rank']==expected
    assert expected-1==record['projective_component_dimension']
    return dict(component_id=record['id'],all_passed=True,
                period_sample=list(map(str,values)),smoothing_sample=[0]*len(smoothings),
                quotient_raw_indices=free,column_order=labels,
                differential_matrix_QQ=[[str(v) for v in row] for row in quotient.tolist()],
                nonzero_minor=witness,
                affine_component_image_dimension_lower_bound=expected,
                projective_component_image_dimension=expected-1,
                s61_component_dimension=expected-10,
                upper_bound_dependency='Existing MPS component dimension and local smoothing-diagram identification; this script proves the matching lower bound only.',
                seconds=time.monotonic()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    records=[]
    for record in json.loads(DIAGRAMS.read_text())['diagrams']:
        if int(record['id'][1:])>38:continue
        report=certificate(record);(args.output/f'{record["id"]}.json').write_text(json.dumps(report,indent=2)+'\n')
        records.append({k:report[k] for k in ['component_id','s61_component_dimension','seconds']})
        print(record['id'],'RANK',report['nonzero_minor']['rank'],'S61',report['s61_component_dimension'],flush=True)
    inputs=[Path(__file__),DIAGRAMS,ROOT/'code/mps_s61_catalogue_20261103_v62.py',ROOT/'code/s61_model.py',
            ROOT/'reports/MPS40_CROSSWALK_20261031.json']
    (args.output/'report.json').write_text(json.dumps(dict(id='MPS_S61_DOMINANCE_20261103_v62',all_passed=True,
        records=records,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        new_components_claimed=0,global_upper_bound_proved=False),indent=2)+'\n')

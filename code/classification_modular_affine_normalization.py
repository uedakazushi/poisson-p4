#!/usr/bin/env python3
"""Verify a supplied all-source affine normalization without dropping charts."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

import sympy as s


def main(path):
    resource.setrlimit(resource.RLIMIT_CPU,(30,31))
    resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    started=time.monotonic();data=json.loads(path.read_text())
    source_path=Path(data['source_report'])
    source=next(v for v in json.loads(source_path.read_text())['branches'] if v['name']==data['source_branch'])
    xx=s.symbols('h x1 x2 x3 x4');nn=s.symbols('n0:15')
    env={str(v):v for v in xx+nn}
    def matrix(values):
        P=s.zeros(5)
        for key,value in values.items():
            i,j=map(int,key.split(','));P[i,j]=s.sympify(value,locals=env)
        return P-P.T
    P=matrix(source['universal_brackets']);C=matrix(data['canonical_brackets'])
    coordinates=s.Matrix([s.sympify(v,locals=env) for v in data['coordinates']])
    A=coordinates.jacobian(xx)
    assert all(s.expand(v)==0 for v in coordinates-A*s.Matrix(xx))
    determinant=s.factor(A.det());assert determinant==1
    assert all(s.expand(v)==0 for v in A*P*A.T-C.subs(dict(zip(xx,coordinates)),simultaneous=True))
    for i,j,k in itertools.combinations(range(5),3):
        assert s.expand(sum(C[i,l]*s.diff(C[j,k],xx[l])+C[j,l]*s.diff(C[k,i],xx[l])+C[k,l]*s.diff(C[i,j],xx[l]) for l in range(5)))==0
    weights=data['weight']
    for (i,j) in itertools.combinations(range(5),2):
        for mon,coef in s.Poly(C[i,j],*xx).terms():
            if coef:assert weights[i]+weights[j]-sum(a*b for a,b in zip(weights,mon))==data['P_weight']
    pf=s.factor(C[1,2]*C[3,4]-C[1,3]*C[2,4]+C[1,4]*C[2,3])
    old_at={env[k]:s.Rational(v) for k,v in source['full_source_specification']['witness'].items()}
    canonical_at={env[k]:s.Rational(v) for k,v in data['canonical_parameter_values'].items()}
    assert C.subs(old_at)==C.subs(canonical_at)
    result=dict(data,universal_all_source_identity=True,affine_determinant=str(determinant),
                universal_raw_jacobi_zero=True,universal_affine_pfaffian=str(pf),
                homogeneous_weight_identity=True,same_frozen_witness_after_affine_change=True,
                source_report_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                producer_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-started,4))
    Path(data['output_report']).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(branch=data['branch'],all_passed=True,pfaffian=str(pf),seconds=result['seconds'])))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('input',type=Path);main(ap.parse_args().input)

#!/usr/bin/env python3
"""Independent homogeneous46-variable expansion of the saved t^4 identity.

Uses no reduction, inverse, scratch cache, or native action evaluator.
It multiplies each saved original-equation cofactor by its actual full quadric
and compares every homogeneous quartic coefficient with t^4-normal_form.
"""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time

import numpy as np

from s61_affine_seventeen_base_lift_20261103_v62 import build
from s61_affine_normal_projection_20261102 import ROOT,P


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--identity',type=Path,default=ROOT/'reports/S61_AFFINE_T_FOURTH_NORMAL_FORM_20261103_v62/identity.json.gz')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();assert not args.output.exists();started=time.monotonic()
    high,cross,base,zpairs,_,_,ypairs,_,inputs,_,_=build()
    pairs=list(itertools.combinations_with_replacement(range(46),2));pi={m:i for i,m in enumerate(pairs)}
    quartics=list(itertools.combinations_with_replacement(range(46),4));qi={m:i for i,m in enumerate(quartics)}
    full=np.zeros((134,len(pairs)),dtype=np.int64)
    for j,mon in enumerate(zpairs):full[:,pi[mon]]=high[:,j]
    for i in range(29):
        for b in range(17):full[:,pi[i,29+b]]=cross[:,i,b]
    for j,(b,c) in enumerate(ypairs):full[:,pi[29+b,29+c]]=base[:,j]
    with gzip.open(args.identity,'rt')as stream:identity=json.load(stream)
    cofactor=np.zeros_like(full)
    for r,i,j,c in identity['cofactor_quadratic']:
        k=pi[i,j];assert cofactor[r,k]==0;cofactor[r,k]=c
    for b,r,k,c in identity['cofactor_base_linear_unknown_linear']:
        j=pi[k,29+b];assert cofactor[r,j]==0;cofactor[r,j]=c
    for b,c,r,value in identity['cofactor_base_quadratic_constant']:
        j=pi[29+b,29+c];assert cofactor[r,j]==0;cofactor[r,j]=value
    expected=np.zeros(len(quartics),dtype=np.int64);expected[qi[(0,)*4]]=1
    normal_indices=[]
    for b,i,j,k,c in identity['normal_base_linear_unknown_cubic']:
        normal_indices.append(qi[tuple(sorted((29+b,i,j,k)))]);expected[normal_indices[-1]]-=c
    for b,c,i,j,value in identity['normal_base_quadratic_unknown_quadratic']:
        normal_indices.append(qi[tuple(sorted((29+b,29+c,i,j)))]);expected[normal_indices[-1]]-=value
    for b,c,d,i,value in identity['normal_base_cubic_unknown_linear']:
        normal_indices.append(qi[tuple(sorted((29+b,29+c,29+d,i)))]);expected[normal_indices[-1]]-=value
    for b,c,d,e,value in identity['normal_base_quartic_constant']:
        normal_indices.append(qi[(29+b,29+c,29+d,29+e)]);expected[normal_indices[-1]]-=value
    assert len(normal_indices)==len(set(normal_indices))==90896
    product=np.array([[qi[tuple(sorted((*a,*b)))] for b in pairs] for a in pairs],dtype=np.int32)
    actual=np.zeros(len(quartics),dtype=np.int64)
    operations=0
    for r in range(134):
        a=np.flatnonzero(cofactor[r]);b=np.flatnonzero(full[r]);operations+=len(a)*len(b)
        np.add.at(actual,product[np.ix_(a,b)].ravel(),(cofactor[r,a][:,None]*full[r,b][None,:]).ravel())
    assert np.array_equal(actual%P,expected%P)
    inputs += [Path(__file__),ROOT/'code/s61_affine_seventeen_base_lift_20261103_v62.py',args.identity]
    report=dict(all_passed=True,scope='independent direct original134 homogeneous46-variable coefficient expansion',
                prime=P,full_variable_count=46,original_quadrics=134,
                independently_verified_quartic_coefficients=len(quartics),
                explicit_normal_terms=len(normal_indices),nonzero_cofactor_terms=int(np.count_nonzero(cofactor)),
                direct_scalar_products=operations,uses_scratch_inverse_cache=False,uses_reduction_algorithm=False,
                full_affine_exhaustion_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(inputs)},seconds=time.monotonic()-started)
    with args.output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print('PASS',len(quartics),operations,report['seconds'])


if __name__=='__main__':main()

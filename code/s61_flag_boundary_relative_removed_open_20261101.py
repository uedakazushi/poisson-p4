#!/usr/bin/env python3
"""Audit exact inputs to the relative dimension argument for a flag removed open."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'reports/S61_FLAG_BOUNDARY_RELATIVE_REMOVED_OPEN_20261101.json')
    args=parser.parse_args();start=time.time()
    names=['reports/S61_FLAG_BOUNDARY_RATIONAL_INITIAL_20261101.json',
        'reports/S61_FLAG_BOUNDARY_ORBIT_SPLIT_FORMAL_20261101.json',
        'reports/S61_FLAG_BOUNDARY_ORBIT_SPLIT_FORMAL_REVIEW_20261101.json',
        'reports/S61_FLAG_UNIFORM_OPEN_20261101.json']
    data=[json.loads((ROOT/name).read_text()) for name in names]
    checked=0
    for r in data:
        assert r.get('all_passed') or r.get('accepted')
        for key in ['input_sha256','source_sha256']:
            for rel,h in r.get(key,{}).items():
                assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==h,(rel,h)
                checked+=1
    initial,split,review,uniform=data
    assert initial['initial_subideal_dimension']==3 and initial['initial_subideal_degree']==5
    assert initial['transverse_length']==49
    assert split['primary_slice_free_rank']==5 and s.Rational(split['orbit_slice_minor'])==4
    geo=initial['geometry'];T=s.Matrix(geo['actual_flag_tangent']).applyfunc(s.Rational)
    L=s.Matrix(geo['b1_derivative']).applyfunc(s.Rational)
    assert T[:3,:].det()==8 and L*T!=s.zeros(1,3)
    assert geo['beta_rank']==16 and geo['joint_rank']==19
    assert geo['boundary_beta_rank']==16 and geo['boundary_joint_rank']==18
    note='research/notes/S61_FLAG_BOUNDARY_RELATIVE_REMOVED_OPEN_20261101.md'
    deps=names+[note,'code/s61_flag_boundary_relative_removed_open_20261101.py',
        'research/notes/S61_FLAG_BOUNDARY_RATIONAL_INITIAL_20261101.md',
        'research/notes/S61_FLAG_BOUNDARY_ORBIT_SPLIT_FORMAL_20261101.md']
    out=dict(all_passed=True,role='Dependency and tangent audit for the separate written proof; no independent six-jet replay.',
        special_formal_rank=5,special_transverse_length=49,base_dimension=15,relative_ring_dimension=18,
        relative_t_invertible_closure_dimension_upper=17,relative_saturation_boundary_dimension_upper=16,
        relative_flag_boundary_dimension=17,group_dimension=10,smooth_relative_dimension=9,
        ambient_flag_boundary_dimension=18,ambient_remaining_boundary_dimension_upper=17,
        specified_point_in_removed_ambient_open_proved=False,entire_general_beta_boundary_removed=False,
        conceptual_review=dict(reviewer='/root/mps_catalog',accepted=True,
            scope='Nakayama at the parameter prime, same-family lower5, dimension of the t=0 boundary, and smooth descent.'),
        input_hashes_checked=checked,input_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in deps},seconds=time.time()-start)
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n')
    print('RELATIVE FLAG REMOVED OPEN INPUT AUDIT PASS',out['seconds'])

if __name__=='__main__':main()

#!/usr/bin/env python3
"""Reproduce the full modular degree-three Macaulay rank with a bounded engine.

Compiles the adjacent C++ source in a temporary directory. No system packages
or libraries are changed. The 45-variable input is generated from immutable
S61 equations and the saved exact toric point.
"""
from __future__ import annotations

import argparse
import hashlib
from itertools import combinations_with_replacement
import json
from pathlib import Path
import subprocess
import tempfile

import sympy as s

from s61_affine_fiber_macaulay import load_fiber

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--seconds',type=int,default=150)
    parser.add_argument('--entry-limit',type=int,default=50000000)
    parser.add_argument('--prime',type=int,default=32003)
    args=parser.parse_args()
    assert s.isprime(args.prime) and args.prime>5
    assert 0<args.seconds<=300 and 0<args.entry_limit<=60000000
    if args.output.exists():raise FileExistsError(args.output)
    equations,columns,base,source,point_source=load_fiber()
    monomials=[mon for degree in (3,2,1,0)
               for mon in combinations_with_replacement(range(45),degree)]
    indices={mon:i for i,mon in enumerate(monomials)}
    modular=[]
    for poly in equations:
        converted={}
        for mon,c in poly.items():
            assert c.denominator%args.prime
            value=c.numerator*pow(c.denominator,-1,args.prime)%args.prime
            if value:converted[mon]=value
        modular.append(converted)
    cpp=ROOT/'code/s61_macaulay_rank.cpp'
    with tempfile.TemporaryDirectory(prefix='s61-degree3-') as temporary:
        directory=Path(temporary);binary=directory/'rank';matrix=directory/'matrix.txt'
        subprocess.run(['clang++','-O3','-std=c++17',str(cpp),'-o',str(binary)],check=True)
        with matrix.open('w') as out:
            print(args.prime,6164,len(monomials),file=out)
            print(*map(len,monomials),file=out)
            for multiplier in [None]+list(range(45)):
                for poly in modular:
                    row=sorted((indices[tuple(sorted(mon+(multiplier,))) if multiplier is not None else mon],c)
                               for mon,c in poly.items())
                    print(len(row),*sum(([key,value] for key,value in row),[]),file=out)
        matrix_hash=hashlib.sha256(matrix.read_bytes()).hexdigest()
        process=subprocess.run([str(binary),str(matrix),str(args.seconds),str(args.entry_limit)],
                               check=True,capture_output=True,text=True,timeout=args.seconds+30)
        result=json.loads(process.stdout)
    result['new_lower_degree_polynomials']=[[[list(monomials[index]),value] for index,value in row]
                                            for row in result.pop('new_lower_degree_rows')]
    complete=result['status']=='complete'
    if complete and args.prime==32003:
        assert result['rows_processed']==6164 and result['rank']==6073
        assert result['pivot_counts']==[0,0,135,5938]
        assert len(result['new_lower_degree_polynomials'])==1
    report={'id':'S61_MACAULAY_DEGREE3_COMPLETE_20261031','field':f'F_{args.prime}',
            'matrix_shape':[6164,len(monomials)],
            'row_order':'the 134 centered equations; then v_j times all 134 equations for j=0,...,44',
            'column_order':'degree 3,2,1,0; combinations_with_replacement(range(45),degree) within each degree',
            'weight_two_u_indices_zero_based':columns,'center_coordinates':[str(x) for x in base],
            'seconds_limit':args.seconds,'entry_limit':args.entry_limit,
            'matrix_sha256':matrix_hash,'computation':result,
            'limitations':['This is one finite-field specialization, not a characteristic-zero fiber degree.',
                           'The extra quadratic is not yet lifted over QQ or over varying first jets.',
                           'No affine finiteness or global exhaustion follows from these Macaulay ranks.'],
            'proof_note':'research/notes/S61_AFFINE_FIBER_REDUCTION_20261031.md',
            'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in (Path(__file__),cpp,ROOT/'code/s61_affine_fiber_macaulay.py',source,point_source)}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report,out,ensure_ascii=False,indent=2);out.write('\n')
    print(json.dumps({key:result[key] for key in ('status','rows_processed','rank','pivot_counts','elapsed_seconds')}))


if __name__=='__main__':main()

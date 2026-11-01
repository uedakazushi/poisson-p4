#!/usr/bin/env python3
"""Bounded resumable degree-three closure of the full 135-quadratic fiber.

The scratch snapshot remains under work; only compact evidence is written to
the requested report. Results concern the chosen prime-field specialization.
"""
import argparse
from fractions import Fraction
import hashlib
from itertools import combinations_with_replacement
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

import sympy as s
from s61_affine_fiber_macaulay import load_fiber

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--seconds',type=int,default=180)
    parser.add_argument('--entry-limit',type=int,default=60000000)
    parser.add_argument('--prime',type=int,default=32003)
    parser.add_argument('--snapshot',type=Path,required=True)
    parser.add_argument('--resume',type=Path)
    args=parser.parse_args()
    assert 5<args.prime<65536 and s.isprime(args.prime)
    assert 0<args.seconds<=300 and 0<args.entry_limit<=100000000
    if args.output.exists() or args.snapshot.exists():raise FileExistsError('Choose fresh report and snapshot paths')
    equations,columns,center,source,point_source=load_fiber()
    extra_source=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    extra=json.loads(extra_source.read_text())
    equations.append({tuple(mon):Fraction(c) for mon,c in extra['sample_centered_new_equation']})
    monomials=[mon for degree in (3,2,1,0) for mon in combinations_with_replacement(range(45),degree)]
    indices={mon:i for i,mon in enumerate(monomials)}
    cpp=ROOT/'code/s61_macaulay_closure_20261101.cpp'
    args.snapshot.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='s61-closure-') as temporary:
        temp=Path(temporary);binary=temp/'closure';matrix=temp/'input.txt'
        subprocess.run(['clang++','-O3','-std=c++17',str(cpp),'-o',str(binary)],check=True)
        with matrix.open('w') as out:
            print(args.prime,len(monomials),45,len(equations),file=out)
            for mon in monomials:print(len(mon),*mon,file=out)
            for poly in equations:
                row=[]
                for mon,c in poly.items():
                    assert c.denominator%args.prime
                    value=c.numerator*pow(c.denominator,-1,args.prime)%args.prime
                    if value:row.append((indices[mon],value))
                row.sort();print(len(row),*sum(([key,value] for key,value in row),[]),file=out)
        input_hash=hashlib.sha256(matrix.read_bytes()).hexdigest()
        if args.resume:
            metadata=json.loads(Path(str(args.resume)+'.meta.json').read_text())
            assert metadata['input_sha256']==input_hash
            assert metadata['engine_sha256']==hashlib.sha256(cpp.read_bytes()).hexdigest()
        temporary_snapshot=temp/'state.bin'
        command=[str(binary),str(matrix),str(args.seconds),str(args.entry_limit),str(temporary_snapshot),
                 str(args.resume.resolve()) if args.resume else '-']
        process=subprocess.run(command,check=True,stdout=subprocess.PIPE,text=True,timeout=args.seconds+60)
        result=json.loads(process.stdout)
        # Persist the result before copying a large optional resume file. Native
        # subprocesses can write their temporary directory in this environment;
        # the Python driver transfers the data into the repository scratch area.
        Path(str(args.snapshot)+'.result.json').write_text(process.stdout)
        shutil.copyfile(temporary_snapshot,args.snapshot)
        Path(str(args.snapshot)+'.meta.json').write_text(json.dumps({
            'input_sha256':input_hash,'engine_sha256':hashlib.sha256(cpp.read_bytes()).hexdigest()}))
    result['new_lower_degree_polynomials']=[[[list(monomials[index]),value] for index,value in row]
                                            for row in result.pop('new_lower_degree_rows')]
    with args.snapshot.open('rb') as snapshot_stream:
        snapshot_hash=hashlib.file_digest(snapshot_stream,'sha256').hexdigest()
    report={'id':'S61_MACAULAY_CLOSURE_20261101','field':f'F_{args.prime}',
            'starting_equations':135,'variables':45,'degree_cap':3,
            'input_sha256':input_hash,'computation':result,
            'snapshot':str(args.snapshot),'resumed_from':str(args.resume) if args.resume else None,
            'snapshot_sha256':snapshot_hash,
            'limits':{'seconds':args.seconds,'stored_entries':args.entry_limit},
            'q_to_u_indices_zero_based':columns,'center_coordinates':[str(v) for v in center],
            'limitations':['One finite-field fiber only; no characteristic-zero degree or exhaustion conclusion.',
                           'Only a complete status certifies closure under the implemented degree-three feedback.',
                           'A scratch binary snapshot is a resume aid, not an independent proof certificate.'],
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in (Path(__file__),cpp,source,point_source,extra_source)}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as out:json.dump(report,out,ensure_ascii=False,indent=2);out.write('\n')
    print(json.dumps({k:result[k] for k in ('status','round','processed','rank','pivot_counts','elapsed_seconds')}))


if __name__=='__main__':main()

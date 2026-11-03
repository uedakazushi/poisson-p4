#!/usr/bin/env python3
"""Check the entire projective line of the exact seventeen-base model.

Restricts all original134 quadrics to Y=s*v+a*w. A finite homogeneous
Macaulay rank computation is used instead of an affine Groebner computation.
The independent five-point evaluation supplies the stopping lower bound.
Explicit coefficients are F32003; generic QQ lifting uses these finite ranks.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import time

import numpy as np

from s61_affine_seventeen_base_lift_20261103_v62 import build
from s61_affine_normal_projection_20261102 import ROOT,P,pivot_minor


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def model(direction_count=2):
    high,cross,base,pairs,_,_,bpairs,v,inputs,_,_=build()
    vectors=[v]
    rng=np.random.default_rng(621127)
    vectors += [rng.integers(-5,6,size=17,dtype=np.int64) for _ in range(direction_count-1)]
    directions=np.column_stack(vectors)
    assert len(pivot_minor(directions)[0])==direction_count
    n=29+direction_count
    order=list(itertools.combinations_with_replacement(range(n),2))
    index={a:i for i,a in enumerate(order)}
    F=np.zeros((134,len(order)),dtype=np.int64)
    for k,pair in enumerate(pairs):F[:,index[pair]]=high[:,k]
    for b,vector in enumerate(vectors):
        linear=cross@vector%P
        for i in range(29):F[:,index[i,29+b]]=linear[:,i]
    for b,c in itertools.combinations_with_replacement(range(direction_count),2):
        if b==c:
            multipliers=[int(vectors[b][i])*int(vectors[b][j])%P for i,j in bpairs]
        else:
            multipliers=[(int(vectors[b][i])*int(vectors[c][j])+int(vectors[c][i])*int(vectors[b][j]))%P for i,j in bpairs]
        F[:,index[29+b,29+c]]=base@np.array(multipliers,dtype=np.int64)%P
    known_path=ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json'
    known=json.loads(known_path.read_text())
    assert known['all_passed'] and known['known_five_degree_two_evaluation_rank']==5
    points=[point+[0]*25+[1]+[0]*(direction_count-1) for point in known['known_five_mod_32003']]
    assert all(len(point)==n for point in points)
    for point in points:
        values=np.array([int(point[i])*int(point[j])%P for i,j in order],dtype=np.int64)
        assert not np.any(F@values%P)
    # t^2 times the four-plane quadratics are among degree-four monomials.
    # Their evaluation rank5 proves the quotient in degree4 has dimension>=5.
    evaluation=np.array([[int(point[0])**2*int(point[i])*int(point[j])%P
                         for i,j in itertools.combinations_with_replacement(range(4),2)]
                        for point in points],dtype=np.int64)
    _,evaluation_columns,evaluation_minor=pivot_minor(evaluation)
    assert len(evaluation_columns)==5
    assert len(pivot_minor(F)[0])==134
    return F,directions,points,evaluation_columns,evaluation_minor,inputs+[known_path]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--direction-count',type=int,default=2,choices=[1,2,3])
    parser.add_argument('--timeout',type=int,default=600)
    parser.add_argument('--reuse-input',type=Path)
    parser.add_argument('--reuse-output',type=Path)
    args=parser.parse_args();assert bool(args.reuse_input)==bool(args.reuse_output)
    args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    F,directions,points,ec,em,inputs=model(args.direction_count)
    n=29+args.direction_count
    path=args.output/'quadrics.in';native=args.output/'native.txt'
    path.write_text(f'{n} 134 32003\n'+' '.join(map(str,F.ravel()))+'\n')
    (args.output/'directions.json').write_text(json.dumps(directions.T.tolist(),indent=2)+'\n')
    source=ROOT/'code/s61_affine_section_overlap_rank_20261103_v62.cpp'
    if args.reuse_input:
        assert sha(path)==sha(args.reuse_input)
        shutil.copyfile(args.reuse_output,native)
        log=args.reuse_output.with_suffix('.log')
        if log.exists():shutil.copyfile(log,args.output/'native.log')
        mode='saved_native_rank_with_matching_input_and_fresh_five_evaluation'
    else:
        executable=args.output/'rank'
        flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],capture_output=True,text=True,check=True).stdout)
        subprocess.run(['c++','-std=c++17','-O3',str(source),*flags,'-o',str(executable)],check=True)
        with (args.output/'native.log').open('w')as stream:
            subprocess.run([str(executable.resolve()),str(path.resolve()),str(native.resolve()),'5'],stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=args.timeout)
        executable.unlink();mode='native_rank_replayed'
    lines=native.read_text().splitlines();tokens=lines[0].split()
    summary={k:int(v) for k,v in zip(tokens[::2],tokens[1::2])}
    assert summary['N']==n and summary['QUOTIENT4_BOUND']==5
    assert summary['PIVOT_PRODUCT']!=0
    assert summary['OVERLAP_RANK']==summary['RESIDUAL4']-5
    assert len(lines[2:])==summary['OVERLAP_RANK']
    overlap_columns=[int(line.split()[2]) for line in lines[2:]]
    assert len(set(overlap_columns))==len(overlap_columns)
    assert summary['RANK4']+5==len(list(itertools.combinations_with_replacement(range(n),4)))
    hf=[1,n,F.shape[1]-134,summary['FREE3'],5]
    inputs += [Path(__file__),source,ROOT/'code/s61_affine_seventeen_base_lift_20261103_v62.py',
               ROOT/'code/s61_affine_finite_normal_module_20261102.py',ROOT/'code/s61_affine_normal_projection_20261102.py',ROOT/'code/s61_affine_fiber_macaulay.py']
    report=dict(id='S61_AFFINE_BASE_LINEAR_SUPPORT_20261103_v62',all_passed=True,prime=P,mode=mode,
                original_quadrics=134,retained_base_variables=17,base_subspace_projective_dimension=args.direction_count-1,
                base_subspace_equation='Y=sum(direction_b*s_b); direction_0=v',
                projective_section_dimension=n-1,unknown_count=n,
                degree_four_row_rank=summary['RANK4'],degree_four_quotient_dimension=5,
                cubic_row_rank=summary['RANK3'],cubic_free_count=summary['FREE3'],
                canonical_quartic_residual_count=summary['RESIDUAL4'],overlap_rank=summary['OVERLAP_RANK'],
                overlap_echelon_pivot_product=summary['PIVOT_PRODUCT'],hilbert_function_through_degree_four=hf,
                known_five_mod_32003=points,known_five_degree_four_evaluation_columns=ec,
                known_five_degree_four_evaluation_minor=em,known_five_are_all_projective_solutions=True,
                affine_parameter_support='s0=1: all remaining s_b=0, t=1; length5 reduced',
                generic_QQ_section_length=5,full_affine_exhaustion_claimed=False,
                limitations=['Explicit field coefficients are frozen-beta F32003.',
                             'Finite rank and known-five ideal inclusion give generic QQ scheme equality in degree4 and after sheafification.',
                             'Only the entire specified projective base subspace is exhausted; the complementary base remains open.'],
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(inputs)},
                evidence_sha256={p.name:sha(p) for p in args.output.iterdir() if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',hf,summary,report['seconds'])


if __name__=='__main__':main()

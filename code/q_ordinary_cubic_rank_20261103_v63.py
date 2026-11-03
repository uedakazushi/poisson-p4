#!/usr/bin/env python3
"""Exact finite-field ranks of the complete true cubic initial-space blocks."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import shlex
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cubic',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    previous=json.loads((args.cubic/'report.json').read_text());assert previous['all_passed']
    matrix=args.cubic/'cubic-blocks.in';assert sha(matrix)==previous['matrix_sha256']
    native=ROOT/'code/q_ordinary_cubic_rank_20261103_v63.cpp';executable=args.output/'rank'
    flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],capture_output=True,text=True,check=True).stdout)
    subprocess.run(['c++','-std=c++17','-O3',str(native),*flags,'-o',str(executable)],check=True)
    evidence=args.output/'ranks.txt'
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(1000000000,1000000000))
        resource.setrlimit(resource.RLIMIT_CPU,(120,121))
    subprocess.run([str(executable.resolve()),str(matrix.resolve()),str(evidence.resolve())],check=True,timeout=120,preexec_fn=limits)
    lines=evidence.read_text().splitlines();blocks=[]
    for line in lines[:-1]:
        weight,columns,nW,nG,rW,rFull,new=map(int,line.split())
        expected=previous['block_shapes'][str(weight)]
        assert columns==expected['distinct_cubic_monomials'] and nW==expected['quadratic_multiples'] and nG==expected['new_cubic_rows']
        blocks.append(dict(old_weight=weight,columns=columns,quadruple_rows=nW,true_cubic_rows=nG,
                           quadratic_linear_multiple_rank=rW,full_cubic_initial_rank=rFull,
                           new_cubic_rank_mod_quadratic_multiples=new))
    total=lines[-1].split();assert total[0]=='TOTAL' and total[4]=='SECONDS'
    rW,rFull,new=map(int,total[1:4]);assert new==rFull-rW
    assert rW==sum(row['quadratic_linear_multiple_rank'] for row in blocks)
    assert rFull==sum(row['full_cubic_initial_rank'] for row in blocks)
    report=dict(id='Q_ORDINARY_CUBIC_RANK_20261103_V63',all_passed=True,coefficient_field='GF(32003)',
                all1995_linear_quadratic_multiples_included=True,
                all135_nonzero_constant_eliminated_cubic_initials_included=True,
                quadratic_linear_multiple_rank=rW,full_cubic_initial_rank=rFull,
                new_cubic_rank_mod_quadratic_multiples=new,
                ordinary_hilbert_values_through3_mod_p=[1,35,630-57,7770-rFull],
                block_ranks=blocks,
                QQ_rank_or_flatness_not_inferred_from_modular_rank=True,
                full_initial_ideal_or_local_dimension_computed=False,
                native_seconds=float(total[5]),
                input_sha256={str(p.resolve().relative_to(ROOT)):sha(p) for p in
                              [Path(__file__),native,args.cubic/'report.json',matrix]},
                ranks_sha256=sha(evidence),seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['input_sha256','block_ranks']},indent=2))


if __name__=='__main__':main()

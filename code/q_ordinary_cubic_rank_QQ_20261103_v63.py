#!/usr/bin/env python3
"""Bounded exact integer ranks of the small old-weight cubic blocks."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import json
from math import gcd,lcm
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
    parser.add_argument('--modular',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--timeout',type=int,default=120)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(3000000000,3000000000))
    packet=json.load(gzip.open(args.cubic/'cubic-initials.json.gz','rt'))
    original=json.load(gzip.open(ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz','rt'))
    model=json.load(gzip.open(ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz','rt'))
    modular=json.loads((args.modular/'report.json').read_text());assert modular['all_passed']
    ranks={row['old_weight']:row for row in modular['block_ranks']}
    blocks=defaultdict(lambda:dict(W=[],G=[]));weights=model['input_weights']
    for source in packet['quadratic_linear_multiples_recipe']['quadratic_original_rows']:
        q={(i,j):Fraction(c) for i,j,c in original['all350_quadratic_rows'][source]}
        for variable in range(35):
            blocks[model['residual_weight_bounds'][source]+weights[variable]]['W'].append(
                {tuple(sorted((*mon,variable))):c for mon,c in q.items()})
    for entry in packet['cubic_initial_relations']:
        if entry['terms']:
            blocks[model['residual_weight_bounds'][entry['source']]]['G'].append(
                {tuple(mon):Fraction(c) for mon,c in entry['terms']})
    matrix=args.output/'integer-blocks.in'
    with matrix.open('w') as stream:
        stream.write(str(len(blocks))+'\n')
        for weight,block in sorted(blocks.items()):
            rows=block['W']+block['G'];columns=sorted({mon for row in rows for mon in row});index={mon:i for i,mon in enumerate(columns)}
            r=ranks[weight]
            stream.write(f"{weight} {len(columns)} {len(block['W'])} {len(block['G'])} {r['quadratic_linear_multiple_rank']} {r['full_cubic_initial_rank']}\n")
            for row in rows:
                denominator=lcm(*(c.denominator for c in row.values()))
                integers={index[mon]:c.numerator*(denominator//c.denominator) for mon,c in row.items()}
                content=gcd(*integers.values())
                stream.write(str(len(integers))+' '+' '.join(f'{j} {c//content}' for j,c in sorted(integers.items()))+'\n')
    native=ROOT/'code/q_ordinary_cubic_rank_QQ_20261103_v63.cpp';executable=args.output/'rank'
    flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],capture_output=True,text=True,check=True).stdout)
    subprocess.run(['c++','-std=c++17','-O3',str(native),*flags,'-o',str(executable)],check=True)
    evidence=args.output/'ranks.txt';status='completed'
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(3000000000,3000000000))
        resource.setrlimit(resource.RLIMIT_CPU,(args.timeout,args.timeout+1))
    with (args.output/'native.log').open('w') as log:
        try:
            result=subprocess.run([str(executable.resolve()),str(matrix.resolve()),str(evidence.resolve())],stdout=log,stderr=subprocess.STDOUT,
                                  timeout=args.timeout,preexec_fn=limits)
            if result.returncode:status=f'exit_{result.returncode}'
        except subprocess.TimeoutExpired:status='timeout'
    lines=evidence.read_text().splitlines() if evidence.exists() else []
    completed=status=='completed' and lines and lines[-1].startswith('TOTAL ')
    rows=[]
    for line in lines:
        if line.startswith('TOTAL '):break
        weight,rW,rFull,new=map(int,line.split());rows.append(dict(old_weight=weight,QQ_quadratic_multiple_rank=rW,QQ_cubic_initial_rank=rFull,new_QQ_cubic_rank=new))
    report=dict(id='Q_ORDINARY_CUBIC_RANK_QQ_20261103_V63',all_blocks_completed=bool(completed),status=status,
                coefficient_field='QQ',method='FLINT exact integer rank after clearing each rational row denominator',
                block_ranks=rows,
                memory_limit_bytes=3000000000,wall_timeout_seconds=args.timeout,
                input_sha256={str(p.resolve().relative_to(ROOT)):sha(p) for p in
                              [Path(__file__),native,args.cubic/'cubic-initials.json.gz',args.modular/'report.json',
                               ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz',
                               ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz']},
                integer_matrix_sha256=sha(matrix),ranks_sha256=sha(evidence) if evidence.exists() else None,
                seconds=time.monotonic()-start)
    if completed:
        tokens=lines[-1].split();rW,rFull,new=map(int,tokens[1:4]);assert tokens[4]=='SECONDS'
        report.update(quadratic_linear_multiple_rank_QQ=rW,full_cubic_initial_rank_QQ=rFull,
                      new_cubic_rank_mod_quadratic_multiples_QQ=new,
                      ordinary_Hilbert_values_through3_QQ=[1,35,573,7770-rFull],
                      all_block_ranks_equal_GF32003=all(row['QQ_quadratic_multiple_rank']==ranks[row['old_weight']]['quadratic_linear_multiple_rank']
                                                     and row['QQ_cubic_initial_rank']==ranks[row['old_weight']]['full_cubic_initial_rank'] for row in rows),
                      native_seconds=float(tokens[5]))
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('input_sha256','block_ranks')},indent=2))


if __name__=='__main__':main()

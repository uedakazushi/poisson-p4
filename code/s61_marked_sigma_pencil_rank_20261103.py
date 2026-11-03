#!/usr/bin/env python3
"""Build bounded Macaulay certificates from the frozen all-h marked-square model.

All computed polynomial row modules in this script are over F_32003[h].
They give rank/support information, not characteristic-zero ambient removal.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
PRIME = 32003
WEIGHTS = [1]*7 + [2]*4 + [3]*2 + [1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@lru_cache(None)
def monomials(weights, degree):
    if not weights:
        return [()] if degree == 0 else []
    return [(k,)+tail for k in range(degree//weights[0]+1)
            for tail in monomials(weights[1:], degree-k*weights[0])]


def read_rows(path):
    reference = json.loads((ROOT/'reports/S61_MARKED_SIGMA_PENCIL_MODEL_20261102.json').read_text())
    assert sha(path) == reference['row_output_sha256']
    records = []
    with path.open() as source:
        assert tuple(map(int, source.readline().split())) == (14,136)
        for row in range(136):
            weight, count = map(int, source.readline().split())
            coefficients = defaultdict(dict)
            for _ in range(count):
                scalar, *powers = source.readline().split()
                powers = tuple(map(int,powers))
                assert len(powers) == 15
                if row == 135:
                    continue  # Observable t is not an additional equation.
                h, exponents = powers[0], powers[1:]
                assert sum(a*b for a,b in zip(exponents,WEIGHTS)) == weight
                scalar = Fraction(scalar)
                value = scalar.numerator * pow(scalar.denominator,-1,PRIME) % PRIME
                coefficients[exponents][h] = (coefficients[exponents].get(h,0)+value)%PRIME
            coefficients = {ex:{h:c for h,c in cs.items() if c} for ex,cs in coefficients.items()}
            records.append((row,weight,{ex:cs for ex,cs in coefficients.items() if cs}))
        assert not source.read().strip()
    return records


def make_module(records, mode, degree):
    kept = (list(range(7,13)) if mode == 'heavy' else
            list(range(13)) if mode == 'zero' else list(range(14)))
    deleted = set(range(14))-set(kept)
    weights = tuple(WEIGHTS[i] for i in kept)
    columns = sorted(monomials(weights,degree))
    index = {mon:i for i,mon in enumerate(columns)}
    rows = []
    for row,weight,source in records:
        if row == 135 or weight > degree:
            continue
        restricted = {tuple(ex[i] for i in kept):cs for ex,cs in source.items()
                      if not any(ex[i] for i in deleted)}
        if not restricted:
            continue
        for shift in monomials(weights,degree-weight):
            value = {index[tuple(a+b for a,b in zip(ex,shift))]:cs for ex,cs in restricted.items()}
            bound = max(max(cs) for cs in value.values())
            rows.append((bound,len(value),row,shift,value))
    rows.sort(key=lambda record:record[:4])
    return columns,[record[-1] for record in rows]


def write_module(path, columns, rows, extra_columns=()):
    with path.open('w') as output:
        print(PRIME,len(rows)+len(extra_columns),len(columns),0,file=output)
        print(0,1,file=output)
        for row in [*rows,*[{column:{0:1}} for column in extra_columns]]:
            print(len(row),file=output)
            for column,cs in sorted(row.items()):
                degree = max(cs)
                print(column,degree+1,*[cs.get(i,0) for i in range(degree+1)],file=output)


def compile_native(source, executable):
    flags = shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],
                                      capture_output=True,text=True,check=True).stdout)
    subprocess.run(['c++','-O3','-std=c++17',str(source),*flags,'-o',str(executable)],check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--tasks',nargs='+',default=['heavy:6','heavy:8','heavy:9','heavy:10'])
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False,parents=True)
    records = read_rows(args.rows)
    span_source = ROOT/'code/s61_marked_sigma_polynomial_rowmodule_20261103.cpp'
    numeric_source = ROOT/'code/s61_marked_sigma_numeric_complement_20261103.cpp'
    span,numeric = args.output/'rowmodule',args.output/'numeric_complement'
    compile_native(span_source,span)
    compile_native(numeric_source,numeric)
    results = []
    for task in args.tasks:
        mode,degree = task.split(':');degree=int(degree)
        assert mode in ['zero','full','heavy']
        columns,rows = make_module(records,mode,degree)
        stem = f'{mode}-{degree}'
        input_path = args.output/f'{stem}.in'
        write_module(input_path,columns,rows)
        extras = []
        if mode == 'full':
            values = list(map(int,subprocess.run([str(numeric),str(input_path),'0'],
                                               capture_output=True,text=True,check=True).stdout.split()))
            rank,extras = values[0],values[1:]
            assert len(extras) == len(columns)-rank == 4
            write_module(input_path,columns,rows,extras)
        output_path,log = args.output/f'{stem}.out',args.output/f'{stem}.log'
        started=time.monotonic()
        with log.open('w') as stream:
            result = subprocess.run([str(span),str(input_path),str(output_path)],stdout=stream,
                                    stderr=subprocess.STDOUT,timeout=185,check=True)
        results.append(dict(task=task,shape=[len(rows)+len(extras),len(columns)],
                            extra_cokernel_monomials=[columns[i] for i in extras],
                            output=output_path.read_text(),seconds=time.monotonic()-started,
                            returncode=result.returncode))
        print(task,results[-1]['shape'],results[-1]['output'],flush=True)
    span.unlink();numeric.unlink()
    inputs=[Path(__file__),span_source,numeric_source,
            ROOT/'reports/S61_MARKED_SIGMA_PENCIL_MODEL_20261102.json']
    report=dict(all_passed=True,field='F32003[h]',global_rows_sha256=sha(args.rows),results=results,
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs},
                evidence_sha256={str(p.relative_to(args.output)):sha(p)
                                 for p in args.output.rglob('*') if p.is_file()},
                full_pencil_removal_proved=False)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()

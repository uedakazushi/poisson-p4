#!/usr/bin/env python3
"""Prove that the 274 explicit AB quadrics generate the full quartic ideal component."""
import argparse
import hashlib
import json
import math
import subprocess
import tempfile
import time
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
NOTE='research/notes/AB_QUADRICS_QUARTIC_MULTIPLICATION_20261102.md'

def verify():
    clock=time.monotonic();pr=32003
    qp=ROOT/'reports/AB_EXPLICIT_QUADRICS_CUBIC_GENERATION_20261102.json';hp=ROOT/'reports/AB_ORIGIN_DEGREE_FOUR_20261101.json';src=ROOT/'code/ab_quadrics_quartic_multiplication_20261102.cpp'
    q=json.loads(qp.read_text());h=json.loads(hp.read_text());assert q['all_passed'] and h['all_passed']
    assert sorted(q['coordinate_weights'])==sorted(h['coordinate_weights']);qq=q['quadratic_relations'];g=h['evaluation_certificate']['groups']
    assert len(qq)==274 and sum(x['rank'] for x in g)==38570
    with tempfile.TemporaryDirectory(prefix='ab-fourth-multiplication-') as tmp:
        wd=Path(tmp);exe=wd/'rank';subprocess.run(['c++','-std=c++17','-O3',str(src),'-o',str(exe)],check=True,timeout=50)
        data=wd/'input';
        with data.open('w') as f:
            print(45,file=f)
            for w in q['coordinate_weights']:print(*w,file=f)
            print(len(qq),file=f)
            for row in qq:
                print(len(row),file=f)
                for i,j,c in row:
                    v=Fraction(c);print(i,j,v.numerator*pow(v.denominator,-1,pr)%pr,file=f)
            print(len(g),file=f)
            for b in g:print(*b['weight'],b['columns'],b['rank'],file=f)
        with data.open() as f:run=subprocess.run([str(exe)],stdin=f,capture_output=True,text=True,check=True,timeout=180)
        sh=hashlib.sha256(data.read_bytes()).hexdigest()
    rows=[list(map(int,line.split())) for line in run.stdout.splitlines()];assert len(rows)==len(g)==249
    assert len({(a,b) for a,b,_,_,_,_ in rows})==249 and all(n==m for a,b,c,k,n,m in rows)
    assert sum(c for a,b,c,k,n,m in rows)==math.comb(48,4)==194580
    assert sum(k for a,b,c,k,n,m in rows)==274*1035==283590
    assert sum(n for a,b,c,k,n,m in rows)==194580-38570==156010
    paths=[Path(__file__),src,qp,hp,ROOT/NOTE]
    return dict(id='AB_QUADRICS_QUARTIC_MULTIPLICATION_20261102',all_passed=True,prime=pr,weight_groups=249,
      quadratic_times_quadratic_monomials=283590,quartic_polynomial_dimension=194580,
      degree_four_of_AB=38570,quartic_ideal_dimension=156010,minimal_quartic_generators=0,
      remaining_possible_new_generator_degrees=[5,6],rank_table_columns=['weight1','weight2','columns','rows','rank_lower_bound','target'],
      rank_table=rows,matrix_description_sha256=sh,input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths},seconds=time.monotonic()-clock)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();result=verify()
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print('PASS',result['seconds'])

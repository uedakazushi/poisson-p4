#!/usr/bin/env python3
"""Certify rho=1 support at the new and old quadratic core values.

All three computations use actual finite fields. Adoption is permitted only
when regenerated input bytes and separately observed exit-status records match.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
PRIME=32003
NAMES=['new','old-plus','old-minus']


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def sources(path):
    ref=json.loads((ROOT/'reports/S61_MARKED_SIGMA_PENCIL_MODEL_20261102.json').read_text())
    assert sha(path)==ref['row_output_sha256']
    arrays={name:[]for name in NAMES}
    with path.open()as stream:
        assert tuple(map(int,stream.readline().split()))==(14,136)
        for i in range(136):
            weight,count=map(int,stream.readline().split());out={name:defaultdict(lambda:[0,0])for name in NAMES}
            for _ in range(count):
                scalar,*powers=stream.readline().split();e=list(map(int,powers));value=Fraction(scalar)
                value=value.numerator*pow(value.denominator,-1,PRIME)%PRIME;mon=tuple(e[1:14]);hp=e[0]
                out['new'][mon][hp%2]=(out['new'][mon][hp%2]+value*pow(-9314,hp//2,PRIME))%PRIME
                for name,root in [('old-plus',3083),('old-minus',28920)]:
                    out[name][mon][0]=(out[name][mon][0]+value*pow(root,hp,PRIME))%PRIME
            for name in NAMES:arrays[name].append({ex:cs for ex,cs in out[name].items()if any(cs)})
        assert not stream.read().strip()
    def polynomial(row):
        return '+'.join('('+str(c0)+'+'+str(c1)+'*a)'+''.join('*v'+str(i)+('^'+str(j)if j>1 else'')
                        for i,j in enumerate(ex)if j)for ex,(c0,c1)in row.items())or'0'
    result={}
    for name in NAMES:
        head=('A=(ZZ/32003)[h];K=GF(A/ideal(h^2+9314),Variable=>a);a=K_0;assert(a^2+9314==0);'
              if name=='new'else'K=ZZ/32003;a=0;')
        text=head+'R=K['+','.join('v'+str(i)for i in range(13))+',MonomialOrder=>{Weights=>{1,1,1,1,1,1,1,2,2,2,2,3,3},GRevLex}];\n'
        text+='I=ideal('+','.join(polynomial(row)for row in arrays[name][:135])+');T='+polynomial(arrays[name][135])+';print "INPUT READY";<<flush;G=gb I;d=dim(R/I);print("DIM="|toString d);if d==0 then(print("LENGTH="|toString numColumns basis(R/I));J=saturate(I,ideal vars R);print("EXTRA_LENGTH="|toString numColumns basis(R/J));print("J="|toString gens gb J));print("T="|toString(T%G));print("T2="|toString(T^2%G));print "FIBER COMPLETE";exit 0;\n'
        result[name]=text
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--rows',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--adopt',type=Path)
    parser.add_argument('--timeout',type=int,default=180);args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    scripts=sources(args.rows);m2=shutil.which('M2');assert m2
    if args.adopt:
        statuses=dict((name,int(code))for name,code in(map(str.split,(args.adopt/'rho-one-special-status.txt').read_text().splitlines())))
        assert all(statuses.get(name)==0 for name in NAMES)
        shutil.copyfile(args.adopt/'rho-one-special-status.txt',args.output/'observed-status.txt')
    results=[]
    for name,text in scripts.items():
        script=args.output/f'rho-one-{name}.m2';script.write_text(text);log=args.output/f'rho-one-{name}.log';started=time.monotonic()
        if args.adopt:
            assert sha(script)==sha(args.adopt/script.name),name
            shutil.copyfile(args.adopt/log.name,log);returncode=statuses[name]
        else:
            with log.open('w')as stream:
                result=subprocess.run([m2,'--script',str(script.resolve())],stdout=stream,stderr=subprocess.STDOUT,
                                      timeout=args.timeout,check=True)
            returncode=result.returncode
        content=log.read_text()
        for line in ['INPUT READY','DIM=0','LENGTH=4','EXTRA_LENGTH=0','J=matrix {{1}}','T=0','T2=0','FIBER COMPLETE']:
            assert line in content.splitlines(),(name,line)
        assert 'error:'not in content
        results.append(dict(case=name,returncode=returncode,dimension=0,length=4,
                            origin_away_saturation_ideal_one=True,only_reduced_point='standard flag origin',
                            t_zero_in_slice=True,execution='adopted exact matching input and observed exit0'if args.adopt else'fresh subprocess',
                            seconds=time.monotonic()-started))
    paths=[Path(__file__),ROOT/'reports/S61_MARKED_SIGMA_PENCIL_MODEL_20261102.json',
           ROOT/'reports/S61_MARKED_SIGMA_STANDARD_FLAG_SPREAD_20261103_v62.json',
           ROOT/'research/notes/S61_MARKED_SIGMA_RHO_ONE_QUADRATICS_20261103_v62.md']
    report=dict(id='S61_MARKED_SIGMA_RHO_ONE_QUADRATICS_20261103_v62',all_passed=True,
                field='GF(32003^2) for h^2+9314; GF(32003) for split roots3083,28920',
                global_rows_sha256=sha(args.rows),cases=results,
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
                evidence_sha256={str(p.relative_to(args.output)):sha(p)for p in args.output.rglob('*')if p.is_file()},
                limitations=['Slice support and slice t=0 alone are not an ambient nilpotence proof.',
                             'Use the separately verified all-h standard flag mixed-characteristic neighborhood.'])
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',report['id'])


if __name__=='__main__':main()

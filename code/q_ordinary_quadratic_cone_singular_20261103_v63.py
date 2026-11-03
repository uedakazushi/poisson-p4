#!/usr/bin/env python3
"""Bounded Singular slimgb replay of the frozen exact ordinary quadratic cone."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
PRIME=32003


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cone',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--runtime',type=Path,default=ROOT/'work/nonsquare-v62/singular-runtime/run-singular.sh')
    parser.add_argument('--timeout',type=int,default=200)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    archive=args.cone/'quadrics.json.gz';data=json.load(gzip.open(archive,'rt'))
    extraction=json.loads((args.cone/'extraction.json').read_text());assert extraction['extraction_completed']
    assert sha(archive)==extraction['quadratic_archive_sha256']
    variables=data['variables'];assert len(variables)==35
    script=args.output/'compute.sing';summary=(args.output/'summary.txt').resolve()
    text=f'ring r={PRIME},('+','.join(variables)+'),dp;\nideal I=\n'
    quadrics=[]
    for row in data['selected_original_rows']:
        terms=[]
        for i,j,c in data['all350_quadratic_rows'][row]:
            value=Fraction(c);assert value.denominator%PRIME
            residue=value.numerator*pow(value.denominator,-1,PRIME)%PRIME
            terms.append(f'{residue}*{variables[i]}*{variables[j]}')
        quadrics.append('+'.join(terms))
    text+=',\n'.join(quadrics)+';\noption(prot);\nideal G=slimgb(I);\n'
    text+='link f=":w '+str(summary)+'";\n'
    text+='write(f,"DIMENSION "+string(dim(G)));\n'
    text+='write(f,"GB_GENERATORS "+string(size(G)));\n'
    text+='write(f,"HILBERT_NUMERATOR "+string(hilb(G,1)));\n'
    text+='close(f);\nquit;\n';script.write_text(text)
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(6000000000,6000000000))
        resource.setrlimit(resource.RLIMIT_CPU,(args.timeout,args.timeout+1))
    status='completed'
    with (args.output/'singular.log').open('w') as log:
        try:
            result=subprocess.run([str(args.runtime.resolve()),'--no-tty','-q',str(script.resolve())],
                                  stdout=log,stderr=subprocess.STDOUT,timeout=args.timeout,preexec_fn=limits)
            if result.returncode:status=f'exit_{result.returncode}'
        except subprocess.TimeoutExpired:status='timeout'
    resulttext=summary.read_text() if summary.exists() else ''
    completed=status=='completed' and 'DIMENSION ' in resulttext
    report=dict(id='Q_ORDINARY_QUADRATIC_CONE_SINGULAR_20261103_V63',cone_GB_completed=completed,status=status,
                prime=PRIME,ordinary_local_variables=35,exact_QQ_quadratic_span_rank=extraction['quadratic_span_rank_QQ'],
                method='Singular slimgb over GF(32003)',wall_timeout_seconds=args.timeout,memory_limit_bytes=6000000000,
                full_formal_model_exact_dimension_claimed=False,
                input_sha256={str(p.resolve().relative_to(ROOT)):sha(p) for p in
                              [Path(__file__),archive,args.cone/'extraction.json',args.runtime]},
                script_sha256=sha(script),seconds=time.monotonic()-start)
    if completed:
        report['finite_field_quadratic_cone_dimension']=int(next(line.split()[1] for line in resulttext.splitlines() if line.startswith('DIMENSION ')))
        report['cone_GB_generators']=int(next(line.split()[1] for line in resulttext.splitlines() if line.startswith('GB_GENERATORS ')))
        report['cone_hilbert_numerator']=next(line[len('HILBERT_NUMERATOR '):] for line in resulttext.splitlines() if line.startswith('HILBERT_NUMERATOR '))
        report['summary_sha256']=sha(summary)
        report['dimension_transfer_requires_proper_projectivization_argument']=True
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='input_sha256'},indent=2))


if __name__=='__main__':main()

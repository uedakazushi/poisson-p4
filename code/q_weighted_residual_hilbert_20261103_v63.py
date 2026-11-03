#!/usr/bin/env python3
"""Bounded weighted-Hilbert computation on an explicitly selected residual ideal.

A weight cutoff below 19 gives a SUBideal of the full residual ideal.  Its
dimension is therefore only an UPPER BOUND for the fixed-a full model.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expansion',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--weight-cut',type=int,default=8)
    parser.add_argument('--timeout',type=int,default=150)
    parser.add_argument('--prepare-only',action='store_true')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    modelpath=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz'
    model=json.load(gzip.open(modelpath,'rt'))
    expansion=json.loads((args.expansion/'report.json').read_text())
    assert expansion['expansion_completed']
    selected=[i for i,(weight,length) in enumerate(zip(model['residual_weight_bounds'],expansion['term_counts']))
              if weight<=args.weight_cut and length]
    selected.sort(key=lambda i:(model['residual_weight_bounds'][i],i))
    selectedset=set(selected)
    script=args.output/'compute.m2'
    with (args.expansion/'residuals.m2').open() as source,script.open('w') as dest:
        dest.write(source.readline())
        count=0
        for line in source:
            if line.startswith('q'):
                row=int(line[1:line.index('=')])
                if row in selectedset:
                    dest.write(line);count+=1
        assert count==len(selected)
        dest.write('I=ideal{'+','.join(f'q{i}' for i in selected)+'};\n')
        dest.write('print "LOADED"; gbTrace=1; G=gb I; print "GB_FINISHED";\n')
        dest.write('print ("GB_GENERATORS "|toString numgens source gens G);\n')
        dest.write('print ("DIMENSION "|toString dim(R/I));\n')
        dest.write('print ("CODIMENSION "|toString codim I);\n')
        dest.write('print ("HILBERT "|toString hilbertSeries(R/I,Reduce=>true));\n')
        dest.write('print ("LEADING_MONOMIALS "|toString flatten entries leadTerm gens G);\nexit 0;\n')
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(6000000000,6000000000))
        resource.setrlimit(resource.RLIMIT_CPU,(args.timeout,args.timeout+1))
    status='prepared' if args.prepare_only else 'completed';returncode=None
    with (args.output/'m2.log').open('w') as log:
        if args.prepare_only:
            log.write('Prepared exact selected input without starting a GB computation.\n')
        else:
            try:
                result=subprocess.run(['M2','--script',str(script)],stdout=log,stderr=subprocess.STDOUT,
                                      timeout=args.timeout,preexec_fn=limits)
                returncode=result.returncode
                if returncode:status=f'exit_{returncode}'
            except subprocess.TimeoutExpired:status='timeout'
    log=(args.output/'m2.log').read_text()
    completed=status=='completed' and 'GB_FINISHED' in log and 'DIMENSION ' in log
    report=dict(id='Q_WEIGHTED_RESIDUAL_HILBERT_20261103_V63',completed=completed,status=status,
                coefficient_field=expansion['coefficient_field'],a_rational=expansion['a_rational'],
                positive_variables=34,weight_cut=args.weight_cut,selected_rows=selected,
                selected_nonzero_polynomials=len(selected),
                selected_total_terms=sum(expansion['term_counts'][i] for i in selected),
                full_residual_ideal_included=args.weight_cut>=19,
                partial_ideal_dimension_is_only_upper_bound=args.weight_cut<19,
                characteristic_zero_dimension_claimed=False,
                memory_limit_bytes=6000000000,wall_timeout_seconds=args.timeout,
                input_sha256={str(p.resolve().relative_to(ROOT)):sha(p)
                              for p in [Path(__file__),modelpath,args.expansion/'report.json',args.expansion/'residuals.m2']},
                script_sha256=sha(script),seconds=time.monotonic()-started)
    if completed:
        for label,key in [('DIMENSION','dimension'),('CODIMENSION','codimension'),('GB_GENERATORS','gb_generators')]:
            report[key]=int(next(line.split()[1] for line in log.splitlines() if line.startswith(label+' ')))
        report['hilbert_series']=next(line[len('HILBERT '):] for line in log.splitlines() if line.startswith('HILBERT '))
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('input_sha256','selected_rows')},indent=2))


if __name__=='__main__':
    main()

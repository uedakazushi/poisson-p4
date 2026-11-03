#!/usr/bin/env python3
"""Safe full-ideal dimension bounds from four Lex leading-generator ideals."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expansion',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--timeout',type=int,default=90)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    modelpath=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz'
    model=json.load(gzip.open(modelpath,'rt'));expansion=json.loads((args.expansion/'report.json').read_text())
    assert expansion['expansion_completed']
    weights=args.output/'weights.in';weights.write_text(' '.join(map(str,model['input_weights'][1:]))+'\n')
    native=ROOT/'code/q_weighted_residual_leading_scan_20261103_v63.cpp';executable=args.output/'scan'
    subprocess.run(['c++','-O3','-std=c++17',str(native),'-o',str(executable)],check=True)
    output=args.output/'leading.txt'
    subprocess.run([str(executable),str(args.expansion/'residuals.sparse'),str(weights),str(output)],check=True,timeout=args.timeout)
    lines=output.read_text().splitlines();assert int(lines[0])==4
    orders=[list(map(int,line.split())) for line in lines[1:5]]
    records=[list(map(int,line.split())) for line in lines[5:]]
    assert len(records)==1400
    polynomials=[[] for _ in range(4)]
    for order,row,coefficient,*exponents in records:
        assert len(exponents)==34
        assert bool(coefficient)==bool(expansion['term_counts'][row])
        if coefficient:
            assert sum(a*b for a,b in zip(exponents,model['input_weights'][1:]))==model['residual_weight_bounds'][row]
            polynomials[order].append('*'.join(f'z{i}^{e}' for i,e in enumerate(exponents) if e))
    summary=(args.output/'summary.txt').resolve();script=args.output/'compute.m2'
    text='R=GF(32003)['+','.join(f'z{i}' for i in range(34))+'];\nf=openOut "'+str(summary)+'";\n'
    for order,polys in enumerate(polynomials):
        text+='J=monomialIdeal('+','.join(polys)+');\n'
        text+=f'f << "ORDER {order} DIMENSION " << dim(R/J) << endl;\n'
        text+='flush f;\n'
    text+='close f;\nexit 0;\n';script.write_text(text)
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(6000000000,6000000000))
        resource.setrlimit(resource.RLIMIT_CPU,(args.timeout,args.timeout+1))
    status='completed'
    with (args.output/'m2.log').open('w') as log:
        try:
            p=subprocess.run(['M2','--script',str(script)],stdout=log,stderr=subprocess.STDOUT,
                             timeout=args.timeout,preexec_fn=limits)
            if p.returncode:status=f'exit_{p.returncode}'
        except subprocess.TimeoutExpired:status='timeout'
    results={}
    if summary.exists():
        for line in summary.read_text().splitlines():
            tokens=line.split();assert tokens[0]=='ORDER' and tokens[2]=='DIMENSION'
            results[int(tokens[1])]=int(tokens[3])
    report=dict(id='Q_WEIGHTED_RESIDUAL_LEADING_SCAN_20261103_V63',all_four_completed=len(results)==4,status=status,
                coefficient_field='GF(32003)',a_rational=expansion['a_rational'],positive_variables=34,
                all350_full_residuals_scanned=True,nonzero_residuals=expansion['nonzero_rows'],
                lex_variable_orders=orders,leading_generator_ideal_dimensions=results,
                full_model_dimension_upper_bound=min(results.values()) if results else None,
                claim='For each order, J=(LM(F_0),...,LM(F_349)) is a subideal of the full initial ideal; dim(full model)<=dim(R/J).',
                full_Groebner_basis_or_exact_dimension_claimed=False,
                characteristic_zero_dimension_claimed=False,
                input_sha256={str(p.resolve().relative_to(ROOT)):sha(p) for p in
                              [Path(__file__),native,modelpath,args.expansion/'report.json',args.expansion/'residuals.sparse']},
                evidence_sha256={p.name:sha(p) for p in [output,script,summary] if p.exists()},
                seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['input_sha256','lex_variable_orders','evidence_sha256']},indent=2))


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Choose explicit17 minimal true cubic initial generators modulo quadrics."""
import argparse
import gzip
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
    parser.add_argument('--cubic',type=Path,default=ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    packet=json.load(gzip.open(args.cubic/'cubic-initials.json.gz','rt'))
    modelpath=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz'
    model=json.load(gzip.open(modelpath,'rt'))
    QQ=json.loads((args.cubic/'exact-rank/report.json').read_text());assert QQ['all_blocks_completed'] and QQ['all_block_ranks_equal_GF32003']
    modular=json.loads((args.cubic/'modular-rank/report.json').read_text());assert modular['all_passed']
    oldrank={row['old_weight']:row for row in modular['block_ranks']}
    groups={}
    for row in packet['cubic_initial_relations']:
        if row['terms']:groups.setdefault(model['residual_weight_bounds'][row['source']],[]).append(row)
    matrix=args.output/'cubic-blocks.in';matrix.write_bytes(gzip.decompress((args.cubic/'cubic-blocks.in.gz').read_bytes()))
    native=ROOT/'code/q_ordinary_cubic_representatives_20261103_v63.cpp';executable=args.output/'representatives'
    flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],capture_output=True,text=True,check=True).stdout)
    subprocess.run(['c++','-O3','-std=c++17',str(native),*flags,'-o',str(executable)],check=True)
    evidence=args.output/'selection.txt'
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(1000000000,1000000000))
        resource.setrlimit(resource.RLIMIT_CPU,(60,61))
    subprocess.run([str(executable.resolve()),str(matrix.resolve()),str(evidence.resolve())],check=True,timeout=60,preexec_fn=limits)
    selected=[];counts={}
    for line in evidence.read_text().splitlines():
        tokens=line.split()
        if tokens[0]=='BLOCK':
            weight,rank,count=map(int,tokens[1:]);assert rank==oldrank[weight]['quadratic_linear_multiple_rank']
            assert count==oldrank[weight]['new_cubic_rank_mod_quadratic_multiples'];counts[weight]=count
        elif tokens[0]=='TOTAL':assert tokens[1]=='17'
        else:
            weight,index,pivot,coefficient=map(int,tokens);assert coefficient%32003
            row=groups[weight][index]
            selected.append(dict(source=row['source'],old_weight=weight,terms=row['terms'],
                                 modular_quotient_pivot_column=pivot,modular_unit_pivot=coefficient))
    assert len(selected)==len({row['source'] for row in selected})==17
    result=dict(variables=packet['variables'],quadratic_generators=packet['quadratic_linear_multiples_recipe']['quadratic_original_rows'],
                cubic_generators=selected)
    artifact=args.output/'generators.json';artifact.write_text(json.dumps(result,indent=2)+'\n')
    report=dict(id='Q_ORDINARY_CUBIC_REPRESENTATIVES_20261103_V63',all_passed=True,
                selected_cubic_source_rows=[row['source'] for row in selected],new_generators_by_old_weight=counts,
                minimal_generators_through3=dict(quadratic=57,cubic=17,total=74),
                independently_modular_selected_classes_form_basis=True,
                QQ_basis_follows_from_exact_QQ_dimension17_and_modular_independence=True,
                V_basis_follows_from_previous_equal_QQ_modp_ranks_and_Nakayama=True,
                exact_ideal_statement='J_(degrees<=3)=(57 selected Q_i,17 selected G_i) over QQ and Z_(32003)',
                next_linearQ_plus_cubic_relation_count=35*57+17-1911,
                full_initial_ideal_or_completed_local_ideal_generators_not_claimed=True,
                input_sha256={str(p.resolve().relative_to(ROOT)):sha(p) for p in
                              [Path(__file__),native,modelpath,args.cubic/'cubic-initials.json.gz',args.cubic/'cubic-blocks.in.gz',
                               args.cubic/'exact-rank/report.json',args.cubic/'modular-rank/report.json']},
                evidence_sha256=dict(generators_json=sha(artifact),selection_txt=sha(evidence)),seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('input_sha256','evidence_sha256')},indent=2))


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Save the completed bounded exact generic-projection attempts and blockers.

This records actual solver failures, not a rank/unit theorem. It checks the
small-parameter inputs against the independently verified reduction producer.
No solver is rerun and no runtime configuration is changed.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import time

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def run(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    work=ROOT/'work/affine-infinity-lift-v62'
    reduction=ROOT/'reports/S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_v62'
    prior=json.loads((reduction/'report.json').read_text());assert prior['all_passed']
    for name,value in prior['input_sha256'].items():assert sha(ROOT/name)==value
    for name,value in prior['evidence_sha256'].items():assert sha(reduction/name)==value
    def ideal(text):return text.split('ideal I=',1)[1].split(';\n',1)[0]
    reference={mode:ideal((reduction/(mode+'-exact.sing')).read_text())for mode in ['secant','wronski']}
    specs=[('generic-pair-QQ',18,43,134,300,124,'time_limit',None),
           ('generic-pair-homogeneous-v2',18,44,134,300,124,'time_limit',None),
           ('secant-nine-QQ',9,45,136,600,14,'address_space_limit','secant'),
           ('wronski-eight-QQ',8,45,136,600,14,'address_space_limit','wronski'),
           ('wronski-eight-sparse-QQ',8,45,136,600,14,'address_space_limit','wronski')]
    attempts=[]
    for name,nparams,nvars,neq,cap,code,status,mode in specs:
        script=work/(name+'.sing');log=work/(name+'.log');text=script.read_text();content=log.read_text()
        header=re.search(r'ring \w+=\(0,([^)]*)\),\(([^)]*)\),dp;',text);assert header
        params=header.group(1).split(',');variables=header.group(2).split(',')
        assert len(params)==nparams and len(variables)==nvars and len(ideal(text).split(','))==neq
        if status=='address_space_limit':assert 'no more memory' in content and 'halt 14' in content
        else:assert 'START' in content and 'DONE' not in content and 'CANDIDATE_READY' not in content
        if mode:assert ideal(text)==reference[mode]
        raw={}
        for path in [script,log]:
            destination=output/(path.name+'.gz')
            with gzip.open(destination,'wb')as stream:stream.write(path.read_bytes())
            assert gzip.decompress(destination.read_bytes())==path.read_bytes()
            raw[path.name]=sha(path)
        attempts.append(dict(name=name,coefficient_parameter_count=nparams,ring_variable_count=nvars,
                             input_generators=neq,algorithm='homogeneous modStd candidate plus exact HF guards'if'homogeneous'in name else'exact direct slimgb',
                             time_cap_seconds=cap,address_space_cap_KiB=6291456,returncode=code,status=status,
                             verified_reduced_system_match=bool(mode),candidate_obtained=False,rank_or_unit_result=False,
                             source_and_log_raw_sha256=raw))
    # Earlier syntax-only attempt is retained separately; it is not counted as
    # an algorithm timeout or as evidence of a mathematical outcome.
    syntax=work/'generic-pair-homogeneous.log'
    if syntax.exists():
        assert 'error occurred' in syntax.read_text()
        with gzip.open(output/'homogeneous-initial-syntax.log.gz','wb')as stream:stream.write(syntax.read_bytes())
    runtime=ROOT/'work/nonsquare-v62/singular-runtime/run-singular.sh';assert runtime.is_file()
    (output/'resume.txt').write_text('Use the already installed runtime, without setup changes:\n'+str(runtime)+' --cpus=0 -q '+str(reduction/'wronski-exact.sing')+'\n\nThe completed direct runs exhausted6GiB address space. An improved coefficient-control method or a separately coordinated larger exclusive-memory budget is needed; no numerical rank bound is known.\n')
    memory={}
    for name in ['memory.max','memory.current','memory.events']:
        path=Path('/sys/fs/cgroup')/name
        if path.exists():memory[name]=path.read_text().strip()
    report=dict(id='S61_AFFINE_GENERIC_PROJECTION_ATTEMPTS_20261103_v62',all_recorded=True,
                runtime='reused same-session Singular4.4.1 wrapper; direct jobs used --cpus=0',
                runtime_wrapper_sha256=sha(runtime),attempts=attempts,cgroup_memory_at_save=memory,
                sparse_order_choice=dict(first_coordinates=[8,24],coefficient_term_counts=[1229,1206],
                                         prior_first_coordinates=[0,1],prior_coefficient_term_counts=[2240,2541],
                                         linear_elimination_minor='16*z^31'),
                generic_rank_upper29_proved=False,generic_unit_proved=False,affine_surfaces_excluded=False,
                unaffected_work='Exact projection reduction, global annihilator proofs, P31 finite-rank section, and independent boundary audits remain valid.',
                next_required_change='Control rational-function coefficient growth or coordinate a larger exclusive-memory run; any modular candidate still requires exact QQ cofactor verification or the homogeneous good-prime Hilbert bridge.',
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),reduction/'report.json',reduction/'secant-exact.sing',reduction/'wronski-exact.sing']},
                evidence_sha256={p.name:sha(p)for p in output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('SAVED',len(attempts),'bounded attempts; no rank/unit theorem')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    run(parser.parse_args().output)

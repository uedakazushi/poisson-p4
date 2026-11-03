#!/usr/bin/env python3
"""Replay cofactors and replace the weak partial-GB bound by a 13-term certificate.

Eleven disjoint singleton supports and two further disjoint supports give
height13 directly.  A covering set of13 variables proves the full329-leading-
term ideal has exactly that height.  No additional GB of the74-generator ideal
is performed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import resource
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    previous=json.loads((args.run/'report.json').read_text());stage=previous['stages'][0]
    assert stage['degree_limit']==3 and stage['basis_count']==329 and stage['cofactors_export_completed']
    source=args.run/'stage-3.witness.m2';summary=args.run/'stage-3.summary'
    assert sha(source)==stage['witness_sha256'] and sha(summary)==stage['summary_sha256']
    line=next(line for line in summary.read_text().splitlines() if line.startswith('LEADING_MONOMIALS '))
    monomials=[mon.strip() for mon in line[len('LEADING_MONOMIALS '):].strip()[1:-1].split(',')]
    def support(mon):return frozenset(0 if name=='u' else int(name[1:])+1 for name in re.findall(r'\bu\b|\bz\d+\b',mon))
    supports=[support(mon) for mon in monomials];assert len(supports)==329 and all(supports)
    forced=sorted({next(iter(s)) for s in supports if len(s)==1});assert len(forced)==11
    selected=[next(i for i,s in enumerate(supports) if s=={variable}) for variable in forced]
    first_pair=next(i for i,s in enumerate(supports) if s=={4,9})
    second_pair=next(i for i,s in enumerate(supports) if s=={19,34})
    selected.extend([first_pair,second_pair]);assert len(selected)==13
    assert all(not(supports[i]&supports[j]) for k,i in enumerate(selected) for j in selected[k+1:])
    covering=set(forced)|{9,34};assert len(covering)==13 and all(s&covering for s in supports)
    # Extract a tiny explicit subcertificate from the already saved full
    # cofactor matrix.  Its13 leading terms have pairwise disjoint supports.
    text=source.read_text().splitlines();assert text[0].startswith('B=matrix{{') and text[0].endswith('}};')
    polynomials=text[0][len('B=matrix{{'):-len('}};')].split(',');assert len(polynomials)==329
    selected_text='B=matrix{{'+','.join(polynomials[i] for i in selected)+'}};\nD=mutableMatrix(R,74,13);\n'
    newindex={index:i for i,index in enumerate(selected)}
    for entry in text[2:]:
        match=re.fullmatch(r'D_\((\d+),(\d+)\)=(.*);',entry)
        if match:
            i,j,polynomial=match.groups();j=int(j)
            if j in newindex:selected_text+=f'D_({i},{newindex[j]})={polynomial};\n'
    selected_text+='C=map(source A,source B,matrix D); assert(A*C==B);\n'
    selected_text+='assert(dim(R/monomialIdeal leadTerm B)==22);\n'
    selectedfile=args.output/'thirteen-cofactors.m2';selectedfile.write_text(selected_text)
    verifier=args.output/'verify.m2'
    verifier.write_text((args.run/'input.m2').read_text()+source.read_text()+
                        'assert(numgens source B==329);\n'+selected_text+
                        'print "ALL329_AND_SELECTED13_COFACTOR_IDENTITIES_PASS"; stdio << flush;\nexit 0;\n')
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
        resource.setrlimit(resource.RLIMIT_CPU,(60,61))
    with (args.output/'replay.log').open('w') as stream:
        result=subprocess.run(['M2','--script',str(verifier)],stdout=stream,stderr=subprocess.STDOUT,timeout=60,preexec_fn=limits)
    log=(args.output/'replay.log').read_text();assert result.returncode==0 and 'ALL329_AND_SELECTED13_COFACTOR_IDENTITIES_PASS' in log
    names=['u']+[f'z{i}' for i in range(34)]
    report=dict(id='Q_THROUGH3_PARTIAL_STRUCTURE_20261103_V63',all_passed=True,coefficient_field='GF(32003)',
                all329_cofactor_identities_independently_replayed=True,
                selected13_cofactor_identities_independently_replayed=True,
                selected_original_basis_column_indices=selected,
                selected_leading_monomials=[monomials[i] for i in selected],
                selected_leading_supports=[[names[j] for j in sorted(supports[i])] for i in selected],
                thirteen_pairwise_disjoint_nonempty_supports=True,
                selected_leading_monomial_ideal_complete_intersection_height=13,
                covering_prime_variables=[names[i] for i in sorted(covering)],
                covering_prime_height=13,
                full329_leading_ideal_dimension=22,
                special74_ideal_dimension_upper_bound=22,
                further_GB_of74_ideal_performed=False,
                QQ_upper_bound_transfer_requires_proper_homogeneous_family=True,
                local_dimension12_or_full_initial_ideal_claimed=False,
                input_sha256={str(p.resolve().relative_to(ROOT)):sha(p) for p in
                              [Path(__file__),args.run/'report.json',source,summary,args.run/'input.m2']},
                evidence_sha256=dict(thirteen_cofactors_m2=sha(selectedfile),verify_m2=sha(verifier),replay_log=sha(args.output/'replay.log')),
                seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['input_sha256','evidence_sha256']},indent=2))


if __name__=='__main__':main()

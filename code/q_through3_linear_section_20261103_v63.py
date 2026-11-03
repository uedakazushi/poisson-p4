#!/usr/bin/env python3
"""One bounded23-variable section of the true74-generator homogeneous ideal."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import random
import resource
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
P=32003
GIB=1024**3


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def expression(terms,names):
    result=[]
    for mon,text in terms:
        c=Fraction(text);assert c.denominator%P
        coefficient=c.numerator*pow(c.denominator,-1,P)%P
        if coefficient:result.append(str(coefficient)+'*'+'*'.join(names[i] for i in mon))
    return '+'.join(result) if result else '0'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--timeout',type=int,default=120)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    minpath=ROOT/'reports/Q_ORDINARY_CUBIC_MINIMAL_GENERATORS_20261103_v63/generators.json.gz'
    quadpath=ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz'
    data=json.load(gzip.open(minpath,'rt'));quad=json.load(gzip.open(quadpath,'rt'))
    oldnames=data['variables'];assert len(oldnames)==35
    names=[f'y{i}' for i in range(23)]
    rng=random.Random(20261103);T=[[rng.choice([-2,-1,1,2]) for j in range(23)] for i in range(12)]
    mapping=['+'.join(f'({coefficient})*{names[j]}' for j,coefficient in enumerate(row)) for row in T]+names
    section=dict(seed=20261103,pivot_matrix='12x12 identity',pivot_determinant_Z=1,
                 original_variables=oldnames,eliminated_variables=oldnames[:12],retained_variables=oldnames[12:],
                 section_variables=names,integer_matrix_T=T,
                 linear_form_convention='x_i-sum_j T_ij*x_(12+j), i=0,...,11',
                 quotient_ring_isomorphism_over_Z=True,genericity_not_assumed=True)
    sectionpath=args.output/'section.json';sectionpath.write_text(json.dumps(section,indent=2)+'\n')
    generators=[]
    for index in data['quadratic_generators']:
        generators.append(expression([([i,j],c) for i,j,c in quad['all350_quadratic_rows'][index]],oldnames))
    for row in data['cubic_generators']:generators.append(expression([(mon,c) for mon,c in row['terms']],oldnames))
    assert len(generators)==74
    header=f'O=GF({P})['+','.join(oldnames)+',MonomialOrder=>GRevLex];\n'
    header+='A35=matrix{{'+',\n'.join(generators)+'}};\n'
    header+=f'R=GF({P})['+','.join(names)+',MonomialOrder=>GRevLex];\n'
    header+='phi=map(R,O,{'+','.join(mapping)+'});\nA=phi A35;\n'
    (args.output/'input.m2').write_text(header)
    text=header+'gbTrace=1; proceed=true;\n'
    for degree in [3,4,5,6]:
        suffix=str(degree);prefix=(args.output/f'stage-{degree}').resolve()
        text+='if proceed then (\n'
        text+=f'G{suffix}=gb(A,ChangeMatrix=>true,DegreeLimit=>{degree},HardDegreeLimit=>6,CodimensionLimit=>23,BasisElementLimit=>2500);\n'
        text+=f'B{suffix}=gens G{suffix}; C{suffix}=getChangeMatrix G{suffix}; assert(A*C{suffix}==B{suffix});\n'
        text+=f'L{suffix}=monomialIdeal leadTerm B{suffix}; bound{suffix}=dim(R/L{suffix}); n{suffix}=numgens source B{suffix};\n'
        text+=f'f=openOut "{prefix}.summary"; f << "DEGREE_LIMIT {degree}" << endl << "BASIS_COUNT " << n{suffix} << endl;\n'
        text+=f'f << "SECTION_DIMENSION_BOUND " << bound{suffix} << endl << "ALL74_SECTION_COFACTOR_IDENTITIES_VERIFIED true" << endl;\n'
        text+=f'f << "LEADING_MONOMIALS " << toString flatten entries leadTerm B{suffix} << endl; close f;\n'
        text+=f'w=openOut "{prefix}.witness.m2"; w << "B=matrix{{{{";\n'
        text+=f'scan(n{suffix},j->(if j>0 then w << ","; w << toString B{suffix}_(0,j)));\n'
        text+=f'w << "}}}};" << endl << "D=mutableMatrix(R,74," << n{suffix} << ");" << endl;\n'
        text+=f'scan(n{suffix},j->scan(74,i->if C{suffix}_(i,j)!=0 then w << "D_(" << i << "," << j << ")=" << toString C{suffix}_(i,j) << ";" << endl));\n'
        text+='w << "C=map(source A,source B,matrix D); assert(A*C==B);" << endl;\n'
        text+=f'w << "assert(dim(R/monomialIdeal leadTerm B)==" << bound{suffix} << ");" << endl; close w;\n'
        text+=f'f=openOut "{prefix}.done"; f << "COFACTORS_EXPORTED true" << endl; close f;\n'
        text+=f'proceed=(bound{suffix}>0 and n{suffix}<2500);\n);\n'
    text+='exit 0;\n';script=args.output/'compute.m2';script.write_text(text)
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(3*GIB,3*GIB))
        resource.setrlimit(resource.RLIMIT_CPU,(args.timeout,args.timeout+1))
    status='completed';maximum_rss=0;maximum_vmsize=0;observations=[]
    with (args.output/'m2.log').open('w') as log:
        process=subprocess.Popen(['M2','--script',str(script)],stdout=log,stderr=subprocess.STDOUT,preexec_fn=limits)
        beginning=time.monotonic()
        while process.poll() is None:
            elapsed=time.monotonic()-beginning
            try:
                state=Path(f'/proc/{process.pid}/status').read_text().splitlines()
                values={line.split(':')[0]:int(line.split()[1]) for line in state if line.startswith(('VmRSS:','VmSize:'))}
                maximum_rss=max(maximum_rss,values.get('VmRSS',0));maximum_vmsize=max(maximum_vmsize,values.get('VmSize',0))
                if not observations or elapsed-observations[-1]['seconds']>=5:observations.append(dict(seconds=elapsed,**values))
                if values.get('VmRSS',0)*1024>int(2.75*GIB):status='memory_growth_stop';process.terminate();break
            except OSError:pass
            if elapsed>args.timeout:status='timeout';process.terminate();break
            time.sleep(0.25)
        try:process.wait(timeout=3)
        except subprocess.TimeoutExpired:process.kill();process.wait()
        if status=='completed' and process.returncode:status=f'exit_{process.returncode}'
    stages=[]
    for degree in [3,4,5,6]:
        summary=args.output/f'stage-{degree}.summary';done=args.output/f'stage-{degree}.done';witness=args.output/f'stage-{degree}.witness.m2'
        if not summary.exists():continue
        lines=summary.read_text().splitlines()
        if not any(line.startswith('SECTION_DIMENSION_BOUND ') for line in lines):continue
        dimension=int(next(line.split()[1] for line in lines if line.startswith('SECTION_DIMENSION_BOUND ')))
        stage=dict(degree_limit=degree,basis_count=int(next(line.split()[1] for line in lines if line.startswith('BASIS_COUNT '))),
                   section_dimension_upper_bound=dimension,original35_dimension_upper_bound=dimension+12,
                   all74_section_cofactor_identities_verified='ALL74_SECTION_COFACTOR_IDENTITIES_VERIFIED true' in lines,
                   cofactors_export_completed=done.exists(),summary_sha256=sha(summary))
        if witness.exists():stage['witness_sha256']=sha(witness)
        stages.append(stage)
    report=dict(id='Q_THROUGH3_LINEAR_SECTION_20261103_V63',status=status,
                new_generators=dict(quadratics=57,true_cubics=17,total=74),
                original_variables=35,linear_forms=12,section_variables=23,section_pivot_determinant_Z=1,
                deterministic_integer_section=True,genericity_assumption_used=False,
                same_cached_GB_degree_checkpoints=[3,4,5,6],hard_degree_limit=6,basis_limit=2500,
                wall_timeout_seconds=args.timeout,memory_limit_bytes=3*GIB,
                maximum_observed_rss_kb=maximum_rss,maximum_observed_vmsize_kb=maximum_vmsize,
                stages=stages,
                Artinian_section_certified=any(s['section_dimension_upper_bound']==0 and s['all74_section_cofactor_identities_verified'] for s in stages),
                QQ_dimension_transfer_requires_explicit_proper_homogeneous_family=True,
                full_initial_ideal_or_branches_claimed=False,
                input_sha256={str(p.resolve().relative_to(ROOT)):sha(p) for p in [Path(__file__),minpath,quadpath]},
                section_sha256=sha(sectionpath),script_sha256=sha(script),seconds=time.monotonic()-start)
    (args.output/'memory-monitor.json').write_text(json.dumps(observations,indent=2)+'\n')
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='input_sha256'},indent=2))


if __name__=='__main__':main()

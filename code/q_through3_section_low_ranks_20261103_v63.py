#!/usr/bin/env python3
"""Small tensor-projection/rank certificate for the failed fixed linear section."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import resource
import shlex
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
P=32003


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def fp(value):
    value=Fraction(value);assert value.denominator%P
    return value.numerator*pow(value.denominator,-1,P)%P


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--section',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    sectionpath=args.section/'section.json';section=json.loads(sectionpath.read_text());T=section['integer_matrix_T']
    minpath=ROOT/'reports/Q_ORDINARY_CUBIC_MINIMAL_GENERATORS_20261103_v63/generators.json.gz'
    quadpath=ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz'
    data=json.load(gzip.open(minpath,'rt'));quad=json.load(gzip.open(quadpath,'rt'))
    original=[[(tuple([i,j]),fp(c)) for i,j,c in quad['all350_quadratic_rows'][source]] for source in data['quadratic_generators']]
    original += [[(tuple(mon),fp(c)) for mon,c in row['terms']] for row in data['cubic_generators']]
    assert len(original)==74
    lines=['57 17']+[' '.join(map(str,row)) for row in T]
    for row in original:lines.append(str(len(row))+' '+' '.join(' '.join(map(str,(*mon,c))) for mon,c in row))
    source=args.output/'projection.in';source.write_text('\n'.join(lines)+'\n')
    native=ROOT/'code/q_through3_section_low_ranks_20261103_v63.cpp';executable=args.output/'rank'
    flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],capture_output=True,text=True,check=True).stdout)
    subprocess.run(['c++','-O3','-std=c++17',str(native),*flags,'-o',str(executable)],check=True)
    evidence=args.output/'projected-polynomials-and-ranks.txt'
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3));resource.setrlimit(resource.RLIMIT_CPU,(60,61))
    subprocess.run([str(executable.resolve()),str(source.resolve()),str(evidence.resolve())],check=True,timeout=60,preexec_fn=limits)
    lines=iter(evidence.read_text().splitlines());assert next(lines)=='POLYNOMIALS';projected=[]
    for expected in range(74):
        index,degree,count=map(int,next(lines).split());assert index==expected
        row=[]
        for _ in range(count):
            coefficient,*mon=map(int,next(lines).split());assert len(mon)==degree and coefficient%P
            row.append((tuple(mon),coefficient))
        projected.append(row)
    tokens=next(lines).split();assert tokens[0]=='RANKS' and tokens[4]=='SECONDS'
    rQ,rW,rAll=map(int,tokens[1:4])
    # Independent scalar checks compare projected coefficient tables with the
    # original exact recipe at three dense points, no tensor expansion reused.
    for seed in [2,11,29]:
        y=[(j*j+seed*j+2*seed+1)%P for j in range(23)]
        x=[sum(c*v for c,v in zip(row,y))%P for row in T]+y
        for i,(before,after) in enumerate(zip(original,projected)):
            def value(poly,point):
                total=0
                for mon,c in poly:
                    term=c
                    for j in mon:term=term*point[j]%P
                    total=(total+term)%P
                return total
            assert value(before,x)==value(after,y),('projection scalar check',seed,i)
    report=dict(id='Q_THROUGH3_SECTION_LOW_RANKS_20261103_V63',all_passed=True,coefficient_field='GF(32003)',
                fixed_integer_section_seed=section['seed'],section_variables=23,
                projected_quadratic_rank=rQ,linear_quadratic_multiple_rank=rW,
                full_degree3_ideal_rank=rAll,new_cubic_class_count=rAll-rW,
                exact_homogeneous_section_Hilbert_values_mod_p_through3=[1,23,276-rQ,2300-rAll],
                Artinian_length_if_finite_at_least=1+23+276-rQ+2300-rAll,
                cached_Sym2_Sym3_tensor_projection=True,all74_projections_match_three_independent_dense_scalar_checks=True,
                exact_full_section_dimension_or_QQ_Hilbert_values_not_claimed=True,
                no_GB_or_higher_jet_performed=True,matrix_memory_estimate_bytes=(57*23*2300+(57*23+17)*2300)*8,
                native_rank_seconds=float(tokens[5]),
                input_sha256={str(p.resolve().relative_to(ROOT)):sha(p) for p in [Path(__file__),native,sectionpath,minpath,quadpath]},
                evidence_sha256=sha(evidence),seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='input_sha256'},indent=2))


if __name__=='__main__':main()

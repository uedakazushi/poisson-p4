#!/usr/bin/env python3
"""Direct full135 central support and prime-unit restoration at [0:0:1]."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize

ROOT=Path(__file__).resolve().parents[1]
PRIME=32003


def verify(output):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    prior_path=ROOT/'reports/S61_MARKED_SIGMA_PENCIL_INFINITY_20261102.json'
    prior=json.loads(prior_path.read_text());assert prior['all_passed']
    for filename,expected in prior['input_sha256'].items():assert hashlib.sha256((ROOT/filename).read_bytes()).hexdigest()==expected,filename
    model=S61();one=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==1]
    chart=[j for j in range(45)if j!=25];cuts=[0,9,27,32];columns=[j for j in chart if j not in cuts]
    marked_path=ROOT/'reports/S61_MARKED_SIGMA_RESIDUAL_FINITE_20261102.json'
    old=json.loads(marked_path.read_text());qw=old['frame']['q_weights']
    rho=s.Symbol('rho');beta=dict.fromkeys(one,s.S.Zero);beta.update({4:rho,11:1,26:1})
    _,forms=specialize(model,beta,True)
    variables=[f'q{j}'for j in columns]+['t','rho'];weights=[qw[chart.index(j)]for j in columns]+[3,1]
    assert len(variables)==42 and min(weights)>0
    polynomials=[];count=0;grades=[]
    for form in forms:
        terms=[];rowgrades=set()
        for mon,coefficient in form.items():
            if any(j in cuts for j in mon):continue
            coefficient*=(-1)**mon.count(25)
            factors=[f'q{j}'for j in mon if j!=25]+['t']*(2-len(mon))
            for (power,),c in s.Poly(coefficient,rho).terms():
                if not c:continue
                assert int(c.q)%PRIME
                rowgrades.add(power+sum(weights[variables.index(x)]for x in factors))
                value=int(c.p)*pow(int(c.q),-1,PRIME)%PRIME
                if value:
                    terms.append(str(value)+('*rho^'+str(power)if power else'')+''.join('*'+x for x in factors));count+=1
        assert len(rowgrades)<=1
        grades.append(next(iter(rowgrades),-1));polynomials.append('+'.join(terms)or'0')
    text='R=(ZZ/32003)['+','.join(variables)+',Degrees=>{'+','.join(map(str,weights))+'}];\n'
    text+='I=ideal('+','.join(polynomials)+');G=gb I;A=ideal('+','.join(variables[:-1])+');\n'
    text+='assert(isSubset(I,A));assert all(flatten entries gens A,f->any(1..20,n->f^n%G==0));\n'
    text+='U=I+ideal(rho);assert(dim(R/U)==0);assert(numColumns basis(R/U)==79);\n'
    text+='P=I+ideal(rho-1);assert(dim(R/P)==0);assert(numColumns basis(R/P)==4);assert(saturate(P,A)==ideal(1_R));\n'
    text+='print("HS="|toString hilbertSeries(R/I,Reduce=>true));print("T="|toString(t%G));print("T2="|toString(t^2%G));\n'
    text+='assert(t%G!=0 and t^2%G==0);print "INFINITY PRIME SUPPORT PASS";exit 0;\n'
    script=output/'infinity_modp.m2';script.write_text(text);log=output/'infinity_modp.log';m2=shutil.which('M2');assert m2
    with log.open('w')as stream:
        result=subprocess.run([m2,'--script',str(script.resolve())],stdout=stream,stderr=subprocess.STDOUT,timeout=150,check=True)
    content=log.read_text();assert'INFINITY PRIME SUPPORT PASS'in content
    line=next(line[3:]for line in content.splitlines()if line.startswith('HS='));T=s.Symbol('T')
    actual=s.sympify(line.replace('^','**'),locals={'T':T})
    expected=(1+7*T+23*T*T+35*T**3-24*T**4-36*T**5-2*T**6)/(1-T)
    assert s.cancel(actual-expected)==0
    # Restore all three group cuts and the H15 normal; q32 stays a base coordinate.
    flag_path=ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json';flag=json.loads(flag_path.read_text())['flag']
    act_path=ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json';act=s.Matrix(json.loads(act_path.read_text())['primitive_action']).applyfunc(s.sympify)
    us=sorted(act.free_symbols,key=lambda x:int(x.name[1:]));at=act.subs(dict(zip(us,[s.sympify(beta[i]).subs(rho,0)for i in one])))
    orbit=s.Matrix(flag['q_orbit_differential']).applyfunc(s.Rational)
    block=orbit[[chart.index(j)for j in [0,9,27]],:].col_join(at[0,:]);det=block[:,[4,5,6,9]].det();assert det==-2
    standard_path=ROOT/'reports/S61_MARKED_SIGMA_STANDARD_FLAG_SPREAD_20261103_v62.json'
    standard=json.loads(standard_path.read_text());assert standard['all_passed']and standard['all_other_beta_parameters_retained']and standard['ambient_t_zero']
    paths=[Path(__file__),prior_path,marked_path,flag_path,act_path,standard_path,
           ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
           ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',
           ROOT/'research/notes/S61_MARKED_SIGMA_PENCIL_INFINITY_PRIME_20261103_v62.md']
    report=dict(id='S61_MARKED_SIGMA_PENCIL_INFINITY_PRIME_20261103_v62',all_passed=True,prime=PRIME,
                core=[0,0,1],raw_variable_count=42,raw135_rows=len(forms),nonzero_terms=count,
                row_weights=grades,all_positive_variable_weights=weights,
                exact_reduced_special_support='standard rho axis',rho_zero_length=79,rho_one_length=4,
                radical_verified_by_two_containments_and_powers_bound=20,
                full_hilbert_series_modp=line,group_restoration_minor=str(det),
                standard_axis_nonzero_mixed_characteristic_t_zero=True,
                actual_saturated_support_to_origin_weight_lemma=True,
                full_ambient_mixed_characteristic_homogenizer_nilpotent=True,
                all_nuisance_at_this_core_removed=True,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
                evidence_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in [script,log]},
                m2_returncode=result.returncode,seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=verify(args.output);print('PASS',result['id'],result['seconds'])

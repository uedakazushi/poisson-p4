#!/usr/bin/env python3
"""Rebuild all135 direct infinity polynomials independently of specialize()."""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path

from s61_model import S61

ROOT=Path(__file__).resolve().parents[1]
P=32003


def verify(output):
    folder=ROOT/'reports/S61_MARKED_SIGMA_PENCIL_INFINITY_PRIME_20261103_v62'
    source=folder/'report.json'
    report=json.loads(source.read_text())
    model=S61()
    one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    qpos=dict(zip(two,range(45)))
    cuts={0,9,27,32}
    retained=[i for i in range(45) if i!=25 and i not in cuts]
    names=[f'q{i}'for i in retained]+['t','rho']
    location=dict(zip(names,range(42)))
    weights=report['all_positive_variable_weights']
    assert len(names)==42 and min(weights)>0
    equations=model.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
    expected=[]
    for row in equations:
        terms=defaultdict(int)
        for mon,n,d in row:
            assert d%P
            coefficient=n*pow(d,-1,P)%P
            powers=[0]*42
            for i in mon:
                if i in qpos:
                    q=qpos[i]
                    if q==25:coefficient=-coefficient
                    elif q in cuts:coefficient=0;break
                    else:powers[location[f'q{q}']]+=1
                elif i==4:powers[location['rho']]+=1
                elif i in {11,26}:pass
                else:coefficient=0;break
            if not coefficient:continue
            powers[location['t']]+=2-sum(i in qpos for i in mon)
            terms[tuple(powers)]=(terms[tuple(powers)]+coefficient)%P
        expected.append({m:c for m,c in terms.items() if c})
    script=(folder/'infinity_modp.m2').read_text()
    raw=script.split('I=ideal(',1)[1].split(');G=gb I;',1)[0].split(',')
    assert len(raw)==135
    for index,expression in enumerate(raw):
        obtained=defaultdict(int)
        if expression!='0':
            for term in expression.split('+'):
                factors=term.split('*')
                coefficient=int(factors[0])
                mon=[0]*42
                for factor in factors[1:]:
                    pieces=factor.split('^')
                    mon[location[pieces[0]]]+=int(pieces[1]) if len(pieces)>1 else 1
                obtained[tuple(mon)]=(obtained[tuple(mon)]+coefficient)%P
        obtained={m:c for m,c in obtained.items() if c}
        assert obtained==expected[index],index
        grades={sum(a*b for a,b in zip(m,weights)) for m in obtained}
        assert len(grades)<=1
    log=(folder/'infinity_modp.log').read_text()
    assert 'INFINITY PRIME SUPPORT PASS' in log and 'T2=0' in log
    assert 'assert(isSubset(I,A))'in script and 'f^n%G==0'in script
    assert 'numColumns basis(R/U)==79'in script and 'numColumns basis(R/P)==4'in script
    assert 'saturate(P,A)==ideal(1_R)'in script
    assert report['group_restoration_minor']=='-2'
    for path,digest in report['input_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    for path,digest in report['evidence_sha256'].items():
        assert hashlib.sha256((folder/path).read_bytes()).hexdigest()==digest
    paths=[Path(__file__),source,folder/'infinity_modp.m2',folder/'infinity_modp.log',
           ROOT/'research/notes/S61_MARKED_SIGMA_PENCIL_INFINITY_PRIME_20261103_v62.md']
    result=dict(id='S61_MARKED_SIGMA_PENCIL_INFINITY_PRIME_20261103_v62_REVIEW',all_passed=True,
                direct_original135_reconstruction_matches_every_modular_coefficient=True,
                positive_variable_count=42,nonpositive_q32_retained_in_base_not_group_cut=True,
                coefficient_denominators_prime_units=True,positive_homogeneity_checked=True,
                saved_M2_full_support_two_inclusions_and_lengths_accepted=True,
                finite_field_Groebner_engine_independently_rerun=False,
                actual_V_t_saturated_quotient_finite_weight_argument_accepted=True,
                saturation_and_specialization_commutation_used=False,
                input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()for path in paths})
    with output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print('PASS DIRECT ORIGINAL135 ALL INFINITY COEFFICIENTS')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();verify(args.output)

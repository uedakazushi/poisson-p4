#!/usr/bin/env python3
"""Retain G6 in the exact a-only ten-variable consequence model.

Replays all33 recorded substitutions against fresh original135 coefficients.
The old partial model kept only rows of weight<=4 and omitted G6; this script
checks every frozen original134 coefficient and adds the exact missing G6.
No completeness of the partial consequence ideal is claimed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile
import time

import sympy as s
from sympy.polys.rings import ring

from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize

ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=ROOT/'research/handoff/20261102/prior_binary_probes/HYPERPLANE_PARTIAL_POLYNOMIAL_ELIMINATION_20261102.tar.gz'
MEMBER='work/generic_hyperplane_polynomial_a_elimination_20261102.json'
EXPECTED='e28e9336d7cee80a5a9150ba7e6aa06fae77fcadd3a443aa4c2dc8f23ff468aa'


def run(output):
    started=time.time()
    output=Path(output).resolve()
    output.mkdir(parents=True,exist_ok=False)
    with tarfile.open(ARCHIVE,'r:gz') as archive:
        raw=archive.extractfile(MEMBER).read()
    assert hashlib.sha256(raw).hexdigest()==EXPECTED
    recorded=json.loads(raw)
    chartpath=ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json'
    oldpath=ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json'
    chart=json.loads(chartpath.read_text())
    old=json.loads(oldpath.read_text())['linear']
    R,a,e,d=ring('a,e,d',s.QQ)
    beta={int(k[1:])-1:s.sympify(v) for k,v in chart['beta'].items()}
    _,source=specialize(S61(),beta,True)
    assert len(source)==135
    weights=chart['local_coordinate_weights']
    rowweights=chart['equation_weights']+[1]
    fixed,columns=old['fixed'],old['cols']
    cuts={19,23}
    equations=[{} for _ in range(136)]
    def add(poly,mon,coefficient):
        if coefficient:
            value=poly.get(mon,R.zero)+coefficient
            if value:
                poly[mon]=value
            else:
                poly.pop(mon,None)
    for row,equation in enumerate(source):
        if rowweights[row]>4:
            continue
        for mon,value in equation.items():
            coefficient=R.from_expr(s.sympify(value))
            if len(mon)==2:
                i,j=mon
                if i!=fixed and j==fixed and columns.index(i) not in cuts:
                    add(equations[row],(columns.index(i),),coefficient)
                if j!=fixed and i==fixed and columns.index(j) not in cuts:
                    add(equations[row],(columns.index(j),),coefficient)
            elif mon==(fixed,):
                add(equations[row],(44,),coefficient)
            if fixed not in mon:
                local=tuple(sorted([columns.index(i) for i in mon]+[44]*(2-len(mon))))
                if not set(local).intersection(cuts):
                    add(equations[row],local,coefficient)
    equations[135]={(5,):R.one}
    target={(44,):R.one}
    target_denominator=0
    g6_history=[]
    def multiply(left,right):
        result={}
        for mon,c in left.items():
            for other,v in right.items():
                add(result,tuple(sorted(mon+other)),c*v)
        return result
    def substitute(poly,variable,replacement,denominator_power):
        maximum=max((mon.count(variable) for mon in poly),default=0)
        if not maximum:
            return poly.copy(),0
        powers=[{():R.one},replacement]
        for _ in range(2,maximum+1):
            powers.append(multiply(powers[-1],replacement))
        result={}
        apowers=[a**(denominator_power*i) for i in range(maximum+1)]
        for mon,c in poly.items():
            count=mon.count(variable)
            rest=tuple(i for i in mon if i!=variable)
            for other,value in powers[count].items():
                add(result,tuple(sorted(rest+other)),c*value*apowers[maximum-count])
        return result,maximum*denominator_power
    assert len(recorded['steps'])==len(recorded['normalization_history'])==33
    for step,(record,history) in enumerate(zip(recorded['steps'],recorded['normalization_history'])):
        row,variable=record['equation'],record['variable']
        equation=equations[row]
        power=record['a_power']
        scalar=s.QQ(record['scalar'])
        assert equation[(variable,)]==scalar*a**power
        assert all(variable not in mon for mon in equation if mon!=(variable,))
        replacement={mon:-value/scalar for mon,value in equation.items() if mon!=(variable,)}
        histories={item['row']:item for item in history}
        updated=[]
        for index,poly in enumerate(equations):
            if index==row:
                updated.append({})
                continue
            value,cleared=substitute(poly,variable,replacement,power)
            if index==134:
                content=min((mon[0] for c in value.values() for mon in c),default=0)
                value={mon:c.exquo(a**content) for mon,c in value.items()}
                g6_history.append(dict(step=step,cleared_a=cleared,removed_a_content=content))
            else:
                previous=histories.get(index,{'cleared_a':0,'unit_powers':[0,0,0,0]})
                assert cleared==previous['cleared_a']
                assert not any(previous['unit_powers'][1:])
                content=previous['unit_powers'][0]
                if content:
                    value={mon:c.exquo(a**content) for mon,c in value.items()}
            updated.append(value)
        equations=updated
        target,cleared=substitute(target,variable,replacement,power)
        target_denominator+=cleared
        for index,poly in enumerate(equations):
            assert all(sum(weights[i] for i in mon)==rowweights[index] for mon in poly)
        if (step+1)%8==0:
            print('replayed',step+1,'G6 normal monomials',len(equations[134]),'seconds',time.time()-started,flush=True)
    def decode(poly):
        return {tuple(mon):R.from_dict({tuple(mm):s.QQ(c) for mm,c in coeffs}) for mon,coeffs in poly}
    saved={item['original_row']:decode(item['polynomial']) for item in recorded['equations_sparse']}
    assert all(equations[i]==saved.get(i,{}) for i in range(134))
    assert not equations[135]
    assert target==decode(recorded['targets_sparse'][-1])
    assert target_denominator==recorded['target_a_denominators'][-1]==9
    assert equations[134] and rowweights[134]==3
    remaining=recorded['remaining_variables']
    assert all(set(mon).issubset(remaining) for poly in equations for mon in poly)
    encoded=lambda poly:[[list(mon),[[list(mm),str(c)] for mm,c in coeff.items()]] for mon,coeff in sorted(poly.items())]
    consequence=dict(G6=encoded(equations[134]),weight=3,t_denominator_a_power=9,
                     t_numerator=encoded(target),remaining_variables=remaining,
                     remaining_weights=[weights[i] for i in remaining])
    datafile=output/'additional-G6.json'
    datafile.write_text(json.dumps(consequence,separators=(',',':'))+'\n')
    files=[Path(__file__),ARCHIVE,chartpath,oldpath,
           ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
           ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    report=dict(id='S61_BINARY_NONSQUARE_PARTIAL135_20261103_v62',all_passed=True,
                field='QQ[a,e,d,a^-1]',recorded_substitutions_replayed=33,
                all_frozen_original134_coefficients_matched=True,exact_t_numerator_matched=True,
                original_G6_retained=True,G6_weight=3,G6_normal_monomials=len(equations[134]),
                G6_core_terms=sum(len(c) for c in equations[134].values()),
                G6_a_normalization_history=g6_history,
                frozen_model_archive_member=MEMBER,frozen_model_member_sha256=EXPECTED,
                limitations=['Only original rows of weight<=4 are represented; this is a consequence model.',
                             'No ideal equality with all135 equations or ambient saturation is asserted.',
                             'No generic or all-core t-nilpotence follows from replay alone.'],
                input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in files},
                generated_sha256={datafile.name:hashlib.sha256(datafile.read_bytes()).hexdigest()},
                seconds=time.time()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'G6_normal_monomials':report['G6_normal_monomials'],'G6_core_terms':report['G6_core_terms'],'seconds':report['seconds']}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)

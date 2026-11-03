#!/usr/bin/env python3
"""All-core coordinate cover for the exact quadratic product pair."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]


def run(output):
    started=time.time()
    output=Path(output).resolve()
    output.mkdir(parents=True,exist_ok=False)
    quadpath=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json'
    pairpath=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_PAIR_20261102.json'
    quad=json.loads(quadpath.read_text())
    pair=json.loads(pairpath.read_text())
    for report in [quad,pair]:
        assert report['all_passed']
        for path,expected in report['input_sha256'].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==expected
    R,a,e,d=ring('a,e,d',s.QQ)
    K=R.to_field()
    known=[a,500*a+90*e**2+3*e-24,
           (50*a+18*e**2-9*e-2)*d+60*e-40,
           100*a+18*e**2-9*e-2]
    def unit(poly):
        if not poly:
            return False
        for factor in known:
            while poly.rem(factor)==0:
                poly=poly.exquo(factor)
        return all(not any(mon) for mon in poly)
    P,*v=ring('v0,v1,v2,v3,v4,v5,v6,v7',R)
    qs=[P.from_dict({tuple(mon):R.from_dict({tuple(mm):s.QQ(c) for mm,c in pp})
                    for mon,pp in row}) for row in quad['quadratic_rows_sparse']]
    Z=3*v[0]-v[1]+v[6]
    W=v[3]
    combo=-(d*e+8)*qs[1]+8*qs[2]+d*qs[3]
    ell=combo.exquo(W)
    assert (-10*a-3*e**2+3*e)*qs[1]+3*(e-1)*qs[3]+10*a*qs[4]==-a*Z*ell
    assert s.expand(ell.as_expr()-s.sympify(pair['linear_form']))==0
    coeff=[]
    for i in range(8):
        mon=tuple(int(i==j) for j in range(8))
        coeff.append(ell.get(mon,R.zero))
    conditions=[coeff[2],coeff[4],coeff[7],coeff[0]+3*coeff[1],coeff[6]+coeff[1]]
    script=output/'ell-cover.m2'
    script.write_text('\n'.join([
        'R=QQ[a,e,d,z];F=500*a+90*e^2+3*e-24;B=(50*a+18*e^2-9*e-2)*d+60*e-40;C=100*a+18*e^2-9*e-2;',
        'I=ideal('+','.join(str(q).replace('**','^') for q in conditions)+',z*a*F*B*C-1);',
        'G=gb I;assert(1%G==0);L=matrix{{1_R}} // gens I;',
        'assert(gens I*L==matrix{{1_R}});print "EXACT_ELL_COVER_AND_BEZOUT_PASS";print toString L;'])+'\n')
    m2=shutil.which('M2')
    assert m2
    cp=subprocess.run([m2,'--script',str(script)],cwd=ROOT,text=True,capture_output=True,timeout=30)
    (output/'ell-cover.log').write_text(cp.stdout+cp.stderr)
    cp.check_returncode()
    assert 'EXACT_ELL_COVER_AND_BEZOUT_PASS' in cp.stdout
    # Quotient by the CONSTANT linear forms Z,W, leaving five light variables.
    oldpositions=[0,1,2,4,7]
    rows=[]
    for index in [0,1,3,5]:
        row={}
        for mon,coefficient in qs[index].items():
            if mon[3]:
                continue
            assert not mon[5]
            k=mon[6]
            for j in range(k+1):
                power=[mon[i] for i in oldpositions]
                power[0]+=j
                power[1]+=k-j
                power=tuple(power)
                row[power]=row.get(power,R.zero)+coefficient*(int(s.binomial(k,j))*(-3)**j)
        rows.append({mon:c for mon,c in row.items() if c})
    monomials=sorted(set().union(*rows))
    matrix=[[K(row.get(mon,R.zero)) for mon in monomials] for row in rows]
    pivots=[]
    used=[]
    for step in range(4):
        candidates=sorted((len(value.numer)+len(value.denom),i,j)
                          for i in range(step,4) for j,value in enumerate(matrix[i])
                          if j not in used and value)
        i,j=next((i,j) for _,i,j in candidates
                 if unit(matrix[i][j].numer) and unit(matrix[i][j].denom))
        matrix[step],matrix[i]=matrix[i],matrix[step]
        value=matrix[step][j]
        pivots.append(dict(row_swap=i,column=j,value=str(value)))
        used.append(j)
        matrix[step]=[x/value for x in matrix[step]]
        for later in range(step+1,4):
            value=matrix[later][j]
            if value:
                matrix[later]=[x-value*y for x,y in zip(matrix[later],matrix[step])]
    paths=[Path(__file__),quadpath,pairpath]
    report=dict(id='S61_BINARY_NONSQUARE_PAIR_COVER_20261103_v62',all_passed=True,field='QQ',
                known_open='a*F*B*C!=0',ell_Z_W_linearly_independent_at_every_core=True,
                ell_dependence_coefficient_conditions=list(map(str,conditions)),
                ell_cover_ideal_is_unit=True,explicit_Bezout_identity_checked=True,
                four_remaining_quadrics_restricted_mod_Z_W_rank=4,
                restriction_variables_zero_based_positions=oldpositions,
                restriction_pivots=pivots,
                limitations=['No assertion of complete intersection at every core.',
                             'The full cubic-product rank and quartic coverage remain separate.',
                             'No new ambient t-nilpotence follows from this coordinate cover alone.'],
                input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
                generated_sha256={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in [script,output/'ell-cover.log']},
                seconds=time.time()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'ell_Z_W_independent':True,'restriction_rank':4,'seconds':report['seconds']}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)

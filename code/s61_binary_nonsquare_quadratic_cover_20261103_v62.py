#!/usr/bin/env python3
"""Exact all-core rank-six cover of the full-S61 weighted quadrics.

The original six rows were independently reconstructed from all135 equations
in the pinned 20261102 certificate.  We retain their exact QQ coefficients,
replay its three known-unit pivots, then compute the ENTIRE residual minor
ideal, rather than only its height-one gcd.
"""
import argparse
import hashlib
from itertools import combinations
import json
from pathlib import Path
import shutil
import subprocess
import time

import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(output):
    started=time.monotonic()
    output=Path(output).resolve()
    output.mkdir(parents=True,exist_ok=False)
    data=json.loads(SOURCE.read_text())
    assert data['all_passed'] and data['quadratic_row_count']==6
    for path,value in data['input_sha256'].items():
        assert sha(ROOT/path)==value,path
    R,a,e,d=ring('a,e,d',s.QQ)
    K=R.to_field()
    known=[a,500*a+90*e**2+3*e-24,
           (50*a+18*e**2-9*e-2)*d+60*e-40,
           100*a+18*e**2-9*e-2]
    def is_unit(poly):
        if not poly:
            return False
        for factor in known:
            while poly.rem(factor)==0:
                poly=poly.exquo(factor)
        return all(not any(mon) for mon in poly)
    rows=[{tuple(mon):R.from_dict({tuple(mm):s.QQ(c) for mm,c in pp})
           for mon,pp in row} for row in data['quadratic_rows_sparse']]
    monomials=sorted(set().union(*rows))
    assert len(monomials)==28
    matrix=[[K(row.get(mon,R.zero)) for mon in monomials] for row in rows]
    pivots=[]
    for row,(swap,column,expected) in enumerate(data['unit_pivots']):
        matrix[row],matrix[swap]=matrix[swap],matrix[row]
        pivot=matrix[row][column]
        assert s.cancel(pivot.as_expr()-s.sympify(expected))==0
        assert is_unit(pivot.numer) and is_unit(pivot.denom)
        pivots.append(column)
        matrix[row]=[z/pivot for z in matrix[row]]
        for later in range(row+1,6):
            value=matrix[later][column]
            if value:
                matrix[later]=[x-value*y for x,y in zip(matrix[later],matrix[row])]
    columns=[j for j in range(28) if j not in pivots and any(matrix[i][j] for i in range(3,6))]
    assert len(columns)==25
    residual=[]
    denominators=[]
    for row in matrix[3:]:
        den=R.one
        for column in columns:
            den=den.lcm(row[column].denom)
        assert is_unit(den)
        vector=[row[column].numer*den.exquo(row[column].denom) for column in columns]
        content=R.zero
        for coefficient in vector:
            content=content.gcd(coefficient)
        assert content and all(not any(mon) for mon in content)
        residual.append([coefficient.exquo(content) for coefficient in vector])
        denominators.append(str(den))
    minors=[]
    selections=[]
    for i,j,k in combinations(range(25),3):
        poly=(residual[0][i]*(residual[1][j]*residual[2][k]-residual[1][k]*residual[2][j])
              -residual[0][j]*(residual[1][i]*residual[2][k]-residual[1][k]*residual[2][i])
              +residual[0][k]*(residual[1][i]*residual[2][j]-residual[1][j]*residual[2][i]))
        if not poly:
            continue
        orders=[]
        for factor in known:
            exponent=0
            while poly.rem(factor)==0:
                poly=poly.exquo(factor)
                exponent+=1
            orders.append(exponent)
        poly=poly.monic()
        if poly not in minors:
            minors.append(poly)
            selections.append(dict(columns=[i,j,k],stripped_known_orders=orders))
    assert len(minors)==1010
    print('exact residual 3x25;1010 distinct nonzero stripped minors',flush=True)
    script=output/'cover.m2'
    lines=['R=QQ[a,e,d,z,MonomialOrder=>GRevLex];',
           'F=500*a+90*e^2+3*e-24;B=(50*a+18*e^2-9*e-2)*d+60*e-40;C=100*a+18*e^2-9*e-2;',
           'I=ideal('+','.join(str(poly).replace('**','^') for poly in minors)+',z*a*F*B*C-1);',
           'g=gb I;assert(1%g==0);assert(dim(R/I)==-1);',
           'print "EXACT_QQ_UNIT_IDEAL_PASS";print toString gens g;']
    script.write_text('\n'.join(lines)+'\n')
    m2=shutil.which('M2')
    assert m2,'Activate the existing /workspace/cloud-setup/activate.sh first.'
    cp=subprocess.run([m2,'--script',str(script)],cwd=ROOT,capture_output=True,text=True,timeout=120)
    log=cp.stdout+cp.stderr
    (output/'cover.log').write_text(log)
    cp.check_returncode()
    assert 'EXACT_QQ_UNIT_IDEAL_PASS' in log and 'matrix {{1}}' in log
    inputs={str(SOURCE.relative_to(ROOT)):sha(SOURCE),**data['input_sha256'],
            str(Path(__file__).relative_to(ROOT)):sha(__file__)}
    report=dict(id='S61_BINARY_NONSQUARE_QUADRATIC_COVER_20261103_v62',all_passed=True,
                field='QQ; no finite-field lifting is used',core_variables=['a','e','d'],
                known_open='a*F*B*C!=0',quadratic_rows=data['quadratic_rows'],
                full_quadratic_matrix_shape=[6,28],replayed_unit_pivots=data['unit_pivots'],
                residual_shape=[3,25],cleared_row_denominators=denominators,
                minor_count=len(minors),minor_selections=selections,
                localized_residual_minor_ideal_is_unit=True,all_core_quadratic_rank=6,
                consequence='No quadratic rank-drop curves or points exist on the stated three-core open.',
                limitations=['Cubic and quartic rank coverage are separate.',
                             'No new ambient t-nilpotence or full marked exhaustion follows from rank6 alone.'],
                input_sha256=inputs,generated_sha256={p.name:sha(p) for p in [script,output/'cover.log']},
                seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'all_core_quadratic_rank':6,'seconds':report['seconds']}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)

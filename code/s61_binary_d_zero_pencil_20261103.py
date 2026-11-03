#!/usr/bin/env python3
"""Rebuild the d=0 finite pencil and its polynomial Smith/Macaulay certificate."""
import argparse
import functools
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

import sympy as s
from sympy.matrices.normalforms import smith_normal_decomp
from sympy.polys.rings import ring
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_binary_endpoint_full_transverse_20261102 import CUTS, monomials, ROOT
from s61_binary_zero_pair_stratum_removal_20261102 import write_quartics, compile_span

P=32003
A=s.Symbol('a')


def integral_polynomial(value):
    poly=s.Poly(value,A,domain=s.QQ)
    assert all(int(c.q)%P for c in poly.coeffs())
    return poly


def unit(value):
    value=s.Rational(value)
    assert int(value.p)%P and int(value.q)%P
    return value


def raw_pencil():
    model=S61()
    chart=json.loads((ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json').read_text())
    linear=json.loads((ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json').read_text())['linear']
    fixed,columns=linear['fixed'],linear['cols']
    weights=chart['local_coordinate_weights']+[1,1]
    grades=chart['equation_weights']
    normal=s.symbols('n0 n1')
    one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    beta=dict.fromkeys(one,s.S.Zero)
    beta.update({3:(-s.Rational(7,5)-A**2)/24,13:(5*A+10)/45,
                 34:1,35:0,48:(5*A+1)/9,49:-5,11:normal[0],12:normal[1]})
    _,equations=specialize(model,beta,True)
    rows=[]
    for grade,equation in zip(grades,equations):
        row={}
        for monomial,coefficient in equation.items():
            ids=tuple(columns.index(i) for i in monomial if i!=fixed)+(44,)*(2-len(monomial))
            if any(i in CUTS for i in ids):
                continue
            for (n0,n1),value in s.Poly(coefficient,*normal).terms():
                mon=tuple(sorted(ids+(45,)*n0+(46,)*n1))
                row[mon]=row.get(mon,0)+value
        rows.append({mon:value for mon,value in row.items() if value})
        assert all(sum(weights[i] for i in mon)==grade for mon in rows[-1])
        for value in rows[-1].values():
            integral_polynomial(value)
    assert all(w>0 for i,w in enumerate(weights) if i not in CUTS)
    return rows,weights,grades


def verify(output):
    started=time.monotonic()
    output.mkdir(parents=True,exist_ok=False)
    raw,weights,grades=raw_pencil()
    blocks=[];free_weights=[];metadata=[]
    for grade in range(1,5):
        rows=[i for i,w in enumerate(grades) if w==grade]
        columns=[i for i,w in enumerate(weights) if w==grade and i not in CUTS]
        matrix=s.Matrix([[raw[i].get((j,),0) for j in columns] for i in rows])
        diagonal,left,right=smith_normal_decomp(matrix,domain=s.QQ.poly_ring(A))
        assert (left*matrix*right-diagonal).applyfunc(s.expand)==s.zeros(*matrix.shape)
        for value in list(left)+list(right)+list(diagonal):
            integral_polynomial(value)
        ld=unit(left.det(method='domain-ge'));rd=unit(right.det(method='domain-ge'))
        constant=[];nonconstant=[]
        for i in range(len(columns)):
            value=diagonal[i,i] if i<min(diagonal.shape) else 0
            if value and value.is_Rational:
                unit(value);constant.append(i)
            else:
                free_weights.append(grade)
                if value:
                    nonconstant.append(str(s.factor(value)))
        metadata.append(dict(grade=grade,shape=list(matrix.shape),
                             left_determinant=str(ld),right_determinant=str(rd),
                             constant_pivots=[str(diagonal[i,i]) for i in constant],
                             nonconstant_diagonal=nonconstant,
                             left=[[str(v) for v in row] for row in left.tolist()],
                             right=[[str(v) for v in row] for row in right.tolist()]))
        blocks.append((rows,columns,diagonal,left,right))
        print('SMITH',grade,matrix.shape,'SECONDS',time.monotonic()-started,flush=True)
    assert free_weights==[1]*7+[2]*2
    (output/'smith.json').write_text(json.dumps(metadata,indent=2)+'\n')
    R,*x=ring('a,'+','.join(f'x{i}' for i in range(len(free_weights))),s.GF(P))
    @functools.lru_cache(None)
    def convert(value):
        result=R.zero
        for (power,),coefficient in integral_polynomial(value).terms():
            residue=int(coefficient.p)*pow(int(coefficient.q),-1,P)%P
            result+=R(residue)*x[0]**power
        return result
    coordinates=[R.zero]*47
    @functools.lru_cache(None)
    def product(monomial):
        result=R.one
        for i in monomial:
            result*=coordinates[i]
        return result
    def evaluate(row):
        result=R.zero
        for monomial,coefficient in row.items():
            if all(coordinates[i] for i in monomial):
                result+=convert(coefficient)*product(monomial)
        return result
    counter=0
    for rows,columns,diagonal,left,right in blocks:
        values=[evaluate(raw[i]) for i in rows]
        transformed=[sum((convert(c)*value for c,value in zip(row,values)),R.zero)
                     for row in left.tolist()]
        changed=[]
        for j in range(len(columns)):
            value=diagonal[j,j] if j<min(diagonal.shape) else 0
            if value and value.is_Rational:
                changed.append(-transformed[j]*convert(1/value))
            else:
                counter+=1;changed.append(x[counter])
        new=[sum((convert(c)*value for c,value in zip(row,changed)),R.zero)
             for row in right.tolist()]
        for j,value in zip(columns,new):
            coordinates[j]=value
        residual=[evaluate(raw[i]) for i in rows]
        check=[sum((convert(c)*value for c,value in zip(row,residual)),R.zero)
               for row in left.tolist()]
        for j in range(min(diagonal.shape)):
            if diagonal[j,j] and diagonal[j,j].is_Rational:
                assert not check[j]
        print('ELIMINATED',weights[columns[0]],'TERMS',sum(len(z) for z in new),
              'SECONDS',time.monotonic()-started,flush=True)
    residual=[evaluate(row) if grade<=4 else R.zero for row,grade in zip(raw,grades)]
    print('RESIDUAL_TERMS',sum(len(z) for z in residual),'SECONDS',time.monotonic()-started,flush=True)
    matrix=output/'weight4.in'
    nrows=write_quartics(matrix,residual,free_weights,grades)
    lines=matrix.read_text().splitlines()
    nr,nc=map(int,lines[0].split()[1:])
    assert nr==nrows==561 and nc==269
    matrix.write_text(lines[0]+' 0\n0 1\n'+'\n'.join(lines[1:])+'\n')
    source=ROOT/'code/s61_binary_hyperplane_vertical_rowmodule_20261102.cpp'
    exe=output/'row-module'
    compile_span(source,exe)
    log=output/'elimination.log';result=output/'determinant.txt'
    with log.open('w') as stream:
        subprocess.run([str(exe.resolve()),str(matrix.resolve()),str(result.resolve())],
                       stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=45)
    assert result.read_text()=='ROWMODULE 269 DEGREE 0 STEPS 1325357\nDETERMINANT_COEFFICIENTS 1 1\nRESIDUAL_SIZE 0 FACTORS 0\n'
    exe.unlink()
    (output/'smith.json').write_text(json.dumps(metadata,indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'code/s61_model.py',
            ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
            ROOT/'code/s61_binary_endpoint_full_transverse_20261102.py',source,
            ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json',
            ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json',
            ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json',
            ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    report=dict(id='S61_BINARY_D_ZERO_PENCIL_20261103',all_passed=True,prime=P,
                shape=[nr,nc],free_weights=free_weights,row_operations=1325357,
                matrix_sha256=hashlib.sha256(matrix.read_bytes()).hexdigest(),
                polynomial_smith_unimodularity_over_V_checked=True,
                nonconstant_smith_entries_not_inverted=True,
                all_finite_mod_p_central_algebras_finite=True,
                endpoint_handled_separately=True,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=verify(args.output)
    print('PASS',result['id'],result['matrix_sha256'],result['seconds'])

#!/usr/bin/env python3
"""Exact small heavy quotients at two exceptional marked-square quadratics.

The model is the frozen all-h full135 residual, with rho and seven lights set
to zero. Positive-dimensional support here is a central quotient calculation;
it does not identify the ambient saturated closure.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from s61_marked_sigma_pencil_rank_20261103 import (
    ROOT, PRIME, monomials, compile_native, sha)


def read_exact(path):
    reference = json.loads((ROOT/'reports/S61_MARKED_SIGMA_PENCIL_MODEL_20261102.json').read_text())
    assert sha(path) == reference['row_output_sha256']
    rows = []
    with path.open() as source:
        assert tuple(map(int,source.readline().split())) == (14,136)
        for i in range(136):
            weight,count = map(int,source.readline().split())
            row = defaultdict(dict)
            for _ in range(count):
                scalar,*powers = source.readline().split()
                powers = tuple(map(int,powers))
                assert len(powers) == 15
                if powers[-1] or any(powers[1:8]):
                    continue
                exponent,h = powers[8:14],powers[0]
                scalar = Fraction(scalar)
                assert scalar.denominator % PRIME
                row[exponent][h] = row[exponent].get(h,0)+scalar
            row = {ex:{h:c for h,c in cs.items() if c} for ex,cs in row.items()}
            rows.append((i,weight,{ex:cs for ex,cs in row.items() if cs}))
        assert not source.read().strip()
    assert rows[-1][2] == {(0,0,0,0,1,0):{0:Fraction(1)}}
    return rows


def minor_certificate(rows):
    weights = (2,2,2,2,3,3)
    columns = sorted(monomials(weights,6))
    index = {mon:i for i,mon in enumerate(columns)}
    records = []
    for i,weight,row in rows[:-1]:
        if weight > 6 or not row:
            continue
        for shift in monomials(weights,6-weight):
            value = {index[tuple(a+b for a,b in zip(ex,shift))]:cs for ex,cs in row.items()}
            records.append((max(max(cs) for cs in value.values()),len(value),i,shift,value))
    records.sort(key=lambda r:r[:4])
    h = sp.Symbol('h')
    matrix = sp.zeros(len(records),len(columns))
    numeric = sp.zeros(*matrix.shape)
    for i,record in enumerate(records):
        for j,cs in record[-1].items():
            matrix[i,j] = sum(sp.Rational(c.numerator,c.denominator)*h**k for k,c in cs.items())
            c = cs.get(0,Fraction(0))
            numeric[i,j] = c.numerator*pow(c.denominator,-1,PRIME)%PRIME
    selected = list(DomainMatrix.from_Matrix(numeric.T).convert_to(sp.GF(PRIME)).rref()[1])
    assert len(selected) == 23 and matrix.shape == (140,23)
    determinant = matrix[selected,:].to_DM().convert_to(sp.QQ.poly_ring(h)).det().as_expr()
    factorization = sp.factor_list(determinant)
    assert sp.Poly(determinant,h).degree() == 26
    assert sp.rem(sp.Poly(determinant,h),sp.Poly(1000*h*h+1127,h)).is_zero
    return dict(shape=list(matrix.shape),selected_rows=selected,
                determinant_QQ=str(determinant),factorization_QQ=str(factorization),
                degree=26,scope='one explicit nonzero minor; not the complete QQ Fitting ideal')


def m2_source(rows,quadratic,prime,root=None):
    base = f'(ZZ/{prime})' if prime else 'QQ'
    polynomial = ('1000*h^2+1127' if quadratic == 'new' else 'h^2+2')
    if prime and quadratic == 'old':
        # h^2+2 splits in F_32003; do not misdeclare its reducible quotient a field.
        assert root is not None and (root*root+2)%prime == 0
        source = f'K={base};a={root};\n'
    else:
        source = f'A={base}[h]; K=toField(A/ideal({polynomial})); a=promote(h,K);\n'
    source += 'R=K[x0,x1,x2,x3,x4,x5,Degrees=>{2,2,2,2,3,3}];\n'
    def encode(row):
        result = []
        for ex,cs in row.items():
            coefficients = []
            for h,c in sorted(cs.items()):
                value = (str(c.numerator*pow(c.denominator,-1,prime)%prime) if prime
                         else f'({c.numerator}/{c.denominator})')
                coefficients.append(value+(f'*a^{h}' if h else ''))
            result.append('('+'+'.join(coefficients)+')'+''.join(
                f'*x{i}^{power}' for i,power in enumerate(ex) if power))
        return '+'.join(result) or '0'
    source += 'I=ideal('+','.join(encode(row) for i,w,row in rows[:-1])+'); G=gb I;\n'
    source += 'T='+encode(rows[-1][2])+'; assert(T==x4); assert(T^2%G==0);\n'
    if quadratic == 'old':
        source += 'J=ideal(x0,x1,x2-a*x3,x4); assert(dim(R/I)==2);\n'
    else:
        source += ('G1=gb(I+ideal(x3-1)); assert(numColumns basis(R/(I+ideal(x3-1)))==1);\n'
                   'c0=x0%G1;c1=x1%G1;c2=x2%G1;\n'
                   'J=ideal(x0-c0*x3,x1-c1*x3,x2-c2*x3,x4,x5);\n'
                   'assert(dim(R/I)==1); print("LINE="|toString{c0,c1,c2});\n')
    # Both J are explicit homogeneous polynomial linear subspaces, so radical.
    source += ('assert(isSubset(I,J));\n'
               'assert all(flatten entries gens J,f->any(1..8,n->f^n%G==0));\n'
               'print("DIM="|toString dim(R/I));\n'
               'print("HS="|toString hilbertSeries(R/I,Reduce=>true));\n'
               'print("RADICAL="|toString gens J); print("GB="|toString gens G);\n'
               'print "HEAVY SUPPORT PASS: t=x4, t^2=0";exit 0;\n')
    return source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--full-module',type=Path)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False,parents=True)
    started=time.monotonic();rows=read_exact(args.rows)
    minor=minor_certificate(rows)
    (args.output/'minor_QQ.json').write_text(json.dumps(minor,indent=2)+'\n')
    m2=shutil.which('M2');assert m2
    results=[]
    cases=[('new',0,None,'new-QQ'),('new',PRIME,None,'new-modp'),
           ('old',0,None,'old-QQ'),('old',PRIME,3083,'old-modp-plus'),
           ('old',PRIME,28920,'old-modp-minus')]
    for quadratic,prime,root,name in cases:
        script=args.output/f'{name}.m2';script.write_text(m2_source(rows,quadratic,prime,root))
        log=args.output/f'{name}.log'
        with log.open('w') as stream:
            result=subprocess.run([m2,'--script',str(script.resolve())],stdout=stream,
                                  stderr=subprocess.STDOUT,timeout=120,check=True)
        content=log.read_text();assert 'HEAVY SUPPORT PASS' in content
        results.append(dict(case=name,returncode=result.returncode,output=content))
        print(name,content,flush=True)
    numeric_results=[]
    if args.full_module:
        lines=args.full_module.read_text().splitlines();header=list(map(int,lines[0].split()))
        assert header == [PRIME,2603,1432,0]
        assert all(lines[-8+2*i]=='1' for i in range(4))
        header[1]-=4;lines[0]=' '.join(map(str,header))
        original=args.output/'full5_original.in';original.write_text('\n'.join(lines[:-8])+'\n')
        executable=args.output/'numeric';compile_native(ROOT/'code/s61_marked_sigma_numeric_complement_20261103.cpp',executable)
        for h,expected in [(0,4),(4499,5),(27504,5),(15938,5),(16065,5),(4,6),(31999,6)]:
            values=list(map(int,subprocess.run([str(executable),str(original),str(h)],
                                             capture_output=True,text=True,check=True,timeout=30).stdout.split()))
            assert len(values[1:])==1432-values[0]==expected
            numeric_results.append(dict(h=h,rank=values[0],cokernel_dimension=expected,nonpivot_columns=values[1:]))
        executable.unlink()
    inputs=[Path(__file__),ROOT/'code/s61_marked_sigma_pencil_rank_20261103.py',
            ROOT/'code/s61_marked_sigma_numeric_complement_20261103.cpp',
            ROOT/'reports/S61_MARKED_SIGMA_PENCIL_MODEL_20261102.json']
    report=dict(id='S61_MARKED_SIGMA_HEAVY_EXCEPTIONS_20261103',all_passed=True,
                global_rows_sha256=sha(args.rows),minor_QQ=minor,cases=results,
                full_S5_numeric_ranks_mod_32003=numeric_results,
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs},
                evidence_sha256={str(p.relative_to(args.output)):sha(p)
                                 for p in args.output.rglob('*') if p.is_file()},
                whole_pencil_removal_proved=False,ambient_nilpotence_proved=False,
                seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()

#!/usr/bin/env python3
"""Certify removal over general beta of the entire c=f=0 nonsquare binary stratum."""
import argparse
import hashlib
import os
import json
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
import sympy as sp
from s61_binary_endpoint_full_transverse_20261102 import (
    ROOT,PRIME,CUTS,A,B,D,E,load_report,original_transversal,reduce_test_line,monomials)

STEM = 'S61_BINARY_ZERO_PAIR_STRATUM_REMOVAL_20261102'


def write_quartics(path, residual, weights, grades):
    columns = {mon:i for i,mon in enumerate(monomials(weights,4))}
    rows = []
    for row,grade in zip(residual,grades):
        if not row or grade > 4:
            continue
        for shift in monomials(weights,4-grade):
            coefficients = {}
            for monomial,c in row.items():
                column = columns[tuple(i+j for i,j in zip(monomial[1:],shift))]
                coefficients.setdefault(column,{})[monomial[0]] = int(c) % PRIME
            rows.append(coefficients)
    with path.open('w') as stream:
        print(PRIME,len(rows),len(columns),file=stream)
        for row in rows:
            print(len(row),file=stream)
            for col,coefficients in sorted(row.items()):
                degree = max(coefficients)
                print(col,degree+1,*[coefficients.get(i,0) for i in range(degree+1)],file=stream)
    return len(rows)


def infinity_program(raw,weights):
    variables = [i for i in range(47) if i not in CUTS]
    substitution = {A:1,B:0,D:12,E:0}
    def expression(row):
        result = []
        for mon,c in row.items():
            c = c.subs(substitution)
            if c:
                result.append(f'({c})'+''.join(f'*q{i}' for i in mon))
        return '+'.join(result) or '0'
    return ('R=(ZZ/32003)['+','.join(f'q{i}' for i in variables)+',Degrees=>{'
            +','.join(str(weights[i]) for i in variables)+'}];I=ideal('
            +','.join(map(expression,raw))+');print toString hilbertSeries(R/I,Reduce=>true);')


def compile_span(source, exe):
    configured = subprocess.run(['pkg-config','--cflags','--libs','flint'],
                                text=True,capture_output=True) if shutil.which('pkg-config') else None
    if configured is not None and configured.returncode == 0:
        import shlex
        flags = shlex.split(configured.stdout)
    else:
        flags = ['-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),
                 '-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),'-I/opt/homebrew/include',
                 str(ROOT/'work/m2-runtime/lib/libflint.24.0.dylib'),'-L/opt/homebrew/lib','-lgmp']
    subprocess.run(['c++','-std=c++17','-O3',str(source),*flags,'-o',str(exe)],check=True,timeout=45)


def verify():
    start=time.monotonic()
    raw,weights,grades,orbit_columns,orbit_minor = original_transversal()
    residual,free_weights,determinants = reduce_test_line(raw,weights,grades)
    source=ROOT/'code/s61_binary_endpoint_polynomial_row_span_20261102.cpp'
    env=os.environ.copy()
    if (ROOT/'work/m2-runtime/lib').is_dir():
        env['DYLD_LIBRARY_PATH']=':'.join([str(ROOT/'work/m2-runtime/lib'),'/opt/homebrew/lib',
                                          '/opt/homebrew/opt/givaro/lib',env.get('DYLD_LIBRARY_PATH','')])
    with tempfile.TemporaryDirectory(prefix='s61-binary-zero-pair-') as tmp:
        tmp=Path(tmp)
        input_path=tmp/'quartics'
        nrows=write_quartics(input_path,residual,free_weights,grades)
        matrix_hash=hashlib.sha256(input_path.read_bytes()).hexdigest()
        exe=tmp/'span'
        compile_span(source,exe)
        run=subprocess.run([str(exe),str(input_path)],env=env,capture_output=True,text=True,check=True,timeout=90)
        nr,nc,operations,det=map(int,run.stdout.split())
        assert nr==nrows and nc==239 and det
        script=tmp/'infinity.m2'
        script.write_text(infinity_program(raw,weights))
        m2=os.getenv('M2_BINARY') or shutil.which('M2-binary') or shutil.which('M2')
        result=subprocess.run([m2,'--script',str(script)],cwd=tmp,env=env,capture_output=True,text=True,check=True,timeout=40)
        answer=sp.sympify(result.stdout.strip().replace('^','**'),locals={'T':sp.Symbol('T')})
        z=sp.Symbol('T')
        assert sp.cancel(answer-(1+7*z+23*z**2+34*z**3+5*z**4+2*z**5))==0
    names = [sp.Symbol('u'+str(i)) for i in [4,5,12,13,14,25,26,27,34,35,36,37,48,49,50,51]]
    constraint=[names.index(sp.Symbol('u'+str(j))) for j in [12,13,35,50,36]]
    group=sp.Matrix(load_report('S61_BINARY_RANK_DROP_FAMILY_20261101')['primitive_action']).applyfunc(sp.sympify)
    constrained=group[constraint,:].subs({names[i]:0 for i in constraint})
    assert constrained.to_DM().rank()==4
    files=[Path(__file__),source,ROOT/'code/s61_binary_endpoint_full_transverse_20261102.py',
           ROOT/'code/toric_closed_jet_coordinates_20261102.py',ROOT/'research/notes'/f'{STEM}.md',
           ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json',
           ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    for s in ['S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102','S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101',
              'S61_BINARY_RANK_DROP_FAMILY_20261101','S61_BINARY_NUISANCE_CONTRACTION_20261102',
              'S61_UNIVERSAL_WEIGHT6_20261031']:
        assert load_report(s)['all_passed']
        files.append(ROOT/'reports'/f'{s}.json')
    return dict(id=STEM,all_passed=True,orbit_cuts=CUTS,orbit_minor_columns=orbit_columns,orbit_minor=str(orbit_minor),
                positive_beta_coordinates=['u12','u13'],test_line={'a':'s^3','b':6,'c':0,'d':12,'e':16,'f':0},
                constant_pivot_determinants=determinants,weight4_row_module_shape=[nr,nc],
                polynomial_row_operations=operations,constant_basis_determinant_mod_32003=det,
                affine_matrix_sha256=matrix_hash,infinity_chart_core=[1,0,12,0],infinity_length_mod_32003=72,
                infinity_hilbert=[1,7,23,34,5,2],d_zero_source_normal_rank=4,d_zero_beta_image_dimension_upper=15,
                entire_c_f_zero_stratum_removed_over_general_beta=True,remaining_marked_c_nonzero_H11=True,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                seconds=time.monotonic()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=verify()
    with args.output.open('x') as stream:
        json.dump(report,stream,indent=2)
        stream.write('\n')
    print('PASS',report['weight4_row_module_shape'],report['constant_basis_determinant_mod_32003'],report['seconds'])

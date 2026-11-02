#!/usr/bin/env python3
"""Exact weighted cubic rank64 off D and seventeen explicitly treated curves."""
import argparse
import hashlib
from itertools import combinations_with_replacement
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
import sympy as s
from sympy.polys.rings import ring
from s61_binary_weighted_seventeen_curves_20261102 import prepare
from s61_binary_deeper_uniform_quadratic_20261101 import PIVOT_ROWS
ROOT=Path(__file__).resolve().parents[1]
NOTE=ROOT/'research/notes/S61_BINARY_WEIGHTED_CUBIC_COVER_20261102.md'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_poly(stream,p):
    stream.write(str(len(p))+'\n')
    for m,c in p.items():stream.write(str(c)+' '+' '.join(map(str,m))+'\n')


def pair_inputs(wd,pair):
    K=s.QQ.frac_field('e','d');C=K.field.ring
    R,*_=ring(','.join('x'+str(i)for i in range(7)),K)
    quadrics=[R.from_expr(s.sympify(q)) for q in pair['new_quadrics']]
    images=[R.from_expr(s.sympify(q)) for q in pair['new_to_old_linear']]
    monomials=[]
    for indices in combinations_with_replacement(range(7),2):
        m=[0]*7
        for i in indices:m[i]+=1
        monomials.append(tuple(m))
    index={m:i for i,m in enumerate(monomials)}
    qfile=wd/'pair_quadrics';cfile=wd/'pair_coordinates'
    with qfile.open('w')as stream:
        stream.write('7 28\n')
        for q in quadrics:
            stream.write(str(len(q))+'\n')
            for m,c in sorted(q.items(),key=lambda v:index[v[0]]):
                stream.write(str(index[m])+'\n');write_poly(stream,c.numer);write_poly(stream,c.denom)
    denominator=C.one
    for q in images:
        for c in q.values():denominator=denominator.lcm(c.denom)
    with cfile.open('w')as stream:
        write_poly(stream,denominator)
        for q in images:
            numerator={}
            for m,c in q.items():
                for n,v in (c.numer*denominator.exquo(c.denom)).items():numerator[n+m]=v
            write_poly(stream,numerator)
    return qfile,cfile


def verify():
    started=time.time()
    finite_path=ROOT/'reports/S61_BINARY_DEEPER_WEIGHTED_FINITE_MODEL_20261102.json'
    blocks_path=ROOT/'reports/S61_BINARY_DEEPER_WEIGHTED_LOW_BLOCKS_20261102.json'
    pair_path=ROOT/'reports/S61_BINARY_DEEPER_PRODUCT_PAIR_20261102.json'
    curves_path=ROOT/'code/data/s61_binary_weighted_seventeen_curves_20261102.json'
    finite=json.loads(finite_path.read_text());blocks=json.loads(blocks_path.read_text());pair=json.loads(pair_path.read_text())
    curves=json.loads(curves_path.read_text())['curves']
    assert finite['all_passed'] and blocks['all_passed'] and pair['all_passed']
    assert pair['old_light_indices']==[0,1,2,3,4,6,8]
    assert set(PIVOT_ROWS[i] for i in pair['old_quadratic_rows'])==set(blocks['blocks']['2']['rows'])
    native=['code/parametric_implicit_jets_20261101.cpp','code/s61_binary_deeper_weighted_extract_20261102.cpp',
            'code/s61_binary_weighted_cubic_schur_20261102.cpp','code/s61_binary_weighted_cubic_pair_20261102.cpp',
            'code/s61_binary_weighted_cubic_basis_20261102.cpp','code/s61_binary_weighted_cubic_dual_20261102.cpp',
            'code/s61_binary_weighted_cubic_dual_check_20261102.cpp','code/s61_binary_weighted_cubic_support_20261102.cpp']
    logs=[]
    with tempfile.TemporaryDirectory(prefix='s61-cubic-cover-')as temp:
        wd=Path(temp);inp,raw,weighted,meta,weights,coeff,known=prepare(wd)
        qfile,cfile=pair_inputs(wd,pair)
        config=wd/'schur_config';knownfile=wd/'known';allfile=wd/'all_factors'
        with config.open('w')as stream:
            stream.write('14\n')
            for row,indices in zip(blocks['blocks']['3']['rows'],blocks['blocks']['3']['monomials']):
                mon=[0]*9
                for j in indices:mon[j]+=1
                stream.write(' '.join(map(str,[row,*mon]))+'\n')
        with knownfile.open('w')as stream:
            stream.write('13\n')
            for p in known:write_poly(stream,p)
        with allfile.open('w')as stream:
            stream.write('30\n')
            for p in known:write_poly(stream,p)
            for curve in curves:write_poly(stream,{tuple(m):s.QQ(c)for m,c in curve['sparse']})
        env=dict(os.environ);env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'
        binaries=[]
        for k,source in enumerate(native):
            binary=wd/f'engine{k}';binaries.append(binary)
            subprocess.run(['c++','-O3','-std=c++17','-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),
                '-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),'-I/opt/homebrew/include',str(ROOT/source),
                str(ROOT/'work/m2-runtime/lib/libflint.24.0.dylib'),'-L/opt/homebrew/lib','-lgmp','-o',str(binary)],check=True)
        schur=wd/'schur';transformed=wd/'transformed';basis=wd/'basis';dual=wd/'dual';support=wd/'support.json'
        commands=[[binaries[0],inp,raw],[binaries[1],raw,coeff,weighted,meta,weights],
            [binaries[2],weighted,config,knownfile,schur],
            [binaries[3],schur,cfile,qfile,transformed],
            [binaries[4],qfile,basis,'3'],
            [binaries[5],basis,transformed,dual,'35','74'],
            [binaries[6],transformed,str(dual)+'.dual'],
            [binaries[7],allfile,str(basis)+'.pivots',dual,str(dual)+'.dual',support]]
        for k,command in enumerate(commands):
            run=subprocess.run(list(map(str,command)),check=True,capture_output=True,text=True,env=env)
            logs.append(run.stderr);print('exact stage',k,'PASS',round(time.time()-started,3),flush=True)
        generated={k:sha(p)for k,p in [('raw_input',inp),('raw_output',raw),('weighted_rows',weighted)]}
        assert generated==finite['generated_sha256']
        factor_support=json.loads(support.read_text());assert factor_support['all_passed']
        artifacts={k:sha(p)for k,p in [('pair_quadrics',qfile),('pair_coordinates',cfile),('schur',schur),
            ('transformed',transformed),('quadratic_basis',basis),('quadratic_pivots',Path(str(basis)+'.pivots')),
            ('dual_minor',dual),('dual_functionals',Path(str(dual)+'.dual'))]}
    paths=[Path(__file__),NOTE,finite_path,blocks_path,pair_path,curves_path,*[ROOT/p for p in native],
        ROOT/'code/s61_binary_deeper_quartic_field_20261101.hpp',
        ROOT/'code/s61_binary_deeper_uniform_quadratic_20261101.py',
        ROOT/'code/s61_binary_weighted_seventeen_curves_20261102.py',
        ROOT/'reports/S61_BINARY_WEIGHTED_SEVENTEEN_CURVES_20261102.json',
        ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json',
        ROOT/'reports/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101.json',
        ROOT/'reports/S61_BINARY_DEEPER_CORE_OPEN_FORMAL_20261101.json',
        ROOT/'research/notes/S61_BINARY_DEEPER_CORE_OPEN_FORMAL_20261101.md',
        ROOT/'research/notes/S61_BINARY_DEEPER_WEIGHTED_FINITE_MODEL_20261102.md']
    return dict(id='S61_BINARY_WEIGHTED_CUBIC_COVER_20261102',all_passed=True,
        nuisance_parameters_zero=True,known_factor_count=13,additional_explicit_curve_count=17,
        heavy_light_invertible_block=14,pure_light_quadratic_products=49,pure_light_product_rank_on_open=48,
        true_pure_cubic_rows=[92,134],dual_free_columns=[35,74],exact_dual_annihilation_identities=98,
        two_by_two_pairing_determinant_unit_on_open=True,weighted_cubic_rank_on_open=64,
        weighted_cubic_quotient_dimension_on_open=34,
        universal_rank_upper64_dependency='Generic free65 from CORE_OPEN_FORMAL, weight2 rank7, and generic weight4 vanishing in the exact weighted model.',
        factor_support=factor_support,generated_sha256=generated,artifact_sha256=artifacts,native_logs=logs,
        proof_note=str(NOTE.relative_to(ROOT)),
        limitations=['The weighted rank64 statement is on nuisance=0 and D*product(f_i) nonzero.',
            'No new uniform rank269 or all-core quartic cover is proved.',
            'Extending the two cubic classes to every nuisance value requires a separate denominator audit.',
            'No affine exhaustion or component upper bound follows.'],
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},seconds=time.time()-started)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args()
    assert not args.output.exists();report=verify()
    with args.output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'])))

#!/usr/bin/env python3
"""Localize the two true cubic classes uniformly over all five nuisance variables."""
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
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]
NOTE=ROOT/'research/notes/S61_BINARY_CUBIC_NUISANCE_LOCALIZATION_20261102.md'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def write_poly(stream,p):
    stream.write(str(len(p))+'\n')
    for m,c in p.items():stream.write(str(c)+' '+' '.join(map(str,m))+'\n')


def specialized_contents(data,known):
    R,*_=ring(','.join(data['variables']),s.QQ);field=R.to_field();C,*_=ring('e,d',s.QQ)
    z=C.zero
    core=lambda p:C.from_dict({m[:2]:c for m,c in p.items()if not any(m[2:])})
    T=[[field.from_expr(s.sympify(v))for v in row]for row in data['T']]
    U=[[field.from_expr(s.sympify(v))for v in row]for row in data['U']]
    Dt=[];Du=R.one
    for j in data['normal']:
        den=R.one
        for row in T:den=den.lcm(row[j].denom)
        Dt.append(den)
    for row in U:
        for v in row:Du=Du.lcm(v.denom)
    # These are the FULL seven-parameter denominators, before specialization.
    powers=[]
    for den in [*Dt,Du]:
        assert all(not any(m[2:]) for m in den)
        p=core(den);record=[]
        for factor in known:
            power=0
            while p.rem(factor)==0:p=p.exquo(factor);power+=1
            record.append(power)
        assert p and all(sum(m)==0 for m in p);powers.append(record)
    Tn=[[core(row[j].numer*Dt[k].exquo(row[j].denom))for k,j in enumerate(data['normal'])]for row in T]
    Un=[[core(v.numer*Du.exquo(v.denom))for v in row]for row in U]
    assert all(v==0 for v in Tn[44]);den=core(Du)
    J=[[core(R.from_expr(s.sympify(v)))for v in row]for row in data['J']]
    nuisance={s.Symbol(v):0 for v in data['variables'][2:]}
    beta={int(k[1:])-1:s.sympify(v).subs(nuisance)for k,v in data['beta'].items()}
    _,equations=specialize(S61(),beta,True)
    mons=list(combinations_with_replacement(range(9),2));index={m:i for i,m in enumerate(mons)};cache={}
    def product(i,j):
        key=tuple(sorted((i,j)))
        if key not in cache:
            row=[z]*45
            for a,v in enumerate(Tn[i]):
                if not v:continue
                for b,w in enumerate(Tn[j]):
                    if w:row[index[tuple(sorted((a,b)))]]+=v*w
            cache[key]=row
        return cache[key]
    selected=data['pivot_rows'];required=set(selected)|{92,134};Q={}
    for i in required:
        row=[z]*45
        for mon,value in equations[i].items():
            if data['fixed'] in mon:continue
            a,b=[data['cols'].index(j)for j in mon]+[44]*(2-len(mon))
            coefficient=C.from_expr(value)
            if not coefficient:continue
            for k,v in enumerate(product(a,b)):
                if v:row[k]+=coefficient*v
        Q[i]=row
    contents={}
    for i in [92,134]:
        projection=[sum((v*Un[j][k]for j,v in enumerate(J[i])if v),z)for k in range(33)]
        A=[den*Q[i][a]-sum((v*Q[selected[k]][a]for k,v in enumerate(projection)if v),z)for a in range(45)]
        g=z
        for v in A:g=g.gcd(v)
        assert g;contents[i]=g
    return contents,powers


def verify():
    started=time.time()
    cert_path=ROOT/'code/certificates/s61_binary_deeper_cubic_20261101.json'
    linear_path=ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json'
    quadratic_path=ROOT/'reports/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101.json'
    certificate=json.loads(cert_path.read_text());data=json.loads(linear_path.read_text())['linear']
    R,e,d=ring('e,d',s.QQ);field=R.to_field()
    known=[R.from_expr(s.sympify(v))for v in json.loads(quadratic_path.read_text())['irreducible_minor_factors']]
    rref=[field.from_expr(s.sympify(v))for row in certificate['quadratic_rref']['rref']for v in row if v!='0']
    relations=[]
    for row in certificate['two_quadratic_row_relations']['relations'].values():
        for polynomial in row:
            for mon,c in polynomial:
                relations.append(tuple(R.from_dict({tuple(m):s.QQ(v)for m,v in c[key]})for key in ['numerator','denominator']))
    assert len(rref)==169 and len(relations)==46
    contents,linear_powers=specialized_contents(data,known)
    native=ROOT/'code/s61_binary_cubic_nuisance_localization_20261102.cpp'
    with tempfile.TemporaryDirectory(prefix='s61-cubic-localization-')as temp:
        wd=Path(temp);factors=wd/'factors';coefficients=wd/'coefficients';out=wd/'native.json';binary=wd/'engine'
        with factors.open('w')as stream:
            stream.write('13\n')
            for p in known:write_poly(stream,p)
        with coefficients.open('w')as stream:
            stream.write('169 46\n')
            for c in rref:write_poly(stream,c.numer);write_poly(stream,c.denom)
            for numerator,denominator in relations:write_poly(stream,numerator);write_poly(stream,denominator)
            stream.write('2\n')
            for i,g in contents.items():stream.write(str(i)+'\n');write_poly(stream,g)
        subprocess.run(['c++','-O3','-std=c++17','-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),
            '-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),'-I/opt/homebrew/include',str(native),
            str(ROOT/'work/m2-runtime/lib/libflint.24.0.dylib'),'-L/opt/homebrew/lib','-lgmp','-o',str(binary)],check=True)
        env=dict(os.environ);env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'
        run=subprocess.run(list(map(str,[binary,factors,coefficients,out])),check=True,capture_output=True,text=True,env=env)
        result=json.loads(out.read_text());assert result['all_passed']
        generated=sha(coefficients)
    paths=[Path(__file__),native,NOTE,cert_path,linear_path,quadratic_path,
        ROOT/'code/s61_binary_deeper_quartic_field_20261101.hpp',ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
        ROOT/'code/s61_binary_deeper_cubic_nuisance_20261101.py',ROOT/'code/s61_binary_deeper_cubic_field_20261101.cpp',
        ROOT/'reports/S61_BINARY_DEEPER_CUBIC_NUISANCE_20261101.json',
        ROOT/'reports/S61_BINARY_DEEPER_WEIGHTED_LOW_BLOCKS_20261102.json',
        ROOT/'reports/S61_BINARY_WEIGHTED_CUBIC_COVER_20261102.json',
        ROOT/'reports/S61_BINARY_WEIGHTED_SEVENTEEN_CURVES_20261102.json',
        ROOT/'research/notes/S61_BINARY_DEEPER_CUBIC_NUISANCE_20261101.md']
    return dict(id='S61_BINARY_CUBIC_NUISANCE_LOCALIZATION_20261102',all_passed=True,
        coefficient_ring='QQ[e,d,D^-1][p,q,r,h,k]',full_seven_parameter_linear_denominators_are_core_only=True,
        full_linear_denominator_D_factor_exponents=linear_powers,
        exact_denominator_audit=result,specialized_raw_A2_contents={str(i):str(g.as_expr())for i,g in contents.items()},
        generated_coefficients_sha256=generated,native_log=run.stderr,
        nuisance_class_invariance_on_entire_D_open=True,
        independent_two_true_cubic_classes_on_DF_open_for_every_nuisance=True,
        quadratic_linear_rank_on_DF_open=129,rank_with_two_true_cubics_on_DF_open=131,
        known_initial_subideal_hilbert_through3=[1,9,21,34],
        limitations=['The ordinary rank131 is for the known two cubic constraints together with I2*S1, not an exhaustion of all true cubic relations.',
            'No uniform fourth-degree vanishing or all-core ambient removal follows.'],
        proof_note=str(NOTE.relative_to(ROOT)),input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},seconds=time.time()-started)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args()
    assert not args.output.exists();report=verify()
    with args.output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'])))

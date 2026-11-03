#!/usr/bin/env python3
"""Independent Fraction audit of two constant column witnesses in frozen V9."""
import argparse
from fractions import Fraction
import hashlib
import io
import json
from pathlib import Path
import tarfile
import sympy as s
from q_smith_global_model_review_20261103_v63 import ROOT, read, sha


def parse(line):
    words=line.split()
    n=int(words[0])
    assert len(words)==1+4*n
    result={}
    for k in range(n):
        mon=tuple(map(int,words[2+4*k:5+4*k]))
        assert mon not in result
        result[mon]=Fraction(words[1+4*k])
    return {m:c for m,c in result.items()if c}


def plus(left,right):
    result=left.copy()
    for mon,c in right.items():result[mon]=result.get(mon,Fraction(0))+c
    return {m:c for m,c in result.items()if c}


def times(left,right):
    result={}
    for mon,c in left.items():
        for other,d in right.items():
            key=tuple(a+b for a,b in zip(mon,other))
            result[key]=result.get(key,Fraction(0))+c*d
    return {m:c for m,c in result.items()if c}


def expr(poly):
    a,e,d=s.symbols('a e d')
    p=s.Poly(s.sympify(poly.replace('^','**'),locals={'a':a,'e':e,'d':d}),a,e,d,domain=s.QQ)
    return {m:Fraction(c)for m,c in p.terms()}


def stream_sha(path):
    result=hashlib.sha256()
    with path.open('rb')as f:
        while chunk:=f.read(1024*1024):result.update(chunk)
    return result.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    directory=ROOT/'research/handoff/20261103_v62/v9_short_unit_pivots'
    manifest=read(directory/'manifest.json')
    status=read(directory/'STATUS.json')
    archive=ROOT/manifest['archive']['path']
    assert sha(archive)==manifest['archive']['sha256']
    assert sha(directory/'STATUS.json')==manifest['status_sha256']
    for filename,digest in status['dependencies'].items():assert sha(ROOT/filename)==digest
    blobs={}
    with tarfile.open(archive)as tar:
        assert set(tar.getnames())=={r['path']for r in manifest['members']}
        for record in manifest['members']:
            data=tar.extractfile(record['path']).read()
            assert len(data)==record['bytes'] and hashlib.sha256(data).hexdigest()==record['sha256']
            blobs[record['path']]=data
    prefix='work/nonsquare-v62/primitive-unit-span-v9/'
    results=[]
    for record in status['cases']:
        name=record['case']
        source=ROOT/record['source']
        assert source.is_file() and stream_sha(source)==record['source_sha256']
        text=blobs[prefix+name+'-target-row.txt'].decode().splitlines()
        n,m,_=map(int,text[0].split())
        assert n==1 and len(text)==1+m
        target=[parse(line)for line in text[1:]]
        with source.open()as f:
            nn,mm,_=map(int,f.readline().split())
            assert mm==m
            if name=='product':f.readline();f.readline()
            for i in range(record['target_row']):
                for _ in range(m):f.readline()
            original=[parse(f.readline())for _ in range(m)]
        assert target==original
        witness_name='product-row10.witness'if name=='product'else'full27-row4.witness'
        witness=blobs[prefix+witness_name]
        assert hashlib.sha256(witness).hexdigest()==record['witness_sha256']
        lines=witness.decode().splitlines()
        assert int(lines[0])==m and len(lines)==m+1
        coefficients=[parse(line)for line in lines[1:]]
        image={}
        for column,coefficient in zip(target,coefficients):image=plus(image,times(column,coefficient))
        if name=='product':
            assert [(j,c)for j,c in enumerate(coefficients)if c]==[(11,{(0,0,0):Fraction(-14985,8)}),(13,{(0,0,0):Fraction(2997,8)})]
            expected=expr('a^3*(5*a*d-6*e-4)^2')
            replacement=11
            determinant=Fraction(-14985,8)
        else:
            assert [(j,c)for j,c in enumerate(coefficients)if c]==[(5,{(0,0,0):Fraction(-1,50)}),(30,{(0,0,0):Fraction(-1,50)})]
            expected=expr('(100*a+18*e^2-9*e-2)*(5*d-6*e-4)')
            replacement=5
            determinant=Fraction(-1,50)
        assert image==expected==parse(record['target_unit_sparse'])
        assert coefficients[replacement]=={(0,0,0):determinant} and determinant
        unitpath='work/nonsquare-product-alt-v62/mod-ell/chart-3-ad-units.txt'if name=='product'else'work/nonsquare-v62/cubic27-v6/P-units'
        units=[expr(line)for line in blobs[unitpath].decode().splitlines()if line]
        assert len(units)==6
        reconstructed={(0,0,0):Fraction(record['target_unit_scalar'])}
        for unit,power in zip(units,record['target_unit_exponents']):
            for _ in range(power):reconstructed=times(reconstructed,unit)
        assert reconstructed==image
        assert units[0]==expr('a') and units[3]==expr('100*a+18*e^2-9*e-2') and units[5]==expr('d')
        assert record['residual_shape']==[nn-1,m-1]
        results.append(dict(case=name,full_source_stream_sha_checked=True,archived_target_row_equals_full_source_row=True,
            exact_Fraction_constant_witness_identity=True,replaced_column=replacement,
            column_replacement_determinant=str(determinant),pivot_polynomial='a^3*J^2'if name=='product'else'C*J',
            complete_declared_six_unit_list_and_factorization_checked=True,
            residual_shape=record['residual_shape'],native_residual_term_count_receipt=record['residual_terms']))
    paths=[Path(__file__),directory/'manifest.json',directory/'STATUS.json',directory/'README.md',archive,
           ROOT/'research/notes/S61_BINARY_NONSQUARE_SHORT_UNIT_PIVOTS_20261103_v62.md',
           ROOT/'research/notes/S61_BINARY_NONSQUARE_SHORT_UNIT_PIVOTS_20261103_v62_REVIEW.md',
           ROOT/'code/q_smith_global_model_review_20261103_v63.py',
           *[ROOT/f for f in status['dependencies']]]
    report=dict(id='S61_BINARY_NONSQUARE_SHORT_UNIT_PIVOTS_20261103_v62_REVIEW',all_passed=True,
        all_archived_member_and_dependency_hashes_checked=True,cases=results,
        same_column_combination_applied_to_every_row_is_invertible_over_QQ=True,
        exact_one_pivot_rank_splitting_on_declared_known_unit_open_accepted=True,
        main_constructive_unit_pivots_lower_bound=31,mod_ell_constructive_unit_pivots_lower_bound=11,
        native_large_residual_recomputed_in_this_review=False,
        all_point_rank41_or57_free_cokernel_quartic65_or_ambient_exhaustion_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
    with args.output.open('x')as f:
        json.dump(report,f,indent=2)
        f.write('\n')
    print('PASS V9: two independent Fraction identities and invertible known-unit one-pivot reductions')


if __name__=='__main__':main()

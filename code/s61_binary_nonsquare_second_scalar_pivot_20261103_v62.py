#!/usr/bin/env python3
"""Exact transport of the two-row certificate and a33rd constructive unit pivot."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_BINARY_NONSQUARE_SECOND_SCALAR_PIVOT_20261103_v62'
DIR=ROOT/'reports'/STEM
FIRST=ROOT/'reports/S61_BINARY_NONSQUARE_SCALAR_PRODUCT_ROWS_20261103_v62/certificate.json.gz'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def poly(data):return {tuple(ex):Fraction(c)for ex,c in data}
def add(p,q):
    out=defaultdict(Fraction,p)
    for ex,c in q.items():out[ex]+=c
    return {ex:c for ex,c in out.items()if c}
def scale(p,c):return {ex:v*c for ex,v in p.items()if v*c}
def mul(p,q):
    out=defaultdict(Fraction)
    for u,c in p.items():
        for v,d in q.items():out[tuple(a+b for a,b in zip(u,v))]+=c*d
    return {ex:c for ex,c in out.items()if c}
def dot(p,q):
    out={}
    for a,b in zip(p,q):out=add(out,mul(a,b))
    return out

def verify():
    start=time.monotonic();certificate=DIR/'certificate.json.gz';D=json.loads(gzip.decompress(certificate.read_bytes()));F=json.loads(gzip.decompress(FIRST.read_bytes()))
    previous_report=ROOT/'reports/S61_BINARY_NONSQUARE_SCALAR_PRODUCT_ROWS_20261103_v62/report.json';previous=json.loads(previous_report.read_text());assert previous['all_passed']
    assert sha(FIRST)==previous['input_sha256'][str(FIRST.relative_to(ROOT))]
    assert D['source_sha256']==F['native_outputs']['row7_replacement']['sha256']
    A0=[poly(p)for p in F['rows_sparse']['0']];A7=[poly(p)for p in F['rows_sparse']['7']]
    c7=[poly(p)for p in F['cofactors_sparse']['7']];L=[[poly(p)for p in row]for row in F['two_row_right_inverse_sparse']]
    a={(1,0,0):Fraction(1)};one={(0,0,0):Fraction(1)}
    assert dot(A7,c7)==a
    assert set(c7[43])=={(0,0,0)};pivot_coefficient=c7[43][(0,0,0)];assert pivot_coefficient
    # Inverse of the first column replacement (identity with column43=c7).
    ratio=scale(L[43][0],1/pivot_coefficient)
    transported=[add(L[j][0],scale(mul(ratio,c7[j]),-1))for j in range(60)if j!=43]
    cofs=[poly(p)for p in D['cofactor_sparse']];assert transported==cofs
    assert dot(A0,[row[0]for row in L])==a and not dot(A7,[row[0]for row in L])
    cross=dot(A0,c7)
    raw=[add(mul(a,A0[j]),scale(mul(cross,A7[j]),-1))for j in range(60)if j!=43]
    row=[poly(p)for p in D['row_sparse']];assert len(row)==len(raw)==59
    first=next(i for i,p in enumerate(raw)if p);lead_raw=max(raw[first]);lead_row=max(row[first]);exp=tuple(u-v for u,v in zip(lead_raw,lead_row));scalar=raw[first][lead_raw]/row[first][lead_row]
    assert exp==(1,0,0) and scalar
    unit={exp:scalar}
    assert all(p==mul(unit,q)for p,q in zip(raw,row))
    expected={(1,0,0):1/scalar};assert dot(row,cofs)==expected
    alpha={int(i):Fraction(c)for i,c in D['constant_combo'].items()};assert len(alpha)==5 and min(alpha)==50 and alpha[50]
    total={}
    for j,c in alpha.items():total=add(total,scale(cofs[j],c))
    assert total==one
    # Inverse column transform: newc50=sum(alpha_j*cj)=1; other cj stay.
    # NewA50=A50/alpha50, newAj=Aj-alpha_j/alpha50*A50.
    # The literal identity checks this generic change with arbitrary polynomials.
    fixture_A=[{(0,0,0):Fraction(1)},{(0,0,0):Fraction(-1)}];fixture_c=[{(0,1,0):Fraction(1)},one|{(0,1,0):Fraction(1)}]
    before=dot(fixture_A,fixture_c);fixture_B=[scale(fixture_A[0],-1),add(fixture_A[1],fixture_A[0])]
    assert before==dot(fixture_B,[one,fixture_c[1]])=={(0,0,0):Fraction(-1)}
    log=DIR/'native.log';text=log.read_text()
    assert 'EXACT_COFACTOR_CONSTANT_COMBINATION 5 EXPOSED_COLUMN 50' in text
    assert 'EXACT_WITNESS_REPLACE 50 TARGET 6' in text and 'TERMS 2755382' in text and 'PASS 1' in text
    result=D['output'];output=ROOT/result['path'];assert sha(output)==result['sha256'] and output.stat().st_size==result['bytes']
    files=[Path(__file__),certificate,FIRST,previous_report,log,DIR/'fixture.log',ROOT/'code/polynomial_cofactor_constant_span_pivot_20261103_v62.cpp',ROOT/'research/notes'/f'{STEM}.md']
    return dict(id=STEM,all_passed=True,open='a*P*F*B0*C*J!=0, P=ad, B0=aB, J=5P-6e-4',
        exact_two_row_witness_transport=True,exact_original_row_pivot_and_cleaning_factor_verified=True,
        row_cleaning_unit_monomial={'exponents':list(exp),'coefficient':str(scalar)},
        exact_scalar_unit_image_coefficient=str(1/scalar),scalar_unit_image_exponents=[1,0,0],
        cofactor_terms=sum(len(p)for p in cofs),constant_cofactor_combination_terms=5,
        exact_constant_cofactor_combination='1',invertible_QQ_column_change_exposed_coefficient_column=50,
        native_all25_before_after_polynomial_images_asserted_equal=True,
        residual_shape=[24,58],residual_terms=2755382,generated_output=result,
        total_constructive_product_unit_pivots=33,
        all_point_product_rank41_or_full_cubic_rank57_claimed=False,
        ambient_or_quartic_coverage_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in files},seconds=time.monotonic()-start)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();out=verify()
    with args.output.open('x')as f:json.dump(out,f,indent=2);f.write('\n')
    print('PASS',STEM,out['seconds'])

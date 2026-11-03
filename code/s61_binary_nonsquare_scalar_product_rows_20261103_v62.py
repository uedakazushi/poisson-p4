#!/usr/bin/env python3
"""Exact scalar unit lifts and a two-row split block of the primary cubic residual.

All accepted coefficient identities use Python Fraction arithmetic. Discovery
GBs and bounded modular searches do not substitute for these QQ identities.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_BINARY_NONSQUARE_SCALAR_PRODUCT_ROWS_20261103_v62'
DIR=ROOT/'reports'/STEM

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def polynomial(data):
    p={tuple(ex):Fraction(c)for ex,c in data}
    assert len(p)==len(data) and all(c and min(ex)>=0 for ex,c in p.items())
    return p

def add_product(out,p,q):
    for u,c in p.items():
        for v,d in q.items():out[tuple(a+b for a,b in zip(u,v))]+=c*d

def dot(A,B):
    out=defaultdict(Fraction)
    for p,q in zip(A,B):add_product(out,p,q)
    return {ex:c for ex,c in out.items()if c}

def verify():
    start=time.monotonic();cert=DIR/'certificate.json.gz';D=json.loads(gzip.decompress(cert.read_bytes()))
    status=ROOT/'research/handoff/20261103_v62/v9_short_unit_pivots/STATUS.json'
    V9=json.loads(status.read_text());main=next(c for c in V9['cases']if c['case']=='main')
    assert D['source_sha256']==main['residual_sha256']=='6d21ebfbeade3c5a0a7bd2dc091bde5ec274622ddf899add905242164d616ea3'
    assert D['source_shape']==[26,60] and V9['main_total_constructive_unit_pivots']==31
    source=ROOT/D['source'];assert source.exists(),'Restore immutable V6 then V9 main replay to reconstruct the pinned26x60 source'
    assert sha(source)==D['source_sha256']
    # Read every original sparse entry, comparing the two certificate rows.
    it=iter(source.read_text().split());n,m,k=[int(next(it))for _ in range(3)];assert(n,m,k)==(26,60,1)
    for _ in range(n+k):next(it)
    rows={}
    for i in range(n):
        row=[]
        for _ in range(m):
            p={}
            for _ in range(int(next(it))):
                c=Fraction(next(it));ex=tuple(int(next(it))for _ in range(3));assert ex not in p;p[ex]=c
            row.append(p)
        if str(i) in D['rows_sparse']:
            assert row==[polynomial(p)for p in D['rows_sparse'][str(i)]];rows[i]=row
    assert next(it,None)is None
    unit={(1,0,0):Fraction(1)};stats={}
    cof={int(i):[polynomial(p)for p in row]for i,row in D['cofactors_sparse'].items()}
    for i in (0,7):
        assert dot(rows[i],cof[i])==unit
        terms=sum(len(p)for p in cof[i]);degree=max(sum(ex)for p in cof[i]for ex in p)
        bits=max(max(abs(c.numerator).bit_length(),c.denominator.bit_length())for p in cof[i]for c in p.values())
        stats[str(i)]=dict(terms=terms,cofactor_degree=degree,max_coefficient_bits=bits)
    assert stats['0']['terms']==64 and stats['7']['terms']==34
    assert set(cof[7][43])=={(0,0,0)} and cof[7][43][(0,0,0)]
    L=[[polynomial(p)for p in row]for row in D['two_row_right_inverse_sparse']]
    assert len(L)==60 and all(len(row)==2 for row in L)
    for i,r in enumerate((0,7)):
        for j in range(2):assert dot(rows[r],[row[j]for row in L])==(unit if i==j else {})
    log=DIR/'full-row7-pivot.log';text=log.read_text()
    assert 'EXACT_WITNESS_REPLACE 43 TARGET 7' in text and 'TERMS 1801170' in text and 'PASS 1' in text
    generated=D['native_outputs']['row7_replacement'];residual=ROOT/generated['path']
    assert residual.exists() and sha(residual)==generated['sha256'] and residual.stat().st_size==generated['bytes']
    files=[Path(__file__),cert,status,log,ROOT/'research/notes'/f'{STEM}.md']
    return dict(id=STEM,all_passed=True,core_coordinates=['a','e','P=ad'],
        open='a*P*F*B0*C*J != 0, J=5P-6e-4, B0=aB',
        source=str(source.relative_to(ROOT)),source_sha256=D['source_sha256'],source_shape=[26,60],
        exact_scalar_identities='source_row0*cofactor0 = source_row7*cofactor7 = a',scalar_witnesses=stats,
        exact_two_original_rows_times_right_inverse='a*identity2',two_row_right_inverse_terms=sum(len(p)for row in L for p in row),
        all_four_two_row_QQ_identities_verified=True,
        preferred_native_operation='column43 replacement using row7 cofactor, followed by one unit pivot',
        residual_shape=[25,59],residual_terms=1801170,generated_residual=generated,
        total_constructive_product_unit_pivots=32,
        all_point_rank41_or_full_cubic_rank57_claimed=False,
        ambient_or_quartic_coverage_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in files},seconds=time.monotonic()-start)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();out=verify()
    with args.output.open('x')as f:json.dump(out,f,indent=2);f.write('\n')
    print('PASS',STEM,out['seconds'])

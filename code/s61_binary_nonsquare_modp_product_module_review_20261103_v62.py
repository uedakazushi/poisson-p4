#!/usr/bin/env python3
"""Independent streaming QQ-prefix binding and all-coefficient GF right-inverse audit."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_BINARY_NONSQUARE_MODP_PRODUCT_MODULE_20261103_v62_REVIEW'
P=32003


def sha(path):
    h=hashlib.sha256()
    with path.open('rb')as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def native(path):
    with (gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:
        def tokens():
            for line in f:
                yield from line.split()
        it=iter(tokens());n,m,k=[int(next(it))for _ in range(3)];assert n>0 and m>0 and k==0
        answer=[]
        for i in range(n):
            row=[]
            for j in range(m):
                count=int(next(it));assert count>=0;poly={}
                for _ in range(count):
                    c=int(next(it))%P;ex=tuple(int(next(it))for _ in range(3))
                    assert c and min(ex)>=0 and ex not in poly;poly[ex]=c
                row.append(poly)
            answer.append(row)
        assert next(it,None)is None
    return answer


def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports/S61_BINARY_NONSQUARE_MODP_PRODUCT_MODULE_20261103_v62'
    author=json.loads((folder/'report.json').read_text());info=json.loads((folder/'source.json').read_text())
    assert author['all_passed'] and info['prime']==P and s.isprime(P)
    checked={}
    for name,want in author['input_sha256'].items():
        actual=sha(ROOT/name);assert actual==want,name;checked[name]=actual
    source=ROOT/info['source'];assert sha(source)==info['source_sha256'];checked[info['source']]=info['source_sha256']
    keep=info['retained_columns'];assert len(keep)==len(set(keep))==53 and keep==sorted(keep)
    assert info['exact_source_rows']==list(range(10))
    M=native(folder/'source_modp.txt.gz');L=native(folder/'right_inverse_modp.txt.gz')
    assert (len(M),len(M[0]),len(L),len(L[0]))==(10,53,53,10)
    locations={col:j for j,col in enumerate(keep)};compared=0;qq_terms=0;zero_columns=[]
    with source.open()as f:
        assert tuple(map(int,f.readline().split()))==(26,60,1)
        metadata=[f.readline().strip(),f.readline().strip()]
        assert metadata==['5','1 2 3 0 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26']
        for i in range(10):
            for j in range(60):
                line=f.readline().split();count=int(line[0]);assert len(line)==1+4*count
                poly={};seen=set()
                for k in range(count):
                    coefficient=Fraction(line[1+4*k]);assert coefficient.denominator%P
                    ex=tuple(map(int,line[2+4*k:5+4*k]));assert min(ex)>=0 and ex not in seen;seen.add(ex)
                    value=coefficient.numerator*pow(coefficient.denominator,-1,P)%P
                    if value:poly[ex]=value
                qq_terms+=count
                if j in locations:assert poly==M[i][locations[j]];compared+=1
                else:assert count==0
    zero_columns=[j for j in range(60)if j not in locations]
    assert len(zero_columns)==7 and compared==530
    assert sum(len(poly)for row in M for poly in row)==author['source_terms']==205156
    assert sum(len(poly)for row in L for poly in row)==author['right_inverse_terms']==183182
    assert max(sum(ex)for row in L for poly in row for ex in poly)==author['right_inverse_max_total_degree']==18
    fixture=native(folder/'native_link_fixture.txt')
    assert fixture==[[{(1,0,0):1},{(0,1,0):1},{(0,0,1):1}],[{(2,0,0):1},{(0,2,0):1},{(0,0,2):1}]]
    try:native(folder/'INVALID_first_row_only.txt.gz')
    except (ValueError,AssertionError,StopIteration):invalid_rejected=True
    else:invalid_rejected=False
    assert invalid_rejected
    failed=(folder/'FAILED_metadata_and_first_row_export.log').read_text()
    assert 'identifier `maxdeg` in use' in failed
    # Decompress with streaming copies; the QQ94MB source is never read_bytes().
    for name,short in [('source_modp.txt.gz','M.txt'),('right_inverse_modp.txt.gz','L.txt')]:
        with gzip.open(folder/name,'rb')as src,(out/short).open('wb')as dest:shutil.copyfileobj(src,dest,1<<20)
    cpp=ROOT/'code/polynomial_modp_right_inverse_review_20261103_v62.cpp';binary=out/'check'
    flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True))
    subprocess.run(['c++','-std=c++17','-O2',str(cpp),'-o',str(binary),*flags],check=True,capture_output=True,text=True)
    result=subprocess.run([str(binary),str(out/'M.txt'),str(out/'L.txt')],check=True,capture_output=True,text=True,timeout=30)
    assert result.stdout.strip()=='INDEPENDENT_OUTER_PRODUCT_ALL100_IDENTITIES_PASS 100'
    (out/'native.log').write_text(result.stdout+result.stderr)
    binary.unlink();(out/'M.txt').unlink();(out/'L.txt').unlink()
    assert not author['original_GF_product41_claimed']
    assert not author['original_prefix_prime_compatibility_audited']
    assert not author['QQ_surjectivity_or_product41_or_full57_claimed']
    note=ROOT/f'research/notes/{STEM}.md'
    paths=[Path(__file__),cpp,note,folder/'report.json',folder/'source.json',source]
    paths += [folder/name for name in ['source_modp.txt.gz','right_inverse_modp.txt.gz','INVALID_first_row_only.txt.gz','FAILED_metadata_and_first_row_export.log','native_link_fixture.txt']]
    report=dict(id=STEM,all_passed=True,field='GF(32003)',
        independent_streaming_source_hash_bytes=source.stat().st_size,
        current_input_hashes_checked=len(checked),
        independently_bound_QQ_prefix_rows=list(range(10)),all530_saved_source_entries_exact_reductions=True,
        all_QQ_prefix_denominators_prime_units=True,omitted7_columns_zero_in_this10row_prefix=zero_columns,
        full_source_shape=[10,53],full_right_inverse_shape=[53,10],source_terms=205156,right_inverse_terms=183182,right_inverse_max_total_degree=18,
        independent_different_monomial_order_and_outer_product_native_all100_polynomial_identities=True,
        exact_identity='M*L=a*P^3*(100a+18e^2-9e-2)*(5P-6e-4)^2*I10',
        no_Singular_or_GB_or_CRT_run=True,
        failed_partial_export_rejected_and_retained=True,explicit_per_entry_native_fixture_checked=True,
        residual_full_image_on='D(a*P*C*J), hence on the stated a*P*F*B0*C*J open',
        all31_original_QQ_prefix_pivots_prime_compatibility_audited=False,
        original_GF_product41_claimed=False,QQ_product41_or_surjectivity_claimed=False,
        full_cubic57_quartic65_or_ambient_removal_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={'native.log':sha(out/'native.log')},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('INDEPENDENT GF RESIDUAL FULL100 IDENTITIES / SCOPE REVIEW PASS',report['seconds'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

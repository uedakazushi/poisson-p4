#!/usr/bin/env python3
"""Exact GF32003 right-inverse certificate for the pinned QQ cubic residual.

No assertion about original characteristic-p products or characteristic0
surjectivity follows without the separate bridges described in the note.
"""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import tempfile
import time
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_BINARY_NONSQUARE_MODP_PRODUCT_MODULE_20261103_v62';DIR=ROOT/'reports'/STEM

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read_matrix(data,p):
    it=iter(data.decode().split());n,m,k=[int(next(it))for _ in range(3)];assert k==0;A=[]
    for i in range(n):
        row=[]
        for j in range(m):
            q={}
            for z in range(int(next(it))):
                c=int(next(it))%p;ex=tuple(int(next(it))for _ in range(3));assert c and ex not in q;q[ex]=c
            row.append(q)
        A.append(row)
    assert next(it,None)is None
    return A

def verify():
    t=time.monotonic();p=32003;info=json.loads((DIR/'source.json').read_text());source=ROOT/info['source'];assert sha(source)==info['source_sha256']
    keep=info['retained_columns'];assert len(keep)==53 and info['exact_source_rows']==list(range(10))
    Mbytes=gzip.decompress((DIR/'source_modp.txt.gz').read_bytes());Lbytes=gzip.decompress((DIR/'right_inverse_modp.txt.gz').read_bytes());M=read_matrix(Mbytes,p);L=read_matrix(Lbytes,p);assert len(M)==10 and len(M[0])==53 and len(L)==53 and len(L[0])==10
    # Independently compare every saved source entry with direct reduction of
    # the pinned ORIGINAL QQ residual. Only identically zero columns are omitted.
    with source.open()as f:
        n,m,k=map(int,f.readline().split());assert(n,m,k)==(26,60,1);f.readline();f.readline()
        for i in range(10):
            for j in range(m):
                tokens=f.readline().split();nt=int(tokens[0]);assert len(tokens)==1+4*nt;poly={}
                for z in range(nt):
                    c=Fraction(tokens[1+4*z]);assert c.denominator%p;ex=tuple(int(v)for v in tokens[2+4*z:5+4*z]);val=c.numerator*pow(c.denominator,-1,p)%p
                    if val:poly[ex]=val
                if j in keep:assert poly==M[i][keep.index(j)]
                else:assert nt==0
    cpp=ROOT/'code/polynomial_modp_right_inverse_check_20261103_v62.cpp'
    with tempfile.TemporaryDirectory(prefix='nonsquare-modp-')as tmp:
        w=Path(tmp);mp=w/'M.txt';lp=w/'L.txt';binary=w/'check';mp.write_bytes(Mbytes);lp.write_bytes(Lbytes)
        flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));subprocess.run(['c++','-O2','-std=c++17',str(cpp),'-o',str(binary),*flags],check=True,capture_output=True,text=True)
        result=subprocess.run([str(binary),str(mp),str(lp),str(p)],check=True,capture_output=True,text=True);assert 'EXACT_ALL_FP_RIGHT_INVERSE_IDENTITIES_PASS 100' in result.stdout
    files=[Path(__file__),cpp,*[DIR/s for s in ['source.json','source_modp.txt.gz','right_inverse_modp.txt.gz','discovery.sing.gz','right_inverse_export.sing.gz','right_inverse_export.log','FAILED_metadata_and_first_row_export.log','INVALID_first_row_only.txt.gz','native_link_fixture.sing','native_link_fixture.txt']],ROOT/'research/notes'/f'{STEM}.md']
    return dict(id=STEM,all_passed=True,field='GF(32003)',source_QQ_residual_sha256=info['source_sha256'],
        source_shape=[10,53],right_inverse_shape=[53,10],source_terms=sum(len(q)for row in M for q in row),right_inverse_terms=sum(len(q)for row in L for q in row),
        right_inverse_max_total_degree=max(sum(ex)for row in L for q in row for ex in q),
        exact_all100_polynomial_identities='M*L=a*P^3*C*J^2*identity10',independent_checker='FLINT nmod_mpoly, separate from Singular discovery',
        residual_surjective_on='a*P*F*B0*C*J!=0, P=ad, B0=aB, J=5P-6e-4',
        original_GF_product41_claimed=False,original_prefix_prime_compatibility_audited=False,
        QQ_surjectivity_or_product41_or_full57_claimed=False,
        failed_partial_export_explicitly_preserved=True,
        input_sha256={str(q.relative_to(ROOT)):sha(q)for q in files},seconds=time.monotonic()-t)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();out=verify()
    with args.output.open('x')as f:json.dump(out,f,indent=2);f.write('\n')
    print('PASS',STEM,out['seconds'])

#!/usr/bin/env python3
"""Independently rebuild all AB quartic evaluation minors and check determinants."""
import argparse,hashlib,json,struct,subprocess,tempfile,time
from fractions import Fraction
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.time()
    source='reports/AB_ORIGIN_DEGREE_FOUR_20261101.json';r=json.loads((ROOT/source).read_text())
    c=r['evaluation_certificate'];prime=32003;assert c['prime']==prime
    points=np.array(c['points'],dtype=np.int64);assert points.shape==(726,10)
    weights=r['variable_weights'];coordinate_weights=r['coordinate_weights']
    power_cache={}
    for f in r['coordinate_polynomials']:
        for ex,_ in f:
            for j,n in enumerate(ex):
                if (j,n)not in power_cache:
                    v=np.ones(726,dtype=np.int64)
                    for k in range(n):v=v*points[:,j]%prime
                    power_cache[j,n]=v
    values=[]
    for i,f in enumerate(r['coordinate_polynomials']):
        v=np.zeros(726,dtype=np.int64)
        for ex,text in f:
            assert [sum(e*w[k]for e,w in zip(ex,weights))for k in range(2)]==coordinate_weights[i]
            co=Fraction(text);assert co.denominator%prime
            term=np.full(726,co.numerator%prime*pow(co.denominator,-1,prime)%prime,dtype=np.int64)
            for j,n in enumerate(ex):term=term*power_cache[j,n]%prime
            v=(v+term)%prime
        values.append(v)
    V=np.array(values,dtype=np.int64).T;total=0;seen=set()
    native='code/prime_square_minor_verifier_20261101.cpp'
    with tempfile.TemporaryDirectory(prefix='ab-minor-independent-')as temp:
        wd=Path(temp);inp=wd/'matrices.bin';binary=wd/'verify'
        with inp.open('wb')as out:
            out.write(struct.pack('<II',prime,len(c['groups'])))
            for g in c['groups']:
                w=tuple(g['weight']);assert w not in seen;seen.add(w)
                n=g['rank'];products=np.array(g['selected_products'],dtype=np.int64);assert products.shape==(n,4)
                assert len(set(map(tuple,products.tolist())))==n
                for inds in products:
                    assert tuple(sum(coordinate_weights[i][k]for i in inds)for k in range(2))==w
                M=np.ones((n,n),dtype=np.int64)
                for k in range(4):M=M*V[:n,products[:,k]]%prime
                out.write(struct.pack('<II',n,g['determinant']));out.write(M.astype('<u4').tobytes(order='C'))
                total+=n
        assert total==38570
        matrix_sha=hashlib.sha256(inp.read_bytes()).hexdigest()
        subprocess.run(['c++','-O3','-std=c++17',str(ROOT/native),'-o',str(binary)],check=True)
        proc=subprocess.run([str(binary),str(inp)],text=True,capture_output=True,check=True)
        result=json.loads(proc.stdout);assert result==dict(all_passed=True,minors=249,total_rank=38570)
    paths=[str(Path(__file__).relative_to(ROOT)),native,source]
    result.update(seconds=time.time()-started,evaluation_matrix_sha256=matrix_sha,
                  scope='Independent Python/NumPy construction of every saved minor and elementary Gaussian determinants; no FLINT matrix code.',
                  limitations=['Uses the primary saved coordinate polynomials; not a second derivation of the incidence map.',
                               'Does not independently prove the cohomological upper bound or all-degree normality.'],
                  proof_note='research/notes/AB_ORIGIN_DEGREE_FOUR_MINORS_REVIEW_20261101.md',
                  input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths})
    with args.output.open('x')as out:json.dump(result,out,indent=2);out.write('\n')
    print('ALL 249 INDEPENDENT MINORS PASS',result['seconds'],flush=True)

if __name__=='__main__':main()

#!/usr/bin/env python3
"""Exact weighted-fourth surjectivity over QQ(e)[d,D^-1].

Candidate denominators and polynomial kernels are verified, not assumed.
Rebuilds all large matrices and uses only 229 QQ chain samples because the
small denominator candidate is supplied. Typical run: about fifteen minutes.
No assertion on all vertical e-values or the full S61 fiber.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile
import time
import numpy as np
import sympy as s
from s61_binary_weighted_seventeen_curves_20261102 import prepare
ROOT=Path(__file__).resolve().parents[1]
NOTE=ROOT/'research/notes/S61_BINARY_WEIGHTED_HORIZONTAL_SUPPORT_20261102.md'
DATA=ROOT/'code/data/s61_binary_weighted_horizontal_support_20261102.tar.gz'
NATIVE=['delta_matrix_export','delta_kernel_peel','delta_content_reduce','delta_transpose_reduce',
        'multikernel_samples','multikernel_interpolation','chain_rational_samples','chain_interpolation',
        'linear_factor_shift','tropical_valuation_normalize','known_factor_valuation',
        'full_modular_rowmodule','determinant_degree_probe','verify_kernels']

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def independent_rank(paths,prime=32003):
    rows=[]
    for path in paths:
        with Path(path).open() as f:
            n,degree,count=map(int,f.readline().split());assert n==269
            for _ in range(count):
                data=f.readline().split();assert len(data)==1+2*int(data[0])
                row=np.zeros(n,dtype=np.int64)
                for j,c in zip(data[1::2],data[2::2]):
                    j=int(j);row[j%n]=(row[j%n]+int(c)%prime*pow(2,j//n,prime))%prime
                rows.append(row)
            assert not f.read().strip()
    basis={};pivots=[]
    for i,row in enumerate(rows):
        for j,b in basis.items():row=(row-int(row[j])*b)%prime
        nz=np.flatnonzero(row)
        if len(nz):
            j=int(nz[0]);row=row*pow(int(row[j]),-1,prime)%prime
            basis[j]=row;pivots.append([i,j])
    return len(basis),pivots

def verify(workdir):
    started=time.time();wd=Path(workdir).resolve();wd.mkdir(parents=True,exist_ok=True)
    dat=wd/'data';dat.mkdir()
    with tarfile.open(DATA) as tf:tf.extractall(dat,filter='data')
    manifest=json.loads((dat/'manifest.json').read_text())
    for name,item in manifest.items():assert sha(dat/name)==item['sha256']
    library=ROOT/'work/m2-runtime/lib';env=dict(os.environ,DYLD_LIBRARY_PATH=str(library)+':/opt/homebrew/lib')
    flags=['clang++','-O3','-std=c++17','-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),
           '-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),'-I/opt/homebrew/include']
    libs=[str(library/'libflint.24.0.dylib'),'-L/opt/homebrew/lib','-lgmp'];exe={};runs=[]
    build=tempfile.TemporaryDirectory(prefix='binary-horizontal-bin-');bindir=Path(build.name)
    sources=[ROOT/('code/s61_binary_weighted_'+name+'_20261102.cpp') for name in NATIVE]
    sources.extend([ROOT/'code/parametric_implicit_jets_20261101.cpp',ROOT/'code/s61_binary_deeper_weighted_extract_20261102.cpp'])
    for src in sources:
        out=bindir/src.stem
        compiled=subprocess.run(flags+[str(src)]+libs+['-o',str(out)],capture_output=True,text=True)
        if compiled.returncode:raise RuntimeError(compiled.stderr)
        exe[src.stem]=out
    def run(name,args,timeout=1200):
        key=name if name in exe else 's61_binary_weighted_'+name+'_20261102'
        log=wd/f'run{len(runs):02d}_{name}.log';begin=time.time()
        with log.open('x') as out:subprocess.run([str(exe[key]),*map(str,args)],env=env,stdout=out,stderr=out,check=True,timeout=timeout)
        runs.append(dict(name=name,seconds=time.time()-begin,log_sha256=sha(log),last_line=log.read_text().splitlines()[-1] if log.stat().st_size else ''))
        return log
    inp,raw,rows,meta,weights,coeff,known=prepare(wd)
    run('parametric_implicit_jets_20261101',[inp,raw])
    run('s61_binary_deeper_weighted_extract_20261102',[raw,coeff,rows,meta,weights])
    delta=wd/'delta';run('delta_matrix_export',[rows,2,delta]);assert list(map(int,delta.read_text().splitlines()[0].split()))[:3]==[2,32003,488]
    original=Path(str(delta)+'.full')
    proof_kernels=['poisson_weighted_leading_kernel_transferred46_20261102.txt',
                   'poisson_weighted_delta_iteration_20261102_s0.exact',
                   'poisson_weighted_delta_iteration_20261102_s1.batch',
                   'poisson_weighted_delta_left_from98_20261102_s0.exact',
                   'poisson_weighted_delta_left_from98_20261102_s1.exact']
    current=original;peels=[]
    for i,(kernel,dimension) in enumerate(zip(proof_kernels,[46,12,40,15,25])):
        output=wd/f'peel{i}';run('delta_kernel_peel',[current,dat/kernel,output]);assert f'PEELED {dimension}' in Path(str(output)+'.meta').read_text()
        peels.append(dimension);current=output
        if i==2:
            for operation in ['delta_content_reduce','delta_transpose_reduce']:
                output=wd/operation;run(operation,[current,output]);assert Path(str(output)+'.meta').read_text().rstrip().endswith('VALUATION 0');current=output
    output=wd/'stage138';run('delta_transpose_reduce',[current,output]);assert Path(str(output)+'.meta').read_text().rstrip().endswith('VALUATION 0');current=output
    assert sum(peels)==138
    pivot14=dat/'free-m0-binary-kernel-degree-probe.txt';sample14=wd/'kernel14';kernel14=wd/'kernel14_exact'
    run('multikernel_samples',[Path(str(current)+'.leading'),pivot14,100,sample14])
    run('multikernel_interpolation',[Path(str(current)+'.leading'),dat/'poisson_weighted_stage138_QQ_denominator_20261102.txt',sample14,97,14,kernel14])
    chain=wd/'chain';chain_result=wd/'chain_result'
    run('chain_rational_samples',[current,pivot14,229,chain])
    run('chain_interpolation',[current,dat/'poisson_chain_QQ_common_D_20261102.txt',chain,226,chain_result,pivot14])
    chaintext=chain_result.read_text();assert 'EXACT_POLYNOMIAL_CHAIN_IDENTITIES 2421' in chaintext and 'FOUR_INITIAL_VECTOR_RANK 4' in chaintext
    free_det=next(x.split(' ',1)[1] for x in chaintext.splitlines() if x.startswith('FOUR_INITIAL_FREE_BLOCK_DETERMINANT_E2 '));assert s.Rational(free_det)
    valuation=138+10+2*3+3;assert valuation==157 and 488-valuation==331
    factors=dat/'poisson_weighted_two_known_factors_reversed_20261102.txt'
    lower=wd/'known_lower';run('known_factor_valuation',[original,factors,lower]);assert lower.read_text().splitlines()[0]=='0 1'
    f=wd/'factor9_shift';run('linear_factor_shift',[original,factors,1,f]);orders=[]
    for operation,expected,label in [('delta_content_reduce',6,'r'),('tropical_valuation_normalize',44,'t'),('delta_content_reduce',0,'f')]:
        out=wd/('factor9_'+label);run(operation,[f,out]);text=Path(str(out)+'.meta').read_text();assert f'VALUATION {expected}' in text;orders.append(expected);f=out
    trans=wd/'factor9_lead_t';run('delta_transpose_reduce',[Path(str(f)+'.leading'),trans]);assert Path(str(trans)+'.meta').read_text().rstrip().endswith('VALUATION 0')
    lead=Path(str(trans)+'.leading');low33=dat/'poisson_factor9_kernel33_exact_20261102.txt';run('verify_kernels',[lead,low33])
    f9samples=wd/'factor9_samples';f9kernels=wd/'factor9_kernels';pivot10=dat/'free-m0-factor9-rational-complement-probe.txt'
    run('multikernel_samples',[lead,pivot10,123,f9samples])
    run('multikernel_interpolation',[lead,dat/'poisson_factor9_QQ_complement_D_20261102.txt',f9samples,120,10,f9kernels])
    rank43,pivots=independent_rank([low33,f9kernels]);assert rank43==43 and sum(orders)+rank43==93
    det=wd/'mod_det';run('determinant_degree_probe',[rows,2,det]);lines=det.read_text().splitlines();assert lines[1]==delta.read_text().splitlines()[1];assert list(map(int,lines[0].split()))==[2,32003,488,331]
    d=s.Symbol('d');v=list(map(int,lines[2].split()));P=s.Poly.from_list(v[1:][::-1],d,modulus=32003)
    h6=s.Poly(456*d+2160,d,modulus=32003);h9=s.Poly(178*d+800,d,modulus=32003)
    Q,rem=s.div(P,h6*h9**93);assert rem.is_zero and Q.degree()==237 and s.gcd(Q,h6*h9).degree()==0
    full=wd/'mod_full';run('full_modular_rowmodule',[rows,2,full]);v=list(map(int,full.read_text().splitlines()[1].split()));F=s.Poly.from_list(v[1:][::-1],d,modulus=32003);assert F.monic()==s.Poly(d+2162,d,modulus=32003)**92
    dependencies=[Path(__file__),NOTE,DATA,*sources,ROOT/'code/s61_binary_univariate_field_20261102.hpp',ROOT/'code/s61_binary_deeper_quartic_field_20261101.hpp',ROOT/'code/s61_binary_weighted_seventeen_curves_20261102.py',ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json',ROOT/'reports/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101.json',ROOT/'reports/S61_BINARY_DEEPER_WEIGHTED_FINITE_MODEL_20261102.json']
    return dict(id='S61_BINARY_WEIGHTED_HORIZONTAL_SUPPORT_20261102',all_passed=True,
        coefficient_ring='QQ(e)[d,D^-1]',matrix_shape=[390,269],row_map_surjective=True,
        exact_peels=peels,short_chain_lengths=[1]*10+[2]*3+[3],delta_order=157,
        tropical_d_degree=488,exact_d_degree=331,known_factor_multiplicities=[1,93],
        factor9_kernel_dimensions=[33,10],factor9_independence_pivots=pivots,
        polynomial_chain_identities=2421,free_four_minor_nonzero=True,
        residual_monic_degree=237,residual_coprime_to_known_factors_mod_prime=True,
        modular_full_row_span_determinant_factor='unit*(d+2162)^92',
        proof_method='Exact QQ polynomial identities, then finite-V CRT and Nakayama; see proof note.',
        limitations=['Generic e only; a finite unspecified set of vertical e-values remains.',
                    'Nuisance-zero weighted matrix only; no uniform all-nuisance same-point upper bound.',
                    'No whole-core saturation or affine S61 exhaustion.'],
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dependencies},
        rebuilt_sha256={'raw_weighted_rows':sha(rows),'original_delta_matrix':sha(original),'stage138_matrix':sha(current)},
        candidate_data_sha256=manifest,runs=runs,seconds=time.time()-started)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--workdir',type=Path)
    args=parser.parse_args();assert not args.output.exists()
    if args.workdir:result=verify(args.workdir)
    else:
        with tempfile.TemporaryDirectory(prefix='binary-horizontal-') as temp:result=verify(Path(temp))
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:result[k] for k in ['all_passed','coefficient_ring','exact_d_degree','known_factor_multiplicities','seconds']}))

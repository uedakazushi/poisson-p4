#!/usr/bin/env python3
"""Construct and independently check a factored global t-action.

The1031 symbols generate the actual full algebra over all17 base coordinates.
The free symbol module need not equal the actual quotient. The circuit gives
a lift T with phi*T=t*phi. A saved dense cofactor identity is independently
expanded here, without recomputing the native cubic/quartic inverses.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import time

import numpy as np

from s61_affine_seventeen_base_lift_20261103_v62 import build
from s61_affine_normal_projection_20261102 import ROOT,P


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ordered_monomials(n,d):
    return sorted(itertools.combinations_with_replacement(range(n),d),
                  key=lambda m:tuple(m.count(v) for v in range(n-1,-1,-1)))


def read_native(path):
    lines=path.read_text().splitlines();tokens=lines[0].split()
    summary={k:int(v) for k,v in zip(tokens[::2],tokens[1::2])}
    piv2=list(map(int,lines[1].split()[1:]));rows3=list(map(int,lines[2].split()[1:]));piv3=list(map(int,lines[3].split()[1:]))
    alternatives=list(map(int,lines[4].split()[1:]));assert len(alternatives)==2*summary['RESIDUAL4']
    trials=[];cursor=5
    for j in range(summary['TRIALS']):
        assert lines[cursor]==f'TRIAL {j} COEFFICIENT_IDENTITY_PASS 1'
        y,u,v=[np.array(list(map(int,line.split()[1:])),dtype=np.int64) for line in lines[cursor+1:cursor+4]]
        assert y.shape==(17,) and u.shape==v.shape==(1031,)
        trials.append((y,u,v));cursor+=4
    assert lines[cursor]==f'SAVED_COFACTOR_TRIAL {len(trials)-1}';cursor+=1
    H=np.zeros((134,435),dtype=np.int64);L=np.zeros((134,29),dtype=np.int64);D=np.zeros(134,dtype=np.int64)
    pairs=list(itertools.combinations_with_replacement(range(29),2));pi={m:i for i,m in enumerate(pairs)}
    for line in lines[cursor:]:
        kind,*a=line.split();a=list(map(int,a))
        if kind=='D':r,c=a;assert D[r]==0;D[r]=c
        elif kind=='L':r,k,c=a;assert L[r,k]==0;L[r,k]=c
        else:
            assert kind=='H';r,i,j,c=a;assert H[r,pi[i,j]]==0;H[r,pi[i,j]]=c
        assert 0<c<P
    return summary,piv2,rows3,piv3,alternatives,trials,H,L,D


def check(high,cross,base,bpairs,piv2,piv3,trials,H,L,D,known_points):
    n=29;pairs=list(itertools.combinations_with_replacement(range(n),2));cubics=ordered_monomials(n,3)
    free2=[m for i,m in enumerate(pairs) if i not in set(piv2)]
    free3=[m for i,m in enumerate(cubics) if i not in set(piv3)]
    basis=[()]+[(i,) for i in range(n)]+free2+free3;assert len(basis)==1031
    monomials=[list(itertools.combinations_with_replacement(range(n),d)) for d in range(5)]
    indices=[{m:i for i,m in enumerate(mons)} for mons in monomials]
    products23=np.array([[indices[3][tuple(sorted((*pair,k)))] for k in range(n)] for pair in pairs])
    products22=np.array([[indices[4][tuple(sorted((*a,*b)))] for b in pairs] for a in pairs])
    products11=np.array([[indices[2][tuple(sorted((i,j)))] for j in range(n)] for i in range(n)])
    Y,u,v=trials[-1];linear=cross@Y%P
    constant=base@np.array([int(Y[i])*int(Y[j])%P for i,j in bpairs],dtype=np.int64)%P
    expansion=[np.zeros(len(mons),dtype=np.int64) for mons in monomials]
    for r in range(134):
        np.add.at(expansion[4],products22.ravel(),(H[r,:,None]*high[r,None,:]).ravel())
        np.add.at(expansion[3],products23.ravel(),(H[r,:,None]*linear[r,None,:]+high[r,:,None]*L[r,None,:]).ravel())
        expansion[2]+=H[r]*constant[r]+D[r]*high[r]
        np.add.at(expansion[2],products11.ravel(),(L[r,:,None]*linear[r,None,:]).ravel())
        expansion[1]+=L[r]*constant[r]+D[r]*linear[r]
        expansion[0][0]+=D[r]*constant[r]
    expected=[np.zeros(len(mons),dtype=np.int64) for mons in monomials]
    for mon,a,b in zip(basis,u,v):
        target=tuple(sorted((0,*mon)))
        expected[len(target)][indices[len(target)][target]]+=a
        expected[len(mon)][indices[len(mon)][mon]]-=b
    assert all(np.array_equal(a%P,b%P) for a,b in zip(expansion,expected))
    for point in known_points:
        Z=np.array(point+[0]*25,dtype=np.int64);assert Z.shape==(29,)
        evaluation=np.array([np.prod(Z[list(mon)],dtype=np.int64)%P if mon else 1 for mon in basis],dtype=np.int64)
        _,source,result=trials[1]
        assert int(evaluation@(result-source))%P==0
    assert trials[0][2][1]==1 and np.count_nonzero(trials[0][2])==1
    return dict(independent_saved_original134_identity=True,
                independently_checked_unknown_monomials_by_degree=list(map(len,monomials)),
                cofactor_nonzero_terms=dict(quadratic=int(np.count_nonzero(H)),linear=int(np.count_nonzero(L)),constant=int(np.count_nonzero(D))),
                known_five_t_action_evaluation=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--timeout',type=int,default=240)
    parser.add_argument('--reuse-input',type=Path)
    parser.add_argument('--reuse-output',type=Path)
    parser.add_argument('--cache-output',type=Path)
    args=parser.parse_args();assert bool(args.reuse_input)==bool(args.reuse_output)
    args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    high,cross,base,pairs,cubics,quartics,bpairs,vector,inputs,_,_=build()
    path=args.output/'quadrics.in';native=args.output/'native.txt'
    path.write_text('29 134 32003 17\n'+'\n'.join(' '.join(map(str,a.ravel())) for a in [high,cross,base])+'\n'+' '.join(map(str,vector%P))+'\n')
    source=ROOT/'code/s61_affine_t_action_circuit_20261103_v62.cpp'
    if args.reuse_input:
        assert sha(path)==sha(args.reuse_input)
        shutil.copyfile(args.reuse_output,native)
        log=args.reuse_output.with_suffix('.log')
        if log.exists():shutil.copyfile(log,args.output/'native.log')
        mode='native_certificate_reused_with_input_hash_and_independent_full_cofactor_check'
    else:
        executable=args.output/'t-action'
        flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],capture_output=True,text=True,check=True).stdout)
        subprocess.run(['c++','-std=c++17','-O3',str(source),*flags,'-o',str(executable)],check=True)
        command=[str(executable.resolve()),str(path.resolve()),str(native.resolve()),'6']
        if args.cache_output:command.append(str(args.cache_output.resolve()))
        with (args.output/'native.log').open('w')as stream:subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=args.timeout)
        executable.unlink();mode='native_circuit_replayed'
    summary,piv2,rows3,piv3,alternatives,trials,H,L,D=read_native(native)
    assert summary['GENERATORS']==1031 and summary['BASE']==17 and summary['TRIALS']==6
    assert summary['DET3']==14008 and summary['DET4']==8797
    assert len(piv2)==134 and len(rows3)==len(piv3)==3795
    known_path=ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json'
    known=json.loads(known_path.read_text());assert known['all_passed']
    checks=check(high,cross,base,bpairs,piv2,piv3,trials,H,L,D,known['known_five_mod_32003'])
    inputs += [Path(__file__),source,ROOT/'code/s61_affine_seventeen_base_lift_20261103_v62.py',known_path,
               ROOT/'code/s61_affine_finite_normal_module_20261102.py',ROOT/'code/s61_affine_normal_projection_20261102.py',ROOT/'code/s61_affine_fiber_macaulay.py']
    report=dict(id='S61_AFFINE_T_ACTION_CIRCUIT_20261103_v62',all_passed=True,prime=P,mode=mode,
                parameter_count=17,symbol_count=1031,symbol_counts_by_unknown_degree=[1,29,301,700],
                full_original_quadrics=134,operator='T_t(Y): B^1031 -> B^1031, phi*T_t=t*phi',
                representation='factored polynomial circuit, linear in arbitrary1031 input coefficients',
                maximum_parameter_degree=4,trial_count=6,
                native_each_trial_original134_coefficient_identity=True,**checks,
                full_dense_matrix_materialized=False,actual_relation_kernel_computed=False,
                free_symbol_module_equals_actual_algebra_claimed=False,full_affine_exhaustion_claimed=False,
                limitations=['Explicit constants are frozen-beta F32003; finite nonzero inverses define the analogous genericQQ circuit.',
                             'Formal all17-parameter correctness follows from exact original-equation subtraction in the circuit; six numerical tests are separate verification.',
                             'The1031 symbols are generators. The free lift can have extra eigenvalues and support; its uncomputed kernel cannot be discarded.'],
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(inputs)},
                evidence_sha256={p.name:sha(p) for p in args.output.iterdir() if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',checks,report['seconds'])


if __name__=='__main__':main()

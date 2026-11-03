#!/usr/bin/env python3
"""Exact bounded fixed-pair degree3 probe, with independent dual certificates.

This checks only the frozen-beta F32003 degree3 Macaulay space. It neither
counts the full affine fiber nor bounds the generic two-coordinate rank.
"""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import time

import numpy as np
import sympy as sp

from s61_affine_finite_normal_module_20261102 import ROOT,P,slice_model,load_fiber,residue


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def model():
    section,Q,axes,inputs=slice_model()
    known_path=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json'
    known=json.loads(known_path.read_text());_,ids,base,_,_=load_fiber();names=sp.symbols('v0:3')
    basis=np.column_stack([np.array(section['fixed_normal_plane_integer_basis']),np.eye(42,dtype=np.int64)[:,axes]])
    indices=[ids.index(i)for i in [5,41]]
    slopes=sp.Matrix([[sp.sympify(known['affine_plane_substitution'][i]).coeff(z)for z in names]for i in indices])
    T=np.array([[residue(base[i])]+[residue(v)for v in slopes[j,:]]+list(basis[i-3])for j,i in enumerate(indices)],dtype=np.int64)%P
    free=[0]+list(range(3,46));inverse=np.array(sp.Matrix(T[:,1:3].tolist()).inv_mod(P)).astype(np.int64)
    substitution=np.zeros((46,44),dtype=np.int64);substitution[free]=np.eye(44,dtype=np.int64)
    substitution[1:3]=-inverse@T[:,free]%P;assert not np.any(T@substitution%P)
    full=substitution.T@Q@substitution%P
    pairs=list(itertools.combinations_with_replacement(range(44),2))
    high=np.column_stack([full[:,i,j]+(full[:,j,i]if i!=j else 0)for i,j in pairs])%P
    return high,pairs,substitution,inputs+[known_path]


def check_dual(path,high,pairs):
    cubes=list(itertools.combinations_with_replacement(range(44),3));ci={m:i for i,m in enumerate(cubes)}
    dual=np.zeros((len(cubes),44),dtype=np.int64)
    lines=path.read_text().splitlines();assert lines[0]=='N 44 DUAL_COLUMNS 44 PRIME 32003'
    for line in lines[1:]:
        a=list(map(int,line.split()));assert len(a)==47
        mon=tuple(a[:3]);assert tuple(sorted(mon))==mon
        row=ci[mon];assert not np.any(dual[row]);dual[row]=a[3:]
    assert np.all((dual>=0)&(dual<P))
    targets=[ci[(0,0,k)]for k in range(44)]
    assert np.array_equal(dual[targets],np.eye(44,dtype=np.int64))
    for k in range(44):
        indices=[ci[tuple(sorted((*pair,k)))]for pair in pairs]
        assert not np.any(high@dual[indices]%P)
    return int(np.count_nonzero(dual)),len(lines)-1


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--timeout',type=int,default=240)
    parser.add_argument('--reuse-input',type=Path)
    parser.add_argument('--reuse-output',type=Path)
    parser.add_argument('--reuse-dual',type=Path)
    args=parser.parse_args();assert bool(args.reuse_input)==bool(args.reuse_output)==bool(args.reuse_dual)
    args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    high,pairs,substitution,inputs=model()
    source=ROOT/'code/s61_affine_fixed_pair_cubic_probe_20261103_v62.cpp'
    path=args.output/'quadrics.in';native=args.output/'native.txt';dual=args.output/'dual-functions.txt'
    path.write_text('44 134 32003\n'+' '.join(map(str,high.ravel()))+'\n')
    if args.reuse_input:
        assert sha(path)==sha(args.reuse_input)
        shutil.copyfile(args.reuse_output,native);shutil.copyfile(args.reuse_dual,dual)
        log=args.reuse_output.with_suffix('.log')
        if log.exists():shutil.copyfile(log,args.output/'native.log')
        mode='same_input_native_certificate_reused_and_all44_duals_independently_checked'
    else:
        executable=args.output/'probe'
        flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],capture_output=True,text=True,check=True).stdout)
        subprocess.run(['c++','-std=c++17','-O3',str(source),*flags,'-o',str(executable)],check=True)
        with (args.output/'native.log').open('w')as stream:
            subprocess.run([str(executable.resolve()),str(path.resolve()),str(native.resolve()),str(dual.resolve())],stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=args.timeout)
        executable.unlink();mode='native_replayed_and_all44_duals_independently_checked'
    tokens=native.read_text().split();summary={k:int(v)for k,v in zip(tokens[::2],tokens[1::2])}
    assert summary==dict(N=44,ORIGINAL134_CUBIC_RANK=5805,CUBIC_FREE=9375,T2_LINEAR_IMAGE_RANK=44,T3_ZERO=0,DUAL_IDENTITIES_VERIFIED=1)
    terms,rows=check_dual(dual,high,pairs)
    with dual.open('rb')as original,gzip.open(args.output/'dual-functions.txt.gz','wb')as compressed:shutil.copyfileobj(original,compressed)
    dual.unlink();(args.output/'substitution.json').write_text(json.dumps(substitution.tolist())+'\n')
    inputs += [Path(__file__),source,ROOT/'code/s61_affine_finite_normal_module_20261102.py',
               ROOT/'code/s61_affine_normal_projection_20261102.py',ROOT/'code/s61_affine_fiber_macaulay.py']
    report=dict(id='S61_AFFINE_FIXED_PAIR_CUBIC_PROBE_20261103_v62',all_passed=True,prime=P,mode=mode,
                fixed_actual_coordinates=['u6','u42'],fixed_values=[0,0],homogeneous_variable_count=44,
                original_equations=134,cubic_matrix_shape=[5896,15180],cubic_row_rank=5805,
                t2_times_all44_linear_classes_rank=44,t_cubed_in_degree3_ideal=False,
                dual_functional_count=44,independent_original_cubic_rows_checked=5896,
                independent_dual_pairings_checked=44*44,dual_nonzero_terms=terms,dual_nonzero_rows=rows,
                full_affine_fiber_length_claimed=False,generic_finite_pair_rank_bound_claimed=False,
                limitations=['This exact frozen-beta GF degree3 probe excludes that particular low-degree rank29 shortcut only.',
                             'Higher degree or t-saturation can change the relevant image;44 is not a lower bound on the actual fiber length.',
                             'No genericQQ upper bound or full affine exhaustion follows.'],
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(inputs)},
                evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS',summary,terms,report['seconds'])


if __name__=='__main__':main()

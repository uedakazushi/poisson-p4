#!/usr/bin/env python3
"""Construct a full17-base homogeneous finite model and check a t^4 lift.

All134 original quadrics and all17 retained coordinates are preserved.
The explicit cofactor/tensor certificate is F32003 arithmetic at the frozen
beta; its nonzero finite matrix ranks give the generic QQ generation bound.
No full S61 exhaustion is claimed.
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

from s61_affine_finite_normal_module_20261102 import ROOT, slice_model
from s61_affine_normal_projection_20261102 import P, pivot_minor


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build():
    section,Q,axes,inputs=slice_model()
    border=ROOT/'reports/S61_AFFINE_BORDER_SUPPORT_20261102.json'
    vector=np.array(json.loads(border.read_text())['retained_sample_direction'][3:],dtype=np.int64)
    assert vector.shape==(17,) and np.any(vector%P)
    change=np.eye(46,dtype=np.int64)
    change[29:,0]=-vector
    full=change.T@Q@change%P
    n,nb=29,17
    pairs=list(itertools.combinations_with_replacement(range(n),2))
    cubics=list(itertools.combinations_with_replacement(range(n),3))
    quartics=list(itertools.combinations_with_replacement(range(n),4))
    base_pairs=list(itertools.combinations_with_replacement(range(nb),2))
    high=np.column_stack([full[:,i,j]+(full[:,j,i] if i!=j else 0) for i,j in pairs])%P
    cross=(full[:,:n,n:]+full[:,n:,:n].transpose(0,2,1))%P
    base=np.column_stack([full[:,n+i,n+j]+(full[:,n+j,n+i] if i!=j else 0) for i,j in base_pairs])%P
    rows,cols,det=pivot_minor(high)
    assert len(rows)==134
    return high,cross,base,pairs,cubics,quartics,base_pairs,vector,inputs+[border],det,cols


def verify_coefficients(high,cross,base,pairs,cubics,quartics,cofactor,target):
    n,nb=cross.shape[1:]
    ci={m:i for i,m in enumerate(cubics)}
    qi={m:i for i,m in enumerate(quartics)}
    products23=np.array([[ci[tuple(sorted((*pair,k)))] for k in range(n)] for pair in pairs])
    products22=np.array([[qi[tuple(sorted((*a,*b)))] for b in pairs] for a in pairs])
    expansion=np.zeros(len(quartics),dtype=np.int64)
    for r in range(134):
        np.add.at(expansion,products22.ravel(),(cofactor[r,:,None]*high[r,None,:]).ravel())
    expected=np.zeros(len(quartics),dtype=np.int64);expected[qi[(target,)*4]]=1
    assert np.array_equal(expansion%P,expected)
    linear=np.zeros((len(cubics),nb),dtype=np.int64)
    product=(cofactor.T@cross.reshape(134,-1)%P).reshape(len(pairs),n,nb)
    np.add.at(linear,products23.ravel(),product.reshape(-1,nb))
    linear%=P
    quadratic=cofactor.T@base%P
    # Independent direct accumulation, rather than using the matrix products above.
    direct_linear=np.zeros_like(linear);direct_quadratic=np.zeros_like(quadratic)
    for r in range(134):
        values=cofactor[r,:,None,None]*cross[r,None,:,:]
        np.add.at(direct_linear,products23.ravel(),values.reshape(-1,nb))
        direct_quadratic+=cofactor[r,:,None]*base[r,None,:]
    assert np.array_equal(linear,direct_linear%P)
    assert np.array_equal(quadratic,direct_quadratic%P)
    return linear,quadratic


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--timeout',type=int,default=180)
    parser.add_argument('--reuse-input',type=Path)
    parser.add_argument('--reuse-output',type=Path)
    args=parser.parse_args();assert bool(args.reuse_input)==bool(args.reuse_output)
    args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    high,cross,base,pairs,cubics,quartics,base_pairs,vector,inputs,det2,cols2=build()
    path=args.output/'quadrics.in'
    path.write_text('29 134 32003\n'+' '.join(map(str,high.ravel()))+'\n')
    native=args.output/'native.txt';source=ROOT/'code/s61_affine_quartic_lift_20261103_v62.cpp'
    mode='native_replayed'
    if args.reuse_input:
        assert sha(path)==sha(args.reuse_input)
        shutil.copyfile(args.reuse_output,native)
        saved_log=args.reuse_output.with_suffix('.log')
        if saved_log.exists(): shutil.copyfile(saved_log,args.output/'native.log')
        mode='native_certificate_reused_with_input_hash_match_and_independent_coefficient_verification'
    else:
        executable=args.output/'lift'
        flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],capture_output=True,text=True,check=True).stdout)
        subprocess.run(['c++','-std=c++17','-O3','-DS61_LIFT_N=29',str(source),*flags,'-o',str(executable)],check=True)
        with (args.output/'native.log').open('w')as stream:
            subprocess.run([str(executable.resolve()),str(path.resolve()),str(native.resolve()),'0'],stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=args.timeout)
        executable.unlink()
    lines=native.read_text().splitlines();tokens=lines[0].split();summary={k:int(v) for k,v in zip(tokens[::2],tokens[1::2])}
    assert summary['TARGET']==0 and summary['RANK4']==len(quartics)
    cofactor=np.zeros_like(high);index={pair:i for i,pair in enumerate(pairs)}
    for line in lines[4:]:
        r,i,j,value=map(int,line.split());assert 0<=r<134 and 0<=value<P
        k=index[i,j];assert cofactor[r,k]==0;cofactor[r,k]=value
    assert np.count_nonzero(cofactor)==summary['TERMS']
    linear,quadratic=verify_coefficients(high,cross,base,pairs,cubics,quartics,cofactor,0)
    l=[];q=[]
    for mon,b in zip(*np.nonzero(linear)):
        l.append([int(b),*cubics[mon],int(linear[mon,b])])
    for mon,b in zip(*np.nonzero(quadratic)):
        q.append([*base_pairs[b],*pairs[mon],int(quadratic[mon,b])])
    identity=dict(convention='t^4=sum(H_i*F_i_full)-sum(Y_j*C_j); in the full algebra t^4=-sum(Y_j*C_j)',
                  unknown_order=['t','plane1','plane2','plane3']+[f'w{i}' for i in range(1,23)]+['old_y1','old_y2','old_y3'],
                  base_order=[f'Y{i}' for i in range(1,18)],shift_vector=list(map(int,vector)),
                  original134_quadratic_cofactors=[[int(r),*pairs[k],int(cofactor[r,k])] for r,k in zip(*np.nonzero(cofactor))],
                  parameter_linear_unknown_cubic_terms=l,parameter_quadratic_unknown_quadratic_terms=q)
    with gzip.open(args.output/'full-parameter-identity.json.gz','wt')as stream:json.dump(identity,stream,separators=(',',':'))
    free2=len(pairs)-134;free3=summary['FREE3'];bound=1+29+free2+free3
    assert summary['RANK3']+free3==len(cubics)
    inputs += [Path(__file__),source,ROOT/'code/s61_affine_finite_normal_module_20261102.py',ROOT/'code/s61_affine_normal_projection_20261102.py',ROOT/'code/s61_affine_fiber_macaulay.py']
    report=dict(id='S61_AFFINE_SEVENTEEN_BASE_LIFT_20261103_v62',all_passed=True,prime=P,mode=mode,
                parameter_base_count=17,unknown_count=29,quadratic_rank=134,quadratic_minor=det2,quadratic_pivot_columns=cols2,
                cubic_rank=summary['RANK3'],cubic_minor=summary['DET3'],cubic_free_count=free3,
                full_quartic_rank=summary['RANK4'],quartic_overlap_count=(len(lines[3].split())-1)//2,quartic_overlap_minor=summary['DET4'],
                homogeneous_fiber_generator_counts=[1,29,free2,free3],generator_bound=bound,
                exact_t4_original134_identity=True,all17_base_coefficient_identity_verified=True,
                expanded_parameter_linear_terms=len(l),expanded_parameter_quadratic_terms=len(q),
                t_nilpotence_global_claimed=False,full_affine_exhaustion_claimed=False,
                limitations=['Explicit coefficients are F32003 at the frozen rational beta.',
                             'Nonzero finite-matrix ranks imply the generic QQ finite-generation bound.',
                             'The full t-action matrix is not materialized; this certifies an actual t4 reduction rule.'],
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(inputs)},
                evidence_sha256={p.name:sha(p) for p in args.output.iterdir() if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',bound,summary,len(l),len(q),report['seconds'])


if __name__=='__main__':
    main()

#!/usr/bin/env python3
"""Certify a nonempty cyclic open of the full seventeen-base finite algebra.

The actual global presentation has a monic degree1031 polynomial and all134
original residual equations; no symbol-freeness assumption is used. Numerical
degree5 gcd evidence concerns only Y=v. The determinant-zero locus is open work.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import time

import numpy as np

from s61_affine_seventeen_base_lift_20261103_v62 import build
from s61_affine_normal_projection_20261102 import ROOT,P


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def evaluate(coefficients,value):
    result=0
    for coefficient in reversed(coefficients):result=(result*int(value)+int(coefficient))%P
    return result


def check_certificate(path,known):
    lines=path.read_text().splitlines()
    first=lines[0].split();summary={k:int(v) for k,v in zip(first[::2],first[1::2])}
    ell=list(map(int,next(line for line in lines if line.startswith('LINEAR_FUNCTIONAL ')).split()[1:]))
    det=int(next(line for line in lines if line.startswith('KRYLOV_DETERMINANT ')).split()[1])
    degree=int(next(line for line in lines if line.startswith('ACTUAL_GCD_DEGREE ')).split()[1])
    gcd=list(map(int,next(line for line in lines if line.startswith('GCD ')).split()[1:]))
    chi=list(map(int,next(line for line in lines if line.startswith('CHI ')).split()[1:]))
    coordinates={int(line.split()[1]):list(map(int,line.split()[2:])) for line in lines if line.startswith('COORDINATE ')}
    assert summary['GENERATORS']==1031 and len(ell)==29 and det!=0
    assert degree==5 and len(gcd)==6 and gcd[-1]==1
    assert len(chi)==1032 and chi[-1]==1 and set(coordinates)==set(range(29))
    assert all(len(row)==1031 for row in coordinates.values())
    combined=np.array(ell,dtype=np.int64)@np.array([coordinates[i]for i in range(29)],dtype=np.int64)%P
    expected=np.zeros(1031,dtype=np.int64);expected[1]=1
    assert np.array_equal(combined,expected)
    roots=[]
    for point in known['known_five_mod_32003']:
        Z=point+[0]*25
        h=sum(int(c)*int(z) for c,z in zip(ell,Z))%P;roots.append(h)
        assert evaluate(chi,h)==evaluate(gcd,h)==0
        assert all(evaluate(coordinates[i],h)==Z[i] for i in range(29))
    assert len(set(roots))==5
    polynomial=[1]
    for root in roots:
        product=[0]*(len(polynomial)+1)
        for i,c in enumerate(polynomial):product[i]=(product[i]-root*c)%P;product[i+1]=(product[i+1]+c)%P
        polynomial=product
    assert polynomial==gcd
    # Save only the compact rank data, one full original cofactor identity,
    # and final univariate certificate; redundant iterative input/output is
    # reproducible scratch output of the native evaluator.
    prefix_end=next(i for i,line in enumerate(lines) if line.startswith('SAVED_COFACTOR_TRIAL '))+1
    while prefix_end<len(lines) and lines[prefix_end].startswith(('H ','L ','D ')):prefix_end+=1
    end_start=next(i for i,line in enumerate(lines) if line.startswith('KRYLOV_DETERMINANT '))
    compact='\n'.join(lines[:prefix_end]+lines[end_start:])+'\n'
    return summary,ell,det,roots,gcd,compact


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--timeout',type=int,default=900)
    parser.add_argument('--reuse-input',type=Path)
    parser.add_argument('--reuse-output',type=Path)
    args=parser.parse_args();assert bool(args.reuse_input)==bool(args.reuse_output)
    args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    high,cross,base,*rest=build();vector=rest[4];inputs=rest[5]
    path=args.output/'quadrics.in'
    path.write_text('29 134 32003 17\n'+'\n'.join(' '.join(map(str,a.ravel())) for a in [high,cross,base])+'\n'+' '.join(map(str,vector%P))+'\n')
    source=ROOT/'code/s61_affine_cyclic_open_20261103_v62.cpp'
    raw=args.output/'raw-native.txt'
    if args.reuse_input:
        assert sha(path)==sha(args.reuse_input)
        raw_path=args.reuse_output
        log=raw_path.with_suffix('.log')
        if log.exists():shutil.copyfile(log,args.output/'native.log')
        mode='native_certificate_reused_with_input_hash_and_independent_coordinate_checks'
    else:
        executable=args.output/'cyclic'
        flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],capture_output=True,text=True,check=True).stdout)
        subprocess.run(['c++','-std=c++17','-O3',str(source),*flags,'-o',str(executable)],check=True)
        with (args.output/'native.log').open('w')as stream:
            subprocess.run([str(executable.resolve()),str(path.resolve()),str(raw.resolve()),'3'],stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=args.timeout)
        executable.unlink();raw_path=raw;mode='native_cyclic_open_replayed'
    known_path=ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json'
    known=json.loads(known_path.read_text());assert known['all_passed']
    summary,ell,det,roots,gcd,compact=check_certificate(raw_path,known)
    raw_digest=sha(raw_path)
    (args.output/'native-certificate.txt').write_text(compact)
    if raw.exists():raw.unlink()
    inputs += [Path(__file__),source,ROOT/'code/s61_affine_seventeen_base_lift_20261103_v62.py',
               ROOT/'code/s61_affine_finite_normal_module_20261102.py',ROOT/'code/s61_affine_normal_projection_20261102.py',ROOT/'code/s61_affine_fiber_macaulay.py',known_path]
    report=dict(id='S61_AFFINE_CYCLIC_OPEN_20261103_v62',all_passed=True,prime=P,mode=mode,
                parameter_count=17,symbol_count=1031,linear_functional=ell,
                sample_base='Y=v',krylov_determinant_mod_prime=det,monic_polynomial_degree=1031,
                native_all134_original_remainders_gcd_degree=5,sample_gcd_coefficients=gcd,
                known_five_h_values=roots,independent_all29_coordinate_values_at_all5_points=True,
                independent_linear_functional_identity=True,
                actual_global_presentation='On detK!=0: A=B[s]/(chi1031(s), F_1(P(s),Y),...,F_134(P(s),Y)). Affine t=1 adds P_0(s)-1.',
                global_gcd_degree5_claimed=False,determinant_zero_locus_resolved=False,full_affine_exhaustion_claimed=False,
                limitations=['Explicit constants and univariate polys are frozen-beta F32003 atY=v.',
                             'The formal global circuit and two inverse maps give the actual presentation on the nonempty cyclic open over genericQQ beta.',
                             'The complementary determinant-zero locus and full17-parameter gcd support remain unresolved.'],
                raw_native_output_sha256=raw_digest,
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(inputs)},
                evidence_sha256={p.name:sha(p) for p in args.output.iterdir() if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',det,roots,gcd,report['seconds'])


if __name__=='__main__':main()

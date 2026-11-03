#!/usr/bin/env python3
"""Independently recompute all134 univariate residuals with Cloud M2.

Does not reconstruct the Krylov matrix. Uses the saved monic polynomial and
coordinate polynomials, then substitutes the actual original134 quadrics.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np

from s61_affine_seventeen_base_lift_20261103_v62 import build
from s61_affine_t_action_circuit_20261103_v62 import ordered_monomials
from s61_affine_normal_projection_20261102 import ROOT,P


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def literal(value):return json.dumps(value,separators=(',',':')).replace('[','{').replace(']','}')


def check_dense_action(lines,ell,high,linear,constant):
    """Direct30-variable homogeneous expansion; no operator reconstruction."""
    piv2=set(map(int,lines[1].split()[1:]));piv3=set(map(int,lines[3].split()[1:]))
    zpairs=list(itertools.combinations_with_replacement(range(29),2))
    gcubics=ordered_monomials(29,3)
    basis=[()]+[(i,)for i in range(29)]+[m for i,m in enumerate(zpairs)if i not in piv2]+[m for i,m in enumerate(gcubics)if i not in piv3]
    assert len(basis)==1031
    trial=next(i for i,line in enumerate(lines)if line.startswith('TRIAL '))
    u=list(map(int,lines[trial+2].split()[1:]));v=list(map(int,lines[trial+3].split()[1:]))
    pairs=list(itertools.combinations_with_replacement(range(30),2));pi={m:i for i,m in enumerate(pairs)}
    quartics=list(itertools.combinations_with_replacement(range(30),4));qi={m:i for i,m in enumerate(quartics)}
    full=np.zeros((134,len(pairs)),dtype=np.int64);cofactor=np.zeros_like(full)
    for j,m in enumerate(zpairs):full[:,pi[m]]=high[:,j]
    for i in range(29):full[:,pi[i,29]]=linear[:,i]
    full[:,pi[29,29]]=constant
    start=next(i for i,line in enumerate(lines)if line.startswith('SAVED_COFACTOR_TRIAL '))+1
    for line in lines[start:]:
        kind,*tokens=line.split()
        if kind not in ['H','L','D']:break
        a=list(map(int,tokens))
        if kind=='H':r,i,j,c=a;k=pi[i,j]
        elif kind=='L':r,i,c=a;k=pi[i,29]
        else:r,c=a;k=pi[29,29]
        assert cofactor[r,k]==0;cofactor[r,k]=c
    expected=np.zeros(len(quartics),dtype=np.int64)
    for mon,a,b in zip(basis,u,v):
        for i,c in enumerate(ell):
            if a and c:
                m=tuple(sorted((*mon,i,*([29]*(3-len(mon))))));expected[qi[m]]+=a*c
        m=tuple(sorted((*mon,*([29]*(4-len(mon))))));expected[qi[m]]-=b
    product=np.array([[qi[tuple(sorted((*a,*b)))]for b in pairs]for a in pairs],dtype=np.int32)
    actual=np.zeros(len(quartics),dtype=np.int64)
    for r in range(134):
        a=np.flatnonzero(cofactor[r]);b=np.flatnonzero(full[r])
        np.add.at(actual,product[np.ix_(a,b)].ravel(),(cofactor[r,a][:,None]*full[r,b][None,:]).ravel())
    assert np.array_equal(actual%P,expected%P)
    return len(quartics)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate',type=Path,default=ROOT/'reports/S61_AFFINE_CYCLIC_OPEN_20261103_v62/native-certificate.txt')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--timeout',type=int,default=120)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    lines=args.certificate.read_text().splitlines()
    ell=list(map(int,next(line for line in lines if line.startswith('LINEAR_FUNCTIONAL ')).split()[1:]))
    chi=list(map(int,next(line for line in lines if line.startswith('CHI ')).split()[1:]))
    gcd=list(map(int,next(line for line in lines if line.startswith('GCD ')).split()[1:]))
    coordinates={int(line.split()[1]):list(map(int,line.split()[2:])) for line in lines if line.startswith('COORDINATE ')}
    high,cross,base,pairs,_,_,bpairs,v,inputs,_,_=build();v=v%P
    linear=cross@v%P;constant=base@np.array([int(v[i])*int(v[j])%P for i,j in bpairs],dtype=np.int64)%P
    action_coefficients=check_dense_action(lines,ell,high,linear,constant)
    text='\n'.join([
        'R=(ZZ/32003)[s];',
        'poly62=l->sum(toList(0..(#l-1)),i->(l#i)*s^i);',
        'chi62=poly62 '+literal(chi)+';',
        'p62=apply('+literal([coordinates[i] for i in range(29)])+',poly62);',
        'ell62='+literal(ell)+';',
        'assert(sum(toList(0..28),i->(ell62#i)*(p62#i))==s);',
        'pairs62='+literal(pairs)+';',
        'products62=apply(pairs62,ij->((p62#(ij#0))*(p62#(ij#1)))%chi62);',
        'high62='+literal(high.tolist())+';',
        'linear62='+literal(linear.tolist())+';',
        'constant62='+literal(constant.tolist())+';',
        'g62=chi62;',
        'scan(toList(0..133),r->(f62:=sum(toList(0..434),j->((high62#r)#j)*(products62#j))+sum(toList(0..28),i->((linear62#r)#i)*(p62#i))+constant62#r;g62=gcd(g62,f62)));',
        'g62=g62*(sub(leadCoefficient g62,coefficientRing R))^(-1);',
        'expected62=poly62 '+literal(gcd)+';',
        'assert(g62==expected62 and numColumns basis(R/ideal g62)==5);',
        'assert(((p62#0)-1)%g62==0);',
        'assert(gcd(g62,diff(s,g62))==1);',
        'print "CYCLIC_OPEN_ALL134_REMAINDERS_PASS";stdio << flush;exit 0;',
    ])+'\n'
    script=args.output/'review.m2';log=args.output/'review.log';script.write_text(text)
    with log.open('w')as stream:subprocess.run([shutil.which('M2'),'--script',str(script.resolve())],stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=args.timeout)
    assert 'CYCLIC_OPEN_ALL134_REMAINDERS_PASS' in log.read_text()
    inputs += [Path(__file__),ROOT/'code/s61_affine_seventeen_base_lift_20261103_v62.py',args.certificate]
    report=dict(all_passed=True,scope='independent CloudM2 substitution of all134 actual quadrics into saved coordinate polynomials',
                prime=P,monic_polynomial_degree=1031,actual_sample_gcd_degree=5,
                t_equals_one=True,gcd_squarefree=True,linear_functional_identity=True,
                independent_dense_h_action_original134_coefficient_identity=True,
                independently_verified_action_coefficients=action_coefficients,
                krylov_determinant_independently_recomputed=False,global_gcd_support_computed=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(inputs)},
                evidence_sha256={p.name:sha(p) for p in args.output.iterdir() if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',report['seconds'])


if __name__=='__main__':main()

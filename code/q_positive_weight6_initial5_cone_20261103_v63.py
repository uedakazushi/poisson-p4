#!/usr/bin/env python3
"""Bounded exactQQ completion and two-plane radical analysis of the fifth data."""
import argparse
from pathlib import Path
import shutil
import subprocess
import time
from q_positive_weight6_quadratic_model_20261103_v63 import read,sha


def polynomial(row,words):
    terms=[]
    for coefficient,word in zip(row,words):
        if coefficient!='0':terms.append('('+coefficient+')*'+'*'.join('z'+str(j)for j in word))
    return '+'.join(terms)or'0'


def main(folder):
    started=time.monotonic();source=read(folder/'fifth-complete-spanning.json.gz')
    selected=read(folder/'selected-fifth.json');words=source['normal_monomials5']
    allpolys=[polynomial(row,words)for row in source['complete_fifth_order_new329']+source['complete_quartic_kernel_fifth_carries48']]
    selectedpoly=[p.replace('**','^')for p in selected['lower21_polynomials']+selected['selected_QQ_fifth_polynomials']]
    text='R=QQ[z0,z1,z2,z3,z4,z5,z6,z7];\nI=ideal('+','.join(selectedpoly)+');\nG=gb I;\n'
    text+='ALL=matrix{{'+','.join(allpolys)+'}};assert(ALL%G==0);\n'
    text+='print("ALL377_EXACT_QQ_REDUCTIONS_ZERO");print("HF5="|toString hilbertFunction(5,R/I));\n'
    text+='A=QQ[u0,u2,u3,u4,u6];fa=map(A,R,{u0,0,u2,u3,u4,0,u6,0});IA=fa I;\n'
    text+='print("A_DIM="|toString dim(A/IA));print("A_GB="|toString gens gb IA);\n'
    text+='print("A_MINIMAL_PRIMES="|toString minimalPrimes IA);\n'
    text+='B=QQ[v0,v2,v3,v5,v6];fb=map(B,R,{v0,0,v2,v3,-2*v6,v5,v6,0});IB=fb I;\n'
    text+='print("B_DIM="|toString dim(B/IB));print("B_GB="|toString gens gb IB);\n'
    text+='print("B_MINIMAL_PRIMES="|toString minimalPrimes IB);\nprint("COMPLETE_FIFTH_AND_RADICAL_SPLIT_DONE");\n'
    (folder/'complete-fifth-cone.m2').write_text(text)
    executable=shutil.which('M2');assert executable
    status='not_run'
    try:
        with (folder/'complete-fifth-cone.log').open('w')as log:
            result=subprocess.run([executable,'--script',str(folder/'complete-fifth-cone.m2')],stdout=log,stderr=subprocess.STDOUT,timeout=50)
        status='completed'if result.returncode==0 else 'exit'+str(result.returncode)
    except subprocess.TimeoutExpired:status='timeout50s'
    print('EXACT FIFTH CONE FINALIZER',status,'SECONDS',time.monotonic()-started,flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--folder',type=Path,required=True);main(ap.parse_args().folder)

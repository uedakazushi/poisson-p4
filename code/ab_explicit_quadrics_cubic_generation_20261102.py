#!/usr/bin/env python3
"""Construct all AB quadrics in the actual 45 S61 coordinates; certify no new cubics."""
import argparse
import hashlib
import itertools
import json
import math
import time
from collections import defaultdict
from pathlib import Path
import numpy as np
import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring
from s61_model import S61, X

ROOT=Path(__file__).resolve().parents[1]
NOTE='research/notes/AB_EXPLICIT_QUADRICS_CUBIC_GENERATION_20261102.md'

def rank_mod(rows,prime):
    cols=sorted(set().union(*(v.keys() for v in rows)));indices={m:i for i,m in enumerate(cols)}
    a=np.zeros((len(rows),len(cols)),dtype=np.int64)
    for i,row in enumerate(rows):
        for k,c in row.items(): a[i,indices[k]]=int(s.numer(c))*pow(int(s.denom(c)),-1,prime)%prime
    r=0;piv=[]
    for j in range(len(cols)):
        nz=np.flatnonzero(a[r:,j]);
        if not nz.size:continue
        b=r+int(nz[0]);a[[b,r]]=a[[r,b]]
        a[r,j:]=(a[r,j:]*pow(int(a[r,j]),-1,prime))%prime
        a[r+1:,j:]=(a[r+1:,j:]-a[r+1:,j,None]*a[r,None,j:])%prime
        piv.append(j);r+=1
        if r==len(rows):break
    return len(cols),r,piv

def verify():
    started=time.monotonic();old=ROOT/'reports/AB_ORIGIN_DEGREE_FOUR_20261101.json';hpath=ROOT/'reports/AB_RATIONAL_SINGULARITIES_HILBERT_20261102.json'
    r=json.loads(old.read_text());h=json.loads(hpath.read_text());assert r['all_passed'] and h['all_passed']
    model=S61();q=s.symbols('q:45');it=iter(q);params=[0 if wt==1 else next(it) for wt in model.reduced['parameter_weights']]
    P=model.reconstruct(params).subs(X[0],1);y=X[1:];J=s.zeros(4);J[0,1]=J[2,3]=1;J-=J.T
    mons2=list(itertools.combinations_with_replacement(range(4),2));coords=[]
    for i,j in itertools.combinations(range(4),2):
        p=s.Poly(P[i+1,j+1]-y[i]*P[0,j+1]+y[j]*P[0,i+1]-J[i,j],*y)
        assert all(sum(m)==2 for m,c in p.terms());coords.extend(p.coeff_monomial(y[a]*y[b]) for a,b in mons2)
    E=s.Matrix(coords).jacobian(q);assert E*s.Matrix(q)==s.Matrix(coords)
    C=E[r['coordinate_indices'],:];det=C.det();assert det;inv=C.inv()
    R,*_=ring(r['variables'],QQ);rp=[R.from_dict({tuple(ex):QQ.convert(s.Rational(co)) for ex,co in p}) for p in r['coordinate_polynomials']]
    qp=[sum((QQ.convert(inv[i,j])*rp[j] for j in range(45) if inv[i,j]),R.zero) for i in range(45)]
    weights=[]
    for f in qp:
        w={tuple(sum(ex[k]*r['variable_weights'][k][j] for k in range(10)) for j in range(2)) for ex in f};assert len(w)==1;weights.append(next(iter(w)))
    groups=defaultdict(list)
    def weight(ind):return tuple(sum(weights[i][j] for i in ind) for j in range(2))
    for ij in itertools.combinations_with_replacement(range(45),2):groups[weight(ij)].append(ij)
    relations=[];qb=[]
    for w,pairs in sorted(groups.items()):
        products=[qp[i]*qp[j] for i,j in pairs];mon=sorted(set().union(*(p.keys() for p in products)))
        ns=s.Matrix([[p.get(m,0) for p in products] for m in mon]).nullspace()
        for v in ns:
            row={ij:c for ij,c in zip(pairs,v) if c};assert sum((QQ.convert(c)*qp[i]*qp[j] for (i,j),c in row.items()),R.zero)==R.zero;relations.append(row)
        qb.append(dict(weight=w,columns=len(pairs),kernel=len(ns)))
    assert len(relations)==274 and sum(x['columns'] for x in qb)==1035
    blocks=defaultdict(list)
    for i in range(45):
        for rel in relations:
            row={tuple(sorted((i,)+ij)):c for ij,c in rel.items()};blocks[weight(next(iter(row)))].append(row)
    prime=32003;cb=[]
    for w,rows in sorted(blocks.items()):
        cols,n,piv=rank_mod(rows,prime);cb.append(dict(weight=w,rows=len(rows),active_columns=cols,rank=n,pivot_columns=piv))
    rank=sum(v['rank'] for v in cb);assert rank==math.comb(47,3)-6691==9524
    paths=[old,hpath,ROOT/'code/s61_model.py',Path(__file__),ROOT/NOTE]
    return dict(id='AB_EXPLICIT_QUADRICS_CUBIC_GENERATION_20261102',all_passed=True,
                coordinate_order='The ordered weight-two parameters of S61, q0..q44.',coordinate_minor=str(det),coordinate_weights=weights,
                quadratic_relations=[[[i,j,str(c)] for (i,j),c in row.items()] for row in relations],quadratic_blocks=qb,
                multiplication_prime=prime,cubic_ranks=cb,quadratic_count=274,cubic_ideal_dimension=9524,minimal_cubic_generators=0,
                full_ideal_quadratically_generated_proved=False,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.monotonic()-started)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--output',required=True,type=Path);arg=a.parse_args();r=verify()
    with arg.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))

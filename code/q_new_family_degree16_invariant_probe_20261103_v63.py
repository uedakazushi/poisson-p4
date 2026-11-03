#!/usr/bin/env python3
"""Bounded degree16 invariant restriction probe; sampled kernels are not identities."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import random
import time
import numpy as np
import sympy as s
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import sp_basis

ROOT=Path(__file__).resolve().parents[1]
PRIME=32003


def read(path):
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:return json.load(f)


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()


def mod(value):
    q=Fraction(str(value));return q.numerator*pow(q.denominator,-1,PRIME)%PRIME


def matrix(rows):return s.Matrix([[s.Rational(x) for x in row] for row in rows])
def modmat(M):return np.array([[mod(c) for c in row] for row in M.tolist()],dtype=np.int64)


def necklaces(n):
    return [w for w in itertools.product(range(2),repeat=n)
            if w==min(w[i:]+w[:i] for i in range(n))]


def echelon(rows):
    basis={};selected=[]
    for index,row in enumerate(rows):
        row=[int(x)%PRIME for x in row]
        for pivot,vector in sorted(basis.items()):
            coefficient=row[pivot]
            if coefficient:row=[(a-coefficient*b)%PRIME for a,b in zip(row,vector)]
        pivot=next((i for i,c in enumerate(row) if c),None)
        if pivot is None:continue
        inverse=pow(row[pivot],-1,PRIME)
        basis[pivot]=[(c*inverse)%PRIME for c in row]
        selected.append(index)
    return selected,sorted(basis),basis


def character(weights,maximum=18):
    dp=[defaultdict(int) for _ in range(maximum+1)];dp[0][0,0]=1
    for a,b in weights:
        for degree in range(1,maximum+1):
            for (i,j),count in list(dp[degree-1].items()):dp[degree][i+a,j+b]+=count
    rho=(2,1);weyl=[]
    for perm in [(0,1),(1,0)]:
        for signs in itertools.product([-1,1],repeat=2):
            wr=tuple(signs[i]*rho[perm[i]] for i in range(2))
            weyl.append((tuple(rho[i]-wr[i] for i in range(2)),signs[0]*signs[1]*(-1 if perm==(1,0) else 1)))
    records=[]
    for d,values in enumerate(dp):
        assert sum(values.values())==int(s.binomial(15+d,d))
        records.append(dict(degree=d,invariant_dimension=sum(sign*values[delta] for delta,sign in weyl),
                            Weyl_terms=[dict(weight=list(delta),sign=sign,multiplicity=values[delta]) for delta,sign in weyl]))
    assert [records[d]['invariant_dimension'] for d in [15,16,17,18]]==[0,20,0,9]
    return records


def main(output,count):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    d8path=ROOT/'reports/Q_NEW_FAMILY_DEGREE8_INVARIANT_EXCLUSION_20261103_v63/certificate.json'
    d12path=ROOT/'reports/Q_NEW_FAMILY_DEGREE12_14_INVARIANT_EXCLUSION_20261103_v63/certificate.json'
    mappath=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/map.json.gz'
    data,d8,d12,payload=map(read,[invpath,d8path,d12path,mappath])
    B=matrix(data['primitive24_by16_basis']);basisrows=list(B.T.rref()[1]);left=B[basisrows,:].inv()
    O=matrix(data['symplectic_form16']);rep=[matrix(R) for R in data['all_ten_representation_matrices']]
    gs=sp_basis();gram=s.Matrix([[s.trace(A*C) for C in gs] for A in gs])
    om,oi,gi=modmat(O),modmat(O.inv()),modmat(gram.inv());reps=[modmat(R) for R in rep]
    w4,w6,w8=necklaces(4),necklaces(6),necklaces(8)
    assert len(w8)==36
    selected8=d8['independent_five_columns'];selected12=d12['degree12_basis_columns']
    productpairs=list(itertools.combinations_with_replacement(range(5),2))
    definitions=[dict(kind='trace8',word=list(w)) for w in w8]
    definitions += [dict(kind='F4_times_trace6',word=list(w6[i])) for i in selected12]
    definitions += [dict(kind='basis8_product',indices=list(pair)) for pair in productpairs]
    assert len(definitions)==60
    def invariants(b):
        b=np.array(b,dtype=np.int64).reshape(16,1)
        V=np.hstack([om@R@b%PRIME for R in reps])
        mu=(b.T@V%PRIME).reshape(10,1)*pow(2,-1,PRIME)%PRIME
        co=gi@mu%PRIME
        theta=sum((int(c)*R for c,R in zip(co[:,0],reps)),np.zeros((16,16),dtype=np.int64))%PRIME
        H=(theta+oi@V%PRIME@gi%PRIME@V.T)%PRIME
        tables={():np.eye(16,dtype=np.int64)}
        for length in range(1,5):
            for w in itertools.product(range(2),repeat=length):tables[w]=tables[w[:-1]]@[theta,H][w[-1]]%PRIME
        def trace(w):
            middle=len(w)//2
            return int(np.sum(tables[w[:middle]]*tables[w[middle:]].T,dtype=np.int64)%PRIME)
        f4=int((mu.T@co)[0,0]%PRIME)
        values8=[trace(w) for w in w4]+[f4*f4%PRIME]
        values8=[values8[j] for j in selected8]
        return [trace(w) for w in w8]+[f4*trace(w6[j])%PRIME for j in selected12]+[values8[i]*values8[j]%PRIME for i,j in productpairs]
    model=S61();z,u,beta=primitive_beta(model)
    ones=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    beta24=[beta[i,j].coeff(t) for i,j in itertools.combinations(range(4),2) for t in z]
    betaU=s.Matrix([[s.expand(c).coeff(v) for v in u] for c in beta24])
    beta16=modmat(left*betaU[basisrows,:])
    def actual(point):
        values=[]
        for node in payload['nodes']:
            op=node[0]
            if op=='q':v=mod(node[1])
            elif op=='v':v=point[payload['input_names'][node[1]]]%PRIME
            elif op=='+':v=sum(values[i] for i in node[1])%PRIME
            elif op=='*':
                v=1
                for i in node[1]:v=v*values[i]%PRIME
            elif op=='inv':v=pow(values[node[1]],-1,PRIME)
            else:raise AssertionError(op)
            values.append(v)
        linear=np.array([values[payload['output_u1_through_u61'][j]] for j in ones],dtype=np.int64)
        return list(map(int,beta16@linear%PRIME))
    rng=random.Random(20261103)
    ambient=[[rng.randrange(PRIME) for _ in range(16)] for _ in range(count)]
    actual_points=[];actual_beta=[]
    while len(actual_beta)<count:
        point={f'n{i}':rng.randint(1,9) for i in range(3)}
        point.update(scale=1);point.update({f'base{i}':rng.randint(-2,2) for i in range(1,5)})
        point.update({f'g{i}':v for i,v in enumerate([1,0,0,1,0,0,0,0,0,0])})
        try:b=actual(point)
        except ValueError:continue
        actual_points.append(point);actual_beta.append(b)
    ambient_values=[invariants(b) for b in ambient]
    actual_values=[invariants(b) for b in actual_beta]
    ar,ac,_=echelon(ambient_values);rr,rc,_=echelon(actual_values)
    chars=character(list(map(tuple,data['Cartan_weights'])))
    result=dict(id='Q_NEW_FAMILY_DEGREE16_INVARIANT_PROBE_20261103_v63',all_passed=True,
        prime=PRIME,number_of_samples=count,number_of_explicit_invariants=len(definitions),
        exact_character_dimensions_15_16_17_18=[0,20,0,9],character_records=chars,
        invariant_definitions=definitions,ambient_modular_rank=len(ar),actual_new_map_modular_rank=len(rr),
        ambient_independent_rows=ar,ambient_independent_columns=ac,
        actual_independent_rows=rr,actual_independent_columns=rc,
        ambient_beta16=ambient,actual_points=actual_points,actual_beta16=actual_beta,
        ambient_evaluation=ambient_values,actual_evaluation=actual_values,
        modular_kernels_not_claimed_as_QQ_identities=True,
        degree16_equation_or_exclusion_proved=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),invpath,d8path,d12path,mappath]},
        seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','number_of_explicit_invariants','ambient_modular_rank','actual_new_map_modular_rank','seconds']},indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--samples',type=int,default=40);args=parser.parse_args();main(args.output.resolve(),args.samples)

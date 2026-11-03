#!/usr/bin/env python3
"""Independent generic-H Jacobian samples using compressed exact beta polynomials."""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import random
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.ddm import DDM
from s61_marked_sigma_point_lift_20261103_v62 import load_model,coerce

ROOT=Path(__file__).resolve().parents[1]


def dm(rows,K):return DomainMatrix.from_rep(DDM(rows,(len(rows),len(rows[0])),K))


def verify(folder,output):
    start=time.monotonic();L=load_model(folder);K=s.QQ.algebraic_field(s.sqrt(-2),alias='h');h=K.unit
    zero,one=K.zero,K.one;v=[zero]*13;v[9]=h;v[10]=one
    point=L.lift(h,v,domain=K);q=point['q'];bi,qi=L.one,L.two
    jac=[defaultdict(dict) for _ in range(135)]
    def add(row,column,mon,value):
        if not value:return
        mon=tuple(sorted(mon));cell=jac[row][column]
        cell[mon]=cell.get(mon,zero)+value
    scalar_cache={}
    for row,form in enumerate(L.original):
        for mon,n,d in form:
            value=s.Rational(n,d)
            if value not in scalar_cache:scalar_cache[value]=coerce(K,value)
            c=scalar_cache[value];qs=[qi.index(i) for i in mon if i in qi];bs=[bi.index(i) for i in mon if i in bi]
            if len(qs)==2:
                for k,index in enumerate(qs):add(row,qi[index],bs,c*q[qs[1-k]])
                qp=c*q[qs[0]]*q[qs[1]]
                for k,index in enumerate(bs):add(row,bi[index],bs[:k]+bs[k+1:],qp)
            elif len(qs)==1:add(row,61,bs,c*q[qs[0]])
    for row in jac:
        for column,cell in row.items():
            for mon in list(cell):
                if not cell[mon]:del cell[mon]
    highest=dm([[jac[i].get(index,{}).get((),zero) for index in bi] for i in range(135)],K)
    H=highest.nullspace().transpose();hrows=list(H.transpose().rref()[1])
    H=H*H.extract(hrows,range(H.shape[1])).inv();assert H.shape==(16,9)
    reference=ROOT/'reports/S61_MARKED_SIGMA_OLD_PLANE_AMBIENT_20261103_v62/report.json'
    old=next(c for c in json.loads(reference.read_text())['cases'] if c['case']=='distinct-u' and c['beta']=='actual')
    frozen={qi[25],*[qi[i] for i in old['q_cut_indices']],*[bi[i] for i in hrows]}
    variables=[i for i in range(62) if i not in frozen];assert len(variables)==44
    def rank(beta):
        products={():one}
        def product(mon):
            if mon not in products:
                value=one
                for i in mon:value*=beta[i]
                products[mon]=value
            return products[mon]
        matrix=[[sum((c*product(mon) for mon,c in jac[i].get(column,{}).items()),zero)
                 for column in variables] for i in range(135)]
        return dm(matrix,K).rank()
    actual_rank=rank(point['beta']);assert actual_rank==43
    samples=[list(range(2,11)),[2,5,11,19,31,43,59,71,89]]
    for seed in [62,620,6200]:samples.append(random.Random(seed).sample(list(range(-97,98)),9))
    results=[]
    for numbers in samples:
        beta=(H*dm([[coerce(K,value)] for value in numbers],K)).to_list();beta=[row[0] for row in beta]
        value=rank(beta)
        results.append(dict(H_coefficients=numbers,full_ambient_Jacobian_rank=value))
        print(numbers,'RANK',value,'SECONDS',time.monotonic()-start,flush=True)
    report=dict(id='S61_MARKED_SIGMA_OLD_PLANE_RANDOM_JACOBIAN_20261103_v62',all_passed=True,
                field='QQ(h), h^2=-2',actual_rank=actual_rank,samples=results,
                compressed_beta_polynomial_terms=sum(len(cell) for row in jac for cell in row.values()),
                generic_reduced_smooth_incidence_proved=any(r['full_ambient_Jacobian_rank']==44 for r in results),
                interpretation='One rank44 sample proves the SAME smooth17 incidence has generic multiplicity1. Rank43 samples alone give no generic upper rank43.',
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in [Path(__file__),ROOT/'code/s61_marked_sigma_point_lift_20261103_v62.py',reference]},
                seconds=time.monotonic()-start)
    with output.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--model-folder',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    r=verify(args.model_folder,args.output);print('PASS',r['id'],r['seconds'])

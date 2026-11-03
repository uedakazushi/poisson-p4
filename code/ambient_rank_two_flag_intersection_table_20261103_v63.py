#!/usr/bin/env python3
"""Exact finite double-flag table for rank2 CA/CB boundary intersections.

The table records actual canonical linear intersections and known-image tests.
An image differential at one point is only a dimension LOWER bound unless an
independent source upper bound matches it. No rank2 exhaustion by N is assumed.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import raw,canonical,vector,X,h,x,y,z,w

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))
BASIS=list(itertools.product(PAIRS,MONS));INDEX={entry:i for i,entry in enumerate(BASIS)}

def permutation(values,perm):
    out=s.zeros(150,1)
    for index,((i,j),mon)in enumerate(BASIS):
        if not values[index]:continue
        p,q=perm[i],perm[j];sign=1 if p<q else -1
        out[INDEX[(tuple(sorted((p,q))),tuple(sorted(perm[k]for k in mon)))]]+=sign*values[index]
    return out

def lie(values,a,b):
    out=s.zeros(150,1)
    for index,((i,j),mon)in enumerate(BASIS):
        c=values[index]
        if not c:continue
        for left,right in ([(a,j)]if i==b else [])+([(i,a)]if j==b else []):
            if left!=right:
                sign=1 if left<right else -1
                out[INDEX[(tuple(sorted((left,right))),mon)]]+=sign*c
        if a in mon:
            m=list(mon);m.remove(a);m.append(b)
            out[INDEX[((i,j),tuple(sorted(m)))]]-=mon.count(a)*c
    return out

def representatives():
    todo=set(itertools.permutations(range(5)));out=[]
    swap=list(range(5));swap[3],swap[4]=4,3
    while todo:
        p=min(todo);orbit={p,tuple(swap[v]for v in p),tuple(p[swap[i]]for i in range(5)),tuple(swap[p[swap[i]]]for i in range(5))}
        todo-=orbit;out.append(p)
    assert len(out)==33
    return out

def main(output):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    ca_path=ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63/report.json'
    cb_path=ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json'
    ca=json.loads(ca_path.read_text());cb=json.loads(cb_path.read_text())
    Aall=s.Matrix.hstack(*(s.Matrix(list(map(s.Rational,v)))for v in ca['basis_canonical150']))
    A=Aall[:,[0,1,2,*range(4,20)]];assert A.rank()==19
    cubics=list(map(s.sympify,ca['cubic_monomials']))
    def rb(f,g,H):
        P=raw(f,0,H);P[3,4]=g;P[4,3]=-g
        P[2,3]-=s.diff(f,x)*z/2;P[3,2]=-P[2,3]
        P[2,4]-=s.diff(f,x)*w/2;P[4,2]=-P[2,4]
        return canonical(P)
    BB={
        'b_zero':s.Matrix.hstack(*(vector(rb(f,0,0))for f in [h*h,h*x,x*x]),*(vector(rb(0,0,H))for H in cubics)),
        'a_zero':s.Matrix.hstack(*(vector(rb(0,f,0))for f in [h*h,h*x,x*x]),*(vector(rb(0,0,H))for H in cubics))}
    assert all(B.rank()==19 for B in BB.values())
    positions={(i,j)for i,j in itertools.product(range(5),repeat=2)
        if not(i==0 and j!=0)and not(i==1 and j not in (0,1))and not(j==2 and i!=2)}
    records=[]
    for perm in representatives():
        pairorbit=25-len(positions&{(perm[i],perm[j])for i,j in positions})
        for name,B0 in BB.items():
            B=s.Matrix.hstack(*(permutation(B0[:,j],perm)for j in range(19)))
            M=A.row_join(-B);kernel=M.nullspace();k=len(kernel)
            KA=s.Matrix.hstack(*(v[:19,0]for v in kernel))if k else s.zeros(19,0)
            L=A*KA;assert A*KA==B*s.Matrix.hstack(*(v[19:,0]for v in kernel))if k else True
            inside=all(KA[2,j]==0 for j in range(k))
            rec=dict(relative_permutation=list(perm),CB_endpoint=name,
                full_flag_pair_orbit_dimension=pairorbit,linear_intersection_dimension=k,
                intersection_basis_in_CA_coordinates=[[str(v)for v in KA[:,j]]for j in range(k)],
                incidence_projective_dimension_upper=pairorbit+k-1 if k else None,
                entire_linear_family_lies_in_N_with_the_CA_flag=inside)
            if k:
                point=L*s.Matrix([j%7+1 for j in range(k)])
                cols=L.row_join(s.Matrix.hstack(*(lie(point,a,b)for a,b in itertools.product(range(5),repeat=2))))
                rank=len(cols.rref()[1]);rec['actual_image_dimension_lower_bound_at_QQ_point']=rank-1
                rec['point_CA_f_x_squared_coefficient']=str((KA*s.Matrix([j%7+1 for j in range(k)]))[2])
            records.append(rec)
            print('FLAG',perm,name,'intersection',k,'sourceupper',rec['incidence_projective_dimension_upper'],'Nsameflag',inside,'imageLower',rec.get('actual_image_dimension_lower_bound_at_QQ_point'),flush=True)
    result=dict(id='AMBIENT_RANK_TWO_FLAG_INTERSECTION_TABLE_20261103_v63',all_passed=True,
        exact_coefficient_field='QQ',parabolic_Weyl_subgroup='swap(z,w)',double_cosets=33,
        endpoint_cases=66,nonempty_cases=sum(bool(v['linear_intersection_dimension'])for v in records),
        entire_N_same_flag_cases=sum(v['entire_linear_family_lies_in_N_with_the_CA_flag']for v in records),
        records=records,global_intersection_exhaustion_by_N_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [Path(__file__),ca_path,cb_path]},
        seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('66-CASE TABLE SAVED',result['seconds'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)

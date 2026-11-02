#!/usr/bin/env python3
"""Raw weights and the parabolic obstruction to quintic/sextic AB generators."""
import argparse
from collections import Counter
import hashlib
from itertools import combinations, combinations_with_replacement, product
import json
from pathlib import Path
import time

import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
STEM='AB_IDEAL_HIGH_GENERATOR_VANISHING_20261102'
PAIRS=list(combinations(range(4),2))
MONOMIALS=list(combinations_with_replacement(range(4),2))
BASIS=list(product(PAIRS,MONOMIALS))
J=sp.diag(sp.Matrix([[0,1],[-1,0]]),sp.Matrix([[0,1],[-1,0]]))
VECTOR_WEIGHTS=[(1,0),(-1,0),(0,1),(0,-1)]


def weight(i,j,r,s):
    return tuple(VECTOR_WEIGHTS[i][k]+VECTOR_WEIGHTS[j][k]
                 -VECTOR_WEIGHTS[r][k]-VECTOR_WEIGHTS[s][k] for k in range(2))


def place(matrix,column,i,j,r,s,value):
    if i==j or not value:
        return
    sign=1 if i<j else -1
    index=BASIS.index((tuple(sorted((i,j))),tuple(sorted((r,s)))))
    matrix[index,column]+=sign*value


def bivector_action(A):
    result=sp.zeros(60,60)
    for column,((i,j),(r,s)) in enumerate(BASIS):
        for h in range(4):
            place(result,column,h,j,r,s,A[h,i])
            place(result,column,i,h,r,s,A[h,j])
            place(result,column,i,j,h,s,-A[r,h])
            place(result,column,i,j,r,h,-A[s,h])
    return result


def raw_modules():
    triples=list(combinations(range(4),3))
    differential=sp.zeros(16,60)
    raw_weights=[]
    for column,((i,j),(r,s)) in enumerate(BASIS):
        raw_weights.append(weight(i,j,r,s))
        for h,(k,l,m) in enumerate(triples):
            for first,left,right in [(k,l,m),(l,m,k),(m,k,l)]:
                if set((left,right))!={i,j}:
                    continue
                sign=1 if (left,right)==(i,j) else -1
                differential[4*h+s,column]+=sign*J[first,r]
                differential[4*h+r,column]+=sign*J[first,s]
    assert differential.rank()==15
    jacobi_spaces={}
    ambient=Counter()
    for w in sorted(set(raw_weights)):
        indices=[i for i,v in enumerate(raw_weights) if v==w]
        local=differential[:,indices].nullspace()
        for vector in local:
            v=sp.zeros(60,1)
            for j,c in zip(indices,vector):
                v[j]=c
            jacobi_spaces.setdefault(sum(w),[]).append(v)
        ambient[w]=len(local)
    assert sum(ambient.values())==45
    fourth=json.loads((ROOT/'reports/AB_ORIGIN_DEGREE_FOUR_20261101.json').read_text())
    selected=fourth['coordinate_indices']
    W=sp.Matrix.hstack(*(v for space in jacobi_spaces.values() for v in space))
    assert W[selected,:].rank()==45
    assert [list(raw_weights[i]) for i in selected]==fourth['coordinate_weights']

    # Fixed vertical L. The two components of F are arbitrary binary cubics.
    embedding=sp.zeros(60,11)
    column=0
    for component in (0,2):
        for k in range(4):
            power=[k,3-k]
            for j,pindex in enumerate((1,3)):
                if not power[j]:
                    continue
                remaining=power.copy()
                remaining[j]-=1
                mon=[0]*remaining[0]+[2]*remaining[1]
                place(embedding,column,component,pindex,*mon,power[j])
            column+=1
    for r,s in combinations_with_replacement((0,2),2):
        place(embedding,column,1,3,r,s,1)
        column+=1
    assert embedding.rank()==11 and differential*embedding==sp.zeros(16,11)
    fiber_weights=[]
    for j in range(11):
        weights={raw_weights[i] for i in range(60) if embedding[i,j]}
        assert len(weights)==1
        fiber_weights.append(weights.pop())
    negate=lambda w:(-w[0],-w[1])
    kernel=Counter({negate(w):n for w,n in ambient.items()})
    kernel.subtract(Counter(map(negate,fiber_weights)))
    assert min(kernel.values())>=0 and sum(kernel.values())==34
    core_weights=[negate(w) for w in fiber_weights[:8]]

    shear=sp.zeros(4)
    shear[1,0]=shear[3,2]=1
    assert shear.T*J+J*shear==sp.zeros(4)
    h=shear.T*shear-shear*shear.T
    assert h==sp.diag(1,-1,1,-1) and h*shear-shear*h==-2*shear
    action=bivector_action(shear)
    injections={}
    for central in (4,2):
        W=sp.Matrix.hstack(*jacobi_spaces[central])
        target=sp.Matrix.hstack(*jacobi_spaces[central-2])
        image=action*W
        assert image.rank()==W.cols and target.row_join(image).rank()==target.cols
        injections[str(central)]={'domain_dimension':W.cols,'target_dimension':target.cols,
                                  'rank_QQ':image.rank()}
    radical_checks=[]
    for values in [(1,0,0),(0,1,0),(0,0,1)]:
        A=sp.zeros(4)
        A[1,0],A[1,2],A[3,0],A[3,2]=values[0],values[1],values[1],values[2]
        image=bivector_action(A)*embedding
        assert image[:,8:]==sp.zeros(60,3)
        assert embedding[:,8:].row_join(image).rank()==3
        radical_checks.append(image.rank())
    return ambient,core_weights,kernel,fiber_weights,injections,radical_checks


def multiply(left,right):
    result=Counter()
    for (a,b),n in left.items():
        for (x,y),m in right.items():
            result[a+x,b+y]+=n*m
    return result


def power_character(weights,n,exterior=False):
    result=Counter()
    choose=combinations if exterior else combinations_with_replacement
    for indices in choose(range(len(weights)),n):
        result[tuple(sum(weights[i][k] for i in indices) for k in range(2))]+=1
    return result


def irreducibles(char):
    char=+char
    result=[]
    while char:
        a,b=max(char)
        mult=char[a,b]
        assert a>=b and mult>0
        result.append([a,b,mult])
        for k in range(a-b+1):
            w=(a-k,b+k)
            char[w]-=mult
            assert char[w]>=0
            if not char[w]:
                del char[w]
    return result


def tor_characters(core_weights):
    binary_square=[(2,0),(1,1),(0,2)]
    dual=[(-a,-b) for a,b in binary_square]
    coordinate_ring=[]
    for n in range(7):
        full=multiply(power_character(binary_square,n),power_character(dual,n))
        if n:
            full.subtract(multiply(power_character(binary_square,n-1),power_character(dual,n-1)))
        shifted=Counter({(a+n,b+n):m for (a,b),m in full.items() if m})
        assert min(shifted.values())>0 and sum(shifted.values())==(n+1)**3
        coordinate_ring.append(shifted)
    assert coordinate_ring[1]==Counter(core_weights)
    shifts={0:0,2:1,3:2,4:3,6:4}
    tors={}
    for n in range(7):
        value=Counter()
        for j in range(n+1):
            term=multiply(power_character(core_weights,j,True),coordinate_ring[n-j])
            if j%2: value.subtract(term)
            else: value.update(term)
        value=Counter({w:m for w,m in value.items() if m})
        if n not in shifts:
            assert not value
            continue
        i=shifts[n];sign=(-1)**i
        char=Counter({w:sign*m for w,m in value.items()})
        assert all(m>0 and sum(w)==2*n for w,m in char.items())
        tors[i]=char
    assert [sum(tors[i].values()) for i in range(5)]==[1,9,16,9,1]
    return tors


def encode(char):
    return [[a,b,n] for (a,b),n in sorted(char.items()) if n]


def verify():
    start=time.time()
    resolution_path=ROOT/'reports/AB_ORIGIN_NORMALITY_20261102.json'
    canonical=json.loads(resolution_path.read_text())
    assert canonical['all_passed'] and canonical['fiber_resolution_shifts']==[0,2,3,4,6]
    hilbert_path=ROOT/'reports/AB_RATIONAL_SINGULARITIES_HILBERT_20261102.json'
    hilbert=json.loads(hilbert_path.read_text())
    assert hilbert['all_passed']
    ambient,core,kernel,fiber,injections,radical=raw_modules()
    tors=tor_characters(core)
    rows=[]
    for i,r,degree in [(1,3,5),(2,2,5),(3,1,5),(4,0,6)]:
        ch=multiply(power_character(list(kernel.elements()),r,True),tors[i])
        highest=irreducibles(ch)
        bott=[]
        for a,b,m in highest:
            if a<=-3:
                bott.append({'Levi':[a,b],'Sp4':[-b-3,-a-3],'multiplicity':m})
        rows.append({'Y_Tor':i,'exterior_linear_degree':r,'internal_degree':degree,
                     'dimension':sum(ch.values()),'Levi_types':highest,'Levi_graded_H3':bott,
                     'maximum_center_after_Serre_dual':-min(sum(w) for w,n in ch.items() if n)-6})
    assert all(not row['Levi_graded_H3'] for row in rows[1:])
    assert [q['maximum_center_after_Serre_dual'] for q in rows]==[2,-4,-10,-18]
    assert injections['4']['domain_dimension']==3
    assert {tuple(item['Sp4']):item['multiplicity'] for item in rows[0]['Levi_graded_H3']}=={
        (0,0):9,(2,0):1,(1,1):1}
    assert irreducibles(tors[4])==[[6,6,1]]
    paths=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),resolution_path,hilbert_path,
           ROOT/'reports/AB_ORIGIN_DEGREE_FOUR_20261101.json',
           ROOT/'research/notes/AB_INCIDENCE_20261101.md',ROOT/'research/notes/AB_ORIGIN_NORMALITY_20261102.md',
           ROOT/'research/notes/AB_ORIGIN_DEGREE_FOUR_COHOMOLOGY_REVIEW_20261101.md']
    return {'id':STEM,'all_passed':True,'original_60_to_45_rank':15,'historical_45_coordinate_restriction_rank':45,
            'ambient_primal_weights':encode(ambient),'fiber_11_primal_weights':fiber,
            'linear_kernel_34_weights':encode(kernel),'linear_kernel_34_Levi_types':irreducibles(kernel),
            'Y_Tor_weights':{str(i):encode(char) for i,char in tors.items()},
            'Y_Tor_Levi_types':{str(i):irreducibles(char) for i,char in tors.items()},
            'Tor4_degree5_and6_Levi_factors':rows,'lowering_injections':injections,
            'three_radical_maps_fiber_into_vertical_ranks':radical,
            'H3_Tor4_degree5_dimension':0,'H3_Tor4_degree6_dimension':0,
            'ideal_minimal_generators_degree_upper':4,
            'argument':'See proof note for parabolic filtration, Serre duality, and highest-vector lowering. No P-equivariant splitting is used.',
            'input_sha256':{str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest() for q in paths},
            'seconds':time.time()-start}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    if args.output.exists(): parser.error(f'output exists: {args.output}')
    result=verify()
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print('PASS',result['Y_Tor_Levi_types'],result['Tor4_degree5_and6_Levi_factors'][0]['Levi_graded_H3'],result['lowering_injections'],result['seconds'])

#!/usr/bin/env python3
"""Complete independent smoothing weights by finite exact Maurer--Cartan steps.

All 150 quadratic bivector monomials are considered. Independent edge weights
give at most one smoothing multidegree to each weight, so there are finitely
many possible corrections. Each correction solves the toric Jacobi differential
over the full rational period field; final projective Jacobi is checked again.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s

from mps_s61_catalogue_20261103_v62 import raw_family,serialize,DIAGRAMS,ROOT
from s61_model import X,jacobi_projective

PAIRS = list(itertools.combinations(range(5),2))
QUADS = list(itertools.combinations(range(5),4))
MONOMIALS = [s.prod(X[i] for i in ids) for ids in itertools.combinations_with_replacement(range(5),4)]


def jacobi_vector(P):
    result=[]
    for quad in QUADS:
        expression=0
        for position,e in enumerate(quad):
            i,j,k=quad[:position]+quad[position+1:]
            jac=sum(P[i,h]*s.diff(P[j,k],X[h])+P[j,h]*s.diff(P[k,i],X[h])
                    +P[k,h]*s.diff(P[i,j],X[h]) for h in range(5))
            expression+=(-1)**position*X[e]*jac
        poly=s.Poly(s.expand(expression),*X)
        result.extend(s.cancel(poly.coeff_monomial(mon)) for mon in MONOMIALS)
    return s.Matrix(result)


def allowed_corrections(record):
    weights=[]
    for edge in record['edges']:
        vector=s.zeros(5,1)
        for i in edge['edge']:vector[i]-=1
        for i,power in edge['orders'].items():vector[int(i)]+=int(power)
        weights.append(vector)
    W=s.Matrix.hstack(*weights)
    if W.rank()!=len(weights):
        raise ValueError('dependent smoothing weights require a separate cycle method')
    rows=list(W.T.rref()[1]);inverse=W[rows,:].inv()
    groups={}
    for i,j in PAIRS:
        for k,l in itertools.combinations_with_replacement(range(5),2):
            weight=s.zeros(5,1)
            weight[k]+=1;weight[l]+=1;weight[i]-=1;weight[j]-=1
            n=inverse*weight[rows,:]
            if W*n!=weight or any(not v.is_Integer or v<0 for v in n) or sum(n)<2:
                continue
            V=s.zeros(5);V[i,j]=X[k]*X[l];V[j,i]=-V[i,j]
            groups.setdefault(tuple(map(int,n)),[]).append(V)
    return groups


def complete(record):
    periods,smoothings,P,pf,A=raw_family(record)
    toric=P.subs(dict.fromkeys(smoothings,0))
    assert not any(jacobi_vector(toric))
    groups=allowed_corrections(record)
    steps=[]
    for n,basis in sorted(groups.items(),key=lambda item:(sum(item[0]),item[0])):
        raw=jacobi_vector(P)
        rhs=s.Matrix([s.Poly(v,*smoothings).coeff_monomial(n) for v in raw])
        if not any(rhs):continue
        columns=[jacobi_vector(toric+V)-jacobi_vector(V) for V in basis]
        matrix=s.Matrix.hstack(*columns)
        active=[i for i in range(350) if rhs[i] or any(matrix[i,j] for j in range(matrix.cols))]
        solution,parameters=matrix[active,:].gauss_jordan_solve(-rhs[active,:])
        solution=solution.subs(dict.fromkeys(parameters,0)).applyfunc(s.cancel)
        assert (matrix*solution+rhs).applyfunc(s.cancel)==s.zeros(350,1)
        term=sum((c*V for c,V in zip(solution,basis)),s.zeros(5))
        coefficient=s.prod(t**power for t,power in zip(smoothings,n))
        P+=coefficient*term
        steps.append(dict(smoothing_multidegree=list(n),weight_space_dimension=len(basis),
                          correction=serialize(term)))
    final=jacobi_vector(P)
    assert not any(final)
    assert all(v==0 for v in jacobi_projective(P))
    return periods,smoothings,P,pf,steps,len(groups)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--families',nargs='+',default=[f'T{i:02d}' for i in range(2,39) if i!=21])
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    diagrams={r['id']:r for r in json.loads(DIAGRAMS.read_text())['diagrams']}
    records=[]
    for name in args.families:
        started=time.monotonic()
        periods,smoothings,P,pf,steps,candidates=complete(diagrams[name])
        record=dict(component_id=name,all_passed=True,period_parameters=list(map(str,periods)),
                    smoothing_parameters=list(map(str,smoothings)),period_pfaffian_open=str(pf),
                    possible_higher_multidegrees=candidates,correction_steps=steps,
                    homogeneous_bivector=serialize(P),
                    exact_projective_Jacobi_checked=True,
                    expected_s61_dimension=diagrams[name]['projective_component_dimension']-9,
                    seconds=time.monotonic()-started)
        (args.output/f'{name}.json').write_text(json.dumps(record,indent=2)+'\n')
        records.append(dict(component_id=name,correction_count=len(steps),seconds=record['seconds']))
        print(name,'PASS','CORRECTIONS',len(steps),'SECONDS',record['seconds'],flush=True)
    inputs=[Path(__file__),Path(__file__).with_name('mps_s61_catalogue_20261103_v62.py'),DIAGRAMS,
            ROOT/'code/s61_model.py']
    report=dict(id='MPS_SMOOTHING_COMPLETION_20261103_v62',all_passed=True,records=records,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                scope='General rational coefficient families with exact finite smoothing corrections. Dominance in S61 and component identification require the separate catalogue check; no global exhaustion.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()

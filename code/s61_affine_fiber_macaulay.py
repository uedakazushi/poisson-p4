#!/usr/bin/env python3
"""Bounded sparse Macaulay probe of the full fixed-first-jet S61 fiber.

The input is the exact T01 point of the 134-equation full S61 model. All 45
weight-two variables remain. QQ constant-row elimination and bounded modular
degree-three row elimination are kept distinct; modular data is not a complex
fiber degree or exhaustion certificate.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import combinations_with_replacement
import json
from pathlib import Path
import time

import sympy as s

from weight1_rank_one_fiber import row_reduce

ROOT = Path(__file__).resolve().parents[1]


def add_term(poly,monomial,coefficient):
    value=poly.get(monomial,Fraction(0))+coefficient
    if value: poly[monomial]=value
    else: poly.pop(monomial,None)


def load_fiber():
    source=ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json'
    certificate=ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json'
    data=json.loads(source.read_text())
    point_record=next(p for p in json.loads(certificate.read_text())['points'] if p['component_id']=='T01')
    point=[Fraction(v) for v in point_record['S61_coordinates_u1_through_u61']]
    columns=[i for i,w in enumerate(data['parameter_weights']) if w==2]
    index={old:new for new,old in enumerate(columns)}
    base=[point[i] for i in columns]
    equations=[]
    for equation in data['equations']:
        poly={}
        for monomial,numerator,denominator in equation:
            coefficient=Fraction(numerator,denominator);free=[]
            for i in monomial:
                if i in index:free.append(index[i])
                else:coefficient*=point[i]
            add_term(poly,tuple(free),coefficient)
        centered={}
        for monomial,coefficient in poly.items():
            assert len(monomial)<=2
            add_term(centered,monomial,coefficient)
            if len(monomial)==1:
                add_term(centered,(),coefficient*base[monomial[0]])
            elif len(monomial)==2:
                i,j=monomial
                add_term(centered,(j,),coefficient*base[i])
                add_term(centered,(i,),coefficient*base[j])
                add_term(centered,(),coefficient*base[i]*base[j])
        assert () not in centered
        equations.append(centered)
    return equations,columns,base,source,certificate


def modular_probe(equations,prime,seconds):
    monomials=[mon for degree in (3,2,1,0)
               for mon in combinations_with_replacement(range(45),degree)]
    index={mon:i for i,mon in enumerate(monomials)}
    modular=[]
    for equation in equations:
        row={}
        for mon,c in equation.items():
            assert c.denominator%prime
            value=c.numerator*pow(c.denominator,-1,prime)%prime
            if value:row[mon]=value
        modular.append(row)
    basis={};started=time.monotonic();rows_processed=0;operations=0
    status='complete'
    generators=[(i,None) for i in range(134)]+[(i,j) for j in range(45) for i in range(134)]
    for equation_index,multiplier in generators:
        if time.monotonic()-started>seconds:
            status='time_limit';break
        row={index[tuple(sorted(mon+(multiplier,))) if multiplier is not None else mon]:c
             for mon,c in modular[equation_index].items()}
        while row:
            pivot=min(row)
            if pivot not in basis:
                inverse=pow(row[pivot],-1,prime)
                basis[pivot]={mon:c*inverse%prime for mon,c in row.items()}
                break
            scale=row[pivot]
            for mon,c in basis[pivot].items():
                value=(row.get(mon,0)-scale*c)%prime
                if value:row[mon]=value
                else:row.pop(mon,None)
            operations+=1
            if operations%100==0 and time.monotonic()-started>seconds:
                status='time_limit';break
        if status!='complete':break
        rows_processed+=1
    degrees=Counter(len(monomials[pivot]) for pivot in basis)
    linear=[[[list(monomials[mon]),c] for mon,c in sorted(row.items())]
            for pivot,row in sorted(basis.items()) if len(monomials[pivot])<=1]
    return {'prime':prime,'status':status,'seconds_limit':seconds,
            'elapsed_seconds':time.monotonic()-started,'generator_rows_total':len(generators),
            'generator_rows_fully_processed':rows_processed,'row_elimination_operations':operations,
            'echelon_rank':len(basis),'pivot_degree_counts':dict(sorted(degrees.items())),
            'basis_nonzero_entries':sum(len(row) for row in basis.values()),
            'linear_consequences':linear,
            'interpretation':'Rank and lower-degree consequences of the processed modular Macaulay rows only. '
                             'A partial run is not the degree-three Macaulay rank; no characteristic-zero fiber degree follows.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--prime',type=int,default=32003)
    parser.add_argument('--seconds',type=float,default=60)
    args=parser.parse_args()
    assert s.isprime(args.prime) and args.prime>5 and 0<args.seconds<=300
    equations,columns,base,source,certificate=load_fiber()
    rational=row_reduce(equations)
    assert len(rational)==134
    assert not any(len(mon)<=1 for mon in rational)
    result=modular_probe(equations,args.prime,args.seconds)
    report={'id':'S61_AFFINE_FIBER_MACAULAY_20261031','input_field':'QQ',
            'equations':134,'variables':45,'fixed_first_jet':'T01 log parameters (2,3,5,7,11,13)',
            'center':'the full rational T01 S61 point in S61_KNOWN_FIBER_AUDIT_20261031.json',
            'weight_two_u_indices_zero_based':columns,'center_coordinates':[str(x) for x in base],
            'QQ_constant_span_rank':134,'QQ_constant_span_linear_consequence_rank':0,
            'modular_degree_three_probe':result,
            'limitations':['No full fiber degree, dimension, radical, or exhaustion claim.',
                           'The rational no-linear statement concerns constant row operations only.',
                           'Modular degree-three data concerns one finite-field specialization.'],
            'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in (Path(__file__),source,certificate)}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report,out,ensure_ascii=False,indent=2);out.write('\n')
    print(json.dumps({k:result[k] for k in ('status','elapsed_seconds','generator_rows_fully_processed',
                                         'echelon_rank','pivot_degree_counts','basis_nonzero_entries')}))


if __name__=='__main__':main()

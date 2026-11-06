#!/usr/bin/env python3
"""Check declared complete modular sources and their full ambient tangents.

The input is a finite mathematical list of solved branches, not sample probes.
Each branch keeps all its parameters and declares one exact witness.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s
from classification_double_quadric import full_tangent_ranks, check_affine_poisson
from classification_multiple_hyperplane_nilpotent import PAIRS, jacobi, pf


def branch(source_report, spec):
    start=time.monotonic()
    chart=source_report['cases'][spec['chart_index']]
    count=chart['full_free_parameter_count']
    nn=s.symbols(f'n0:{count}')
    xx4=s.symbols('x1:5')
    local={str(x):x for x in (*nn,*xx4)}
    sub={local[k]:s.sympify(v,locals=local) for k,v in spec['substitution'].items()}
    params=tuple(local[k] for k in spec['free_parameters'])
    point_sub={local[k]:s.Rational(v) for k,v in spec['witness'].items()}
    assert set(point_sub)==set(params)
    entries=[s.factor(s.sympify(f,locals=local).subs(sub,simultaneous=True))
             for f in chart['full_six_brackets']]
    assert set().union(*(f.free_symbols for f in entries)) <= set(params+xx4)
    for f in chart['all_distinct_jacobi_coefficients']:
        assert s.factor(s.sympify(f,locals=local).subs(sub,simultaneous=True))==0
    Pa=s.zeros(4)
    for (i,j), f in zip(PAIRS,entries):
        Pa[i,j],Pa[j,i]=f,-f
    assert all(s.factor(f)==0 for f in jacobi(Pa,xx4))
    pfaffian=s.factor(pf(Pa))
    assert pfaffian != 0
    assert s.factor(pfaffian.subs(point_sub)) != 0
    h=s.symbols('h');xx=(h,*xx4)
    scale=s.ilcm(*[c.q for f in entries for c in s.Poly(f,*params,*xx4).coeffs()])
    P=s.zeros(5)
    for (i,j),f in zip(PAIRS,entries):
        P[i+1,j+1]=s.expand(scale*sum(a*h**(2-sum(e))*s.prod(x**v for x,v in zip(xx4,e))
                                          for e,a in s.Poly(f,*xx4).terms()))
        P[j+1,i+1]=-P[i+1,j+1]
    point=P.subs(point_sub)
    check_affine_poisson(point,xx)
    directions=[P.diff(p).subs(point_sub) for p in params]
    for tensor in (point,*directions):
        assert all(c.q==1 for f in tensor for c in s.Poly(f,*xx).coeffs())
    result=full_tangent_ranks(point,xx,directions)
    result.update(name=spec['name'],full_source_specification=spec,
                  universal_brackets={f'{i+1},{j+1}':str(P[i+1,j+1]) for i,j in PAIRS},
                  universal_pfaffian=str(pfaffian),all_parameter_raw_jacobi_zero=True,
                  all_source_equations_zero=True,integer_scale=int(scale),
                  witness_brackets={f'{i+1},{j+1}':str(point[i+1,j+1]) for i,j in PAIRS},
                  seconds=round(time.monotonic()-start,4))
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('input',type=Path);ap.add_argument('output',type=Path)
    args=ap.parse_args();data=json.loads(args.input.read_text())
    source_path=Path(data['source_report']);source=json.loads(source_path.read_text())
    results=[branch(source,spec) for spec in data['branches']]
    report={'status':'PASS','input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),
            'source_report_sha256':hashlib.sha256(source_path.read_bytes()).hexdigest(),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'branches':results,'completeness_scope':data['completeness_scope']}
    args.output.write_text(json.dumps(report,indent=2,default=int)+'\n')
    print(json.dumps({'status':'PASS','branches':[{k:r[k]for k in (
        'name','matrix_ranks_mod_prime','projective_tangent_dimension_upper',
        'projective_family_dimension_lower','seconds')}for r in results]}))

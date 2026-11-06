#!/usr/bin/env python3
"""Whole-modulus faithful S61 graph test for the p0 two-hyperplane cubic source."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time
import sympy as s
from s61_model import S61, X, pushforward
from classification_degenerate_quadric_q4_nil31_I_graph import darboux

ROOT=Path(__file__).resolve().parents[1]
SOURCE='reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_COMPONENTS_20261107.json'

def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    begun=time.monotonic();source=json.loads((ROOT/SOURCE).read_text())['branches'][0]
    old=s.symbols('h x1:5');lam=s.Symbol('lam');P=s.zeros(5)
    section={s.Symbol(k):(lam if k=='n3' else 0) for k in source['full_source_specification']['free_parameters']}
    for key,value in source['universal_brackets'].items():
        i,j=map(int,key.split(','));P[i,j]=s.sympify(value).subs(section,simultaneous=True).subs(dict(zip(old,X)),simultaneous=True);P[j,i]=-P[i,j]
    translation=s.eye(5);translation[1,0]=translation[2,0]=-1
    translated=pushforward(P,translation)
    value=translated.subs(dict(zip(X,(1,0,0,0,0))))[1:5,1:5]
    assert s.factor(value.det())==16
    frame=s.eye(5);frame[1:5,1:5]=darboux(value)
    framed=pushforward(translated,frame)
    model=S61();normalized=model.normalize(framed);model.verify_point(normalized,check_equations=True)
    point=normalized['parameters'];us=s.symbols('u1:62');results={}
    for name,path in [('T39','reports/Q51_UNIVERSAL_SECTION_20261101.json'),('T40','reports/Q52_UNIVERSAL_SECTION_20261101.json')]:
        graph=json.loads((ROOT/path).read_text());indices=graph['weight2_u_indices_zero_based']
        residual=[s.cancel(point[k]-s.sympify(q).subs(dict(zip(us,point)),simultaneous=True)) for k,q in zip(indices,graph['quadratic_section_q'])]
        results[name]=dict(whole_section_in_component=all(r==0 for r in residual),
            all45_residuals=list(map(str,residual)),nonzero_residuals={str(k):str(r) for k,r in zip(indices,residual) if r})
    paths=[SOURCE,'code/classification_two_hyperplanes_cubic_nilpotent_graph.py','code/s61_model.py',
           'code/classification_degenerate_quadric_q4_nil31_I_graph.py',
           'reports/Q51_UNIVERSAL_SECTION_20261101.json','reports/Q52_UNIVERSAL_SECTION_20261101.json']
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,source='j21_hx_h_p0',section=section.__repr__(),
        all61_parameters=list(map(str,point)),all134_original_equations_zero=True,full150_quotient_reconstruction_checked=True,
        translation=[[str(a) for a in row] for row in translation.tolist()],Darboux_frame=[[str(a) for a in row] for row in frame.tolist()],
        perspective=list(map(str,normalized['perspective'])),graph_membership=results,
        whole_source_coverage_requires_centralizer_dominance_proof=True,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun)
    with output.open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps(dict(graph_membership={k:v['whole_section_in_component']for k,v in results.items()},seconds=result['seconds'])))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output)

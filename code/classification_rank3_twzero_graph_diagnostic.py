#!/usr/bin/env python3
"""Faithful S61 normalization and whole-graph diagnostic at the Tw0 witness."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

import sympy as s
from s61_model import S61,X,pushforward

ROOT=Path(__file__).resolve().parents[1]
SOURCE='reports/CLASSIFICATION_RANK3_TWZERO_GENERIC_REVIEW_20261107.json'


def darboux(K):
    unit=[s.eye(4)[:,i]for i in range(4)]
    i,j=next((i,j)for i,j in itertools.combinations(range(4),2)if K[i,j])
    r0=unit[i];r1=unit[j]/K[i,j]
    remaining=[]
    for r in unit:
        q=r-(r.T*K*r1)[0]*r0+(r.T*K*r0)[0]*r1
        if q!=s.zeros(4,1) and (not remaining or s.Matrix.hstack(*remaining,q).rank()>len(remaining)):
            remaining.append(q)
    assert len(remaining)==2
    r2,r3=remaining;r3=r3/(r2.T*K*r3)[0]
    S=s.Matrix.vstack(*(r.T for r in (r0,r1,r2,r3)))
    J=s.zeros(4);J[0,1]=J[2,3]=1;J[1,0]=J[3,2]=-1
    assert S*K*S.T==J
    return S


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic();source=json.loads((ROOT/SOURCE).read_text())
    old=s.symbols('h u v z w');P=s.zeros(5)
    for ij in itertools.combinations(range(5),2):
        f=s.sympify(source['point_brackets'].get(str(ij),'0')).subs(dict(zip(old,X)),simultaneous=True)
        P[ij]=f;P[ij[::-1]]=-f
    translation=s.eye(5);translation[4,0]=-1
    translated=pushforward(P,translation)
    value=translated.subs(dict(zip(X,(1,0,0,0,0))))[1:5,1:5]
    assert value.det()!=0
    frame=s.eye(5);frame[1:5,1:5]=darboux(value)
    framed=pushforward(translated,frame)
    model=S61();normalized=model.normalize(framed)
    model.verify_point(normalized,check_equations=True)
    point=normalized['parameters'];us=s.symbols('u1:62')
    results={}
    for name,path in [('T39','reports/Q51_UNIVERSAL_SECTION_20261101.json'),
                      ('T40','reports/Q52_UNIVERSAL_SECTION_20261101.json')]:
        graph=json.loads((ROOT/path).read_text())
        two=graph['weight2_u_indices_zero_based']
        residual=[s.cancel(point[k]-s.sympify(q).subs(dict(zip(us,point)),simultaneous=True))
                  for k,q in zip(two,graph['quadratic_section_q'])]
        results[name]=dict(in_entire_reduced_component_graph=all(r==0 for r in residual),
                           all45_residuals=list(map(str,residual)),
                           nonzero_residuals={str(k):str(r)for k,r in zip(two,residual)if r})
    paths=[SOURCE,'code/classification_rank3_twzero_graph_diagnostic.py','code/s61_model.py',
           'reports/Q51_UNIVERSAL_SECTION_20261101.json','reports/Q52_UNIVERSAL_SECTION_20261101.json']
    report=dict(task_id='E_quadric3',all_passed=True,
        translation=[[str(t)for t in row]for row in translation.tolist()],
        nondegenerate_evaluation=[[str(t)for t in row]for row in value.tolist()],
        Darboux_frame=[[str(t)for t in row]for row in frame.tolist()],
        perspective=list(map(str,normalized['perspective'])),
        all61_parameters=list(map(str,point)),
        full150_quotient_reconstruction_checked=True,all134_equations_zero=True,
        graph_membership=results,
        whole_source_membership_not_claimed=True,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        seconds=time.monotonic()-begun)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(dict(graph_membership={k:v['in_entire_reduced_component_graph']for k,v in results.items()},seconds=report['seconds'])))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output)

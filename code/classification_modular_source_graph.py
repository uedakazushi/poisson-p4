#!/usr/bin/env python3
"""Whole Q51/Q52 graph recognition of a declared full modular-source witness."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import sympy as S
from s61_model import S61, X, pushforward
from classification_rank3_twzero_graph_diagnostic import darboux

ROOT=Path(__file__).resolve().parents[1]


def main(source_path,branch_name,output):
    began=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    signal.alarm(60)
    branch=next(b for b in json.loads(source_path.read_text())['branches'] if b['name']==branch_name)
    old=S.symbols('h x1 x2 x3 x4')
    env={str(v):v for v in old}
    P=S.zeros(5)
    for ij in itertools.combinations(range(5),2):
        entry=S.sympify(branch['witness_brackets'].get(f'{ij[0]},{ij[1]}','0'),locals=env)
        P[ij]=entry.subs(dict(zip(old,X)),simultaneous=True)
        P[ij[::-1]]=-P[ij]
    value=P.subs(dict(zip(X,(1,0,0,0,0))))[1:5,1:5]
    assert value.det()!=0
    frame=S.eye(5)
    frame[1:5,1:5]=darboux(value)
    framed=pushforward(P,frame)
    model=S61()
    normalized=model.normalize(framed)
    model.verify_point(normalized,check_equations=True)
    point=normalized['parameters']
    us=S.symbols('u1:62')
    graph_paths=[ROOT/'reports/Q51_UNIVERSAL_SECTION_20261101.json',
                 ROOT/'reports/Q52_UNIVERSAL_SECTION_20261101.json']
    results={}
    for name,path in zip(('T39','T40'),graph_paths):
        graph=json.loads(path.read_text())
        indices=graph['weight2_u_indices_zero_based']
        residual=[S.cancel(point[k]-S.sympify(q).subs(dict(zip(us,point)),simultaneous=True))
                  for k,q in zip(indices,graph['quadratic_section_q'])]
        assert len(residual)==45
        results[name]=dict(in_entire_reduced_component_graph=all(v==0 for v in residual),
                           all45_residuals=list(map(str,residual)),
                           nonzero_residuals={str(k):str(v) for k,v in zip(indices,residual) if v!=0})
    paths=[Path(__file__).resolve(),source_path.resolve(),ROOT/'code/s61_model.py',
           ROOT/'code/classification_rank3_twzero_graph_diagnostic.py',*graph_paths]
    result=dict(all_passed=True,task_id='E_multiple_hyperplane',branch=branch_name,
        original_origin_nondegenerate=True,nondegenerate_evaluation=[[str(v) for v in row] for row in value.tolist()],
        Darboux_frame=[[str(v) for v in row] for row in frame.tolist()],
        perspective=list(map(str,normalized['perspective'])),all61_parameters=list(map(str,point)),
        all150_reconstruction_checked=True,all134_equations_zero=True,
        graph_membership=results,point_result_requires_smoothness_to_transport_to_whole_source=True,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.monotonic()-began)
    with output.open('x') as handle:
        json.dump(result,handle,indent=2)
        handle.write('\n')
    signal.alarm(0)
    print(json.dumps(dict(all_passed=True,graphs={k:v['in_entire_reduced_component_graph'] for k,v in results.items()},
                          seconds=result['seconds'])),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_report',type=Path)
    parser.add_argument('branch')
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    main(args.source_report,args.branch,args.output)

#!/usr/bin/env python3
"""Actual projective flag fibers of the closed CA/CB resolutions at Q.

Generic affine fibers and resolutions are proved separately. Here the full
rank-drop flag fiber is computed, including its nonreduced scheme structure.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from s61_model import S61,X,pushforward
from weight1_structure import coefficient_map
ROOT=Path(__file__).resolve().parents[1]

def length_of_gb(gb):
    leads=[p.LM(order=gb.order).exponents for p in gb.polys];bounds=[]
    for i in range(3):
        pure=[e[i] for e in leads if e[i] and all(not e[j] for j in range(3) if j!=i)]
        if not pure:return None
        bounds.append(min(pure))
    return sum(not any(all(a>=b for a,b in zip(exp,lm)) for lm in leads)
               for exp in itertools.product(*(range(v) for v in bounds)))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic();path=ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json'
    body=json.loads(path.read_text());base=list(map(s.Symbol,body['base_parameters']));fiber=list(map(s.Symbol,body['fiber_parameters']))
    beta_indices=body['beta_indices_zero_based'];q_indices=body['q_indices_zero_based'];model=S61()
    uCA=s.Matrix(list(map(s.sympify,body['families']['CA']['normalized_61_coordinates'])))
    # The distinguished v24 Q, transported exactly to constant value J.
    point=dict(zip(base,[0,0,s.Rational(4,3),s.Rational(3,2),0,0]))|dict.fromkeys(fiber,0);point[fiber[1]]=s.Rational(9,2)
    uQ=uCA.subs(point);P=model.reconstruct(uQ);adapter=s.Matrix(body['beta16_to_intrinsic24_adapter'])
    beta24=adapter*uQ[beta_indices,0];K=coefficient_map(beta24);assert K.rank()==2
    assert K.nullspace()==[s.Matrix([0,1,0,0]),s.Matrix([0,0,0,1])]
    r,z,w=s.symbols('r z w');g=s.Matrix([[1,r,-w,z],[0,1,0,0],[0,z,1,0],[0,w,0,1]])
    full=s.eye(5);full[1:5,1:5]=g.inv();transformed=model.normalize(pushforward(P,full));u=s.Matrix(transformed['parameters'])
    constraints={};sources={}
    for name,family in body['families'].items():
        bleft=s.Matrix(family['beta_left_inverse']);beta_constraints=s.Matrix(family['beta_closed_constraints'])*u[beta_indices,0]
        params=bleft*u[beta_indices,0];shift=s.Matrix(list(map(s.sympify,family['q_shift']))).subs(dict(zip(base,params)),simultaneous=True)
        q_constraints=s.Matrix(family['fiber_closed_constraints'])*(u[q_indices,0]-shift)
        rows=[s.expand(v) for v in [*beta_constraints,*q_constraints] if s.expand(v)!=0]
        constraints[name]=rows;sources[name]=dict(unsimplified_constraints=list(map(str,rows)))
    constraints['N']=constraints['CA']+constraints['CB']
    for name,rows in constraints.items():
        gb=s.groebner(rows,r,z,w,order='lex',domain=s.QQ);length=length_of_gb(gb)
        nilpowers={}
        if length is not None:
            for variable in (r,z,w):
                nilpowers[str(variable)]=next(n for n in range(1,length+1) if gb.reduce(variable**n)[1]==0)
        sources[name]=sources.get(name,{})|dict(groebner_basis=list(map(str,gb)),quotient_length=length,
            reference_line_support_nilpotence_powers=nilpowers,flag_fiber_supported_at_reference_line=len(nilpowers)==3)
        print(name,'FLAG FIBER GB',list(gb),'LENGTH',length,flush=True)
    # K_Q forces a line into span(e_y,e_w); the sole point outside v_y!=0
    # is e_w. It fails even the base S6 condition, so no off-chart flag exists.
    swap=s.Matrix([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]]);off=s.eye(5);off[1:5,1:5]=swap
    uoff=s.Matrix(model.normalize(pushforward(P,off))['parameters']);bcon=s.Matrix(body['families']['CA']['beta_closed_constraints'])*uoff[beta_indices,0]
    assert any(v!=0 for v in bcon)
    # The exact rank3 open on the reference six-plane.
    B=s.Matrix(body['families']['CA']['beta_embedding']);Kb=coefficient_map(adapter*B*s.Matrix(base))
    minor=s.factor(Kb.extract([0,1,2],[0,2,3]).det());tau,_,_,C,D,F=base
    assert s.expand(minor-tau*(4*C*F-D*D+tau*tau/4))==0
    inputs=[Path(__file__),path,ROOT/'code/s61_model.py',ROOT/'code/weight1_structure.py']
    result=dict(id='AMBIENT_RANKDROP_FLAG_FIBERS_V2_20261103_v63',all_passed=True,Q_coordinates=list(map(str,uQ)),
        base_point={str(k):str(v) for k,v in point.items()},intrinsic_coefficient_matrix=K.tolist(),
        flag_chart_variables=list(map(str,[r,z,w])),flag_fibers=sources,no_off_chart_flags=True,
        off_chart_base_constraint_witness=list(map(str,bcon)),rank3_reference_minor=str(minor),
        rankdrop_scope='Exact projective source flag fibers at Q; generic component fibers/resolutions are established in the separate accepted note.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},seconds=time.monotonic()-started)
    def serial(v):
        if isinstance(v,s.Rational):return str(v)
        raise TypeError(type(v))
    (args.output/'report.json').write_text(json.dumps(result,indent=2,default=serial)+'\n');print('RANKDROP FLAGS SAVED',result['seconds'],flush=True)

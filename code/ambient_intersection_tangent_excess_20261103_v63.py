#!/usr/bin/env python3
"""Embedded component cotangents from exact pinching rings, not normalization df.

The full normalizer maps are frozen inputs. Regular source coordinates, p²,p³,
and the I*p generators give the entire cotangent space of each image ring.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    source=ROOT/'reports/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63/report.json'
    body=json.loads(source.read_text());base=json.loads((ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json').read_text())
    qi=base['q_indices_zero_based'];p,sv=s.symbols('p s');c0,c2,c5=s.symbols('c0 c2 c5');free=[c0,c2,c5]
    matrices={};records={}
    for name,record in body['families'].items():
        pars=list(map(s.Symbol,record['source_fiber_parameters']));u=s.Matrix(list(map(s.sympify,record['full_61_coordinates'])))
        q=u[qi,0];M=q.jacobian(pars);b=q.subs(dict.fromkeys(pars,0));selected=record['actual_q_selected_rows'];L=M[selected,:]
        Li=[L.diff(p,n).subs({p:0,sv:0})/s.factorial(n) for n in range(4)]
        bi=[b[selected,0].diff(p,n).subs({p:0,sv:0})/s.factorial(n) for n in range(4)]
        point={v:0 for v in pars};point[s.Symbol('k0')]=c0;point[s.Symbol('k2')]=c2;point[s.Symbol('k5')]=c5
        w0=s.Matrix([point[v] for v in pars]);ws=[w0]
        for n in range(1,4):
            ws.append(-Li[0].inv()*(bi[n]+sum((Li[i]*ws[n-i] for i in range(1,n+1)),s.zeros(len(pars),1))))
        curve=u.subs(dict(zip(pars,list(sum((p**n*ws[n] for n in range(4)),s.zeros(len(pars),1))))),simultaneous=True).subs(sv,0)
        columns=[curve.diff(p,n).subs(p,0)/s.factorial(n) for n in [2,3]]
        zero={p:0,sv:0}|point;regular_K=u.jacobian(pars).subs(zero)*Li[0].inv()
        w_s=-Li[0].inv()*(L.diff(sv).subs({p:0,sv:0})*w0+b[selected,0].diff(sv).subs({p:0,sv:0}))
        regular_s=u.diff(sv).subs(zero)+u.jacobian(pars).subs(zero)*w_s
        columns += [regular_K[:,i] for i in range(regular_K.cols)]+[regular_s]
        odd=s.Matrix(record['first_odd_coefficient_matrix'])*Li[0].inv();embedded_odd=s.zeros(61,len(pars))
        for i,row in enumerate(qi):embedded_odd[row,:]=odd[i,:]
        columns += [embedded_odd[:,i] for i in range(embedded_odd.cols)]
        wp=ws[1]
        w_sp=-Li[0].inv()*(L.diff(sv).diff(p).subs({p:0,sv:0})*w0+b[selected,0].diff(sv).diff(p).subs({p:0,sv:0})
            +L.diff(p).subs({p:0,sv:0})*w_s+L.diff(sv).subs({p:0,sv:0})*wp)
        source_curve=w0+wp*p+w_s*sv+w_sp*p*sv
        mixed=u.subs(dict(zip(pars,list(source_curve))),simultaneous=True).diff(sv).diff(p).subs({p:0,sv:0})
        columns.append(mixed);matrix=s.Matrix.hstack(*columns).applyfunc(s.expand)
        rank=matrix.rank();assert rank==(19 if name=='N' else 21)
        assert all(s.Poly(v,*free).total_degree()<=1 for v in matrix)
        matrices[name]=matrix
        records[name]=dict(beta_slice_embedded_tangent_matrix=[[str(v) for v in row] for row in matrix.tolist()],
            beta_slice_tangent_dimension=rank,whole_S61_tangent_dimension=rank+7,
            source_point={str(k):str(v) for k,v in point.items()},
            construction='Target regular coordinates +p²,+p³,+all mixed regular*p classes; these generate the exact pinching-ring cotangent, whereas normalization df alone omits the latter classes.')
        print(name,'EXACT SYMBOLIC EMBEDDED TANGENT',rank,time.monotonic()-started,flush=True)
    union=matrices['CA'].row_join(matrices['CB']);union_rank=union.rank();assert union_rank==22
    # Polynomial matrices: rank22 overQQ(c) implies all23-minors identically0
    # at EVERY specialization. The3 separate component ranks are fixed by
    # their proved complete local rings, including each specialfree value.
    assert matrices['N'].row_join(matrices['CA']).rank()==21
    assert matrices['N'].row_join(matrices['CB']).rank()==21
    for point in [(0,0,0),(1,0,0),(2,3,5)]:
        spoint=dict(zip(free,point));assert union.subs(spoint).rank()==22
    paths=[Path(__file__),source,ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json',
        ROOT/'research/notes/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63.md',
        ROOT/'reports/AMBIENT_WHOLE_REDUCED_INTERSECTION_20261103_v63/report.json']
    result=dict(id='AMBIENT_INTERSECTION_TANGENT_EXCESS_20261103_v63',all_passed=True,field='QQ(c0,c2,c5)',
        free_nonnormal_parameters=list(map(str,free)),families=records,
        tangent_union_symbolic_rank=union_rank,tangent_union_every_specialization_upper=22,
        intersection_tangent_every_specialization_lower=20,intersection_tangent_whole_S61_lower=27,
        reduced_N_tangent_whole_S61_dimension=26,
        conclusion='The actual scheme CA intersection CB is NONREDUCED along the common nonnormal locus over the criticalrank2 beta orbit, in particular at distinguished Q(a=1). This detects nonreducedness without identifying the whole intersection ideal or its nilpotence order.',
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('INTERSECTION TANGENT EXCESS SAVED',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)

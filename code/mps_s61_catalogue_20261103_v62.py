#!/usr/bin/env python3
"""Build exact general MPS families and rational S61 charts from stored diagrams.

Period coordinates are first converted by the linear four-dimensional Pfaffian
duality; they are not confused with the original biresidue entries. The initial
smoothing sum is checked against projective Jacobi before it is used. Diagrams
requiring higher corrections are retained as unresolved, not silently accepted.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s

from s61_model import S61, X, J, pushforward, jacobi_projective, coefficients

ROOT = Path(__file__).resolve().parents[1]
DIAGRAMS = ROOT/'legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json'
PAIRS4 = list(itertools.combinations(range(1, 5), 2))


def raw_family(record):
    periods = s.symbols(f'p0:{record["period_dimension"]}')
    smoothings = s.symbols(f't0:{record["edge_count"]}')
    six = [sum(p*s.Rational(row[i]) for p,row in zip(periods,record['period_basis']))
           for i in range(6)]
    B = s.zeros(4)
    for (i,j),value in zip(itertools.combinations(range(4),2),six):
        B[i,j],B[j,i] = value,-value
    pf = s.expand(B[0,1]*B[2,3]-B[0,2]*B[1,3]+B[0,3]*B[1,2])
    assert pf != 0
    A = (-pf*B.inv()).applyfunc(s.cancel)
    assert all(s.Poly(v,*periods).total_degree() <= 1 for v in A)
    P = s.zeros(5)
    for (i,j),value in zip(PAIRS4,[A[i,j] for i,j in itertools.combinations(range(4),2)]):
        P[i,j],P[j,i] = value*X[i]*X[j],-value*X[i]*X[j]
    for value,edge in zip(smoothings,record['edges']):
        i,j = edge['edge']
        monomial = s.prod(X[int(k)]**int(power) for k,power in edge['orders'].items())
        P[i,j] += value*monomial
        P[j,i] -= value*monomial
    return periods,smoothings,P,pf,A


def normalize_family(P):
    # Send (1:1:1:1:1) to e0, then use a rational symplectic pivot chart.
    move = s.eye(5)
    for i in range(1,5):
        move[i,0] = -1
    shifted = pushforward(P, move)
    point = {X[0]:1, **{x:0 for x in X[1:]}}
    C = shifted[1:5,1:5].subs(point).applyfunc(s.cancel)
    pivot = next((i,j) for i,j in itertools.combinations(range(4),2) if C[i,j] != 0)
    perm = [*pivot,*[i for i in range(4) if i not in pivot]]
    reorder = s.eye(4)[perm,:]
    changed = reorder*C*reorder.T
    a,b,c,d,e,f = [changed[i,j] for i,j in itertools.combinations(range(4),2)]
    delta = s.factor(a*f-b*e+c*d)
    frame = s.Matrix([[1,0,0,0],[0,1/a,0,0],
                      [d/a,-b/a,1,0],[e/delta,-c/delta,0,a/delta]])*reorder
    frame = frame.applyfunc(s.cancel)
    assert (frame*C*frame.T-J[1:5,1:5]).applyfunc(s.cancel) == s.zeros(4)
    homogeneous = s.eye(5)
    homogeneous[1:5,1:5] = frame
    based = pushforward(shifted,homogeneous)
    model = S61()
    normalized = model.normalize(based)
    return model,normalized,dict(point=[1]*5,pivot_indices=list(pivot),
                                value_pivot=str(s.factor(a)),
                                value_pfaffian=str(delta),
                                frame=[[str(v) for v in row] for row in frame.tolist()])


def serialize(P):
    return {f'{i},{j}':str(s.expand(P[i,j])) for i,j in itertools.combinations(range(5),2)
            if P[i,j] != 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--families',nargs='+',default=[f'T{i:02d}' for i in range(2,39)])
    parser.add_argument('--maps',action='store_true',help='Construct/verify all61 general rational coordinates')
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    diagrams = {r['id']:r for r in json.loads(DIAGRAMS.read_text())['diagrams']}
    records = []
    for name in args.families:
        start = time.monotonic()
        rec = diagrams[name]
        periods,smoothings,P,pf,A = raw_family(rec)
        jacobi = jacobi_projective(P)
        nonzero = [v for v in jacobi if v != 0]
        report = dict(component_id=name,graph_type=rec['graph_type'],
                      period_parameters=list(map(str,periods)),smoothing_parameters=list(map(str,smoothings)),
                      biresidue_pfaffian_open=str(pf),homogeneous_bivector=serialize(P),
                      projective_Jacobi_identically_zero=not nonzero,
                      nonzero_Jacobi_count=len(nonzero),
                      first_nonzero_Jacobi=str(nonzero[0]) if nonzero else None,
                      proposed_s61_dimension=10+len(periods)+len(smoothings),
                      expected_s61_dimension=rec['projective_component_dimension']-9,
                      scope='General period/edge family; component identification uses the stored MPS diagram and published local deformation theorem. Residual Sp4 action is necessary for a full S61 parametrization.')
        if args.maps and not nonzero:
            model,norm,frame = normalize_family(P)
            report['normalization'] = frame
            report['s61_parameters_u1_through_u61'] = list(map(str,norm['parameters']))
            report['perspective_parameters'] = list(map(str,norm['perspective']))
            report['verification'] = model.verify_point(norm)
        report['seconds'] = time.monotonic()-start
        (args.output/f'{name}.json').write_text(json.dumps(report,indent=2)+'\n')
        records.append(report)
        print(name,'RAW',not nonzero,'MAP',report.get('verification'),report['seconds'],flush=True)
    inputs=[Path(__file__),DIAGRAMS,ROOT/'code/s61_model.py']
    out=dict(id='MPS_S61_CATALOGUE_20261103_v62',
             all_attempted_raw_families_passed=all(r['projective_Jacobi_identically_zero'] for r in records),
             records=[dict(component_id=r['component_id'],raw_Jacobi=r['projective_Jacobi_identically_zero'],
                           s61_map_verified='verification' in r,seconds=r['seconds']) for r in records],
             input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
             new_components_claimed=0,global_upper_bound_proved=False)
    (args.output/'report.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    main()

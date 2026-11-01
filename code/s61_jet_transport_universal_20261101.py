#!/usr/bin/env python3
"""Certify universal affine-q, quadratic-beta primitive-jet transport.

This is an exact coordinate formula. It does not assert constancy of any
elliptic-curve invariant or membership of additional equations in S61.
"""
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_jet_basepoint_transport_20261101 import transport
ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    model=S61();weights=model.reduced['parameter_weights']
    point=[s.Symbol('u'+str(i+1)) if w==1 else 0 for i,w in enumerate(weights)]
    one,two,q,basis,result=transport(model,point);u=[point[i] for i in one]
    assert all(not e['coefficients'].free_symbols for e in result)
    assert all(s.Poly(v,*u).total_degree()<=2 for e in result for v in e['constant'])
    combined=s.Matrix.vstack(*(e['coefficients'] for e in result));assert combined.rank()==45
    report={'id':'S61_JET_TRANSPORT_UNIVERSAL_20261101','all_passed':True,
        'field':'QQ[u_weight1,q_weight2]','weight1_u_order':[i+1 for i in one],
        'weight2_u_indices_zero_based':two,'q_convention':'uncentered original weight2 coordinates',
        'directions':[{'direction':e['direction'],
                       'q_coefficient_matrix':[[str(v) for v in row] for row in e['coefficients'].tolist()],
                       'quadratic_beta_constant':list(map(str,e['constant']))} for e in result],
        'combined_constant_q_rank':45,
        'scope':'Universal basepoint-transport identity only; no proposed invariant-conservation equation is assumed or proved.',
        'source_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in (Path(__file__),ROOT/'code/s61_jet_basepoint_transport_20261101.py',ROOT/'code/s61_model.py',ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json')}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'combined_q_rank':45,'output':str(args.output)}))

if __name__=='__main__':main()

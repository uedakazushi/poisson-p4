#!/usr/bin/env python3
"""Explicit projective identity between CA's f-in-hW boundary and CB squares."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,raw as rawCA,canonical,vector
from ambient_CB_global_flag_model_20261103_v63 import raw as rawCB

ROOT=Path(__file__).resolve().parents[1]

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic();h,x,y,z,w=X
    alpha,beta=s.symbols('alpha beta',nonzero=True);coeff=s.symbols('k0:16')
    mon=[s.prod((h,x,z,w)[i] for i in ids) for ids in itertools.combinations_with_replacement(range(4),3) if any(i>=2 for i in ids)]
    H=sum(a*b for a,b in zip(coeff,mon));P=rawCA(h*(alpha*h+beta*x),0,H)
    change=s.Matrix([[alpha,beta,0,0,0],[1,0,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]])
    subs=dict(zip(X,change.inv()*s.Matrix(X)));transformed=(change*P.subs(subs,simultaneous=True)*change.T).applyfunc(s.cancel)
    new_H=s.cancel(H.subs(subs,simultaneous=True));expected=rawCB(-beta*x*x,0,new_H)
    Euler=s.Matrix(X);field=s.Matrix([0,0,beta*x,0,0]);vertical=Euler*field.T-field*Euler.T
    assert (transformed-expected-vertical).applyfunc(s.cancel)==s.zeros(5)
    assert canonical(transformed-expected).applyfunc(s.cancel)==s.zeros(5)
    paths=[Path(__file__),ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ROOT/'code/ambient_CB_global_flag_model_20261103_v63.py']
    body=dict(id='AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63',all_passed=True,field='QQ(alpha,beta)',
        source_CA='f=h*(alpha*h+beta*x),g=0,H arbitrary16',
        coordinate_change=[[str(v) for v in row] for row in change.tolist()],
        target_CB='f=-beta*x^2,g=0,H transformed by inverse coordinatechange',
        exact_Euler_vertical_difference='E wedge (beta*x*d_y)',
        all_raw25_difference_entries_identically_zero=True,all_canonical150_difference_entries_identically_zero=True,
        cubic_monomials=list(map(str,mon)),transformed_H=str(new_H),
        conclusion='Global N rank<=2 boundary equals the square-binary-F images in CB endpoints: g=0 N boundary = Bb_square, and f=0 CA boundary = Ba_square. Proper closure includes coincident/zero cases. This does not exhaust all globalCA/CB rank2 intersections.',
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(body,indent=2)+'\n');print('GLOBAL BOUNDARY SQUARE IDENTITY PASS',body['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)

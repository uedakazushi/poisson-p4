#!/usr/bin/env python3
"""Proof audit and independent linear algebra for the discriminant ruling."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
MAIN='TORIC_DISCRIMINANT_LINEAR_RULING_20261102'
STEM='TORIC_DISCRIMINANT_LINEAR_RULING_REVIEW_20261102'


def verify():
    primary=ROOT/'reports'/(MAIN+'.json')
    r=json.loads(primary.read_text())
    for path,digest in r['input_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    assert r['all_passed'] and r['ambient_beta_jacobian_generic_rank']==10
    assert r['ordinary_node_hessian_value']=='-78265'
    u=s.symbols('u0:52')
    fixed=[u[4],u[5],u[12],u[26],u[27],u[49]-2*u[14]]
    jet=[u[4],u[5],-2*u[26],-u[12],-u[27],-2*u[14]+u[49]]
    assert s.Matrix(fixed).jacobian(u).rref()[0]==s.Matrix(jet).jacobian(u).rref()[0]
    A,B,C,D,E=s.symbols('A B C D E')
    recovery=[2*A*u[13]+3*D*u[14], 2*A*u[35]-3*C*u[14],
              2*B*u[25]-D*u[51], B*(2*u[37]+u[48])+E*u[51]]
    M=s.Matrix(recovery).jacobian([A,B,C,D,E])
    assert M.rank()==4
    relation=M.nullspace()[0]
    normalized=(relation/relation[3]).applyfunc(s.cancel)
    expected=s.Matrix([-3*u[14]/(2*u[13]),u[51]/(2*u[25]),-u[35]/u[13],1,
                       -(2*u[37]+u[48])/(2*u[25])])
    assert normalized==expected
    paths=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),primary,
           ROOT/'research/notes'/(MAIN+'.md'),ROOT/'code'/(MAIN.lower()+'.py')]
    return dict(id=STEM,all_passed=True,reviewed=MAIN,
                independent_checks=['six fixed conditions have the critical-jet rowspace',
                                    'four recovery rows have a unique normalized projective kernel'],
                scope='Source/proof audit; raw Jacobian and Hessian sample use primary certificates.',
                primary_input_hashes_verified=len(r['input_sha256']),
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in paths})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=verify()
    with args.output.open('x') as stream:
        json.dump(report,stream,indent=2)
        stream.write('\n')
    print('PASS')

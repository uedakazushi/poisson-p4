#!/usr/bin/env python3
"""Replay a spatial-degree obstruction for the fixed SS20 contraction.

This does not rule out another Poisson family or component membership.
The input's complete linear elimination is adopted; the Jacobi obstruction
for its entire two-parameter A projection is recomputed without its producer.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();begun=time.monotonic()
    root=Path(__file__).resolve().parents[1]
    input_path=root/'reports/CLASSIFICATION_NIL22_QUADRATIC_SS_LIMIT_20261107.json'
    original=json.loads(input_path.read_text())
    q1,q2,x,y,lam,r19,r21=s.symbols('q1 q2 x y lam r19 r21')
    coords=(q1,q2,x,y);env={str(z):z for z in (*coords,lam,r19,r21)}
    A=s.zeros(4)
    for pair,text in original['linear_reduction']['exact_A'].items():
        i,j=map(int,pair.split(','));A[i,j]=s.sympify(text,locals=env);A[j,i]=-A[i,j]
    assert set().union(*(z.free_symbols for z in A))<={*coords,lam,r19,r21}
    H=s.diag(-1,-1,1,1)
    divergence=s.Matrix([sum(s.diff(A[i,j],coords[i])for i in range(4))for j in range(4)])
    assert (divergence-H*s.Matrix(coords)).applyfunc(s.expand)==s.zeros(4,1)
    J=[];records=[]
    for triple in itertools.combinations(range(4),3):
        i,j,k=triple
        expression=s.expand(sum(A[i,l]*s.diff(A[j,k],coords[l])+A[j,l]*s.diff(A[k,i],coords[l])+A[k,l]*s.diff(A[i,j],coords[l])for l in range(4)))
        J.append(expression)
        for exponent,coefficient in s.Poly(expression,*coords).terms():
            if coefficient!=0:records.append((triple,exponent,s.factor(coefficient)))
    selected=[]
    targets=(r19,r21,r19*r21+s.Rational(1,16))
    for target in targets:
        for triple,exponent,coefficient in records:
            quotient=s.cancel(coefficient/target)
            if quotient!=0 and quotient.free_symbols<={lam}:
                selected.append(dict(triple=list(triple),exponent=list(exponent),raw_coefficient=str(coefficient),nonzero_generic_multiplier=str(quotient),normalized_condition=str(target)))
                break
        else:raise AssertionError('Expected direct Jacobi coefficient was absent.')
    assert s.expand(16*targets[2]-16*r21*targets[0])==1
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,
        adopted_input_scope='Complete 76x44 linear elimination, whose A projection is exactly the displayed two-parameter family.',
        exact_full_A=original['linear_reduction']['exact_A'],
        recomputed_full_Jacobi_A=list(map(str,J)),direct_three_coefficient_obstruction=selected,
        unit_ideal_identity='1=16*(r19*r21+1/16)-16*r21*r19',
        spatial_degree_argument='For homogeneous quadratic A,B and constant C, the spatial-degree-three part of Jac(C+tau*A+tau^2*B) at tau^2 is Jac(A). The spatial-degree-one term [C,B] cannot cancel it.',
        no_prescribed_actual_SS20_contraction=True,
        scope='Only C+tau*A+tau^2*B with the specified H-weight-zero/divergence and Laurent endpoint conditions is ruled out. No nonmembership or higher-degree claim.',
        input_sha256={str(input_path.relative_to(root)):hashlib.sha256(input_path.read_bytes()).hexdigest(),str(Path(__file__).relative_to(root)):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        seconds=round(time.monotonic()-begun,6))
    with args.output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))


if __name__=='__main__':main()

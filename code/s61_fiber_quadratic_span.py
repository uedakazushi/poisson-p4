#!/usr/bin/env python3
"""Certify a generic obstruction to constant-row affine-linear elimination.

The 134 quadratic-leading parts of the full fixed-jet equations are linearly
independent at an exact toric jet. Their nonzero minor therefore gives a
nonempty open subset of W16 where no nonzero constant row combination is
affine-linear. This is not a statement about nonlinear ideal consequences.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import sympy as s

from s61_affine_fiber_macaulay import load_fiber
from weight1_rank_one_fiber import row_reduce

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    equations,columns,base,source,point_source=load_fiber()
    basis=row_reduce(equations)
    pivots=list(basis)
    assert len(pivots)==134 and all(len(mon)==2 for mon in pivots)
    matrix=s.Matrix([[poly.get(mon,0) for mon in pivots] for poly in equations])
    determinant=matrix.to_DM().det();assert determinant
    data=json.loads(source.read_text());weights=data['parameter_weights']
    equation_weights=[]
    for equation in data['equations']:
        degrees={sum(weights[i] for i in mon) for mon,_,_ in equation}
        assert len(degrees)==1
        equation_weights.append(degrees.pop())
    assert Counter(equation_weights)=={4:70,5:64}
    report={
        'id':'S61_FIBER_QUADRATIC_SPAN_20261031','all_passed':True,'field':'QQ',
        'equations':134,'weight_two_variables':45,
        'fixed_first_jet':'T01 log Poisson coefficients (2,3,5,7,11,13)',
        'pivot_monomials_in_45_variables_zero_based':[list(mon) for mon in pivots],
        'corresponding_u_indices_zero_based':[[columns[i] for i in mon] for mon in pivots],
        'minor_matrix':[[str(x) for x in row] for row in matrix.tolist()],
        'minor_determinant':str(determinant),
        'equation_weight_distribution':{'4':70,'5':64},
        'polynomial_degree_bound_for_minor_in_first_jet':64,
        'claim':'On the nonempty first-jet open set where this 134 by 134 quadratic-coefficient minor is nonzero, the 134 quadratic leading parts are independent. No nonzero constant combination of the equations is affine-linear.',
        'limitations':['This does not rule out affine-linear consequences with polynomial multipliers.',
                       'No fiber dimension, degree, radical, or global component count is determined.',
                       'This rank statement is compatible with the positive-dimensional unsaturated infinity locus.'],
        'proof_note':'research/notes/S61_AFFINE_FIBER_REDUCTION_20261031.md',
        'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in (Path(__file__),ROOT/'code/s61_affine_fiber_macaulay.py',source,point_source)},
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report,out,ensure_ascii=False,indent=2);out.write('\n')
    print(json.dumps({'all_passed':True,'quadratic_span_rank':134,'generic_constant_linear_consequences':0,
                      'output':str(args.output)}))


if __name__=='__main__':main()

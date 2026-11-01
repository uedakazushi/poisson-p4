#!/usr/bin/env python3
"""Two rank certificates excluding C0 minus R0 at generic beta.

The accompanying note proves the commuting-pair Jordan reduction. These rank
certificates concern original134 leading equations, so no extra G6 assumption
is needed for the exclusion.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
import sympy as s
from s61_infinity_flags import leading_equations_matrix
from s61_model import S61,X,J as J5
from weight1_structure import J,minor_witness,sp_basis

ROOT=Path(__file__).resolve().parents[1]


def run():
    model=S61();z=s.Matrix(X[1:]);cases=[]
    cartan1=s.diag(1,-1,0,0);cartan2=s.diag(0,0,1,-1)
    H=s.diag(-1,1,1,-1);N=s.zeros(4);N[1,2]=N[3,0]=1
    for label,A,B in [('cartan',cartan1,cartan2),('paired_jordan',H,N)]:
        assert A.T*J+J*A==s.zeros(4) and B.T*J+J*B==s.zeros(4)
        assert A*B==B*A
        aa=s.Matrix([0,*list(A*z)]);bb=s.Matrix([0,*list(B*z)])
        Q=aa*bb.T-bb*aa.T
        normalized=model.normalize(X[0]**2*J5+Q)
        point=normalized['parameters']
        leading=leading_equations_matrix(model,point)
        assert leading[:,0]==s.zeros(134,1)
        linear=leading[:,1:];witness=minor_witness(linear)
        assert witness['rank']==12
        # Both representatives really lie outside the constant-factor locus.
        vv=s.Matrix(s.symbols('v0:4'));affine=Q[1:,1:]
        tangent=[]
        for i,j,k in combinations(range(4),3):
            expr=s.expand(vv[i]*affine[j,k]-vv[j]*affine[i,k]+vv[k]*affine[i,j])
            tangent.extend(s.Poly(expr,*z).coeffs())
        factor_matrix=s.Matrix(tangent).jacobian(vv)
        factor_witness=minor_witness(factor_matrix);assert factor_witness['rank']==4
        regular=A+2*B if label=='cartan' else A+B
        centralizer=s.Matrix.hstack(*[s.Matrix(regular*C-C*regular).reshape(16,1) for C in sp_basis()])
        centralizer_witness=minor_witness(centralizer)
        assert centralizer_witness['rank']==8
        cases.append({'label':label,'A':[[str(v) for v in row] for row in A.tolist()],
            'B':[[str(v) for v in row] for row in B.tolist()],
            'S61_coordinates':list(map(str,point)),
            'beta_leading_matrix':[[str(v) for v in row] for row in linear.tolist()],
            'beta_rank_minor':witness,'constant_factor_minor':factor_witness,
            'regular_element':[[str(v) for v in row] for row in regular.tolist()],
            'centralizer_rank_minor':centralizer_witness,
            'crude_incidence_dimension_bound':14})
    inputs=[Path(__file__),ROOT/'code/s61_infinity_flags.py',ROOT/'code/s61_model.py',
            model.source/'symplectic_slice_reduced_equations.json',
            ROOT/'legacy/v07/proofs/zero_fiber_residual_reduction.txt']
    return {'id':'S61_C0_INFINITY_20261101','all_passed':True,'field':'QQ','cases':cases,
        'scope':'For sufficiently general beta, the C0 part of original134 unsaturated '
                'infinity is contained in R0, conditional only on the stated proper '
                'commuting-pair description and the elementary Jordan proof in the note. '
                'No exclusion of L0 or nonzero-modular Q45 strata is claimed.',
        'proof_note':'research/notes/S61_C0_INFINITY_20261101.md',
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();report=run();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'beta_ranks':[12,12],'incidence_bounds':[14,14],'output':str(args.output)}))


if __name__=='__main__':main()

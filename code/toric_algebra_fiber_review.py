#!/usr/bin/env python3
"""Independent tensor-level review of the 20-variable toric algebra model.

Constructs general symmetric products directly, without importing the model's
splitting or associativity functions. Checks polynomial identities proving the
omitted upper commutator equations follow from the 96 lower equations.
"""
from __future__ import annotations

import argparse
import hashlib
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
J = s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
PAIRS = list(combinations(range(4),2))
SYM_PAIRS = list(combinations_with_replacement(range(4),2))
I = s.eye(4)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    variables = s.symbols('a0:40')
    tensor = {(i,j,k): variables[4*n+k] for n,(i,j) in enumerate(SYM_PAIRS) for k in range(4)}
    matrices = [s.Matrix(4,4,lambda k,j:tensor[(*sorted((i,j)),k)]) for i in range(4)]
    trace = s.Matrix([s.trace(m) for m in matrices])
    trace_map = trace.jacobian(variables)
    assert trace_map.rank() == 4
    beta = [J*m+m.T*J for m in matrices]
    flat_beta = s.Matrix([beta[k][i,j] for i,j in PAIRS for k in range(4)])
    linear_map = flat_beta.jacobian(variables)
    assert linear_map.rank() == 20
    assert all(s.expand(beta[i][j,k]-beta[j][i,k]+beta[k][i,j]) == 0
               for i,j,k in combinations(range(4),3))
    assert s.Matrix([b[0,1]+b[2,3] for b in beta]) == trace
    trace_free_basis = s.Matrix.hstack(*trace_map.nullspace())
    assert trace_free_basis.cols == 36
    assert (linear_map*trace_free_basis).rank() == 16
    kernel = linear_map.nullspace()
    assert len(kernel) == 20
    # Lower with J; every element of the kernel is totally symmetric.
    for vector in kernel:
        mm = [m.subs(dict(zip(variables,vector,strict=True))) for m in matrices]
        tt = lambda i,j,k:(mm[i][:,j].T*J)[k]
        assert all(tt(i,j,k)==tt(i,k,j) for i in range(4) for j in range(4) for k in range(4))

    # Check a universal syzygy before imposing trace-free or q=tr(mm)/3.
    qvariables = s.symbols('q0:10')
    qentries = dict(zip(SYM_PAIRS,qvariables,strict=True))
    q = s.Matrix(4,4,lambda i,j:qentries[tuple(sorted((i,j)))])
    c = lambda i,j: matrices[i]*matrices[j]-matrices[j]*matrices[i]+I[:,i]*q[j,:]-I[:,j]*q[i,:]
    d = lambda i,j: q[i,:]*matrices[j]-q[j,:]*matrices[i]
    identities = 0
    for i,j,k in combinations(range(4),3):
        lhs = I[:,i]*d(j,k)+I[:,j]*d(k,i)+I[:,k]*d(i,j)
        rhs = s.zeros(4)
        for a,b,e in ((i,j,k),(j,k,i),(k,i,j)):
            cc = c(b,e)
            rhs -= matrices[a]*cc-cc*matrices[a]
        assert (lhs-rhs).applyfunc(s.expand) == s.zeros(4)
        identities += 16
    # General trace-free projection, still parameterized by all 40 inputs.
    mm = [(matrices[i]-(trace[i]*I+I[:,i]*trace.T)/5).applyfunc(s.expand) for i in range(4)]
    assert all(s.trace(m)==0 for m in mm)
    assert all(mm[i][:,j]==mm[j][:,i] for i in range(4) for j in range(4))
    qt = s.Matrix(4,4,lambda i,j:s.expand(s.trace(mm[i]*mm[j])/3))
    regular = [s.Matrix.vstack(s.Matrix([[0]]).row_join(qt[i,:]),I[:,i].row_join(mm[i])) for i in range(4)]
    assert all(s.trace(m)==0 for m in regular)
    assert all(s.expand(s.trace(regular[i]*regular[j])-5*qt[i,j])==0
               for i in range(4) for j in range(4))

    # Independently build the split algebra by five character evaluations.
    # This uses no coframe or model implementation functions.
    a = s.Matrix([[1,2,0,0],[0,1,3,0],[0,0,1,2],[0,0,0,1]])
    ell0 = -s.ones(1,4)*a/5
    ell = s.Matrix.vstack(ell0,*(ell0+a[i,:] for i in range(4)))
    character = s.ones(5,1).row_join(ell)
    assert character.det() != 0
    regular_split = [character.inv()*s.diag(*ell[:,i])*character for i in range(4)]
    split_m = [r[1:,1:] for r in regular_split]
    split_q = s.Matrix.vstack(*(r[0,1:] for r in regular_split))
    assert split_q == s.Matrix(4,4,lambda i,j:s.trace(split_m[i]*split_m[j])/3)
    assert split_q.det() != 0
    assert all(r[0,0]==0 and r[1:,0]==I[:,i] for i,r in enumerate(regular_split))
    y = s.Matrix(s.symbols('y0:4'))
    frame = I+sum((y[i]*split_m[i] for i in range(4)),s.zeros(4))-y*(y.T*split_q)
    z = a*y
    h0 = 1-sum(z)/5
    coframe = a.inv()*(I-z*s.ones(1,4)/5)*s.diag(*(h0+v for v in z))*a
    assert (frame-coframe).applyfunc(s.expand)==s.zeros(4)
    norm = (s.eye(5)+sum((y[i]*r for i,r in enumerate(regular_split)),s.zeros(5))).det(method='domain-ge')
    assert s.expand(norm-s.prod(1+v for v in ell*y))==0
    assert s.expand(frame.det(method='domain-ge')-norm)==0

    report = {
        'all_passed':True, 'field':'QQ', 'id':'TORIC_ALGEBRA_FIBER_REVIEW_20261031',
        'general_symmetric_product_variables':40,
        'trace_constraint_rank':4, 'trace_free_dimension':36,
        'closed_jet_map_rank':20, 'primitive_jet_map_rank':16,
        'totally_symmetric_kernel_dimension':20,
        'universal_lower_to_upper_syzygy_entries_checked':identities,
        'syzygy':'e_i D_jk + e_j D_ki + e_k D_ij = -([m_i,C_jk]+[m_j,C_ki]+[m_k,C_ij]); distinct i,j,k imply D=0 when C=0',
        'regular_trace_gram':'diag(5,5q)',
        'character_constructed_split_algebra_checked':True,
        'full_frame_reconstruction_checked':True,
        'norm_and_five_hyperplanes_agree':True,
        'review_note':'research/notes/TORIC_ALGEBRA_FIBER_REVIEW_20261031.md',
        'limitations':['The semisimple-to-toric bijection is justified in the proof note, not inferred from the single split algebra example.',
                       'No general solution count or classification of det(q)=0.',
                       'No classification of non-toric S61 extensions.'],
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'reviewed_model_sha256':hashlib.sha256((ROOT/'code/toric_algebra_fiber.py').read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report,stream,ensure_ascii=False,indent=2);stream.write('\n')
    print(json.dumps({'all_passed':True,'universal_syzygy_entries':identities,'output':str(args.output)}))


if __name__=='__main__':
    main()

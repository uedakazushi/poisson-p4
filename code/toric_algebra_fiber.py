#!/usr/bin/env python3
"""A label-free quadratic inverse problem for the transverse toric first jet.

For fixed beta, trace-free commutative multiplications inducing beta form an
affine space modelled on Sym^3(V*) (20 variables). Associativity reduces to 96
quadrics. The open determinant(q)!=0 is exactly the semisimple five-character
case. The generic number of solutions is not computed by this program.
"""
from __future__ import annotations

import argparse
import hashlib
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path

import sympy as s

from weight1_structure import COORDS, J, constraints, form_matrices, minor_witness
from toric_coframe_fiber import Y, affine_bivector, jet, linear_frame

ROOT = Path(__file__).resolve().parents[1]
TRIPLES = list(combinations_with_replacement(range(4),3))


def multiplication_from_jet(vector, cubic):
    """m_i has j-th column m(e_i,e_j); T_ijk=<m(e_i,e_j),e_k>_J."""
    beta = form_matrices(vector)
    t = lambda i,j,k: (beta[i][j,k]+beta[j][i,k])/3 + cubic[tuple(sorted((i,j,k)))]
    return [s.Matrix.hstack(*(J*s.Matrix([t(i,j,k) for k in range(4)])
                              for j in range(4))) for i in range(4)]


def jet_from_multiplication(matrices):
    return s.Matrix([(J*matrices[k]+matrices[k].T*J)[i,j] for (i,j),k in COORDS])


def trace_form(matrices):
    return s.Matrix(4,4,lambda i,j:s.expand(s.trace(matrices[i]*matrices[j])/3))


def associativity_quadrics(matrices):
    q=trace_form(matrices)
    equations=[]
    for i,j in combinations(range(4),2):
        commutator = matrices[i]*matrices[j]-matrices[j]*matrices[i]
        commutator += s.eye(4)[:,i]*q[j,:]-s.eye(4)[:,j]*q[i,:]
        equations.extend(s.expand(x) for x in commutator)
    return q,equations


def coframe_multiplication(r):
    linear=r.inv()*linear_frame(r*Y)*r
    return [linear.applyfunc(lambda x:s.expand(x).coeff(y)) for y in Y]


def full_multiplication_matrices(matrices,q):
    # L_(e_i) for the five-dimensional algebra C*1 + V.
    return [s.Matrix.vstack(s.Matrix([[0]]).row_join(q[i,:]),
                             s.eye(4)[:,i].row_join(matrices[i])) for i in range(4)]


def verify_linear_reconstruction():
    constraint=s.Matrix.vstack(*constraints())
    basis=constraint.nullspace()
    zero={t:s.S.Zero for t in TRIPLES}
    for vector in basis:
        matrices=multiplication_from_jet(vector,zero)
        assert jet_from_multiplication(matrices)==vector
        assert all(s.trace(m)==0 for m in matrices)
        assert all(matrices[i][:,j]==matrices[j][:,i] for i in range(4) for j in range(4))
    kernel=[]
    for triple in TRIPLES:
        cubic=zero|{triple:s.S.One}
        matrices=multiplication_from_jet(s.zeros(24,1),cubic)
        assert jet_from_multiplication(matrices)==s.zeros(24,1)
        assert all(s.trace(m)==0 for m in matrices)
        kernel.append(s.Matrix.vstack(*(s.Matrix(list(m)) for m in matrices)))
    assert s.Matrix.hstack(*kernel).rank()==20
    return {'primitive_jet_basis_vectors_checked':16,'symmetric_cubic_kernel_rank':20,
            'trace_free_commutative_multiplication_dimension':36}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    linear_checks=verify_linear_reconstruction()
    # The split algebra at the canonical coframe is verified before any sample.
    canonical=coframe_multiplication(s.eye(4))
    q0,equations0=associativity_quadrics(canonical)
    assert equations0==[0]*96
    assert q0==(5*s.eye(4)-s.ones(4))/25
    assert q0.det()==s.Rational(1,3125)
    regular=full_multiplication_matrices(canonical,q0)
    assert all(a*b-b*a==s.zeros(5) for a,b in combinations(regular,2))
    # Fixed jet from the already-certified coframe witness.
    coframe_report=json.loads((ROOT/'reports/TORIC_COFRAME_FIBER_20261031.json').read_text())
    r=s.Matrix([[s.Rational(x) for x in row] for row in coframe_report['sample_coframe']])
    vector=jet(r)
    known=coframe_multiplication(r)
    coefficients=s.symbols('c0:20')
    cubic=dict(zip(TRIPLES,coefficients,strict=True))
    matrices=multiplication_from_jet(vector,cubic)
    q,equations=associativity_quadrics(matrices)
    assert len(equations)==96
    assert all(s.Poly(eq,*coefficients).total_degree()<=2 for eq in equations)
    base=multiplication_from_jet(vector,{t:s.S.Zero for t in TRIPLES})
    known_cubic={}
    for i,j,k in TRIPLES:
        known_cubic[tuple((i,j,k))]=((known[i][:,j]-base[i][:,j]).T*J)[k]
    known_sub=dict(zip(coefficients,[known_cubic[t] for t in TRIPLES],strict=True))
    assert all(m.subs(known_sub)==n for m,n in zip(matrices,known,strict=True))
    assert all(eq.subs(known_sub)==0 for eq in equations)
    known_q=q.subs(known_sub)
    assert known_q==r.T*q0*r
    assert known_q.det()==r.det()**2/3125 != 0
    my=sum((Y[i]*known[i] for i in range(4)),s.zeros(4))
    frame=s.eye(4)+my-Y*(Y.T*known_q)
    p=(frame*J*frame.T).applyfunc(s.expand)
    assert p==affine_bivector(r)
    derivative=s.Matrix(equations).jacobian(coefficients).subs(known_sub)
    witness=minor_witness(derivative)
    assert witness['rank']==20
    # Determinant of regular multiplication is the actual affine quintic.
    full=s.eye(5)+sum((Y[i]*a for i,a in enumerate(full_multiplication_matrices(known,known_q))),s.zeros(5))
    z=r*Y;h0=1-sum(z)/5
    assert s.expand(full.det(method='domain-ge')-h0*s.prod(h0+x for x in z))==0
    report={
        'id':'TORIC_ALGEBRA_FIBER_20261031','all_passed':True,'field':'QQ',
        'linear_checks':linear_checks,'number_of_unknowns':20,'number_of_equations':96,
        'equation_degrees':sorted({s.Poly(eq,*coefficients).total_degree() for eq in equations}),
        'unknown_order':[{'variable':str(c),'symmetric_tensor_index_zero_based':list(t)}
                         for t,c in zip(TRIPLES,coefficients,strict=True)],
        'fixed_first_jet':[str(x) for x in vector],
        'multiplication_matrices':[[[str(x) for x in row] for row in m.tolist()] for m in matrices],
        'trace_form_matrix':[[str(x) for x in row] for row in q.tolist()],
        'quadratic_equations':[str(eq) for eq in equations],
        'open_condition':'det(trace_form_matrix) != 0',
        'known_solution':[str(known_sub[c]) for c in coefficients],
        'known_trace_form_determinant':str(known_q.det()),
        'known_solution_jacobian_witness':witness,
        'canonical_trace_form':[[str(x) for x in row] for row in q0.tolist()],
        'norm_and_coframe_formula_agree':True,
        'claim':'For the fixed primitive jet, the solutions of the 96 quadrics with det(q)!=0 are in bijection with transverse toric Poisson extensions; there is no remaining S5 label quotient.',
        'limitations':['The number of solutions and generic degree are not computed.',
                       'No claim that all full S61 extensions are toric.',
                       'The det(q)=0 locus is not classified or identified with the toric boundary.'],
        'proof_note':'research/notes/TORIC_ALGEBRA_FIBER_20261031.md',
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in (Path(__file__),ROOT/'code/toric_coframe_fiber.py',ROOT/'code/weight1_structure.py')},
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report,out,ensure_ascii=False,indent=2);out.write('\n')
    print(json.dumps({'all_passed':True,'quadrics':96,'variables':20,'rank':20,
                      'minor':witness['determinant'],'output':str(args.output)}))


if __name__=='__main__':
    main()

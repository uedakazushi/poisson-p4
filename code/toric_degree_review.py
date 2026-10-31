#!/usr/bin/env python3
"""Independent exact elimination audit of the proposed toric generic degree.

Starts with the previously audited multiplication matrices, reconstructs all
commutators, and performs its own nullspace elimination and Groebner basis.
No routines from the degree certificate generator are imported.
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
I = s.eye(4)


def equations(m):
    q = s.Matrix(4,4,lambda i,j:s.expand(s.trace(m[i]*m[j])/3))
    eq = []
    for i,j in combinations(range(4),2):
        eq.extend(s.expand(v) for v in m[i]*m[j]-m[j]*m[i]+I[:,i]*q[j,:]-I[:,j]*q[i,:])
    return q,eq


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    args = parser.parse_args()
    source = ROOT/'reports/TORIC_ALGEBRA_FIBER_20261031.json'
    data = json.loads(source.read_text())
    c = s.symbols('c0:20')
    m = [s.Matrix([[s.sympify(v) for v in row] for row in block])
         for block in data['multiplication_matrices']]
    q,eq = equations(m)
    assert eq == [s.sympify(v) for v in data['quadratic_equations']]
    assert q == s.Matrix(data['trace_form_matrix'])
    # Rebuild the beta-independent homogeneous part from a symmetric tensor.
    triples = list(combinations_with_replacement(range(4),3))
    cubic = dict(zip(triples,c,strict=True))
    mc = [s.Matrix.hstack(*(J*s.Matrix([cubic[tuple(sorted((i,j,k)))] for k in range(4)])
                             for j in range(4))) for i in range(4)]
    assert all((a-b).jacobian(c)==s.zeros(16,20)
               for a,b in ((s.Matrix(list(mi)),s.Matrix(list(ni))) for mi,ni in zip(m,mc,strict=True)))
    _,eqc = equations(mc)
    quadratic_monomials = [c[i]*c[j] for i,j in combinations_with_replacement(range(20),2)]
    monomials = quadratic_monomials+list(c)+[s.S.One]
    polynomials = [s.Poly(v,*c) for v in eq]
    coefficient = s.Matrix([[p.coeff_monomial(t) for t in monomials] for p in polynomials])
    homogeneous = coefficient[:,:210]
    assert homogeneous == s.Matrix([[s.Poly(v,*c).coeff_monomial(t) for t in quadratic_monomials] for v in eqc])
    print('reconstructed homogeneous coefficient matrix',flush=True)
    nullspace = homogeneous.T.nullspace()
    assert len(nullspace)==51
    nullrows = s.Matrix.vstack(*(v.T for v in nullspace))
    lower = nullrows*coefficient[:,210:]
    row_indices = list(lower[:,:20].T.rref()[1])
    assert len(row_indices)==18
    linear = lower[row_indices,:]
    columns = list(linear[:,:20].rref()[1])
    assert len(columns)==18
    free_indices = [i for i in range(20) if i not in columns]
    free = [c[i] for i in free_indices]
    pivot = linear[:,columns]
    pivot_determinant = pivot.det(method='domain-ge')
    assert pivot_determinant!=0
    solutions = -pivot.inv()*(linear[:,free_indices]*s.Matrix(free)+linear[:,20])
    substitution = dict(zip((c[i] for i in columns),solutions,strict=True))
    assert all(s.expand(v)==0 for v in linear*s.Matrix([substitution.get(ci,ci) for ci in c]+[1]))
    print('eliminated eighteen variables',free_indices,flush=True)
    residual = [s.expand(v.subs(substitution,simultaneous=True)) for v in eq]
    x,y = free
    residual_coeff = s.Matrix([[s.Poly(v,x,y).coeff_monomial(t) for t in (x*x,x*y,y*y,x,y,1)] for v in residual])
    conic_indices = list(residual_coeff[:,:3].T.rref()[1])
    assert len(conic_indices)==3
    conic_det = residual_coeff[conic_indices,:3].det()
    assert conic_det!=0
    assert residual_coeff.rank()==3
    # Row-span equality makes the selected three conics generate all residuals.
    selected = [residual[i] for i in conic_indices]
    selected_coeff = residual_coeff[conic_indices,:]
    reconstruction = residual_coeff[:,:3]*selected_coeff[:,:3].inv()
    assert reconstruction*selected_coeff==residual_coeff
    basis = s.groebner(selected,x,y,order='lex',domain=s.QQ)
    assert all(basis.reduce(v)[1]==0 for v in residual)
    assert len(basis.polys)==2
    first,second = [p.as_expr() for p in basis.polys]
    assert s.Poly(first,x).degree()==1 and first.coeff(x)==1 and not second.has(x)
    assert s.Poly(second,y).degree()==3
    assert s.gcd(second,s.diff(second,y))==1
    q_reduced = q.applyfunc(lambda v:s.expand(v.subs(substitution,simultaneous=True)))
    qdet = q_reduced.det(method='domain-ge')
    qdet_mod = basis.reduce(qdet)[1]
    assert s.gcd(qdet_mod,second)==1
    # Upper-block equations are also zero in the resulting exact quotient.
    m_reduced = [mi.applyfunc(lambda v:s.expand(v.subs(substitution,simultaneous=True))) for mi in m]
    assert all(basis.reduce(s.expand(v))[1]==0 for i,j in combinations(range(4),2)
               for v in q_reduced[i,:]*m_reduced[j]-q_reduced[j,:]*m_reduced[i])
    # Audit the other implementation's saved certificate, rather than only
    # compare its final answer with this independent dense-matrix computation.
    certificate_path = ROOT/'reports/TORIC_ALGEBRA_DEGREE_20261031.json'
    certificate = json.loads(certificate_path.read_text())
    combinations_matrix = s.zeros(18,96)
    for row,item in enumerate(certificate['universal_linear_consequences']):
        for index,value in item['constant_combination_of_96_equations']:
            combinations_matrix[row,index] = s.Rational(value)
    assert combinations_matrix*homogeneous == s.zeros(18,210)
    certified_lower = combinations_matrix*coefficient[:,210:]
    certified_equations = certified_lower*s.Matrix(list(c)+[1])
    assert all(s.expand(value-s.sympify(item['equation_at_sample']))==0
               for value,item in zip(certified_equations,certificate['universal_linear_consequences'],strict=True))
    certificate_columns = certificate['linear_pivot_columns_zero_based']
    assert certified_lower[:,certificate_columns] == s.Matrix(certificate['linear_pivot_matrix_at_sample'])
    assert certified_lower[:,certificate_columns].det(method='domain-ge') == s.Rational(certificate['linear_pivot_determinant_at_sample'])
    assert all(s.expand(substitution.get(ci,ci)-s.sympify(value))==0
               for ci,value in zip(c,certificate['eliminated_c0_through_c19'],strict=True))
    assert residual == [s.sympify(v) for v in certificate['residual_quadrics']]
    assert conic_det == s.Rational(certificate['residual_leading_determinant'])
    assert [first,second] == [s.sympify(v) for v in certificate['lex_groebner_basis']]
    assert qdet_mod == s.sympify(certificate['trace_determinant_mod_cubic'])
    bezout_q,bezout_cubic = [s.sympify(v) for v in certificate['trace_determinant_bezout']]
    assert s.expand(bezout_q*qdet_mod+bezout_cubic*second)==1
    normal = (-selected_coeff[:,:3].inv()*selected_coeff[:,3:]*s.Matrix([x,y,1])).applyfunc(s.expand)
    assert list(normal) == [s.sympify(v) for v in certificate['quadratic_normal_forms']]
    vector = lambda v:s.Matrix([s.Poly(v,x,y).coeff_monomial(t) for t in (1,x,y)])
    mx = s.Matrix.hstack(s.Matrix([0,1,0]),vector(normal[0]),vector(normal[1]))
    my = s.Matrix.hstack(s.Matrix([0,0,1]),vector(normal[1]),vector(normal[2]))
    assert mx == s.Matrix(certificate['multiplication_by_x'])
    assert my == s.Matrix(certificate['multiplication_by_y'])
    assert mx*my == my*mx
    print('verified reduced length-three sample, all semisimple',flush=True)
    report = {
        'id':'TORIC_DEGREE_REVIEW_20261031','all_passed':True,'field':'QQ',
        'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'reviewed_certificate_sha256':hashlib.sha256(certificate_path.read_bytes()).hexdigest(),
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'homogeneous_quadratic_rank':45,'constant_left_kernel_dimension':51,
        'independent_linear_consequences':18,
        'linear_consequence_nullspace_indices_zero_based':row_indices,
        'eliminated_variable_indices_zero_based':columns,
        'free_variable_indices_zero_based':free_indices,
        'linear_pivot_determinant':str(pivot_determinant),
        'selected_conic_equation_indices_zero_based':conic_indices,
        'conic_leading_determinant':str(conic_det),
        'all_residuals_in_selected_conic_row_span':True,
        'sample_lex_groebner_basis':[str(p) for p in (first,second)],
        'sample_length':3,'sample_reduced':True,
        'trace_form_determinant_mod_groebner':str(qdet_mod),
        'trace_form_determinant_coprime_to_eliminant':True,
        'upper_commutator_blocks_vanish_in_quotient':True,
        'stored_linear_combination_certificate_independently_checked':True,
        'stored_quotient_multiplication_matrices_independently_checked':True,
        'stored_trace_form_bezout_identity_independently_checked':True,
        'generic_degree_deduction':'Constant quadratic-part annihilators give universal linear consequences. Localize the nonzero pivot and conic-leading determinants to obtain a finite module generated by 1,x,y. The length-three reduced semisimple sample and smooth GL4/S5 structure give etaleness over a smaller nonempty base open, hence degree three.',
        'limitations':['Degree three is for the transverse toric open, not the full S61 fiber.',
                       'The generic elliptic degrees and full component exhaustion are not computed.'],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report,stream,ensure_ascii=False,indent=2);stream.write('\n')
    print(json.dumps({'all_passed':True,'sample_length':3,'free_variables':[str(v) for v in free],
                      'output':str(args.output)}))


if __name__=='__main__':
    main()

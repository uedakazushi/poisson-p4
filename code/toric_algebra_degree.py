#!/usr/bin/env python3
"""Certify generic toric first-jet degree three via a finite local model.

Tracks constant row operations giving universal linear consequences, then
certifies an exact three-dimensional special-fiber algebra with no discarded
points. The accompanying note supplies the finite/etale spreading argument;
a special-fiber degree alone would not prove the generic degree.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as s

from toric_algebra_fiber import TRIPLES, associativity_quadrics, multiplication_from_jet

ROOT = Path(__file__).resolve().parents[1]


def sparse(expression, variables):
    return {tuple(i for i,n in enumerate(powers) for _ in range(n)):Fraction(value)
            for powers,value in s.Poly(expression,*variables).terms() if value}


def add_multiple(target, source, scalar):
    for key,value in source.items():
        updated=target.get(key,Fraction(0))+scalar*value
        if updated: target[key]=updated
        else: target.pop(key,None)


def tracked_rows(polynomials):
    basis={}
    for index,polynomial in enumerate(polynomials):
        row=dict(polynomial); combination={index:Fraction(1)}
        while row:
            pivot=min(row,key=lambda mon:(-len(mon),mon))
            if pivot not in basis:
                value=row[pivot]
                basis[pivot]=({mon:c/value for mon,c in row.items()},
                              {j:c/value for j,c in combination.items()})
                break
            old,certificate=basis[pivot]
            value=-row[pivot]
            add_multiple(row,old,value); add_multiple(combination,certificate,value)
    return basis


def expression(row,variables):
    return sum(s.Rational(c.numerator,c.denominator)*s.prod(variables[i] for i in mon)
               for mon,c in row.items())


def serialize(matrix):
    return [[str(x) for x in row] for row in matrix.tolist()]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    source=ROOT/'reports/TORIC_ALGEBRA_FIBER_20261031.json'
    data=json.loads(source.read_text())
    variables=s.symbols('c0:20')
    original=[s.sympify(eq) for eq in data['quadratic_equations']]
    polynomials=[sparse(eq,variables) for eq in original]
    basis=tracked_rows(polynomials)
    linear=[(row,combination) for pivot,(row,combination) in basis.items() if len(pivot)<=1]
    assert len(basis)==63 and len(linear)==18
    # Verify the combination certificate and independence from beta in degree 2.
    zero_m=multiplication_from_jet(s.zeros(24,1),dict(zip(TRIPLES,variables,strict=True)))
    _,zero_equations=associativity_quadrics(zero_m)
    zero_polys=[sparse(eq,variables) for eq in zero_equations]
    for sample,zero in zip(polynomials,zero_polys,strict=True):
        assert {mon:v for mon,v in sample.items() if len(mon)==2}==zero
    for row,combination in linear:
        recovered={}; zero_recovered={}
        for index,multiplier in combination.items():
            add_multiple(recovered,polynomials[index],multiplier)
            add_multiple(zero_recovered,zero_polys[index],multiplier)
        assert recovered==row and not zero_recovered
    equations=[expression(row,variables) for row,_ in linear]
    coefficients=s.Matrix([[row.get((i,),Fraction(0)) for i in range(20)] for row,_ in linear])
    constants=s.Matrix([row.get((),Fraction(0)) for row,_ in linear])
    rref,pivots=coefficients.row_join(constants).rref()
    assert len(pivots)==18 and all(i<20 for i in pivots)
    free_indices=[i for i in range(20) if i not in pivots]
    assert free_indices==[18,19]
    free=[variables[i] for i in free_indices]; x,y=free
    solution=list(variables)
    for row,pivot in enumerate(pivots):
        solution[pivot]=-rref[row,20]-sum(rref[row,j]*variables[j] for j in free_indices)
    substitution=dict(zip(variables,solution,strict=True))
    assert all(s.expand(eq.subs(substitution,simultaneous=True))==0 for eq in equations)
    linear_minor=coefficients[:,list(pivots)]
    linear_determinant=linear_minor.to_DM().det()
    assert linear_determinant
    residual=[s.expand(eq.subs(substitution,simultaneous=True)) for eq in original]
    monomials=[x*x,x*y,y*y]
    leading=s.Matrix([[s.Poly(eq,x,y).coeff_monomial(mon) for mon in monomials] for eq in residual])
    rows=list(leading.T.rref()[1]); assert len(rows)==3
    leading_minor=leading[rows,:]
    leading_determinant=leading_minor.det();assert leading_determinant
    selected=s.Matrix([residual[i] for i in rows])
    lower=selected-leading_minor*s.Matrix(monomials)
    normal_forms=(-leading_minor.inv()*lower).applyfunc(s.expand)
    assert all(s.Poly(eq,x,y).total_degree()<=1 for eq in normal_forms)
    normalized_quadrics=[s.expand(mon-nf) for mon,nf in zip(monomials,normal_forms,strict=True)]
    # Every residual equation is already in the span of these three quadrics.
    for row,eq in enumerate(residual):
        assert s.expand(eq-sum(leading[row,j]*normalized_quadrics[j] for j in range(3)))==0
    def normal_vector(poly):
        p=s.Poly(poly,x,y)
        return s.Matrix([p.coeff_monomial(1),p.coeff_monomial(x),p.coeff_monomial(y)])
    nx, nxy, ny=[normal_vector(eq) for eq in normal_forms]
    mx=s.Matrix.hstack(s.Matrix([0,1,0]),nx,nxy)
    my=s.Matrix.hstack(s.Matrix([0,0,1]),nxy,ny)
    assert mx*my==my*mx
    identity=s.eye(3)
    assert all((mx*mx-nx[0]*identity-nx[1]*mx-nx[2]*my)==s.zeros(3)
               for _ in [0])
    assert mx*my-nxy[0]*identity-nxy[1]*mx-nxy[2]*my==s.zeros(3)
    assert my*my-ny[0]*identity-ny[1]*mx-ny[2]*my==s.zeros(3)
    # Cyclic y explicitly identifies the entire algebra with a separable cubic.
    cubic=y**3-s.Rational(24010,13)*y
    relation=x-s.Rational(13,3430)*y*y+y/7+s.Rational(14,3)
    assert my**3-s.Rational(24010,13)*my==s.zeros(3)
    assert mx-s.Rational(13,3430)*my**2+my/7+s.Rational(14,3)*identity==s.zeros(3)
    unit=s.Matrix([1,0,0])
    cyclic=s.Matrix.hstack(unit,my*unit,my**2*unit)
    assert cyclic.det()!=0
    discriminant=s.discriminant(cubic,y);assert discriminant!=0
    # Independent polynomial ideal equality with the displayed lex generators.
    groebner=s.groebner([relation,cubic],x,y,order='lex')
    assert all(groebner.reduce(eq)[1]==0 for eq in normalized_quadrics)
    check=s.groebner(normalized_quadrics,x,y,order='lex')
    assert check.reduce(relation)[1]==check.reduce(cubic)[1]==0
    # Compute det(q) inside the three-dimensional algebra, without huge expansion.
    q=s.Matrix(data['trace_form_matrix']).subs(substitution,simultaneous=True)
    q=q.subs(x,s.Rational(13,3430)*y*y-y/7-s.Rational(14,3))
    q=q.applyfunc(lambda eq:s.rem(s.expand(eq),cubic,y))
    determinant_remainder=s.rem(s.expand(q.det()),cubic,y)
    bezout_q,bezout_cubic,gcd=s.gcdex(determinant_remainder,cubic,y)
    assert gcd==1
    assert s.expand(bezout_q*determinant_remainder+bezout_cubic*cubic)==1
    report={
        'id':'TORIC_ALGEBRA_DEGREE_20261031','all_passed':True,'field':'QQ',
        'original_equations':96,'quadratic_span_rank':45,'total_span_rank':63,
        'universal_linear_consequences':[
            {'equation_at_sample':str(expression(row,variables)),
             'constant_combination_of_96_equations':[[index,str(value)] for index,value in sorted(combination.items())]}
            for row,combination in linear],
        'quadratic_coefficients_equal_beta_zero_model':True,
        'linear_pivot_columns_zero_based':list(pivots),
        'linear_pivot_matrix_at_sample':serialize(linear_minor),
        'linear_pivot_determinant_at_sample':str(linear_determinant),
        'free_variables':[str(v) for v in free],
        'eliminated_c0_through_c19':[str(v) for v in solution],
        'residual_quadrics':[str(v) for v in residual],
        'residual_leading_monomials':[str(v) for v in monomials],
        'residual_leading_pivot_rows_zero_based':rows,
        'residual_leading_minor':serialize(leading_minor),
        'residual_leading_determinant':str(leading_determinant),
        'quadratic_normal_forms':[str(v) for v in normal_forms],
        'quotient_basis':['1',str(x),str(y)],
        'multiplication_by_x':serialize(mx),'multiplication_by_y':serialize(my),
        'lex_groebner_basis':[str(relation),str(cubic)],
        'sample_algebra_dimension':3,'sample_fiber_reduced':True,
        'cubic_discriminant':str(discriminant),
        'trace_determinant_mod_cubic':str(determinant_remainder),
        'trace_determinant_bezout':[str(bezout_q),str(bezout_cubic)],
        'all_three_sample_points_have_invertible_trace_form':True,
        'generic_toric_first_jet_degree':3,
        'generic_ordered_coframe_degree':360,
        'generic_degree_justification':'Universal linear consequences and the invertible leading-conic matrix give a finite model near this jet; after removing finite images of discriminant-zero and nonetale loci, the semisimple model GL4/S5 is finite etale with three-point fiber. See proof note.',
        'limitations':['The explicit open base equations and branch divisor are not expanded.',
                       'No degree for T39/T40 and no full S61 exhaustion.',
                       'Boundary and exceptional fibers remain unclassified.'],
        'proof_note':'research/notes/TORIC_GENERIC_DEGREE_20261031.md',
        'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in (Path(__file__),ROOT/'code/toric_algebra_fiber.py',source)},
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report,out,ensure_ascii=False,indent=2);out.write('\n')
    print(json.dumps({'all_passed':True,'linear_rank':18,'sample_degree':3,
                      'generic_toric_degree':3,'coframe_degree':360,
                      'leading_determinant':str(leading_determinant),'output':str(args.output)}))


if __name__=='__main__':
    main()

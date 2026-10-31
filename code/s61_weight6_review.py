#!/usr/bin/env python3
"""Independent certificate checker for the universal weight-six consequence.

No generator or model routines are imported. Exact rational kernel identities
give an upper rank bound; a separate modular calculation gives the matching
lower bound, thereby certifying the rank over QQ without lifting guesses.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def add(target,key,value):
    updated = target.get(key,F(0))+value
    if updated:
        target[key]=updated
    else:
        target.pop(key,None)


def add_polynomial(target,source,scalar=F(1),variable=None):
    for mon,value in source.items():
        key = tuple(sorted((*mon,variable))) if variable is not None else mon
        add(target,key,scalar*value)


def decode(encoded):
    out = {}
    for mon,numerator,denominator in encoded:
        add(out,tuple(sorted(mon)),F(numerator,denominator))
    return out


def relation(encoded):
    out = {index:F(value) for index,value in encoded}
    assert len(out)==len(encoded)
    return out


def modular_rank(rows,prime=32003):
    basis = {}
    for original in rows:
        row = {}
        for key,value in original.items():
            value = F(value)
            assert value.denominator%prime
            coefficient = value.numerator*pow(value.denominator,-1,prime)%prime
            if coefficient:
                row[key]=coefficient
        while row:
            pivot = min(row)
            if pivot not in basis:
                inverse = pow(row[pivot],-1,prime)
                basis[pivot]={key:value*inverse%prime for key,value in row.items()}
                break
            scalar = row[pivot]
            for key,value in basis[pivot].items():
                updated = (row.get(key,0)-scalar*value)%prime
                if updated:
                    row[key]=updated
                else:
                    row.pop(key,None)
    return len(basis)


def leading_identity(coefficients,leading):
    out = {}
    for index,scalar in coefficients.items():
        variable,equation = divmod(index,70)
        assert 0<=variable<45 and 0<=equation<70
        add_polynomial(out,leading[equation],scalar,variable)
    return out


def centered_consequence(coefficients,equations,center):
    out = {}
    for index,scalar in coefficients.items():
        variable,equation = divmod(index,70)
        add_polynomial(out,equations[equation],scalar,variable)
        add_polynomial(out,equations[equation],scalar*center[variable])
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    report_path = ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    source_path = ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json'
    point_path = ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json'
    flag_path = ROOT/'reports/S61_INFINITY_FLAGS_20261031_final.json'
    certificate = json.loads(report_path.read_text())
    source = json.loads(source_path.read_text())
    weights = source['parameter_weights']
    columns = [i for i,w in enumerate(weights) if w==2]
    index = {old:new for new,old in enumerate(columns)}
    assert columns==certificate['q_to_u_indices_zero_based']
    equations = [decode(eq) for eq in source['equations']]
    assert [{sum(weights[i] for i in mon) for mon in eq} for eq in equations]==[{4}]*70+[{5}]*64
    leading = []
    for equation in equations[:70]:
        row = {}
        for mon,value in equation.items():
            if all(i in index for i in mon):
                assert len(mon)==2
                add(row,tuple(index[i] for i in mon),value)
        leading.append(row)
    kernels = [relation(row) for row in certificate['all_36_kernel_vectors']]
    bianchi = [relation(row) for row in certificate['canonical_35_Bianchi_vectors']]
    selected = relation(certificate['selected_non_Bianchi_syzygy'])
    assert len(kernels)==36 and len(bianchi)==35 and len(selected)==125
    assert all(not leading_identity(row,leading) for row in kernels+bianchi+[selected])
    assert modular_rank(kernels)==36
    assert modular_rank(bianchi)==35 and modular_rank(bianchi+[selected])==36
    cubic_rows = ({tuple(sorted((*mon,j))):value for mon,value in equation.items()}
                  for j in range(45) for equation in leading)
    assert modular_rank(cubic_rows)==3114
    # 36 independent QQ kernel vectors prove rank <= 3114; reduction mod32003
    # proves rank >= 3114. These two exact inequalities prove rank_QQ=3114.
    print('all rational syzygies checked; exact QQ rank3114 follows from matching bounds',flush=True)

    universal = {}
    for key,value in selected.items():
        j,i = divmod(key,70)
        add_polynomial(universal,equations[i],value,columns[j])
    assert universal==decode(certificate['universal_new_equation'])
    assert len(universal)==9257
    assert {sum(weights[i] for i in mon) for mon in universal}=={6}
    assert {sum(weights[i]==2 for i in mon) for mon in universal}=={1,2}
    flag = {}
    for mon,value in universal.items():
        if [i for i in mon if weights[i]==2]==[37,37]:
            add(flag,tuple(i for i in mon if weights[i]==1),value)
    assert not flag
    # A second flag distinguishes incidence vanishing from vanishing for all
    # beta: here the restriction is a nonzero multiple of the incidence form.
    extra_flag = json.loads(flag_path.read_text())['additional_flag_check_all_first_jets']
    flag_direction = {int(name[1:])-1:F(value) for name,value in extra_flag['quadratic_direction_S61'].items()}
    one = [i for i,w in enumerate(weights) if w==1]
    beta_variables = s.symbols('b0:16')
    beta_index = dict(zip(one,beta_variables,strict=True))
    restriction = s.S.Zero
    for mon,value in universal.items():
        qpart = [i for i in mon if weights[i]==2]
        if len(qpart)==2:
            restriction += s.Rational(value)*s.prod(flag_direction.get(i,0) for i in qpart)*s.prod(beta_index[i] for i in mon if weights[i]==1)
    restriction = s.expand(restriction)
    incidence = sum(s.Rational(value)*b for value,b in zip(extra_flag['incidence_coefficients_in_old_weight_one_order'],beta_variables,strict=True))
    quotient,remainder = s.div(restriction,incidence,*beta_variables)
    assert remainder==0 and quotient!=0

    points = json.loads(point_path.read_text())['points']
    point = [F(v) for v in next(item for item in points if item['component_id']=='T01')['S61_coordinates_u1_through_u61']]
    center = [point[i] for i in columns]
    def specialize_and_center(polynomial):
        out = {}
        for mon,value in polynomial.items():
            free = []
            for i in mon:
                if i in index:
                    free.append(index[i])
                else:
                    value *= point[i]
            expanded = {():value} if value else {}
            for j in free:
                nxt = {}
                add_polynomial(nxt,expanded,F(1),j)
                add_polynomial(nxt,expanded,center[j])
                expanded = nxt
            add_polynomial(out,expanded)
        return out
    full = [specialize_and_center(eq) for eq in equations]
    assert all(() not in eq for eq in full)
    new = specialize_and_center(universal)
    saved_new = {tuple(mon):F(value) for mon,value in certificate['sample_centered_new_equation']}
    assert new==saved_new==centered_consequence(selected,full,center)
    assert () not in new and all(len(mon)<=2 for mon in new)
    assert len(new)==790
    for row,saved in zip(bianchi,certificate['sample_Bianchi_constant_redundancies'],strict=True):
        actual = centered_consequence(row,full,center)
        predicted = {}
        for equation,value in relation(saved).items():
            assert 0<=equation<134
            add_polynomial(predicted,full[equation],value)
        assert actual==predicted
    monomials = [tuple(mon) for mon in certificate['sample_augmented_minor_monomials']]
    assert len(monomials)==len(set(monomials))==135 and all(len(mon)==2 for mon in monomials)
    matrix = s.Matrix([[s.Rational(row.get(mon,0)) for mon in monomials] for row in full+[new]])
    determinant = matrix.to_DM().det()
    assert determinant and str(determinant)==certificate['sample_augmented_minor_determinant']
    print('universal9257-term G6, all35 redundancies and nonzero135minor independently verified',flush=True)
    report = {
        'id':'S61_UNIVERSAL_WEIGHT6_REVIEW_20261031','all_passed':True,
        'field':'QQ, with exact mod32003 rank lower bounds paired with rational kernel upper bounds',
        'linear_syzygy_matrix_shape':[3150,16215],'exact_QQ_rank':3114,
        'verified_rational_kernel_vectors':36,'canonical_Bianchi_rank':35,
        'selected_non_Bianchi_multiplier_terms':125,'universal_consequence_terms':9257,
        'universal_weight':6,'degrees_in_q':[1,2],
        'sample_new_equation_terms':790,'sample_Bianchi_redundancies_checked':35,
        'sample_augmented_quadratic_rank':135,
        'sample_augmented_minor_determinant':str(determinant),
        'augmented_minor_degree_bound_in_beta':66,
        'specified_u38_flag_leading_restriction_is_zero':True,
        'additional_flag_incidence_divisibility_checked':True,
        'additional_flag_quotient_in_old_weight_one_order':str(s.factor(quotient)),
        'not_proved':['No full affine-fiber dimension or degree bound.',
                      'No new Poisson condition: G6 belongs to the existing ideal.',
                      'The all-flags extension is not part of this certificate.'],
        'proof_note':'research/notes/S61_UNIVERSAL_WEIGHT6_REVIEW_20261031.md',
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                        for path in (report_path,source_path,point_path,flag_path,Path(__file__))},
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report,stream,ensure_ascii=False,indent=2);stream.write('\n')
    print(json.dumps({'all_passed':True,'kernel_dimension':36,'augmented_rank':135,'output':str(args.output)}))


if __name__=='__main__':
    main()

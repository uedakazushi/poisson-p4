#!/usr/bin/env python3
"""Independent review of the unsaturated infinity-flag obstruction.

Uses the immutable sparse data directly, without importing s61_model or the
flag generator. Reconstructs directional bivectors and checks their Schouten
brackets, both flag restrictions, all Bianchi identities, and primitive trace.
"""
from __future__ import annotations

import argparse
import hashlib
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'legacy/v02/data'
J = s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
PAIRS = list(combinations(range(4),2))


def polynomial(encoded):
    out = {}
    for mon,numerator,denominator in encoded:
        key = tuple(sorted(mon))
        out[key] = out.get(key,0)+s.Rational(numerator,denominator)
    return {key:value for key,value in out.items() if value}


def subtract_multiple(target,source,coefficient,variable=None):
    for mon,value in source.items():
        key = tuple(sorted(mon+(variable,))) if variable is not None else mon
        target[key] = target.get(key,0)-coefficient*value
        if not target[key]:
            del target[key]


def affine_from_raw(raw,z):
    xx = (s.S.One,*z)
    monomials = [xx[i]*xx[j] for i,j in combinations_with_replacement(range(5),2)]
    ambient = s.zeros(5)
    for block,(i,j) in enumerate(combinations(range(5),2)):
        ambient[i,j] = s.expand(sum(raw[15*block+k]*mon for k,mon in enumerate(monomials)))
        ambient[j,i] = -ambient[i,j]
    return s.Matrix(4,4,lambda i,j:s.expand(ambient[i+1,j+1]-z[i]*ambient[0,j+1]+z[j]*ambient[0,i+1]))


def jacobi(matrix,z):
    return [s.expand(sum(matrix[i,h]*s.diff(matrix[j,k],z[h])
                         +matrix[j,h]*s.diff(matrix[k,i],z[h])
                         +matrix[k,h]*s.diff(matrix[i,j],z[h]) for h in range(4)))
            for i,j,k in combinations(range(4),3)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    paths = {name:SOURCE/name for name in ('symplectic_slice_certificate.json',
        'symplectic_slice_equations.json','symplectic_slice_reduced_equations.json',
        'symplectic_bianchi_identities.json')}
    cert,raw,reduced,bianchi = [json.loads(path.read_text()) for path in paths.values()]
    reviewed_path = ROOT/'reports/S61_INFINITY_FLAGS_20261031_final.json'
    reviewed = json.loads(reviewed_path.read_text())
    weights = reduced['parameter_weights']
    one = [i for i,w in enumerate(weights) if w==1]
    two = [i for i,w in enumerate(weights) if w==2]
    one_index = {old:k for k,old in enumerate(one)}
    u = s.symbols('b0:16')
    assert len(one)==16 and len(two)==45
    assert len(raw['eliminated_raw_variables'])==55
    assert max(sum(weights[i]==2 for i in mon) for item in raw['eliminated_raw_variables']
               for mon,_,_ in item['polynomial'])==1
    assert all(sum(weights[i]==2 for i in mon)<=2 for equation in reduced['equations'] for mon,_,_ in equation)
    old = [polynomial(eq) for eq in raw['equations']]
    assert [raw['equations'][i] for i in reduced['old_generator_indices']]==reduced['equations']
    for identity in bianchi['generator_redundancy_identities']:
        difference = dict(old[identity['eliminated_old_generator']])
        for index,numerator,denominator in identity['constant_coefficients_on_retained_old_generators']:
            subtract_multiple(difference,old[index],s.Rational(numerator,denominator))
        for index,parameter,numerator,denominator in identity['linear_coefficients_on_weight4_old_generators']:
            assert weights[parameter]==1
            subtract_multiple(difference,old[index],s.Rational(numerator,denominator),parameter)
        assert not difference
    assert len(bianchi['generator_redundancy_identities'])==20

    z = s.Matrix(s.symbols('z1:5'))
    flag_reports = []
    for key in ('representative','additional_flag_check_all_first_jets'):
        flag = reviewed[key]
        v,w = (s.Matrix([s.Rational(a) for a in flag[name]]) for name in ('v','w'))
        assert (v.T*J*w)[0]==0 and v.row_join(w).rank()==2
        direction = [s.S.Zero]*61
        for name,value in flag['quadratic_direction_S61'].items():
            direction[int(name[1:])-1] = s.Rational(value)
        assert all(direction[i]==0 for i in one)
        directional_raw = [s.S.Zero]*150
        for index,value in zip(cert['raw_free_indices'],direction,strict=True):
            directional_raw[index] = value
        for item in raw['eliminated_raw_variables']:
            value = 0
            for mon,numerator,denominator in item['polynomial']:
                qindices = [i for i in mon if i in two]
                if len(qindices)!=1:
                    continue
                value += s.Rational(numerator,denominator)*direction[qindices[0]]*s.prod(u[one_index[i]] for i in mon if i in one_index)
            directional_raw[item['raw_index']] = s.expand(value)
        b = affine_from_raw(directional_raw,z)
        q = b.subs(dict.fromkeys(u,0))
        ell = -J*v
        assert q == (ell.dot(z)**2*(v*w.T-w*v.T)).applyfunc(s.expand)
        assert jacobi(q,z)==[0]*4
        assert jacobi(J+q,z)==[0]*4
        # Directly check the quadratic-in-direction Jacobiator [B,B]/2.
        direct_rows = []
        for expression in jacobi(b,z):
            for value in s.Poly(expression,*z).coeffs():
                p = s.Poly(value,*u)
                assert p.total_degree()<=1 and p.coeff_monomial(1)==0
                direct_rows.append([p.coeff_monomial(ui) for ui in u])
        direct = s.Matrix(direct_rows)
        incidence = s.Matrix([[s.Rational(value) for value in flag['incidence_coefficients_in_old_weight_one_order']]])
        assert direct.rank()==direct.col_join(incidence).rank()==1
        # Restrict the 134 leading quadrics without using the generator's helper.
        leading = []
        for equation in reduced['equations']:
            row = [s.S.Zero]*17
            for mon,numerator,denominator in equation:
                qindices = [i for i in mon if i in two]
                if len(qindices)!=2:
                    continue
                linear = [i for i in mon if i in one_index]
                assert len(linear)<=1
                column = 1+one_index[linear[0]] if linear else 0
                row[column] += s.Rational(numerator,denominator)*s.prod(direction[i] for i in qindices)
            leading.append(row)
        leading = s.Matrix(leading)
        assert leading == s.Matrix(flag['leading_equation_matrix_constant_then_16_weight_one'])
        assert leading[:,0]==s.zeros(134,1)
        assert leading[:,1:].rank()==leading[:,1:].col_join(incidence).rank()==1
        flag_reports.append({'flag':key,'direct_directional_Jacobi_rank':1,
                             'leading_134_restriction_rank':1,'incidence_row_agrees':True})
        print('direct Schouten and all134 verified:',key,flush=True)

    # Check the displayed primitive trace from its differential-form formula.
    fcoeff = s.Matrix(s.symbols('f1:5'))
    f = fcoeff.dot(z)
    contraction = -J*z
    beta = (fcoeff*contraction.T-contraction*fcoeff.T+2*f*J)/3
    assert s.expand(beta[0,1]+beta[2,3]-5*f/3)==0
    a = s.Matrix(s.symbols('a1:5'))
    vv,ww = s.Matrix(s.symbols('v1:5')),s.Matrix(s.symbols('w1:5'))
    tangent = -2*a.dot(vv)*J+vv*(J*a).T-(J*a)*vv.T
    assert s.expand((vv.T*(-J*tangent*J)*ww)[0]+3*a.dot(vv)*(vv.T*J*ww)[0])==0
    injection = reviewed['incidence_injectivity']
    basis = s.Matrix(injection['primitive_closed_basis'])
    assert basis.shape==(24,16) and basis.rank()==16
    forms = []
    for column in range(16):
        matrices = [s.zeros(4) for _ in range(4)]
        for block,(i,j) in enumerate(PAIRS):
            for k in range(4):
                matrices[k][i,j] = basis[4*block+k,column]
                matrices[k][j,i] = -basis[4*block+k,column]
        assert all(m[0,1]+m[2,3]==0 for m in matrices)
        assert all(matrices[i][j,k]-matrices[j][i,k]+matrices[k][i,j]==0 for i,j,k in combinations(range(4),3))
        beta_v = sum((vv[k]*m for k,m in enumerate(matrices)),s.zeros(4))
        forms.append(s.Poly(s.expand((vv.T*beta_v*ww)[0]),*vv,*ww))
    omega = (vv.T*J*ww)[0]
    forms += [s.Poly(s.expand(omega*vi),*vv,*ww) for vi in vv]
    matrix = s.Matrix([[p.coeff_monomial(tuple(mon)) for p in forms] for mon in injection['coefficient_monomials']])
    assert matrix==s.Matrix(injection['coefficient_matrix'])
    witness = injection['rank20_witness']
    assert matrix[witness['rows_zero_based'],witness['columns_zero_based']].det()==625
    report = {
        'id':'S61_INFINITY_FLAGS_REVIEW_20261031','all_passed':True,'field':'QQ',
        'affine_direction_reconstruction_entries':55,'exact_Bianchi_identities_checked':20,
        'flag_checks':flag_reports,'incidence_injectivity_minor':'625',
        'primitive_trace':'Lambda((df wedge i_E Omega + 2 f Omega)/3) = 5 f/3',
        'perspective_incidence_change':'-3 a(v) Omega(v,w)',
        'accepted_claim':'The UNSATURATED generator-wise quadratic homogenization contains a closed flag incidence of pure dimension3 for primitive beta nonzero, or dimension4 for beta zero.',
        'not_proved':['Membership in the true projective closure after t-saturation.',
                      'Positive-dimensional affine fibers or any full generic affine-fiber bound.',
                      'Reducedness of the flag incidence divisor or equality with the entire infinity scheme.'],
        'proof_review_note':'research/notes/S61_INFINITY_FLAGS_REVIEW_20261031.md',
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                        for path in (*paths.values(),reviewed_path,Path(__file__))},
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report,stream,ensure_ascii=False,indent=2);stream.write('\n')
    print(json.dumps({'all_passed':True,'direct_flags':2,'Bianchi_identities':20,'output':str(args.output)}))


if __name__=='__main__':
    main()

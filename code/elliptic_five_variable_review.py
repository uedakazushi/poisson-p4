#!/usr/bin/env python3
"""Independent tensor, normalization, and finite-subgroup audit.

No functions from elliptic_six_variable_fiber are imported. Universal Lie
covariance is checked with 24 arbitrary alternating-tensor coefficients.
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


def tensors(vector):
    out = [s.zeros(4) for _ in range(4)]
    for block,(i,j) in enumerate(PAIRS):
        for k in range(4):
            out[k][i,j] = vector[4*block+k]
            out[k][j,i] = -vector[4*block+k]
    return out


def covariants(t):
    # Matrix traces give a separate implementation of the four-index sum.
    m = s.Matrix(4,4,lambda i,a:s.trace(t[i]*J*t[a]*J))
    a = J*m
    h = J*m*J.T
    v = -J*sum((t[k]*h[:,k] for k in range(4)),s.zeros(4,1))
    return m,a,v,s.Matrix.hstack(v,a*v,a*a*v,a*a*a*v)


def signature(t):
    frame = covariants(t)[-1]
    gram = frame.T*J*frame
    transformed = [sum((frame[a,k]*frame.T*t[a]*frame for a in range(4)),s.zeros(4)) for k in range(4)]
    return s.Matrix([gram[i,j] for i,j in PAIRS]+[transformed[k][i,j] for i,j in PAIRS for k in range(4)])


def scalar_invariants(t):
    _,a,v,_ = covariants(t)
    moved = a*v
    at_moved = sum((moved[k]*t[k] for k in range(4)),s.zeros(4))
    i4 = s.trace(a*a)
    i14 = (v.T*at_moved*moved)[0]
    return i4,i14,i14/i4**3


def normalized_tensor(homogeneous,coordinates,lam,value,point,scalar):
    z = s.Matrix(s.symbols('z1:5'))
    subst = dict(zip(coordinates,(1,*z),strict=True))
    affine = s.Matrix(4,4,lambda i,j:s.expand((homogeneous[i+1,j+1]
        -coordinates[i+1]*homogeneous[0,j+1]+coordinates[j+1]*homogeneous[0,i+1]).subs(subst)))
    point_sub = dict(zip(z,point,strict=True)) | {lam:s.Rational(value)}
    k = scalar*affine.subs(point_sub)
    derivative = [scalar*affine.diff(zi).subs(point_sub) for zi in z]
    # Symplectic Gram--Schmidt, rather than the saved rational frame formula.
    e = s.eye(4)
    f0 = e[0,:]
    f1 = e[3,:]/k[0,3]
    project = lambda row: row-(row*k*f1.T)[0]*f0+(row*k*f0.T)[0]*f1
    f2 = project(e[1,:])
    f3 = project(e[2,:])
    f3 /= (f2*k*f3.T)[0]
    frame = s.Matrix.vstack(f0,f1,f2,f3)
    assert frame*k*frame.T==J
    inverse = frame.inv()
    linear = [frame*sum((inverse[h,j]*derivative[h] for h in range(4)),s.zeros(4))*frame.T for j in range(4)]
    lowered = [-J*l*J for l in linear]
    tau = sum(z[k]*(lowered[k][0,1]+lowered[k][2,3]) for k in range(4))
    dtau = s.Matrix([s.diff(tau,zi) for zi in z])
    contraction = -J*z
    removed = (2*tau*J+dtau*contraction.T-contraction*dtau.T)/5
    beta = sum((z[k]*lowered[k] for k in range(4)),s.zeros(4))-removed
    result = [beta.diff(zi) for zi in z]
    assert all(t[0,1]+t[2,3]==0 for t in result)
    assert all(result[i][j,k]-result[j][i,k]+result[k][i,j]==0 for i,j,k in combinations(range(4),3))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    tensor_variables = s.symbols('t0:24')
    t = tensors(tensor_variables)
    m,a,v,_ = covariants(t)
    m = m.applyfunc(s.expand)
    v = v.applyfunc(s.expand)
    h = J*m*J.T
    assert m==m.T
    assert a*J+J*a.T==s.zeros(4)
    for i,j in combinations_with_replacement(range(4),2):
        symmetric = s.zeros(4)
        symmetric[i,j]=symmetric[j,i]=1
        x = J*symmetric
        delta_t = [-sum((x[b,k]*t[b] for b in range(4)),s.zeros(4))-x.T*t[k]-t[k]*x for k in range(4)]
        delta_m = s.Matrix(4,4,lambda k,l:s.trace(delta_t[k]*J*t[l]*J+t[k]*J*delta_t[l]*J))
        assert (delta_m+x.T*m+m*x).applyfunc(s.expand)==s.zeros(4)
        delta_h = J*delta_m*J.T
        delta_v = -J*sum((delta_t[k]*h[:,k]+t[k]*delta_h[:,k] for k in range(4)),s.zeros(4,1))
        assert (delta_v-x*v).applyfunc(s.expand)==s.zeros(4,1)
    print('universal covariance checked on all10 sp4 generators',flush=True)

    report_path = ROOT/'reports/ELLIPTIC_SIX_VARIABLE_FIBER_20261031.json'
    source_path = ROOT/'reports/ELLIPTIC_WEIGHT1_DOMINANCE_20261031.json'
    reviewed = json.loads(report_path.read_text())
    source = json.loads(source_path.read_text())
    family_results = []
    xcoords = s.symbols('x0:5')
    lam = s.Symbol('lam')
    zeta = s.Symbol('zeta')
    cyclotomic = sum(zeta**i for i in range(5))
    shift,reflection = s.zeros(5),s.zeros(5)
    for i in range(5):
        shift[i,(i+1)%5]=1
        reflection[i,(-i)%5]=1
    phase = s.diag(*(zeta**i for i in range(5)))
    phase_inverse = s.diag(*(zeta**((-i)%5) for i in range(5)))
    for index,saved in enumerate(reviewed['families']):
        homogeneous = s.zeros(5)
        for (i,j),value in zip(combinations(range(5),2),source['families'][index]['homogeneous_bivector_upper_triangle'],strict=True):
            homogeneous[i,j]=s.sympify(value)
            homogeneous[j,i]=-homogeneous[i,j]
        tensor = normalized_tensor(homogeneous,xcoords,lam,2,[1,2,3,4],1)
        vector = s.Matrix([tensor[k][i,j] for i,j in PAIRS for k in range(4)])
        assert vector==s.Matrix(saved['primitive_jet'])
        mm,aa,vv,rr = covariants(tensor)
        assert mm==s.Matrix(saved['M']) and aa==s.Matrix(saved['A'])
        assert vv==s.Matrix(saved['v']) and rr==s.Matrix(saved['R'])
        determinant = rr.det(method='domain-ge')
        assert determinant==s.Rational(saved['det_R']) and determinant!=0
        i4,i14,h = scalar_invariants(tensor)
        assert i4==s.Rational(saved['I4'])!=0
        assert i14==s.Rational(saved['I14'])!=0
        for name,g,inverse,character in (('shift',shift,shift.T,1),('phase',phase,phase_inverse,1),('reflection',reflection,reflection,-1)):
            pushed = g*homogeneous.subs(dict(zip(xcoords,inverse*s.Matrix(xcoords),strict=True)),simultaneous=True)*g.T
            assert all(s.expand(s.rem(s.expand(value),cyclotomic,zeta))==0 for value in pushed-character*homogeneous)
        sig = signature(tensor)
        column_degrees = [3,5,7,9]
        signature_degrees = ([column_degrees[i]+column_degrees[j] for i,j in PAIRS]
            +[1+column_degrees[k]+column_degrees[i]+column_degrees[j] for i,j in PAIRS for k in range(4)])
        assert all(degree%2==0 for degree in signature_degrees)
        for scalar in (s.Integer(-2),s.Integer(3)):
            scaled_tensor = normalized_tensor(homogeneous,xcoords,lam,2,[1,2,3,4],scalar)
            scaled_i4,scaled_i14,scaled_h = scalar_invariants(scaled_tensor)
            assert scaled_i4==scalar**2*i4 and scaled_i14==scalar**7*i14
            assert scaled_h==scalar*h
            assert covariants(scaled_tensor)[-1].det(method='domain-ge')==scalar**12*determinant
            recovered_scalar = scaled_h/h
            assert recovered_scalar==scalar
            weighted_signature = s.Matrix([recovered_scalar**(degree//2)*entry
                for degree,entry in zip(signature_degrees,sig,strict=True)])
            assert signature(scaled_tensor)==weighted_signature
        # Actual chart parameters and scalar signs, not just ambient identities.
        for g,character in ((shift,1),(reflection,-1)):
            moved = g*s.Matrix([1,1,2,3,4])
            moved_point = [moved[i]/moved[0] for i in range(1,5)]
            moved_tensor = normalized_tensor(homogeneous,xcoords,lam,2,moved_point,character)
            assert signature(moved_tensor)==sig
            moved_frame = covariants(moved_tensor)[-1]
            recovered = moved_frame*rr.inv()
            assert recovered*J*recovered.T==J
            inverse = recovered.inv()
            expected = [sum((inverse[b,k]*inverse.T*tensor[b]*inverse for b in range(4)),s.zeros(4)) for k in range(4)]
            assert expected==moved_tensor
        family_results.append({'component_id':saved['component_id'],'det_R':str(determinant),
            'I4':str(i4),'I14':str(i14),
            'independent_Gram_Schmidt_normalization_agrees':True,
            'generic_lambda_shift_phase_reflection_identities':True,
            'actual_shift_reflection_parameter_signature_checks':True,
            'scalar_recovery_and_all30_weighted_signatures_checked_at':['-2','3']})
        print('independent normalization and50 subgroup checked:',saved['component_id'],flush=True)
    report = {
        'id':'ELLIPTIC_FIVE_VARIABLE_REVIEW_20261031','all_passed':True,'field':'QQ',
        'universal_tensor_variables':24,'infinitesimal_covariance_generators':10,
        'covariance':'M -> S^-T M S^-1, A -> S A S^-1, v -> S v, R -> S R',
        'families':family_results,'scheme_equivalence_on_localized_domain':True,
        'parameter_variables':5,'signature_equations':30,
        'signature_homogeneous_degrees':signature_degrees,
        'open_conditions':'lambda*a*Delta*det(R_source_at_c1)*det(R_target)*I4(source_at_c1)*I14(source_at_c1)*I4(target)*I14(target) != 0',
        'recovered_symplectic_element':'R_target R_source^-1',
        'recovered_scalar':'h(target)/h(source_at_c1), h=I14/I4^3',
        'verified_redundancy_subgroup_order':50,
        'limitations':['The det(R)=0 locus and normalization-chart complements are excluded.',
                      'Thirty signature equalities are a redundant system, not a claim of thirty independent equations.',
                      'The order50 subgroup is not a determination of all parameter redundancies.',
                      'Neither elliptic component degree nor a full S61 generic-fiber count is computed.'],
        'proof_note':'research/notes/ELLIPTIC_FIVE_VARIABLE_REVIEW_20261031.md',
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in (source_path,report_path,ROOT/'code/elliptic_six_variable_fiber.py',Path(__file__))},
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report,stream,ensure_ascii=False,indent=2);stream.write('\n')
    print(json.dumps({'all_passed':True,'families':2,'universal_sp4_generators':10,'output':str(args.output)}))


if __name__=='__main__':
    main()

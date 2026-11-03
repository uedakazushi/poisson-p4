#!/usr/bin/env python3
"""Global necessary T01 constraints onQa, with correct closed-jet coordinates.

These vanishing tests are NOT the membership proof.  That proof is the separate
commuting-matrix construction, valid even at the singular canonical beta.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import sympy as s

from s61_model import S61,X,pushforward
from q51_universal_section_20261101 import affine
from toric_closed_jet_coordinates_20261102 import coords,pairs,m2
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from elliptic_curve_invariants_20261101 import curve_invariants_from_beta
from toric_coframe_fiber import affine_bivector,Y
import toric_zero_fiber_coframe_quadrics_20261102 as lifts

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    a=s.Symbol('a');model=S61();local=ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63/certificate.json'
    parameters=list(map(s.sympify,json.loads(local.read_text())['scalar_section_normal61']))
    Q=s.zeros(5)
    for(i,j),value in {(1,4):a*X[0]**2,(2,3):-s.Rational(2,9)*X[0]**2,
                       (2,4):s.Rational(4,3)*X[0]*X[1],(3,4):s.Rational(2,3)*X[0]*X[2]+X[1]**2}.items():Q[i,j]=value;Q[j,i]=-value
    frame=s.eye(5);frame[1:,1:]=s.Matrix([[1,0,0,0],[0,0,0,1/a],[0,1,0,0],[0,0,-s.Rational(9,2),0]])
    normalized=model.normalize(pushforward(Q,frame));assert normalized['parameters']==parameters and normalized['perspective']==[0]*4
    model.verify_point(normalized)
    z=s.symbols('y0:4');D=coords();P=affine(model.reconstruct(parameters),s.Matrix(z))
    linear=s.Matrix(4,4,lambda i,j:sum(s.Poly(P[i,j],*z).coeff_monomial(x)*x for x in z))
    alpha=[s.diff(linear[0,1]+linear[2,3],x)/5 for x in z];assert alpha==[0]*4
    def closedjet(P):
        ll=s.Matrix([s.Poly(P[i,j],*z).coeff_monomial(x) for i,j in pairs for x in z])
        qq=s.Matrix([s.Poly(P[i,j],*z).coeff_monomial(z[k]*z[l]) for i,j in pairs for k,l in m2])
        return list(ll[D['li'],:])+list(D['Pi']*qq)
    values=closedjet(P)
    assert [(i,s.factor(v))for i,v in enumerate(values)if v]==[(11,-s.Rational(4,3)/a),(15,3/a),(41,s.Rational(9,2)/a)]
    lifts.W=D['weights'];lifts.ACT=[lifts.action_lists(matrix)for matrix in D['acts']]
    lifts.sel=D['li'];lifts.P4=np.array((4*D['Pi']).tolist(),dtype=object)
    # Independently calibrate the weight1/weight2 scaling in the frozen coframe
    # relations.  This avoids confusing Pi with the producer's4Pi after its
    # (10detR)^2 denominator clearance.
    actual=affine_bivector(s.eye(4)).subs(dict(zip(Y,z)),simultaneous=True)
    actualjet=closedjet(actual)
    eye=np.eye(4,dtype=int).astype(object);scaled=lifts.jets(eye,eye,1)
    assert list(scaled)==[10*x for x in actualjet[:16]]+[100*x for x in actualjet[16:]]
    relationpath=ROOT/'reports/TORIC_ZERO_FIBER_COFRAME_QUADRICS_20261102.json'
    stored=json.loads(relationpath.read_text());assert stored['all_passed']
    highest={weight:[{tuple(mon):s.Rational(c)for mon,c in poly}for poly in polys]
             for weight,polys in stored['highest_lifts_in_projected_coordinates'].items()}
    relations=lifts.orbit(highest);assert len(relations)==274
    results=[]
    for polynomial in relations:
        value=s.factor(sum(c*s.prod(values[i]for i in mon) for mon,c in polynomial.items()if all(values[i]for i in mon)))
        assert value==0;results.append(str(value))
    e2path=ROOT/'reports/TORIC_POLYNOMIAL_E2_COORDINATE_20261101.json';e2=json.loads(e2path.read_text())
    us=s.symbols('u1:62');sub=dict(zip(us,parameters))
    lam=[s.factor(s.sympify(value).subs(sub))for value in e2['lambda_coefficients_q']]
    assert [(i,value)for i,value in enumerate(lam)if value]==[(6,288/a**2)]
    center=s.factor(s.sympify(e2['lambda_center']).subs(sub));assert center==0
    two=e2['weight2_u_indices_zero_based'];Xpol=s.factor(sum(value*parameters[index]for index,value in zip(two,lam))-center);assert Xpol==0
    zz,uu,beta=primitive_beta(model);bs=beta.subs(sub).subs(a,1)
    vector=[s.diff(bs[i,j],x)for i,j in pairs for x in zz]
    c4,c6,delta,j=curve_invariants_from_beta(vector);assert c4==c6==delta==0
    report=dict(id='Q_T01_GLOBAL_CONSTRAINT_CHECKS_20261103_V63',all_passed=True,coefficient_field='QQ(a)',
                full150_S61_normalization_134equations_and_projective_Jacobi_rechecked=True,
                normalization_denominators_only_a=True,primitive_perspective_alpha_zero=True,
                closedjet_nonzero_zero_based={'11':'-4/(3*a)','15':'3/a','41':'9/(2*a)'},
                single_coframe_scaling10_and100_independently_calibrated=True,
                all274_globally_valid_T01_lifts_identically_zero=results,
                polynomial_root_coordinate_X_identically_zero=True,
                lambda_nonzero_qordinal6='288/a^2',lambda_center_zero=True,
                canonical_beta_c4_c6_Delta_identically_zero=True,
                invariants_checked_at_a1_then_homogeneous_scaling=True,
                membership_NOT_inferred_from_necessary_constraint_vanishing=True,
                Delta_open_cover_formulas_not_used=True,
                input_sha256={str(path.resolve().relative_to(ROOT)):sha(path)for path in
                              [Path(__file__),local,relationpath,e2path,ROOT/'code/toric_closed_jet_coordinates_20261102.py',
                               ROOT/'code/toric_zero_fiber_coframe_quadrics_20261102.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
                               ROOT/'code/elliptic_curve_invariants_20261101.py',ROOT/'code/s61_model.py']},
                seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS all274globalconstraints/vacuousE2, singularbeta; seconds',report['seconds'])


if __name__=='__main__':main()

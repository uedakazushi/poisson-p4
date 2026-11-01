#!/usr/bin/env python3
"""Universal Pfaffian lift and orthogonal duality for the null-correlation curve.

Polynomial identities hold for all16 primitive first jets. The genus/Picard
interpretation uses the bundle proof and Shinder--Zhang Proposition2.8, not the
finite invariant samples alone.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib,json,time
from itertools import combinations
from pathlib import Path
import sympy as s
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from q51_universal_section_20261101 import exterior_derivative
from s61_model import S61
from weight1_structure import J
from elliptic_curve_invariants_20261101 import pfaff4,pfaffOmega,invOmega,curve_invariants_from_beta
ROOT=Path(__file__).resolve().parents[1]


def verify():
    started=time.monotonic();model=S61();z,u,beta=primitive_beta(model);y=s.Matrix(z)
    alpha=-beta*y
    assert (exterior_derivative(alpha,y)-3*beta).applyfunc(s.expand)==s.zeros(4)
    assert s.expand(beta[0,1]+beta[2,3])==0
    f=sum(a*y[i] for i,a in enumerate(s.symbols('f0:4')))
    gauge=exterior_derivative(f*J*y,y)/3
    assert s.expand(gauge[0,1]+gauge[2,3]+5*f/3)==0

    t=s.Symbol('t');block=beta+t*J
    phi=block.row_join(J*y).col_join((-(J*y).T).row_join(s.zeros(1)))
    pfaffians=s.Matrix([(-1)**i*pfaff4(phi.extract([j for j in range(5) if j!=i],
        [j for j in range(5) if j!=i])) for i in range(5)])
    quadratic=(beta-t*J)*y
    pfaffian=s.expand(pfaff4(beta))
    assert (pfaffians[:4,0]+J*quadratic).applyfunc(s.expand)==s.zeros(4,1)
    assert s.expand(pfaffians[4]-pfaffian-t*t)==0
    assert ((J*beta)**2+pfaffian*s.eye(4)).applyfunc(s.expand)==s.zeros(4)
    # The fifth quadric is implied by the first four on every y_i != 0 chart.
    assert ((J*beta-t*s.eye(4))*J*quadratic+(pfaffian+t*t)*y).applyfunc(s.expand)==s.zeros(4,1)
    assert pfaffians[4].subs(dict.fromkeys(y,0))==t*t

    r=s.symbols('r0:5');pairs=list(combinations(range(4),2));Z=s.zeros(4)
    for (i,j),value in zip(pairs,[*r,-r[0]],strict=True):Z[i,j]=value;Z[j,i]=-value
    assert sum(J[i,j]*Z[i,j] for i,j in pairs)==0
    ell=s.Matrix([sum(beta[i,j].diff(z[k])*Z[i,j] for i,j in pairs) for k in range(4)])
    w=-J*ell;K=Z.row_join(w).col_join((-w.T).row_join(s.zeros(1)))
    pairing=s.expand(sum(K[i,j]*phi[i,j] for i,j in combinations(range(5),2)))
    assert pairing==0
    star=J*Z*J
    canonical=star.row_join(ell).col_join((-ell.T).row_join(s.zeros(1)))
    congruence=s.diag(J,s.ones(1))
    assert (congruence*canonical*congruence.T+K).applyfunc(s.expand)==s.zeros(5)

    # Optional numeric normalization cross-check, not the duality proof.
    samples=[]
    for values in [list(range(1,17)),[(-1)**i*(i*i+1) for i in range(16)]]:
        substitution=dict(zip(u,values,strict=True));actual=phi.subs(substitution)
        denominator=s.ilcm(*(s.denom(coefficient) for value in actual for coefficient in s.Poly(value,*y,t).coeffs()))
        omega,_=pfaffOmega(denominator*actual,[*y,t]);raw=invOmega(omega,[*y,t])
        # Scaling all entries of a degree5 Pfaffian matrix has c4,c6 weights20,30.
        result=(raw[0]/denominator**20,raw[1]/denominator**30,raw[2]/denominator**60,raw[3])
        coefficients=s.Matrix([beta[i,j].diff(z[k]).subs(substitution) for i,j in pairs for k in range(4)])
        expected=curve_invariants_from_beta(coefficients)
        assert result==expected and result[2]!=0
        samples.append({'first_jet_u_values':values,'c4':str(result[0]),'c6':str(result[1]),
                        'Delta':str(result[2]),'j':str(result[3]),'equals_canonical_model':True})
    inputs=[Path(__file__),ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
            ROOT/'code/q51_universal_section_20261101.py',ROOT/'code/s61_model.py',
            ROOT/'code/weight1_structure.py',ROOT/'code/elliptic_curve_invariants_20261101.py',
            ROOT/'reports/PRIMITIVE_JET_ELLIPTIC_CURVE_20261101.json',
            ROOT/'research/notes/S61_UNIVERSAL_WEIGHT6_REVIEW_20261031.md']
    return {'id':'NULL_CORRELATION_DUALITY_20261101','all_passed':True,
        'generated_utc':datetime.now(timezone.utc).isoformat(),'field':'QQ[all16 primitive first jets]',
        'null_correlation_sequence':'0 -> O(1) -> Omega1_P3(3) -> N(2) -> 0',
        'section_isomorphism':'beta -> class(i_E beta), H0 dimension16',
        'Cartan_identity_d_iE_beta_equals3beta':True,'primitive_kernel_trace_factor':'-5/3',
        'Phi_E_formula':'[[beta(y)+t*J,J*y],[-(J*y)^T,0]]',
        'Pfaffian_vector_first4':'-J*(beta(y)-t*J)*y',
        'Pfaffian_last':'t^2+Pf(beta(y))','all_Pfaffian_identities_checked':True,
        'projection_center_absent':True,'projection_scheme_inverse_on_each_chart':
            'On (J*y)_i != 0, t=(beta(y)*y)_i/(J*y)_i; the fifth quadric follows from (J*beta)^2=-Pf(beta)*I.',
        'K_graph_formula':'[[Z,-J*ell_beta(Z)],[(-J*ell_beta(Z))^T*(-1),0]]; <J,Z>=0',
        'universal_annihilator_pairing_zero':True,'both_linear_spaces_dimension':5,
        'canonical_Phi_congruence':'diag(J,1)*Phi_C*diag(J,1)^T = -K_graph',
        'generic_zero_curve':{'degree':5,'genus':1,'connected':True,'smooth':True,
            'ambient':'P3','linearly_normal':False,'normal_embedding':'Phi_E in P4'},
        'duality_theorem':'Shinder-Zhang, Proposition2.8, BLMS52(2020),395-409',
        'primary_url':'https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/blms.12339',
        'duality_conclusion':'On the smooth transverse locus, E_beta is Pic^3(C_beta), equivalently Pic^2 via the degree5 polarization. Their Jacobians agree; no canonical generic-field isomorphism E_beta=C_beta is asserted.',
        'invariant_samples':samples,
        'G6_scope':'Existing weight6 review section6 proves stability of the entire135-leading zero scheme modulo134, before its flag application; no new saturation claim.',
        'proof_note':'research/notes/NULL_CORRELATION_DUALITY_20261101.md',
        'elapsed_seconds':round(time.monotonic()-started,3),
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as handle:json.dump(report,handle,indent=2);handle.write('\n')
    print(json.dumps({'all_passed':True,'duality':'Pic^3','seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()

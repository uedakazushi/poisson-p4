#!/usr/bin/env python3
"""Construct the ENTIRE Qa family inside actualT01 by commuting matrix powers.

The proof uses a polynomial canonical Poisson arc, explicit diagonalization
away fromt=0, and the definition of the toric component.  No Delta-open formula
or specialization of an unsaturated generic ideal is used.
"""
import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as sp
from sympy.polys.rings import ring

from q_weighted_finite_model_validation_20261103_v63 import physical_equations

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(5),2))
MONS=list(itertools.combinations_with_replacement(range(5),2))


def matrix_zero(zero):return [[zero for _ in range(5)]for _ in range(5)]
def identity(zero,one):return [[one if i==j else zero for j in range(5)]for i in range(5)]
def multiply(A,B,zero):return [[sum((A[i][k]*B[k][j]for k in range(5)),zero)for j in range(5)]for i in range(5)]
def plus(A,B):return [[A[i][j]+B[i][j]for j in range(5)]for i in range(5)]
def scale(A,c):return [[c*value for value in row]for row in A]
def equal(A,B):return all(A[i][j]==B[i][j]for i in range(5)for j in range(5))


def wedge_coefficients(A,B):
    result=[]
    for i,j in PAIRS:
        for p,q in MONS:
            value=A[i][p]*B[j][q]-B[i][p]*A[j][q]
            if p!=q:value+=A[i][q]*B[j][p]-B[i][q]*A[j][p]
            result.append(value)
    return result


def combine(left,right,a,b):return [a*x+b*y for x,y in zip(left,right)]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    a,s=sp.symbols('a s');K=sp.QQ.frac_field(a,s);R,t=ring('t',K)
    ground=lambda expression:K.from_sympy(sp.cancel(expression))
    lift=lambda expression:R(ground(expression))
    class Domain:
        zero=R.zero;one=R.one
        @staticmethod
        def q(value):
            value=Fraction(value)
            return R(K.convert(sp.Rational(value.numerator,value.denominator)))
        @staticmethod
        def add(left,right):return left+right
        @staticmethod
        def mul(left,right):return left*right
        @staticmethod
        def neg(value):return -value
    I=identity(R.zero,R.one);links=[sp.S.One,s,(s+3)/2,-9/(s+3)]
    N=matrix_zero(R.zero)
    for i,value in enumerate(links):N[i+1][i]=lift(value)
    powersN=[I]
    for _ in range(5):powersN.append(multiply(powersN[-1],N,R.zero))
    assert all(not value for row in powersN[5]for value in row)
    assert powersN[4][4][0]==lift(-9*s/2)
    Nt=[[N[i][j]+(i*t if i==j else R.zero)for j in range(5)]for i in range(5)]
    powers=[I]
    for _ in range(4):powers.append(multiply(powers[-1],Nt,R.zero))
    tau=[sum((powers[k][i][i]for i in range(5)),R.zero)*lift(sp.Rational(1,5))for k in range(1,5)]
    expected_tau=[2*t,6*t*t,20*t**3,lift(sp.Rational(354,5))*t**4]
    assert tau==expected_tau
    W=[plus(powers[k],scale(I,-tau[k-1]))for k in range(1,5)]
    assert all(sum((matrix[i][i]for i in range(5)),R.zero)==0 for matrix in W)
    for i in range(4):
        for j in range(4):assert equal(multiply(W[i],W[j],R.zero),multiply(W[j],W[i],R.zero))
    lam=lift(-2*a/(9*s));mu=lift(-4/(9*s*s*(s+3)))
    raw=combine(wedge_coefficients(W[0],W[3]),wedge_coefficients(W[1],W[2]),lam,mu)
    assert len(raw)==150
    # Independent physical differentiation/multiplication in x0,...,x4, with
    # full QQ(a,s)[t] coefficients: no truncation in either x ort.
    divergence,jacobi=physical_equations(raw,Domain)
    assert not any(divergence)and not any(jacobi)
    plain=combine(wedge_coefficients(powers[1],powers[4]),wedge_coefficients(powers[2],powers[3]),lam,mu)
    Z=plus(scale(plus(scale(powers[1],tau[3]),scale(powers[4],-tau[0])),lam),
           scale(plus(scale(powers[2],tau[2]),scale(powers[3],-tau[1])),mu))
    euler=wedge_coefficients(I,Z)
    assert all(c-p-e==0 for c,p,e in zip(raw,plain,euler))
    # All150 limit errors before imposing the quadratic cover relation.
    target=[K.zero]*150
    entries={((1,4),(0,0)):a,((2,3),(0,0)):-sp.Rational(2,9),
             ((2,4),(0,1)):sp.Rational(4,3),((3,4),(0,2)):sp.Rational(2,3),((3,4),(1,1)):1}
    basis=list(itertools.product(PAIRS,MONS))
    for key,value in entries.items():target[basis.index(key)]=ground(value)
    errors=[poly.get((0,),K.zero)-expected for poly,expected in zip(raw,target)]
    F=3*a*s*s+(9*a-4)*s-6
    assert {i for i,value in enumerate(errors)if value}=={121,137}
    assert errors[121]==ground(F/(3*(s+3)))and errors[137]==ground(F/(6*s))
    a_of_s=(4*s+6)/(3*s*(s+3))
    assert all(sp.cancel(value.as_expr().subs(a,a_of_s))==0 for value in errors)
    assert sp.expand(F.subs(s,0))==-6 and sp.expand(F.subs(s,-3))==6
    assert sp.expand(F.subs(a,0))==-4*s-6 and sp.cancel(F.subs({a:0,s:-sp.Rational(3,2)}))==0
    # Fixed e0 remains a rank4 point whenevera!=0, independently oft.
    const=lambda i,j:raw[basis.index(((i,j),(0,0)))]
    pf_e0=const(1,2)*const(3,4)-const(1,3)*const(2,4)+const(1,4)*const(2,3)
    assert pf_e0==lift(-2*a/9)
    # Explicit eigenframe.  It is allowed to have t-poles; rawarc above has
    # polynomialt coefficients and is regular atzero.
    tt=sp.Symbol('t');L=sp.QQ.frac_field(a,s,tt)
    cv=lambda expression:L.from_sympy(sp.cancel(expression))
    NtL=[[cv(value.as_expr())for value in row]for row in Nt]
    SL=identity(L.zero,L.one)
    for i in range(5):
        for j in range(i):SL[i][j]=cv((-1)**(i-j)*sp.prod(links[k]for k in range(j,i))/(sp.factorial(i-j)*tt**(i-j)))
    diag=matrix_zero(L.zero)
    for i in range(5):diag[i][i]=cv(i*tt)
    assert equal(multiply(NtL,SL,L.zero),multiply(SL,diag,L.zero))
    U=plus(SL,scale(identity(L.zero,L.one),-L.one));Sinv=identity(L.zero,L.one);power=identity(L.zero,L.one)
    for k in range(1,5):power=multiply(power,U,L.zero);Sinv=plus(Sinv,scale(power,L((-1)**k)))
    assert equal(multiply(Sinv,SL,L.zero),identity(L.zero,L.one))
    Wdiag=[]
    for k in range(1,5):
        actual=multiply(multiply(Sinv,[[cv(value.as_expr())for value in row]for row in W[k-1]],L.zero),SL,L.zero)
        expected=matrix_zero(L.zero);mean=sp.Rational(sum(i**k for i in range(5)),5)
        for i in range(5):expected[i][i]=cv(tt**k*(i**k-mean))
        assert equal(actual,expected);Wdiag.append(actual)
    diagonal_raw=combine(wedge_coefficients(Wdiag[0],Wdiag[3]),wedge_coefficients(Wdiag[1],Wdiag[2]),cv(-2*a/(9*s)),cv(-4/(9*s*s*(s+3))))
    assert len(diagonal_raw)==150
    gamma={}
    for i,j in itertools.combinations(range(1,5),2):
        gamma[i,j]=lam*(i*j**4-i**4*j)+mu*(i*i*j**3-i**3*j*j)
    pf_gamma=gamma[1,2]*gamma[3,4]-gamma[1,3]*gamma[2,4]+gamma[1,4]*gamma[2,3]
    assert pf_gamma==288*lam*mu
    c12=sp.factor(gamma[1,2].as_expr().subs(a,a_of_s));generic_pf=sp.factor(pf_gamma.as_expr().subs(a,a_of_s))
    assert c12!=0 and generic_pf!=0
    # The untraced diagonal coefficients are exactly the standard6 log
    # coefficients onindices1,...,4; their projective class is the canonicalarc.
    Dpowers=[identity(L.zero,L.one)]
    for _ in range(4):Dpowers.append(multiply(Dpowers[-1],diag,L.zero))
    log_raw=combine(wedge_coefficients(Dpowers[1],Dpowers[4]),wedge_coefficients(Dpowers[2],Dpowers[3]),cv(-2*a/(9*s)),cv(-4/(9*s*s*(s+3))))
    expected_log=[L.zero]*150
    for (i,j),coefficient in gamma.items():expected_log[basis.index(((i,j),(i,j)))]=cv(coefficient.as_expr())*cv(tt**5)
    assert log_raw==expected_log
    # Existing normalizedzero-section coordinates are only used to identify
    # the same Qa inS61; a=0 is retained as a globalrank2 boundary point.
    qlocal=ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63/certificate.json'
    normal=json.loads(qlocal.read_text())['scalar_section_normal61']
    assert [(i+1,text)for i,text in enumerate(normal)if text!='0']==[(25,'-4/(3*a)'),(36,'3/a'),(38,'9/(2*a)')]
    catalogue=ROOT/'legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json'
    toric=next(record for record in json.loads(catalogue.read_text())['diagrams']if record['id']=='T01')
    assert toric['edge_count']==0 and toric['period_dimension']==6
    assert sp.Matrix(toric['period_basis']).rank()==6 and toric['affine_component_dimension']==26
    weights=json.loads((ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json').read_text())['parameter_weights']
    weighted=[sp.cancel(a**weight*sp.sympify(value))for value,weight in zip(normal,weights)]
    assert [(i+1,value)for i,value in enumerate(weighted)if value]!=[]
    assert [(i+1,value)for i,value in enumerate(weighted)if value]==[(25,-sp.Rational(4,3)),(36,3),(38,9*a/2)]
    artifacts=dict(cover_relation=str(F),a_of_s=str(a_of_s),N=[[str(value.as_expr())for value in row]for row in N],
                   lambda_parameter=str(lam.as_expr()),mu_parameter=str(mu.as_expr()),
                   canonical_arc_raw150=[[[mon[0],str(c.as_expr())]for mon,c in sorted(value.items())]for value in raw],
                   limit_error_raw150={str(i):str(value.as_expr())for i,value in enumerate(errors)if value},
                   eigenframe=[[str(value.as_expr())for value in row]for row in SL],
                   eigenframe_inverse=[[str(value.as_expr())for value in row]for row in Sinv],
                   toric_periods={f'{i},{j}':str(value.as_expr())for (i,j),value in gamma.items()},
                   canonical_trace_corrections=list(map(str,[value.as_expr()for value in tau])),
                   fixed_nonzero_beta_line=dict(u25='-4/3',u36='3',u38='c',all_other_u='0'))
    certificate=args.output/'certificate.json';certificate.write_text(json.dumps(artifacts,indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'code/q_weighted_finite_model_validation_20261103_v63.py',ROOT/'code/s61_model.py',
            qlocal,catalogue,ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json',
            ROOT/'research/notes/TORIC_AFFINE_MODIFICATION_20261101.md',ROOT/'research/notes/TORIC_STRICT_STABILIZER_20261101.md']
    report=dict(id='Q_T01_MEMBERSHIP_20261103_V63',all_passed=True,
                coefficient_ring='QQ(a,s)[t], coverF=0; equivalently QQ[s,1/s,1/(s+3)][t]',
                complete_Qa_family_in_actual_T01=True,all_complex_a_including_zero=True,
                nonzero_a_normalized_S61_zero_section_in_actual_T01=True,
                entire_affine_qline_over_fixed_nonzero_beta_in_actual_T01=True,
                proof_method='Explicit regularcanonicalcommuting-matrixarc andstandard6-parametertoricclosure',
                saturation_not_specialized=True,Delta_inverse_or_smooth_beta_assumption_used=False,
                N_regular_nilpotent_N4_entry='-9*s/2',commuting_tracefree_powers=True,
                full150_limit_identity_verified_with_explicit_cover_errors=True,
                physical_divergence25_and_rawJacobi350_identically_zero=True,
                full150_canonical_vs_plain_Euler_wedge_identity=True,
                explicit_unitriangular_eigenframe_and_inverse_verified=True,
                full150_diagonal_toric_coefficient_identity=True,
                toric_period_pfaffian='288*lambda*mu*t^10',generic_toric_chart_c12=c12.__str__(),
                generic_toric_chart_pfaffian=str(generic_pf),
                cover_excludes_s_zero_and_minus3_for_every_a=True,
                a_zero_admissible_cover_point='s=-3/2',
                fixed_e0_projective_value_pfaffian='-2*a/9, independentoft',
                generic_cover_degree2=True,discriminant_in_s='81*a^2+16',
                full_T01_image_or_local_branch_exhaustion_claimed=False,
                certificate_sha256=sha(certificate),
                input_sha256={str(path.resolve().relative_to(ROOT)):sha(path)for path in inputs},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS entireQa inactualT01; full150+350 canonicalarc, eigenframe; seconds',report['seconds'])


if __name__=='__main__':main()

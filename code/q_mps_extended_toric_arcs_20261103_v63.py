#!/usr/bin/env python3
"""Four actual MPS components through ALL Qa, using vanishing arc terms."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as sp
from sympy.polys.rings import ring

from q_t01_membership_20261103_v63 import PAIRS,MONS,identity,matrix_zero,multiply,plus,scale,equal,wedge_coefficients
from q_weighted_finite_model_validation_20261103_v63 import physical_equations
from q_mps_toric_arc_memberships_20261103_v63 import PAIRS4,annihilator,period_permutations

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def arc(record):
    parameter=sp.Symbol(record['free']);K=sp.QQ.frac_field(parameter);R,t=ring('t',K)
    convert=lambda expr:K.from_sympy(sp.cancel(expr));lift=lambda expr:R(convert(expr))
    r=sp.sympify(record['rho']);alpha=sp.sympify(record['alpha']);beta=sp.sympify(record['beta'])
    s=sp.Rational(3,2)*(r-1);a=8*r/(9*(r*r-1));mu=-4/(9*s*s*(s+3));lam=r*mu
    assert sp.cancel(sp.diff(r,parameter))!=0 and sp.cancel(sp.diff(a,parameter))!=0
    assert all(sp.cancel(v)!=0 for v in[s,s+3,r,r-1,r+1])
    di=list(map(sp.sympify,record['eigenvalues']));assert all(sp.cancel(di[i]-di[j])!=0 for i in range(5)for j in range(i))
    links=[sp.S.One,s,(s+3)/2,-9/(s+3)];I=identity(R.zero,R.one);N=matrix_zero(R.zero)
    for i,value in enumerate(links):N[i+1][i]=lift(value)
    power=I
    for k in range(1,6):
        power=multiply(power,N,R.zero)
        if k==4:assert power[4][0]==lift(-9*s/2)
    assert all(not value for row in power for value in row)
    Nt=[[N[i][j]+(t*lift(di[i])if i==j else R.zero)for j in range(5)]for i in range(5)]
    powers=[I]
    for k in range(4):powers.append(multiply(powers[-1],Nt,R.zero))
    tau=[sum((powers[k][i][i]for i in range(5)),R.zero)*lift(sp.Rational(1,5))for k in range(1,5)]
    W=[plus(powers[k],scale(I,-tau[k-1]))for k in range(1,5)]
    for A in W:assert sum((A[i][i]for i in range(5)),R.zero)==0
    for i in range(4):
        for j in range(i):assert equal(multiply(W[i],W[j],R.zero),multiply(W[j],W[i],R.zero))
    terms=[(0,1,lift(mu*alpha)*t*t),(0,2,lift(mu*beta)*t),(0,3,lift(lam)),(1,2,lift(mu))]
    raw=[R.zero]*150;plain=[R.zero]*150;Z=matrix_zero(R.zero)
    for i,j,c in terms:
        raw=[a+c*b for a,b in zip(raw,wedge_coefficients(W[i],W[j]))]
        plain=[a+c*b for a,b in zip(plain,wedge_coefficients(powers[i+1],powers[j+1]))]
        Z=plus(Z,scale(plus(scale(powers[i+1],tau[j]),scale(powers[j+1],-tau[i])),c))
    assert all(a-b-c==0 for a,b,c in zip(raw,plain,wedge_coefficients(I,Z)))
    target=[K.zero]*150;basis=list(itertools.product(PAIRS,MONS))
    for key,value in {((1,4),(0,0)):a,((2,3),(0,0)):-sp.Rational(2,9),((2,4),(0,1)):sp.Rational(4,3),((3,4),(0,2)):sp.Rational(2,3),((3,4),(1,1)):1}.items():target[basis.index(key)]=convert(value)
    assert all(poly.get((0,),K.zero)==value for poly,value in zip(raw,target))
    class Domain:
        zero=R.zero;one=R.one
        @staticmethod
        def q(v):return lift(sp.Rational(v.numerator,v.denominator))
        @staticmethod
        def add(a,b):return a+b
        @staticmethod
        def mul(a,b):return a*b
        @staticmethod
        def neg(a):return -a
    divergence,jacobi=physical_equations(raw,Domain)
    assert not any(divergence)and not any(jacobi)
    # Explicit eigenframe recurrence; its unit diagonal gives determinant1.
    S=[[sp.S.One if i==j else sp.S.Zero for j in range(5)]for i in range(5)]
    for i in range(5):
        for j in range(i):
            S[i][j]=sp.factor((-1)**(i-j)*sp.prod(links[k]for k in range(j,i))/sp.prod(di[k]-di[j]for k in range(j+1,i+1)))
            assert sp.cancel((di[i]-di[j])*S[i][j]+links[i-1]*S[i-1][j])==0
    return dict(component=record['component'],parameter=str(parameter),a=str(sp.factor(a)),
        a_derivative=str(sp.factor(sp.diff(a,parameter))),rho=str(r),alpha=str(alpha),beta=str(beta),
        N_links=list(map(str,links)),diagonal=list(map(str,di)),
        coefficient_wedges='mu*(alpha*t^2 W1^W2+beta*t W1^W3+rho W1^W4+W2^W3)',
        canonical_raw150=[[[int(mon[0]),str(coefficient.as_expr())]for mon,coefficient in sorted(poly.items())]for poly in raw],
        eigenframe_coefficients=[[str(c)for c in row]for row in S],
        eigenframe_rule='Sij=coefficient/t^(i-j)',
        all150_Qa_limit=True,all150_Euler_identity=True,all25_divergence_and350_raw_Jacobi=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    witnesspath=ROOT/'code/data/q_mps_extended_toric_arc_witnesses_20261103_v63.json';witnesses=json.loads(witnesspath.read_text())
    catalog=ROOT/'legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json';records=json.loads(catalog.read_text())['diagrams']
    M=sp.Matrix(period_permutations()[0].tolist());results=[];arcs=[]
    for witness in witnesses:
        rec=next(v for v in records if v['id']==witness['component']);assert witness['permutation']==[0,1,2,3,4]
        di=list(map(sp.sympify,witness['eigenvalues']));rho=sp.sympify(witness['rho']);alpha=sp.sympify(witness['alpha']);beta=sp.sympify(witness['beta'])
        G=sp.Matrix([[di[i]*di[j]**2-di[i]**2*di[j],di[i]*di[j]**3-di[i]**3*di[j],di[i]*di[j]**4-di[i]**4*di[j],di[i]**2*di[j]**3-di[i]**3*di[j]**2]for i,j in PAIRS4])
        L=sp.Matrix(annihilator(rec).tolist());K=L*M*G;coefficients=sp.Matrix([alpha,beta,rho,1])
        assert (K*coefficients).applyfunc(sp.cancel)==sp.zeros(K.rows,1)
        gamma=G*coefficients;period=M*gamma;basis=sp.Matrix(rec['period_basis']).applyfunc(sp.Rational).T
        selected=list(basis.T.rref()[1]);period_coordinates=basis[selected,:].inv()*period[selected,:]
        assert (basis*period_coordinates-period).applyfunc(sp.cancel)==sp.zeros(6,1)
        pf=gamma[0]*gamma[5]-gamma[1]*gamma[4]+gamma[2]*gamma[3]
        vandermonde=sp.prod(di[i]for i in range(1,5))*sp.prod(di[j]-di[i]for i in range(1,5)for j in range(i+1,5))
        assert sp.cancel(pf-rho*vandermonde)==0 and sp.cancel(pf)!=0
        size=K.rows;assert size in[2,3];det=sp.factor(K[:,:size].det());assert det!=0
        result=dict(component=rec['id'],MPS_internal_component_id=True,projective_component_dimension=rec['projective_component_dimension'],
            original_period_basis=rec['period_basis'],identity_label_permutation=True,
            exact_period_constraint_matrix=[[str(sp.factor(c))for c in row]for row in K.tolist()],
            coefficient_solve_determinant=str(det),period_coordinates=list(map(lambda c:str(sp.factor(c)),period_coordinates)),
            all_period_equations_exact_zero=True,pfaffian=str(sp.factor(pf)),
            all_finite_complex_Qa_membership_by_rational_density_and_closedness=True)
        results.append(result);arcs.append(arc(witness))
    assert [v['component']for v in results]==['T02','T03','T16','T17']
    (args.output/'period-certificates.json').write_text(json.dumps(results,indent=2)+'\n')
    (args.output/'full150-arcs.json.gz').write_bytes(gzip.compress(json.dumps(arcs,separators=(',',':')).encode(),mtime=0))
    inputs=[Path(__file__),witnesspath,catalog,ROOT/'code/q_mps_toric_arc_memberships_20261103_v63.py',
        ROOT/'code/q_t01_membership_20261103_v63.py',ROOT/'code/q_weighted_finite_model_validation_20261103_v63.py',
        ROOT/'research/notes/Q_MPS_EXTENDED_TORIC_ARCS_20261103_v63.md']
    report=dict(id='Q_MPS_EXTENDED_TORIC_ARCS_20261103_v63',all_passed=True,
        all_finite_complex_Qa_in_actual_MPS_components=['T02','T03','T16','T17'],
        includes_a_zero_by_closedness=True,all150_Qa_limits_and_Euler_identities=True,
        all25_divergence_and350_rawJacobi_per_arc_exact=True,
        explicit_rational_period_parameters_verified=True,
        nonconstant_rho_and_a_maps_verified=True,
        proof_uses_no_modular_specialization=True,
        MPS_T16_T17_not_R19_boundary_stratum_labels=True,
        component_or_local_branch_exhaustion_claimed=False,known_global_component_lower_bound_unchanged=57,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS ALLQa in actualMPS T02,T03,T16,T17;full150+350 each;noexhaustion;',report['seconds'])


if __name__=='__main__':main()

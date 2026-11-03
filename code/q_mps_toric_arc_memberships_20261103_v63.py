#!/usr/bin/env python3
"""Actual ALL-Qa membership in T02,T03 via period-constrained toric arcs.

Fixed-diagonal and one-free-eigenvalue scans use all40*120 labelings.  The
positive two-eigenvalue curve witnesses are sufficient, not an exhaustion.
"""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np
import sympy as sp
from sympy.polys.rings import ring

from q_t01_membership_20261103_v63 import (
    PAIRS,MONS,identity,matrix_zero,multiply,plus,scale,equal,wedge_coefficients,combine,
)

ROOT=Path(__file__).resolve().parents[1]
PAIRS4=list(itertools.combinations(range(1,5),2))
PERMS=list(itertools.permutations(range(5)))


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def dual(v):return [v[5],-v[4],v[3],v[2],-v[1],v[0]]


def extend(v):
    B=np.zeros((5,5),dtype=np.int64)
    for(i,j),c in zip(PAIRS4,v):B[i,j]=c;B[j,i]=-c
    for i in range(1,5):B[0,i]=sum(B[i,1:]);B[i,0]=-B[0,i]
    assert not np.any(B@np.ones(5,dtype=np.int64));return B


def period_permutations():
    result=[]
    for perm in PERMS:
        columns=[];sign=(-1)**sum(perm[i]>perm[j]for i in range(5)for j in range(i+1,5))
        for k in range(6):
            v=[int(i==k)for i in range(6)];B=extend(dual(v));columns.append([int(B[perm[i],perm[j]])for i,j in PAIRS4])
            A=np.zeros((5,5),dtype=np.int64)
            for(i,j),c in zip(PAIRS4,v):A[i,j]=c;A[j,i]=-c
            chart=[int(A[perm[i],perm[j]]-A[perm[i],perm[0]]-A[perm[0],perm[j]])for i,j in PAIRS4]
            assert dual(chart)==[sign*c for c in columns[-1]]
        result.append(np.array(columns,dtype=np.int64).T)
    return result


def annihilator(record):
    basis=sp.Matrix(record['period_basis']).applyfunc(sp.Rational)
    assert basis.rank()==record['period_dimension']
    rows=[]
    for column in basis.nullspace():
        den=sp.ilcm(*[q.q for q in column]);values=[int(q*den)for q in column]
        divisor=math.gcd(*values);rows.append([v//divisor for v in values])
    return np.array(rows,dtype=np.int64).reshape(6-basis.rows,6)


def convolution(a,b):
    answer=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):answer[i+j]+=int(x)*int(y)
    return answer


def modpoly(a,p):
    a=[int(c)%p for c in a]
    while a and not a[-1]:a.pop()
    return a


def modrem(a,b,p):
    a,b=modpoly(a,p),modpoly(b,p);inv=pow(b[-1],-1,p)
    for i in range(len(a)-len(b),-1,-1):
        c=a[i+len(b)-1]*inv%p
        for j,v in enumerate(b):a[i+j]=(a[i+j]-c*v)%p
    return modpoly(a,p)


def modpow(a,n,f,p):
    answer=[1]
    while n:
        if n&1:answer=modrem(convolution(answer,a),f,p)
        n//=2
        if n:a=modrem(convolution(a,a),f,p)
    return answer


def irred_specialization(f,x,y,witness):
    assert sp.gcd_list(sp.Poly(f,x).all_coeffs())==1
    p,yy=witness['prime'],witness['y_value'];degree=sp.degree(f,x);assert degree in[2,5]
    value=sp.Poly(f.subs(y,yy),x,modulus=p);assert value.degree()==degree
    g=[int(value.nth(i))%p for i in range(degree+1)];g=[v*pow(g[-1],-1,p)%p for v in g]
    power=[0,1];first=None
    for k in range(degree):
        power=modpow(power,p,g,p)
        if k==0:first=power
    assert power==[0,1]
    difference=list(first)+[0]*max(0,2-len(first));difference[1]=(difference[1]-1)%p
    a,b=g,modpoly(difference,p)
    while b:a,b=b,modrem(a,b,p)
    assert len(a)==1
    return dict(prime=p,y_value=yy,monic_modular_coefficients=g,
                Rabin_x_p_to_degree_equals_x=True,Rabin_gcd_x_p_minus_x_is_one=True)


def universal_arc():
    rho,x,y=sp.symbols('rho x y');K=sp.QQ.frac_field(rho,x,y);R,t=ring('t',K)
    lift=lambda value:R(K.from_sympy(sp.cancel(value)))
    ss=sp.Rational(3,2)*(rho-1);aa=8*rho/(9*(rho**2-1))
    lam=-2*aa/(9*ss);mu=-4/(9*ss**2*(ss+3));assert sp.cancel(lam/mu-rho)==0
    links=[1,ss,(ss+3)/2,-9/(ss+3)];di=[0,1,2,x,y];I=identity(R.zero,R.one);N=matrix_zero(R.zero)
    for i,c in enumerate(links):N[i+1][i]=lift(c)
    Nt=[[N[i][j]+(lift(di[i])*t if i==j else R.zero)for j in range(5)]for i in range(5)]
    powers=[I]
    for k in range(4):powers.append(multiply(powers[-1],Nt,R.zero))
    tau=[sum((powers[k][i][i]for i in range(5)),R.zero)*lift(sp.Rational(1,5))for k in range(1,5)]
    W=[plus(powers[k],scale(I,-tau[k-1]))for k in range(1,5)]
    for A in W:assert sum((A[i][i]for i in range(5)),R.zero)==0
    for i in range(4):
        for j in range(i):assert equal(multiply(W[i],W[j],R.zero),multiply(W[j],W[i],R.zero))
    raw=combine(wedge_coefficients(W[0],W[3]),wedge_coefficients(W[1],W[2]),lift(lam),lift(mu))
    target=[K.zero]*150;basis=list(itertools.product(PAIRS,MONS))
    for key,c in {((1,4),(0,0)):aa,((2,3),(0,0)):-sp.Rational(2,9),((2,4),(0,1)):sp.Rational(4,3),((3,4),(0,2)):sp.Rational(2,3),((3,4),(1,1)):1}.items():target[basis.index(key)]=K.from_sympy(sp.cancel(c))
    assert all(value.get((0,),K.zero)==c for value,c in zip(raw,target))
    plain=combine(wedge_coefficients(powers[1],powers[4]),wedge_coefficients(powers[2],powers[3]),lift(lam),lift(mu))
    Z=plus(scale(plus(scale(powers[1],tau[3]),scale(powers[4],-tau[0])),lift(lam)),scale(plus(scale(powers[2],tau[2]),scale(powers[3],-tau[1])),lift(mu)))
    assert all(a-b-c==0 for a,b,c in zip(raw,plain,wedge_coefficients(I,Z)))
    S=[[sp.S.One if i==j else sp.S.Zero for j in range(5)]for i in range(5)]
    for i in range(5):
        for j in range(i):
            S[i][j]=sp.factor((-1)**(i-j)*sp.prod(links[k]for k in range(j,i))/sp.prod(di[k]-di[j]for k in range(j+1,i+1)))
            assert sp.cancel((di[i]-di[j])*S[i][j]+links[i-1]*S[i-1][j])==0
    return dict(coefficient_field='QQ(rho,x,y)[t]',s_of_rho=str(ss),a_of_rho=str(aa),lambda_parameter=str(sp.factor(lam)),mu_parameter=str(sp.factor(mu)),
                N_links=list(map(str,links)),diagonal=list(map(str,di)),eigenframe_t_denominator_rule='Sij=listed_coefficient/t^(i-j)',
                eigenframe_coefficients=[[str(v)for v in row]for row in S],
                all150_limit_coefficients_equal_Qa=True,all150_Euler_wedge_difference_checked=True,
                commuting_tracefree_power_identities=True,
                canonical_arc_raw150=[[[int(mon[0]),str(c.as_expr())]for mon,c in sorted(poly.items())]for poly in raw])


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    catalog=ROOT/'legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json';records=json.loads(catalog.read_text())['diagrams'];assert len(records)==40
    maps=period_permutations();L={r['id']:annihilator(r)for r in records}
    # Build the one-variable coefficient arrays; d=4 gives the original scan.
    d=sp.Symbol('d');ds=[0,1,2,3,d]
    coefficients=np.array([[[int(sp.expand(expression).coeff(d,k))for k in range(5)]for expression in
        [ds[i]*ds[j]**4-ds[i]**4*ds[j],ds[i]**2*ds[j]**3-ds[i]**3*ds[j]**2]]for i,j in PAIRS4],dtype=np.int64)
    fixed=coefficients@np.array([4**k for k in range(5)],dtype=np.int64)
    cases=[];fixed_hits=[];generic_nonconstant=[]
    for record in records:
        for perm,M in zip(PERMS,maps):
            K=L[record['id']]@M;eq=K@fixed;kind='none';rho=None
            non=[(int(a),int(b))for a,b in eq if a or b]
            if not non:kind='all_rho'
            else:
                lead=next(((a,b)for a,b in non if a),None)
                if lead:
                    rr=Fraction(-lead[1],lead[0])
                    if all(a*rr+b==0 for a,b in non):kind='finite_rho';rho=rr
                elif all(a==0 for a,b in non):kind='infinity_only'
            case=dict(component=record['id'],permutation=list(perm),kind=kind,constraint_pairs=[list(v)for v in non])
            if rho is not None:
                case['rho']=str(rho);case['allowed_Qa_arc']=rho not in[-1,1]
                if rho not in[-1,1]:
                    aa=8*rho/(9*(rho*rho-1));case['a']=str(aa);case['s']=str(Fraction(3,2)*(rho-1));fixed_hits.append(case)
            cases.append(case)
            A=K@coefficients[:,0,:];B=K@coefficients[:,1,:]
            rows=[(a.tolist(),b.tolist())for a,b in zip(A,B)if np.any(a)or np.any(b)]
            if not rows:generic_nonconstant.append(dict(component=record['id'],permutation=list(perm),kind='free_rho'))
            else:
                lead=next(((a,b)for a,b in rows if any(a)),None)
                if lead and all(convolution(a,lead[1])==convolution(b,lead[0])for a,b in rows):
                    a,b=lead;j=next(i for i,c in enumerate(a)if c)
                    if any(b[k]*a[j]!=b[j]*a[k]for k in range(5)):
                        generic_nonconstant.append(dict(component=record['id'],permutation=list(perm),kind='nonconstant_rational_rho'))
    assert len(cases)==4800 and len(generic_nonconstant)==120 and {v['component']for v in generic_nonconstant}=={'T01'}
    assert {c['component']for c in fixed_hits}=={'T02','T03'}
    assert len([c for c in fixed_hits if c['component']=='T02'])==len([c for c in fixed_hits if c['component']=='T03'])==8
    (args.output/'fixed120scan.json.gz').write_bytes(gzip.compress(json.dumps(cases,separators=(',',':')).encode(),mtime=0))
    (args.output/'fixed-witnesses.json').write_text(json.dumps(fixed_hits,indent=2)+'\n')
    (args.output/'one-variable-generic-scan.json').write_text(json.dumps(dict(cases_checked=4800,
        nonconstant_or_free_rho_cases=generic_nonconstant,
        scope='Generic d over QQ(d) only; no global nonmembership or exhaustion.'),indent=2)+'\n')
    x,y,rho=sp.symbols('x y rho');ds=[0,1,2,x,y]
    gamma=[rho*(ds[i]*ds[j]**4-ds[i]**4*ds[j])+(ds[i]**2*ds[j]**3-ds[i]**3*ds[j]**2)for i,j in PAIRS4]
    period=sp.Matrix(dual(gamma));vandermonde=2*x*y*(x-1)*(x-2)*(y-1)*(y-2)*(y-x)
    pf=gamma[0]*gamma[5]-gamma[1]*gamma[4]+gamma[2]*gamma[3]
    assert sp.expand(pf-rho*vandermonde)==0
    curvepath=ROOT/'code/data/q_mps_toric_arc_curves_20261103_v63.json';witnesses=json.loads(curvepath.read_text());curve_results=[]
    for witness in witnesses:
        rec=next(r for r in records if r['id']==witness['component']);assert witness['permutation']==[0,1,2,3,4]
        f=sp.sympify(witness['curve']);F=sp.Poly(f,x,y);rr=sp.cancel(sp.sympify(witness['rho']));N,D=sp.fraction(rr)
        assert all(c.is_Integer for c in sp.Poly(N,x,y).coeffs()+sp.Poly(D,x,y).coeffs())
        reduction=lambda expression:sp.Poly(sp.expand(expression),x,y).rem(F).as_expr()
        irreducibility=irred_specialization(f,x,y,witness['irreducibility_witness'])
        values=period.applyfunc(lambda value:sp.expand(N*sp.expand(value).coeff(rho)+D*value.subs(rho,0)))
        constraints=sp.Matrix(L[rec['id']].tolist())*values
        assert all(reduction(value)==0 for value in constraints)
        basis=sp.Matrix(rec['period_basis']).applyfunc(sp.Rational).T;chosen=list(basis.T.rref()[1]);parameters=basis[chosen,:].inv()*values[chosen,:]
        assert all(reduction(value)==0 for value in basis*parameters-values)
        tangent=sp.expand((sp.diff(N,x)*D-N*sp.diff(D,x))*sp.diff(f,y)-(sp.diff(N,y)*D-N*sp.diff(D,y))*sp.diff(f,x))
        factors=[vandermonde,D,N,N-D,N+D,tangent,N*N+D*D]
        remainders=[reduction(value)for value in factors];assert all(value!=0 for value in remainders)
        # The verified QQ-irreducibility makes the quotient a domain: each
        # factor nonzero implies their PRODUCT nonzero, with no huge expansion.
        curve_results.append(dict(component=rec['id'],curve=str(f),rho_numerator=str(N),rho_denominator=str(D),
            exact_period_parameters=[str(sp.factor(p/D))for p in parameters],
            irreducibility=irreducibility,period_constraint_remainders_all_zero=True,
            admissibility_and_tangent_nonzero_remainders=list(map(str,remainders)),
            product_nonzero_by_irreducible_domain=True,a_image_dense=True,
            all_finite_complex_Qa_membership_by_closedness=True))
    (args.output/'dense-curve-certificates.json').write_text(json.dumps(curve_results,indent=2)+'\n')
    arc=universal_arc()
    (args.output/'universal-two-eigen-arc.json.gz').write_bytes(gzip.compress(json.dumps(arc,separators=(',',':')).encode(),mtime=0))
    # Qa1 provenance: the previously proved T37 family P(b)+t*tau(b)
    # has tau(1)=Qa1. Its projective t=infinity limit proves this membership.
    legacy=ROOT/'legacy/v23/proofs/literature/ss_geometry_mps_certificate.json';old=json.loads(legacy.read_text())
    assert old['all_passed']and old['MPS_all_brackets_match']and old['catalogue_id']=='T37'
    assert old['checks']['full_MPS_identification_all_a_t']and old['checks']['old_catalog_T37_biresidue']
    (args.output/'scope-and-provenance.json').write_text(json.dumps(dict(
        new_membership='ALL finite complexQa in actual T02 and T03',
        fixed_diagonal_witnesses_are_not_exhaustive_memberships=True,
        Qa1_T37_provenance='v23 fullP(b)+t*tau(b) MPSExample5.6 identity; tau(1)=Qa1; projective t=infinity limit',
        Qa1_T37_is_not_from_the_current_fixed_diagonal_test=True,
        no_other_component_nonmembership_claim=True,no_full_local_branch_claim=True,
        no_new_component_count=True),indent=2)+'\n')
    inputs=[Path(__file__),curvepath,catalog,ROOT/'code/mps_s61_catalogue_20261103_v62.py',
        ROOT/'code/q_t01_membership_20261103_v63.py',ROOT/'reports/Q_T01_MEMBERSHIP_20261103_v63/report.json',legacy,
        ROOT/'legacy/v23/proofs/literature/SS_geometry_and_MPS_identification.txt',
        ROOT/'research/notes/Q_MPS_TORIC_ARC_MEMBERSHIPS_20261103_v63.md']
    report=dict(id='Q_MPS_TORIC_ARC_MEMBERSHIPS_20261103_v63',all_passed=True,
        actual_all_finite_complex_Qa_in_T02=True,actual_all_finite_complex_Qa_in_T03=True,
        all_a_zero_included_by_closedness=True,
        proof='Actual semisimple commuting-matrix arcs, exact period-space membership on nonconstant algebraic curves, density and closedness',
        fixed_diagonal_label_tests=4800,one_variable_generic_label_tests=4800,
        all120_permutation_duality_signs_checked=True,
        period_biresidue_not_confused_with_log_bivector=True,
        full150_generic_two_eigen_arc_limit_and_Euler_identity_checked=True,
        curve_irreducibility_exact_mod5_Rabin_certificates=True,
        all_curve_period_remainders_exact_QQ=True,
        nonmembership_or_component_exhaustion_claimed=False,local_branch_exhaustion_claimed=False,
        known_component_count_lower_bound_unchanged=57,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS ALLQa in actualT02,T03;40x120 fixed andgeneric-one-eigen scans;noexhaustion;',report['seconds'])


if __name__=='__main__':main()

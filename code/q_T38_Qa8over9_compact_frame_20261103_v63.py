#!/usr/bin/env python3
"""Exact quartic-field T38 toric frame and full canonical150 arc to Qa8/9."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.rings import ring
from q_t01_membership_20261103_v63 import identity,matrix_zero,multiply,plus,scale,equal,wedge_coefficients,PAIRS,MONS
from q_weighted_finite_model_validation_20261103_v63 import physical_equations

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_T38_QA8OVER9_COMPACT_FRAME_20261103_v63'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    K=s.QQ.algebraic_field(s.sqrt(25-8*s.sqrt(5)),alias='theta')
    theta=K.unit;q=K.convert
    assert list(map(s.Rational,K.mod.to_list()))==[1,0,-50,0,305]
    # theta^4-50theta²+305 is 5-Eisenstein, so this is a field, not a product algebra.
    assert 50%5==305%5==0 and 305%25!=0
    rho=(q(33)-theta**2)/q(16)
    z=(theta-q(7))*(rho-K.one)/q(4)
    x=-(3*rho+K.one)*z-(2*rho+K.one)
    y=-(2*rho+K.one)*z-(3*rho+K.one)
    assert rho*rho-rho-K.one==K.zero
    assert 2*rho*(z*z+K.one)+7*z==K.zero
    assert y==x+rho*(z-K.one)
    alpha=x*x+rho*z*(x+y)+rho*x*y
    beta=-rho*(z+x+y)-x
    ds=[K.zero,K.one,z,x,y]
    def enc(v):
        coeff=list(reversed(v.to_list()));return [str(c)for c in coeff]+['0']*(4-len(coeff))
    def norm(v):
        columns=[enc(v*theta**i)for i in range(4)]
        return s.Matrix([[s.Rational(columns[j][i])for j in range(4)]for i in range(4)]).det()
    units=[]
    for name,value in [(f'd{i}-d{j}',ds[i]-ds[j])for i in range(5)for j in range(i)]:
        assert value and norm(value)!=0
        inverse=K.one/value;assert value*inverse==K.one
        units.append(dict(name=name,value=enc(value),inverse=enc(inverse),norm=str(norm(value))))
    gamma={}
    for i,j in itertools.combinations(range(1,5),2):
        u,v=ds[i],ds[j]
        gamma[i,j]=u*v*(v-u)*(alpha+beta*(u+v)+rho*(u*u+u*v+v*v)+u*v)
    p0=gamma[1,2]
    assert p0 and gamma[1,3]==p0/2 and gamma[2,4]==-p0/2
    assert gamma[1,4]==gamma[2,3]==gamma[3,4]==K.zero
    pf=gamma[1,2]*gamma[3,4]-gamma[1,3]*gamma[2,4]+gamma[1,4]*gamma[2,3]
    vand=s.prod(ds[i]for i in range(1,5))*s.prod(ds[j]-ds[i]for i in range(1,5)for j in range(i+1,5))
    assert pf==p0*p0/4==rho*vand and pf
    for name,value in [('gamma12',p0),('gamma_Pfaffian',pf)]:
        units.append(dict(name=name,value=enc(value),inverse=enc(K.one/value),norm=str(norm(value))))
    ss=3*(rho-K.one)/2
    mu=q(-4)/(9*ss*ss*(ss+3))
    aa=8*rho/(9*(rho*rho-K.one));assert aa==q(s.Rational(8,9))
    links=[K.one,ss,(ss+3)/2,q(-9)/(ss+3)]
    R,t=ring('t',K);lift=lambda v:R(v)
    I=identity(R.zero,R.one);N=matrix_zero(R.zero)
    for i,c in enumerate(links):N[i+1][i]=lift(c)
    Nt=[[N[i][j]+(t*lift(ds[i])if i==j else R.zero)for j in range(5)]for i in range(5)]
    powers=[I]
    for _ in range(4):powers.append(multiply(powers[-1],Nt,R.zero))
    trace=[sum((powers[k][i][i]for i in range(5)),R.zero)/5 for k in range(1,5)]
    W=[plus(powers[k],scale(I,-trace[k-1]))for k in range(1,5)]
    assert all(sum((A[i][i]for i in range(5)),R.zero)==R.zero for A in W)
    for i in range(4):
        for j in range(i):assert equal(multiply(W[i],W[j],R.zero),multiply(W[j],W[i],R.zero))
    raw=[R.zero]*150
    for i,j,c in [(0,1,lift(mu*alpha)*t*t),(0,2,lift(mu*beta)*t),(0,3,lift(mu*rho)),(1,2,lift(mu))]:
        raw=[a+c*b for a,b in zip(raw,wedge_coefficients(W[i],W[j]))]
    target=[K.zero]*150;basis=list(itertools.product(PAIRS,MONS))
    for key,value in {((1,4),(0,0)):aa,((2,3),(0,0)):q(s.Rational(-2,9)),((2,4),(0,1)):q(s.Rational(4,3)),((3,4),(0,2)):q(s.Rational(2,3)),((3,4),(1,1)):K.one}.items():target[basis.index(key)]=value
    assert all(poly.get((0,),K.zero)==value for poly,value in zip(raw,target))
    class Domain:
        zero=R.zero;one=R.one
        @staticmethod
        def q(v):return lift(K.convert(s.Rational(v.numerator,v.denominator)))
        @staticmethod
        def add(a,b):return a+b
        @staticmethod
        def mul(a,b):return a*b
        @staticmethod
        def neg(a):return -a
    divergence,jacobi=physical_equations(raw,Domain)
    assert not any(divergence)and not any(jacobi)
    # Laurent eigenframe Sij=S0ij*t^(j-i), exact over K(t).
    S0=identity(K.zero,K.one)
    for i in range(5):
        for j in range(i):S0[i][j]=-links[i-1]*S0[i-1][j]/(ds[i]-ds[j])
    F=R.to_field();tt=F.gens[0]
    S=[[F(S0[i][j])*tt**(j-i)if i>=j else F.zero for j in range(5)]for i in range(5)]
    IF=identity(F.zero,F.one);U=plus(S,scale(IF,-F.one));inv=IF;power=IF
    for k in range(1,5):power=multiply(power,U,F.zero);inv=plus(inv,scale(power,F((-1)**k)))
    assert equal(multiply(S,inv,F.zero),IF)and equal(multiply(inv,S,F.zero),IF)
    diag=matrix_zero(F.zero)
    for i in range(5):diag[i][i]=tt*F(ds[i])
    assert equal(multiply([[F(v)for v in row]for row in Nt],S,F.zero),multiply(S,diag,F.zero))
    diagonal=[]
    for k in range(1,5):
        transformed=multiply(multiply(inv,[[F(v)for v in row]for row in W[k-1]],F.zero),S,F.zero)
        expected=matrix_zero(F.zero);mean=sum((v**k for v in ds),K.zero)/5
        for i in range(5):expected[i][i]=tt**k*F(ds[i]**k-mean)
        assert equal(transformed,expected);diagonal.append(expected)
    # Canonical toric representative in the eigenframe; Euler quotient is checked below.
    diagraw=[F.zero]*150
    for i,j,c in [(0,1,F(mu*alpha)*tt*tt),(0,2,F(mu*beta)*tt),(0,3,F(mu*rho)),(1,2,F(mu))]:
        diagraw=[a+c*b for a,b in zip(diagraw,wedge_coefficients(diagonal[i],diagonal[j]))]
    canonical_log=[v/(F(mu)*tt**5)for v in diagraw]
    assert all(all(mon==(0,)for mon in v.numer)and all(mon==(0,)for mon in v.denom)for v in canonical_log)
    canonical_log=[v.numer.get((0,),K.zero)/v.denom.get((0,),K.zero)for v in canonical_log]
    # Reconstruct the three-entry T38 log representative and its canonical projection.
    log=[K.zero]*150
    for pair,value in [((1,2),p0),((1,3),p0/2),((2,4),-p0/2)]:log[basis.index((pair,pair))]=value
    div=[[K.zero]*5 for _ in range(5)]
    for ((i,j),(u,v)),value in zip(basis,log):
        for a,b,c in [(i,j,value),(j,i,-value)]:
            if a==u:div[b][v]+=c
            if a==v:div[b][u]+=c
    vertical=wedge_coefficients(identity(K.zero,K.one),div)
    projected=[v-c/5 for v,c in zip(log,vertical)]
    assert canonical_log==projected
    def encpoly(v):return [[int(mon[0]),enc(c)]for mon,c in sorted(v.items())]
    payload=dict(field_polynomial_theta=[305,0,-50,0,1],field_basis='1,theta,theta^2,theta^3',
        theta_real_embedding='sqrt(25-8*sqrt(5)); identities valid in all four embeddings',
        field_irreducibility='5-Eisenstein',rho=enc(rho),z=enc(z),x=enc(x),y=enc(y),
        alpha=enc(alpha),beta=enc(beta),mu=enc(mu),s=enc(ss),a=enc(aa),
        diagonal_eigenvalues=[enc(v)for v in ds],N_subdiagonal_links=[enc(v)for v in links],
        unit_certificates=units,gamma={f'{i},{j}':enc(v)for (i,j),v in gamma.items()},
        period_vector_Pfaffian_dual=[enc(gamma[3,4]),enc(-gamma[2,4]),enc(gamma[2,3]),enc(gamma[1,4]),enc(-gamma[1,3]),enc(gamma[1,2])],
        T38_period_basis=[['0'],['1/2'],['0'],['0'],['-1/2'],['1']],
        eigenframe_coefficient_S0=[[enc(v)for v in row]for row in S0],
        eigenframe_rule='S_ij=S0_ij*t^(j-i), i>=j; x=S*y',
        eigenframe_inverse_rule='I-U+U^2-U^3+U^4, U=S-I',
        canonical150_arc=[encpoly(v)for v in raw],canonical150_limit=[enc(v)for v in target],
        canonical_T38_log150=[enc(v)for v in canonical_log],
        full_arc_formula='mu*(alpha*t^2*W1^W2+beta*t*W1^W3+rho*W1^W4+W2^W3)',
        W_definition='Wk=(N+t*diag(d))^k-tr((N+t*diag(d))^k)*I/5')
    with (out/'certificate.json.gz').open('wb')as file:
        with gzip.GzipFile(fileobj=file,mode='wb',filename='',mtime=0)as f:f.write((json.dumps(payload,separators=(',',':'))+'\n').encode())
    note=ROOT/f'research/notes/{STEM}.md'
    inputs=[Path(__file__),note,ROOT/'code/q_t01_membership_20261103_v63.py',ROOT/'code/q_weighted_finite_model_validation_20261103_v63.py',
            ROOT/'reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T38.json',ROOT/'reports/Q_ALL_TORIC_MPS_COMPONENTS_20261103_v63.json']
    result=dict(id=STEM,all_passed=True,exact_field_degree=4,field_irreducibility_5_Eisenstein=True,
        exact_rho_equation='rho^2-rho-1=0',exact_z_equation='2*rho*(z^2+1)+7*z=0',
        all5_eigenvalues_distinct=True,all10_difference_inverses_and_norms_saved=True,
        nonzero_period_scale_gamma12=True,nonzero_period_Pfaffian=True,
        exact_identity_label_T38_period_equations=True,
        actual_polynomial_full150_canonical_arc_regular_at_t0=True,
        full150_t0_equal_Qa8over9=True,all25_physical_divergence_and350_strict_Jacobi_zero=True,
        exact_actual_Laurent_eigenframe_and_inverse=True,all4_tracefree_powers_diagonalized=True,
        full150_eigenframe_tensor_equals_canonical_T38_log_with_scale='mu*t^5',
        T38_source_scope='actual three-entry toric boundary of the frozen full T38 smoothing source',
        smoothing_edges_enabled=False,positive_axis_membership_or_P2_actual_branch_claimed=False,
        all_Qa_abstract_membership_already_known_not_recounted=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
        evidence_sha256={'certificate.json.gz':sha(out/'certificate.json.gz')},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('EXACT T38 QUARTIC FRAME / ACTUAL150 QA8OVER9 ARC PASS',result['seconds'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

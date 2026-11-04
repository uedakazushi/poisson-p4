#!/usr/bin/env python3
"""Support/denominator preflight only: no degree-26 polynomial or rank job."""
import argparse,gzip,hashlib,json,math,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_V_ZERO_DEGREE26_SUPPORT_PREFLIGHT_20261104_v63'
ZERO=(0,0,0,0)
def read(p):
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else p.open()) as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(a,b):
    c=dict(a)
    for m,x in b.items():
        c[m]=c.get(m,Q(0))+x
        if not c[m]:del c[m]
    return c
def scale(a,c):return {m:x*c for m,x in a.items() if x*c}
def mul(a,b):
    c={}
    for m,x in a.items():
        for n,y in b.items():
            k=tuple(u+v for u,v in zip(m,n));c[k]=c.get(k,Q(0))+x*y
    return {m:x for m,x in c.items() if x}
def poly(e,vs):return {tuple(m):Q(str(c)) for m,c in s.Poly(e,*vs).terms() if c}
def rows(a):return [[[list(m),str(c)]for m,c in sorted(p.items())]for p in a]
def lcm(xs):return math.lcm(*xs)
def l1(p):return sum(map(abs,p.values()),Q(0))
@lru_cache(None)
def support(d,w):
    """Outer support from homogeneous ambient degree and both Cartan weights."""
    w1,w2=w;ans=[]
    for b in range(d+1):
        if (d+w2-3*b)%2 or (d-w2+b)%2:continue
        np=(d+w2-3*b)//2;nm=(d-w2+b)//2
        if min(np,nm)<0 or b%2!=w1%2:continue
        for a in range(d-b+1):
            if (2*a+3*b+w1-2*w2)%4:continue
            for c in range(np+1):
                for e in range(np-c+1):ans.append((a,b,c,e))
    return ans
def ceilq(x):return (x.numerator+x.denominator-1)//x.denominator
def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    ip=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    sp=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    vp=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63/certificate.json'
    dp=ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63/definitions.json'
    gates=[ROOT/'reports/Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63_REVIEW/report.json',ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed']for p in gates)
    info,spin,vzero,defs=map(read,[ip,sp,vp,dp]);weights=[tuple(w)for w in info['Cartan_weights']]
    xi,r,ss,q,B,T,R,H,L=s.symbols('xi0 r0 s0 q0 B0 T R H L')
    raw=list(map(s.sympify,vzero['beta_bar0_16']));beta=[];torus=[]
    for i,(z,(w1,w2))in enumerate(zip(raw,weights)):
        assert s.expand(2*xi*s.diff(z,xi)-(1-w1-w2)*z)==0
        transformed=s.cancel(z.subs({xi:1,B:T*T,r:R**4})*T**(w1+2*w2)*R**(-w1+2*w2))
        if transformed:
            p=s.Poly(transformed,R,T,ss,q);terms={}
            for (er,et,es,eq),c in p.terms():
                assert er>=3*et and (er-3*et)%2==0
                terms[((er-3*et)//2,et,es,eq)]=Q(str(c))
            beta.append(terms)
        else:beta.append({})
        torus.append(str(transformed))
    expected=[-45*ss/16,-30*L,2400,0,0,H*(5*q+ss+1)/240,ss/54000,0,-32*H/3,-s.Rational(32,675),0,0,0,0,H*(ss+1)/120,-64*H/3]
    assert beta==[poly(z,[H,L,ss,q])for z in expected]
    assert sum(bool(p)for p in beta)==9 and sum(map(len,beta))==12
    for i,p in enumerate(beta):assert set(p)<=set(support(1,weights[i]))
    Omega=s.Matrix(info['symplectic_form16']);Oi=Omega.inv();rho=[s.Matrix(z)for z in info['all_ten_representation_matrices']]
    J=s.diag(s.Matrix([[0,1],[-1,0]]),s.Matrix([[0,1],[-1,0]]));gs=[]
    for i in range(4):
        for j in range(i,4):
            S=s.zeros(4);S[i,j]=S[j,i]=1;gs.append(J*S)
    for h,ww in [(s.diag(1,-1,0,0),[w[0]for w in weights]),(s.diag(0,0,1,-1),[w[1]for w in weights])]:
        sym=J.inv()*h;co=[sym[i,j]for i in range(4)for j in range(i,4)]
        assert sum((c*A for c,A in zip(co,rho)),s.zeros(16))==s.diag(*ww)
    mon=[tuple(z)for z in spin['quadratic_monomials']];mi={m:i for i,m in enumerate(mon)}
    coeff=[s.Rational(z)for key in ['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']for row in spin[key]for z in row]
    CD=s.Matrix(32,136,coeff);quadratic=[];derivatives=[];cdweights=[]
    spinweights=[(1,0),(-1,0),(0,1),(0,-1)]
    def tensor_weight_check(A,n,ww):
        for row in range(n*n):
            i,j=divmod(row,n);w=tuple(x-y for x,y in zip(ww[i],ww[j]))
            assert all(not A[row,k]or tuple(weights[u][h]+weights[v][h]for h in range(2))==w for k,(u,v)in enumerate(mon))
    tensor_weight_check(CD[:16,:],4,spinweights);tensor_weight_check(CD[16:,:],4,spinweights)
    for row in range(32):
        value={};grad=[{}for _ in range(16)]
        for k,(i,j)in enumerate(mon):
            c=Q(str(CD[row,k]))
            if c:
                value=add(value,scale(mul(beta[i],beta[j]),c))
                grad[i]=add(grad[i],scale(beta[j],c));grad[j]=add(grad[j],scale(beta[i],c))
        rr,cc=divmod(row%16,4);w=tuple(x-y for x,y in zip(spinweights[rr],spinweights[cc]));cdweights.append(w)
        assert set(value)<=set(support(2,w))
        check={}
        for k,p in enumerate(grad):
            assert set(p)<=set(support(1,tuple(w[h]-weights[k][h]for h in range(2))))
            check=add(check,mul(beta[k],p))
        assert check==scale(value,2)
        quadratic.append(value);derivatives.extend(grad)
    # Full global 2theta and 2H16 from the accepted quartic, only degree TWO.
    Gi=s.Matrix([[s.trace(A*B)for B in gs]for A in gs]).inv();mu=[]
    for A in rho:
        S=Omega*A;assert S==S.T;mu.append([S[i,j]/(2 if i==j else 1)for i,j in mon])
    rflat=s.Matrix([[A[i,j]for A in rho]for i in range(16)for j in range(16)])
    Theta=2*rflat*Gi*s.Matrix(mu);hess=s.zeros(256,136)
    for indices,c in spin['quartic_trace_moment_square']:
        counts=Counter(indices);c=s.Rational(c)/2
        for i,ei in counts.items():
            for j,ej in counts.items():
                m=ei*(ej-(i==j))
                if m:
                    rest=list(indices);rest.remove(i);rest.remove(j);hess[16*i+j,mi[tuple(sorted(rest))]]+=c*m
    Hess=s.zeros(256,136)
    for i in range(16):
        for k in range(16):
            if Oi[i,k]:
                for j in range(16):Hess[16*i+j,:]+=2*Oi[i,k]*hess[16*k+j,:]
    assert all(c.is_Integer for c in Theta)and all(c.is_Integer for c in Hess)
    tensor_weight_check(Theta,16,weights);tensor_weight_check(Hess,16,weights)
    hvalues=[]
    for A in [Theta,Hess]:
        for row in range(256):
            p={}
            for k,(i,j)in enumerate(mon):
                if A[row,k]:p=add(p,scale(mul(beta[i],beta[j]),Q(str(A[row,k]))/2))
            i,j=divmod(row,16);w=tuple(x-y for x,y in zip(weights[i],weights[j]))
            assert set(p)<=set(support(2,w));hvalues.append(p)
    assert len(defs['Inv26_basis'])==45
    # Only combinatorial support envelopes are formed after quadratic entrance.
    table={str(d):len(support(d,(0,0)))for d in [4,8,12,14,18,26]}
    assert table['26']==2874
    gradient_bounds={str(d):[len(support(d-1,tuple(-x for x in w)))for w in weights]for d in[8,12]}
    matrixbounds=[]
    for n,ww,degree_range in[(4,spinweights,range(2,13,2)),(16,weights,range(2,13,2))]:
        for d in degree_range:
            matrixbounds.append(dict(size=n,degree=d,max_entry_support=max(len(support(d,tuple(x-y for x,y in zip(w,v))))for w in ww for v in ww)))
    Kbeta=lcm([c.denominator for p in beta for c in p.values()]);Kcd=lcm([Q(str(c)).denominator for c in CD]);Ko=lcm([Q(str(c)).denominator for c in Oi])
    M=max(map(l1,quadratic));DM=max(map(l1,derivatives));MH=max(map(l1,hvalues));OO=sum(abs(Q(str(c)))for c in Oi)
    trace=lambda ell:4*(4*M)**ell
    gradient=lambda ell:4*ell*(4*M)**(ell-1)*(4*DM)
    F=4*(12*M)**2;T4=trace(4);T6=trace(6);G4=gradient(4);G6=gradient(6)
    N=trace(6)+trace(2)*trace(4)/2+24*M**4*trace(2)
    H6=16*(16*MH)**6;Q14=OO*G4*G4;Q18=OO*G4*G6
    inv12_norms=[F*T4]*5+[N,T6,T6,H6]
    qden=Ko*Kcd**8*Kbeta**14;pden=Ko*Kcd**10*Kbeta**18
    inv12_dens=[Kcd**6*Kbeta**12]*5+[2*Kcd**6*Kbeta**12,Kcd**6*Kbeta**12,Kcd**6*Kbeta**12,2**6*Kbeta**12]
    norms=[];dens=[]
    for z in defs['Inv26_basis']:
        if z['kind']=='Inv14_times_Inv12':norms.append(Q14*inv12_norms[z['Inv12_basis_index']]);dens.append(qden*inv12_dens[z['Inv12_basis_index']])
        else:
            assert z['kind']=='primitive18_times_Inv8';norms.append(Q18*T4);dens.append(pden*Kcd**4*Kbeta**8)
    bits=[ceilq(n*d).bit_length()for n,d in zip(norms,dens)]
    counts=[2874]*45;binary=sum(n*((b+7)//8+16)for n,b in zip(counts,bits));decimal=sum(n*(math.ceil(b*math.log10(2))+32)for n,b in zip(counts,bits))
    # All exponent data is small. No invariant degree >=4 has been expanded.
    cert=dict(variables_polynomial=['H','L','s0','q0'],coordinates_original=['xi0','r0','s0','q0','B0'],
        gauges=dict(T_squared='B0',R_fourth='r0',positive_pullback='diag(T/R,R/T,T^2*R^2,T^-2*R^-2)',H='R^2',L='R^3*T'),
        beta16=[str(z)for z in expected],beta16_sparse=rows(beta),beta_nonzero9_terms12=True,
        xi_homogeneity='I_degree_d(beta0)=xi0^(d/2)*I_degree_d(beta0|xi0=1)',
        finite_descent='H^a L^b: 2a+3b=0 mod4; r0^((2a+3b)/4)*B0^(b/2)',
        support_proof=dict(ambient_degree='d',ambient_Cartan_weight=['w1','w2'],L_exponent='b',plus_count='(d+w2-3b)/2',minus_count='(d-w2+b)/2',constraints=['counts nonnegative integers','b=w1 mod2','0<=a<=d-b','2a+3b=-w1+2w2 mod4','s_degree+q_degree<=plus_count']),
        invariant_support_counts=table,Inv26_support_by_L_exponent={str(b):sum(m[1]==b for m in support(26,(0,0)))for b in range(9)},
        gradient_support_bounds=gradient_bounds,matrix_entry_support_bounds=matrixbounds,
        CD_quadratic32=rows(quadratic),CD_ambient_derivative512=rows(derivatives),CD_actual_support_counts=list(map(len,quadratic)),CD_derivative_max_terms=max(map(len,derivatives)),
        theta_H_global_quadratic_weights_checked=True,theta_H_source_quadratic_max_terms=max(map(len,hvalues)),
        denominator_constants=dict(beta=Kbeta,CD=Kcd,Omega_inverse=Ko),Inv26_constant_denominators=list(map(str,dens)),
        Inv26_coefficient_numerator_bit_bounds=bits,Inv26_support_bound_per_column=2874,
        rigorous_uncompressed_binary_coefficient_exponent_bytes_upper_bound=binary,conservative_decimal_json_bytes_upper_bound=decimal,
        source_variable_denominators_after_descent=False,restored_scalar_power=13,
        actual_full_45_basis=defs['Inv26_basis'],
        no_degree26_polynomial_expansion_or_rank_or_new_point=True,
        proposed_next_job='One native FLINT exact low trace/ambient-gradient restriction, then ONE GF32003 coefficient rank of selected45; requires separate parent allocation. Stream outputs, hard60s/2GB/32MiB, intermediate150k/output100k caps; no retry.',
        warning='Support and coefficient-height upper bounds do not predict cancellation or certify runtime. No Inv26 restriction or new equation is asserted.')
    data=json.dumps(cert,indent=2)+'\n';assert len(data.encode())<32*1024*1024;(out/'certificate.json').write_text(data)
    paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',ip,sp,vp,dp,*gates]
    report=dict(id=ID,all_passed=True,only_support_denominator_parity_preflight=True,all16_scalar_and_two_Cartan_identities=True,ambient16_CD_derivatives_before_substitution=True,
        actual_low_quadratic_and_Euler_checks=True,no_degree26_polynomial_expansion_rank_or_source_point=True,
        inv26_support_bound2874=True,variable_denominators_removed_by_Sp4_and_finite_descent=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=report['seconds'],support2874=True,coefficient_bits_max=max(bits),binary_bytes=binary,decimal_bytes=decimal)))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)

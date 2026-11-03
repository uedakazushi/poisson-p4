#!/usr/bin/env python3
"""Exact no-proportionality certificate for Q0,Q3 on B=W=0, a*P*J!=0."""
from pathlib import Path
import hashlib
import json
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'work/nonsquare-v62/gcd-proportional-certificate-v58'

def main():
    started=time.monotonic();OUT.mkdir(parents=True,exist_ok=True)
    source=ROOT/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json'
    assert hashlib.sha256(source.read_bytes()).hexdigest()=='80f7c46b1bd4ffef4a78f492ad4f7a2bb40b2ffa27cddcd23d5b5d660db4a7e0'
    data=json.loads(source.read_text());a,e,P=s.symbols('a e P');x=s.symbols('x0:6')
    r,U,V,T,R=s.symbols('r U V T R');ys=(r,U,V,T,R);quads=[]
    for i,row in enumerate(data['cleared_qfour_sparse']):
        q=0
        for mon,terms in row:
            q+=sum(s.Rational(co)*a**(ex[0]-ex[2]+[0,2,2,1][i])*e**ex[1]*P**ex[2]
                   for ex,co in terms)*s.prod(v**k for v,k in zip(x,mon[:6]))
        quads.append(q)
    J=5*P-6*e-4;Rc=15*P**2-6*P*e+6*P+80*a;Tc=-5*P+6*e-2
    Sc=21*P*e+16*P+140*a;Vc=9*P*e+4*P+60*a;Kc=85*P**2-12*P*e+58*P+600*a
    B=30*a*Rc*x[0]-9*P*Sc*x[1]+90*P*a*Tc*x[2]-18*P*Vc*x[3]-2*a*Kc*x[4]
    A=(B+10*a*P*J*x[4])/3;Q2=s.cancel((quads[2]-e*quads[1])/(a*a*J))
    assert s.expand(quads[1].subs(x[5],0)-A*B/6075)==0
    assert s.expand(Q2.subs(x[5],0)-40*x[0]*B/6075)==0
    assert s.expand(3*Sc-7*Vc-20*P)==0
    inverse={x[0]:r,x[1]:V,x[2]:R/a-U-6*(e+1)*V/(5*a)-(9*e+4)*T/(10*a),
             x[3]:(T-V)/2,x[4]:5*U+5*r-4*R/a+24*(e+1)*V/(5*a)+2*(9*e+4)*T/(5*a),x[5]:0}
    assert s.cancel(4*quads[0].subs(inverse,simultaneous=True)-U*R+5*V*T)==0
    Q=s.Poly(s.cancel(a*a*quads[3].subs(inverse,simultaneous=True)),*ys)
    for c in Q.coeffs():s.Poly(c,a,e,P,domain=s.QQ)
    H=s.hessian(Q.as_expr(),ys)/2
    BN=s.Poly(s.cancel(B.subs(inverse,simultaneous=True)),*ys)
    br,bU,bV,bT,bR=[BN.coeff_monomial(y)for y in ys];bs=[bU,bV,bT,bR]
    DB=20*P**2+3*P*e+20*P+180*a;assert s.expand(br+20*a*DB)==0
    cleared_U2=H[0,0]*bU*bU-2*br*H[0,1]*bU+br*br*H[1,1]
    assert s.expand(cleared_U2+a**5*br**2*J/12)==0
    # Hence on DB!=0 the U^2 coefficient of unscaled Q3|B is -a^3*J/12.
    Lp=20*P+3*e+20;ad=-P*Lp/180
    assert s.cancel(H[0,0].subs(a,ad))==0
    h=[s.Poly(s.cancel(H[0,i+1].subs(a,ad)),P,e,domain=s.QQ)for i in range(4)]
    b=[s.Poly(s.cancel(v.subs(a,ad)),P,e,domain=s.QQ)for v in bs]
    residual=[];powers=[];indices=[]
    for i in range(4):
        for j in range(i+1,4):
            poly=h[i]*b[j]-h[j]*b[i];vals=[]
            for unit in (P,Lp,J):
                unitpoly=s.Poly(unit,P,e,domain=s.QQ);n=0
                while not poly.is_zero:
                    quotient,rem=poly.div(unitpoly)
                    if rem:break
                    poly=quotient;n+=1
                vals.append(n)
            residual.append(poly.as_expr());powers.append(vals);indices.append([i,j])
    delta=-5*P+6*e-5
    f1=-100*P**2+120*P*e-25*P-36*e-60
    f2=-50*P**2+60*P*e-25*P-18*e-15
    f3=-60*P*e-40*P+72*e**2+6*e-29
    scales=[s.Integer(12754584000000),s.Integer(8503056000000),-s.Integer(141717600000),
            s.Integer(93312000000),s.Integer(15552000000)]
    G=[s.expand(c*q)for c,q in zip(scales,residual)]
    for got,want in zip(G[:3],[delta*f1,delta*f2,delta*f3]):assert s.expand(got-want)==0
    p0=s.Rational(6,5);e0=s.Rational(13,6)
    k2=s.cancel((f2-f2.subs(P,p0))/(P-p0))
    k3=s.cancel((f3-f3.subs(P,p0))/(P-p0))
    kE=s.cancel((f3.subs(P,p0)-118)/(e-e0))
    alpha=-k3/25+kE*k2/1350;beta=-kE/54-2*alpha
    assert s.expand(alpha*f1+beta*f2+f3-118)==0
    assert s.expand(alpha*G[0]+beta*G[1]+G[2]-118*delta)==0
    q1=25*P**2+35*P-6;q2=675*P**2+2445*P+1594
    assert s.expand(G[3].subs(e,(5*P+5)/6)-q1)==0
    assert s.expand(G[4].subs(e,(5*P+5)/6)-q2)==0
    d4=s.cancel((G[3]-q1)/delta);d5=s.cancel((G[4]-q2)/delta)
    s.Poly(d4,P,e,domain=s.QQ);s.Poly(d5,P,e,domain=s.QQ)
    u,v,gcd=s.gcdex(q1,q2,P);assert gcd==1
    weight=s.expand(u*d4+v*d5)
    cof=[-weight*alpha/118,-weight*beta/118,-weight/118,u,v]
    assert s.expand(sum(c*q for c,q in zip(cof,G))-1)==0
    for c in cof:s.Poly(c,P,e,domain=s.QQ)
    encode=lambda q,vs:[[list(ex),str(c)]for ex,c in s.Poly(q,*vs,domain=s.QQ).terms()]
    cert=dict(normalized_minors=[encode(g,(P,e))for g in G],Bezout_coefficients=[encode(c,(P,e))for c in cof],
              corresponding_cross_minor_indices=indices[:5],known_unit_powers_P_L_J=powers[:5],normalization_QQ_scalars=list(map(str,scales)),
              global_H_a2Q3=[[encode(H[i,j],(a,e,P))for j in range(5)]for i in range(5)],
              global_B_normal_coefficients=[encode(v,(a,e,P))for v in (br,*bs)])
    (OUT/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    result=dict(all_exact_checks_passed=True,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                DB=str(DB),known_open='a*P*J!=0',DB_nonzero_Q3_U2_coefficient='-a^3*J/12',
                DB_zero_a_elimination=str(ad),DB_zero_extra_known_unit='20P+3e+20=-180a/P',
                exact_two_parameter_Bezout_identity=True,no_proportionality_uniform_on_known_open=True,
                combined_with_rank2_packet_gcd_gate_failure_zero_dimensional=True,
                five_quadric_rank_gate_not_proved=True,pair_gate_not_proved=True,raw_rank24_or_original_product41_not_proved=True,
                seconds=time.monotonic()-started)
    (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)

if __name__=='__main__':main()

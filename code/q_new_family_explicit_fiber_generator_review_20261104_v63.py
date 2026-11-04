#!/usr/bin/env python3
"""Independent sparse rational arithmetic and actual-field/deck/descent proof audit."""
import argparse
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.fields import field

ROOT=Path(__file__).resolve().parents[1]
BASE='Q_NEW_FAMILY_EXPLICIT_FIBER_GENERATOR_20261104_v63';STEM=BASE+'_REVIEW'
DOUBLE='Q_NEW_FAMILY_FIBER_DOUBLE_COVER_20261104_v63';DREVIEW=DOUBLE+'_REVIEW'
A=256
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def pathkey(p):
    try:return str(p.relative_to(ROOT))
    except ValueError:return str(p)
def determinant(M,one):
    answer=one*0
    for p in it.permutations(range(len(M))):
        term=one*((-1)**sum(p[i]>p[j] for i in range(len(M)) for j in range(i+1,len(M))))
        for i in range(len(M)):term*=M[i][p[i]]
        answer+=term
    return answer
def Jformula(u,K,Iv,w):
    L=Iv*A*A*K*((K+1)*w-(K-1)*u)+8*u*(K-1)*((K+1)*u-(K-1)*w)
    return L*L/(K*u*u*(K-1)**2*(u*u-w*w))


def main(out,dout):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    folder=ROOT/'reports'/BASE;primary=read(folder/'report.json');cert=read(folder/'certificate.json')
    dfolder=ROOT/'reports'/DOUBLE;dprimary=read(dfolder/'report.json');dcert=read(dfolder/'certificate.json')
    assert primary['all_passed'] and dprimary['all_passed']
    gate_names=['Q_NEW_FAMILY_GENERIC_STABILIZER_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_SIGMA_FINITE_EQUIVALENCE_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_FULL_SOURCE_DEGREE_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_FULL_TENSOR_FIBER_PREFLIGHT_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW',
        'Q_NEW_FAMILY_RATIONAL_GROUP_INVERSE_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63_REVIEW']
    gates=[ROOT/'reports'/n/'report.json' for n in gate_names]
    gd=[read(p) for p in gates];assert all(x['all_passed'] for x in gd)
    assert gd[2]['original_full_source_generic_degree_exact8']
    assert gd[2]['descended_full_source_fiber_generic_degree_exact4_QQ_and_C']
    assert gd[3]['same_K_B15_field_inclusion_into_K_z_reviewed']
    assert gd[6]['new_F0_actual_quotient_index1'] and gd[6]['actual_full_boundary_Y15_to_B15_birational_section']
    # Independently use rational PolyElement fields instead of expression cancel.
    F,u,v,c,X,W,B=field('u,v,c,X,W,B',QQ)
    epath=ROOT/'reports/Q_NEW_FAMILY_GENERIC_STABILIZER_20261104_v63/certificate.json'
    spath=ROOT/'reports/Q_NEW_FAMILY_SIGMA_FINITE_EQUIVALENCE_20261104_v63/certificate.json'
    E=[[F.from_expr(s.sympify(x)) for x in row] for row in read(epath)['E']]
    branch=next(b for b in read(spath)['branches'] if b['sign']==1)
    M=[[F.from_expr(s.sympify(x)) for x in row] for row in branch['homogeneous_matrix']]
    pt=[F.one,X,F.zero,F.zero,W]
    image=lambda N:[sum((a*b for a,b in zip(row,pt)),F.zero) for row in N]
    rr=lambda vv,xx,ww:A+vv*(384*xx+3072*ww-1)
    ss=lambda vv,xx,ww:A+vv*(384*xx-3072*ww+1)
    ep=image(E);mp=image(M);vp=-32*u-v
    assert E[1][3]==E[4][3]==M[1][3]==M[4][3]==0
    assert rr(v,ep[1],ep[4])==-rr(v,X,W) and ss(v,ep[1],ep[4])==-ss(v,X,W)
    assert rr(vp,mp[1],mp[4])==-rr(v,X,W) and ss(vp,mp[1],mp[4])==ss(v,X,W)
    detE=determinant([row[1:] for row in E[1:]],F.one)
    detM=determinant([row[1:] for row in M[1:]],F.one)
    assert detE==1
    kap=F.from_expr(s.sympify(branch['Poisson_scalar']))
    cp=c/kap;Bp=kap**2*detM*B
    assert cp==-c*v*v/(vp*vp) and Bp==B*v**3/vp**3
    geo=read(ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json')
    substitution={s.Symbol('n0'):s.Symbol('u'),s.Symbol('n1'):1,s.Symbol('n2'):s.Symbol('v')}
    modular=[[F.from_expr(s.sympify(z).subs(substitution)) for z in row] for row in geo['homogeneous_modular_N']]
    assert not any(modular[1]) and not any(modular[4]) and modular[2][0]==4*u/3

    C,u,v,c,B,r,sf=field('u,v,c,B,r,sf',QQ)
    def constants(uu,vv,cc,bb,rf,ss):
        h=1-A*A/(rf*rf)
        return [uu,h/(1-A*A/(ss*ss)),
            (vv-8*uu*(rf-ss)**2/(rf*ss))/(rf*ss*h),
            cc*vv*vv/(rf*ss*h*h),bb*vv**3/(rf*rf*ss*ss*h**3)]
    invariants=constants(u,v,c,B,r,sf);vp=-32*u-v
    assert constants(u,v,c,B,-r,-sf)==invariants
    assert constants(u,vp,-c*v*v/vp**2,B*v**3/vp**3,-r,sf)==invariants
    tau=constants(-u,-v,-c,B,sf,r);K0,I0,C0,B0=invariants[1:]
    assert tau==[-u,1/K0,-K0*I0,-K0**2*C0,-K0**3*B0]

    T,K,t=field('K,t',QQ)
    def r_of(t0):
        D=K*t0*t0-K+1;N=K*t0-K+1
        return A*(D-2*N)/D
    def s_of(t0):
        D=K*t0*t0-K+1;N=K*t0-K+1
        return A*(D-2*N)/(D-2*t0*N)
    def m1(t0):return (K*t0-K+1)/(K*(t0-1))
    def m2(t0):return (K-1)/(K*t0)
    rt,st=r_of(t),s_of(t)
    assert m1(m1(t))==m2(m2(t))==t and m1(m2(t))==m2(m1(t))
    assert m1(t)!=t and m2(t)!=t and m1(t)!=m2(t)
    assert r_of(m1(t))==-rt and s_of(m1(t))==-st
    assert r_of(m2(t))==-rt and s_of(m2(t))==st
    assert (1-A*A/(rt*rt))/(1-A*A/(st*st))==K
    q=(rt*rt-A*A)/(st*st-A*A)
    assert q==1+(K-1)*rt*rt/(A*A)
    phi=(q-1)/(q+1);tpoly=T.ring.gens[1]
    assert max(phi.numer.degree(tpoly),phi.denom.degree(tpoly))==4
    assert phi.numer.gcd(phi.denom).degree(tpoly)==0
    assert phi== (K-1)*rt*rt/(2*A*A+(K-1)*rt*rt)

    Fw,u,K,Iv,w=field('u,K,Iv,w_fiber',QQ)
    R=2*A*A*w/((K-1)*(u-w));S=A*A*K*R/((K-1)*R+A*A)
    assert u*(K-1)*R/(2*A*A+(K-1)*R)==w
    jw=Jformula(u,K,Iv,w)
    assert jw==(Iv*S*(R-A*A)+8*u*(R+S))**2/(u*u*R*S)
    assert jw==Fw.from_expr(s.sympify(cert['J_in_w']))
    assert Jformula(-u,1/K,-K*Iv,w)==jw
    wp=Fw.ring.gens[3]
    assert max(jw.numer.degree(wp),jw.denom.degree(wp))==2
    assert jw.numer.gcd(jw.denom).degree(wp)==0
    # Compact coprimeness proof: L(+u) and L(-u) are generic nonzero.
    L0=u*(K-1)*(-Iv*A*A*K+8*u*(K+1))
    L1=Iv*A*A*K*(K+1)-8*u*(K-1)**2
    assert L0+u*L1==2*u*(A*A*Iv*K+8*u*(K-1))!=0
    assert L0-u*L1==2*u*K*(-A*A*Iv*K+8*u*(K-1))!=0
    Ft,ut,Kt,Ivt,tt=field('u,K,Iv,t',QQ)
    D=Kt*tt*tt-Kt+1;N=Kt*tt-Kt+1
    rf=A*(D-2*N)/D;sf=A*(D-2*N)/(D-2*tt*N)
    wt=ut*(Kt-1)*rf*rf/(2*A*A+(Kt-1)*rf*rf)
    vt=Ivt*rf*sf*(1-A*A/(rf*rf))+8*ut*(rf-sf)**2/(rf*sf)
    jt=(vt/ut+16)**2
    assert jt==Jformula(ut,Kt,Ivt,wt)==Ft.from_expr(s.sympify(cert['J_in_t']))
    tp=Ft.ring.gens[3]
    assert max(jt.numer.degree(tp),jt.denom.degree(tp))==8
    assert jt.numer.gcd(jt.denom).degree(tp)==0

    # New v0 coordinates, explicitly kept distinct from the old t2 coordinates.
    Fd,u,r0,s0,q0,w,J=field('u,r0,s0,q0,w_fiber,J',QQ)
    Ksub=(8*s0+u*q0)/(8*s0-u*q0);Ivsub=u/(512*(8*s0+u*q0))
    h=s0*s0-25*r0*q0*q0;AA=200*r0*q0*(s0-1);BB=s0-25*r0*q0*q0
    DD=r0*q0*q0*h;UU=1600*r0
    expected=256*(AA+BB*w)**2/(25*DD*(UU-w*w))
    difference=Jformula(u,Ksub,Ivsub,w)-expected
    def square_reduce(poly):
        out={}
        for ex,co in poly.items():
            e=list(ex);power=e[0]//2;e[0]%=2;e[1]+=power
            e=tuple(e);out[e]=out.get(e,QQ.zero)+co*1600**power
        return Fd.ring.from_dict({e:c for e,c in out.items() if c})
    assert square_reduce(difference.numer)==0 and square_reduce(difference.denom)!=0
    assert expected==Fd.from_expr(s.sympify(dcert['J_in_w']))
    jcrit=256-256/(25*r0*q0*q0)
    assert 256*(AA*AA-BB*BB*UU)/(25*DD*UU)==jcrit
    ca=-J*DD-256*BB*BB/25;cb=-512*AA*BB/25;cc=J*DD*UU-256*AA*AA/25
    disc=cb*cb-4*ca*cc
    assert disc==(80*DD)**2*r0*J*(J-jcrit)
    assert disc==Fd.from_expr(s.sympify(dcert['discriminant']))
    def j_at(z):return 256*(AA+BB*z)**2/(25*DD*(UU-z*z))
    w0=-AA/BB;wc=-BB*UU/AA;boundary=200*r0*q0/s0
    assert j_at(w0)==0 and j_at(wc)==jcrit and j_at(boundary)==256
    assert square_reduce((u*(Ksub-1)/(Ksub+1)-boundary).numer)==0
    assert jcrit!=0 and w0!=wc

    checked={}
    for p in [folder/'report.json',dfolder/'report.json',*gates]:
        d=read(p)
        for key in ['input_sha256','checked_dependency_sha256']:
            for name,want in d.get(key,{}).items():
                target=ROOT/name;assert sha(target)==want,str(target);checked[name]=want
        for name,want in d.get('evidence_sha256',{}).items():
            target=p.parent/name;assert sha(target)==want,str(target);checked[str(target.relative_to(ROOT))]=want
    out.mkdir(parents=True,exist_ok=False);dout.mkdir(parents=True,exist_ok=False)
    proof={'independent_sparse_rational_matrix_E_sigma_and_all5_constant_identities':True,
        'faithful_Mobius_V4_and_r_squared_fixed_field_degree':4,
        'tau_fixed_w_and_inverse_R_formula_verified':True,
        'same_actual_field_equality_uses_full_source_degree8_not_universal_m_degree':True,
        'actual_generator':'w=u*(K-1)*r_f^2/(2*256^2+(K-1)*r_f^2)',
        'J_compact_formula':'L(w)^2/[K*u^2*(K-1)^2*(u^2-w^2)]',
        'L0':'u*(K-1)*(-Iv*256^2*K+8*u*(K+1))',
        'L1':'Iv*256^2*K*(K+1)-8*u*(K-1)^2',
        'J_actual_degree':2,'J_source_t_degree':8,
        'generic_function_fields_QQ_and_C_only':True,'new_full150_expansion_point_rank_or_invariant_job':False}
    ep=out/'independent-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    dproof={'new_v0_coordinate_substitution_verified_by_sparse_quadratic_reduction':True,
        'J_discriminant':'(80*D)^2*r0*J*(J-Jcrit)','Jcrit':'256-256/(25*r0*q0^2)',
        'actual_field':'K_B(w)=K_B(J)(sqrt(r0*J*(J-Jcrit)))',
        'branch_values':['0','Jcrit'],'branch_w_points':['-A/B','-1600*r0*B/A'],
        'v0_section_w':'200*r0*q0/s0','v0_section_J':256,
        'all_special_fibers_global_morphism_or_NEW_equation_claim':False}
    dp=dout/'independent-double-cover-proof.json';dp.write_text(json.dumps(dproof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),ROOT/'research/notes'/(DREVIEW+'.md'),
        folder/'report.json',folder/'certificate.json',dfolder/'report.json',dfolder/'certificate.json',*gates]
    for p in paths:checked[str(p.relative_to(ROOT))]=sha(p)
    common={'all_passed':True,'input_sha256':checked,'seconds':time.monotonic()-start,
        'no_author_rerun_full150_newpoint_rank_high_invariant_or_source26':True,'frozen':True}
    report=dict(common,id=STEM,
        full_actual_V4_deck_lifts_constants_and_unique_beta_frame_reviewed=True,
        faithful_V4_fixed_field_and_tau_descent_verified=True,
        actual_same_beta_field_equals_KB_w_by_accepted_full_source_degree8=True,
        actual_J_map_degree2_and_source_t_degree8=True,
        old_parameter_V4_not_confused_with_source_deck_V4=True,
        generic_QQ_and_C_field_scope_no_specialfiber_or_NEW_equation_claim=True,
        evidence_sha256={'independent-proof.json':sha(ep)})
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    dinputs=dict(checked);dinputs[pathkey(out/'report.json')]=sha(out/'report.json')
    dreport=dict(common,id=DREVIEW,input_sha256=dinputs,
        independent_new_v0_J_formula_and_quadratic_discriminant_verified=True,
        actual_quadratic_field_sqrt_rJ_JminusJcrit_and_branch_values_verified=True,
        generic_two_rational_branch_points_and_v0_J256_boundary_verified=True,
        new_v0_coordinates_distinguished_from_old_t2=True,
        no_all_specialfiber_classification_NEW_equation_or_spectral_degree_claim=True,
        evidence_sha256={'independent-double-cover-proof.json':sha(dp)})
    (dout/'report.json').write_text(json.dumps(dreport,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':report['seconds']}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--double-output',type=Path,required=True)
    args=p.parse_args();main(args.output.resolve(),args.double_output.resolve())

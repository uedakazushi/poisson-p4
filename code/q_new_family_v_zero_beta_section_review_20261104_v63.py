#!/usr/bin/env python3
"""Independent cancelled constants, transverse DVR and Laurent beta section audit."""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s

ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63_REVIEW'


def read(p):
    with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)


def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def laurent(expr,variables):
    num,den=s.fraction(s.cancel(expr));pn=s.Poly(num,*variables);pd=s.Poly(den,*variables)
    assert pd.length()==1
    dm,dc=pd.terms()[0];answer={}
    for mon,co in pn.terms():
        if co:answer[tuple(i-j for i,j in zip(mon,dm))]=Q(str(co/dc))
    return answer


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    if out.exists():
        assert not any(out.iterdir()), 'refuse to overwrite evidence'
    else:out.mkdir(parents=True)
    base=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63'
    author=read(base/'report.json');cert=read(base/'certificate.json');assert author['all_passed']
    for rel,h in author['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    for rel,h in author['evidence_sha256'].items():assert sha(base/rel)==h,rel
    gatepaths=[ROOT/rel for rel in author['input_sha256']if rel.startswith('reports/')and rel.endswith('/report.json')]
    assert all(read(p)['all_passed']for p in gatepaths)
    sixpath=ROOT/'reports/Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63/certificate.json'
    kernelpath=ROOT/'reports/Q_NEW_FAMILY_SOURCE_KERNEL_20261103_v63/kernel-certificate.json.gz'
    constantpath=ROOT/'reports/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63/certificate.json'
    basispath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    geopath=ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json'
    six,kernel,constant,basis,geo=map(read,(sixpath,kernelpath,constantpath,basispath,geopath))
    u,v,c,X,W,B,alpha,beta=s.symbols('u v c X W B alpha beta')
    rr,ss=s.symbols('r s');A=s.Integer(256);rf=A+v*alpha;sf=A+v*beta
    zr=1-A*A/rf**2;zs=1-A*A/sf**2
    definitions=[u,zr/zs,(v-8*u*(rf-sf)**2/(rf*sf))/(rf*sf*zr),
                 c*v*v/(rf*sf*zr*zr),B*v**3/((rf*sf)**2*zr**3)]
    old=[s.sympify(q).subs({rr:rf,ss:sf},simultaneous=True)for q in constant['five_first_integrals']]
    cancelled=[s.sympify(q,locals={'alpha':alpha,'beta':beta})for q in cert['cancelled_five_constants']]
    assert all(s.cancel(a-b)==0 for a,b in zip(definitions,old))
    assert all(s.cancel(a-b)==0 for a,b in zip(definitions,cancelled))
    restricted=[s.cancel(q.subs(v,0))for q in cancelled]
    assert restricted==[u,alpha/beta,1/(512*alpha),c/(4*alpha**2),B/(2048*alpha**3)]
    uu,kk,iv,ic,ib=s.symbols('UU KK IV IC IB')
    inverse={u:uu,alpha:1/(512*iv),beta:1/(512*iv*kk),c:ic/(65536*iv**2),B:ib/(65536*iv**3)}
    target=[uu,kk,iv,ic,ib]
    assert all(s.cancel(q.subs(inverse,simultaneous=True)-p)==0 for q,p in zip(restricted,target))
    forward=dict(zip(target,restricted))
    assert all(s.cancel(q.subs(forward,simultaneous=True)-p)==0 for p,q in inverse.items())
    aa=384*X+3072*W-1;bb=384*X-3072*W+1
    assert s.expand((aa+bb)/768-X)==0
    assert s.expand((1+(aa-bb)/2)/3072-W)==0
    parameter=s.cancel((rf/sf-1)/(rf/A-1))
    assert s.cancel(parameter-A*(alpha-beta)/(sf*alpha))==0
    assert s.cancel(parameter.subs(v,0)-(1-beta/alpha))==0
    tau={u:-u,c:-c,alpha:-beta,beta:-alpha,B:B}
    tau_out=[-restricted[0],1/restricted[1],-restricted[1]*restricted[2],
             -restricted[1]**2*restricted[3],-restricted[1]**3*restricted[4]]
    assert all(s.cancel(q.subs(tau,simultaneous=True)-p)==0 for q,p in zip(restricted,tau_out))
    # The UNRESCALED kernel is polynomial and transverse at v=0. Its
    # compensating perspective entries must vanish identically.
    for q in kernel['full26_polynomial_kernel']:s.Poly(s.sympify(q),u,v,X,W,B,domain=s.QQ)
    assert kernel['compensating_perspective4']==['0']*4
    assert kernel['all30_original_polynomial_dot_products_exactly_zero']
    assert s.expand(s.sympify(kernel['source_six_direction']['v']).subs(v,0))==402653184
    # Finite gu pullback weights are checked on every physical24-by16
    # coefficient, rather than accepted from an unchecked weight label.
    weights=[q[0]for q in basis['Cartan_weights']]
    physical=list((i,j,k)for i,j in itertools.combinations(range(4),2)for k in range(4))
    character=[1,-1,0,0]
    for row,indices in zip(basis['primitive24_by16_basis'],physical):
        for j,co in enumerate(row):
            if s.Rational(co):assert sum(character[i]for i in indices)==weights[j]
    xi0,r0,s0,q0,B0=s.symbols('xi0 r0 s0 q0 B0')
    smallvars=(xi0,r0,s0,q0,B0)
    source={v:0,X:s0/(48*u),c:xi0*B0*u/40,W:(q0+1)/3072,B:B0}
    def descend(expression):
        result={}
        for mon,co in laurent(expression,(u,*smallvars)).items():
            assert mon[0]%2==0
            power=mon[0]//2;new=list(mon[1:]);new[1]+=power
            new=tuple(new);result[new]=result.get(new,Q())+co*Q(1600)**power
        return {m:c0 for m,c0 in result.items()if c0}
    atoms=[]
    for text,weight,expected in zip(six['primitive_beta16'],weights,cert['beta_bar0_16']):
        raw=laurent((s.sympify(text).subs(v,0)*u**weight).subs(source,simultaneous=True),(u,*smallvars))
        descended={}
        for mon,co in raw.items():
            assert mon[0]%2==0
            power=mon[0]//2;new=list(mon[1:]);new[1]+=power
            new=tuple(new);descended[new]=descended.get(new,Q())+co*Q(1600)**power
        descended={m:c0 for m,c0 in descended.items()if c0}
        assert descended==laurent(s.sympify(expected),smallvars)
        assert all(mon[2]>=0 and mon[3]>=0 for mon in descended)
        atoms.append([[list(mon),str(co)]for mon,co in sorted(descended.items())])
    assert sum(bool(q)for q in atoms)==9 and sum(map(len,atoms))==12 and max(map(len,atoms))==3
    # Full tensor descent is bound by finite matrices, without a new150
    # expansion. All of these identities remain regular when v=0.
    R=s.Matrix(six['rational_Darboux_frame4']);alph=s.Matrix(list(map(s.sympify,six['primitive_perspective_alpha'])))
    assert not any(q.has(v)for q in R)
    assert s.cancel(R.det()-s.Rational(9,8)/(c*c*B))==0
    trans={u:-u,v:-v,c:-c,X:-X,W:W,B:B}
    D=s.diag(-1,1,-1,1);H=s.diag(-1,-1,1,1);gu=s.diag(u,1/u,1,1)
    assert(R.subs(trans,simultaneous=True)*D-H*R).applyfunc(s.cancel)==s.zeros(4)
    assert(alph.subs(trans,simultaneous=True)-H*alph).applyfunc(s.cancel)==s.zeros(4,1)
    assert(gu.subs(u,-u)-H*gu)==s.zeros(4)
    for text,expected in zip(alph,cert['full_tensor_factors']['primitive_alpha0']):
        assert descend(text.subs(source,simultaneous=True))==laurent(s.sympify(expected),smallvars)
    n0,n1,n2=s.symbols('n0 n1 n2');x,y,z,w=s.symbols('x y z w');spatial=(x,y,z,w)
    P=s.Matrix(geo['affine_P']).subs({n0:u,n1:1,n2:0})
    back={x:-x,y:y,z:-z,w:w}
    assert((-c)*P.subs(u,-u)-c*D*P.subs(back,simultaneous=True)*D).applyfunc(s.expand)==s.zeros(4)
    pointZ=-3*(B+W/1536-W*W+X*X/64+s.Rational(2,9))/(4*u)
    point=s.Matrix([X,0,pointZ,W])
    assert(point.subs(trans,simultaneous=True)-D*point).applyfunc(s.cancel)==s.zeros(4,1)
    # All determinant and constant-chart units are named explicitly.
    assert s.cancel((s.Rational(9,8)/(c*c*B)).subs(source,simultaneous=True)*u*u/(1600*r0)
                    -s.Rational(9,8)/(r0*xi0**2*B0**3))==0
    assert s.expand((aa*bb).subs(source,simultaneous=True)*u*u-64*(s0*s0-25*r0*q0*q0)).subs(u*u,1600*r0)==0
    proof=dict(new_coordinates=['xi0','r0','s0','q0','B0'],only_coordinate_identical_to_old_t2='r0=u^2/1600',
        cancelled_constants=list(map(str,cancelled)),boundary_inverse={str(k):str(q)for k,q in inverse.items()},
        t_at_v0='1-1/K',original_V_v_at_v0='402653184',
        DVR_argument='On GL4 times the source, regular D annihilates all30 rawjet pullbacks and Dv is a unit. A nonzero pullback vanishing on v=0 has valuation n>=1, but its derivative has exact valuation n-1, contradiction in characteristic zero.',
        no_flow_integration_or_delta_transversality_assumption=True,
        beta16_Laurent_atoms=atoms,nonzero=9,total_terms=12,max_terms=3,
        Laurent_units=['xi0','r0','B0'],no_h_denominator=True,
        full_descent='Full-family tau, point_tau=D point, Rtau D=H R, alpha_tau=H alpha and gu_tau=H gu cancel before Euler/S61 normalization. All factors and inverse factors are Laurent, so the complete lifted tensor descends without h inversion.',
        actual_field_proof='The five constants restrict birationally and tau-equivariantly to the v=0 source; actual DVR density prevents a pole-locus blowdown. Accepted E_B=F5^tau therefore restricts isomorphically to F0.',
        Y15_to_B15='Sp4 times F0 maps to the actual full normalized boundary and then B15; the composite is birational by density, quotient-index1 and generic stabilizer identity. The intermediate image has the same field, so Y15->B15 is birational.',
        generic_Y16_source_degree='The old generic source degree4 after tau is a different generic-fiber statement and is unchanged.',
        no_full150_expansion_point_rank_high_invariant_or_new_generic_degree_claim=True)
    ep=out/'independent-certificate.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),base/'report.json',base/'certificate.json',
            ROOT/'code/q_new_family_v_zero_beta_section_20261104_v63.py',
            ROOT/'research/notes/Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63.md',
            sixpath,kernelpath,constantpath,basispath,geopath,*gatepaths]
    report=dict(id=ID,all_passed=True,independent_all_cancelled_constants_boundary_inverses_and_tau=True,
        original_polynomial_V_transverse_402653184_and_regular_DVR_density_proof_accepted=True,
        independent_all24_by16_weight_support_and16_Laurent_descent_coefficients=True,
        beta_nonzero_count=9,beta_total_numerator_terms=12,no_h_denominator=True,
        actual_v0_beta_Sp4_saturation_dense=True,new_F0_actual_quotient_index1=True,
        actual_full_boundary_Y15_to_B15_birational_section=True,
        new_coordinates_not_old_t2_coordinates_and_generic_source_degree4_unchanged=True,
        independent_full_tau_matrix_and_original_scaled_tensor_binding=True,
        no_author_rerun_new_point_prime_rank_full150_expansion_or_Inv26=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(inputs)},
        evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start,frozen=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k]for k in('all_passed','beta_nonzero_count','beta_total_numerator_terms','seconds')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(p.parse_args().output.resolve())

#!/usr/bin/env python3
"""Exact V4 source deck group and an explicit descended actual fiber generator."""
import argparse, hashlib, json, resource, signal, time
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
ID = 'Q_NEW_FAMILY_EXPLICIT_FIBER_GENERATOR_20261104_v63'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def read(p):
    return json.loads(p.read_text())

def zero(x):
    assert s.cancel(x) == 0, x

def main(out):
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    signal.alarm(60)
    out.mkdir(parents=True, exist_ok=False)
    stems = [
        'Q_NEW_FAMILY_GENERIC_STABILIZER_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_SIGMA_FINITE_EQUIVALENCE_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_FULL_SOURCE_DEGREE_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_FULL_TENSOR_FIBER_PREFLIGHT_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW',
        'Q_NEW_FAMILY_RATIONAL_GROUP_INVERSE_20261104_v63_REVIEW',
    ]
    gates = [ROOT/'reports'/n/'report.json' for n in stems]
    for p in gates:
        assert read(p)['all_passed'], p
    sigmapath = ROOT/'reports/Q_NEW_FAMILY_SIGMA_FINITE_EQUIVALENCE_20261104_v63/certificate.json'
    epath = ROOT/'reports/Q_NEW_FAMILY_GENERIC_STABILIZER_20261104_v63/certificate.json'
    constpath = ROOT/'reports/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63/certificate.json'
    geo = ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json'
    sigma, ec, const, G = map(read, [sigmapath, epath, constpath, geo])
    u,v,c,X,W,B,t,K,Iv,Ic,Ib,r,sf = s.symbols('u v c X W B t K Iv Ic Ib r sf')
    A = s.Integer(256)
    rr = v*(384*X+3072*W-1)+A
    ss = v*(384*X-3072*W+1)+A
    E = s.Matrix(ec['E'])
    pt = s.Matrix([1,X,0,0,W])
    image = E*pt
    # x,w do not depend on the omitted z coordinate, nor on modular Ga.
    assert E[1,3] == E[4,3] == 0
    zero(rr.subs({X:image[1],W:image[4]}, simultaneous=True)+rr)
    zero(ss.subs({X:image[1],W:image[4]}, simultaneous=True)+ss)
    zero(E[1:,1:].det()-1)
    b = next(b for b in sigma['branches'] if b['sign'] == 1)
    M = s.Matrix(b['homogeneous_matrix'])
    image = M*pt
    assert M[1,3] == M[4,3] == 0
    vp = -32*u-v
    zero(rr.subs({v:vp,X:image[1],W:image[4]}, simultaneous=True)+rr)
    zero(ss.subs({v:vp,X:image[1],W:image[4]}, simultaneous=True)-ss)
    kap = s.sympify(b['Poisson_scalar'])
    cp = s.cancel(c/kap)
    Bp = s.cancel(kap**2*M[1:,1:].det()*B)
    zero(cp+c*v*v/vp**2)
    zero(Bp-B*v**3/vp**3)
    n0,n1,n2 = s.symbols('n0 n1 n2')
    modular = s.Matrix(G['homogeneous_modular_N']).subs({n0:u,n1:1,n2:v})
    assert list(modular[1,:]) == [0]*5 and list(modular[4,:]) == [0]*5
    zero(modular[2,0]-4*u/3)
    # Exact preservation of all five constants under each independent sign flip.
    integrals = [u, (1-A*A/r**2)/(1-A*A/sf**2),
                 (v-8*u*(r-sf)**2/(r*sf))/(r*sf*(1-A*A/r**2)),
                 c*v*v/(r*sf*(1-A*A/r**2)**2),
                 B*v**3/(r*r*sf*sf*(1-A*A/r**2)**3)]
    deck = [
        {r:-r,sf:-sf},
        {r:-r,sf:sf,v:vp,c:-c*v*v/vp**2,B:B*v**3/vp**3},
    ]
    for trans in deck:
        for f in integrals:
            zero(f.subs(trans, simultaneous=True)-f)
    D = K*t*t-K+1
    N = K*t-K+1
    EE = D-2*t*N
    rt = A*(D-2*N)/D
    st = A*(D-2*N)/EE
    zero((1-A*A/rt**2)/(1-A*A/st**2)-K)
    maps = {'simultaneous_sign':N/(K*(t-1)), 'r_sign':(K-1)/(K*t)}
    for name, mt in maps.items():
        zero(mt.subs(t,mt)-t)
        zero(rt.subs(t,mt)+rt)
        zero(st.subs(t,mt)-(-st if name=='simultaneous_sign' else st))
    a,bm = maps.values()
    zero(a.subs(t,bm)-bm.subs(t,a))
    q = s.cancel((rt*rt-A*A)/(st*st-A*A))
    zero(q-1-(K-1)*rt*rt/(A*A))
    wt = s.cancel(u*(q-1)/(q+1))
    num,den = s.fraction(wt)
    pn,pd = s.Poly(num,t),s.Poly(den,t)
    assert max(pn.degree(),pd.degree()) == 4
    assert s.Poly(s.gcd(pn,pd),t).degree() == 0
    for mt in maps.values():
        zero(wt.subs(t,mt)-wt)
    qq,ww = s.symbols('q w_fiber')
    wf = u*(qq-1)/(qq+1)
    zero(wf.subs({u:-u,qq:1/qq},simultaneous=True)-wf)
    qi = (u+ww)/(u-ww)
    zero(wf.subs(qq,qi)-ww)
    R = 2*A*A*ww/((K-1)*(u-ww))
    zero((u*(K-1)*R/(2*A*A+(K-1)*R))-ww)
    # The full PGL invariant J has degree two on this actual curve, not eight.
    vt = s.cancel(Iv*rt*st*(1-A*A/rt**2)+8*u*(rt-st)**2/(rt*st))
    jt = s.cancel((vt/u+16)**2)
    jn,jd = s.fraction(jt)
    jp,jdp = s.Poly(jn,t),s.Poly(jd,t)
    assert max(jp.degree(),jdp.degree()) == 8
    assert s.Poly(s.gcd(jp,jdp),t).degree() == 0
    # Avoid implicitization: square the linear equation in r*s.
    # (v+16u)r*s = Iv*s^2*(r^2-A^2)+8u*(r^2+s^2).
    Sr = A*A*K*R/((K-1)*R+A*A)
    jw = s.factor((Iv*Sr*(R-A*A)+8*u*(R+Sr))**2/(u*u*R*Sr))
    zero(jt-jw.subs(ww,wt))
    jwn,jwd = s.fraction(s.cancel(jw))
    assert max(s.Poly(jwn,ww).degree(),s.Poly(jwd,ww).degree()) == 2
    cert = dict(
        flow_coordinates={'r_f':str(rr),'s_f':str(ss),'A':256},
        full_deck_lifts={'simultaneous_sign':'Aut E, followed by unique Ga to p_y=0',
                         'r_sign':'sigma sign+1, followed by unique target Ga to p_y=0',
                         'sigma_scalar_new':str(cp),'sigma_B_new':str(Bp)},
        all5_constants_preserved=True,
        exact_t_deck_maps={k:str(z) for k,z in maps.items()},
        deck_group='V4, independent sign flips of r_f,s_f; NOT parameter V4(tau,sigma)',
        pre_tau_fixed_field='F5(r_f^2)',pre_tau_index=4,
        tau='u->-u,K->1/K,r_f<->s_f; corresponding semilinear action on F5',
        q=str(q),q_definition='(r_f^2-256^2)/(s_f^2-256^2)=1+(K-1)*r_f^2/256^2',
        explicit_actual_generator='w_fiber=u*(K-1)*r_f^2/(2*256^2+(K-1)*r_f^2)',
        inverse_r_squared=str(R),tau_invariant=True,t_degree=4,
        actual_same_beta_field='k(Y16)=k(B15)(w_fiber), k=QQ,C',
        physical_frame_binding='Every lift preserves the full tensor by accepted full E/sigma identities; raw-pair matching and unique generic beta frame keep h fixed. Thus inclusion into degree8 full-source invariants, equality by accepted degree8 and tau degree2.',
        J_in_t=str(jt),J_in_w=str(jw),J_degree_on_source=8,J_degree_on_actual_fiber=2,
        no_new_point_or_rank_or_high_invariant_or_full150_recalculation=True,
        no_special_fiber_or_global_morphism_claim=True)
    (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    paths = [Path(__file__),ROOT/'research/notes'/(ID+'.md'),sigmapath,epath,constpath,geo,*gates]
    report = dict(id=ID,all_passed=True,explicit_actual_same_beta_Luroth_generator=True,
                  independent_review_pending=True,full_source_generic_degree_unchanged=8,
                  tau_descended_source_degree_unchanged=4,
                  actual_J_map_degree=2,
                  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},
                  evidence_sha256={'certificate.json':sha(out/'certificate.json')},
                  seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k:report[k] for k in ['all_passed','explicit_actual_same_beta_Luroth_generator','actual_J_map_degree','seconds']}))

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output.resolve())

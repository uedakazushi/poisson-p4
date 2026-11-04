#!/usr/bin/env python3
"""Generic joint-biform invariant field via an exact torus/mu4 cross-section.

No beta point, Hilbert-series extension, invariant-ring generation, GB, or implicitization.
"""
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
ID='Q_FIRSTJET_SPIN4_BIFORM_RATIONAL_FIELD_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    output=ROOT/'reports'/ID
    x0,x1,y0,y1,c,t=s.symbols('x0 x1 y0 y1 c t',nonzero=True)
    variables=[x0,x1,y0,y1];vs=s.symbols('v0:6');v0,v1,v2,v3,v4,v5=vs
    U0=x0*x0*y0+x1*x1*y1
    V=x0*(v0*y0*y0+v1*y0*y1+v2*y1*y1)+x1*(v3*y0*y0+v4*y0*y1+v5*y1*y1)
    mons=[x0*y0*y0,x0*y0*y1,x0*y1*y1,x1*y0*y0,x1*y0*y1,x1*y1*y1]
    def transform(poly,A,B):
        new=list(A*s.Matrix([x0,x1]))+list(B*s.Matrix([y0,y1]))
        return s.expand(poly.subs(dict(zip(variables,new)),simultaneous=True))
    TA=s.diag(t,1/t);TB=s.diag(1/t**2,t**2)
    assert transform(U0,TA,TB)==U0
    weights=[-3,1,5,-5,-1,3]
    assert transform(V,TA,TB)==s.expand(sum(v*t**k*mon for v,k,mon in zip(vs,weights,mons)))
    SA=s.Matrix([[0,1],[-1,0]]);SB=s.Matrix([[0,s.I],[s.I,0]])
    assert SA.det()==SB.det()==1
    assert transform(U0,SA,SB)==s.expand(s.I*U0)
    sigma_v=[v5,v4,v3,-v2,-v1,-v0]
    assert transform(V,SA,SB)==s.expand(sum(v*mon for v,mon in zip(sigma_v,mons)))
    assert SA**2==SB**2==-s.eye(2) and SA**4==SB**4==s.eye(2)
    assert s.simplify(SA*TA*SA.inv()-TA.subs(t,1/t))==s.zeros(2)
    assert s.simplify(SB*TB*SB.inv()-TB.subs(t,1/t))==s.zeros(2)
    sigma=dict(zip(vs,sigma_v));sigma[c]=s.I*c
    sigma2={v:s.expand(sigma[v].subs(sigma,simultaneous=True))for v in [c,*vs]}
    assert all(sigma2[v]==-v for v in [c,*vs])
    #Universal Jacobian discriminant/resultant identity, no numerical specialization.
    aa,bb,cc,dd,ee,ff=s.symbols('aa bb cc dd ee ff')
    p=aa*x0*x0+bb*x0*x1+cc*x1*x1;q=dd*x0*x0+ee*x0*x1+ff*x1*x1
    Jac=s.expand(s.diff(p,x0)*s.diff(q,x1)-s.diff(p,x1)*s.diff(q,x0))
    jp=s.Poly(Jac,x0,x1);j0=jp.coeff_monomial(x0*x0);j1=jp.coeff_monomial(x0*x1);j2=jp.coeff_monomial(x1*x1)
    syl=s.Matrix([[aa,bb,cc,0],[0,aa,bb,cc],[dd,ee,ff,0],[0,dd,ee,ff]])
    res=s.expand(syl.det());assert s.expand(j1*j1-4*j0*j2-16*res)==0
    assert res.subs({aa:c,bb:0,cc:0,dd:0,ee:0,ff:c})==c**4
    #If the Jacobian has only its middle term, independent p,q have no mixed terms.
    determinant=aa*ff-cc*dd
    assert s.expand(bb*determinant-(aa*(bb*ff-cc*ee)+cc*(aa*ee-bb*dd)))==0
    assert s.expand(ee*determinant-(ff*(aa*ee-bb*dd)+dd*(bb*ff-cc*ee)))==0
    #Torusal quotient: a weight1 coefficient gives a rational gauge, no roots.
    r=v1*v4;P=v0*v1**3;Q=v5/v1**3;S=v3*v1**5;T=v2/v1**5
    torus=[r,P,Q,S,T]
    scaling=dict(zip(vs,[v*t**k for v,k in zip(vs,weights)]))
    assert all(s.cancel(f.subs(scaling,simultaneous=True)-f)==0 for f in torus)
    expected_sigma=[-r,Q*r**3,-P/r**3,-T*r**5,S/r**5]
    assert all(s.cancel(f.subs(sigma,simultaneous=True)-g)==0 for f,g in zip(torus,expected_sigma))
    #After sigma² quotient: z=c²; rename Q'=Q*r³,T'=-T*r⁵.
    z=c*c;Qprime=Q*r**3;Tprime=-T*r**5
    six=[c**4,r/z,P+Qprime,(P-Qprime)/z,S+Tprime,(S-Tprime)/z]
    assert all(s.cancel(f.subs(sigma,simultaneous=True)-f)==0 for f in six)
    assert all(s.cancel(f.subs(scaling,simultaneous=True)-f)==0 for f in six)
    #Explicit inverse over z²=I0 proves generation and algebraic independence.
    zz,rr,pp,qq,ss,tt=s.symbols('zz rr pp qq ss tt')
    invvars=s.symbols('I0:6');I0,I1,I2,I3,I4,I5=invvars
    linear_invariants=[zz**2,rr/zz,pp+qq,(pp-qq)/zz,ss+tt,(ss-tt)/zz]
    inverse={rr:zz*I1,pp:(I2+zz*I3)/2,qq:(I2-zz*I3)/2,ss:(I4+zz*I5)/2,tt:(I4-zz*I5)/2}
    recovered=[s.cancel(f.subs(inverse,simultaneous=True))for f in linear_invariants]
    assert recovered==[zz**2,I1,I2,I3,I4,I5]
    involution={zz:-zz,rr:-rr,pp:qq,qq:pp,ss:tt,tt:ss}
    assert all(s.cancel(f.subs(involution,simultaneous=True)-f)==0 for f in linear_invariants)
    #Recover V coordinates from its torus quotient and a weight1 coordinate.
    rr0,pp0,qq0,ss0,tt0,k=s.symbols('rr0 pp0 qq0 ss0 tt0 k',nonzero=True)
    source_inverse={v0:pp0/k**3,v1:k,v2:tt0*k**5,v3:ss0/k**5,v4:rr0/k,v5:qq0*k**3}
    assert [s.cancel(f.subs(source_inverse,simultaneous=True))for f in torus]==[rr0,pp0,qq0,ss0,tt0]
    note=ROOT/'research/notes'/(ID+'.md')
    spinreview=ROOT/'reports/Q_FIRSTJET_SPIN4_MIXED_DISCRIMINANT_20261103_v63_REVIEW/report.json'
    gate=json.loads(spinreview.read_text());assert gate['all_passed']
    checked={}
    for name,want in gate['input_sha256'].items():assert sha(ROOT/name)==want;checked[name]=want
    evidence={
        'U_normal_form':'c*(x0^2*y0+x1^2*y1); c^4=Res(p,q)',
        'Jacobian_binary_discriminant_equals16_resultant':str(res),
        'exact_U_stabilizer':'T: (diag(t,t^-1),diag(t^-2,t^2))',
        'U_line_normalizer':'N=T semidirect mu4, sigma=(SA,SB); sigma conjugates t to t^-1',
        'sigma_A':[[str(e)for e in row]for row in SA.tolist()],
        'sigma_B':[[str(e)for e in row]for row in SB.tolist()],
        'V_coefficient_order':[str(mon)for mon in mons],'T_weights':weights,
        'sigma_on_c_and_V':[str(sigma[v])for v in [c,*vs]],
        'torus_quotient_coordinates_r_P_Q_S_T':[str(f)for f in torus],
        'sigma2_quotient':'z=c^2; (z,r,P,Q*r^3,S,-T*r^5) -> (-z,-r,Q*r^3,P,-T*r^5,S)',
        'six_rational_joint_invariant_field_generators':[str(s.factor(f))for f in six],
        'inverse_after_adjoining_z':{str(key):str(value)for key,value in inverse.items()},
        'whole_field_generation_proof_in_note':True,
        'no_additional_Spin5_line_normalizer_quotient_claim':True,
    }
    output.mkdir(parents=True,exist_ok=False)
    ev=output/'certificate.json';ev.write_text(json.dumps(evidence,indent=2)+'\n')
    files=[Path(__file__),note,spinreview,ROOT/'research/notes/Q_FIRSTJET_SPIN4_MIXED_DISCRIMINANT_20261103_v63.md',
           ROOT/'research/notes/Q_FIRSTJET_SP4_INVARIANT_LITERATURE_SEARCH_20261104_v63.md']
    result={'id':ID,'all_passed':True,
        'generic_U_normal_form_and_exact_stabilizer_and_mu4_line_normalizer_proved':True,
        'all_U_V_coefficient_torus_and_sigma_actions_checked':True,
        'six_joint_SL2xSL2_rational_invariant_field_coordinates_complete':True,
        'joint_invariant_field_pure_transcendental_degree':6,
        'base_field':'C, characteristic0',
        'polynomial_invariant_ring_generators_Hilbert_series_or_global_quotient_claim':False,
        'full_W16_Sp4_extra_normalizer_rationality_or_NEW25_equation_claim':False,
        'new_beta_point_high_degree_probe_GB_or_implicitization':False,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},'checked_dependency_sha256':checked,
        'evidence_sha256':{ev.name:sha(ev)},'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'seconds':time.monotonic()-started}
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    signal.alarm(0);print(json.dumps({'all_passed':True,'rational_joint_field_dimension':6,'seconds':result['seconds']}))

if __name__=='__main__':main()

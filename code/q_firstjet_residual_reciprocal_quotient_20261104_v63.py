#!/usr/bin/env python3
"""Generic reciprocal involution and quartic quotient; no beta evaluation."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_FIRSTJET_RESIDUAL_RECIPROCAL_QUOTIENT_20261104_v63'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out.mkdir(parents=True,exist_ok=False)
    # One small universal coefficient computation, with no actual beta/source point.
    c,d,f,u,v,R,T,k,b,x,y,r,t,z,D=s.symbols('c d f u v R T k b x y r t z D')
    mm=u-v
    rr=-v*R+f*u/R;tt=-c*T+u*v/T
    zz=(v*R+f*u/R)*(c*T+u*v/T)
    L1=x*rr+y*tt-k
    aa=s.Symbol('a')
    L2=b+c*R*T+f*u*u/(R*T)-c*f*T/R-v*v*R/T-2*c*d*T+2*aa*f*u/R-2*aa*v*R+2*d*u*v/T
    Ri=-f*u/(v*R);Ti=-u*v/(c*T)
    assert s.cancel(L1.subs({R:Ri,T:Ti},simultaneous=True)-L1)==0
    assert s.cancel(L2.subs({R:Ri,T:Ti},simultaneous=True)-L2)==0
    assert s.cancel(zz*zz-(rr*rr+4*f*u*v)*(tt*tt+4*c*u*v))==0
    assert s.cancel(2*u*v*(L2-b-2*aa*rr-2*d*tt)-(u+v)*rr*tt-mm*zz)==0
    Psi=r*r*t*t+(u+v)*r*t*D+u*v*D*D-mm*mm*(c*r*r+f*t*t+4*c*f*u*v)
    W=2*u*v*D+(u+v)*r*t
    assert s.expand(W*W-mm*mm*(r*r+4*f*u*v)*(t*t+4*c*u*v)-4*u*v*Psi)==0
    Psi1=s.Poly(s.expand(Psi.subs({t:(k-x*r)/y,D:b+2*aa*r+2*d*(k-x*r)/y})),r)
    assert Psi1.degree()==4 and s.cancel(Psi1.LC()-x*x/(y*y))==0
    # Actual generic coefficients, cleared with named x,y denominators only.
    A,*gens=ring('x,y,a,k1,k2,r,m',s.QQ)
    xx,yy,a,k1,k2,r,m=gens
    hh=4*xx*yy;S=xx*xx+yy*yy;delta=xx*xx-yy*yy
    ell=(a+xx)**2-yy*yy
    dnum=-a*xx-S;cnum=m-ell;fnum=yy*yy*m-xx*xx*ell
    B2=-4*a*a*xx*xx-2*a*xx**3-6*a*xx*yy*yy+3*xx**4+yy**4
    U=k2-B2-(xx-yy)**2*m;V=k2-B2-(xx+yy)**2*m
    abar=(hh*a*dnum*(a*yy*yy+dnum*xx)+yy*(a*yy*yy+dnum*xx)*(U+V)
          +hh*yy*yy*cnum*dnum+hh*a*xx*fnum)
    Kbar=hh*yy*yy*k1/2-abar
    Tbar=Kbar-hh*yy*yy*xx*r
    Bbar=hh*(a*a+cnum)*(dnum*dnum+fnum)+hh*yy*yy*(m*m-delta*delta)+2*a*dnum*yy*(U+V)
    Dbar=yy*yy*Bbar+2*a*hh*yy**4*r+2*dnum*Tbar
    plane=(hh*hh*yy*yy*r*r*Tbar*Tbar+hh*yy*(U+V)*r*Tbar*Dbar+U*V*Dbar*Dbar
           -m*m*(hh**4*yy**8*cnum*r*r+hh*hh*fnum*Tbar*Tbar+4*hh*hh*yy**6*cnum*fnum*U*V))
    assert plane
    # divide only a literal common x/y monomial and rational content; no factor/GB.
    exmin=[min(ex[j] for ex in plane) for j in range(len(gens))]
    divisor=A.one
    for j in (0,1):divisor*=gens[j]**exmin[j]
    primitive=plane.exquo(divisor)
    content=s.gcd_list([s.Rational(q.numerator,q.denominator) for q in primitive.values()])
    primitive=primitive.exquo(A(content))
    degree_r=max(ex[5] for ex in primitive);degree_m=max(ex[6] for ex in primitive)
    total=max(ex[5]+ex[6] for ex in primitive)
    assert degree_r==4
    # The nonzero B0^2+16a^2fuv norm has m^4 leading coefficient4.
    norm_numerator=Bbar*Bbar+16*a*a*fnum*U*V*yy*yy
    assert max(ex[6] for ex in norm_numerator)==4
    leading_norm={ex[:6]:co for ex,co in norm_numerator.items() if ex[6]==4}
    expect=4*hh*hh*yy**4
    assert A.from_dict({ex+(0,):co for ex,co in leading_norm.items()})==expect
    # Store a compact polynomial table grouped by r,m, not an invariant expansion.
    grouped={}
    for ex,co in primitive.items():
        key=f'{ex[5]},{ex[6]}'
        grouped.setdefault(key,[]).append([list(ex[:5]),str(co)])
    certificate={
        'universal_involution':{'R':'-f*u/(v*R)','T':'-u*v/(c*T)','m':'m'},
        'quotient_coordinates':{'r':str(rr),'t':str(tt),'z':str(zz)},
        'universal_quartic':str(Psi),'t_recovery':'(K1-x*r)/y',
        'z_recovery':'-(2*u*v*(B0+2*a*r+2*d*t)+(u+v)*r*t)/m',
        'double_cover':'h^2=r^2+4*f*u*v',
        'R_recovery':'(h-r)/(2*v)','T_recovery':'(z/h-t)/(2*c)',
        'actual_plane_clearing_factor':'(4*x*y)^4*y^8',
        'removed_literal_monomial':str(divisor),'removed_rational_content':str(content),
        'actual_primitive_plane_coefficients_r_m':grouped,
        'coefficient_variable_order':['x','y','a','k1','k2'],
        'plane_bidegree_r_m':[degree_r,degree_m],'plane_total_degree_r_m':total,
        'plane_number_of_nonzero_coefficients_r_m':len(grouped),
        'plane_total_monomials_all_variables':len(primitive),
        'generic_h_unit_norm_B0_squared_plus_16a2fuv_m_degree':4,
        'generic_h_unit_norm_m_leading':4,
    }
    (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    paths=[Path(__file__),
           ROOT/'reports/Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63_REVIEW/report.json',
           ROOT/'reports/Q_FIRSTJET_RESIDUAL_OCTIC_PRESENTATION_20261103_v63_REVIEW/report.json']
    for p in paths[1:]:assert json.loads(p.read_text())['all_passed']
    report={'id':STEM,'all_passed':True,
            'reciprocal_involution_preserves_both_original_Laurent_equations':True,
            'universal_quartic_quotient_and_quadratic_cover_identities_verified':True,
            'generic_total_intersection_quartic_times_quadratic_presentation':True,
            'actual_plane_bidegree_r_m':[degree_r,degree_m],
            'actual_plane_total_degree_r_m':total,
            'curve_genus_irreducibility_selected_component_degree_claimed':False,
            'no_beta_point_new_invariant_GB_or_octic_factorization':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
            'evidence_sha256':{'certificate.json':sha(out/'certificate.json')},
            'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','actual_plane_bidegree_r_m','actual_plane_total_degree_r_m','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)

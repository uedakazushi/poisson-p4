#!/usr/bin/env python3
"""Exact generic octic presentation of the existing fixed-base total intersection.

Small universal coefficient identities only: no new beta point, GB, or Inv24 probe.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
ID='Q_FIRSTJET_RESIDUAL_OCTIC_PRESENTATION_20261103_v63'

def read(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main(output):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    output.mkdir(parents=True,exist_ok=False)
    review=ROOT/'reports/Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63_REVIEW/report.json'
    primary=ROOT/'reports/Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63/report.json'
    certpath=primary.parent/'certificate.json.gz'
    checked={}
    for path in [review,primary]:
        report=read(path);assert report['all_passed']
        for key in ['input_sha256','evidence_sha256']:
            for name,want in report.get(key,{}).items():
                p=next((p for p in [ROOT/name,path.parent/name]if p.is_file()and sha(p)==want),None)
                assert p is not None,(path,name);checked[str(p.relative_to(ROOT))]=want
    actual=read(certpath)
    x,y,c,d,f,u,v,m,R,T,K,B=s.symbols('x y c d f u v m R T K B')
    A=R*R-2*d*R-f
    b=-x*v*R*R+f*x*u-K*R
    q=v*R*R+f*u
    a1=-c*y*R;c1=y*u*v*R
    a=s.Symbol('a')
    b2=B*R+2*a*f*u-2*a*v*R*R
    Q1=a1*T*T+b*T+c1
    Q2=c*A*T*T+b2*T-v*v*R*R+f*u*u+2*d*u*v*R
    L=y*R*b2+A*b
    M=y*m*R*q
    P=u*v*L*L-m*q*b*L-c*y*y*m*m*R*R*q*q
    def reduce(expression):return s.expand(expression.subs(m,u-v))
    assert reduce(y*R*Q2+A*Q1-L*T-M)==0
    E=L*T+M
    #Entire ideal identities, rather than testing only the common zero set.
    assert reduce(L*L*Q1-a1*E*E-(b*L-2*a1*M)*E-y*R*P)==0
    assert s.expand(a1*M*M-b*M*L+c1*L*L-y*R*P)==0
    c2=-v*v*R*R+f*u*u+2*d*u*v*R
    assert reduce(y*R*(c*A*M*M-b2*M*L+c2*L*L)+A*y*R*P)==0
    octic=s.Poly(reduce(P),R)
    assert octic.degree()==8
    assert s.expand(octic.LC()-x*x*v**4)==0
    assert s.expand(octic.TC()-x*x*f**4*u**4)==0
    #The obstruction to a common root of A and q is exactly this quadratic resultant.
    resultant=f*f*(u+v)**2+4*d*d*f*u*v
    #A*q reduction has linear remainder 2dvR+f(u+v).
    assert s.expand(q-v*A-2*d*v*R-f*(u+v))==0
    l=2*d*v;constant=f*(u+v)
    assert s.expand(constant*constant+2*d*l*constant-f*l*l-resultant)==0
    #Bind every compact coefficient to the actual two-adjoint fixed-base model.
    ax,ay,aa,k2,am=s.symbols('ax ay aa k2 am')
    delta=ax*ax-ay*ay;S=ax*ax+ay*ay;ell=(aa+ax)**2-ay*ay
    dd=-(aa*ax+S)/ay;cc=am-ell;ff=am-ax*ax*ell/ay**2
    B2=-4*aa*aa*ax*ax-2*aa*ax**3-6*aa*ax*ay*ay+3*ax**4+ay**4
    kap=k2-B2;uu=(kap-(ax-ay)**2*am)/(4*ax*ay);vv=(kap-(ax+ay)**2*am)/(4*ax*ay)
    assert s.cancel(uu-vv-am)==0
    actual_result=s.Poly(s.cancel(resultant.subs({f:ff,u:uu,v:vv,d:dd})),am)
    assert actual_result.degree()==4
    lead_expected=S*S/(4*ax*ax*ay*ay)
    assert s.cancel(actual_result.LC()-lead_expected)==0
    #Literal identity with both frozen polynomial22 equations; K and B are retained
    #as short finite circuits, not expanded as high-degree beta invariants.
    ax_to_original={ax:s.Symbol('x'),ay:s.Symbol('y'),aa:s.Symbol('a'),am:s.Symbol('m')}
    compact=actual['compact_coefficients']
    substitutions={x:s.Symbol('x'),y:s.Symbol('y'),a:s.Symbol('a'),c:cc.subs(ax_to_original),
                   d:dd.subs(ax_to_original),f:ff.subs(ax_to_original),u:uu.subs(ax_to_original),
                   v:vv.subs(ax_to_original),K:s.sympify(compact['K1']),B:s.sympify(compact['constant'])}
    assert s.cancel(Q1.subs(substitutions,simultaneous=True)-s.sympify(actual['polynomial22_equation1']))==0
    assert s.cancel(Q2.subs(substitutions,simultaneous=True)-s.sympify(actual['polynomial22_equation2']))==0
    #The k1-independent coefficient A0 makes K=k1/2-A0 transcendental over the
    #other generic base/m coordinates. Verify the exact source coefficient.
    k1=s.Symbol('k1');Kactual=s.sympify(compact['K1'])
    assert s.diff(Kactual,k1)==s.Rational(1,2)
    assert s.diff(s.sympify(compact['constant']),k1)==0
    note=ROOT/'research/notes'/(ID+'.md')
    evidence={
        'generic_field':'QQ(x,y,a,k2,m)(K), K=k1/2-A0; x*y*(x^2-y^2)*m*c*f*u*v nonzero generically',
        'auxiliary_polynomials':{name:str(expression)for name,expression in [('A',A),('b',b),('q',q),('L',L),('M',M),('Q1',Q1),('Q2',Q2),('P8',P)]},
        'relation':'m=u-v',
        'linear_T_identity':'y*R*Q2+A*Q1=L*T+M',
        'whole_ideal_identity':'y*R*P8=L^2*Q1-a1*(L*T+M)^2-(b*L-2*a1*M)*(L*T+M); a1=-c*y*R',
        'T_recovery':'T=-y*m*R*(v*R^2+f*u)/L',
        'octic_coefficients_R0_through_R8':[str(s.factor(octic.nth(ii)))for ii in range(9)],
        'octic_degree':8,'R_leading':'x^2*v^4','R_constant':'x^2*f^4*u^4',
        'A_q_resultant':str(resultant),'actual_resultant_m_degree':4,
        'actual_resultant_m_leading':str(lead_expected),
        'generic_L_unit_proof_in_note':True,
        'inverse_L_finite_recipe':'adjugate of multiplication-by-L on K[R]/(P8) divided by its nonzero determinant',
        'no_coefficient_expansion_as_beta_polynomial':True,
    }
    ev=output/'certificate.json';ev.write_text(json.dumps(evidence,indent=2)+'\n')
    files=[Path(__file__),note,review,primary,certpath,
           ROOT/'research/notes/Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63.md',
           ROOT/'research/notes/Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63_REVIEW.md']
    report={
        'id':ID,'all_passed':True,
        'small_universal_linear_T_identity_and_whole_ideal_equalities':True,
        'literal_binding_to_frozen_actual_fixed_base_equations':True,
        'generic_L_and_R_and_T_units_proved_without_point_or_GB':True,
        'total_generic_intersection_has_monogenic_degree8_presentation':True,
        'selected_actual_component_degree_exact8_or_genus_claim':False,
        'total_intersection_reduced_irreducible_or_eight_geometric_points_claim':False,
        'NEW25_multisection_degree_or_equation_claim':False,
        'all_fixed_base_specializations_or_all_m_L_unit_claim':False,
        'no_beta_point_or_Inv24_evaluation_or_GB':True,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},'checked_dependency_sha256':checked,
        'evidence_sha256':{ev.name:sha(ev)},'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'seconds':time.monotonic()-started,
    }
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    signal.alarm(0);print(json.dumps({'all_passed':True,'total_algebra_dimension':8,'seconds':report['seconds']}))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(parser.parse_args().output.resolve())

#!/usr/bin/env python3
"""Fixed spectral-base curve: two explicit (2,2) equations and one old-point rank6."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import resource
import signal
import sys
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63'
sys.set_int_max_str_digits(0)


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def rows(M):
    return [[str(x) for x in row] for row in M.tolist()]


def piv(M):
    return list(DomainMatrix.from_Matrix(M).convert_to(s.QQ).rref()[1])


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    out.mkdir(parents=True, exist_ok=False)
    normalpath = ROOT / 'reports/Q_FIRSTJET_TWO_ADJOINT_NORMAL_FORM_20261103_v63/certificate.json.gz'
    pairpath = ROOT / 'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63/certificate.json.gz'
    dompath = ROOT / 'reports/Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63/certificate.json.gz'
    invpath = ROOT / 'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    normal, pair, dom, info = map(read, [normalpath, pairpath, dompath, invpath])
    x,y,a,g,h,i,j,t,z,k1,k2,m,R,T = s.symbols('x y a g h i j t z k1 k2 m R T')
    Delta = x*x-y*y
    S = x*x+y*y
    ell = (a+x)**2-y*y
    d = -(a*x+S)/y
    c, f = m-ell, m-x*x*ell/y**2
    C = s.Matrix([[a,1,g,h],[c,-a,i,j],[-j,h,d,1],[i,-g,f,-d]])
    D = s.diag(x,-x,y,-y)
    original = C.subs(m,g*j-h*i)
    detpencil = s.Poly(s.cancel((t*original+z*D).det(method='domain-ge')),t,z)
    p = [s.factor(detpencil.coeff_monomial(t**(4-k)*z**k)) for k in range(5)]
    p2constant = -4*a*a*x*x-2*a*x**3-6*a*x*y*y+3*x**4+y**4
    assert s.expand(p[2]-p2constant-(x+y)**2*g*j+(x-y)**2*h*i)==0
    assert s.expand(p[3]+2*a*x*Delta+2*x*x*S)==0 and p[4]==x*x*y*y
    vark3=s.Symbol('k3')
    a_from_p3=-(vark3+2*x*x*S)/(2*x*Delta)
    assert s.cancel(p[3].subs(a,a_from_p3)-vark3)==0
    N1,den1=s.fraction(s.cancel(p[1]-k1))
    N0,den0=s.fraction(s.cancel(p[0]-Delta**2))
    assert den1==den0==y*y
    assert s.expand(N0-s.sympify(normal['final_polynomial']))==0
    K2=k2-p2constant
    u=(K2-(x-y)**2*m)/(4*x*y)
    v=(K2-(x+y)**2*m)/(4*x*y)
    assert s.cancel(u-v-m)==0 and s.cancel((x+y)**2*u-(x-y)**2*v-K2)==0
    A0=a*d*(a*y+d*x)+(a*y+d*x)*(u+v)+c*d*y+f*a*x
    K1=k1/2-A0
    constant=(a*a+c)*(d*d+f)+m*m+2*a*d*(u+v)-Delta**2
    L1=-c*y*T-x*v*R+f*x*u/R+y*u*v/T-K1
    L2=constant+c*R*T+f*u*u/(R*T)-c*f*T/R-v*v*R/T-2*c*d*T+2*a*f*u/R-2*a*v*R+2*d*u*v/T
    Q1=s.Poly(s.expand(R*T*L1),R,T)
    Q2=s.Poly(s.expand(R*T*L2),R,T)
    assert Q1.degree(R)==Q1.degree(T)==Q2.degree(R)==Q2.degree(T)==2
    # Direct quotient homomorphism on Laurent monomials; no square-root expansion.
    def reduce_even(expression):
        expanded=s.expand(expression.subs({j:u/g,i:v/h},simultaneous=True))
        result=0
        for term in s.Add.make_args(expanded):
            powers=term.as_powers_dict()
            eg,eh=int(powers.get(g,0)),int(powers.get(h,0))
            assert (eg-eh)%2==0 and (eg+eh)%2==0
            coefficient=s.cancel(term/(g**eg*h**eh))
            assert not coefficient.has(g,h)
            result+=coefficient*R**((eg-eh)//2)*T**((eg+eh)//2)
        return s.cancel(result)
    generic_p1=2*a*a*d*y+2*a*d*d*x+2*a*g*j*y+2*a*h*i*y+c*(2*d*y-2*g*h*y)\
               +2*d*g*j*x+2*d*h*i*x+f*(2*a*x+2*h*j*x)-2*g*i*x+2*i*j*y
    assert s.cancel(reduce_even(generic_p1-k1)-2*L1)==0
    assert s.cancel(reduce_even(C.det(method='domain-ge')-Delta**2)-L2)==0
    # Homogeneous boundary checks on P1_R x P1_T.
    R0,R1,T0,T1=s.symbols('R0 R1 T0 T1')
    def homogenize(P):
        return sum(co*R1**er*R0**(2-er)*T1**et*T0**(2-et) for (er,et),co in P.terms())
    H1,H2=map(homogenize,[Q1,Q2])
    corners=[{R0:1,R1:0,T0:1,T1:0},{R0:1,R1:0,T0:0,T1:1},
             {R0:0,R1:1,T0:1,T1:0},{R0:0,R1:1,T0:0,T1:1}]
    expected=[f*u*u,-c*f,-v*v,c]
    assert all(s.cancel(H1.subs(corner))==0 for corner in corners)
    assert all(s.cancel(H2.subs(corner)-e)==0 for corner,e in zip(corners,expected))
    boundary=[s.factor(H1.subs({R0:1,R1:0})),s.factor(H1.subs({R0:0,R1:1})),
              s.factor(H1.subs({T0:1,T1:0})),s.factor(H1.subs({T0:0,T1:1}))]
    # Intrinsic rational m and its SAME-point universal-W derivative, not NEW-source rank6.
    C0,D0=s.Matrix(pair['actual_C4']),s.Matrix(pair['actual_D4'])
    basis=[s.Matrix(A) for A in pair['pair_coordinate_basis4']]
    def trword(word,variation=None):
        A=s.eye(4)
        for k,ch in enumerate(word):
            A=A*(variation[1] if variation and variation[0]==k else {'C':C0,'D':D0}[ch])
        return s.trace(A)
    def derivative(word,letter,E):
        return sum(trword(word,(k,E)) for k,ch in enumerate(word) if ch==letter)
    qC,qD,U4,U6,V4=[trword(w) for w in ['CC','DD','CCDD','CDDCDD','DDDD']]
    S0=qD/2
    P0=qD*qD/8-V4/4
    dis=S0*S0-4*P0
    numerator=U6-S0*U4+P0*qC
    intrinsic_m=numerator/(2*dis)
    assert P0==D0.det() and dis
    # Verify intrinsic expression identically in the full diagonal-D symbolic chart.
    qCs=s.trace(C*C);qDs=s.trace(D*D);Ssymbol=qDs/2;Psymbol=D.det();ds=Ssymbol*Ssymbol-4*Psymbol
    intrinsic_symbolic=(s.trace(C*D*D*C*D*D)-Ssymbol*s.trace(C*C*D*D)+Psymbol*qCs)/(2*ds)
    assert s.cancel(intrinsic_symbolic-(g*j-h*i))==0
    gradient=s.zeros(1,20)
    for column,letter in enumerate(['C']*10+['D']*10):
        E=basis[column%10]
        dqC,dqD,dU4,dU6,dV4=[derivative(w,letter,E) for w in ['CC','DD','CCDD','CDDCDD','DDDD']]
        dS=dqD/2;dP=qD*dqD/4-dV4/4
        ddis=2*S0*dS-4*dP
        dnum=dU6-dS*U4-S0*dU4+dP*qC+P0*dqC
        gradient[0,column]=s.cancel((dnum*dis-numerator*ddis)/(2*dis*dis))
    gradbeta=gradient*s.Matrix(pair['jacobian20x16'])
    six=s.Matrix(dom['trace8_gradient5x16_beta']).col_join(gradbeta)
    columns=piv(six)
    rank=len(columns)
    assert rank<=6
    minor=six[:,columns] if rank==6 else None
    determinant=minor.det(method='domain-ge') if rank==6 else None
    if rank==6:assert determinant
    b0=s.Matrix(dom['beta16'])
    orbit=s.Matrix.hstack(*(s.Matrix(A)*b0 for A in info['all_ten_representation_matrices']))
    assert six*orbit==s.zeros(6,10)
    # The same NEW source still has five invariant directions; no rank6 claim there.
    new_six=six*s.Matrix(dom['full_beta16x18_source_Jacobian'])
    selected=dom['selected_source_columns']
    alpha=new_six[5,selected]*new_six[:5,selected].inv()
    assert alpha*new_six[:5,:]==new_six[5,:]
    cert=dict(fixed_base_parameters=['x','y','a','k1=p1','k2=p2'],
              spectral_base_dictionary=dict(F4=str(20*S),p0=str(Delta**2),p3=str(p[3]),p4=str(p[4])),
              a_from_p3=str(a_from_p3), full_pencil_coefficients=list(map(str,p)),
              original_four_B_variables=['g','h','i','j'],three_fixed_base_equations=[str(N1),str(p[2]-k2),str(N0)],
              equation_denominators=[str(den1),'1',str(den0)],
              auxiliary_m_dictionary=dict(c=str(c),f=str(f),d=str(d),u_gj=str(u),v_hi=str(v)),
              relative_sign_quotient_coordinates=dict(R='g/h',T='g*h',g_squared='R*T',h='T/g',i='v*g/T',j='u/g'),
              Laurent_equation1=str(L1),Laurent_equation2=str(L2),
              polynomial22_equation1=str(Q1.as_expr()),polynomial22_equation2=str(Q2.as_expr()),
              compact_coefficients=dict(A0=str(A0),K1=str(K1),constant=str(constant)),
              generic_known_units='x*y*(x^2-y^2)*g*h*u*v*c*f*R*T',
              homogeneous_boundary_L1=list(map(str,boundary)),
              homogeneous_L2_at_four_corners=list(map(str,expected)),
              intrinsic_m_formula='(tr(CD^2CD^2)-S tr(C^2D^2)+P tr(C^2))/(2(S^2-4P)); S=tr(D^2)/2,P=detD',
              actual_same_point_intrinsic_m=str(intrinsic_m),intrinsic_m_gradient1x20=rows(gradient),
              intrinsic_m_gradient1x16_beta=rows(gradbeta),universal_six_gradient6x16=rows(six),
              universal_gradient_rank=rank,selected_beta_columns=columns,
              exact6_minor=rows(minor) if minor is not None else None,
              exact6_minor_determinant=str(determinant) if determinant is not None else None,
              same_NEW_source_six_gradient6x18=rows(new_six),NEW_source_m_differential_linear_combination=rows(alpha))
    with gzip.open(out/'certificate.json.gz','wt') as stream:
        json.dump(cert,stream,indent=2);stream.write('\n')
    result=dict(id=STEM,all_passed=True,full_three_fixed_base_equations_verified=True,
                relative_sign_quotient_two_Laurent_equations_verified=True,
                both_polynomial_equations_bidegree=[2,2],
                generic_proper_P1xP1_intersection_length=8,
                actual_selected_component_degree_bound_conditional_on_rank6=8,
                same_existing_point_only=True,universal_W16_six_invariant_gradient_rank=rank,
                intrinsic_m_nonconstant_on_generic_spectral_fiber=rank==6,
                NEW_source_invariant_differential_rank_stays5=True,
                actual_selected_component_degree_exact8_or_genus_claimed=False,
                no_new_beta_points_no_degree24_probe_no_GB=True,
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),normalpath,pairpath,dompath,invpath]},
                seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['all_passed','universal_W16_six_invariant_gradient_rank',
          'intrinsic_m_nonconstant_on_generic_spectral_fiber','seconds']},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output.resolve())

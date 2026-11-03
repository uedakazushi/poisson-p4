#!/usr/bin/env python3
"""Certify the precise gap after a=0 and d=0 marked-stratum removals.

The constructed finite equivariant incidence is NOT a Jacobi solution family.
It shows why the existing intrinsic/finite-anchor data alone do not exclude
remaining marked nonsquare points.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def verify():
    started=time.time()
    A,B,C,D,h,k=s.symbols('A B C D h k')
    variables=[A,B,C,D,h,k]
    r,t=s.symbols('r t')
    P=A*r**3+B*r*r*t+C*r*t*t+D*t**3
    ell=h*r+k*t
    marked=s.Matrix([-k,h])
    at=dict(zip([r,t],marked))
    gradient=s.Matrix([P.diff(r),P.diff(t)]).subs(at,simultaneous=True)
    other=s.Matrix([gradient[1],-gradient[0]])
    resultant=s.expand(P.subs(at,simultaneous=True))
    hessian=-s.expand(s.hessian(P,[r,t]).det().subs(at,simultaneous=True))/4
    discriminant=s.discriminant(P.subs(t,1),r)
    polar_value=s.expand(P.subs(dict(zip([r,t],other)),simultaneous=True))
    sextic,remainder=s.div(polar_value,resultant,*variables)
    assert remainder==0
    sextic=s.expand(sextic)
    assert s.expand(sextic**2-4*hessian**3+27*resultant**2*discriminant)==0
    assert s.expand(s.det(s.Matrix.hstack(marked,other))+3*resultant)==0
    # The unshifted polar lies on the previously removed c=f=0 stratum.
    # Shift along the marked vector with matching homogeneous degree to
    # obtain an example avoiding that stratum as well.
    shifted=sextic*other-resultant**2*marked
    shifted_value=s.expand(P.subs(dict(zip([r,t],shifted)),simultaneous=True))
    assert s.expand(shifted_value-resultant*(sextic**4+3*resultant**2*hessian*sextic**2-resultant**6))==0
    assert s.expand(gradient.dot(shifted)+3*resultant**3)==0
    assert s.expand(s.det(s.Matrix.hstack(marked,shifted))+3*resultant*sextic)==0
    f=s.expand(hessian-resultant)
    for polynomial,degree in [(resultant,4),(hessian,4),(discriminant,4),(sextic,6),(f,4)]:
        assert s.Poly(polynomial,*variables).total_degree()==degree
        assert all(sum(mon)==degree for mon in s.Poly(polynomial,*variables).monoms())
    actions=[]
    for dr,dt in [(t,0),(0,r)]:
        dp=s.Poly(dr*P.diff(r)+dt*P.diff(t),r,t)
        dl=s.Poly(dr*ell.diff(r)+dt*ell.diff(t),r,t)
        actions.append([dp.coeff_monomial(r**(3-i)*t**i) for i in range(4)]+
                       [dl.coeff_monomial(r),dl.coeff_monomial(t)])
    for vector in actions:
        for invariant in [resultant,hessian,discriminant,sextic,f]:
            assert s.expand(sum(v*invariant.diff(u) for u,v in zip(variables,vector)))==0
    # Absolute irreducibility: f is primitive and linear in A.  Its A
    # coefficient is k*L, where L is primitive linear in C over any char0
    # constant field.  Neither k nor L divides the A-free remainder.
    coefficient=s.diff(f,A)
    constant=f.subs(A,0)
    L=-3*C*k+9*D*h+k*k
    assert s.expand(coefficient-k*L)==0
    assert constant.subs({k:0,h:1,C:0,D:1,B:0})==-1
    assert L.subs({k:1,h:0,C:s.Rational(1,3)})==0
    assert constant.subs({k:1,h:0,C:s.Rational(1,3),B:1})==1
    assert s.gcd(coefficient,constant)==1
    sample={A:1,B:0,C:1,D:1,h:1,k:0}
    assert f.subs(sample)==0
    assert resultant.subs(sample)==1 and hessian.subs(sample)==1
    assert sextic.subs(sample)==29 and discriminant.subs(sample)==-31
    assert list(marked.subs(sample))==[0,1]
    assert list(other.subs(sample))==[3,-1]
    assert polar_value.subs(sample)==29
    assert list(shifted.subs(sample))==[87,-30]
    assert shifted_value.subs(sample)==709803
    assert gradient.dot(shifted).subs(sample)==-3
    sample_marked=marked.subs(sample)
    sample_other=shifted.subs(sample)
    determinant=s.det(s.Matrix.hstack(sample_other,sample_marked))
    frame=s.Matrix.hstack(sample_other/determinant,sample_marked)
    assert frame.det()==1
    transported=s.Poly(P.subs(sample).subs(dict(zip([r,t],frame*s.Matrix([r,t]))),simultaneous=True),r,t)
    coefficients=[transported.coeff_monomial(r**(3-i)*t**i) for i in range(4)]
    assert ell.subs(sample).subs(dict(zip([r,t],frame*s.Matrix([r,t]))),simultaneous=True)==r
    aa,BB,CC,DD=coefficients
    LL=-s.S.One
    cc=-CC/6
    bb=(BB+2*LL)/9
    ee=(5*BB+LL)/9
    dd=-DD
    ff=-5*cc
    assert aa*cc*dd*LL !=0
    deeper=600*bb**2-255*bb*ee-63*ee**2+100*aa*ff
    assert deeper !=0
    old_path=ROOT/'reports/S61_AFFINE_ANCHOR_RESIDUE_DEGREE_20261102.json'
    old=json.loads(old_path.read_text())
    assert old['all_passed']
    rename=dict(zip(s.symbols('x0:6'),variables))
    quartics=[s.sympify(x).subs(rename) for x in old['binary_pair_SL2_invariant_polynomials']['4']['polynomials']]
    old_sextic=s.sympify(old['binary_pair_SL2_invariant_polynomials']['6']['polynomials'][0]).subs(rename)
    assert all(s.expand(x-y)==0 for x,y in zip(quartics,[discriminant,hessian,resultant]))
    assert s.expand(old_sextic-sextic)==0
    # SL2 covariance of both the marked vector and the polar vector follows
    # from J*g^T=adj(g)*J and the chain rule.  Verify its matrix identity.
    ga,gb,gc,gd=s.symbols('ga gb gc gd')
    g=s.Matrix([[ga,gb],[gc,gd]])
    J=s.Matrix([[0,1],[-1,0]])
    assert J*g.T==g.adjugate()*J
    files=[Path(__file__),old_path,
           ROOT/'research/notes/S61_AFFINE_ANCHOR_RESIDUE_DEGREE_20261102.md',
           ROOT/'research/notes/S61_AFFINE_BOUNDARY_FINITE_ANCHORS_20261102.md']
    return dict(id='S61_NONSQUARE_MARKED_POLAR_OBSTRUCTION_20261103',all_passed=True,field='QQ',
                cubic=str(P),linear=str(ell),marked_vector=list(map(str,marked)),polar_vector=list(map(str,other)),
                shifted_polar_vector=list(map(str,shifted)),
                invariants=dict(resultant=str(resultant),hessian=str(hessian),discriminant=str(discriminant),sextic=str(sextic)),
                verified_identities=['P(polar)=resultant*sextic','det(marked,polar)=-3*resultant',
                                     'sextic^2=4*hessian^3-27*resultant^2*discriminant',
                                     'P(shifted)=resultant*(sextic^4+3*resultant^2*hessian*sextic^2-resultant^6)',
                                     'gradientP(marked).shifted=-3*resultant^3',
                                     'det(marked,shifted)=-3*resultant*sextic'],
                invariant_countermodel_support_divisor=str(f),support_divisor_absolutely_irreducible=True,
                sample={str(v):str(z) for v,z in sample.items()},sample_polar_value='29',
                sample_shifted_polar_value='709803',sample_core=dict(a=str(aa),b=str(bb),c=str(cc),d=str(dd),e=str(ee),f=str(ff)),
                sample_deeper_condition=str(deeper),
                model_generic_support_plane_degree=10,
                conclusion='Intrinsic finite-anchor geometry plus a=0,d=0 removal does not alone force all nonsquare marked anchors empty.',
                limitations=['This is a geometric countermodel to a proposed inference, NOT a Jacobi solution.',
                             'The unshifted polar is c=f=0; the shifted polar is used to avoid that already removed stratum.',
                             'No containment in the actual residual weighted-rank defect locus is claimed.',
                             'It establishes no actual remaining boundary point or new Poisson component.',
                             'The actual localized saturated Jacobi ideal still requires exclusion.'],
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                seconds=time.time()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=verify()
    with args.output.open('x') as output:
        json.dump(report,output,indent=2)
        output.write('\n')
    print(json.dumps({'all_passed':True,'polar_identity':True,'countermodel_is_not_Jacobi':True,'seconds':report['seconds']}))

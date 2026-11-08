#!/usr/bin/env python3
"""Exact independent check of the constant-TRIPLE -> Ga9 bridge (v1).

Only SymPy is required. No archived producer, report, or numerical sampling is
used. This checks the displayed family, not global coverage or componenthood.
Run: python code/verify_final_bridge.py --output reports/final_bridge.json
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import platform
import sympy as s


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    checks: dict[str, int] = {}

    def zero(name: str, expressions) -> None:
        expressions = list(expressions)
        for i, expression in enumerate(expressions):
            remainder = s.cancel(s.expand(expression))
            if remainder != 0:
                raise AssertionError(f'{name}[{i}]: {remainder}')
        checks[name] = len(expressions)

    def bracket(f, g, coordinates, matrix):
        return s.expand(sum(s.diff(f, coordinates[i])*matrix[i,j]*s.diff(g, coordinates[j])
                            for i in range(len(coordinates)) for j in range(len(coordinates))))

    def skew(size, upper):
        matrix = s.zeros(size)
        for (i,j), value in upper.items():
            matrix[i,j] = value
            matrix[j,i] = -value
        return matrix

    def jacobi(coordinates, matrix):
        return [sum(bracket(coordinates[i], matrix[j,k], coordinates, matrix)
                    for i,j,k in [(a,b,c),(b,c,a),(c,a,b)])
                for a,b,c in itertools.combinations(range(len(coordinates)),3)]

    e,T,x,y,S,t,r = s.symbols('e T x y S t r')
    H=x**3+3*x*y-s.Rational(3,2)*y
    k=x*y/3+y**2/2-x**4/12
    A=e/12
    B=e*(15-e)/72
    M=A*x*x*y+B*y*y
    old=(T,x,y,S)
    new=(t,x,y,r)
    forward=((1-e*x/3)*T-x*S/3-k+M,x,y,S+e*T-H)
    inverse=(t+x*(r+H)/3+k-M,x,y,r+H-e*(t+x*(r+H)/3+k-M))
    zero('polynomial_inverse',
         [v.subs(dict(zip(old,inverse)), simultaneous=True)-w for v,w in zip(forward,new)]+
         [v.subs(dict(zip(new,forward)), simultaneous=True)-w for v,w in zip(inverse,old)])
    zero('jacobian_determinant', [s.Matrix(forward).jacobian(old).det()-1])
    product=skew(4,{(0,3):S,(1,2):1})
    actual=skew(4,{(i,j):bracket(forward[i],forward[j],old,product).subs(dict(zip(old,inverse)),simultaneous=True)
                   for i,j in itertools.combinations(range(4),2)})
    actual=actual.applyfunc(s.expand)
    claimed=skew(4,{
        (0,1):x/3+(1-2*B)*y-A*x*x,
        (0,2):(-1+2*A)*x*y+y/6-r/3,
        (0,3):-e*t+(1-e/3)*x*r+(-3+3*e/4+e**2/8-e**3/72)*y*y+5*e*x*y/12-5*y/4+r/2,
        (1,2):1,(1,3):-3*x+s.Rational(3,2),(2,3):3*x*x+3*y})
    zero('six_brackets_from_product',[actual[i,j]-claimed[i,j] for i,j in itertools.combinations(range(4),2)])
    degrees=[s.Poly(actual[i,j],*new).total_degree() for i,j in itertools.combinations(range(4),2)]
    if max(degrees)>2:
        raise AssertionError('Nonquadratic bracket')
    checks['quadratic_degree_bounds']=6
    zero('affine_jacobi',jacobi(new,actual))
    pf=actual[0,1]*actual[2,3]-actual[0,2]*actual[1,3]+actual[0,3]*actual[1,2]
    zero('pfaffian',[pf-inverse[3]])
    div=[sum(s.diff(actual[i,j],new[i]) for i in range(4)) for j in range(4)]
    zero('modular_field',[div[i]-v for i,v in enumerate((-1+e*x/3,0,0,-e))])

    # A homogeneous quadratic Poisson tensor on A^5 with h central descends
    # to P^4. Verify strict Jacobi and dehomogenization, independently of the
    # repository's 150/350 convention; all coefficients are exact in Q[e].
    h=s.Symbol('h')
    homogeneous=skew(5,{(i+1,j+1):s.expand(h*h*actual[i,j].subs({v:v/h for v in new},simultaneous=True))
                        for i,j in itertools.combinations(range(4),2)})
    if any(not s.Poly(v,h,*new,e).is_zero and s.Poly(v,h,*new).total_degree()!=2 for v in homogeneous):
        raise AssertionError('Homogenization failed')
    zero('homogeneous_strict_jacobi',jacobi((h,)+new,homogeneous))
    zero('dehomogenization',[homogeneous[i+1,j+1].subs(h,1)-actual[i,j] for i,j in itertools.combinations(range(4),2)])
    D=s.Poly(s.expand(h**5*pf.subs({v:v/h for v in new},simultaneous=True)),h,*new)
    zero('projective_pfaffian',[s.expand(h*(homogeneous[1,2]*homogeneous[3,4]-homogeneous[1,3]*homogeneous[2,4]+homogeneous[1,4]*homogeneous[2,3]))-D.as_expr()])
    if s.expand(actual[1,2]) != 1:
        raise AssertionError('Family may pass through zero tensor')
    checks['nonzero_all_parameters']=1
    zero('special_triple_divisor',[D.as_expr().subs(e,0)-h*h*(h*h*r+x**3+3*h*x*y-s.Rational(3,2)*h*h*y)])
    F4=s.cancel(D.as_expr()/h)
    zero('generic_quartic_coefficients',[s.diff(F4,t)+e*h**3,
        s.Poly(F4,h,t,x,y,r).coeff_monomial(x**4)+e/4])

    # Ga9 is reconstructed from its potential on (u,Y,Z,w), not from the
    # printed six specialized brackets.
    X,Y,Z,w,u=s.symbols('X Y Z w u')
    transformed=(t-r/e,e*x/3,-r/e,3*y/e)
    back=(X-Z,3*Y/e,e*w/3,-e*Z)
    zero('affine_Ga9_inverse',
         [v.subs(dict(zip(new,back)),simultaneous=True)-q for v,q in zip(transformed,(X,Y,Z,w))]+
         [v.subs(dict(zip((X,Y,Z,w),transformed)),simultaneous=True)-q for v,q in zip(back,new)])
    zero('Ga9_frame_determinant',[s.Matrix(transformed).jacobian(new).det()-1/e])
    pushed=skew(4,{(i,j):s.cancel(bracket(transformed[i],transformed[j],new,actual).subs(dict(zip(new,back)),simultaneous=True))
                  for i,j in itertools.combinations(range(4),2)})
    rho=81/e**4
    b=3/e
    alpha=-s.Rational(1,2)
    aa=-27/e**4
    bb=b-s.Rational(1,4)
    ee=alpha-(1-3*b)*aa/rho
    ff=-(1-b)*bb/(2*rho)
    phi=rho*Y**4/4+aa*Y**3+bb*Y**2*w+ee*Y*w+ff*w*w+w/2
    Vy=-alpha-b*Y
    Vw=b*w+rho*Y*Y
    A2=s.expand(phi+Vy*s.diff(phi,Y)+Vw*s.diff(phi,w))
    Vu=-u-A2
    base=(u,Y,Z,w)
    base_bracket=skew(4,{(0,1):s.diff(phi,w),(0,2):-Vu,(0,3):-s.diff(phi,Y),
                         (1,2):-Vy,(1,3):1,(2,3):Vw})
    canonical=(u+Y*Z,Y,Z,w)
    target=skew(4,{(i,j):s.cancel(bracket(canonical[i],canonical[j],base,base_bracket).subs(u,X-Y*Z))
                   for i,j in itertools.combinations(range(4),2)})
    zero('Ga9_potential_full_six',[pushed[i,j]-target[i,j] for i,j in itertools.combinations(range(4),2)])
    zero('Ga9_parameter_identities',[ee+b+s.Rational(1,6),ff-e*e*(e-3)*(e-12)/648])
    # The exceptional F2 zeros are NOT excluded by the r != 0 source domain.
    zero('Ga9_F2_zero_parameters',[expr.subs(e,value) for value in (3,12)
                                  for expr in [pushed[i,j]-target[i,j] for i,j in itertools.combinations(range(4),2)]])
    if any(s.denom(s.cancel(v)).subs(e,3)==0 or s.denom(s.cancel(v)).subs(e,12)==0 for v in target):
        raise AssertionError('Unexpected pole at e=3 or e=12')
    if s.cancel(rho*e**4) != 81:
        raise AssertionError('Ga9 source condition rho != 0 failed')
    for value in list(target)+[rho,b,alpha,aa,ee,ff]:
        denominator = s.Poly(s.denom(s.cancel(value)),e)
        if len(denominator.terms()) != 1:
            raise AssertionError('Unexpected denominator away from e=0')
    checks['Ga9_source_nonzero_rho_and_only_e_poles']=1
    record={
        'date':'2026-11-09','version':'v1','status':'PASS',
        'scope':'Exact symbolic verification of one polynomial family and its Ga9 source identification. NOT a verification of global 147-component classification, source coverage, or Ga9 componenthood.',
        'base_commit':'ffc178bff0a69fb0edcc99a4d8d7de98c8a2a222',
        'field':'Q(e); polynomial identities checked identically, not by samples',
        'python':platform.python_version(),'sympy':s.__version__,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'checks':checks,'check_count':sum(checks.values()),'coefficient_degrees':degrees,
        'pfaffian':str(s.expand(pf)),
        'homogeneous_divisor':str(D.as_expr()),
        'Ga9_parameters':{'rho':str(rho),'b':str(b),'alpha':str(alpha),'A':str(aa),'H':'1/2','p,q,D,G,k':'0','E':str(s.cancel(ee)),'F2':str(s.factor(ff))},
        'not_checked':['Coverage by the four original TRIPLE sources','Ga9 ten-source exhaustiveness and componenthood','Prior 142 components','T37 dense-orbit residue identities','Normal and nonnormal quintic recognition','Full classification manuscript and PDF rebuild']}
    text=json.dumps(record,ensure_ascii=False,indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text,encoding='utf-8')
    print(text)


if __name__=='__main__':
    main()

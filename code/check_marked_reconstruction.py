#!/usr/bin/env python3
"""Exact checks of the shared modular primitive identity and two lifting families.

No old producers are imported. This is not a verifier of the 147 components.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sympy as s

PAIRS = list(itertools.combinations(range(4), 2))
TRIPLES = list(itertools.combinations(range(4), 3))


def skew(values):
    p = s.zeros(4)
    for (i, j), value in zip(PAIRS, values):
        p[i, j], p[j, i] = value, -value
    return p


def pf(p):
    return s.expand(p[0, 1]*p[2, 3]-p[0, 2]*p[1, 3]+p[0, 3]*p[1, 2])


def dual(p):
    return skew([p[2, 3], -p[1, 3], p[1, 2], p[0, 3], -p[0, 2], p[0, 1]])


def div(p, coordinates):
    return s.Matrix([sum(s.diff(p[i, j], coordinates[i]) for i in range(4)) for j in range(4)])


def jac(p, coordinates):
    return s.Matrix([s.expand(sum(p[i, l]*s.diff(p[j, k], coordinates[l])
        + p[j, l]*s.diff(p[k, i], coordinates[l])
        + p[k, l]*s.diff(p[i, j], coordinates[l]) for l in range(4)))
        for i, j, k in TRIPLES])


def run():
    checks = []

    def check(name, values):
        for value in values:
            residual = s.cancel(s.expand(value))
            if residual != 0:
                raise AssertionError(f'{name}: {residual}')
        checks.append(name)

    # Universal independent first jets, including Pfaffian zero.
    vals = s.symbols('p0:6')
    dp = {v: s.symbols(f'd{j}_0:4') for j, v in enumerate(vals)}
    p = skew(vals)
    def derivative(g, k):
        return sum(s.diff(g, v)*dp[v][k] for v in vals)
    z = s.Matrix([sum(derivative(p[i,j], i) for i in range(4)) for j in range(4)])
    jj = {idx: s.expand(sum(p[i,l]*derivative(p[j,k],l)+p[j,l]*derivative(p[k,i],l)
        +p[k,l]*derivative(p[i,j],l) for l in range(4))) for idx in TRIPLES for i,j,k in [idx]}
    theta = s.Matrix([(-1)**(3-i)*jj[tuple(k for k in range(4) if k != i)] for i in range(4)])
    # Contraction i_Z alpha has column alpha.T*Z.
    check('universal_primitive_jacobi', s.Matrix([derivative(pf(p),k) for k in range(4)])+dual(p).T*z-theta)
    x = s.symbols('x0:4')
    ranktwo = skew([1, x[0], 0, 0, 0, 0])
    check('rank_two_primitive_excludes_counterexample', [pf(ranktwo),
        (dual(ranktwo).T*div(ranktwo,x))[3]-1, jac(ranktwo,x)[0]-1])

    t,x,y,r,e = s.symbols('t x y r e')
    coordinates=(t,x,y,r)
    H=x**3+3*x*y-s.Rational(3,2)*y
    k=x*y/3+y*y/2-x**4/12
    K=k+x*(r+H)/3
    p0=skew([x/3+y,-x*y+y/6-r/3,-3*y*y+x*r-5*y/4+r/2,1,-3*x+s.Rational(3,2),3*x*x+3*y])
    p1=skew([-5*y/12-x*x/12,x*y/6,-t-x*r/3+3*y*y/4+5*x*y/12,0,0,0])
    p2=skew([y/36,0,y*y/8,0,0,0])
    p3=skew([0,0,-y*y/72,0,0,0])
    pe=p0+e*p1+e*e*p2+e**3*p3
    f0=r+H
    phi=-t-K
    psi=x*x*y/12+5*y*y/24
    chi=-y*y/72
    fe=f0+e*phi+e*e*psi+e**3*chi
    check('triple_all_parameter_jacobi', jac(pe,coordinates))
    check('triple_prescribed_pfaffian', [pf(pe)-fe])
    z0=s.Matrix([-1,0,0,0]); z1=s.Matrix([x/3,0,0,-1])
    check('triple_modular',div(pe,coordinates)-z0-e*z1)
    check('triple_first_primitive_lift',s.Matrix([s.diff(phi,v) for v in coordinates])+dual(p1).T*z0+dual(p0).T*z1)
    check('triple_second_primitive_lift',s.Matrix([s.diff(psi,v) for v in coordinates])+dual(p2).T*z0+dual(p1).T*z1)
    check('triple_second_pfaffian_lift',[pf(p1),s.expand(pf(p0+e*p2)).coeff(e,1)-psi])
    trunc=s.expand(jac(p0+e*p1,coordinates)[1]).coeff(e,2)
    check('triple_omitted_second_correction',[trunc-x*x/12-5*y/12])
    if psi == 0 or trunc == 0:
        raise AssertionError('The required correction disappeared')
    A,B=s.symbols('A B')
    M=A*x*x*y+B*y*y
    high=s.expand(p0[0,3]-e*t-e*K+e*M-(s.diff(M,x)*s.diff(H,y)-s.diff(M,y)*s.diff(H,x)))
    check('triple_forced_degree_correction',[high.coeff(x,4)-(-e/4+3*A),
        high.coeff(x,2).coeff(y,1)-(-e+e*A-3*A+6*B)])
    check('triple_high_degree_cancellation',[high.subs({A:e/12,B:e*(15-e)/72})-pe[0,3]])
    T,ss=s.symbols('T ss')
    Me=M.subs({A:e/12,B:e*(15-e)/72})
    forward=s.Matrix([(1-e*x/3)*T-x*ss/3-k+Me,x,y,ss+e*T-H])
    invT=t+x*(r+H)/3+k-Me; invs=r+H-e*invT
    check('triple_polynomial_inverse_forward',forward.subs({T:invT,ss:invs},simultaneous=True)-s.Matrix(coordinates))
    inv=s.Matrix([invT,x,y,invs])
    check('triple_polynomial_inverse_backward',inv.subs(dict(zip(coordinates,forward)),simultaneous=True)-s.Matrix([T,x,y,ss]))
    check('triple_polynomial_jacobian',[forward.jacobian([T,x,y,ss]).det()-1])
    product=skew([0,0,ss,1,0,0])
    dforward=forward.jacobian([T,x,y,ss])
    pushed=(dforward*product*dforward.T).subs({T:invT,ss:invs},simultaneous=True)
    check('triple_product_pushforward',pushed-pe)
    # Signed cofactor convention: coefficient tx*xr-... checked above directly.
    X,Y,Z,w=s.symbols('X Y Z w')
    affine=s.Matrix([t-r/e,e*x/3,-r/e,3*y/e])
    inverse=s.Matrix([X-Z,3*Y/e,e*w/3,-e*Z])
    check('triple_actual_affine_inverse',affine.subs(dict(zip(coordinates,inverse)),simultaneous=True)-s.Matrix([X,Y,Z,w]))
    check('triple_actual_affine_jacobian',[affine.jacobian(coordinates).det()-1/e])
    mapped=(affine.jacobian(coordinates)*(z0+e*z1)).subs(dict(zip(coordinates,inverse)),simultaneous=True)
    check('triple_actual_affine_modular',mapped-s.Matrix([Y,0,1,0]))
    for exceptional in [0,3,12]:
        check(f'triple_parameter_{exceptional}',[pe[1,2].subs(e,exceptional)-1, fe.subs({e:exceptional,t:0,x:0,y:0,r:1})-1])

    x,y,z,w=s.symbols('x y z w'); coordinates=(x,y,z,w)
    pc=skew([-2*x*y,2*x*z,0,-2*y*z,0,z*(1+e*w)])
    fc=-2*x*y*z*(1+e*w)
    check('collision_all_parameter_jacobi',jac(pc,coordinates))
    check('collision_pfaffian',[pf(pc)-fc])
    direction=pc.diff(e)
    check('collision_no_second_correction',[pf(direction)]+list(jac(direction,coordinates)))
    check('collision_highest_degree',[s.Poly(fc,x,y,z,w).coeff_monomial(x*y*z*w)+2*e])
    # Recover the whole simple-H source (tau,a,b) and the degree-drop boundary.
    tau,a,b=s.symbols('tau a b')
    pq=skew([36*tau*w+36*a*z+36*b*w*z+6*x*z-8*y*y,
       36*a*w+36*b*w*w+6*w*x-10*y*z,-10*w*y-z*z,6*w*y-9*z*z,-12*w*z,-12*w*w])
    fq=9*((2*w*y-z*z)**2-48*tau*w**3)
    check('simple_H_whole_three_parameter_source',list(jac(pq,coordinates))+[pf(pq)-fq])
    check('simple_H_tau_zero_complement',[fq.subs(tau,0)-9*(2*w*y-z*z)**2])

    return {'date':'2026-11-09','version':'v68-marked-reconstruction','status':'PASS',
        'scope':'Universal modular primitive identity; entire displayed lifting curves and simple-H source. No 147-component verification.',
        'checks':checks,'sampling_used':False,'old_producers_imported':False,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=run()
    output=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(output)
    print(output)


if __name__=='__main__':
    main()

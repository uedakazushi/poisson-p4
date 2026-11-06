#!/usr/bin/env python3
"""Full irreducible inverse-form source for the tangent rank3 T(w)=0 branch.

The SL2 relation is imposed by the explicit chart in projective_chart().
No unknown Poisson parameter or residual inverse-form class is specialized.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def full_forms():
    u, v, z, w = x = s.symbols('u v z w')
    a, b, c, d = s.symbols('a b c d')
    delta, D, D0, S, C, g, T0, T = s.symbols('delta D D0 S C g T0 T')
    params = (a,b,c,d,delta,D,D0,S,C,g,T0,T)
    U, V, W = a*b+c*d*w, a*a+c*c*w, b*b+d*d*w
    ss = s.Matrix([2*a*b,2*a*a,S])
    dd = s.Matrix([D,delta,D0])
    by, ce = (dd-ss)/2, (dd+ss)/2
    gv, gz = (delta*W-D*U+2*w)/2, (delta*U-D*V)/2
    gw = (s.diff(gv,w)-1)*v+s.diff(gz,w)*z+g*w
    theta = T0+T*u+(by.T*s.Matrix([v,z,w]))[0]
    beta = s.Matrix([theta, u*ce[0]+gv, u*ce[1]+gz, u*ce[2]+gw])
    r = w*(c*d*v+c*c*z+C*w)
    gamma = beta+s.Matrix([u*s.diff(r,xx)-(r if i==0 else 0)
                           for i,xx in enumerate(x)])
    Q = w+u*u
    alpha = s.zeros(4)
    for i,j in itertools.combinations(range(4),2):
        alpha[i,j] = s.expand(Q*(s.diff(gamma[j],x[i])-s.diff(gamma[i],x[j]))
                            -s.diff(Q,x[i])*gamma[j]+s.diff(Q,x[j])*gamma[i])
        alpha[j,i] = -alpha[i,j]
    return dict(x=x,params=params,beta=beta,r=r,Q=Q,alpha=alpha,detSL2=a*d-b*c)


def projective_chart():
    """Actual polynomial Poisson family over a*delta != 0; all 11 directions.

    Return (P, homogeneous_coordinates, independent_parameters).  At the
    B!=0 witness set (a,b,c,delta,C)=(1,0,1,2,1), other parameters zero.
    The a!=0 chart is only for derivatives; full_forms and the proof retain a=0.
    """
    data = full_forms()
    x=data['x']; a,b,c,d,*remaining=data['params']
    alpha=data['alpha'].subs(d,(1+b*c)/a)
    coeffs=(alpha[2,3],-alpha[1,3],alpha[1,2],
            alpha[0,3],-alpha[0,2],alpha[0,1])
    h=s.symbols('h'); xx=(h,*x); P=s.zeros(5)
    for (i,j),f in zip(itertools.combinations(range(1,5),2),coeffs):
        P[i,j]=sum(aa*h**(2-sum(e))*s.prod(yy**ee for yy,ee in zip(x,e))
                   for e,aa in s.Poly(s.expand(f),*x).terms())
        P[j,i]=-P[i,j]
    return P,xx,(a,b,c,*remaining)


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    start=time.monotonic(); data=full_forms(); x=data['x']; Q=data['Q']
    alpha=data['alpha']; delta=data['params'][4]
    pf=s.expand(alpha[0,1]*alpha[2,3]-alpha[0,2]*alpha[1,3]+alpha[0,3]*alpha[1,2])
    target=-delta*Q*(x[0]**2+data['detSL2']**2*x[3])
    assert s.expand(pf-target)==0
    assert all(s.Poly(alpha[i,j],*x).total_degree()<=2
               for i,j in itertools.combinations(range(4),2))
    for i,j,k in itertools.combinations(range(4),3):
        da=s.diff(alpha[j,k],x[i])-s.diff(alpha[i,k],x[j])+s.diff(alpha[i,j],x[k])
        wedge=s.diff(Q,x[i])*alpha[j,k]-s.diff(Q,x[j])*alpha[i,k]+s.diff(Q,x[k])*alpha[i,j]
        assert s.expand(Q*da-2*wedge)==0
    P,xx,params=projective_chart()
    vals={p:0 for p in params}; vals.update(dict(zip(params[:3],(1,0,1))))
    vals[delta]=2; vals[s.Symbol('C')]=1
    witness=P.subs(vals)
    jac=[s.expand(sum(witness[i,l]*s.diff(witness[j,k],xx[l])
                      +witness[j,l]*s.diff(witness[k,i],xx[l])
                      +witness[k,l]*s.diff(witness[i,j],xx[l]) for l in range(5)))
         for i,j,k in itertools.combinations(range(5),3)]
    assert jac==[0]*10
    report=dict(task_id='E_quadric3',all_passed=True,
        full_source='SL2 x Gm x A7 after the beta+dQ kernel gauge',
        full_source_surjectivity='Analytic rank-one symmetric-pencil proof in the companion note',
        universal_alpha_degree_at_most_two=True,universal_closed_inverse_identity=True,
        universal_pfaffian=str(target),SL2_relation=str(data['detSL2']-1),
        beta=list(map(str,data['beta'])),r=str(data['r']),
        alpha={f'{i},{j}':str(alpha[i,j]) for i,j in itertools.combinations(range(4),2)},
        derivative_chart_parameters=list(map(str,params)),
        chart_open='a != 0; complement a=0 retained by the full SL2 source',
        nonzero_B_witness={f'{i},{j}':str(witness[i,j]) for i,j in itertools.combinations(range(5),2)},
        witness_direct_raw_Jacobi_zero=True,full_ambient_tangent_rank_not_computed=True,
        known_component_coverage_not_proved=True,
        input_sha256={'code/classification_rank3_twzero_source.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        elapsed_seconds=time.monotonic()-start)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=report['elapsed_seconds'],output=str(output))))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output)

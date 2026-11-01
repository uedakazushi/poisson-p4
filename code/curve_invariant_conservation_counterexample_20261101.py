#!/usr/bin/env python3
"""A global quadratic Poisson counterexample to unrestricted c4/c6 conservation.

Its canonical curve is singular (Delta=0), so this does not contradict a
possible conservation theorem after localizing at Delta.
"""
import argparse
import hashlib
import json
from itertools import product
from pathlib import Path

import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

from elliptic_curve_invariants_20261101 import pfaffOmega
from elliptic_weight1_dominance import primitive_vector
from s61_model import S61, X, jacobi_projective, pushforward

ROOT=Path(__file__).resolve().parents[1]


def univariate_invariants(omega,x,parameter):
    """Fisher's contractions over QQ[t], not numerical interpolation."""
    polynomial_ring,_=ring(str(parameter),QQ)
    rr=range(5)
    w={key:polynomial_ring.from_expr(s.diff(omega[key[0],key[1]],x[key[2]],x[key[3]]))
       for key in product(rr,repeat=4)}
    z={(i,j,a,b):sum((w[i,r,u,a]*w[j,u,r,b]+w[i,r,u,b]*w[j,u,r,a]
                     for r in rr for u in rr),polynomial_ring.zero)
       for i,j,a,b in w}
    tensor={(i,j,k,a,b,c):sum((z[i,j,r,a]*w[r,k,b,c]+z[i,j,r,b]*w[r,k,a,c]
                             +z[i,j,r,c]*w[r,k,a,b] for r in rr),polynomial_ring.zero)
            for i,j,k,a,b,c in product(rr,repeat=6)}
    c4=sum((z[i,j,a,b]*z[a,b,i,j] for i,j,a,b in w),polynomial_ring.zero)/QQ(13440)
    c6=-sum((v*tensor[a,b,c,i,j,k] for (i,j,k,a,b,c),v in tensor.items()),
            polynomial_ring.zero)/QQ(1036800)
    return c4.as_expr(),c6.as_expr()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    z=s.Matrix(s.symbols('z0:4'));y=s.Matrix(s.symbols('y0:4'));t=s.Symbol('t')
    source=s.zeros(4);source[0,1]=z[0];source[1,0]=-z[0]
    source[2,3]=1;source[3,2]=-1
    phi=s.Matrix([z[0]+z[1]**2+z[2]**2+z[3]**2,z[1],z[2],z[3]])
    inverse=[y[0]-y[1]**2-y[2]**2-y[3]**2,y[1],y[2],y[3]]
    assert phi.subs(dict(zip(z,inverse)),simultaneous=True)==y
    jacobian=phi.jacobian(z);assert jacobian.det()==1
    affine=(jacobian*source*jacobian.T).subs(dict(zip(z,inverse)),simultaneous=True)
    affine=affine.applyfunc(s.expand)
    assert affine[0,1]==y[0]-y[1]**2-y[2]**2-y[3]**2
    assert affine[0,2]==-2*y[3] and affine[0,3]==2*y[2] and affine[2,3]==1
    homogeneous=s.zeros(5)
    for i in range(4):
        for j in range(4):
            homogeneous[i+1,j+1]=s.expand(X[0]**2*affine[i,j].subs(
                dict(zip(y,[X[k+1]/X[0] for k in range(4)])),simultaneous=True))
    assert all(value==0 for value in jacobi_projective(homogeneous))

    point={y[0]:t,y[1]:0,y[2]:0,y[3]:0}
    constant=affine.subs(point);frame=s.diag(1,1/t,1,1);frame_inverse=s.diag(1,t,1,1)
    standard=s.diag(s.Matrix([[0,1],[-1,0]]),s.Matrix([[0,1],[-1,0]]))
    assert frame*constant*frame.T==standard
    derivatives=[affine.diff(v).subs(point) for v in y]
    linear=[frame*sum((frame_inverse[h,j]*derivatives[h] for h in range(4)),s.zeros(4))*frame.T
            for j in range(4)]
    beta=primitive_vector(linear).applyfunc(s.cancel)

    # Scale only the last row/column, by 5t, to work over QQ[t].
    variables=s.symbols('r0:5');a,b,c,d,e=variables;pl=[a,b,c,d,e,-a]
    star=s.Matrix([[0,-a,-e,d],[a,0,c,-b],[e,-c,0,a],[-d,b,-a,0]])
    ell=s.Matrix([sum(5*t*beta[4*h+i]*pl[h] for h in range(6)) for i in range(4)])
    model_phi=star.row_join(ell).col_join((-ell.T).row_join(s.zeros(1)))
    omega,_=pfaffOmega(model_phi,variables)
    scaled4,scaled6=univariate_invariants(omega,variables,t)
    c4=s.cancel(scaled4/(5*t)**8);c6=s.cancel(scaled6/(5*t)**12)
    assert c4==256/t**4 and c6==4096/t**6
    assert s.cancel(c4**3-c6**2)==0
    assert s.diff(c4,t)==-1024/t**5 and s.diff(c6,t)==-24576/t**7

    # Independently place the entire one-parameter family in the exact S61 model.
    projective_frame=s.eye(5);projective_frame[1:5,1:5]=frame
    projective_frame[1,0]=-t
    adapter=S61();normalized=adapter.normalize(pushforward(homogeneous,projective_frame))
    verification=adapter.verify_point(normalized)
    report={
        'id':'CURVE_INVARIANT_CONSERVATION_COUNTEREXAMPLE_20261101',
        'all_passed':True,'field':'QQ(t)','open_condition':'t != 0',
        'polynomial_automorphism':list(map(str,phi)),
        'affine_poisson_matrix':[[str(v) for v in row] for row in affine.tolist()],
        'homogeneous_poisson_matrix':[[str(v) for v in row] for row in homogeneous.tolist()],
        'moving_point':['t','0','0','0'],
        'primitive_beta_coefficients':list(map(str,beta)),
        'c4':str(c4),'c6':str(c6),'Delta':'0',
        'derivative_c4':str(s.diff(c4,t)),'derivative_c6':str(s.diff(c6,t)),
        'parameters_u1_through_u61':list(map(str,normalized['parameters'])),
        'perspective':list(map(str,normalized['perspective'])),
        'verification':verification,
        'scope':[
            'A global quadratic projective Poisson tensor, verified in full S61 over QQ(t).',
            'Refutes unrestricted c4/c6 basepoint constancy and unrestricted full-ideal membership.',
            'The canonical curve has Delta=0; no counterexample on the smooth-curve open is claimed.',
            'No new component, component identification, or generic-fiber conclusion is asserted.'],
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (Path(__file__),ROOT/'code/elliptic_curve_invariants_20261101.py',
                      ROOT/'code/elliptic_weight1_dominance.py',ROOT/'code/s61_model.py')}
    }
    with args.output.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps({'all_passed':True,'c4':str(c4),'c6':str(c6),'Delta':0,
                      'full_S61_verified_over_QQ_t':True}))


if __name__=='__main__':main()

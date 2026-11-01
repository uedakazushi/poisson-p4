#!/usr/bin/env python3
"""Global quadratic Poisson counterexample even on the smooth-curve open.

The canonical c4,c6 and j vary with the base point while Delta is nonzero.
This does not assert that this component dominates the full first-jet space.
"""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import sympy as s

from curve_invariant_conservation_counterexample_20261101 import univariate_invariants
from elliptic_curve_invariants_20261101 import pfaffOmega,curve_invariants_from_beta
from elliptic_weight1_dominance import projective_affine,symplectic_frame,primitive_vector
from s61_model import S61,X,jacobi_projective,pushforward

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    bivector=s.zeros(5);bivector[0,1]=-2*X[2]**2;bivector[1,0]=-bivector[0,1]
    for (i,j),a in zip(combinations(range(1,5),2),[2,8,12,3,5,7]):
        bivector[i,j]=a*X[i]*X[j];bivector[j,i]=-bivector[i,j]
    assert all(v==0 for v in jacobi_projective(bivector))
    y=s.symbols('y0:4');t=s.Symbol('t');affine=projective_affine(bivector,X,y)
    point=dict(zip(y,[t,1,1,1]));constant=affine.subs(point)
    frame=symplectic_frame(constant);inverse=frame.inv()
    derivatives=[affine.diff(v).subs(point) for v in y]
    linear=[frame*sum((inverse[h,j]*derivatives[h] for h in range(4)),s.zeros(4))*frame.T
            for j in range(4)]
    beta=primitive_vector(linear).applyfunc(s.cancel)
    scale=s.lcm([s.denom(v) for v in beta])
    r=s.symbols('r0:5');a,b,c,d,e=r;pl=[a,b,c,d,e,-a]
    star=s.Matrix([[0,-a,-e,d],[a,0,c,-b],[e,-c,0,a],[-d,b,-a,0]])
    ell=s.Matrix([sum(s.cancel(scale*beta[4*h+i])*pl[h] for h in range(6)) for i in range(4)])
    phi=star.row_join(ell).col_join((-ell.T).row_join(s.zeros(1)))
    omega,_=pfaffOmega(phi,r);scaled4,scaled6=univariate_invariants(omega,r,t)
    c4=s.factor(scaled4/scale**8);c6=s.factor(scaled6/scale**12)
    delta=s.factor((c4**3-c6**2)/1728);j=s.factor(c4**3/delta)
    assert delta!=0 and s.diff(c4,t)!=0 and s.diff(c6,t)!=0 and s.diff(j,t)!=0
    samples=[];adapter=S61()
    for value in [2,3]:
        actual=curve_invariants_from_beta(beta.subs(t,value))
        expected=tuple(s.cancel(v.subs(t,value)) for v in [c4,c6,delta,j])
        assert actual==expected and actual[2]!=0
        projective_frame=s.eye(5);sample_frame=frame.subs(t,value)
        projective_frame[1:5,1:5]=sample_frame
        projective_frame[1:5,0]=-sample_frame*s.Matrix([value,1,1,1])
        normalized=adapter.normalize(pushforward(bivector,projective_frame))
        verification=adapter.verify_point(normalized)
        samples.append({'point':[value,1,1,1],'c4':str(actual[0]),'c6':str(actual[1]),
                        'Delta':str(actual[2]),'j':str(actual[3]),
                        'parameters_u1_through_u61':list(map(str,normalized['parameters'])),
                        'perspective':list(map(str,normalized['perspective'])),
                        'verification':verification})
    assert samples[0]['j']!=samples[1]['j']

    report={
        'id':'CURVE_INVARIANT_SMOOTH_COUNTEREXAMPLE_20261101','all_passed':True,'field':'QQ(t)',
        'homogeneous_poisson_matrix':[[str(v) for v in row] for row in bivector.tolist()],
        'affine_poisson_matrix':[[str(v) for v in row] for row in affine.tolist()],
        'moving_point':['t','1','1','1'],'canonical_matrix_denominator_scale':str(s.factor(scale)),
        'primitive_beta_coefficients':list(map(str,beta)),
        'c4':str(c4),'c6':str(c6),'Delta':str(delta),'j':str(j),
        'derivative_c4':str(s.factor(s.diff(c4,t))),
        'derivative_c6':str(s.factor(s.diff(c6,t))),
        'derivative_j':str(s.factor(s.diff(j,t))),
        'exact_rational_samples':samples,
        'verification_scope':'Global projective Jacobi and invariant formulas exact over QQ(t); full S61 reconstruction and134equations checked separately at both rational points.',
        'scope':[
            'Refutes c4,c6 and j conservation on the full Delta-nonzero S61 locus.',
            'Thus the proposed transport equations are not in radical(I_S61:Delta^infinity).',
            'Does not refute conservation over a possibly smaller generic first-jet open.',
            'No new component or generic-fiber exhaustiveness claim is made.'],
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (Path(__file__),ROOT/'code/curve_invariant_conservation_counterexample_20261101.py',
                      ROOT/'code/elliptic_curve_invariants_20261101.py',
                      ROOT/'code/elliptic_weight1_dominance.py',ROOT/'code/s61_model.py')}
    }
    with args.output.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps({'all_passed':True,'c4':str(c4),'Delta':str(delta),
                      'j_varies_on_Delta_nonzero_locus':True,'full_S61_at_two_rational_points':True}))


if __name__=='__main__':main()

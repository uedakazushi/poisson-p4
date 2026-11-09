#!/usr/bin/env python3
"""Exact full-magnetic normal forms for the remaining nil22 image sources.

No point sampling or jets. Fiber affine translations are solved as a complete
six-by-six action, with inverse maps and all tensor coefficients retained.
"""
import argparse
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time

import sympy as s
from classification_nil22_image_source_review import build_P,coords,p,F0,mag,x,y,q1,q2

ROOT=Path(__file__).resolve().parents[1]


def normal_form(index,S,V0,V1):
    P,M,D=build_P(x*x*y,S,V0,V1)
    monomials=[s.Integer(1),x,y,x*x,x*y,y*y]
    alpha=s.diff(P[0,1],q1);beta=s.diff(P[0,1],q2)
    def deriv(field,value):return sum(field[j]*s.diff(value,p[j])for j in range(2))
    effects=[s.expand(-alpha*z-deriv(M[1,:],z))for z in [1,x,y]]+[s.expand(-beta*z+deriv(M[0,:],z))for z in [1,x,y]]
    T=s.Matrix([[s.Poly(f,x,y).coeff_monomial(mon)for f in effects]for mon in monomials])
    invariants=s.Matrix.vstack(*[a.T for a in T.T.nullspace()])
    assert invariants*T==s.zeros(invariants.rows,6)
    pivots=list(invariants.rref()[1]);assert len(pivots)==invariants.rows
    R=invariants[:,pivots].inv()*invariants
    canonical=s.zeros(6,1)
    for j,position in enumerate(pivots):canonical[position]=(R*s.Matrix(mag))[j]
    solution,free=T.gauss_jordan_solve(canonical-s.Matrix(mag))
    solution=solution.subs({a:0 for a in free})
    assert (T*solution-canonical+s.Matrix(mag)).applyfunc(s.expand)==s.zeros(6,1)
    L1=sum(solution[j]*mon for j,mon in enumerate([1,x,y]))
    L2=sum(solution[j+3]*mon for j,mon in enumerate([1,x,y]))
    newcoordinates=s.Matrix([q1+L1,q2+L2,x,y]);G=newcoordinates.jacobian(s.Matrix(coords))
    assert G.det()==1
    transported=(G*P*G.T).subs({q1:q1-L1,q2:q2-L2},simultaneous=True).applyfunc(s.expand)
    canonical_F=sum(canonical[j]*mon for j,mon in enumerate(monomials))
    expected=P.copy();expected[0,1]+=canonical_F-F0;expected[1,0]=-expected[0,1]
    assert (transported-expected).applyfunc(s.expand)==s.zeros(4)
    record=dict(source_case_index=index,full_magnetic_action_matrix=[[str(a)for a in row]for row in T.tolist()],
        translation_rank=T.rank(),complete_translation_kernel=[[str(a)for a in col]for col in T.nullspace()],
        full_magnetic_invariant_map=[[str(a)for a in row]for row in R.tolist()],
        residual_monomials=list(map(str,[monomials[j]for j in pivots])),residual_invariants=list(map(str,R*s.Matrix(mag))),
        exact_fiber_translation=[str(s.factor(L1)),str(s.factor(L2))],
        all_six_brackets_and_all90_coefficients_conjugacy_checked=True,translation_determinant='1',
        no_magnetic_coefficient_inverted=True,canonical_residual_F0=str(s.factor(canonical_F)),
        actual_affine_determinant=str(s.factor(D)),ambient_component_membership_claimed=False)
    if index==0:
        assert R.rows==1 and pivots==[1]
        invariant=s.expand((R*s.Matrix(mag))[0])
        assert invariant==mag[1]-mag[3]-4*mag[4]/3
        I=s.Symbol('I',nonzero=True);scale=1/I
        canonical_P=P.subs(dict(zip(mag,[0,I,0,0,0,0])))
        scaling=s.diag(scale,scale,1,1)
        actual=(scaling*canonical_P*scaling.T).subs({q1:q1/scale,q2:q2/scale},simultaneous=True)/scale
        target=P.subs(dict(zip(mag,[0,1,0,0,0,0])))
        assert (actual-target).applyfunc(s.cancel)==s.zeros(4)
        tau=s.Symbol('tau');curve=P.subs(dict(zip(mag,[0,tau,0,0,0,0])))
        for i,j,k in it.combinations(range(4),3):
            assert s.expand(sum(curve[i,l]*s.diff(curve[j,k],coords[l])+curve[j,l]*s.diff(curve[k,i],coords[l])+curve[k,l]*s.diff(curve[i,j],coords[l])for l in range(4)))==0
        record.update(entire_nonzero_invariant_source_in_one_actual_projective_orbit=True,
            invariant_nonzero_normalizer='q_new=q/I, tensor_scale=I',
            zero_invariant_complement='All zero-invariant points translate exactly to F0=0; P_tau with F0=tau*x is an actual polynomial Poisson curve from the nonzero orbit to it.',
            entire_six_parameter_source_in_single_actual_orbit_closure=True,
            all_points_in_single_actual_orbit_claimed=False,
            generic_canonical_F0='x',boundary_canonical_F0='0',
            generic_point_torus_from_q_scaling_available=False,
            boundary_point_q_scaling_torus_available=True)
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    begun=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(110,115));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    specifications=[(0,s.Matrix([[s.Rational(8,3),1],[1,-s.Rational(3,64)]]),s.Matrix([-s.Rational(1,8),-s.Rational(5,3)]),s.diag(-s.Rational(1,4),-s.Rational(1,2))),
        (1,s.diag(-8,1),s.Matrix([0,1]),s.diag(-s.Rational(1,4),-s.Rational(1,2))),
        (2,s.Matrix([[2,1],[1,0]]),s.Matrix([0,-2]),s.diag(0,-1))]
    records=[normal_form(*entry)for entry in specifications]
    inputs=[Path(__file__),ROOT/'code/classification_nil22_image_source_review.py',ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_SOURCE_REVIEW_V2_20261107.json']
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,no_producer_import=True,
        input_sha256={str(z.relative_to(ROOT)):hashlib.sha256(z.read_bytes()).hexdigest()for z in inputs},cases=records,
        implication='double_A entire6source is reduced by explicit affine/projective maps to its F0=x actualorbit and the F0=0 boundary. This reduces recognition to one actual orbit, without proving any ambient component membership.',
        reviewer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps(dict(all_passed=True,translation_ranks=[a['translation_rank']for a in records],seconds=result['seconds'])))


if __name__=='__main__':main()

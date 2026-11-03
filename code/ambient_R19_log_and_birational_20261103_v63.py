#!/usr/bin/env python3
"""Exact Log1112 containment, rational Casimirs, and quadratic Poisson map."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,canonical,vector

ROOT=Path(__file__).resolve().parents[1]

def sign(seq):return (-1)**sum(seq[i]>seq[j] for i in range(len(seq)) for j in range(i+1,len(seq)))

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic();h,x,y,z,w=X;A,B,F=s.symbols('A B F')
    Q=A*h*h+B*h*z+F*z*z;V=s.Matrix([0,4*x*x,0,-2*x*z,2*x*w+Q]);ey=s.Matrix([0,0,1,0,0]);P=ey*V.T-V*ey.T;E=s.Matrix(X)
    K=120*x*w+20*A*h*h+15*B*h*z+12*F*z*z
    alpha=s.Matrix([-3/h,1/x,0,2/z,0]);beta=s.Matrix([s.diff(K,u)/K for u in X])+s.Matrix([-2/h,-1/x,0,1/z,0])
    product=h*x*z*K;omega={};dual={}
    for i,j in itertools.combinations(range(5),2):
        comp=[a for a in range(5) if a not in [i,j]];a,b,c=comp
        dual[i,j]=sign(comp+[i,j])*(E[a]*P[b,c]-E[b]*P[a,c]+E[c]*P[a,b])
        omega[i,j]=s.cancel(product*(alpha[i]*beta[j]-alpha[j]*beta[i]))
        assert s.expand(omega[i,j]+60*dual[i,j])==0
    dform={}
    for a,b,c in itertools.combinations(range(5),3):dform[a,b,c]=s.expand(s.diff(omega[b,c],X[a])-s.diff(omega[a,c],X[b])+s.diff(omega[a,b],X[c]))
    logP=s.zeros(5)
    for i,j in itertools.combinations(range(5),2):
        comp=tuple(a for a in range(5) if a not in [i,j]);logP[i,j]=sign([i,j,*comp])*dform[comp];logP[j,i]=-logP[i,j]
    assert (logP+300*canonical(P)).applyfunc(s.expand)==s.zeros(5)
    ra=s.Matrix([-3,1,2,0]);rb=s.Matrix([-2,-1,1,1]);degree=s.Matrix([1,1,1,2]);assert ra.dot(degree)==rb.dot(degree)==0
    residue=ra*rb.T-rb*ra.T;old_a,old_b,old_c=s.Rational(5,4),s.Rational(1,4),s.Rational(3,4)
    expected=[4*old_a,4*old_b,-2*(old_a+old_b),4*old_c,2*(old_a-old_c),2*(old_b+old_c)]
    assert [residue[i,j] for i,j in itertools.combinations(range(4),2)]==expected
    phi=x*z*z/h**3;psi=z*K/(120*x*h*h)
    assert s.cancel(sum(V[i]*s.diff(phi,X[i]) for i in range(5)))==0
    assert s.cancel(sum(V[i]*s.diff(psi,X[i]) for i in range(5)))==0
    neww=K/(120*x);change=s.Matrix([h,x,y,z,neww]);jacobian=change.jacobian(X)
    assert jacobian.det()==1 and (jacobian*E-change).applyfunc(s.cancel)==s.zeros(5,1)
    pushed=(jacobian*P*jacobian.T).subs(w,w-A*h*h/(6*x)-B*h*z/(8*x)-F*z*z/(10*x),simultaneous=True).applyfunc(s.cancel)
    Vzero=s.Matrix([0,4*x*x,0,-2*x*z,2*x*w]);assert pushed==ey*Vzero.T-Vzero*ey.T
    r,t,u,Y=s.symbols('r t u Y');inverse={h:1,x:r/(u*u),y:Y,z:u,w:t/u-A*u*u/(6*r)-B*u**3/(8*r)-F*u**4/(10*r)}
    assert s.cancel(phi.subs(inverse,simultaneous=True))==r
    assert s.cancel(psi.subs(inverse,simultaneous=True))==t
    old_path=ROOT/'legacy/v14/data/six_rank2_components.json';old=json.loads(old_path.read_text());cases=old if isinstance(old,list) else old['cases'] if 'cases' in old else old['components'];logold=next(v for v in cases if v.get('name')=='Log1112')
    assert logold['projective_component_dimension']==28 and logold['ambient_affine_tangent_dimension']==logold['family_affine_differential_rank']==29
    paths=[Path(__file__),ROOT/'legacy/v14/code/rank2_lower_bounds/verify_rank2_components.py',old_path,
        ROOT/'research/notes/AMBIENT_R19_COMPLEMENT_ORBITS_20261103_v63.md']
    body=dict(id='AMBIENT_R19_LOG_AND_BIRATIONAL_20261103_v63',all_passed=True,field='QQ(A,B,F)',
        quadratic_divisor=str(K),logarithmic_divisors=['h','x','z',str(K)],residue_vectors=[list(map(int,ra)),list(map(int,rb))],
        legacy_Log1112_residue_parameters=['5/4','1/4','3/4'],log_2form_to_Euler_dual_P_factor='-60',
        legacy_domega_dual_canonical_P_factor='-300',all10_form_and150_canonical_coefficients_exact=True,
        entire_R_containment='R19 subset existing Log1112_28 component by exactnormalform identity and closedPGLsaturation',
        existing_Log1112_component_certificate=logold,
        rational_Casimirs=[str(phi),str(psi)],birational_Poisson_coordinate_change=list(map(str,change)),
        homogeneous_quadratic_Cremona_coordinates=['120*h*x','120*x²','120*x*y','120*x*z',str(K)],
        actual_birational_Poisson_push_equals_Qzero=True,rational_generic_leaf_inverse={str(k):str(v) for k,v in inverse.items()},
        source_scope='EveryR minusN point hasoneoftheseVQforms; birationalPoissonequivalence andrationalleaves asserted there. ClosedLogcontainment extendsENTIRER.',
        known_full_Poisson_components_through_entire_R=['CA28','CB28','Log1112_28','LPB39'],
        old_component_certificate_reused_not_rerun=True,component_count_change=0,
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(body,indent=2)+'\n');print('ENTIRE R LOG + BIRATIONAL POISSON PASS',body['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)

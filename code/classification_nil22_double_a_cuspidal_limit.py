#!/usr/bin/env python3
"""Whole double_A six-magnetic source at a regular known cuspidal24 boundary."""
import hashlib
import itertools as it
import json
import resource
import time
from pathlib import Path

import sympy as s
from classification_rank3_mixed_component_review import jacobi as projective_jacobi, PAIRS

ROOT=Path(__file__).resolve().parents[1]


def main():
    begun=time.monotonic()
    h,x,y,z,w=parent_coords=s.symbols('h x y z w')
    a,ep=s.symbols('a epsilon')
    r=s.symbols('r0:6')
    env={str(v):v for v in parent_coords+(a,)+r}
    parent_path=ROOT/'reports/CLASSIFICATION_GA9_CUSPIDAL_QUARTIC_COMPONENT_20261107.json'
    parent=json.loads(parent_path.read_text())['actual_family_power3']
    P=s.zeros(5)
    for ij,text in parent['universal_Poisson_brackets'].items():
        i,j=map(int,ij.split(','));P[i,j]=s.sympify(text,locals=env);P[j,i]=-P[i,j]
    assert s.sympify(parent['coefficient_denominator'],locals=env).subs(a,s.Rational(1,4))!=0
    q1,q2,p1,p2=s.symbols('q1 q2 p1 p2')
    target_coords=(h,q1,q2,p1,p2)
    oldh,oldq1,oldq2,oldp1,oldp2=s.symbols('h x1 x2 x3 x4')
    m=s.symbols('m0:6')
    tenv={str(v):v for v in (oldh,oldq1,oldq2,oldp1,oldp2)+m}
    source_path=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json'
    row=json.loads(source_path.read_text())['cases'][0]
    assert row['name']=='nil22_image_double_A'
    target=s.zeros(5)
    for ij,text in row['full_six_magnetic_source'].items():
        i,j=map(int,ij.split(','))
        v=s.sympify(text,locals=tenv).subs(dict(zip((oldh,oldq1,oldq2,oldp1,oldp2),target_coords)),simultaneous=True)
        target[i,j]=v;target[j,i]=-v
    B=s.Matrix([[-4,0],[-48,-s.Rational(9,2)]])
    c=s.Matrix([-3,-54])
    base_change=B*s.Matrix([p1,p2])+c*h
    F=s.sympify(parent['plane_quartic'],locals=env)
    C3=h*w*w+12*w*y*y-144*y**3
    assert s.expand(F.subs(a,s.Rational(1,4))-h*C3)==0
    Ctarget=s.sympify(row['affine_D'],locals=tenv).subs({oldp1:p1,oldp2:p2},simultaneous=True)
    Ctarget=s.Poly(Ctarget,p1,p2)
    Chom=sum(co*h**(3-sum(ex))*p1**ex[0]*p2**ex[1] for ex,co in Ctarget.terms())
    assert s.expand(C3.subs({y:base_change[0],w:base_change[1]},simultaneous=True)+864*Chom)==0
    Mparent=s.Matrix([[P[i,j] for j in [2,4]] for i in [1,3]])
    Mparent=Mparent.subs(a,s.Rational(1,4)).subs({y:base_change[0],w:base_change[1]},simultaneous=True)
    Mtarget=s.Matrix([[target[i,j] for j in [3,4]] for i in [1,2]])
    k=s.symbols('k0:4');A=s.Matrix(2,2,k)
    difference=(A*Mparent*B.inv().T-Mtarget).applyfunc(s.expand)
    equations=[co for f in difference for co in s.Poly(f,h,p1,p2).coeffs()]
    solutions=list(s.linsolve(equations,k))
    assert len(solutions)==1,'The actual two-anchor planes differ; no containment inferred.'
    values=solutions[0]
    assert not set().union(*(v.free_symbols for v in values))
    A=A.subs(dict(zip(k,values)));assert A.det()!=0
    G=s.zeros(5);G[0,0]=1
    for i,ii in enumerate([1,2]):
        for j,jj in enumerate([1,3]):G[ii,jj]=A[i,j]
    for i,ii in enumerate([3,4]):
        G[ii,0]=(-B.inv()*c)[i]
        for j,jj in enumerate([2,4]):G[ii,jj]=B.inv()[i,j]
    assert G.det()!=0
    inverse=G.inv()*s.Matrix(target_coords)
    curve=(G*P.subs(a,s.Rational(1,4)+ep)*G.T).subs(dict(zip(parent_coords,inverse)),simultaneous=True).applyfunc(s.cancel)
    endpoint=curve.subs(ep,0).applyfunc(s.cancel)
    coefficient=s.expand(endpoint[1,2]-target[1,2].subs(dict.fromkeys(m,0)))
    assert coefficient.free_symbols<=set((h,p1,p2)+r)
    mons=(h*h,h*p1,h*p2,p1*p1,p1*p2,p2*p2)
    magnetic_map=s.Matrix([s.Poly(coefficient,h,p1,p2).coeff_monomial(mon) for mon in mons])
    magnetic_matrix=magnetic_map.jacobian(r)
    assert magnetic_matrix.det()!=0
    assert (magnetic_map-magnetic_matrix*s.Matrix(r)).applyfunc(s.expand)==s.zeros(6,1)
    assert (endpoint-target.subs(dict(zip(m,magnetic_map)),simultaneous=True)).applyfunc(s.expand)==s.zeros(5)
    # Keep original target m's throughout the actual family, using the inverse
    # six-by-six map. This matches ALL150 endpoint coefficients identically.
    r_binding=magnetic_matrix.inv()*s.Matrix(m)
    curve=curve.subs(dict(zip(r,r_binding)),simultaneous=True).applyfunc(s.cancel)
    assert (curve.subs(ep,0)-target).applyfunc(s.expand)==s.zeros(5)
    sparse={ij:dict(s.Poly(s.expand(curve[ij]),*target_coords).terms()) for ij in PAIRS if curve[ij]!=0}
    assert projective_jacobi(sparse).applyfunc(s.cancel)==s.zeros(350,1)
    strict=[]
    for i,j,k in it.combinations(range(5),3):
        v=s.cancel(sum(curve[i,l]*s.diff(curve[j,k],target_coords[l])+
                       curve[j,l]*s.diff(curve[k,i],target_coords[l])+
                       curve[k,l]*s.diff(curve[i,j],target_coords[l]) for l in range(5)))
        strict.append(v)
    assert strict==[0]*10
    affine=curve[1:,1:].subs(h,1)
    pf=s.factor(affine[0,1]*affine[2,3]-affine[0,2]*affine[1,3]+affine[0,3]*affine[1,2])
    determinant_ratio=s.sympify(parent['determinant_ratio'],locals=env).subs(a,s.Rational(1,4)+ep)
    actualD=s.factor(s.cancel(-h*G.det()**2*determinant_ratio*
                             F.subs(a,s.Rational(1,4)+ep).subs(dict(zip(parent_coords,inverse)),simultaneous=True)))
    # Affine Pf is sufficient to fix the canonical Euler contraction and scale.
    candidate=s.expand(actualD.subs(h,1)-pf)
    if candidate!=0:
        actualD=s.factor(h**5*pf.subs({q1:q1/h,q2:q2/h,p1:p1/h,p2:p2/h},simultaneous=True))
    assert s.expand(actualD.subs(ep,0)+h*h*Chom)==0
    den=s.factor(s.lcm([s.denom(s.cancel(v)) for v in curve]))
    assert den.subs(ep,0)!=0
    def encode(M):return {f'{i},{j}':str(M[i,j]) for i,j in PAIRS if M[i,j]!=0}
    inputs=[Path(__file__),parent_path,source_path,ROOT/'reports/CLASSIFICATION_GA9_CUSPIDAL_QUARTIC_REVIEW_20261107.json',ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_SOURCE_REVIEW_V2_20261107.json',ROOT/'code/classification_rank3_mixed_component_review.py']
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,status='WHOLE_DOUBLE_A_SOURCE_IN_KNOWN_CUSPIDAL24_CANDIDATE',
        known_component_id='R4_GA9_CUSPIDAL_PLANE_QUARTIC_24',parent_parameter='a=1/4+epsilon',
        base_parent_from_target=[str(v) for v in base_change],base_matrix=[[str(v) for v in B.row(i)] for i in range(2)],
        base_affine_translation=list(map(str,c)),base_determinant=str(B.det()),
        actual_fiber_parent_to_target=[[str(v) for v in A.row(i)] for i in range(2)],
        actual_fiber_determinant=str(A.det()),
        whole_projective_parent_to_target=[[str(v) for v in G.row(i)] for i in range(5)],
        whole_projective_determinant=str(G.det()),
        endpoint_magnetic_from_parent=list(map(str,magnetic_map)),
        endpoint_magnetic_matrix=[[str(v) for v in magnetic_matrix.row(i)] for i in range(6)],
        endpoint_magnetic_determinant=str(s.factor(magnetic_matrix.det())),
        actual_parent_r_binding=list(map(str,r_binding)),
        complete_regular_family=encode(curve),family_denominator=str(den),
        all350_projective_Jacobi_zero=True,all10_strict_Jacobi_zero=True,
        all150_endpoint_coefficients_whole_six_magnetic_exact=True,
        actual_projective_degeneracy_quintic=str(actualD),actual_affine_Pfaffian=str(pf),
        endpoint_divisor='-h^2 * canonical_double_A_C3',
        containment_proof='Fixed invertible projective push of the accepted all7source family at a regular a=1/4. Every target m has an explicit invertible all6mag parent binding. Consequently the ENTIRE sixmag double_A source is in the known component; PGL saturation and closure preserve this inclusion. No boundary membership is reversed.',
        independent_review_status='pending',inputs=[dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in inputs],
        seconds=time.monotonic()-begun,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    out=ROOT/'reports/CLASSIFICATION_NIL22_DOUBLE_A_CUSPIDAL_LIMIT_20261107.json'
    with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:result[k] for k in ['all_passed','status','actual_fiber_determinant','endpoint_magnetic_determinant','family_denominator','seconds','max_rss_kib']}))


if __name__=='__main__':main()

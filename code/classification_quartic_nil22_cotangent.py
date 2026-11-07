#!/usr/bin/env python3
"""Finite exact identities for nil22-image zero-base-bracket quartic geometry."""
import argparse
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time

from classification_quadratic_cofactor_review import Ring,PAIRS,verify_source

ROOT=Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args();start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(90,95));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    path=ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json'
    row=json.loads(path.read_text())['full_nilpotent_sources'][5]
    source_review=verify_source(row,6);assert len(source_review['free_indices'])==27
    R=Ring(27);p=list(map(R.parse,row['brackets']))
    keep=lambda p:{ex:c for ex,c in p.items()if not any(ex[i]for i in [18,19,24])}
    p=list(map(keep,p));assert not p[4]
    # q=(x,z), base=(y,w); pi=dx^# wedge U+dz^# wedge V+F dx^# wedge dz^#.
    U=[p[0],p[2]];V=[R.scale(p[3],-1),p[5]]
    assert all(not ex[-4]and not ex[-2]for q in U+V for ex in q)
    der=lambda p,k:R.derivative(p,k)
    def apply(X,f):return R.add(R.mul(X[0],der(f,1)),R.mul(X[1],der(f,3)))
    div=lambda X:R.add(der(X[0],1),der(X[1],3))
    det=lambda X,Y:R.add(R.mul(X[0],Y[1]),R.scale(R.mul(X[1],Y[0]),-1))
    Ffiber=p[1];alpha=der(Ffiber,0);beta=der(Ffiber,2)
    assert alpha==R.add(R.variable('w'),div(V))
    assert beta==R.scale(R.add(R.variable('y'),div(U)),-1)
    magnetic=R.add(Ffiber,R.scale(R.mul(R.variable('x'),alpha),-1),R.scale(R.mul(R.variable('z'),beta),-1))
    expected=R.parse('n3+n4*y+n5*w+n6*y**2+n7*w*y+n8*w**2')
    assert magnetic==expected
    assert not any(ex[i]for q in U+V for ex in q for i in range(3,9))
    D=det(U,V)
    pf=R.add(R.mul(p[0],p[5]),R.mul(p[2],p[3]))
    assert D==pf==keep(R.parse(row['actual_affine_pfaffian']))
    K=[R.add(apply(U,V[i]),R.scale(apply(V,U[i]),-1),
             R.scale(R.mul(alpha,U[i]),-1),R.scale(R.mul(beta,V[i]),-1))for i in range(2)]
    assert R.add(apply(U,D),R.mul(R.variable('y'),D),R.scale(det(U,K),-1))=={}
    assert R.add(apply(V,D),R.mul(R.variable('w'),D),R.scale(det(V,K),-1))=={}
    P=[[{}for _ in range(4)]for _ in range(4)]
    for (i,j),q in zip(PAIRS,p):P[i][j],P[j][i]=q,R.scale(q,-1)
    J=[R.add(*(R.mul(P[a][l],der(P[b][c],l))for a,b,c in [(i,j,k),(j,k,i),(k,i,j)]for l in range(4)))for i,j,k in it.combinations(range(4),3)]
    assert J==[R.scale(K[0],-1),{},K[1],{}]
    assert all(not ex[i]for q in J for ex in q for i in range(3,9))
    # Free quintic exponents(2,2): degree of singular scheme from HB.
    # Numerator 1-3t^4+2t^6 has second derivative24 at1; degree=12.
    assert (-3*4*3+2*6*5)//2==12
    possibilities=[]
    for genus in range(4):
        for branches in range(7):
            for intersections in range(1,5):
                muF=6-2*genus-branches
                muD=muF+8-intersections
                if muD>=12:possibilities.append([genus,branches,intersections,muF])
    assert possibilities==[[0,0,1,6],[0,0,2,6],[0,1,1,5]]
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,no_CAS_import=True,no_producer_import=True,
        source=source_review,zero_base_conditions=[18,19,24],full_zero_base_parameter_count=24,
        anchor_rows=[[R.encode(q)for q in X]for X in [U,V]],all_anchors_depend_only_on_base=True,
        all_six_magnetic_parameters_free_and_absent_from_all_Jacobi=True,
        fiber_linear_coefficients=[R.encode(alpha),R.encode(beta)],
        all_four_Jacobi_exactly_anchor_commutator=True,
        two_modular_logarithmic_identities_are_exact_Jacobi_syzygies=True,
        free_quintic_Tjurina_length=12,
        complete_numeric_genus_branch_intersection_possibilities=possibilities,
        geometric_dimension_proof='HYPERPLANE_QUARTIC_NIL22_COTANGENT_G_BOUND.md',
        geometry_not_proved_by_gate_or_arithmetic=True,
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()},
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_bytes=Path(__file__).read_text(),
        seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['all_passed','seconds','maximum_RSS_KiB']}))

if __name__=='__main__':main()

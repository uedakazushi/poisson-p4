#!/usr/bin/env python3
"""Exact Fisher gradients and exploratory S61 basepoint-transport rows.

This computes candidate affine-linear q equations at one fixed beta. It DOES
NOT assert their membership in the full Poisson ideal or invariant constancy
on unknown components.
"""
import argparse,hashlib,json
from itertools import product
from pathlib import Path
import sympy as s
from elliptic_curve_invariants_20261101 import canonical_phi,pfaffOmega,curve_invariants_from_beta

ROOT=Path(__file__).resolve().parents[1]
RANGE=range(5)
KEY4=list(product(RANGE,repeat=4));KEY6=list(product(RANGE,repeat=6))

def omega_tensor(beta):
    x,A,den=canonical_phi(beta);assert den==1
    O,_=pfaffOmega(A,x)
    polys={(i,j):s.Poly(O[i,j],*x) for i in RANGE for j in RANGE}
    return {(i,j,a,b):int(polys[i,j].coeff_monomial(x[a]*x[b]))*(2 if a==b else 1)
            for i,j,a,b in KEY4}

def secondary(W):
    Z={(i,j,a,b):sum(W[i,r,u,a]*W[j,u,r,b]+W[i,r,u,b]*W[j,u,r,a]
                    for r in RANGE for u in RANGE) for i,j,a,b in KEY4}
    T={(i,j,k,a,b,c):sum(Z[i,j,r,a]*W[r,k,b,c]+Z[i,j,r,b]*W[r,k,a,c]+Z[i,j,r,c]*W[r,k,a,b]
                        for r in RANGE) for i,j,k,a,b,c in KEY6}
    return Z,T

def differential(W,Z,T,V):
    U={(i,j,a,b):sum(V[i,r,u,a]*W[j,u,r,b]+W[i,r,u,a]*V[j,u,r,b]
                     +V[i,r,u,b]*W[j,u,r,a]+W[i,r,u,b]*V[j,u,r,a]
                     for r in RANGE for u in RANGE) for i,j,a,b in KEY4}
    dc4=s.Rational(2*sum(U[i,j,a,b]*Z[a,b,i,j] for i,j,a,b in KEY4),13440)
    total=0
    for i,j,k,a,b,c in KEY6:
        dt=sum(U[i,j,r,a]*W[r,k,b,c]+Z[i,j,r,a]*V[r,k,b,c]
              +U[i,j,r,b]*W[r,k,a,c]+Z[i,j,r,b]*V[r,k,a,c]
              +U[i,j,r,c]*W[r,k,a,b]+Z[i,j,r,c]*V[r,k,a,b] for r in RANGE)
        total+=dt*T[a,b,c,i,j,k]
    return dc4,-s.Rational(2*total,1036800)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    source=ROOT/'reports/S61_JET_BASEPOINT_TRANSPORT_20261101.json'
    data=json.loads(source.read_text());basis=s.Matrix(data['primitive_coefficient_basis'])
    u=s.Matrix(data['fixed_u_values']);beta=basis*u
    scale=s.ilcm(*(s.denom(v) for v in list(beta)+list(basis)))
    ib=scale*beta;W=omega_tensor(ib);Z,T=secondary(W)
    c4,c6,disc,j=curve_invariants_from_beta(beta)
    assert s.Rational(sum(Z[i,j,a,b]*Z[a,b,i,j] for i,j,a,b in KEY4),13440)==c4*scale**8
    assert -s.Rational(sum(T[i,j,k,a,b,c]*T[a,b,c,i,j,k] for i,j,k,a,b,c in KEY6),1036800)==c6*scale**12
    grad4=[];grad6=[]
    for column in range(16):
        direction=scale*basis[:,column]
        plus=omega_tensor(ib+direction);minus=omega_tensor(ib-direction)
        # Omega is polynomial homogeneous of beta-degree2, so centraldifference
        # gives its exact derivative, not a finite-difference approximation.
        V={key:(plus[key]-minus[key])//2 for key in KEY4}
        assert all(2*V[key]==plus[key]-minus[key] for key in KEY4)
        d4,d6=differential(W,Z,T,V)
        grad4.append(d4/scale**8);grad6.append(d6/scale**12)
    g4=s.Matrix([grad4]);g6=s.Matrix([grad6]);gj=2*c4*g6-3*c6*g4
    assert (g4*u)[0]==8*c4 and (g6*u)[0]==12*c6
    # An independent direction check uses the exact degree12 interpolation.
    direction=basis*s.ones(16,1)
    points=list(range(-6,7));vals4=[];vals6=[]
    for t in points:
        a,b,_,_=curve_invariants_from_beta(beta+t*direction);vals4.append(a);vals6.append(b)
    weights=s.finite_diff_weights(1,points,0)[1][-1]
    assert sum(w*v for w,v in zip(weights,vals4))==sum(grad4)
    assert sum(w*v for w,v in zip(weights,vals6))==sum(grad6)
    rows4=[];rows6=[];rowsj=[]
    for entry in data['derivative_by_basepoint_direction']:
        matrix=s.Matrix(entry['q_coefficient_matrix']).row_join(s.Matrix(entry['constant']))
        rows4.append(g4*matrix);rows6.append(g6*matrix);rowsj.append(gj*matrix)
    m4=s.Matrix.vstack(*rows4);m6=s.Matrix.vstack(*rows6);mj=s.Matrix.vstack(*rowsj)
    together=m4.col_join(m6)
    known=json.loads((ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json').read_text())
    point=next(p for p in known['points'] if p['component_id']=='T01')['S61_coordinates_u1_through_u61']
    q=s.Matrix([s.Rational(point[i]) for i in data['q_to_u_indices_zero_based']]+[1])
    assert together*q==s.zeros(8,1) and mj*q==s.zeros(4,1)
    serial=lambda mat:[[str(v) for v in row] for row in mat.tolist()]
    report={'id':'PRIMITIVE_CURVE_INVARIANT_GRADIENTS_20261101','all_passed':True,'field':'QQ',
        'fixed_u_order':data['fixed_primitive_u_order'],'fixed_u_values':data['fixed_u_values'],
        'c4':str(c4),'c6':str(c6),'j':str(j),'gradient_c4':list(map(str,grad4)),
        'gradient_c6':list(map(str,grad6)),'gradient_j_numerator':list(map(str,gj)),
        'exact_degree12_direction_interpolation_checked':True,'euler_degree_identities_checked':True,
        'candidate_c4_transport_affine_rows':serial(m4),'candidate_c6_transport_affine_rows':serial(m6),
        'candidate_j_transport_affine_rows':serial(mj),
        'candidate_linear_ranks':{'c4':m4[:,:45].rank(),'c6':m6[:,:45].rank(),
                                  'combined':together[:,:45].rank(),'j':mj[:,:45].rank()},
        'candidate_augmented_ranks':{'combined':together.rank(),'j':mj.rank()},
        'known_toric_center_annihilated':True,'q_coordinate_convention':data['q_coordinates'],
        'scope':['Exploratory candidate equations at one beta, not proven Poisson ideal consequences.',
                 'Constancy of these invariants on arbitrary unknown Poisson structures remains unproved.',
                 'The toriccenter check and exact ranks do not justify global elimination or exhaustion.'],
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in (Path(__file__),source,ROOT/'code/elliptic_curve_invariants_20261101.py')}}
    with args.output.open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps({'all_passed':True,'candidate_ranks':report['candidate_linear_ranks'],'not_ideal_membership':True}))

if __name__=='__main__':main()

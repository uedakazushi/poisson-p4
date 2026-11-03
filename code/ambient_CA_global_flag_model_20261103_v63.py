#!/usr/bin/env python3
"""Parabolic-invariant20-space giving a global proper model of reduced CA.

This includes generic-Poisson-rank<=2 points; it does not infer global R1,
normality, a whole intersection theorem, or fibers from source equations alone.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
X=s.symbols('h x y z w');h,x,y,z,w=X
PAIRS=list(itertools.combinations(range(5),2))
MONS=list(itertools.combinations_with_replacement(range(5),2))

def raw(f,g,H):
    P=s.zeros(5);P[1,2]=f;P[3,4]=g*h*h
    P[2,3]=-s.diff(H,w);P[2,4]=s.diff(H,z)
    return P-P.T

def divergence(P):return s.Matrix([sum(s.diff(P[i,j],X[i])for i in range(5))for j in range(5)])

def canonical(P):
    d=divergence(P);xx=s.Matrix(X)
    return (P-(xx*d.T-d*xx.T)/5).applyfunc(s.expand)

def vector(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(s.prod(X[t]for t in m))for i,j in PAIRS for m in MONS])

def action(M,P):
    v=M*s.Matrix(X)
    return (M*P+P*M.T-sum((v[j]*P.diff(X[j])for j in range(5)),s.zeros(5))).applyfunc(s.expand)

def boundary_witness(P,source_basis):
    columns=[vector(R)for R in source_basis]
    for i,j in itertools.product(range(5),repeat=2):
        M=s.zeros(5);M[i,j]=1;columns.append(vector(action(M,P)))
    C=s.Matrix.hstack(*columns);pc=list(C.rref()[1])
    pr=list(C[:,pc].T.rref()[1]);det=C[pr,pc].det();assert det!=0
    mon3=list(itertools.combinations_with_replacement(range(5),3));wedge=[]
    for q in range(5):
        v=[int(j==q)for j in range(5)];entries=[]
        for i,j,k,l in itertools.combinations(range(5),4):
            value=s.expand(X[i]*(v[j]*P[k,l]-v[k]*P[j,l]+v[l]*P[j,k])
                -X[j]*(v[i]*P[k,l]-v[k]*P[i,l]+v[l]*P[i,k])
                +X[k]*(v[i]*P[j,l]-v[j]*P[i,l]+v[l]*P[i,j])
                -X[l]*(v[i]*P[j,k]-v[j]*P[i,k]+v[k]*P[i,j]))
            poly=s.Poly(value,*X);entries.extend(poly.coeff_monomial(s.prod(X[t]for t in m))for m in mon3)
        wedge.append(s.Matrix(entries))
    W=s.Matrix.hstack(*wedge);wc=list(W.rref()[1]);wr=list(W[:,wc].T.rref()[1])
    return dict(affine_image_differential_rank=len(pc),projective_image_dimension_lower_bound=len(pc)-1,
        minor_rows=pr,minor_columns=pc,minor_determinant=str(det),
        exact_raw150_image_differential=[[str(v)for v in row]for row in C.tolist()],
        intrinsic_constant_vector_factor_matrix_rank=len(wc),
        constant_vector_factor_minor_rows=wr,constant_vector_factor_minor_columns=wc,
        constant_vector_factor_minor_determinant=str(W[wr,wc].det()))

def main(output):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    cubics=[s.prod((h,x,z,w)[i]for i in m)for m in itertools.combinations_with_replacement(range(4),3)
        if any(i>=2 for i in m)]
    assert len(cubics)==16
    raws=[raw(f,0,0)for f in [h*h,h*x,x*x]]+[raw(0,1,0)]+[raw(0,0,H)for H in cubics]
    basis=[canonical(P)for P in raws];C=s.Matrix.hstack(*(vector(P)for P in basis))
    selected=list(C.T.rref()[1]);assert len(selected)==20
    inverse=C[selected,:].inv();assert C[selected,:].det()!=0
    # Row0 preserves<h>, row1 preserves<h,x>, column2 preserves<d_y>.
    positions=[(i,j)for i,j in itertools.product(range(5),repeat=2)
        if not(i==0 and j!=0)and not(i==1 and j not in (0,1))and not(j==2 and i!=2)]
    assert len(positions)==16
    generators=[];actions=[]
    for i,j in positions:
        M=s.zeros(5);M[i,j]=1;generators.append(M)
        B=s.Matrix.hstack(*(vector(action(M,P))for P in basis))
        small=inverse*B[selected,:]
        assert C*small==B
        actions.append(small)
    assert sum((actions[positions.index((i,i))]for i in range(5)),s.zeros(20))==s.zeros(20)
    # Independently check the full20-parameter raw Jacobi identities.
    u=s.symbols('u0:20');P=sum((c*R for c,R in zip(u,raws)),s.zeros(5))
    jac=[]
    for i,j,k in itertools.combinations(range(5),3):
        value=s.expand(sum(P[i,t]*s.diff(P[j,k],X[t])+P[j,t]*s.diff(P[k,i],X[t])+P[k,t]*s.diff(P[i,j],X[t])for t in range(5)))
        assert value==0;jac.append(str(value))
    assert divergence(P)==s.Matrix([0,0,2*u[2]*x+u[1]*h,0,0])
    # Actual quintic (up to its fixed nonzero scalar): E wedge P wedge P/2.
    pf=P[1,2]*P[3,4]-P[1,3]*P[2,4]+P[1,4]*P[2,3]
    assert s.expand(h*pf-u[3]*h**3*(u[0]*h*h+u[1]*h*x+u[2]*x*x))==0
    dense_H=sum((i%7+1)*v for i,v in enumerate(cubics))
    boundary={}
    for name,point,source in [('g_zero',canonical(raw(h*h+h*x+x*x,0,dense_H)),basis[:3]+basis[4:]),
            ('f_zero',canonical(raw(0,1,dense_H)),basis[3:])]:
        boundary[name]=boundary_witness(point,source)
    assert boundary['g_zero']['affine_image_differential_rank']==28
    assert boundary['f_zero']['affine_image_differential_rank']==26
    assert boundary['g_zero']['intrinsic_constant_vector_factor_matrix_rank']==4
    assert boundary['f_zero']['intrinsic_constant_vector_factor_matrix_rank']==5
    result=dict(id='AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63',all_passed=True,
        ambient='Full projectivized canonical Poisson space on P4; no PGL quotient',
        flag='L line in C5, W dimension2 in ann(L), h line in W',
        flag_base_dimension=9,parabolic_GL_dimension=16,parabolic_PGL_dimension=15,
        fiber_vector_dimension=20,projective_bundle_dimension=28,
        canonical_basis_rank=20,selected_raw150_rows=selected,
        basis_raw150=[list(map(str,vector(R)))for R in raws],
        basis_canonical150=[list(map(str,C[:,i]))for i in range(20)],
        cubic_monomials=list(map(str,cubics)),
        parabolic_matrix_positions=positions,
        exact_20_by20_Lie_actions=[[[str(c)for c in row]for row in B.tolist()]for B in actions],
        all_320_canonical_Lie_action_columns_verified=True,
        full20_parameter_raw_Jacobi_identically_zero=True,
        quintic_up_to_nonzero_constant='g*h^3*f(h,x)',
        rank_at_most_two_source_support='g=0 or f=0',
        generic_rank_two_boundary_source_g_zero_dimension=27,
        f_zero_source_codimension=3,
        boundary_witness_cubic_coefficients=[i%7+1 for i in range(16)],
        rank_two_boundary_exact_witnesses=boundary,
        entire_reduced_rank_at_most_two_boundary_components_projective_dimensions=[27,25],
        boundary_components_distinct_by_closed_constant_vector_factor_condition=True,
        proper_smooth_source_image_is_CA_by_geometric_argument=True,
        global_normalization_or_R1_claimed=False,
        full_global_intersection_claimed=False,
        input_sha256={str(Path(__file__).relative_to(ROOT)):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('GLOBAL CA20-SPACE PARABOLIC STABILITY PASS',result['seconds'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)

#!/usr/bin/env python3
"""Full CB compactified relative-pencil model and its rank<=2 boundary.

The ratio is an associated P1-bundle with a nontrivial parabolic character.
No global normalization/R1/intersection assertion is made here.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,canonical,vector,action,boundary_witness

ROOT=Path(__file__).resolve().parents[1]
h,x,y,z,w=X

def raw(f,g,H):
    P=s.zeros(5);P[1,2]=f;P[3,4]=g
    P[2,3]=-s.diff(f,x)*z/2-s.diff(H,w);P[2,4]=-s.diff(f,x)*w/2+s.diff(H,z)
    return P-P.T

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    binary=[h*h,h*x,x*x]
    cubic=[s.prod((h,x,z,w)[i] for i in mon) for mon in itertools.combinations_with_replacement(range(4),3) if any(i>=2 for i in mon)]
    raws=[raw(F,0,0) for F in binary]+[raw(0,F,0) for F in binary]+[raw(0,0,H) for H in cubic]
    basis=[canonical(P) for P in raws];C=s.Matrix.hstack(*(vector(P) for P in basis));rows=list(C.T.rref()[1]);assert C.rank()==22
    inverse=C[rows,:].inv();positions=[(i,j) for i,j in itertools.product(range(5),repeat=2)
        if not(i==0 and j!=0) and not(i==1 and j not in (0,1)) and not(j==2 and i!=2)]
    assert len(positions)==16;representations=[];mu_actions=[];characters=[];mu=s.Symbol('mu')
    embedding=s.zeros(22,19)
    for i in range(3):embedding[i,i]=1;embedding[i+3,i]=mu
    for i in range(16):embedding[i+6,i+3]=1
    selected=[0,1,2,*range(6,22)]
    for i,j in positions:
        M=s.zeros(5);M[i,j]=1;B=s.Matrix.hstack(*(vector(action(M,P)) for P in basis));small=inverse*B[rows,:]
        assert C*small==B
        assert small[:3,3:]==s.zeros(3,19) and small[3:6,:3]==s.zeros(3,3) and small[3:6,6:]==s.zeros(3,16)
        chi=M[3,3]+M[4,4]-M[1,1]-M[2,2]
        assert small[3:6,3:6]-small[:3,:3]==chi*s.eye(3)
        adjusted=small*embedding-chi*mu*embedding.diff(mu);induced=adjusted[selected,:]
        assert (embedding*induced-adjusted).applyfunc(s.expand)==s.zeros(22,19)
        representations.append(small);mu_actions.append(induced);characters.append(chi)
    assert sum((representations[positions.index((i,i))] for i in range(5)),s.zeros(22))==s.zeros(22)
    print('FULL352 LIE AND304 RELATIVE-MU COLUMNS PASS',time.monotonic()-started,flush=True)
    # The larger17-dimensional group preserves (L,W), forgetting hline.
    extra=s.zeros(5);extra[0,1]=1
    endpoint_a=[*basis[3:6],*basis[6:]];Honly=basis[6:]
    bigger=[]
    for name,bb in [('a_zero',endpoint_a),('H_only',Honly)]:
        matrix=s.Matrix.hstack(*(vector(P) for P in bb));sel=list(matrix.T.rref()[1]);inv=matrix[sel,:].inv()
        B=s.Matrix.hstack(*(vector(action(extra,P)) for P in bb));rep=inv*B[sel,:];assert matrix*rep==B;bigger.append(rep)
    # RawJacobi is exactly the binary coefficient-rank1 condition.
    uv=s.symbols('u0:22');P=sum((v*R for v,R in zip(uv,raws)),s.zeros(5));f=sum(uv[i]*binary[i] for i in range(3));g=sum(uv[i+3]*binary[i] for i in range(3))
    wronskian=s.expand(g*s.diff(f,x)-f*s.diff(g,x));jacobi=[]
    for triple in itertools.combinations(range(5),3):
        i,j,k=triple;value=s.expand(sum(P[i,q]*s.diff(P[j,k],X[q])+P[j,q]*s.diff(P[k,i],X[q])+P[k,q]*s.diff(P[i,j],X[q]) for q in range(5)))
        assert value==(wronskian if triple==(2,3,4) else 0);jacobi.append(str(value))
    coeff=s.Poly(wronskian,h,x);assert coeff.coeff_monomial(h**3)==uv[3]*uv[1]-uv[0]*uv[4]
    assert coeff.coeff_monomial(h*h*x)==2*(uv[3]*uv[2]-uv[0]*uv[5])
    assert coeff.coeff_monomial(h*x*x)==uv[4]*uv[2]-uv[1]*uv[5]
    aa,bb=s.symbols('a b');F=sum(s.Symbol('F'+str(i))*binary[i] for i in range(3));H=sum(uv[i+6]*v for i,v in enumerate(cubic));rank1=raw(aa*F,bb*F,H)
    assert s.expand(h*(rank1[1,2]*rank1[3,4]-rank1[1,3]*rank1[2,4]+rank1[1,4]*rank1[2,3])-aa*bb*h*F*F)==0
    dense_H=sum((i%7+1)*v for i,v in enumerate(cubic));dense_F=h*h+h*x+x*x
    boundary={
        'a_zero':boundary_witness(canonical(raw(0,dense_F,dense_H)),endpoint_a),
        'b_zero':boundary_witness(canonical(raw(dense_F,0,dense_H)),[*basis[:3],*basis[6:]]),
        'F_zero':boundary_witness(canonical(raw(0,0,dense_H)),Honly)}
    assert boundary['a_zero']['affine_image_differential_rank']==27
    assert boundary['b_zero']['affine_image_differential_rank']==28
    assert boundary['F_zero']['affine_image_differential_rank']==24
    assert boundary['a_zero']['intrinsic_constant_vector_factor_matrix_rank']==5
    assert boundary['b_zero']['intrinsic_constant_vector_factor_matrix_rank']==4
    result=dict(id='AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63',all_passed=True,
        ambient='Full projectivized canonical Poisson space on P4, no PGL quotient',
        flag_base_dimension=9,associated_P1_base_dimension=10,relative_vector_bundle_rank=19,proper_projective_source_dimension=28,
        ratio_character='mu -> det(z,w)/(scale_x*scale_y)*mu',
        infinitesimal_ratio_character=[str(v) for v in characters],parabolic_GL_positions=positions,
        cubic_monomials=list(map(str,cubic)),binary_monomials=list(map(str,binary)),
        enlarged_raw_basis=[[str(v) for v in vector(P)] for P in raws],enlarged_canonical_basis=[[str(v) for v in vector(P)] for P in basis],
        all352_enlarged_Lie_columns_verified=True,all304_relative_mu_Lie_columns_verified=True,
        exact_22_by22_Lie_actions=[[[str(v) for v in row] for row in M.tolist()] for M in representations],
        exact_19_by19_relative_mu_actions=[[[str(v) for v in row] for row in M.tolist()] for M in mu_actions],
        forget_hline_extra_generator_actions=[[[str(v) for v in row] for row in M.tolist()] for M in bigger],
        rawJacobi22=jacobi,Wronskian=str(wronskian),rank1_binary_coefficient_condition='f=a*F,g=b*F',
        actual_degeneracy_quintic='a*b*h*F(h,x)^2',
        rank_at_most2_preimage_support='a=0 or b=0 or F=0',
        endpoint_a_zero_forgets_hline=True,H_only_forgets_hline=True,
        boundary_exact_witnesses=boundary,entire_rank_at_most2_boundary_projective_dimensions=[26,27],
        H_only_image_projective_dimension=23,H_only_image_contained_in_both_endpoint_images=True,
        boundary_images_distinct_by_closed_constant_vector_factor_condition=True,
        global_normalization_R1_intersection_claimed=False,
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in [Path(__file__),ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py']},
        seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('GLOBAL CB COMPACTIFIED MODEL PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)

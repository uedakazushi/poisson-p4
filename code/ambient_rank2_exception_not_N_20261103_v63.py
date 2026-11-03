#!/usr/bin/env python3
"""Exact intrinsic obstruction for the exceptional rank2 intersection witness.

This does not itself prove its CB membership, family dimension, or global
intersection exhaustion; those are the separate Bruhat-table calculation.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,canonical,raw,vector

ROOT=Path(__file__).resolve().parents[1]

def witness(M):
    cols=list(M.rref()[1]);rows=list(M[:,cols].T.rref()[1]);return dict(rank=len(cols),rows=rows,columns=cols,determinant=str(M[rows,cols].det()))

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic();h,x,y,z,w=X;U=s.Matrix([h,x,z,w]);V=s.Matrix([0,4*x*x,0,-2*x*z,3*h*h+z*z+2*x*w])
    ey=s.Matrix([0,0,1,0,0]);P=ey*V.T-V*ey.T
    assert P==raw(-4*x*x,0,3*h*h*z+2*w*x*z+z**3/3)
    Q=canonical(P);C=s.Matrix([Q[2,i] for i in [0,1,3,4]])
    expected=s.Matrix([-8*h*x/5,12*x*x/5,-18*x*z/5,3*h*h+z*z+2*x*w/5]);assert C==expected
    assert sum(s.diff(C[i],U[i]) for i in range(4))==0
    # No constant covector annihilates the canonical bivector.
    ell5=s.Matrix(s.symbols('ell0:5'));mon2=[s.prod(X[i] for i in ids) for ids in itertools.combinations_with_replacement(range(5),2)]
    rows=[s.Poly(v,*X).coeff_monomial(m) for v in Q*ell5 for m in mon2];kernel_matrix,_=s.linear_eq_to_matrix(rows,list(ell5));assert kernel_matrix.rank()==5
    # Full intrinsic condition E wedge v wedge Q=0 has unique factorL=ey.
    E=s.Matrix(X);vv=s.Matrix(s.symbols('v0:5'));mon3=[s.prod(X[i] for i in ids) for ids in itertools.combinations_with_replacement(range(5),3)];eq=[]
    for i,j,k,l in itertools.combinations(range(5),4):
        value=s.expand(E[i]*(vv[j]*Q[k,l]-vv[k]*Q[j,l]+vv[l]*Q[j,k])-E[j]*(vv[i]*Q[k,l]-vv[k]*Q[i,l]+vv[l]*Q[i,k])
            +E[k]*(vv[i]*Q[j,l]-vv[j]*Q[i,l]+vv[l]*Q[i,j])-E[l]*(vv[i]*Q[j,k]-vv[j]*Q[i,k]+vv[k]*Q[i,j]))
        eq += [s.Poly(value,*X).coeff_monomial(m) for m in mon3]
    factor_matrix,_=s.linear_eq_to_matrix(eq,list(vv));assert factor_matrix.rank()==4 and factor_matrix.nullspace()==[ey]
    a,b,c,d,lam=s.symbols('a b c d lambda');ell=a*h+b*x+c*z+d*w
    Cell=sum(s.diff(ell,U[i])*C[i] for i in range(4));difference=s.Poly(Cell-lam*ell*ell,*U)
    assert difference.coeff_monomial(w*w)==-lam*d*d
    assert difference.coeff_monomial(h*h)==3*d-lam*a*a
    assert difference.coeff_monomial(z*z)==d-lam*c*c
    # lambda=0 permits no nonzero ell because the four C components are linearlyindependent.
    linear_coeff=[s.Poly(Cell,*U).coeff_monomial(m) for m in [h*h,h*x,h*z,h*w,x*x,x*z,x*w,z*z,z*w,w*w]]
    first_matrix,_=s.linear_eq_to_matrix(linear_coeff,[a,b,c,d]);assert first_matrix.rank()==4
    alpha=s.Symbol('alpha');k=a*h+b*x+c*z+d*w
    condition=s.expand(sum(s.diff(k,U[i])*C[i] for i in range(4))+s.Rational(48,5)*x*k+alpha*x*x)
    expected_condition=8*a*h*x+(12*b+alpha)*x*x+6*c*x*z+d*(3*h*h+z*z+10*x*w)
    assert condition==s.expand(expected_condition)
    equations=s.Poly(condition,*U).coeffs();second_matrix,_=s.linear_eq_to_matrix(equations,[a,b,c,d,alpha])
    assert second_matrix.nullspace()==[s.Matrix([0,-s.Rational(1,12),0,0,1])]
    paths=[Path(__file__),ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',
        ROOT/'research/notes/AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63.md']
    result=dict(id='AMBIENT_RANK2_EXCEPTION_NOT_N_20261103_v63',all_passed=True,field='QQ',raw_vector_field=list(map(str,V)),
        canonical_quotient_vector_field=list(map(str,C)),canonical_covector_kernel_matrix=kernel_matrix.tolist(),canonical_covector_kernel_witness=witness(kernel_matrix),
        intrinsic_constant_factor_matrix=factor_matrix.tolist(),intrinsic_factor_witness=witness(factor_matrix),intrinsic_factor_kernel=[0,0,1,0,0],
        self_cofactor_equation_key_coefficients=['-lambda*d^2','3*d-lambda*a^2','d-lambda*c^2'],
        lambda_zero_matrix=first_matrix.tolist(),lambda_zero_rank=4,
        necessary_N_first_flag_conclusion='ell proportionalx,lambda*ell=12x/5',
        necessary_N_second_flag_equation=str(condition),second_flag_matrix=second_matrix.tolist(),second_flag_kernel=list(map(str,second_matrix.nullspace()[0])),
        conclusion='Witness is NOT in global N: no canonical covectorCasimirs excludesA_f; uniqueintrinsicfactor+twoflagcoefficient identities excludeN_g0. This source does not establish CB membership or new-family dimension/exhaustion.',
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(result,indent=2,default=str)+'\n');print('EXCEPTIONAL RANK2 WITNESS NOT N PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)

#!/usr/bin/env python3
"""Independent finite matrix/Jacobi checks for the regular-pencil G proof.

The global Jordan/common-image/geometry argument is a separate hand proof.
These exact checks cover all three single-eigenvalue normal forms.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time

import classification_degenerate_quadric_q4_nil31_I_review as q

ROOT=Path(__file__).resolve().parents[1]

def zero():return [[F()for _ in range(5)]for _ in range(5)]
def eye():return [[F(i==j)for j in range(5)]for i in range(5)]
def scale(A,s):return [[s*x for x in r]for r in A]
def add(*As):return [[sum(A[i][j]for A in As)for j in range(5)]for i in range(5)]
def comm(A,B):return add(q.matmul(A,B),scale(q.matmul(B,A),-1))
def flat(A):return [x for r in A for x in r]
def powmat(A,k):
    out=eye()
    for _ in range(k):out=q.matmul(out,A)
    return out
def encode(A):return [[str(x)for x in r]for r in A]

def rank_certificate(A):
    rank,cols=q.rank(A)
    rows=q.rank(q.transpose([[r[i]for i in cols]for r in A]))[1]
    determinant=q.determinant([[A[i][j]for j in cols]for i in rows])
    assert determinant and len(rows)==len(cols)==rank
    return dict(rank=rank,rows=rows,columns=cols,determinant=str(determinant))

def quadric(Q):
    return q.add(*(q.atom(Q[i][j],i,j)for i,j in it.product(range(5),repeat=2)))

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args();started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(25,26));resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,)*2)
    shift=[[F(j==i+1)for j in range(5)]for i in range(5)]
    cases=[]
    for name in ['5','4+1','3+2']:
        if name=='5':
            Q=[[F(i+j==4)for j in range(5)]for i in range(5)]
            N=shift;S=[[F(2-i)if i==j else F()for j in range(5)]for i in range(5)]
            basis=[S]
        elif name=='4+1':
            Q=[[F(i+j==3 and i<4 and j<4 or i==j==4)for j in range(5)]for i in range(5)]
            N=[[F(j==i+1 and j<4)for j in range(5)]for i in range(5)]
            diag=[F(3,2),F(1,2),F(-1,2),F(-3,2),F()]
            S=[[diag[i]if i==j else F()for j in range(5)]for i in range(5)]
            R=zero();R[0][4]=F(1);R[4][3]=F(-1);basis=[S,R]
            assert q.rank(R)[0]==2
            C=q.matmul(R,q.inverse(Q))
            assert all(C[i][j]*C[k][l]-C[i][k]*C[j][l]+C[i][l]*C[j][k]==0 for i,j,k,l in it.combinations(range(5),4))
            assert powmat(R,2)!=zero() # Matrix square is not the exterior square.
        else:
            Q=[[F((-1)**i)if i+j==4 else F()for j in range(5)]for i in range(5)]
            N=powmat(shift,2);diag=[F(1),F(1,2),F(),F(-1,2),F(-1)]
            S=[[diag[i]if i==j else F()for j in range(5)]for i in range(5)]
            basis=[S,shift,powmat(shift,3)]
        assert q.determinant(Q) and q.matmul(q.transpose(N),Q)==q.matmul(Q,N)
        assert comm(S,N)==N
        for T in basis:assert add(q.matmul(q.transpose(T),Q),q.matmul(Q,T))==zero()
        for T in basis[1:]:assert comm(T,N)==zero()
        assert q.rank(q.transpose(list(map(flat,basis))))[0]==len(basis)
        # All10 skew A unknowns plus alpha,beta in [AQ,N]=alpha I+beta N.
        columns=[]
        for i,j in it.combinations(range(5),2):
            A=zero();A[i][j]=F(1);A[j][i]=F(-1)
            columns.append(flat(comm(q.matmul(A,Q),N)))
        columns.extend([flat(scale(eye(),-1)),flat(scale(N,-1))])
        matrix=q.transpose(columns);certificate=rank_certificate(matrix)
        assert 12-certificate['rank']==len(basis)
        record=dict(Jordan_partition=name,Q=encode(Q),N=encode(N),grading=encode(S),
            full_skew_normalizer_matrix=encode(matrix),normalizer_rank_certificate=certificate,
            full_normalizer_dimension=len(basis),normalizer_basis=[encode(T)for T in basis])
        if name=='4+1':record.update(centralizer_matrix_rank=2,associated_bivector_exterior_square_zero=True,
            centralizer_matrix_square_is_not_zero=True)
        cases.append(record)
    # Explicit preserving automorphisms in the only residual 3+2 case.
    R=shift;R3=powmat(R,3)
    K=scale(add(q.matmul(N,S),q.matmul(S,N)),F(1,2))
    L=add(scale(N,F(1,2)),K)
    assert add(q.matmul(q.transpose(K),Q),q.matmul(Q,K))==zero()
    assert comm(K,N)==powmat(N,2)
    QN=q.matmul(Q,N)
    delta=lambda T,B:add(q.matmul(q.transpose(T),B),q.matmul(B,T))
    assert delta(L,Q)==QN and delta(L,QN)==zero()
    assert delta(S,Q)==zero() and delta(S,QN)==scale(QN,-1)
    for T in [R,R3]:assert delta(T,Q)==delta(T,QN)==zero()
    preserving_matrix=q.transpose([flat(T)for T in [eye(),R,R3,S,L]])
    preserving_certificate=rank_certificate(preserving_matrix);assert preserving_certificate['rank']==5
    # All six actual source basis directions, and a literal nonPoisson member.
    q1,q2=quadric(Q),quadric(QN);Qinv=q.inverse(Q)
    As=[q.matmul(T,Qinv)for T in [S,R,R3]]
    assert q.rank(As[0])[0]==4
    source=[]
    for polynomial in [q1,q2]:
        for A in As:source.append(q.antisymmetric({ij:q.scale(polynomial,A[ij[0]][ij[1]])for ij in q.PAIRS}))
    source_matrix=q.transpose([q.vec(P)for P in source]);source_certificate=rank_certificate(source_matrix)
    assert source_certificate['rank']==6
    actualJ=q.jacobi(source[0]);vector=q.jvec(actualJ)
    assert len(vector)==350 and any(vector)
    witness_index=next(i for i,v in enumerate(vector)if v)
    paths=[ROOT/'research/classification/REGULAR_PENCIL_CI22_ABELIAN_CONORMAL_NONDOMINANCE.md',
           ROOT/'research/classification/RANK_ZERO_SURFACE_GENERIC_FIRST_JET.md',
           ROOT/'code/classification_degenerate_quadric_q4_nil31_I_review.py']
    result=dict(task_id='G',all_passed=True,no_producer_import=True,no_CAS_import=True,
        all_three_single_eigenvalue_normalizer_cases=cases,
        explicit_Jordan32_preserving_fields={name:encode(T)for name,T in zip(['R','R3','S','L'],[R,R3,S,L])},
        preserving_fields_mod_scalar_independence_certificate=preserving_certificate,
        all_four_preserving_pencil_derivatives_exact=True,whole_six_source_rank_certificate=source_certificate,
        nonPoisson_rank4_constant_A_times_nondegenerate_q1_full350=[str(v)for v in vector],
        explicit_nonzero_projective_Jacobi_witness=dict(index=witness_index,value=str(vector[witness_index]),row_label=str(q.ROWS[witness_index])),
        mathematical_scope='Finite arithmetic verifies the canonical matrices/all6 source and a proper Jacobi equation. The all-eigenvalue common-image factor, complete normal-form coverage, conormal necessity and dimension24 image argument are separately reviewed mathematical proofs, not inferred from these canonical computations.',
        fixed_residual_Poisson_source_dimension_upper=5,actual_preserving_group_dimension_lower=4,
        saturated_projective_dimension_upper=24,all_member_degenerate_pencils_retained=True,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
        input_note_snapshots={str(p.relative_to(ROOT)):p.read_text()for p in paths[:2]},
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_bytes=Path(__file__).read_text(),
        seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['all_passed','saturated_projective_dimension_upper','seconds','maximum_RSS_KiB']}))

if __name__=='__main__':main()

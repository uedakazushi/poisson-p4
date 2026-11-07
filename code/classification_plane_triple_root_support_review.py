#!/usr/bin/env python3
"""Independent full four-dimensional F1 canonical source support.

Fraction only; reuses previously independent polynomial arithmetic. All ten
quadratic polarizations of all350 Jacobi and all5 raw Pfaffians are retained.
"""
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement
from pathlib import Path
import hashlib,json,resource,signal,time
import classification_degenerate_quadric_q4_nil31_I_review as q

ROOT=Path(__file__).resolve().parents[1]

def tensor(v):
    return q.antisymmetric({ij:{e+(0,0,0):c for (pair,e),c in zip(q.COLS,v)if pair==ij and c}for ij in q.PAIRS})

def sumtensor(A,B):
    return q.antisymmetric({ij:q.add(A.get(ij,{}),B.get(ij,{}))for ij in q.PAIRS})

def pfaffians(P):
    out=[]
    for missing in range(5):
        i,j,k,l=[a for a in range(5)if a!=missing]
        out.append(q.add(q.mul(P.get((i,j),{}),P.get((k,l),{})),q.scale(q.mul(P.get((i,k),{}),P.get((j,l),{})),-1),q.mul(P.get((i,l),{}),P.get((j,k),{}))))
    return out

def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,)*2)
    job=json.loads((ROOT/'research/classification/jobs/plane_triple_root_F1_support_review.json').read_text());out=ROOT/'reports/CLASSIFICATION_PLANE_TRIPLE_ROOT_F1_SUPPORT_REVIEW_20261107.json'
    if out.exists():raise FileExistsError(out)
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    raw=json.loads((ROOT/'reports/CLASSIFICATION_PLANE_M4_TRIPLE_ROOT_QUINTICS_20261107.json').read_text());row=raw['cases'][1];assert row['epsilon']==1
    assert raw['full150_column_labels']==[dict(pair=list(ij),monomial=list(e))for ij,e in q.COLS]
    K=[[F(c)for c in v]for v in row['full150_kernel']];assert len(K)==4 and q.rank(K)[0]==4
    M=[[F(c)for c in v]for v in row['full1075x150_strict_Casimir_matrix']]
    assert all(not sum(a*b for a,b in zip(r,v))for r in M for v in K)
    P=list(map(tensor,K));pairs=list(combinations_with_replacement(range(4),2));Jac=[];Pf=[]
    for i,j in pairs:
        J=q.jvec(q.jacobi(P[i]if i==j else sumtensor(P[i],P[j])));pf=pfaffians(P[i]if i==j else sumtensor(P[i],P[j]))
        if i!=j:
            J=[a-b-c for a,b,c in zip(J,q.jvec(q.jacobi(P[i])),q.jvec(q.jacobi(P[j])))];pi,pj=pfaffians(P[i]),pfaffians(P[j]);pf=[q.add(a,q.scale(b,-1),q.scale(c,-1))for a,b,c in zip(pf,pi,pj)]
        Jac.append(J);Pf.append(pf)
    Jrows=q.transpose(Jac);pfrows=[[poly.get(e+(0,0,0),F(0))for poly in blocks]for blocks in zip(*Pf)for e in q.monomials(4)]
    # The transpose above has five blocks, each consisting of ten polynomials.
    D=[]
    for pf in Pf:D.append(q.add(*(q.scale(q.mul(q.atom(1,k),v),(-1)**k)for k,v in enumerate(pf))))
    Drows=[[p.get(e+(0,0,0),F(0))for p in D]for e in q.monomials(5)]
    basis=[];inds=[]
    for i,r in enumerate(Jrows):
        if q.rank(basis+[r])[0]>len(basis):basis.append(r);inds.append(i)
    # Complete rowspace rank is checked, no leading-subideal conclusion.
    assert q.rank(Jrows)[0]==len(basis)
    rr=len(basis);D_in_rowspace=all(not any(r) for r in Drows) if not basis else all(q.rank(basis+[r])[0]==rr for r in Drows)
    result=dict(task_id='E_reduced_nonnormal',all_arithmetic_passed=True,all_four_source_parameters_retained=True,all1075_saved_rows_times_all4_independent_kernel_zero=True,independent_source_completeness_reused='PLANE_M4_REPEATED_SEPARATE_REVIEW.md full1075matrix independently reconstructed there.',parameter_quadratic_pairs=[list(p)for p in pairs],full350_Jacobi_quadratic_matrix=[[str(c)for c in r]for r in Jrows],independent_original350_rows=inds,Jacobi_rowspace_basis=[[str(c)for c in r]for r in basis],full5_by70_Pfaffian_quadratic_matrix=[[str(c)for c in r]for r in pfrows],full126_actualD_quadratic_matrix=[[str(c)for c in r]for r in Drows],all_raw_Pfaffians_identically_zero=not any(any(r)for r in pfrows),all_actualD_identically_zero=not any(any(r)for r in Drows),all_actualD_rows_in_full_Jacobi_rowspace=D_in_rowspace,scope='Complete F1 canonical necessarysource, all350 equations and actualD. Empty only if fullDzero on every Poissonpoint is certified, not from tangent dimension or an arbitraryslice.',inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in ['all_arithmetic_passed','Jacobi_rowspace_basis','all_raw_Pfaffians_identically_zero','all_actualD_identically_zero','all_actualD_rows_in_full_Jacobi_rowspace','seconds']}))

if __name__=='__main__':main()

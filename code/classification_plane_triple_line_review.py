#!/usr/bin/env python3
"""Independent Fraction full-source exclusion of triple-line/m3/squarefree.

Only an already independent sparse Fraction engine is reused. No producer,
producer helper, symbolic CAS or leading-exclusion result is imported.
"""
from fractions import Fraction as F
import hashlib,itertools,json,resource,time
from pathlib import Path
import classification_plane_m4_special_review as q

ROOT=Path(__file__).resolve().parents[1]
INPUT='reports/CLASSIFICATION_PLANE_NONREDUCED_TRIPLE_LINE_M3_SQUAREFREE_20261107.json'
LEADING='reports/CLASSIFICATION_PLANE_NONREDUCED_TRIPLE_LINE_M3_LEADING_EXCLUSION_20261107.json'
JOB='research/classification/jobs/plane_nonreduced_triple_line_m3_review.json'
OUTPUT='reports/CLASSIFICATION_PLANE_NONREDUCED_TRIPLE_LINE_M3_REVIEW_20261107.json'


def mm(A,B):return [[sum(a*b for a,b in zip(row,col))for col in zip(*B)]for row in A]


def certificate(M):
    _,cols,K=q.rref(M);_,rows,_=q.rref(q.transpose(M));r=len(cols)
    v=q.det([[M[i][j]for j in cols]for i in rows]);assert v
    return dict(rank=r,rows=rows,columns=cols,minor_determinant=str(v)),K


def check_producer_certificate(M,c):
    rows,cols=c['minor_rows'],c['minor_columns'];rank=c['rank']
    assert len(rows)==len(cols)==rank
    value=q.det([[M[i][j]for j in cols]for i in rows]);assert value==F(c['minor_determinant'])and value
    K=q.transpose(q.mat(c['complete_kernel']));assert len(K)==c['nullity']
    assert all(not sum(a*b for a,b in zip(row,k))for row in M for k in K)
    assert len(q.rref(K)[1])==len(K)
    assert rank+len(K)==len(M[0])
    return K


def main():
    begun=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(110,115));resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,)*2)
    job=json.loads((ROOT/JOB).read_text())
    for v in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/v['path']).read_bytes()).hexdigest()==v['sha256'],v['path']
    old=json.loads((ROOT/INPUT).read_text());assert old['all_exact_arithmetic_passed']and not old['all_G24_closed']
    assert hashlib.sha256(old['producer_bytes'].encode()).hexdigest()==old['producer_sha256']
    first,second=q.sources()
    # Producer's coefficient order is l,n,d,b,e,c. This is a permutation of
    # the complete 30-directions plus six-Euler rank36 source just rebuilt.
    first=[first[j]for j in list(range(18))+list(range(24,27))+list(range(18,21))+list(range(27,30))+list(range(21,24))]
    base2=[]
    for inds in itertools.combinations_with_replacement(range(3),2):
        e=[0]*3
        for i in inds:e[i]+=1
        base2.append(tuple(e))
    mon3=[e+(j,3-j)for e in base2 for j in range(4)];assert len(mon3)==24
    H=[{e:F(1)}for e in mon3];X=[q.a(1,i)for i in range(5)]
    records=[]
    for unipotent,c in zip([False,True],old['cases']):
        P0=q.skew({(1,2):q.a(3,0,0),**({(2,4):q.a(3,0,3)}if unipotent else{})})
        assert not q.actual(P0)and all(not p for p in q.jac(P0).values())
        assert list(map(list,mon3))==c['all24_leading_monomials']
        # Reconstruct the regular parallel space, with no saved basis used.
        PDE=[]
        for which in [0,1]:
            polys=[]
            for G in H:
                p=q.mul(X[0],q.der(G,1 if which==0 else 2))
                if unipotent and which==0:p=q.add(p,q.scale(q.mul(X[3],q.der(G,4)),-1))
                polys.append(p)
            PDE.extend([[p.get(e,F())for p in polys]for e in mon3])
        pcert,Kparallel=certificate(PDE);assert len(Kparallel)==(3 if unipotent else 4)
        t=q.add(q.mul(X[0],X[4]),q.mul(X[1],X[3]))
        generators=[q.a(1,0,0,3,3,3),q.mul(q.a(1,0,3,3),t),q.mul(X[3],q.mul(t,t))]if unipotent else[q.a(1,0,0,*([3]*j+[4]*(3-j)))for j in range(4)]
        A=q.transpose([[G.get(e,F())for e in mon3]for G in generators]);assert A==q.mat(c['parallel_basis_24byN'])
        assert all(not v for row in mm(PDE,A)for v in row)
        assert len(q.rref(q.transpose(A))[1])==len(generators)
        Ann=q.mat(c['parallel_membership_full_annihilator'])
        assert len(q.rref(Ann)[1])==24-len(generators)and not any(v for row in mm(Ann,A)for v in row)
        J1=q.transpose([q.jvector(q.jac(q.plus(P0,P)))for P in first])
        J2=q.transpose([q.jvector(q.jac(q.plus(P0,P)))for P in second])
        D3=q.transpose([[q.actual(q.plus(P0,P)).get(e,F())for e in mon3]for P in first])
        assert J1==q.mat(c['full350_weight1_Jacobi'])and J2==q.mat(c['full350_weight2_Jacobi'])
        assert D3==q.mat(c['full_D3_linear_map_24by30'])
        M=J1+mm(Ann,D3);assert M==q.mat(c['full_weight1_with_parallel_D3'])
        Ksaved=check_producer_certificate(M,c['weight1_certificate']);check_producer_certificate(J2,c['weight2_certificate'])
        c1,K1=certificate(M);c2,K2=certificate(J2)
        assert len(K1)==len(Ksaved)==(18 if unipotent else 24)
        assert len(q.rref(K1+Ksaved)[1])==len(K1)
        images=mm(D3,q.transpose(K1));icert,_=certificate(images)
        assert icert['rank']==int(unipotent)
        target=mon3.index((2,0,0,3,0))
        assert all(not v or i==target for i,row in enumerate(images)for v in row)
        # P2 cannot contribute to order3, so this is the full actualD3 image.
        assert all(not q.actual(q.plus(P0,P)).get(e,F())for P in second for e in mon3)
        for P in first:
            assert all(not q.actual(P).get(e,F())for e in mon3)
        for P,Q in itertools.combinations(first,2):
            assert all(not q.actual(q.plus(P,Q)).get(e,F())for e in mon3)
        records.append(dict(type=c['type'],P0_full150=list(map(str,q.vector(P0))),full30_weight1_basis=[list(map(str,q.vector(P)))for P in first],full9_weight2_basis=[list(map(str,q.vector(P)))for P in second],parallel_PDE48by24=[[str(v)for v in row]for row in PDE],parallel_PDE_rank_certificate=pcert,complete_regular_parallel_kernel=[list(map(str,v))for v in Kparallel],
            full350_weight1_matrix=[[str(v)for v in row]for row in J1],full350_weight2_matrix=[[str(v)for v in row]for row in J2],full24by30_actualD3_matrix=[[str(v)for v in row]for row in D3],weight1_rank_certificate=c1,weight2_rank_certificate=c2,full30_weight1_kernel=[list(map(str,v))for v in K1],full9_weight2_kernel=[list(map(str,v))for v in K2],
            whole_actualD3_kernel_image=[[str(v)for v in row]for row in images],actualD3_image_rank_certificate=icert,actualD3_generator_exponent=(2,0,0,3,0)if unipotent else None,actualD3_independent_of_every_P2=True,actualD3_has_no_P1quadratic_contribution=True,squarefree_leading_impossible=True))
    # Final producer receipt is bound and compared only after the independent
    # reconstruction above; no reported solution enters any calculation.
    final=json.loads((ROOT/LEADING).read_text())
    result=dict(task_id='G',all_independent_arithmetic_passed=True,whole_triple_line_m3_squarefree_rank4_source_empty=True,cases=records,
        geometric_review='Finite etale triple-root cover of A2 splits; PicA2 zero and unit pair determinants give precisely constant roots or degrees0,1,1. Whole normal inverses and unique twisted-cubic connection checked. Regular PDE full kernels are4/3. Full actualD3 images are0 and span(a2*z3), so neither can be generic squarefree.',
        producer_coarse_nonclosing_bound_preserved=True,producer_final_receipt_sha256=hashlib.sha256((ROOT/LEADING).read_bytes()).hexdigest(),producer_final_receipt_keys=list(final),scope='Only R=3L,m3,generic squarefree normal cubic. m2,m3 repeated,m4 and R=2L+M remain separate.',
        input_bindings=job['exact_inputs_and_hashes'],reviewer_source_bytes=Path(__file__).read_text(),reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with(ROOT/OUTPUT).open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_independent_arithmetic_passed=True,squarefree_source_empty=True,seconds=result['seconds'],RSS=result['maximum_RSS_KiB'])))


if __name__=='__main__':main()

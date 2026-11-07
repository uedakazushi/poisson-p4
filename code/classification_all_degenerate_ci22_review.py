#!/usr/bin/env python3
"""Independent all-singular-pencil CI22 model review using Fraction only.

Completeness of symmetric-pencil canonical forms is a stated primary input.
This checker verifies all20 model directions, necessary conormal equations,
the whole four-parameter actual divisor, and the single-Jordan4 normalizer.
"""
import argparse,hashlib,itertools as it,json,resource,time
from pathlib import Path
from classification_cubic_scroll_abelian_review import Q,PAIRS,add,scale,mul,diff,parse,rref,monomials,jacobi_and_divisor

ROOT=Path(__file__).resolve().parents[1]


def matrix_kernel(rows,n):
    R=rref(rows);free=[i for i in range(n)if i not in R];basis=[]
    for j in free:
        v=[Q(i==j)for i in range(n)]
        for p in R:v[p]=-R[p].get(j,Q(0))
        basis.append(v)
    return R,free,basis


def mm(A,B):return [[sum(a*b for a,b in zip(row,col))for col in zip(*B)]for row in A]


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(55,60));resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,512*1024**2))
    names=[f'x{i}'for i in range(5)];q=[parse(t,names)for t in ['x0*x1+x3**2','x0*x2+x4**2']];mons2=monomials(5,2);mons5=monomials(5,5)
    bases=[]
    for g in range(2):
        for i,j in PAIRS:
            P=[[{}for _ in range(5)]for _ in range(5)];P[i][j]=q[g];P[j][i]=scale(q[g],-1);bases.append(P)
    def vector(P):return [P[i][j].get(e,Q(0))for i,j in PAIRS for e in mons2]
    source=[vector(P)for P in bases];assert len(rref(source))==20
    vertical=[]
    for a,b in it.product(range(5),repeat=2):
        P=[[{}for _ in range(5)]for _ in range(5)]
        for i,j in PAIRS:
            expression={}
            if j==a:
                e=[0]*5;e[i]+=1;e[b]+=1;expression[tuple(e)]=Q(1)
            if i==a:
                e=[0]*5;e[j]+=1;e[b]+=1;expression=add(expression,{tuple(e):Q(-1)})
            P[i][j]=expression
        vertical.append(vector(P))
    assert len(rref(vertical))==24 and len(rref(vertical+source))==44
    hA=[]
    for i,j in PAIRS:hA.append(add(mul(diff(q[0],i),diff(q[1],j)),scale(mul(diff(q[0],j),diff(q[1],i)),-1)))
    linear=[[p.get(e,Q(0))for p in hA]+[-q[0].get(e,Q(0)),-q[1].get(e,Q(0))]for e in mons2]
    R,free,K=matrix_kernel(linear,12);assert len(R)==10 and free==[PAIRS.index((1,3)),PAIRS.index((2,4))]
    assert all(all(sum(c*v[i]for i,c in enumerate(row))==0 for row in linear)for v in K)
    doubled=[row+[Q(0)]*12 for row in linear]+[[Q(0)]*12+row for row in linear];assert len(rref(doubled))==20
    selected=[bases[PAIRS.index((1,3))],bases[PAIRS.index((2,4))],bases[10+PAIRS.index((1,3))],bases[10+PAIRS.index((2,4))]]
    _,D=jacobi_and_divisor(selected);actual={};columns=list(it.combinations_with_replacement(range(4),2))
    for exponent,row in zip(mons5,D):
        for index,c in row.items():
            i,j=columns[index];parameter=tuple((i,j).count(k)for k in range(4));actual[parameter+exponent]=c
    expected=parse('-x0*(a*(x0*x1+x3**2)+c*(x0*x2+x4**2))*(b*(x0*x1+x3**2)+d*(x0*x2+x4**2))',['a','b','c','d',*names]);assert actual==expected
    # All possible singular blocks in size5; regular subpart is arbitrary.
    partitions=[]
    for count in range(1,6):
        for indices in it.combinations_with_replacement(range(3),count):
            total=sum(2*d+1 for d in indices)
            if total>5:continue
            if indices.count(0)>=2:case='common_radical_at_least2'
            elif indices==(2,):case='M2_contains_plane'
            elif indices==(1,):case='M1_regular2'
            elif indices==(0,):case='M0_regular4'
            elif indices==(0,1):case='M0_M1_regular1_rank2_member'
            else:raise AssertionError(indices)
            partitions.append(dict(minimal_indices=list(indices),regular_part_size=5-total,case=case))
    assert len(partitions)==9
    for texts,zero_coordinates in [(['x0*x3+x1*x4','x1*x3+x2*x4'],[3,4]),(['x0*x1+2*x3*x4','x0*x2+x4**2'],[0,4])]:
        assert all(not{e:c for e,c in parse(t,names).items()if all(e[i]==0 for i in zero_coordinates)}for t in texts)
    # M0+regular4, one eigenvalue: the complete skew normalizer is one-dimensional.
    n=4;pairs4=list(it.combinations(range(n),2));Qmat=[[Q(i+j==3)for j in range(n)]for i in range(n)];N=[[Q(j==i+1)for j in range(n)]for i in range(n)];blocks=[]
    for i,j in pairs4:
        A=[[Q(0)for _ in range(n)]for _ in range(n)];A[i][j]=1;A[j][i]=-1;T=mm(A,Qmat);TN=mm(T,N);NT=mm(N,T);blocks.append([TN[i][j]-NT[i][j]for i in range(n)for j in range(n)])
    blocks +=[[-Q(i==j)for i in range(n)for j in range(n)],[-N[i][j]for i in range(n)for j in range(n)]]
    normalizer=list(map(list,zip(*blocks)));RR,ff,kk=matrix_kernel(normalizer,8);assert len(RR)==7 and len(kk)==1 and kk[0][-2]==0
    A=[[Q(0)for _ in range(n)]for _ in range(n)]
    for value,(i,j)in zip(kk[0],pairs4):A[i][j]=value;A[j][i]=-value
    assert len(rref(A))==4
    proof=ROOT/'research/classification/ALL_DEGENERATE_CI22_ABELIAN_CONORMAL_EXCLUSION.md';regular=ROOT/'reports/CLASSIFICATION_REGULAR_PENCIL_CI22_REVIEW_V2_20261107.json'
    inputs=[Path(__file__),ROOT/'code/classification_cubic_scroll_abelian_review.py',proof,ROOT/'research/classification/REGULAR_PENCIL_CI22_ABELIAN_CONORMAL_REVIEW.md',regular]
    result=dict(task_id='G',all_passed=True,primary_canonical_form=dict(url='https://gauss.uc3m.es/fdopico/papers/jst2018.pdf',location='Section2 Theorem2.1 printedp5, blockdefinitionsp4',scope='Every complex symmetric pencil under congruence, including singular Md and finite/infinite regular blocks; not just generic bounded-rank pencils.'),all_size5_singular_block_partitions=partitions,
        complete_unconstrained_model_source_dimension=20,vertical_rank=24,source_plus_vertical_rank=44,all_fifteen_quadratic_conormal_rows=[[str(c)for c in row]for row in linear],conormal_rank=10,complete_A_alpha_beta_kernel=[[str(c)for c in v]for v in K],surviving_skew_entries=[[1,3],[2,4]],doubled_full_source_constraint_shape=[30,24],doubled_constraint_rank=20,whole_necessary_source_dimension=4,actual_full_four_parameter_divisor='-x0*(a*q1+c*q2)*(b*q1+d*q2)',all_actual_divisor_coefficients_verified=True,
        M0_Jordan4_full_normalizer_matrix=[[str(c)for c in row]for row in normalizer],M0_Jordan4_normalizer_rank=7,M0_Jordan4_normalizer_kernel=[[str(c)for c in v]for v in kk],M0_Jordan4_skew_generator_rank=4,
        reviewed_root_proof_snapshot=proof.read_text(),source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},reviewer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        conclusion='Every all-degenerate defining pencil of an integral nondegenerate CI22 has no rankzero source with irreducible reduced quintic D. Together with the adopted regular-pencil review, all integral nondegenerate CI22 rankzero G images have dimension<=24.',not_claimed=['All Poisson sources on CI22 zero','E coverage of reducible divisor sources','GlobalGcomplete'])
    with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,whole_necessary_source_dimension=4,seconds=result['seconds'])))


if __name__=='__main__':main()

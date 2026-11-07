#!/usr/bin/env python3
"""Independent whole four concurrent special-j sources, Fraction only.
Reuses only the earlier independent polynomial implementation, no producer/CAS.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,resource,signal,time
import classification_plane_m4_special_review as S
ROOT=Path(__file__).resolve().parents[1]
SOURCE='reports/CLASSIFICATION_PLANE_CONCURRENT_M4_SPECIAL_20261107.json'
JOB='research/classification/jobs/plane_concurrent_special_review.json'
OUT='reports/CLASSIFICATION_PLANE_CONCURRENT_SPECIAL_REVIEW_20261107.json'
MANIFEST='research/classification/inputs/plane_concurrent_m4_special_manifest_20261107.json'
def mm(A,B):return [[sum((a*b for a,b in zip(row,col)),F(0))for col in zip(*B)]for row in A]
def mat(data):return [[F(v)for v in row]for row in data]
def decode(poly):return {tuple(t['exponent']):F(t['coefficient'])for t in poly}
def kernels(M,K):
    assert all(not any(r)for r in mm(M,S.transpose(K)))
    assert len(S.rref(K)[1])==len(K)
def model(name):
    a=S.a;add=S.add;mul=S.mul
    if name=='J0_secant':v=[F(-1,8),F(5,8),0,0,0,0,F(5,8),F(-1,8)];G=[a(1,0,3,3,3,3),a(1,1,4,4,4,4)]
    elif name=='J0_compression':v=[F(-1,2),1,0,0,0,0,1,F(-1,2)];G=[a(1,0,3,3,3,4),a(1,1,3,4,4,4)]
    elif name=='I0_mixed':v=[0,F(1,2),0,0,0,0,1,F(-1,2)];G=[a(1,0,3,3,3,3),a(1,1,3,4,4,4)]
    else:
        assert name=='I0_three_distinct';v=[F(-1,2),F(1,4),0,F(1,2),0,F(3,8),F(1,2),F(-1,4)]
        G=[add(mul(add(a(2,0),a(-1,1)),a(4,3,3,3,4)),mul(a(1,1),add(a(F(3,4),3,3,3,3),a(6,3,3,4,4),a(-4,4,4,4,4))))]
    M=[[add(a(v[4*i+2*j],0),a(v[4*i+2*j+1],1))for j in range(2)]for i in range(2)]
    return list(map(F,v)),M,G
def main():
    started=time.monotonic();signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    assert not(ROOT/OUT).exists();job=json.loads((ROOT/JOB).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    man=json.loads((ROOT/MANIFEST).read_text())
    for b in man['files']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    packet=json.loads((ROOT/SOURCE).read_text());assert packet['all_arithmetic_passed']and packet['all_concurrent_special_rank2_G24_closed']is False
    first,second=S.sources();assert [S.vector(P)for P in first]==mat(packet['full30_weight1_basis'])and[S.vector(P)for P in second]==mat(packet['full9_weight2_basis'])
    f=S.add(S.a(1,0,0,1),S.a(-1,0,1,1));Wbase=[S.der(f,1),S.scale(S.der(f,0),-1),{},{},{}]
    six=S.monomials(6);D3mons=[e for e in S.monomials(5)if e[3]+e[4]==3];assert len(D3mons)==24
    output=[]
    for c in packet['cases']:
        values,M,G=model(c['type']);assert values==list(map(F,c['unique_raw_normal_matrix_linear_coefficients']))and G==[decode(P)for P in c['whole_leading_generators']]
        labels=[(i,j,k)for i in (3,4)for j in (3,4)for k in (0,1)];assert [list(t)for t in labels]==c['full_connection_basis_labels']
        con=[];rhs=[]
        for g in G:
            cols=[S.add(S.mul(S.a(1,k,j),S.der(g,i)),S.scale(S.mul(S.a(1,k),g),-int(i==j)))for i,j,k in labels]
            base=S.add(*(S.mul(Wbase[i],S.der(g,i))for i in range(2)))
            con.extend([[P.get(e,F(0))for P in cols]for e in six]);rhs.extend([-base.get(e,F(0))for e in six])
        assert con==mat(c['full_connection_equations'])and rhs==list(map(F,c['full_connection_rhs']))
        assert all(sum(a*b for a,b in zip(row,values))==v for row,v in zip(con,rhs));cm=S.minor(con,c['connection_rank_certificate'])
        normal=[S.add(*(S.mul(M[i][j],S.a(1,3+j))for j in range(2)))for i in range(2)]
        W=Wbase[:3]+normal;P0=S.wedge([{}, {},S.a(1),{},{}],W)
        assert S.vector(P0)==list(map(F,c['full150_P0']))and not any(S.jac(P0).values())and not S.actual(P0)
        for mon in ((3,3),(3,4),(4,4)):
            assert S.actual(S.plus(P0,S.skew({(3,4):S.a(1,*mon)})))==S.scale(S.mul(f,S.a(1,*mon)),3)
        lead=[]
        for i in range(3):
            for k in range(4,-1,-1):lead.append(S.a(1,i,*([3]*k+[4]*(4-k))))
        assert [list(next(iter(P)))for P in lead]==c['all15_leading_monomials']
        tr=S.add(M[0][0],M[1][1]);parallel=[S.add(*(S.mul(W[i],S.der(g,i))for i in (0,1,3,4)),S.scale(S.mul(tr,g),-1))for g in lead]
        L=S.transpose([[P.get(e,F(0))for e in six]for P in parallel]);assert L==mat(c['full_parallel_leading_matrix'])
        K=mat(c['full_parallel_leading_kernel']);kernels(L,K);pm=S.minor(L,c['parallel_leading_rank_certificate']);assert len(K)==(1 if c['type']=='I0_three_distinct'else 2)
        assert len(K)==15-c['parallel_leading_rank_certificate']['rank']
        J1=S.transpose([S.jvector(S.jac(S.plus(P0,P)))for P in first]);D3=[[S.actual(S.plus(P0,P)).get(e,F(0))for P in first]for e in D3mons]
        assert J1==mat(c['complete350_weight1_matrix'])and D3==mat(c['all24_D3_rows'])
        A=J1+D3;K1=mat(c['full30_weight1_kernel']);kernels(A,K1);m1=S.minor(A,c['weight1_rank_certificate']);assert len(K1)==12 and c['weight1_rank_certificate']['rank']==18
        fields=[S.skew({ij:S.add(*(S.scale(P.get(ij,{}),v)for P,v in zip(first,row)))for ij in S.PAIRS})for row in K1]
        J2=S.transpose([S.jvector(S.jac(S.plus(P0,P)))for P in second]);assert J2==mat(c['full350_weight2_matrix'])
        assert all(not v for row in J2 for j,v in enumerate(row)if j>=3)
        m2=S.minor(J2,c['full_weight2_rank_certificate']);assert c['full_weight2_rank_certificate']['rank']==3
        fields+=second[3:];assert len(fields)==18
        assert all(not poly or 2 in ij for P in [P0]+fields for ij,poly in P.items())
        assert all(not e[2]for P in [P0]+fields for poly in P.values()for e in poly)
        # Any linear combination is dz wedge z-independent W. This proves all
        # original quadratic compatibility identicallyzero and rank<=2 at once.
        assert all(not any(row)for row in mat(c['full_original_by_quadratic_monomials_matrix']))
        cent=[]
        for base in (0,1):
            A=[[values[4*i+2*j+base]for j in range(2)]for i in range(2)];cols=[]
            for at in range(4):
                C=[[F(int(2*i+j==at))for j in range(2)]for i in range(2)];CA=mm(C,A);AC=mm(A,C);cols.append([CA[i][j]-AC[i][j]for i in range(2)for j in range(2)])
            cent+=S.transpose(cols)
        assert cent==mat(c['normal_centralizer_matrix']);KC=mat(c['full_normal_centralizer_kernel']);kernels(cent,KC);S.minor(cent,c['normal_centralizer_rank_certificate'])
        group=3+6+len(KC);assert group==c['actual_connected_group_fiber_dimension']==(10 if c['type']=='I0_three_distinct'else 11)
        output.append(dict(type=c['type'],full_connection_rank8_minor=cm,all15_parallel_kernel_dimension=len(K),parallel_minor=pm,
            full374_by30_rank18_minor=m1,complete_P1_dimension12=True,full350_by9_rank3_minor=m2,fullP2_kernel_six_standard_columns=True,
            all18_source_fields_common_partial_z_and_z_independent=True,whole_rank_at_most2_actualDzero=True,correct_group_fiber=group,old_coarse_upper=c['whole_projective_image_dimension_upper']))
        print(c['type'],'whole18 rank2 PASS',flush=True)
    result=dict(task_id='G',all_passed=True,all_four_whole_sources=output,
        theorem='Concurrent reduced-cubic ranktwo-plane m4 generic-squarefree I0/J0 rank2 sources are empty for rank4. Combined accepted genericj/rank1 emptiness gives allconcurrent m4 squarefree empty.',
        no_CAS_or_producer_imports=True,old_false31_32_preserved=True,
        complements=['repeatedleadingroots','m2/m3','nonreducedR','other ranktwo surfaces','wholeE'],
        inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (ROOT/OUT).write_text(json.dumps(result,indent=2)+'\n');print('PASS',OUT,round(result['seconds'],3),flush=True)
if __name__=='__main__':main()

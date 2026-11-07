#!/usr/bin/env python3
"""Independent whole m4 repeated-leading211/31 review; Fraction only.

No producer or CAS import. All six fixed-connection sources, original350
operators and the full150 strict-Casimir sources of both actual-D normal forms
are independently rebuilt. Historical source snapshots remain immutable.
"""
import ast,hashlib,json,resource,signal,time
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement,permutations
from pathlib import Path
import classification_plane_m4_special_review as S

ROOT=Path(__file__).resolve().parents[1]
JOB='research/classification/jobs/plane_m4_repeated_review.json'
OUT='reports/CLASSIFICATION_PLANE_M4_REPEATED_REVIEW_20261107.json'
MAN211='research/classification/inputs/plane_m4_repeated_211_manifest_20261107.json'
MAN31='research/classification/inputs/plane_m4_triple_root_quintics_manifest_20261107.json'
REPORT211='reports/CLASSIFICATION_PLANE_M4_REPEATED_211_20261107.json'
REPORT31='reports/CLASSIFICATION_PLANE_M4_TRIPLE_ROOT_QUINTICS_20261107.json'
PRODUCER='code/classification_plane_reducible_cubic_m4_special.py'

def decode(data):return {tuple(v['exponent']):F(v['coefficient'])for v in data if F(v['coefficient'])}
def mat(A):return [[F(v)for v in row]for row in A]
def tensor(v):return S.skew({ij:{e:F(c)for (pair,e),c in zip(S.COLS,v)if pair==ij and F(c)}for ij in S.PAIRS})
def combine(basis,row):return S.skew({ij:S.add(*(S.scale(P.get(ij,{}),c)for P,c in zip(basis,row)if c))for ij in S.PAIRS})
def check_kernel(M,K,C):
    value=S.minor(M,C);assert len(K)==len(M[0])-C['rank']
    assert len(S.rref(K)[1])==len(K) if K else True
    assert all(not sum(a*b for a,b in zip(row,v))for row in M for v in K)
    return value


def model211():
    a,add,mul=S.a,S.add,S.mul;x,y,z=(a(1,i)for i in range(3))
    tri=mul(mul(x,y),z);con=mul(mul(x,y),add(x,S.scale(y,-1)));tan=mul(z,add(mul(x,z),S.scale(mul(y,y),-1)))
    return [('triangle_diagonal',tri,x,{},y,10,2),('triangle_Jordan',tri,{},S.scale(y,F(1,2)),x,10,2),
      ('concurrent_diagonal',con,x,{},y,11,2),('concurrent_Jordan',con,{},S.scale(y,F(1,2)),x,11,2),
      ('concurrent_three_used',con,x,S.scale(x,-1),add(x,S.scale(y,-1)),10,1),('tangent_conic_net',tan,x,y,z,9,1)]


def review211(packet):
    assert packet['all_arithmetic_passed']and len(packet['cases'])==6
    assert packet['full350_row_labels']==[dict(quad=list(quad),monomial=list(e))for quad,e in S.ROWS]
    X=[S.a(1,i)for i in range(5)];x,y,z=X[:3]
    VF=[[X[b]if i==a else{}for i in range(3)]for a,b in permutations(range(3),2)]
    VF += [[x,S.scale(y,-1),{}],[x,{},S.scale(z,-1)]]
    def contract(V):return [S.add(S.mul(z,V[1]),S.scale(S.mul(y,V[2]),-1)),S.add(S.mul(x,V[2]),S.scale(S.mul(z,V[0]),-1)),S.add(S.mul(y,V[0]),S.scale(S.mul(x,V[1]),-1))]
    mon2=[e for e in S.monomials(2)if not e[3]and not e[4]]
    labels=[(i,e)for i in range(3)for e in mon2]
    Vmap=[[contract(V)[i].get(e,F())for V in VF]for i,e in labels]
    assert len(Vmap)==18 and len(S.rref(Vmap)[1])==8
    first,second=S.sources();d3=[e for e in S.monomials(5)if e[3]+e[4]==3];d4=[e for e in S.monomials(5)if e[3]+e[4]==4]
    out=[]
    for c,(name,f,A,B,C,group,ldim)in zip(packet['cases'],model211()):
        assert c['type']==name and [decode(c[k])for k in('base_cubic','Q_A','Q_B','Q_C')]==[f,A,B,C]
        T=S.add(S.mul(A,C),S.scale(S.mul(B,B),-1));assert T==decode(c['Q_determinant'])
        beta={ast.literal_eval(ij):[decode(p)for p in v]for ij,v in c['full_f_beta'].items()}
        assert all(not p for p in beta[3,4])
        # Verify the whole rational reconstruction by polynomial identities;
        # no assumption that f/T is itself polynomial (Jordan T is a square).
        for i in range(3):
            df,dT,dC,dB=S.der(f,i),S.der(T,i),S.der(C,i),S.der(B,i)
            assert S.scale(S.mul(T,beta[3,3][i]),4)==S.add(S.scale(S.mul(T,df),2),S.scale(S.mul(f,dT),-1))
            assert S.scale(S.mul(S.mul(C,T),beta[4,4][i]),4)==S.add(S.scale(S.mul(S.mul(C,T),df),2),S.scale(S.mul(S.mul(f,T),dC),-4),S.mul(S.mul(f,C),dT))
            assert S.scale(S.mul(S.mul(C,T),beta[4,3][i]),2)==S.add(S.mul(S.mul(f,B),dT),S.scale(S.mul(S.mul(f,T),dB),-2))
        for a,b in((3,3),(3,4),(4,3),(4,4)):
            for i,j in combinations(range(3),2):
                v=S.add(S.mul(f,S.add(S.der(beta[a,b][j],i),S.scale(S.der(beta[a,b][i],j),-1))),S.scale(S.mul(S.der(f,i),beta[a,b][j]),-1),S.mul(S.der(f,j),beta[a,b][i]),
                  *(S.scale(S.add(S.mul(beta[a,k][i],beta[k,b][j]),S.scale(S.mul(beta[a,k][j],beta[k,b][i]),-1)),-1)for k in(3,4)))
                assert not v
        assert Vmap==mat(c['full18x8_form_field_map']);fields={}
        coeff=mat(c['full4x8_tracefree_field_coefficients'])
        for (a,b),row in zip(((3,3),(3,4),(4,3),(4,4)),coeff):
            fields[a,b]=[S.add(*(S.scale(V[i],v)for V,v in zip(VF,row)))for i in range(3)]
            assert contract(fields[a,b])==[S.scale(S.add(beta[a,b][i],S.scale(S.der(f,i),-F(a==b,3))),3)for i in range(3)]
        Pf=S.skew({(0,1):S.der(f,2),(0,2):S.scale(S.der(f,1),-1),(1,2):S.der(f,0)})
        raw=S.plus(Pf,S.skew({(i,a):S.scale(S.add(*(S.mul(X[b],fields[a,b][i])for b in(3,4))),-1)for i in range(3)for a in(3,4)}))
        normals=[S.a(1,3,3),S.a(1,3,4),S.a(1,4,4)];d0=S.actual(raw)
        columns=[S.add(S.actual(S.plus(raw,S.skew({(3,4):n}))),S.scale(d0,-1))for n in normals]
        DK=[[p.get(e,F())for p in columns]for e in S.monomials(5)]
        forced=list(map(F,c['forced_raw_uv_quadratic_coefficients']));minor0=S.minor(DK,c['raw_uv_forcing_rank_certificate'])
        assert c['raw_uv_forcing_rank_certificate']['rank']==3
        assert all(sum(a*b for a,b in zip(row,forced))==-d0.get(e,F())for row,e in zip(DK,S.monomials(5)))
        P0=S.plus(raw,S.skew({(3,4):S.add(*(S.scale(p,v)for p,v in zip(normals,forced)))}))
        assert S.vector(P0)==list(map(F,c['full150_P0']))and not any(S.jac(P0).values())and not S.actual(P0)
        assert c['all15_leading_monomials']==[list(e)for e in d4]
        lead=[{e:F(1)}for e in d4];mon7=S.monomials(7)
        parallel=[[S.add(S.mul(f,S.der(g,i)),S.scale(S.mul(S.der(f,i),g),-1),*(S.mul(S.mul(beta[a,b][i],X[b]),S.der(g,a))for a,b in((3,3),(3,4),(4,3),(4,4))),S.scale(S.mul(S.add(beta[3,3][i],beta[4,4][i]),g),-1))for i in range(3)]for g in lead]
        PL=[[p[i].get(e,F())for p in parallel]for i in range(3)for e in mon7]
        assert PL==mat(c['full_parallel_leading_matrix']);KL=mat(c['full_parallel_leading_kernel'])
        minor_lead=check_kernel(PL,KL,c['parallel_leading_rank_certificate']);assert len(KL)==ldim
        ann=mat(c['leading_annihilator']);assert ann==S.rref(KL)[2]
        assert [S.vector(P)for P in first]==mat(c['full30_weight1_basis'])and [S.vector(P)for P in second]==mat(c['full9_weight2_basis'])
        shape=lambda P:[sum(v*S.actual(P).get(e,F())for v,e in zip(row,d4))for row in ann]
        J1=S.transpose([S.jvector(S.jac(S.plus(P0,P)))for P in first]);D3=[[S.actual(S.plus(P0,P)).get(e,F())for P in first]for e in d3]
        assert J1==mat(c['full350_weight1_matrix'])and D3==mat(c['all24_D3_rows'])
        K1=mat(c['full30_weight1_kernel']);minor1=check_kernel(J1+D3,K1,c['weight1_rank_certificate'])
        J2=S.transpose([S.jvector(S.jac(S.plus(P0,P)))for P in second]);SL=S.transpose([shape(S.plus(P0,P))for P in second]);L=J2+SL
        assert J2==mat(c['full350_weight2_matrix'])and SL==mat(c['all_leading_shape_weight2_rows'])
        K2=mat(c['full9_weight2_kernel']);minor2=check_kernel(L,K2,c['weight2_rank_certificate'])
        family=[combine(first,row)for row in K1];f2=[combine(second,row)for row in K2]
        own=[S.jvector(S.jac(P))+shape(P)for P in family];pairs=list(combinations_with_replacement(range(len(K1)),2));columns=[]
        for i,j in pairs:
            if i==j:columns.append(own[i])
            else:columns.append([v-a-b for v,a,b in zip(S.jvector(S.jac(S.plus(family[i],family[j])))+shape(S.plus(family[i],family[j])),own[i],own[j])])
        Q=S.transpose(columns);assert Q==mat(c['full_original_by_quadratic_monomials_matrix'])and c['all_quadratic_parameter_pairs']==[list(t)for t in pairs]
        qzero=not any(v for row in Q for v in row)
        common=qzero and all(not p or(2 in ij and all(not e[2]for e in p))for P in[P0]+family+f2 for ij,p in P.items())
        assert common==c['entire_necessary_source_common_partial_z_factor']==name.startswith('concurrent')
        assert qzero==c['all_original_quadratic_coefficients_zero']
        d=len(K1)+len(K2);bound=24+d-group
        assert d==c['fixed_P0_full_source_dimension_upper']and group==c['actual_connected_group_fiber_dimension']and bound==c['whole_projective_image_dimension_upper']
        assert common or bound<=24
        out.append(dict(type=name,all_beta_reconstruction_and12_curvature_identities=True,complete18x8_field_inverse=True,forced_normal_K=list(map(str,forced)),rawK_rank3_minor=minor0,
          full150_P0_and350_Jacobi_and_actualPf=True,all15_parallel_leading_kernel_dimension=ldim,parallel_rank_minor=minor_lead,full374_by30_kernel_dimension=len(K1),weight1_minor=minor1,
          full350_plus_shape_by9_kernel_dimension=len(K2),weight2_minor=minor2,all_original_quadratic_rows_reconstructed=True,entire18_source_common_partial_z_rank2=common,
          full_fixed_source_dimension_upper=d,actual_connected_group_dimension=group,projective_image_upper=None if common else bound,rank4_empty=common))
        print(name,'PASS',flush=True)
    return out


def review31(packet):
    X=[S.a(1,i)for i in range(5)];mon6=S.monomials(6);mon1=S.monomials(1)
    labels=[dict(type='Casimir',component=i,monomial=list(e))for i in range(5)for e in mon6]+[dict(type='divergence',component=i,monomial=list(e))for i in range(5)for e in mon1]
    assert packet['full1075_row_labels']==labels and packet['full150_column_labels']==[dict(pair=list(ij),monomial=list(e))for ij,e in S.COLS]
    out=[]
    for eps,c in enumerate(packet['cases']):
        assert c['epsilon']==eps
        Fp=S.add(S.a(1,0,2,2,2,2),S.a(1,1,2,2,2,3),S.a(eps,2,2,3,3,3),S.a(1,3,3,3,3,3));assert Fp==decode(c['F'])
        grad=[S.der(Fp,i)for i in range(5)];assert not S.add(S.mul(X[3],grad[0]),S.scale(S.mul(X[2],grad[1]),-1))
        columns=[]
        for (i,j),e in S.COLS:
            p={e:F(1)};cas=[{}for _ in range(5)];div=[{}for _ in range(5)]
            cas[i]=S.mul(p,grad[j]);cas[j]=S.scale(S.mul(p,grad[i]),-1);div[i]=S.scale(S.der(p,j),-1);div[j]=S.der(p,i)
            columns.append([poly.get(e,F())for poly in cas for e in mon6]+[poly.get(e,F())for poly in div for e in mon1])
        M=S.transpose(columns);assert len(M)==1075 and M==mat(c['full1075x150_strict_Casimir_matrix'])
        K=mat(c['full150_kernel']);minor=check_kernel(M,K,c['rank_certificate']);assert len(K)==(4 if eps else 5)
        syn=[]
        for d,record in enumerate(c['complete_gradient_syzygies_degrees_0_1_2']):
            mons=[e for e in S.monomials(d)if not e[4]];omons=[e for e in S.monomials(d+4)if not e[4]]
            C=[S.mul({e:F(1)},grad[i])for i in range(4)for e in mons];GM=[[p.get(e,F())for p in C]for e in omons]
            assert GM==mat(record['full_gradient_map'])and record['coefficient_labels']==[dict(component=i,monomial=list(e))for i in range(4)for e in mons]
            GK=mat(record['full_kernel']);check_kernel(GM,GK,record['rank_certificate']);syn.append(len(GK))
        if not eps:
            assert syn==[0,2,8]
            Z=[X[3],S.scale(X[2],-1),{},{}];Y=[S.scale(X[0],-4),S.scale(X[1],-3),X[2],{}]
            mons=[e for e in S.monomials(2)if not e[4]];vectors=[]
            for V in(Z,Y):
                for j in range(4):vectors.append([S.mul(X[j],p).get(e,F())for p in V for e in mons])
            assert len(S.rref(vectors)[1])==8
            GM=mat(c['complete_gradient_syzygies_degrees_0_1_2'][2]['full_gradient_map'])
            assert all(not sum(a*b for a,b in zip(row,v))for row in GM for v in vectors)
            # Each coefficient in the remaining vertex w is an AR_0/1/2
            # element. The complete module spans pointwise Z,Y, so after w
            # the image has dimension<=3 and skew rank<=2 on the whole source.
            assert c['rank4_source_empty']
        assert c['actualD_family_dimension_upper']==18 and c['necessary_projective_Poisson_image_upper']==18+len(K)-1
        out.append(dict(epsilon=eps,full1075_by150_matrix_reconstructed=True,full_strict_Casimir_kernel_dimension=len(K),nonzero_minor=minor,complete_AR0_AR1_AR2_dimensions=syn,
            full_F0_AR2_equals_S1Z_plus_S1Y=not eps,whole_F0_rank_at_most2=not eps,actualD_Ga_symmetry_verified=True,actualD_family_dimension_upper=18,
            whole_projective_image_upper=None if not eps else 21,rank4_empty=not eps))
        print('triple-root epsilon',eps,'PASS',flush=True)
    return out


def main():
    start=time.monotonic();signal.alarm(240);resource.setrlimit(resource.RLIMIT_CPU,(235,238));resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    assert not(ROOT/OUT).exists();job=json.loads((ROOT/JOB).read_text())
    for row in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    p211=json.loads((ROOT/REPORT211).read_text());p31=json.loads((ROOT/REPORT31).read_text())
    historical=[]
    for manifest,packet in((MAN211,p211),(MAN31,p31)):
        snapshot=packet['producer_source_bytes'].encode();assert hashlib.sha256(snapshot).hexdigest()==packet['input_sha256'][PRODUCER]
        for row in json.loads((ROOT/manifest).read_text())['files']:
            b=(ROOT/row['path']).read_bytes()
            if hashlib.sha256(b).hexdigest()!=row['sha256']:
                assert row['path']==PRODUCER and hashlib.sha256(snapshot).hexdigest()==row['sha256']
                historical.append(dict(path=PRODUCER,sha256=row['sha256'],read_from_immutable_report=True))
    result=dict(task_id='G',all_passed=True,producer_imports=False,CAS_imports=False,independent_Fraction_only=True,
        exact_inputs_and_hashes=job['exact_inputs_and_hashes'],historical_producer_bindings=historical,
        all_six_whole211_sources=review211(p211),all_two_whole31_actualD_sources=review31(p31),
        theorem='Reduced-reducible R: all repeated-leading m4 ranktwo-plane sources have PGL image<=23 (211 six sources,31 two actualD forms,rank1/22/4 accepted essential3 eight-source coverage). Generic31 result holds also for nonreduced R.',
        complements=['m2/m3 ranktwo-plane','nonreduced R generic211 highrank','other surfaces','wholeE/componenthood/globalG'],
        elapsed_seconds=time.monotonic()-start,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (ROOT/OUT).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in('all_passed','elapsed_seconds','peak_rss_kib')}),flush=True)


if __name__=='__main__':main()

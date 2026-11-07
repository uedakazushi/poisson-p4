#!/usr/bin/env python3
"""Separate Fraction reconstruction of the three entire m4 special-j sources.

No arithmetic, producer or CAS module is imported. Degree-five divisors,
projective Schouten operators, full Euler-quotient sources and rank certificates
are rebuilt from their definitions. The geometric whole-pencil proof is audited
in the accompanying note; these are whole sources, not point jets.
"""
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement
import hashlib,json,resource,signal,time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
INPUT='reports/CLASSIFICATION_PLANE_REDUCIBLE_CUBIC_M4_SPECIAL_20261107.json'
OUTPUT='reports/CLASSIFICATION_PLANE_REDUCIBLE_CUBIC_M4_SPECIAL_REVIEW_20261107.json'
JOB='reports/CLASSIFICATION_PLANE_REDUCIBLE_CUBIC_M4_SPECIAL_REVIEW_JOB_20261107.json'
PAIRS=list(combinations(range(5),2));TRIPLES=list(combinations(range(5),3));QUADS=list(combinations(range(5),4))


def a(c=1,*inds):
    e=[0]*5
    for k in inds:e[k]+=1
    return {tuple(e):F(c)} if c else {}
def add(*Ps):
    out={}
    for P in Ps:
        for e,c in P.items():out[e]=out.get(e,F(0))+c
    return {e:c for e,c in out.items() if c}
def scale(P,c):return {e:c*v for e,v in P.items() if c*v}
def mul(P,Q):
    out={}
    for e,c in P.items():
        for f,d in Q.items():
            g=tuple(v+w for v,w in zip(e,f));out[g]=out.get(g,F(0))+c*d
    return {e:c for e,c in out.items() if c}
def der(P,k):
    out={}
    for e,c in P.items():
        if e[k]:f=list(e);f[k]-=1;out[tuple(f)]=c*e[k]
    return out

def skew(data):
    out={}
    for (i,j),P in data.items():out[i,j]=P;out[j,i]=scale(P,-1)
    return out

def plus(*Ps):return skew({ij:add(*(P.get(ij,{}) for P in Ps)) for ij in PAIRS})
def wedge(X,Y):return skew({ij:add(mul(X[ij[0]],Y[ij[1]]),scale(mul(X[ij[1]],Y[ij[0]]),-1)) for ij in PAIRS})
def monomials(degree):
    out=[]
    for indices in combinations_with_replacement(range(5),degree):
        e=[0]*5
        for i in indices:e[i]+=1
        out.append(tuple(e))
    return out
COLS=[(ij,e) for ij in PAIRS for e in monomials(2)]
ROWS=[(quad,e) for quad in QUADS for e in monomials(4)]
def vector(P):return [P.get(ij,{}).get(e,F(0)) for ij,e in COLS]
def jvector(J):return [J[quad].get(e,F(0)) for quad,e in ROWS]
def jac(P):
    strict={ijk:add(*(mul(P.get((i,l),{}),der(P.get((j,k),{}),l)) for i,j,k in (ijk,(ijk[1],ijk[2],ijk[0]),(ijk[2],ijk[0],ijk[1])) for l in range(5))) for ijk in TRIPLES}
    return {quad:add(*(scale(mul(a(1,quad[k]),strict[quad[:k]+quad[k+1:]]),(-1)**k) for k in range(4))) for quad in QUADS}
def actual(P):
    out={}
    for at in range(5):
        i,j,k,l=[n for n in range(5) if n!=at]
        pf=add(mul(P.get((i,j),{}),P.get((k,l),{})),scale(mul(P.get((i,k),{}),P.get((j,l),{})),-1),mul(P.get((i,l),{}),P.get((j,k),{})))
        out=add(out,scale(mul(a(1,at),pf),(-1)**at))
    return out

def transpose(M):return list(map(list,zip(*M)))
def mat(M):return [[F(v) for v in row] for row in M]
def rref(M):
    A=[list(map(F,row)) for row in M];n=len(A[0]);piv=[]
    for j in range(n):
        r=len(piv);k=next((k for k in range(r,len(A)) if A[k][j]),None)
        if k is None:continue
        A[r],A[k]=A[k],A[r];d=A[r][j];A[r]=[v/d for v in A[r]]
        for k in range(len(A)):
            if k!=r and A[k][j]:
                d=A[k][j];A[k]=[v-d*w for v,w in zip(A[k],A[r])]
        piv.append(j)
    free=[j for j in range(n) if j not in piv];K=[]
    for j in free:
        v=[F(0)]*n;v[j]=F(1)
        for row,p in zip(A,piv):v[p]=-row[j]
        assert all(not sum(a*b for a,b in zip(row,v)) for row in M);K.append(v)
    return A[:len(piv)],piv,K

def det(M):
    B=[list(map(F,row)) for row in M];out=F(1)
    for j in range(len(B)):
        k=next((k for k in range(j,len(B)) if B[k][j]),None)
        if k is None:return F(0)
        if k!=j:B[k],B[j]=B[j],B[k];out=-out
        pivot=B[j][j];out*=pivot
        for k in range(j+1,len(B)):
            c=B[k][j]/pivot
            if c:
                for l in range(j+1,len(B)):B[k][l]-=c*B[j][l]
    return out

def minor(M,C):
    rr,cc=C['rows'],C['columns'];r=C['rank']
    assert len(rr)==len(cc)==r and len(set(rr))==len(set(cc))==r
    value=det([[M[i][j] for j in cc] for i in rr]);assert value and value==F(C['minor_determinant'])
    return str(value)


def sources():
    P1=[]
    for normal in (3,4):
        for ij in combinations(range(3),2):
            for k in range(3):P1.append(skew({ij:a(1,normal,k)}))
    for normal,ex in [(3,(3,4)),(3,(4,4)),(4,(3,3)),(4,(3,4))]:
        for i in range(3):P1.append(skew({(i,normal):a(1,*ex)}))
    P2=[skew({ij:a(1,*ex)}) for ij in combinations(range(3),2) for ex in [(3,3),(3,4),(4,4)]]
    X=[a(1,i) for i in range(5)];gauges=[]
    for normal in (3,4):
        for i in range(3):
            Y=[{} for _ in range(5)];Y[i]=a(1,normal);gauges.append(wedge(X,Y))
    _,piv,_=rref(transpose([vector(P) for P in P1+gauges]));assert len(piv)==36
    raw=[(ij,e) for ij,e in COLS if e[3]+e[4]-sum(i>=3 for i in ij)==1]
    assert len(raw)==36
    assert all(not v for P in P1+gauges for (ij,e),v in zip(COLS,vector(P)) if (ij,e) not in raw)
    assert len(P2)==9
    return P1,P2


def main():
    begun=time.monotonic();signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,)*2)
    dest=ROOT/OUTPUT
    if dest.exists():raise FileExistsError(dest)
    job=json.loads((ROOT/JOB).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    packet=json.loads((ROOT/INPUT).read_text());assert packet['all_arithmetic_passed']
    P1,P2=sources()
    assert [vector(P) for P in P1]==mat(packet['full30_weight1_basis'])
    assert [vector(P) for P in P2]==mat(packet['full9_weight2_basis'])
    assert packet['full350_row_labels']==[dict(quad=list(quad),monomial=list(e)) for quad,e in ROWS]
    f=a(1,0,1,2);Pf=skew({(0,1):der(f,2),(0,2):scale(der(f,1),-1),(1,2):der(f,0)})
    D3mons=[e for e in monomials(5) if e[3]+e[4]==3]
    D4mons=[e for e in monomials(5) if e[3]+e[4]==4]
    assert len(D3mons)==24 and len(D4mons)==15
    models=[('J0_secant',(F(1,8),F(3,8),F(1,2)),(F(3,8),F(1,8),F(1,2)),[(4,0),(0,1)]),
            ('J0_compression',(F(0),F(1,2),F(1,2)),(F(1,2),F(0),F(1,2)),[(3,0),(1,1)]),
            ('I0_mixed',(F(1,6),F(1,3),F(1,2)),(F(1,2),F(0),F(1,2)),[(4,0),(1,1)])]
    output=[]
    for record,(name,bu,bv,leading) in zip(packet['cases'],models):
        assert record['type']==name and list(bu)==list(map(F,record['beta_u_logarithmic_coefficients'])) and list(bv)==list(map(F,record['beta_v_logarithmic_coefficients']))
        assert sum(bu)==sum(bv)==1
        tests=[]
        for k in range(5):
            ex=tuple(1-(k-1)*u-(3-k)*v for u,v in zip(bu,bv))
            coord=next((i for i in range(3) if ex==tuple(F(i==j) for j in range(3))),None)
            tests.append(dict(normal_u_degree=k,required_base_monomial_exponents=list(map(str,ex)),linear_coefficient_coordinate=coord))
        assert tests==record['all_five_parallel_monomial_tests']
        assert [(r['normal_u_degree'],r['linear_coefficient_coordinate']) for r in tests if r['linear_coefficient_coordinate'] is not None]==sorted(leading)
        fields=[]
        for beta in (bu,bv):
            tau=[v-F(1,3) for v in beta];fields.append((tau[2]-tau[1],tau[0]-tau[2],tau[1]-tau[0]))
        assert [list(v) for v in fields]==mat(record['raw_mixed_diagonal_base_fields'])
        raw=plus(Pf,skew({(i,k):a(-fields[j][i],i,k) for j,k in enumerate((3,4)) for i in range(3)}))
        D=actual(raw);assert set(D)<={(1,1,1,1,1)}
        kuv=-D.get((1,1,1,1,1),F(0))/3
        assert kuv==F(record['forced_raw_normal_uv_coefficient'])
        P0=plus(raw,skew({(3,4):a(kuv,3,4)}));assert vector(P0)==list(map(F,record['full150_P0']))
        assert not any(jac(P0).values()) and not actual(P0)
        # Independent affine-chart derivation of beta and of the forced uv term.
        X=[a(1,i) for i in range(5)];Y=[{},{},{},a(fields[0][0],3),a(fields[1][0],4)]
        gauge=plus(P0,wedge(X,Y))
        assert not gauge[0,3] and not gauge[0,4]
        for beta,k in ((bu,3),(bv,4)):
            assert gauge[1,k]==a(3*beta[2]-1,1,k)
            assert gauge[2,k]==a(1-3*beta[1],2,k)
        affine_uv=3*(bu[1]*bv[2]-bu[2]*bv[1]);assert gauge[3,4]==a(affine_uv,3,4)
        J1=transpose([jvector(jac(plus(P0,P))) for P in P1]);D3=[[actual(plus(P0,P)).get(e,F(0)) for P in P1] for e in D3mons]
        assert J1==mat(record['complete350_weight1_matrix']) and D3==mat(record['all24_D3_rows'])
        M1=J1+D3;K1=mat(record['full30_weight1_kernel']);minor1=minor(M1,record['weight1_rank_certificate'])
        assert len(K1)==6 and all(len(v)==30 and all(not sum(a*b for a,b in zip(row,v)) for row in M1) for v in K1)
        _,kpiv,_=rref(K1);assert len(kpiv)==6
        assert 30-record['weight1_rank_certificate']['rank']==len(K1)
        allowed=[]
        for k,i in leading:
            e=[0]*5;e[i]=1;e[3]=k;e[4]=4-k;allowed.append(tuple(e))
        shapes=[e for e in D4mons if e not in allowed]
        assert [list(e) for e in allowed]==record['allowed_two_leading_degree5_monomials']
        assert [list(e) for e in shapes]==record['all13_excluded_leading_monomials'] and len(shapes)==13
        shape=lambda P:[actual(P).get(e,F(0)) for e in shapes]
        J2=transpose([jvector(jac(plus(P0,P))) for P in P2]);S2=transpose([shape(plus(P0,P)) for P in P2]);L=J2+S2
        assert J2==mat(record['full350_weight2_matrix']) and S2==mat(record['all13_leading_shape_weight2_rows'])
        minor2=minor(L,record['full_weight2_rank_certificate']);_,piv2,K2=rref(L)
        assert len(piv2)==record['full_weight2_rank_certificate']['rank']
        family=[skew({ij:add(*(scale(P.get(ij,{}),c) for P,c in zip(P1,row) if c)) for ij in PAIRS}) for row in K1]
        quadratic_pairs=list(combinations_with_replacement(range(6),2));assert record['all_quadratic_parameter_pairs']==[list(ij) for ij in quadratic_pairs]
        own=[jvector(jac(P))+shape(P) for P in family];quadratic=[]
        for i,j in quadratic_pairs:
            if i==j:quadratic.append(own[i])
            else:quadratic.append([v-a-b for v,a,b in zip(jvector(jac(plus(family[i],family[j])))+shape(plus(family[i],family[j])),own[i],own[j])])
        Q=transpose(quadratic);assert Q==mat(record['full363_by_quadratic_monomials_matrix'])
        # A literal inhomogeneous system L*P2=-Q(P1) is all that the proof
        # requires; extra nonlinear equations can only shrink its fibers.
        k2=len(K2);fixed=len(K1)+k2;image=24+fixed-10
        assert k2==record['weight2_fiber_dimension_upper'] and fixed==record['fixed_P0_full_source_dimension_upper'] and image==record['whole_projective_image_dimension_upper'] and image<=24
        output.append(dict(type=name,all150_P0_and350_Jacobi_and_Pf=True,
            independent_affine_chart_beta_check=True,forced_raw_uv=str(kuv),chart_uv=str(affine_uv),
            full_weight_one_matrix_shape=[374,30],complete_weight_one_kernel=record['full30_weight1_kernel'],rank1=24,nonzero_minor1=minor1,
            full_weight_two_matrix_shape=[363,9],complete_weight_two_kernel=[[str(v) for v in row] for row in K2],rank2=len(piv2),nonzero_minor2=minor2,
            all363_quadratic_rows_and21_monomials=True,quadratic_restriction_not_needed=True,
            fixed_source_dimension_upper=fixed,actual_group_fiber_dimension=10,projective_image_dimension_upper=image))
        print(name,'PASS fixedsource',fixed,'whole image upper',image,flush=True)
    result=dict(task_id='G',all_passed=True,cases=output,complete_weight_one_raw36_Euler6_quotient30=True,
        no_arithmetic_producer_or_CAS_imports=True,
        theorem='All nonconcurrent reduced-cubic ranktwo-plane m4 squarefree special-j strata have whole projective image<=24: rankone/conic branches by accepted scalar source; triangle ranktwo special pencils20/22/21.',
        complements=['concurrent reducedcubic','generic repeated binaryroot','m2/m3','nonreduced cubic','whole E'],
        no_full_nonlinear_source_equality_claim=True,producer_snapshot=packet['producer_source_bytes'],
        inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        seconds=time.monotonic()-begun,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    dest.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',OUTPUT,round(result['seconds'],3),'seconds',result['maximum_RSS_KiB'],'KiB',flush=True)

if __name__=='__main__':main()

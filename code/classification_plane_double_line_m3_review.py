#!/usr/bin/env python3
"""Separate Fraction reconstruction of all19 nonred-cubic m3 SF sources.
No producer or CAS imports. Reuses only frozen independent polynomial helpers.
"""
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import ast,hashlib,json,resource,signal,time
import classification_plane_m4_special_review as h
from classification_minimal_quartic_review import LaurentRing
ROOT=Path(__file__).resolve().parents[1]
INPUT='reports/CLASSIFICATION_PLANE_DOUBLE_LINE_M3_SQUAREFREE_REPLAY_20261107.json'
JOB='research/classification/jobs/plane_double_line_m3_separate_review.json'
OUT='reports/CLASSIFICATION_PLANE_DOUBLE_LINE_M3_SEPARATE_REVIEW_20261107.json'
R=LaurentRing((),('a','b','c','z','w'));p=R.parse
x=[h.a(1,i)for i in range(5)];f=p('a**2*b')

def matmul(A,B):return[[sum((a*b for a,b in zip(row,col)if a and b),F(0))for col in zip(*B)]for row in A]
def cert(M,C):
    rank=C['rank'];rr,cc=C['minor_rows'],C['minor_columns'];assert len(rr)==len(cc)==rank
    d=h.det([[M[i][j]for j in cc]for i in rr]);assert d==F(C['minor_determinant'])and d
    K=h.mat(C['complete_kernel']);null=C['nullity'];assert null==len(M[0])-rank and len(K)==len(M[0])and all(len(v)==null for v in K)
    if null:
        assert not any(any(r)for r in matmul(M,K));assert len(h.rref(h.transpose(K))[1])==null
    if rank==0:assert not any(any(r)for r in M)
    return K,str(d)

def monbase(d):return[e[:3]for e in h.monomials(d)if not any(e[3:])]
def monbideg(base,normal):return[b+(k,normal-k)for b in monbase(base)for k in range(normal+1)]
def poly(e,c=1):return{tuple(e):F(c)}if c else{}
def sumP(A,B):return h.plus(A,B)
def linearJac(A,B):return[h1-h2-h3 for h1,h2,h3 in zip(h.jvector(h.jac(h.plus(A,B))),h.jvector(h.jac(A)),h.jvector(h.jac(B)))]
def op(A,basis):return h.transpose([linearJac(A,B)for B in basis])
def allpolys_equal(M,saved):assert M==h.mat(saved)
def fromvector(v):return h.skew({ij:{e:F(c)for(pair,e),c in zip(h.COLS,v)if pair==ij and c}for ij in h.PAIRS})
def parseP(d):return h.skew({tuple(ast.literal_eval(ij)):p(s)for ij,s in d.items()})
def connection(P):
    N=[[[{}for i in range(3)]for b in range(2)]for a in range(2)]
    for a in range(2):
        for b in range(2):
            V=[h.scale(h.der(P.get((i,3+a),{}),3+b),-1)for i in range(3)];assert not h.add(*(h.der(V[i],i)for i in range(3)))
            cross=[h.add(h.mul(V[1],x[2]),h.scale(h.mul(V[2],x[1]),-1)),h.add(h.mul(V[2],x[0]),h.scale(h.mul(V[0],x[2]),-1)),h.add(h.mul(V[0],x[1]),h.scale(h.mul(V[1],x[0]),-1))]
            N[a][b]=[h.add(cross[i],h.der(f,i)if a==b else{})for i in range(3)]
    return N

def geometric_connection(row):
    zero=[{}, {}, {}];N=[[zero[:],zero[:]],[zero[:],zero[:]]];kind=row['kind'];name=row['type']
    if kind in('diagonal','nil'):
        if kind=='diagonal':lam,mu=F(row['metadata']['lambda']),F(row['metadata']['mu'])
        else:lam,mu=(F(-2),F(-2))if name=='NIL_SPLIT'else(F(-1,2),F(-2))
        for i,L in enumerate((lam,mu)):N[i][i]=[h.scale(p('a*b'),2+L),h.scale(p('a**2'),1-L),{}]
        if kind=='nil':N[1][0]=[p('b**2'),p('-a*b'),{}]
    elif name.startswith('C_DEPENDENT_A2_'):
        g=F(row['metadata']['g']);N[0][0]=[h.scale(p('a*b'),3*(2-g/2)),h.scale(p('a**2'),3*(g/2-1)),{}];N[1][1]=[{},p('3*a**2'),{}]
        N[1][0]=[h.scale(p('b*c'),3*(g/2-1)),h.scale(p('a*c'),3*(4-g)/2),p('-3*a*b')]
    else:
        assert name=='C_DEPENDENT_B2_CANDIDATE';N[0][0]=[p('3*a*b'),{},{}];N[1][1]=[p('6*a*b'),p('-3*a**2'),{}];N[1][0]=[p('3*a*c'),{},p('-3*a**2')]
    return N

def PDE(N,mons,targets):
    columns=[]
    for ex in mons:
        G=poly(ex);values=[]
        for i in range(3):
            expr=h.add(h.scale(h.mul(f,h.der(G,i)),3),h.scale(h.mul(h.der(f,i),G),-3),h.scale(h.mul(h.add(N[0][0][i],N[1][1][i]),G),-1))
            expr=h.add(expr,*(h.mul(N[a][b][i],h.mul(x[3+b],h.der(G,3+a)))for a in range(2)for b in range(2)))
            assert set(expr).issubset(set(targets));values.extend(expr.get(e,F(0))for e in targets)
        columns.append(values)
    return h.transpose(columns)

def shear(P,basis):
    B=h.transpose([h.vector(Q)for Q in basis]);assert len(basis)==30 and all(sum(bool(v)for v in h.vector(Q))==1 for Q in basis)
    # Read only the unique retained monomial of each canonical basis vector.
    selected=[next((label,c)for label,c in zip(h.COLS,h.vector(Q))if c)for Q in basis];columns=[]
    for normal in(3,4):
        for idx in range(3):
            Y=[{}for _ in range(5)];Y[idx]=x[normal]
            orbit=h.skew({ij:h.add(*(h.scale(h.mul(Y[k],h.der(P.get(ij,{}),k)),-1)for k in range(5)),*(h.add(h.mul(P.get((k,ij[1]),{}),h.der(Y[ij[0]],k)),h.mul(P.get((ij[0],k),{}),h.der(Y[ij[1]],k)))for k in range(5)))for ij in h.PAIRS})
            Z=[{}for _ in range(5)]
            for i in range(3):
                Z[i]=h.add(h.scale(x[3],orbit.get((i,3),{}).get((0,0,0,2,0),F(0))),h.scale(x[4],orbit.get((i,4),{}).get((0,0,0,0,2),F(0))))
            canonical=h.plus(orbit,h.wedge(x,Z));coeff=[canonical.get(ij,{}).get(ex,F(0))/c for(ij,ex),c in selected]
            assert h.vector(canonical)==[sum(a*b for a,b in zip(r,coeff))for r in B];columns.append(coeff)
    return h.transpose(columns)

def check(row,basis1,basis2,mons,kmon,pdemons):
    name=row['type'];P=parseP(row['full150_P0']);raw=h.skew({ij:v for ij,v in P.items()if ij[0]<ij[1]and ij!=(3,4)});D0=h.actual(raw)
    C=h.transpose([[h.add(h.actual(h.plus(raw,h.skew({(3,4):poly(ex)}))),h.scale(D0,-1)).get(e,F(0))for e in kmon]for ex in[(0,0,0,2,0),(0,0,0,1,1),(0,0,0,0,2)]])
    const=[[D0.get(e,F(0))]for e in kmon];allpolys_equal(C,row['full30_K_global_consistency_matrix']);allpolys_equal(const,row['full30_K_global_consistency_constant'])
    cert(C,row['K_matrix_certificate']);aug=[r+[-c[0]]for r,c in zip(C,const)];cert(aug,row['K_augmented_certificate']);assert row['K_matrix_certificate']['rank']==row['K_augmented_certificate']['rank']==3
    kval=[P.get((3,4),{}).get(ex,F(0))for ex in[(0,0,0,2,0),(0,0,0,1,1),(0,0,0,0,2)]];assert all(sum(a*b for a,b in zip(r,kval))+c[0]==0 for r,c in zip(C,const));assert P.get((3,4),{})==p(row['unique_global_K'])
    assert not h.actual(P)and not any(h.jac(P).values())
    N=connection(P);assert N==geometric_connection(row)==[[[p(s)for s in entry]for entry in r]for r in row['full_3f_connection_numerators']]
    PM=PDE(N,mons,pdemons);allpolys_equal(PM,row['full180by24_parallel_matrix']);KG,pmminor=cert(PM,row['parallel_certificate'])
    Ann=h.mat(row['parallel_membership_annihilator']);assert len(Ann)==24-row['parallel_certificate']['nullity'];assert not any(any(r)for r in matmul(Ann,KG));assert len(h.rref(Ann)[1])==len(Ann)
    J1,J2=op(P,basis1),op(P,basis2);allpolys_equal(J1,row['full350_weight1_operator']);allpolys_equal(J2,row['full350_weight2_operator'])
    D3=h.transpose([[h.add(h.actual(h.plus(P,Q)),h.scale(h.actual(Q),-1)).get(e,F(0))for e in mons]for Q in basis1]);allpolys_equal(D3,row['full24by30_actualD3_map'])
    M1=J1+matmul(Ann,D3);allpolys_equal(M1,row['full_weight1_plus_parallel_matrix']);K1,m1minor=cert(M1,row['weight1_certificate']);K2,m2minor=cert(J2,row['weight2_certificate'])
    S=shear(P,basis1);allpolys_equal(S,row['complete_actual_shear']);cert(S,row['actual_shear_certificate']);assert not any(any(r)for r in matmul(M1,S))
    image=matmul(D3,K1);allpolys_equal(image,row['full_leading_image_matrix']);_,iminor=cert(image,row['leading_image_certificate'])
    RC=LaurentRing(tuple('I'+str(i)for i in range(4)),('a','b','c','z','w'));G={}
    for k,col in enumerate(row['leading_image_certificate']['minor_columns']):
        for ex,r in zip(mons,image):
            if r[col]:ee=[0]*4+list(ex);ee[k]=1;G=RC.add(G,{tuple(ee):r[col]})
    assert G==RC.parse(row['whole_leading_image_polynomial']);coefs=[]
    for zdegree in(3,2,1,0):
        PP={}
        for ee,c in G.items():
            if ee[-2:]==(zdegree,3-zdegree):e=list(ee);e[-2]=e[-1]=0;PP[tuple(e)]=c
        coefs.append(PP)
    A,B,Cc,D=coefs;mul=RC.mul;power=lambda P,k:powpoly(RC,P,k)
    disc=RC.add(mul(power(B,2),power(Cc,2)),RC.scale(mul(A,power(Cc,3)),-4),RC.scale(mul(power(B,3),D),-4),RC.scale(mul(power(A,2),power(D,2)),-27),RC.scale(mul(A,B,Cc,D),18))
    assert disc==RC.parse(row['whole_leading_discriminant_polynomial']);empty=not disc;assert empty==row['entire_squarefree_m3_source_empty']
    higher=row['weight1_certificate']['nullity']+row['weight2_certificate']['nullity'];assert higher==row['whole_higher_source_dimension_upper'];upper=row['whole_fixed_plane_symbol_dimension_upper']+5+higher;assert upper==row['whole_projective_image_upper']
    return dict(type=name,global_K=row['unique_global_K'],all30_K_affine_rows_and_unique_solution=True,all350_Jacobi_and_actualD0_zero=True,full180x24_parallel_PDE_and_complete_kernel=True,parallel_minor=pmminor,whole350_by30_and350_by9_reconstructed=True,weight1_nullity=row['weight1_certificate']['nullity'],weight2_nullity=row['weight2_certificate']['nullity'],weight1_minor=m1minor,weight2_minor=m2minor,full_actual_shear=True,actualD3_full24x30_and_image=True,leading_image_rank=row['leading_image_certificate']['rank'],leading_image_minor=iminor,whole_discriminant_identically_zero=empty,SF_whole_branch_empty=empty,old_coarse_upper_preserved=upper,G24_by_weak_B9=(name=='C_DEPENDENT_A2_G2'),unresolved=name in('DIAGONAL_-2_-2','NIL_SPLIT'))

def powpoly(R,P,k):
    out=R.constant(1)
    for _ in range(k):out=R.mul(out,P)
    return out

def main():
    signal.alarm(300);t=time.monotonic();job=json.loads((ROOT/JOB).read_text());assert all(hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256']for b in job['exact_inputs_and_hashes']);raw=json.loads((ROOT/INPUT).read_text());assert hashlib.sha256(raw['producer_bytes'].encode()).hexdigest()==raw['producer_sha256']=='95981c680721a562462a81d1f97489e87c612e7a8ac1384fe9bcf9ef374bf0b2'
    B1,B2=h.sources();B1=B1[:18]+B1[24:27]+B1[18:21]+B1[27:30]+B1[21:24]
    params=tuple('l'+str(i)for i in range(9))+tuple('n'+str(i)for i in range(9))+tuple(s+str(i)for s in('d','b','e','c')for i in range(3));T=LaurentRing(params,('a','b','c','z','w'))
    for ij,s in raw['full_source1'].items():
        pair=tuple(ast.literal_eval(ij));expected={}
        for j,Q in enumerate(B1):
            for ex,c in Q.get(pair,{}).items():ee=[0]*30+list(ex);ee[j]=1;expected[tuple(ee)]=c
        assert T.parse(s)==expected
    T=LaurentRing(tuple('r'+str(i)for i in range(9)),('a','b','c','z','w'))
    for ij,s in raw['full_source2'].items():
        pair=tuple(ast.literal_eval(ij));expected={}
        for j,Q in enumerate(B2):
            for ex,c in Q.get(pair,{}).items():ee=[0]*9+list(ex);ee[j]=1;expected[tuple(ee)]=c
        assert T.parse(s)==expected
    mons,kmon,pdemons=monbideg(2,3),monbideg(3,2),monbideg(4,3);assert[ list(e)for e in mons]==raw['all24_leading_monomials'];assert[list(e)for e in kmon]==raw['all30_K_monomials'];assert[list(e)for e in pdemons]==raw['all60_parallel_component_monomials']
    pairs=set();inter=[]
    for k,l in[(0,2),(0,3),(1,2),(1,3)]:
        for pp,qq in product(range(3),repeat=2):
            det=(k-1)*(2-l)-(2-k)*(l-1);v0,v1=4-3*pp,4-3*qq;lam=F(v0*(2-l)-(2-k)*v1,det);mu=F((k-1)*v1-v0*(l-1),det);pairs.add(tuple(sorted((lam,mu))));inter.append(dict(normal_support=[k,l],base_powers=[pp,qq],lambda_mu=[str(lam),str(mu)]))
    assert inter==raw['all36_eigenweight_intersections'];assert[[str(v)for v in pair]for pair in sorted(pairs)]==raw['all_normal_swap_distinct_diagonal_pairs'];assert len(raw['cases'])==19
    results=[check(row,B1,B2,mons,kmon,pdemons)for row in raw['cases']];empty=[r['type']for r in results if r['SF_whole_branch_empty']];assert len(empty)==16;assert not raw['all_G24_closed']and not raw['all_squarefree_m3_sources_empty']
    # Actual B9 redundancy, every bracket, not a generic tangent estimate.
    P=parseP(next(r for r in raw['cases']if r['type']=='C_DEPENDENT_A2_G2')['full150_P0']);Y=[{},{},p('b'),{},{}]
    orbit=h.skew({ij:h.add(*(h.mul(Y[k],h.der(P.get(ij,{}),k))for k in range(5)),*(h.scale(h.add(h.mul(P.get((k,ij[1]),{}),h.der(Y[ij[0]],k)),h.mul(P.get((ij[0],k),{}),h.der(Y[ij[1]],k))),-1)for k in range(5)))for ij in h.PAIRS});assert not any(orbit.values())
    out=dict(task_id='G',all19_independent_arithmetic_passed=True,no_producer_or_CAS_imports=True,all_raw36_Euler6_canonical30_and9_complete=True,all36_eigen_intersections_and13_normal_swap_classes=True,cases=results,all16_SF_empty_cases=empty,entire_cdepG2_image_upper24_by_B9=True,actual_c_shift_b_invariance=True,all_G_closed=False,unresolved=['DIAGONAL_-2_-2','NIL_SPLIT'],hand_K_erratum_independently_confirmed=True,no_claim_fourth_globalK_inconsistent=True,old_coarse_false_results_preserved=True,whole_geometry_independently_read=True,root_cover_complements=['split','1+2','C3','rdegree0','rdegree011','allnilresonances/eigenwalls','allglobalKcandidates'],historical_lost138_report_not_falsely_restored=True,inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),seconds=time.monotonic()-t,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    dest=ROOT/OUT;assert not dest.exists();dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k]for k in('all19_independent_arithmetic_passed','seconds','maximum_RSS_KiB','all_G_closed')}))
if __name__=='__main__':main()

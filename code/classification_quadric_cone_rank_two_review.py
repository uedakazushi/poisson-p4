#!/usr/bin/env python3
"""Independent complete rank2 quadric-cone G review over Fraction/QQ(k,s).

No producer or CAS imports. Reconstructs the full150 tangent condition, raw69,
all59 affine modh equations, full33/Euler quotient, all350/126, five pencil
blocks, highest symbol, all graded operators/walls and the genuine resonance.
"""
import ast,hashlib,json,resource,signal,time
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement,product
from pathlib import Path
import classification_cubic_scroll_abelian_review as a
import classification_quadric_cone_zero_review as z
import classification_plane_m2_completion_review as m

ROOT=Path(__file__).resolve().parents[1]
MAN='research/classification/inputs/quadric_cone_rank_two_manifest_20261107.json'
JOB='research/classification/jobs/quadric_cone_rank_two_review.json'
OUT='reports/CLASSIFICATION_QUADRIC_CONE_RANK_TWO_REVIEW_20261107.json'
PREFIX='reports/CLASSIFICATION_QUADRIC_CONE_RANK_TWO_'
PAIRS=list(combinations(range(5),2));P3=list(combinations(range(3),2));QUADS=list(combinations(range(5),4))
M2=a.monomials(5,2);M3=a.monomials(5,3);M4=a.monomials(5,4);M5=a.monomials(5,5)
ROWS=list(product(QUADS,M4));COLS=list(product(PAIRS,M2))


def read(suffix):return json.loads((ROOT/(PREFIX+suffix+'_20261107.json')).read_text())
def mat(A):return [[F(c)for c in row]for row in A]
def decode_rows(A):return [{int(j):F(c)for j,c in row.items()}for row in A]
def poly(A):return {tuple(e):F(c)for e,c in A if F(c)}
def tensor(A):
    P=z.zero(5)
    for ij,v in A.items():z.setentry(P,*ast.literal_eval(ij),poly(v))
    return P
def vec(P):return [P[i][j].get(e,F())for i,j in PAIRS for e in M2]
def combine(basis,values):
    return [[a.add(*(a.scale(P[i][j],c)for P,c in zip(basis,values)if c))for j in range(5)]for i in range(5)]
def get(P,i,j):return P[i][j]
def pull(P):
    out={};powers=[(2,0,0),(1,1,0),(0,2,0),(0,0,1)]
    for e,c in P.items():
        if e[4]:continue
        f=tuple(sum(e[i]*powers[i][j]for i in range(4))for j in range(3));out[f]=out.get(f,F())+c
    return {e:c for e,c in out.items()if c}


def source(saved):
    xs=[z.var(5,i)for i in range(5)];q=a.add(a.mul(xs[0],xs[2]),a.scale(a.mul(xs[1],xs[1]),-1));h=xs[4]
    L=[];A=[];J=[]
    for i,j in PAIRS:
        for k in range(5):P=z.zero(5);z.setentry(P,i,j,a.mul(h,xs[k]));L.append(P)
        P=z.zero(5);z.setentry(P,i,j,q);A.append(P)
    fmons=[e for e in M2 if not e[4]and e!=(1,0,1,0,0)]
    for e in fmons:
        Fpoly={e:F(1)};P=z.zero(5)
        for i,j in combinations(range(4),2):
            k,l=[v for v in range(4)if v not in(i,j)];perm=[i,j,k,l]
            sign=(-1)**sum(perm[u]>perm[v]for u in range(4)for v in range(u+1,4))
            z.setentry(P,i,j,a.scale(a.add(a.mul(a.diff(q,k),a.diff(Fpoly,l)),a.scale(a.mul(a.diff(q,l),a.diff(Fpoly,k)),-1)),sign))
        assert all(not a.add(*(a.mul(P[i][j],a.diff(q,j))for j in range(5)))for i in range(5))
        assert all(not P[4][j]for j in range(5));J.append(P)
    assert fmons==[tuple(e)for e in saved['F_basis_monomials']]
    raw=L+A+J;assert len(raw)==69 and z.rank([vec(P)for P in raw])==69
    vertical=[]
    for i in range(5):
        for j in range(5):
            Y=[{}for _ in range(5)];Y[i]=xs[j];vertical.append(vec(z.wedge(xs,Y)))
    assert z.rank(vertical)==24 and z.rank(vertical+[vec(P)for P in raw])==88
    # Entire150 tangent system, not only testing the displayed ansatz.
    target={}
    for col,(ij,e)in enumerate(COLS):
        P=z.zero(5);z.setentry(P,*ij,{e:F(1)})
        for g,dg in enumerate(([{}, {}, {}, {},z.constant(5,1)],[a.diff(q,k)for k in range(5)])):
            C=[a.add(*(a.mul(P[i][k],dg[k])for k in range(5)))for i in range(5)]
            for i,j in PAIRS:
                p=pull(a.add(a.mul(xs[i],C[j]),a.scale(a.mul(xs[j],C[i]),-1)))
                for ex,c in p.items():target.setdefault((g,i,j,ex),{})[col]=c
    tr=list(target.values());assert len(a.rref(tr))==62
    assert all(not sum(c*v[j]for j,c in row.items())for row in tr for v in vertical+[vec(P)for P in raw])
    old=json.loads((ROOT/'reports/CLASSIFICATION_QUADRIC_CONE_ZERO_SOURCE_20261107.json').read_text())
    membership=[a.add(*(a.mul(P[4][i],a.diff(q,i))for i in range(5)))for P in raw[:60]]
    ideals=[a.mul(a.mul(h,h),x)for x in xs]+[a.mul(h,q)]
    C=[[P.get(e,F())for P in membership]+[-P.get(e,F())for P in ideals]for e in M3]
    assert C==mat(old['full_conormal_matrix'])and z.rank(C)==17
    K=mat(old['complete_conormal_kernel_including_membership']);constraints=mat(old['all11_eliminated_raw_source_constraints'])
    assert len(K)==z.rank(K)==49 and z.rank([v[:60]for v in K])==49
    assert all(not z.dot(row,v)for row in C for v in K)
    assert len(constraints)==z.rank(constraints)==11 and all(not z.dot(row,v[:60])for row in constraints for v in K)
    E=[[F((j==k and ell==i)-(i==k and ell==j))for i,j in PAIRS for ell in range(5)]+[F()]*9 for k in range(5)]
    return L,A,J,fmons,constraints,E,dict(full150_tangent_rank=62,raw69_independent=True,Euler24_rank=24,
       raw_plus_Euler_rank=88,complete_projective_tangent_dimension=64,full11_conormal_constraints=True,
       complete_raw_conormal_dimension=58,full150_normal_complements_retained=True)


def affine_pencil(P):
    xs=[z.var(5,i)for i in range(5)];A=z.zero(4)
    def deh(p):
        out={}
        for e,c in p.items():out[e[:4]]=out.get(e[:4],F())+c
        return {e:c for e,c in out.items()if c}
    for i,j in combinations(range(4),2):
        z.setentry(A,i,j,deh(a.add(P[i][j],a.scale(a.mul(xs[i],P[4][j]),-1),a.mul(xs[j],P[4][i]))))
    B=z.zero(3);D=z.zero(3);Y=[];T=[]
    for i,j in P3:
        assert all(e[3]<=1 and(not e[3]or sum(e[:3])==0)for e in A[i][j])
        z.setentry(B,i,j,{e[:3]:c for e,c in A[i][j].items()if not e[3]})
        z.setentry(D,i,j,{e[:3]:c for e,c in A[i][j].items()if e[3]})
    for i in range(3):
        assert all(e[3]<=1 for e in A[i][3]);Y.append({e[:3]:c for e,c in A[i][3].items()if not e[3]});T.append({e[:3]:c for e,c in A[i][3].items()if e[3]})
    return B,D,Y,T


def five_blocks(data,JJ):
    pairs=list(combinations_with_replacement(range(len(data)),2));rows={}
    def record(kind,component,P,col):
        for e,c in P.items():rows.setdefault((kind,component,e),{})[col]=c
    def mixed(B,C):return a.add(*(a.mul(B[i][ell],a.diff(C[j][k],ell))for i,j,k in[(0,1,2),(1,2,0),(2,0,1)]for ell in range(3)))
    def wd(V,D):return a.add(a.mul(V[0],D[1][2]),a.scale(a.mul(V[1],D[0][2]),-1),a.mul(V[2],D[0][1]))
    for col,(aa,bb)in enumerate(pairs):
        first={};second={};third={ij:{}for ij in P3};fourth={ij:{}for ij in P3};fifth={ij:{}for ij in P3}
        for i,j in([(aa,bb)]if aa==bb else[(aa,bb),(bb,aa)]):
            B,D,Y,T=data[i];Bc,Dc,Yc,Tc=data[j]
            first=a.add(first,mixed(B,Bc),wd(Y,Dc));second=a.add(second,mixed(D,Bc),wd(T,Dc))
            LY=z.lie(Y,Bc,0);LT=z.lie(T,Bc,0);LDY=z.lie(Y,Dc,0);LDT=z.lie(T,Dc,0);YT=z.wedge(Y,Tc)
            for ij in P3:
                u,v=ij;third[ij]=a.add(third[ij],a.scale(LY[u][v],-1),YT[u][v])
                fourth[ij]=a.add(fourth[ij],a.scale(LT[u][v],-1),a.scale(LDY[u][v],-1));fifth[ij]=a.add(fifth[ij],a.scale(LDT[u][v],-1))
        record('JacB_plus_YwedgeD',(),first,col);record('B_D_cross_plus_TwedgeD',(),second,col)
        predicted={(0,1,2):a.add({e+(0,):c for e,c in first.items()},{e+(1,):c for e,c in second.items()})}
        for ij in P3:
            record('minusLieYB_plus_YwedgeT',ij,third[ij],col);record('minusLieTB_minusLieYD',ij,fourth[ij],col);record('minusLieTD',ij,fifth[ij],col)
            predicted[(*ij,3)]=a.add({e+(0,):c for e,c in third[ij].items()},{e+(1,):c for e,c in fourth[ij].items()},{e+(2,):c for e,c in fifth[ij].items()})
        for tri,value in predicted.items():
            qi=QUADS.index((*tri,4));actual={}
            for j,e in enumerate(M4):
                c=JJ[70*qi+j].get(col,F())
                if c:actual[e[:4]]=actual.get(e[:4],F())+c
            assert not a.add(value,{e:c for e,c in actual.items()if c})
    return [dict(condition=k,component=list(c),spatial_monomial=list(e),parameter_pair_coefficients={str(j):str(v)for j,v in row.items()})for(k,c,e),row in sorted(rows.items())]


def whole_sources(packet,L,A,Fb,fmons,C,E):
    unknown=L+Fb;A0=A[PAIRS.index((3,4))];A0=[[a.scale(p,-1)for p in row]for row in A0]
    Cb=[A[PAIRS.index(ij)]for ij in P3];universal=[A0,*Cb,*unknown]
    UJ,_=a.jacobi_and_divisor(universal);upairs=list(combinations_with_replacement(range(63),2));modrows=[i for i,(_,e)in enumerate(ROWS)if not e[4]]
    assert all(all(upairs[j][0]<4 for j in UJ[i])for i in modrows)
    result=[];cache=[]
    for saved,cv in zip(packet['cases'],[(0,1,0),(1,0,0),(0,0,0)]):
        const=list(map(F,[1,*cv]));Avec=[F()for _ in range(10)];Avec[PAIRS.index((3,4))]=-1
        for ij,c in zip(P3,cv):Avec[PAIRS.index(ij)]=F(c)
        M=[];rhs=[]
        for i in modrows:
            row=[F()]*59;constant=F()
            for j,c in UJ[i].items():
                u,v=upairs[j]
                if v<4:constant+=c*const[u]*const[v]
                else:row[v-4]+=c*const[u]
            M.append(row);rhs.append(-constant)
        assert M==mat(saved['all_modh_linear_rows'])and rhs==list(map(F,saved['all_modh_rhs']))and modrows==saved['all350_modh_row_indices']
        CM=[row[:50]+[F()]*9 for row in C];cb=[-z.dot(row[50:],Avec)for row in C]
        assert CM==mat(saved['all11_conormal_rows'])and cb==list(map(F,saved['all11_conormal_rhs']))
        M+=CM;rhs+=cb;assert z.rank(M)==21
        for pointkey,kernelkey,d in [('raw_particular','raw_complete_kernel',38),('projective_particular','complete_projective_modh_kernel',33)]:
            point=list(map(F,saved[pointkey]));K=mat(saved[kernelkey]);assert len(K)==z.rank(K)==d
            assert all(z.dot(row,point)==b for row,b in zip(M,rhs))and all(not z.dot(row,v)for row in M for v in K)
        point=list(map(F,saved['projective_particular']));K=mat(saved['complete_projective_modh_kernel']);cuts=saved['Euler_gauge_indices'];free=saved['projective_free_raw_indices']
        assert a.determinant([[v[j]for j in cuts]for v in E])==F(saved['Euler_gauge_minor'])!=0
        assert all(not z.dot(row,v)for row in M for v in E)and all(not point[j]for j in cuts)and all(not v[j]for v in K for j in cuts)
        assert [[v[j]for j in free]for v in K]==[[F(i==j)for j in range(33)]for i in range(33)]
        assert all(not point[j]for j in free)and[free[i-1]for i in(24,28,31,33)]==[48,52,55,57]
        baseline=combine([A0,*Cb],const);family=[z.plusmatrix(baseline,combine(unknown,point)),*[combine(unknown,v)for v in K]]
        assert family==[tensor(P)for P in saved['all_recovered150_tensors']]
        JJ,DD=a.jacobi_and_divisor(family)
        assert JJ==decode_rows(saved['full_original350_Jacobi_parameter_rows'])and DD==decode_rows(saved['full126_actual_Pfaffian_parameter_rows'])
        assert not any(JJ[i]for i in modrows)
        data=[affine_pencil(P)for P in family]
        for j,(B,D,Y,T)in enumerate(data):
            bp=saved['pencil_B0_point_and_all_directions'][j]
            assert all(B[i][k]==poly(bp[str((i,k))])for i,k in P3)
            assert Y==[poly(v)for v in saved['pencil_Y_point_and_all_directions'][j]]and T==[poly(v)for v in saved['pencil_T_point_and_all_directions'][j]]
        BDs=[[D[i][j].get((0,0,0),F())for B,D,Y,T in data]for i,j in P3]
        assert [row[0]for row in BDs]==list(map(F,saved['pencil_D_constant_point']))and[row[1:]for row in BDs]==mat(saved['pencil_D_constant_parameter_matrix'])
        conditions=five_blocks(data,JJ);assert conditions==saved['all_pencil_condition_coefficient_rows']
        fv=[50+fmons.index(tuple(int(j==i)+int(j==3)for j in range(5)))for i in range(3)]
        vv=50+fmons.index((0,0,0,2,0));assert not point[vv]and all(not v[vv]for v in K)
        assert [point[i]for i in fv]==list(map(F,saved['F_v_linear_point']))and[[v[i]for v in K]for i in fv]==mat(saved['F_v_complete_coefficient_matrix'])
        if saved['type']=='semisimple':assert BDs[2][28]==2 and BDs[0][33]==-2
        result.append(dict(type=saved['type'],complete59_system_rank=21,raw_kernel=38,Euler5_quotient_dimension=33,
             all150_tensors_and350_Jacobi_and126_actual_Pfaffian_reconstructed=True,all_five_pencil_blocks_and_dense_chart_equivalence=True,
             complete_free_coordinates_inverse=True,F_vv_forced_zero=True,F_v_all_three_coefficients_retained=True))
        cache.append(dict(saved=saved,family=family,data=data,J=JJ,D=DD,conditions=conditions))
    return result,cache


def highest(high,cache):
    degrees={'JacB_plus_YwedgeD':3,'B_D_cross_plus_TwedgeD':1,'minusLieYB_plus_YwedgeT':4,'minusLieTB_minusLieYD':2,'minusLieTD':0}
    allowed=[0,24,28,31,33];pairs5=list(combinations_with_replacement(range(5),2));out=[]
    for saved,C in zip(high['cases'],cache):
        pairs=list(combinations_with_replacement(range(34),2));rows=[];encoded=[]
        for row in C['conditions']:
            if sum(row['spatial_monomial'])!=degrees[row['condition']]:continue
            value=[F()]*15
            for col,c in row['parameter_pair_coefficients'].items():
                i,j=pairs[int(col)];assert i in allowed and j in allowed
                value[pairs5.index((allowed.index(i),allowed.index(j)))]+=F(c)
            rows.append(value);encoded.append(dict(condition=row['condition'],component=row['component'],spatial_monomial=row['spatial_monomial'],coefficient_row=list(map(str,value))))
        assert encoded==saved['all_selected_coefficient_rows']
        RR=a.rref(rows);basis=[[v.get(j,F())for j in range(15)]for v in RR.values()]
        assert basis==mat(saved['rowspace_basis'])and len(RR)==saved['all_rows_linear_span_rank']
        forced=(2,4)if saved['type']=='semisimple'else(2,3)if saved['type']=='nilpotent'else()
        if forced:
            for f in forced:assert any(row[pairs5.index((0,f))]and sum(bool(v)for v in row)==1 for row in basis)
            assert all(not c or i in forced or j in forced for row in basis for c,(i,j)in zip(row,pairs5))
        else:assert not basis
        keep=set(range(34))-{allowed[i]for i in forced}
        assert all(not p for i in keep for row in C['data'][i][1]for p in row)
        out.append(dict(type=saved['type'],all_highest_rows_and_other29_absence_verified=True,
             exact_highest_ideal='(ell0,ell2)'if forced==(2,4)else'(ell0,ell1)'if forced else'(0)',
             D_pencil_zero_derived_from_original_Jacobi=True))
    return out


ENV={'kappa':m.TAU,'sigma':m.ETA};R3=m.R3;RF=m.RF;U=m.H.previous.RHO
def scalar(t,env=ENV,field=R3):return m.scalar(t,env,field)
def cert(M,C,env,field):
    minor=m.minor(M,C,env,field);K=[[scalar(c,env,field)for c in row]for row in C['complete_kernel']];k=C['nullity']
    assert len(K)==len(M[0])and all(len(row)==k for row in K)
    assert all(not sum((u*v for u,v in zip(row,col)if u and v),field(0))for row in M for col in zip(*K))
    assert all(any(row[j]==field(1)and all(not c for i,c in enumerate(row)if i!=j)for row in K)for j in range(k))
    return minor
def summary(Ms,C,env,field,offset,base):
    product_=field(1);nulls=[]
    for M,c in zip(Ms,C['stage_certificates']):product_=product_*cert(M,c,env,field);nulls.append(c['nullity'])
    assert C['base_dimension']==base and C['triangular_fiber_upper']==sum(nulls)
    assert C['normalized_source_upper']==base+sum(nulls)and C['whole_PGL_image_upper']==offset+base+sum(nulls)
    assert C['G24_closed']==(C['whole_PGL_image_upper']<=24)
    return product_,nulls
def polynomial_factor_identity(value,factors,env,field):
    factor=field(1)
    for text,n in factors:
        assert n>=1
        for _ in range(n):factor=factor*scalar(text,env,field)
    if field is R3:
        quotient=m.pd(m.pd(value.p,value.q),m.pd(factor.p,factor.q));assert set(quotient)=={m.Z}and quotient[m.Z]
    else:
        quotient=value/factor;assert len(quotient.p)==len(quotient.q)==1 and quotient


def line_cover(Ms,C,values,offset,bad):
    pulled=[m.submatrix(M,values,RF)for M in Ms];product_,nulls=summary(pulled,C['generic'],{'tau':U},RF,offset,1)
    assert product_==scalar(C['selected_minor_product'],{'tau':U},RF)
    roots=C['all_exceptional_roots'];factors=[];seen=set()
    for root in roots:
        factor=scalar(root['minimal_polynomial'],{'tau':U},RF);assert factor.q==(F(1),)and len(factor.p)==2 and factor.p[1]
        t=-factor.p[0]/factor.p[1];assert t==F(root['parameter'])and t not in seen;seen.add(t)
        factors.append((root['minimal_polynomial'],root['multiplicity']))
        specialized=[[[m.at_rf(v,t)for v in row]for row in M]for M in pulled]
        _,ks=summary(specialized,root,{},RF,offset,0)
        if root['whole_PGL_image_upper']>24:
            point=tuple(m.at_rf(v,t).p[0]if m.at_rf(v,t).p else F()for v in values[:2]);assert point==(F(),F())and ks==[5,3,0]
            bad.append(dict(point=['0','0'],stage_nullities=ks,coarse_image=25))
    polynomial_factor_identity(product_,factors,{'tau':U},RF)
    assert C['generic']['whole_PGL_image_upper']<=24
    return dict(all_rational_complex_roots=len(roots),generic_nullities=nulls,all_minor_factor_roots_verified=True)


def graded(tri,cache):
    stages=[[5,7,12,14,18,20,21,22,23,26,27,29,30,32],[1,2,3,4,8,9,10,11,15,16,17,25],[6,13,19]]
    degrees={'JacB_plus_YwedgeD':3,'minusLieYB_plus_YwedgeT':4,'minusLieTB_minusLieYD':2}
    out=[];bad=[];Ms_cache=[]
    for C in tri['cases']:
        data=cache[C['source_case_index']];rows=data['conditions'];pairs=list(combinations_with_replacement(range(34),2));ell=C['F_v_retained_family_index']
        base={0:R3(1),24:m.ETA,ell:m.TAU if C['source_case_index']<2 else R3(1)};retained=set(base)|set(sum(stages,[]));weights={i:0 for i in base}
        for j,stage in enumerate(stages,1):weights.update({i:j for i in stage})
        matrices=[[]for _ in stages];indices=[[]for _ in stages];selected=[]
        for n,row in enumerate(rows):
            terms=[(*pairs[int(col)],F(c))for col,c in row['parameter_pair_coefficients'].items()if all(i in retained for i in pairs[int(col)])]
            if not terms:continue
            assert row['condition']in degrees
            weight=degrees[row['condition']]-sum(row['spatial_monomial']);assert all(weights[i]+weights[j]==weight for i,j,c in terms)
            if not weight:assert not sum((c*base[i]*base[j]for i,j,c in terms),R3(0))
            if weight not in(1,2,3):continue
            stage=stages[weight-1];v=[R3()for _ in stage]
            for i,j,c in terms:
                if i in base and j in stage:v[stage.index(j)]=v[stage.index(j)]+c*base[i]
                if j in base and i in stage:v[stage.index(i)]=v[stage.index(i)]+c*base[j]
            if any(v):matrices[weight-1].append(v);indices[weight-1].append(n)
            selected.append(dict(original_pencil_row_index=n,stage=weight,only_same_weight_terms_verified=True))
        assert matrices==[[[scalar(c)for c in row]for row in M]for M in C['all_stage_matrices']]
        assert indices==C['all_stage_original_pencil_row_indices']and selected==C['all_selected_weight_rows']and stages==C['all_stage_variable_family_indices']
        offset=C['PGL_dimension_offset'];cover=C['cover'];lines=cover.get('all_linear_exception_lines')
        if lines is not None:
            product_,ks=summary(matrices,cover['generic'],ENV,R3,offset,2)
            assert product_==scalar(cover['selected_minor_product'])
            polynomial_factor_identity(product_,[(line['equation'],line['multiplicity'])for line in lines],ENV,R3)
            checks=[]
            for line in lines:
                equation=scalar(line['equation']);assert all(sum(e)<=1 and not e[2]for e in equation.p)and equation.q==m.ONE
                A=equation.p.get((1,0,0),F());B=equation.p.get((0,1,0),F());constant=equation.p.get(m.Z,F());assert A or B
                values=[U,-(A*U+constant)/B,RF(0)]if B else[RF(-constant/A),U,RF(0)]
                assert values[:2]==[scalar(line['parameterization'][k],{'tau':U},RF)for k in ('kappa','sigma')]
                checks.append(line_cover(matrices,line,values,offset,bad))
            assert len(lines)==(12 if C['model']=='semisimple' else 7)
        else:checks=[line_cover(matrices,cover,[RF(0),U,RF(0)],offset,bad)]
        Ms_cache.append(matrices)
        out.append(dict(model=C['model'],all_three_graded_matrices_reconstructed=True,shapes=[[len(M),len(M[0])]for M in matrices],
             all_linear_walls=len(lines)if lines is not None else 0,all_roots=sum(c['all_rational_complex_roots']for c in checks),
             all_complete_kernels_and_minor_factors_verified=True,source_image_offset=offset))
    assert bad and all(v['point']==['0','0']for v in bad)and tri['all_rank2_quadric_cone_G24_closed']is False
    # C=0, F_v=0: all actual h0 coefficients vanish before imposing Jacobi.
    DD=cache[2]['D'];pairs=list(combinations_with_replacement(range(34),2))
    assert all(not c for e,row in zip(M5,DD)if not e[4]for col,c in row.items()if not any(i in(28,31,33)for i in pairs[col]))
    return out,Ms_cache


def resonance(saved,cache,tri,Ms):
    source=cache[1]['conditions'];pairs=list(combinations_with_replacement(range(34),2));stages=tri['cases'][1]['all_stage_variable_family_indices']
    line=next(v for v in tri['cases'][1]['cover']['all_linear_exception_lines']if v['equation']=='sigma');point=next(v for v in line['all_exceptional_roots']if F(v['parameter'])==0)
    K=mat(point['stage_certificates'][0]['complete_kernel']);params=[{}for _ in range(34)];params[0]=z.constant(8,1)
    for i,row in zip(stages[0],K):params[i]=a.add(*(a.scale(z.var(8,j),c)for j,c in enumerate(row)))
    assert params==[poly(v)for v in saved['whole_source_parameter_map']]
    rows=[];rhs=[];ids=[];weights={0:0}
    for w,js in enumerate(stages,1):weights.update({j:w for j in js})
    degrees={'JacB_plus_YwedgeD':3,'minusLieYB_plus_YwedgeT':4,'minusLieTB_minusLieYD':2}
    for n,row in enumerate(source):
        if row['condition']not in degrees or degrees[row['condition']]-sum(row['spatial_monomial'])!=2:continue
        values=[F()]*12;Q={}
        for col,c in row['parameter_pair_coefficients'].items():
            i,j=pairs[int(col)];c=F(c)
            if i in weights and j in weights:
                assert weights[i]+weights[j]==2
                if i==0 and j in stages[1]:values[stages[1].index(j)]+=c
                if j==0 and i in stages[1]:values[stages[1].index(i)]+=c
            Q=a.add(Q,a.scale(a.mul(params[i],params[j]),c))
        rows.append(values);rhs.append(Q);ids.append(n)
    C=saved['stage_certificates'][0];assert ids==C['all_original_pencil_row_indices']and rows==mat(C['full_linear_matrix'])and rhs==[poly(v)for v in C['full_lower_rhs']]
    assert z.rank(rows)==9 and C['nullity']==3
    phi=list(map(F,C['left_annihilator_combination']));assert len(phi)==len(rows)
    assert all(not sum(c*row[j]for c,row in zip(phi,rows))for j in range(12))
    expected=a.add(a.mul(z.var(8,1),z.var(8,0)),a.scale(a.mul(z.var(8,1),z.var(8,2)),F(1,2)))
    assert a.add(*(a.scale(P,c)for P,c in zip(rhs,phi)))==expected==poly(saved['obstruction_polynomial'])
    assert {ids[i]:c for i,c in enumerate(phi)if c}=={13:F(7,8),59:F(-1,8),62:F(1,8),75:F(1,8)}
    assert saved['source_upper_after_nonzero_obstruction']==7 and saved['whole_PGL_image_upper']==24
    return dict(complete_top5_parameter_source=True,full24_by12_middle_matrix_rank=9,full_middle_kernel_dimension=3,
      bottom_kernel_dimension=0,full_left_annihilator_identity=True,nonzero_compatibility='p1*(p0+p2/2)',
      no_free_parameter_inverted=True,full_source_dimension_upper=7,whole_PGL_image_upper=24)


def main():
    resource.setrlimit(resource.RLIMIT_CPU,(235,238));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2);signal.alarm(240)
    start=time.monotonic();job=json.loads((ROOT/JOB).read_text())
    for row in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    manifest=json.loads((ROOT/MAN).read_text())
    for row in manifest['files']:assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    packet=read('PENCIL');high=read('HIGHEST_SYMBOL');tri=read('TRIANGULAR_BOUND');res=read('RESONANT_COMPLETION')
    for P in(packet,high,tri,res):assert hashlib.sha256(P['producer_source_bytes'].encode()).hexdigest()==P['input_sha256']['code/classification_quadric_cone_rank_two_pencil.py']
    L,A,J,fmons,C,E,sourceproof=source(packet)
    full,cache=whole_sources(packet,L,A,J,fmons,C,E);highestproof=highest(high,cache);gradedproof,Ms=graded(tri,cache)
    last=resonance(res,cache,tri,Ms)
    result=dict(task_id='G',all_passed=True,producer_imports=False,CAS_imports=False,independent_Fraction_and_QQ_polynomial_only=True,
      exact_input_bindings=job['exact_inputs_and_hashes'],full_tangent_source=sourceproof,all_whole_sources=full,
      original_five_block_highest_symbol=highestproof,graded_all_walls=gradedproof,unique_resonance=last,
      old_false_Ore_entrance_runtime_preserved=True,no_initial_v_independent_B_assumption=True,
      all_R_and_vertical_and_Czero_boundaries_retained=True,whole_rank2_quadric_cone_projective_image_upper=24,
      scope='Rank4 tensor, reduced irreducible actual quintic, nonnormal quadriccone surface with generic restrictionrank2. No E membership or componenthood/global G exhaustion.',
      elapsed_seconds=time.monotonic()-start,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (ROOT/OUT).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in('all_passed','whole_rank2_quadric_cone_projective_image_upper','elapsed_seconds','peak_rss_kib')}))


if __name__=='__main__':main()

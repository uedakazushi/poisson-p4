#!/usr/bin/env python3
"""Independent whole seventeen-pencil concurrent M2 double-root review.

Only frozen earlier independent Fraction/sparse helpers are imported. No
producer or CAS is imported. Algebraic walls are checked over the full quotient
algebra with explicit unit minors, so all complex roots are covered without
requiring a factorization or irreducibility oracle.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json,resource,signal,time
import classification_plane_m2_double_tangent_review as D
from classification_unimodular_logarithmic_review import RF,ua,um,ud
H=D.H
ROOT=Path(__file__).resolve().parents[1]
INPUT='reports/CLASSIFICATION_PLANE_M2_CONCURRENT_DOUBLE_20261108.json'
JOB='research/classification/jobs/plane_m2_concurrent_double_separate_review.json'
OUT='reports/CLASSIFICATION_PLANE_M2_CONCURRENT_DOUBLE_SEPARATE_REVIEW_20261108.json'
A=RF((F(0),F(1)))

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def load(p):return json.loads((ROOT/p).read_text())
def trim(p):
    p=list(p)
    while p and not p[-1]:p.pop()
    return tuple(p)

class Quotient:
    """QQ[a]/g with independent coefficient convolution and unit inversion."""
    modulus=()
    def __init__(self,v=0):
        if isinstance(v,Quotient):self.p=v.p;return
        p=list(map(F,v))if isinstance(v,(tuple,list))else[F(v)]
        n=len(self.modulus)-1
        for j in range(len(p)-1,n-1,-1):
            c=p[j]
            if c:
                for k in range(n):p[j-n+k]-=c*self.modulus[k]
        self.p=trim(p[:n])
    def __add__(self,v):
        v=Quotient(v);return Quotient(tuple((self.p[j]if j<len(self.p)else F())+(v.p[j]if j<len(v.p)else F())for j in range(max(len(self.p),len(v.p)))))
    __radd__=__add__
    def __neg__(self):return Quotient(tuple(-c for c in self.p))
    def __sub__(self,v):return self+-Quotient(v)
    def __rsub__(self,v):return Quotient(v)+-self
    def __mul__(self,v):
        v=Quotient(v);p=[F()]*max(0,len(self.p)+len(v.p)-1)
        for i,c in enumerate(self.p):
            for j,d in enumerate(v.p):p[i+j]+=c*d
        return Quotient(p)
    __rmul__=__mul__
    def inv(self):
        # Solve the multiplication matrix for 1. This tests that the element
        # is a unit even if a supplied modulus were reducible.
        n=len(self.modulus)-1
        cols=[(self*Quotient((F(),)*j+(F(1),))).p for j in range(n)]
        M=[[(col[i]if i<len(col)else F())for col in cols]+[F(i==0)]for i in range(n)]
        for j in range(n):
            k=next((k for k in range(j,n)if M[k][j]),None)
            assert k is not None,('nonunit at an algebraic root',self.p,self.modulus)
            M[j],M[k]=M[k],M[j];d=M[j][j];M[j]=[v/d for v in M[j]]
            for k in range(n):
                if k!=j and M[k][j]:
                    c=M[k][j];M[k]=[v-c*w for v,w in zip(M[k],M[j])]
        result=Quotient(tuple(row[-1]for row in M));assert self*result==1;return result
    def __truediv__(self,v):return self*Quotient(v).inv()
    def __rtruediv__(self,v):return Quotient(v)*self.inv()
    def __bool__(self):return bool(self.p)
    def __eq__(self,v):return self.p==Quotient(v).p
    def __repr__(self):return str(tuple(map(str,self.p)))

def compact(v):
    if isinstance(v,RF)and v.q==(F(1),)and len(v.p)<=1:return v.p[0]if v.p else F()
    if isinstance(v,Quotient)and len(v.p)<=1:return v.p[0]if v.p else F()
    return v
def scalar(v):
    if isinstance(v,dict):
        if 'residue'in v:
            assert tuple(map(F,v['modulus']))==Quotient.modulus
            return compact(Quotient(tuple(map(F,v['residue']))))
        return compact(RF(tuple(map(F,v['numerator'])),tuple(map(F,v['denominator']))))
    return F(v)
def matrix(M):return[[scalar(v)for v in row]for row in M]
def actualpoly(v):return{tuple(row['exponent']):scalar(row['coefficient'])for row in v}
def determinant(M):
    B=[list(row)for row in M];out=F(1)
    for j in range(len(B)):
        k=next((k for k in range(j,len(B))if B[k][j]),None)
        if k is None:return F()
        if k!=j:B[k],B[j]=B[j],B[k];out=-out
        d=B[j][j];out=compact(out*d)
        for k in range(j+1,len(B)):
            if B[k][j]:
                c=B[k][j]/d
                for l in range(j+1,len(B)):B[k][l]=compact(B[k][l]-c*B[j][l])
    return out
def certificate(M,C):
    rr,cc=C['rows'],C['columns'];r=C['rank'];n=len(M[0]);assert len(rr)==len(cc)==r
    assert len(set(rr))==len(set(cc))==r
    value=determinant([[M[i][j]for j in cc]for i in rr]);assert value and value==scalar(C['minor_determinant'])
    if isinstance(value,Quotient):assert value*value.inv()==1
    K=matrix(C['complete_kernel']);assert len(K)==C['nullity']==n-r
    free=[j for j in range(n)if j not in cc]
    assert all(len(row)==n for row in K)
    assert [[row[j]for j in free]for row in K]==[[F(i==j)for j in range(len(K))]for i in range(len(K))]
    for row in M:
        nz=[(i,c)for i,c in enumerate(row)if c]
        assert all(not sum((c*k[i]for i,c in nz if k[i]),F())for k in K)
    return K,value

def models(kind,b,c):
    x,y,z,u,v=[H.a(1,i)for i in range(5)]
    f=H.mul(H.mul(x,y),H.add(x,H.scale(y,-1)))
    Pf=H.skew({(0,2):H.scale(H.der(f,1),-1),(1,2):H.der(f,0)})
    Y=[{},{},H.add(H.scale(x,A),H.scale(y,b))]
    if kind=='h21':
        X=[{},{},H.add(H.scale(x,A),H.scale(y,b+3))]
        Z=[H.scale(x,-1),H.scale(y,2),H.add(H.scale(z,-1),H.scale(x,c))];k=b+4
        target=H.scale(H.a(1,0,0,1,3,3),3*(A+b+4))
        excluded={-b-4}
    else:
        E=-(b+8)/(3*(A+b+4));k=E*(b+8)
        X=[{},{},H.add(H.scale(x,A-3),H.scale(y,b+6))]
        Z=[{},H.add(x,H.scale(y,3*E)),H.scale(z,-3*E)]
        target=H.scale(H.a(1,0,0,0,3,3),A-4)
        excluded={F(4),-b-4}
    P=H.plus(Pf,H.skew({(i,3):H.mul(u,X[i])for i in range(3)}),H.skew({(i,4):H.add(H.mul(u,Z[i]),H.mul(v,Y[i]))for i in range(3)}),H.skew({(3,4):H.scale(H.mul(u,u),k)}))
    return P,target,excluded

def hom(v,point):
    if not isinstance(v,RF):return F(v)
    if len(v.p)<=1 and len(v.q)==1:return compact(v)
    def value(poly):
        out=F()
        for c in reversed(poly):out=compact(out*point+c)
        return out
    return compact(value(v.p)/value(v.q))
def specialize_poly(P,point):return{e:hom(c,point)for e,c in P.items()if hom(c,point)}
def specialize_tensor(P,point):return H.skew({ij:specialize_poly(P.get(ij,{}),point)for ij in H.PAIRS})
def monic(g):g=tuple(map(F,g));assert g and g[-1];return tuple(c/g[-1]for c in g)
def supported(poly,factors):
    p=tuple(poly);found=set()
    for g in factors:
        while len(p)>=len(g):
            q,rem=ud(p,g)
            if rem:break
            p=q;found.add(g)
    assert len(p)<=1,('uncovered complex parameter roots',p)
    return found

def check(row,P,target,operators,B1,B2):
    assert H.vector(P)==list(map(scalar,row['full150_P0']))
    assert not any(H.jac(P).values())and H.actual(P)==target==actualpoly(row['actualG2'])and target
    J1,J2,G=operators
    assert J1==matrix(row['full350by30_J1'])and J2==matrix(row['full350by9_J2'])and G==matrix(row['full150_Euler_Levi_scale_matrix'])
    K1,d1=certificate(J1,row['J1_certificate']);K2,d2=certificate(J2,row['J2_certificate']);KG,dg=certificate(G,row['normalizer_certificate'])
    bad=[[sum((Q.get((i,3),{}).get((0,0,0,0,2),F())*k[j]for j,Q in enumerate(B1)),F())for k in K1]for i in range(3)]
    assert bad==matrix(row['all_bad_urow_maps'])and not any(bad[0])and not any(bad[1])
    ufactor=not any(bad[2]);assert ufactor==row['actual_u_factor_for_every_solution']
    assert all(e[3]for i in range(5)if i!=3 for e in P.get((i,3),{}))
    assert all(not Q.get((i,3),{})for Q in B2 for i in range(3))
    s=37-row['normalizer_certificate']['rank'];fiber=s+6
    upper=24-fiber+row['source_parameter_dimension']+len(K1)+len(K2)
    assert s==row['graded_Levi_projective_stabilizer_dimension']and fiber==row['actual_filtered_group_fiber_dimension']and upper==row['whole_projective_image_dimension_upper']
    assert row['G24_closed']==(ufactor or upper<=24)
    assert len(row['all_full39_actualD_linear_terms'])==39
    for Q,saved in zip(B1+B2,row['all_full39_actualD_linear_terms']):
        assert H.add(H.actual(H.plus(P,Q)),H.scale(target,-1))==actualpoly(saved)
    result=dict(type=row['type'],field=row['field'],K1=len(K1),K2=len(K2),filtered_group_fiber=fiber,parameter_dimension=row['source_parameter_dimension'],image_upper=upper,entire_actual_u_factor=ufactor,all_full350_150_maps=True,all_rank_kernel_unit_certificates=True)
    return result,(d1,d2,dg),K1+K2+KG

def main():
    begun=time.monotonic();signal.alarm(420);resource.setrlimit(resource.RLIMIT_AS,(3072*1024**2,)*2)
    job=load(JOB);assert all(sha(r['path'])==r['sha256']for r in job['exact_inputs_and_hashes'])
    raw=load(INPUT);assert hashlib.sha256(raw['producer_source_bytes'].encode()).hexdigest()==raw['producer_sha256']=='448c081c5110dc7f99e53f29631c4d655d16b29cd5597a35a66e490b601e70b8'
    B1,B2=H.sources();assert[H.vector(P)for P in B1]==matrix(raw['full30_weight1_basis'])and[H.vector(P)for P in B2]==matrix(raw['full9_weight2_basis'])
    assert raw['full350_row_labels']==[dict(quad=list(q),monomial=list(e))for q,e in H.ROWS]
    # Previously independently reconstructed full 39-polarization table is
    # reused by exact equality; no old suite is rerun.
    prior=load('reports/CLASSIFICATION_PLANE_M2_DOUBLE_TANGENT_20261107.json')
    assert raw['universal_full39_quadratic_Jacobi_and_actualD']==prior['universal_full39_quadratic_Jacobi_and_actualD']
    records=[];all_fields=0;nonrational=0;residual=[]
    expected=[('h21',F(b),F(c))for b in(-2,-5,-8,F(-7,2),-4)for c in(0,1)]+[('h3',F(b),None)for b in(-14,-11,F(-19,2),-8,F(-13,2),-5,-2)]
    assert len(raw['cases'])==len(expected)==17
    for row,(kind,b,c)in zip(raw['cases'],expected):
        assert row['kind']==kind and F(row['b'])==b and(None if row['c']is None else F(row['c']))==c
        P,target,excluded=models(kind,b,c)
        assert set(map(F,row['source_domain_excluded_points']))==excluded
        ops=(D.op(P,B1),D.op(P,B2),D.group(P))
        result,dets,K=check(row,P,target,ops,B1,B2)
        assert result['entire_actual_u_factor']
        euler=[[F(v)for v in line[:25]]for line in ops[2]]
        assert len(H.rref(euler)[1])==24
        factors={monic(g)for g in row['all_generic_minor_and_kernel_factors']}
        used=set()
        for value in list(dets)+[v for k in K for v in k]:
            z=RF(value)
            used|=supported(z.p,factors);used|=supported(z.q,factors)
        assert used==factors
        domain_factors=factors|{(-a,F(1))for a in excluded}
        for O in ops:
            for line in O:
                for value in line:supported(RF(value).q,domain_factors)
        walls=[];wall_polys=set()
        for wall in row['all_allowed_parameter_walls']:
            g=monic(wall['minimal_polynomial']);assert g in factors and g not in wall_polys;wall_polys.add(g)
            if len(g)==2:
                point=-g[0];assert point==F(wall['a'])and point not in excluded
            else:
                Quotient.modulus=g;point=Quotient((0,1));nonrational+=1
                assert wall['all_complex_embeddings_covered']
            W=specialize_tensor(P,point);WT=specialize_poly(target,point)
            os=tuple([[hom(v,point)for v in line]for line in O]for O in ops)
            checked,_,_=check(wall,W,WT,os,B1,B2);checked['modulus']=list(map(str,g));checked['all_complex_roots_verified_by_unit_minors']=True
            walls.append(checked);all_fields+=1
            if not checked['entire_actual_u_factor']:residual.append(checked)
        skip=set()
        for omitted in row['excluded_normalizer_or_m3_points']:
            g=monic(omitted['minimal_polynomial']);assert len(g)==2 and -g[0]in excluded and F(omitted['a'])==-g[0];skip.add(g)
        assert wall_polys|skip==factors and not wall_polys&skip
        result.update(all_walls=walls,all_parameter_factors=list(map(lambda g:list(map(str,g)),sorted(factors))),whole_complex_parameter_line_covered=True)
        records.append(result)
        print(json.dumps(dict(pencil=row['type'],walls=len(walls),complete=True)),flush=True)
    assert(all_fields,nonrational,len(residual))==(84,5,3)
    assert all(r['image_upper']==23 and r['parameter_dimension']==0 for r in residual)
    assert {(r['type'].split('_wall_')[0],r['K1'],r['K2'],r['filtered_group_fiber'])for r in residual}=={('h21_b-2_c0',9,0,10),('h3_b-2',8,0,9)}
    out=dict(task_id='G',all_independent_arithmetic_passed=True,no_producer_or_CAS_imports=True,generic_pencils=17,allowed_wall_quotient_algebras=all_fields,nonrational_wall_quotient_algebras=nonrational,all_complex_roots_covered_by_unit_minors=True,cases=records,irreducible_residual_sources=residual,whole_concurrent_m2_double_image_upper=23,individual_E_membership_or_existence_not_claimed=True,universal39_quadratic_table_reused_by_exact_equality=True,inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),seconds=time.monotonic()-begun,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    dest=ROOT/OUT;assert not dest.exists();dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in('all_independent_arithmetic_passed','seconds','maximum_RSS_KiB')}))
def triple_double_eight_main():
    """Eight original full operators, corrected V2, and whole-P2 Ore proof."""
    started=time.monotonic();signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    jobp='research/classification/jobs/plane_triple_line_m2_double_eight_separate_review.json'
    job=load(jobp);assert all(sha(r['path'])==r['sha256']for r in job['exact_inputs_and_hashes'])
    old=load('reports/CLASSIFICATION_PLANE_TRIPLE_LINE_M2_DOUBLE_EIGHT_FULLHIGHER_20261108.json')
    new=load('reports/CLASSIFICATION_PLANE_TRIPLE_LINE_M2_DOUBLE_EIGHT_FULLHIGHER_V2_20261108.json')
    assert hashlib.sha256(old['producer_source'].encode()).hexdigest()==old['producer_sha256']=='67c58613e9180447c76692e62e11323fcd52591b89d515f3535e56a4b6e3322e'
    assert hashlib.sha256(new['producer_source'].encode()).hexdigest()==new['producer_sha256']
    import ast
    names=('a','b','c','u','v')
    def parse(text):
        def go(n):
            if isinstance(n,ast.Constant):return H.a(F(n.value))
            if isinstance(n,ast.Name):return H.a(1,names.index(n.id))
            if isinstance(n,ast.UnaryOp):return H.scale(go(n.operand),-1 if isinstance(n.op,ast.USub)else 1)
            if isinstance(n,ast.BinOp):
                if isinstance(n.op,ast.Add):return H.add(go(n.left),go(n.right))
                if isinstance(n.op,ast.Sub):return H.add(go(n.left),H.scale(go(n.right),-1))
                if isinstance(n.op,ast.Mult):return H.mul(go(n.left),go(n.right))
                if isinstance(n.op,ast.Div):return H.scale(go(n.left),1/go(n.right)[(0,)*5])
                if isinstance(n.op,ast.Pow):
                    z=H.a(1)
                    for _ in range(n.right.value):z=H.mul(z,go(n.left))
                    return z
            raise ValueError(ast.dump(n))
        return go(ast.parse(text,mode='eval').body)
    def combo(bs,ks):return H.skew({ij:H.add(*(H.scale(P.get(ij,{}),k)for P,k in zip(bs,ks)if k))for ij in H.PAIRS})
    def qremove(P):
        B=H.der(P.get((1,2),{}),2)
        assert not B or all(e[0]+e[1]+e[2]==0 and e[3]+e[4]==1 for e in B)
        E=[H.a(1,i)for i in range(5)];T=[{},B,{},{},{}]
        U=H.plus(P,H.wedge(E,T));assert all(not e[2]for ij in H.PAIRS for e in U.get(ij,{}))
        return H.skew({ij:U.get(ij,{})for ij in H.PAIRS if 2 not in ij}),B
    source=load('research/classification/inputs/plane_triple_line_m2_constant_double_eight_first_symbols_20261108.json')
    B1,B2=H.sources();a,b,c,u,v=[H.a(1,i)for i in range(5)];records=[];ore=[]
    assert len(old['cases'])==len(new['cases'])==len(source['profiles'])==8
    for row,v2,profile in zip(old['cases'],new['cases'],source['profiles']):
        assert row['name']==v2['name']==profile['name']
        for key in row:
            if key not in{'actual_P0_zero','actual_D_has_proper_a_factor','G24_closed'}:assert row[key]==v2[key],(row['name'],key)
        assert row['actual_P0_zero']is True and v2['actual_P0_zero']is False
        P=H.skew({tuple(names.index(s)for s in ij.split(',')):parse(text)for ij,text in profile['raw_brackets'].items()})
        assert H.vector(P)==list(map(F,row['full150_P0']))and not any(H.jac(P).values())
        actual=H.actual(P);assert actual==parse(row['actual_G2'])==parse(profile['G2'])and actual
        M1,M2,G=D.op(P,B1),D.op(P,B2),D.group(P)
        assert M1==matrix(row['complete350x30'])and M2==matrix(row['complete350x9'])and G==matrix(row['full150x39_Euler_Levi_scale'])
        K1,_=certificate(M1,row['weight1_unconstrained_certificate']);K2,_=certificate(M2,row['weight2_certificate']);KG,_=certificate(G,row['filtered_group_certificate'])
        assert row['constrained_weight1_certificate']==row['weight1_unconstrained_certificate']
        higher1=[combo(B1,k)for k in K1];hom2=[combo(B2,k)for k in K2]
        def factor(families,at):return all(e[at]>0 for Q in families for i in range(5)if i!=at for e in Q.get((at,i),{}))
        uf=factor([P]+higher1+B2,3);af=factor([P]+higher1+B2,0)
        assert uf==row['actual_D_has_proper_u_factor']==v2['actual_D_has_proper_u_factor']
        assert factor([P]+higher1+hom2,0)==row['actual_D_has_proper_a_factor']
        assert af==v2['actual_D_has_proper_a_factor']
        B=row['filtered_group_certificate']['rank']-24;bound=B+len(K1)+len(K2)+5
        assert B==row['first_symbol_dimension']and bound==row['whole_projective_image_upper']
        assert v2['G24_closed']==(uf or af or bound<=24)
        records.append(dict(name=row['name'],K1=len(K1),K2=len(K2),first_symbol_dimension=B,coarse_image_upper=bound,whole_u_factor=uf,whole_a_factor_with_FULL9=af,old_metadata_rejected=True,all_original_operators_kernels_minors=True))
        if row['name'].startswith('NIL'):
            Q0,_=qremove(P);assert H.vector(Q0)==H.vector(H.skew({(1,4):H.mul(a,u)}))
            epsilon=int(row['name'][-1]);checks=[]
            for j,Q in enumerate(higher1):
                extracted,gotB=qremove(Q);t=[F(i==j)for i in range(len(K1))]
                if epsilon==0:
                    L=(t[2]+t[4]+3*t[6]+t[10])/2;V=(t[5]-t[7])/2
                    predicted=H.skew({(0,1):H.mul(a,H.add(H.scale(u,L),H.scale(v,V))),(1,3):H.add(H.scale(H.mul(u,u),-t[2]),H.scale(H.mul(u,v),3*t[7]-t[5])),(1,4):H.add(H.scale(H.mul(u,u),t[8]),H.scale(H.mul(u,v),t[10]-t[2]),H.scale(H.mul(v,v),-t[5]))})
                    targetB=H.add(H.scale(u,t[2]),H.scale(v,t[5]))
                else:
                    aa=t[3]/3+t[4]-2*t[6]-t[9]/3;tt=t[6];jj=t[9];rr=t[7];m=t[5];L=2*aa+3*tt+jj
                    predicted=H.skew({(0,1):H.mul(a,H.add(H.scale(u,L),H.scale(v,10*m/3))),(0,4):H.scale(H.mul(u,u),tt),(1,3):H.add(H.scale(H.mul(u,u),-aa),H.scale(H.mul(u,v),-14*m/3)),(1,4):H.add(H.scale(H.mul(u,u),rr),H.scale(H.mul(u,v),jj-aa),H.scale(H.mul(v,v),-23*m/3))})
                    targetB=H.add(H.scale(u,aa),H.scale(v,23*m/3))
                assert H.vector(extracted)==H.vector(predicted)and gotB==targetB
                if epsilon:
                    bad=Q.get((2,3),{}).get((0,0,0,0,2),F());assert bad==t[5]
                checks.append(True)
            projected2=[qremove(Q)[0]for Q in B2]
            assert all(H.vector(Q)==H.vector(H.skew({(0,1):Q.get((0,1),{})}))for Q in projected2)
            assert len([Q for Q in projected2 if any(Q.values())])==3
            ore.append(dict(name=row['name'],all_K1_basis_Q_formulas_and_Euler_gauge=True,all_FULL9_P2_retained=True,arbitrary_F2_in_Qab=True,original_bad_m_coefficient_verified=bool(epsilon),constant_factor_Q_when_Tzero=True))
    assert sum(x['whole_u_factor']for x in records)==6 and not any(x['whole_a_factor_with_FULL9']for x in records)
    assert records[-2]['coarse_image_upper']==records[-1]['coarse_image_upper']==28
    out=dict(task_id='G',all_independent_arithmetic_passed=True,no_producer_or_CAS_imports=True,all8_full350_30_9_150_group_certificates=records,V1_V2_all_math_arrays_identical=True,historical_defects_not_adopted=['V1 actual_P0_zero=true contradicts verified nonzero G2','V1 NIL0 a-factor used only homogeneous K2, missing affine particular; withdrawn'],whole_Ore_Q_extraction=ore,independent_pureproof_scope='c coefficient gives strict Q and module; constant factor gives F independent of b,c; NIL1 strict abu/abv give Tm=0 and T(A+3T)=0 with arbitrary F2. Full9 retained. See separate proof note.',accepted_result='Six semisimple sources empty; both nil whole solutions empty or in accepted essential-three eight components; no eight-source coarse linear bound called a proof.',inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    dest=ROOT/'reports/CLASSIFICATION_PLANE_TRIPLE_LINE_M2_DOUBLE_EIGHT_SEPARATE_REVIEW_20261108.json';assert not dest.exists();dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in('all_independent_arithmetic_passed','seconds','maximum_RSS_KiB')}))

if __name__=='__main__':
    import sys
    if '--triple-line-double-eight'in sys.argv:triple_double_eight_main()
    else:main()

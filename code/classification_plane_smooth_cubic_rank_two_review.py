#!/usr/bin/env python3
"""Independent universal plane/smooth-cubic weight reduction.

Fraction sparse arithmetic only; no producer, CAS or repository imports.
No cubic or nonzero binary-quadratic coefficient is specialized.
"""
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(it.combinations(range(5),2))
TRIPLES=list(it.combinations(range(5),3))
QUADS=list(it.combinations(range(5),4))
BP=list(it.combinations(range(3),2))

def add(*args):
    d={}
    for p in args:
        for k,v in p.items():d[k]=d.get(k,F(0))+v
    return {k:v for k,v in d.items()if v}

def scale(p,a):return {k:a*v for k,v in p.items()if a*v}

def mul(*args):
    if not args:return {}
    p=args[0]
    for q in args[1:]:
        d={}
        for e,a in p.items():
            for f,b in q.items():
                k=tuple(x+y for x,y in zip(e,f));d[k]=d.get(k,F(0))+a*b
        p={k:v for k,v in d.items()if v}
    return p

def der(p,i):
    d={}
    for e,a in p.items():
        if e[i]:
            f=list(e);f[i]-=1;d[tuple(f)]=a*e[i]
    return d

def mons(n,d):
    out=[]
    for ss in it.combinations_with_replacement(range(n),d):
        e=[0]*n
        for i in ss:e[i]+=1
        out.append(tuple(e))
    return sorted(out,reverse=True)

def rref(rows):
    b={};chosen=[]
    for i,row in enumerate(rows):
        v={k:F(a)for k,a in row.items()if a}
        for p in sorted(b):
            if p in v:v=add(v,scale(b[p],-v[p]))
        if not v:continue
        p=min(v);v=scale(v,1/v[p])
        for q in list(b):
            if p in b[q]:b[q]=add(b[q],scale(v,-b[q][p]))
        b[p]=v;chosen.append(i)
    return b,chosen

class Ring:
    def __init__(self,names):self.names=names;self.n=len(names);self.zero=(0,)*self.n
    def one(self,a=1):return {self.zero:F(a)}if a else {}
    def v(self,s):
        e=list(self.zero);e[self.names.index(s)]=1;return {tuple(e):F(1)}
    def power(self,p,k):
        q=self.one()
        for _ in range(k):q=mul(q,p)
        return q

def antisym(P):
    P={ij:f for ij,f in P.items()if f}
    for i,j in list(P):
        if i<j:P[j,i]=scale(P[i,j],-1)
    return P

def tensor_add(*values):return antisym({ij:add(*(p.get(ij,{})for p in values))for ij in PAIRS})

def tensor_scale(P,a):return antisym({ij:scale(P.get(ij,{}),a)for ij in PAIRS})

def main():
    resource.setrlimit(resource.RLIMIT_CPU,(169,170))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    start=time.monotonic()
    inputs=['research/classification/PLANE_SMOOTH_CUBIC_RANK_TWO_EXCLUSION.md',
            'research/classification/PLANE_SMOOTH_CUBIC_CONORMAL_CONNECTION_REVIEW.md']
    names=[f'f{i}'for i in range(10)]+['a']+[f'l{i}'for i in range(9)]+[f'm{i}'for i in range(9)]
    names += [f'{v}{i}'for v in ['b','c','d','e']for i in range(3)]
    names += [f'q{i}'for i in range(9)]+[f'{v}{i}'for v in ['Y','Z']for i in range(3)]+['x0','x1','x2','z','w']
    R=Ring(names);cx=list(range(R.n-5,R.n));x=[R.v(s)for s in names[-5:]];z,w=x[3:];a=R.v('a')
    def xp(i):return x[i]
    def diff(p,i):return der(p,cx[i])
    def mono(e):return {tuple((0,)*(R.n-5)+tuple(e)):F(1)}
    f=add(*(mul(R.v(f'f{k}'),mono(tuple(e)+(0,0)))for k,e in enumerate(mons(3,3))))
    Pf=antisym({(0,1):diff(f,2),(0,2):scale(diff(f,1),-1),(1,2):diff(f,0)})
    def vec(s):return [R.v(f'{s}{i}')for i in range(3)]
    B,C,D,E,Y,Z=map(vec,['b','c','d','e','Y','Z'])
    def wedge(V,W):return antisym({(i,j):add(mul(V[i],W[j]),scale(mul(V[j],W[i]),-1))for i,j in BP})
    def U(V):return antisym({ij:scale(add(*(mul(V[k],diff(Pf.get(ij,{}),k))for k in range(3))),-1)for ij in BP})
    L=antisym({ij:add(*(mul(R.v(f'l{3*k+i}'),x[i])for i in range(3)))for k,ij in enumerate(BP)})
    M=antisym({ij:add(*(mul(R.v(f'm{3*k+i}'),x[i])for i in range(3)))for k,ij in enumerate(BP)})
    def first(L,M,B,C,D,E):
        out={ij:add(mul(z,L.get(ij,{})),mul(w,M.get(ij,{})))for ij in BP}
        for i in range(3):
            out[i,3]=add(mul(z,w,B[i]),mul(w,w,C[i]))
            out[i,4]=add(mul(z,z,D[i]),mul(z,w,E[i]))
        return antisym(out)
    P1=first(L,M,B,C,D,E)
    P2=antisym({ij:add(*(mul(R.v(f'q{3*k+j}'),q)for j,q in enumerate([mul(z,z),mul(z,w),mul(w,w)])))for k,ij in enumerate(BP)})
    def strict(P,Q):
        return {tri:add(*(add(mul(P.get((i,l),{}),diff(Q.get((j,k),{}),l)),
                              mul(Q.get((i,l),{}),diff(P.get((j,k),{}),l)))
                         for i,j,k in [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]for l in range(5)))for tri in TRIPLES}
    def euler(S):return {quad:add(*(scale(mul(x[quad[i]],S[quad[:i]+quad[i+1:]]),(-1)**i)for i in range(4)))for quad in QUADS}
    def bracket(P,Q):return euler(strict(P,Q))
    def ambient(P):
        return {(ij,e):p for ij in PAIRS for e,p in P.get(ij,{}).items()}
    def normal_part(p,dz,dw):
        ans={}
        for ee,v in p.items():
            if ee[cx[3]]==dz and ee[cx[4]]==dw:
                ee=list(ee);ee[cx[3]]=ee[cx[4]]=0;ans[tuple(ee)]=v
        return ans
    def scalar3(S):return S.get((0,1,2),{})
    def base_euler(T):return add(mul(x[0],T.get((1,2),{})),scale(mul(x[1],T.get((0,2),{})),-1),mul(x[2],T.get((0,1),{})))
    # Complete homogeneous coefficient source and exact Euler quotient.
    cols=[(ij,e)for ij in PAIRS for e in mons(5,2)]
    def basis_tensor(j):
        ij,e=cols[j];return antisym({ij:mono(e)})
    basis=[basis_tensor(j)for j in range(150)]
    def flat(P):return {j:P.get(ij,{}).get(tuple((0,)*(R.n-5)+e),F(0))for j,(ij,e)in enumerate(cols)if P.get(ij,{}).get(tuple((0,)*(R.n-5)+e),F(0))}
    weight=[sum(e[3:])-sum(i>=3 for i in ij)for ij,e in cols]
    raw={d:[j for j,v in enumerate(weight)if v==d]for d in range(-2,3)}
    assert [len(raw[d])for d in range(-2,3)]==[6,42,57,36,9]
    V=[];vw=[]
    for i,j in it.product(range(5),repeat=2):
        vec5=[{}for _ in range(5)];vec5[i]=x[j]
        T=antisym({ij:add(mul(x[ij[0]],vec5[ij[1]]),scale(mul(x[ij[1]],vec5[ij[0]]),-1))for ij in PAIRS})
        V.append(flat(T));vw.append((j>=3)-(i>=3))
    assert len(rref(V)[0])==24
    def chart(P,i,j):
        value=add(P.get((i,j),{}),scale(mul(x[j],P.get((i,0),{})),-1),scale(mul(x[i],P.get((0,j),{})),-1))
        ans={}
        for ee,c in value.items():
            e=list(ee);e[cx[0]]=0;e=tuple(e);ans[e]=ans.get(e,F(0))+c
        return {e:c for e,c in ans.items()if c}
    def constraints(source,zero_connection=False):
        rows={}
        for column,P in enumerate(source):
            tests=[]
            if zero_connection:
                tests.append(((1,2),0));tests += [((i,j),1)for i in [1,2]for j in [3,4]]
            else:tests += [((i,j),0)for i in [1,2]for j in [3,4]]+[((3,4),0),((3,4),1)]
            for ij,d in tests:
                for ee,c in chart(P,*ij).items():
                    if ee[cx[3]]+ee[cx[4]]==d:rows.setdefault((ij,ee),{})[column]=c
        return list(rows.values())
    invrows=constraints(basis);invred,invchosen=rref(invrows)
    assert len(invred)==42
    allowed=[flat(basis[j])for d in [0,1,2]for j in raw[d]]+[v for v,d in zip(V,vw)if d==-1]
    assert len(rref(allowed)[0])==108
    assert all(not sum(c*v.get(j,F(0))for j,c in row.items())for row in invrows for v in allowed)
    c0=constraints([basis[j]for j in raw[0]],True);c0red,c0chosen=rref(c0)
    assert len(c0red)==42
    zero0=[v for v,d in zip(V,vw)if d==0]
    Kbasis=[flat(antisym({(3,4):q}))for q in [mul(z,z),mul(z,w),mul(w,w)]]
    assert len(rref(zero0)[0])==12 and len(rref(zero0+Kbasis)[0])==15
    assert all(not sum(c*v.get(raw[0][j],F(0))for j,c in row.items())for row in c0 for v in zero0+Kbasis)
    p1vars=[f'l{i}'for i in range(9)]+[f'm{i}'for i in range(9)]+[f'{s}{i}'for s in ['b','c','d','e']for i in range(3)]
    def parameter_vectors(P,parameters):
        result=[]
        for v in parameters:
            k=names.index(v);T={}
            for ij in PAIRS:
                value={}
                for ee,c in P.get(ij,{}).items():
                    if ee[k]:
                        e=list(ee);assert e[k]==1;e[k]=0;value[tuple(e)]=c
                T[ij]=value
            result.append(flat(T))
        return result
    p1columns=parameter_vectors(P1,p1vars);g1=[v for v,d in zip(V,vw)if d==1]
    assert len(p1columns)==30 and len(rref(p1columns)[0])==30
    assert len(rref(g1)[0])==6 and len(rref(p1columns+g1)[0])==36
    assert set().union(*(set(v)for v in p1columns+g1))==set(raw[1])
    assert len(parameter_vectors(P2,[f'q{i}'for i in range(9)]))==len(raw[2])==9
    # Universal exact identities in all ten coefficients of f and every a.
    modes={};A00=wedge(x[:3],E);A0B=wedge(x[:3],B);A0C=wedge(x[:3],C);A0D=wedge(x[:3],D)
    for mode,K in [('zero',{}),('distinct',mul(a,z,w)),('square',mul(a,z,z))]:
        P0=tensor_add(Pf,antisym({(3,4):K}))
        assert all(not v for v in strict(P0,P0).values())
        S=strict(P0,P1);J=euler(S);assert len(J)==5
        if mode=='zero':expected=[U(D),U(E),tensor_scale(U(B),-1),tensor_scale(U(C),-1)]
        elif mode=='distinct':expected=[tensor_add(U(D),{ij:mul(a,v)for ij,v in A0D.items()}),
                                      tensor_add(U(E),{ij:scale(mul(a,v),-1)for ij,v in L.items()}),
                                      tensor_add(tensor_scale(U(B),-1),{ij:scale(mul(a,v),-1)for ij,v in M.items()}),
                                      tensor_add(tensor_scale(U(C),-1),{ij:mul(a,v)for ij,v in A0C.items()})]
        else:expected=[tensor_add(U(D),{ij:scale(mul(a,v),-1)for ij,v in tensor_add(L,A00).items()}),
                       tensor_add(U(E),{ij:scale(mul(a,v),-1)for ij,v in M.items()},{ij:mul(a,v)for ij,v in A0B.items()}),
                       tensor_add(tensor_scale(U(B),-1),{ij:scale(mul(a,v),2)for ij,v in A0C.items()}),
                       tensor_scale(U(C),-1)]
        for k,(dz,dw)in enumerate([(3,0),(2,1),(1,2),(0,3)]):
            assert all(normal_part(J[(i,j,3,4)],dz,dw)==expected[k].get((i,j),{})for i,j in BP)
        if mode=='distinct':
            assert normal_part(J[(0,1,2,4)],2,0)==base_euler(U(D))
            assert normal_part(J[(0,1,2,3)],0,2)==base_euler(U(C))
        if mode!='zero':
            if mode=='distinct':closed=first(U(Y),U(Z),[scale(mul(a,v),-1)for v in Z],[{}]*3,[{}]*3,[mul(a,v)for v in Y])
            else:closed=first(tensor_add(U(Y),{ij:scale(mul(a,v),-1)for ij,v in wedge(x[:3],Z).items()}),U(Z),[{}]*3,[{}]*3,[mul(a,v)for v in Y],[mul(a,v)for v in Z])
            assert all(not v for v in bracket(P0,closed).values())
            # Actual infinitesimal coordinate shear, followed only by Euler gauge.
            X=[add(mul(z,Y[i]),mul(w,Z[i]))for i in range(3)]+[{},{}]
            orbit=antisym({(i,j):add(*(add(mul(P0.get((l,j),{}),diff(X[i],l)),
                                        mul(P0.get((i,l),{}),diff(X[j],l)),
                                        scale(mul(X[l],diff(P0.get((i,j),{}),l)),-1))for l in range(5)))for i,j in PAIRS})
            if mode=='square':
                VG=[mul(a,z,v)for v in Z]+[{},{}]
                gauge=antisym({(i,j):add(mul(x[i],VG[j]),scale(mul(x[j],VG[i]),-1))for i,j in PAIRS})
                orbit=tensor_add(orbit,tensor_scale(gauge,-1))
            assert all(orbit.get(ij,{})==closed.get(ij,{})for ij in PAIRS)
            # The shear is an actual global linear map, not merely a tangent.
            # Its inverse substitutes x -> x-zY-wZ; both fix z,w.
            substitutions=[add(x[i],scale(X[i],-1))for i in range(5)]
            def subst(p):
                ans={}
                for ee,c in p.items():
                    pp={tuple(ee[:-5]+(0,)*5):c}
                    for i,n in enumerate(ee[-5:]):pp=mul(pp,R.power(substitutions[i],n))
                    ans=add(ans,pp)
                return ans
            T=[[R.one(i==j)for j in range(5)]for i in range(5)]
            for i in range(3):T[i][3]=Y[i];T[i][4]=Z[i]
            actual=antisym({(i,j):add(*(mul(T[i][aa],T[j][bb],subst(P0.get((aa,bb),{})))
                                              for aa,bb in it.product(range(5),repeat=2)))for i,j in PAIRS})
            raworbit=orbit if mode=='distinct' else tensor_add(orbit,gauge)
            for ij in PAIRS:
                parts={d:{}for d in range(3)}
                for ee,c in actual.get(ij,{}).items():
                    d=ee[cx[3]]+ee[cx[4]]-sum(i>=3 for i in ij)
                    assert d in parts;parts[d][ee]=c
                assert parts[0]==P0.get(ij,{}) and parts[1]==raworbit.get(ij,{})
        J2=bracket(P0,P2)
        assert all(J2[(i,j,3,4)]==scale(mul(K,P2.get((i,j),{})),-2)for i,j in BP)
        pf={quad:add(mul(P0.get((quad[0],quad[1]),{}),P0.get((quad[2],quad[3]),{})),
                     scale(mul(P0.get((quad[0],quad[2]),{}),P0.get((quad[1],quad[3]),{})),-1),
                     mul(P0.get((quad[0],quad[3]),{}),P0.get((quad[1],quad[2]),{})))for quad in QUADS}
        divisor=add(*(scale(mul(x[i],pf[tuple(j for j in range(5)if j!=i)]),(-1)**i)for i in range(5)))
        assert divisor==scale(mul(f,K),3)
        # Literal all350 row packet; parameter polynomials remain universal.
        packet=[]
        for quad in QUADS:
            for ee in mons(5,4):
                packet.append({e[:-5]:c for e,c in J[quad].items()if e[-5:]==ee})
        assert len(packet)==350
        modes[mode]={'all350_universal_linear_Jacobi_reconstructed':True,
                     'normal_cubic_four_bivector_coefficients_verified':True,
                     'full_weight2_forcing_minus2KQ_verified':True,
                     'actual_Euler_Pfaffian_is_3fK':True,
                     'nonzero_full350_rows':sum(bool(p)for p in packet),
                     'full350_sha256':hashlib.sha256(repr(packet).encode()).hexdigest()}
    # Complete K=0 kernel: base linear bivector closed iff its curl is zero.
    # Universal identity uses all 9 coefficients, not a chosen Hessian cubic.
    l=[L.get((1,2),{}),scale(L.get((0,2),{}),-1),L.get((0,1),{})]
    curl=[add(diff(l[2],1),scale(diff(l[1],2),-1)),add(diff(l[0],2),scale(diff(l[2],0),-1)),add(diff(l[1],0),scale(diff(l[0],1),-1))]
    SL=strict(Pf,L)[(0,1,2)]
    expected=add(*(mul(diff(f,i),curl[i])for i in range(3)))
    assert SL==expected or SL==scale(expected,-1)
    curlrows=[]
    for p in curl:
        curlrows.append({i:p.get(tuple(1 if j==names.index(v)else 0 for j in range(R.n)),F(0))for i,v in enumerate([f'l{i}'for i in range(9)])})
    assert len(rref(curlrows)[0])==3
    for Vv in [B,C,D,E,Y,Z]:
        vf=add(*(mul(Vv[i],diff(f,i))for i in range(3)))
        assert base_euler(U(Vv))==scale(vf,-2)
    # Constant-base bivectors commute with the Jacobian bivector universally.
    assert all(not v for v in strict(Pf,P2).values())
    out={'task_id':'E_reduced_nonnormal','all_passed':True,
         'scope':'Every smooth ternary cubic f and every a!=0 in both binary quadratic types; K=0 separately retained. Connection rigidity is an accepted geometric input.',
         'coefficient_arithmetic':'Independent standard-library Fraction sparse polynomials; no producer/CAS/repository arithmetic imports',
         'all10_ternary_cubic_coefficients_kept':True,'all30_weight1_parameters_kept':True,'all9_weight2_parameters_kept':True,
         'raw_weight_dimensions':[6,42,57,36,9],
         'invariant_and_abelian_normal_source':{'ambient150_constraint_rank':42,'complete_kernel':108,'nonnegative_raw':102,'negative_Euler_kernel':6},
         'weight0_fixed_restriction_and_connection':{'raw57_constraint_rank':42,'kernel15':True,'Euler12_plus_binaryK3_complete':True},
         'weight1_quotient':{'raw36':True,'Euler6':True,'unique_canonical30':True,'full_direct_sum':True},
         'weight2_quotient9_complete':True,'modes':modes,
         'actual_finite_shears_verified_with_complete_quadratic_weight_bound':True,
         'linear_bivector_curl_map_rank3_kernel6':True,
         'smoothness_only_used_to_deduce_constant_Xf_zero_implies_X_zero':True,
         'K0_gradient_quadratic_kernel12_and_basebase_weight2_give_rank_at_most2':True,
         'finite_shear_argument':'X=zY+wZ has square-zero linear matrix; expX=Id+X. Weight1 is its actual orbit derivative modulo only explicit Euler gauge. Acting by the inverse shear removes weight1, and global quadratic degree bounds the remainder at weight2; injectivity then removes it.',
         'no_component_count_or_global_exhaustion_claim':True,
         'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in inputs},
         'reviewer_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    target=ROOT/'reports/CLASSIFICATION_PLANE_SMOOTH_CUBIC_RANK_TWO_REVIEW_20261107.json'
    assert not target.exists();target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()

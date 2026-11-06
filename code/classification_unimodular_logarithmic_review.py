#!/usr/bin/env python3
"""Independent complete L112 certificate review over Fraction(u), without CAS."""
import ast
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
import hashlib,json,resource,time
from pathlib import Path

def trim(p):
    p=list(p)
    while p and not p[-1]:p.pop()
    return tuple(p)
def ua(p,q):return trim([(p[i] if i<len(p) else F(0))+(q[i] if i<len(q) else F(0)) for i in range(max(len(p),len(q)))])
def um(p,q):
    a=[F(0)]*max(0,len(p)+len(q)-1)
    for i,c in enumerate(p):
        for j,d in enumerate(q):a[i+j]+=c*d
    return trim(a)
def ud(p,q):
    p=list(p);out=[F(0)]*max(0,len(p)-len(q)+1)
    while p and len(p)>=len(q):
        i=len(p)-len(q);c=p[-1]/q[-1];out[i]=c
        for j,d in enumerate(q):p[i+j]-=c*d
        p=list(trim(p))
    return trim(out),tuple(p)
def ug(p,q):
    while q:p,q=q,ud(p,q)[1]
    return tuple(v/p[-1] for v in p) if p else ()
class RF:
    def __init__(self,p=0,q=(F(1),)):
        if isinstance(p,RF):self.p,self.q=p.p,p.q;return
        p=trim(p if isinstance(p,(tuple,list)) else (F(p),));q=trim(q)
        assert q
        if not p:self.p,self.q=(),(F(1),);return
        g=ug(p,q);p=ud(p,g)[0];q=ud(q,g)[0];a=q[-1]
        self.p,self.q=tuple(v/a for v in p),tuple(v/a for v in q)
    def __add__(self,v):
        v=RF(v);return RF(ua(um(self.p,v.q),um(v.p,self.q)),um(self.q,v.q))
    __radd__=__add__
    def __neg__(self):return RF(tuple(-v for v in self.p),self.q)
    def __sub__(self,v):return self+-RF(v)
    def __rsub__(self,v):return RF(v)+-self
    def __mul__(self,v):
        v=RF(v);return RF(um(self.p,v.p),um(self.q,v.q))
    __rmul__=__mul__
    def __truediv__(self,v):
        v=RF(v);assert v.p;return RF(um(self.p,v.q),um(self.q,v.p))
    def __rtruediv__(self,v):return RF(v)/self
    def __bool__(self):return bool(self.p)
    def __eq__(self,v):
        v=RF(v);return self.p==v.p and self.q==v.q
    def __hash__(self):return hash((self.p,self.q))
    def __repr__(self):return str((tuple(map(str,self.p)),tuple(map(str,self.q))))

N=11;ZERO=(0,)*N;NAMES=('x','y','z','w')+tuple('n'+str(i) for i in range(5))+('a0','a1')
U=RF((F(0),F(1)))
def atom(v=1,k=None):
    e=[0]*N
    if k is not None:e[k]=1
    return {tuple(e):RF(v)} if v else {}
def add(*ps):
    a={}
    for p in ps:
        for e,c in p.items():a[e]=a.get(e,RF())+c
    return {e:c for e,c in a.items() if c}
def scale(p,c):return {e:c*v for e,v in p.items() if c*v}
def mul(p,q):
    a={}
    for e,c in p.items():
        for f,d in q.items():
            g=tuple(v+w for v,w in zip(e,f));a[g]=a.get(g,RF())+c*d
    return {e:c for e,c in a.items() if c}
def der(p,k):
    a={}
    for e,c in p.items():
        if e[k]:
            f=list(e);f[k]-=1;a[tuple(f)]=c*e[k]
    return a
def parse(s):
    def walk(v):
        if isinstance(v,ast.Constant):return atom(F(v.value))
        if isinstance(v,ast.Name):return atom(U) if v.id=='u' else atom(1,NAMES.index(v.id))
        if isinstance(v,ast.UnaryOp):return scale(walk(v.operand),-1 if isinstance(v.op,ast.USub) else 1)
        if isinstance(v,ast.BinOp):
            a=walk(v.left)
            if isinstance(v.op,ast.Pow):
                b=atom()
                for _ in range(v.right.value):b=mul(b,a)
                return b
            b=walk(v.right)
            if isinstance(v.op,ast.Add):return add(a,b)
            if isinstance(v.op,ast.Sub):return add(a,scale(b,-1))
            if isinstance(v.op,ast.Mult):return mul(a,b)
            if isinstance(v.op,ast.Div):assert set(b)=={ZERO};return scale(a,1/b[ZERO])
        raise ValueError(ast.dump(v))
    return walk(ast.parse(str(s),mode='eval').body)
def scalar(s):
    p=parse(s);assert not p or set(p)=={ZERO};return p.get(ZERO,RF())
def matrix(m):return [[scalar(v) for v in r] for r in m]
def transpose(m):return list(map(list,zip(*m)))
def rank(m):
    b=[list(r) for r in m];r=0
    for j in range(len(b[0])):
        k=next((k for k in range(r,len(b)) if b[k][j]),None)
        if k is None:continue
        b[k],b[r]=b[r],b[k];a=b[r][j];b[r]=[v/a for v in b[r]]
        for k in range(r+1,len(b)):
            a=b[k][j]
            if a:b[k]=[v-a*w for v,w in zip(b[k],b[r])]
        r+=1
    return r
def det(m):
    b=[list(r) for r in m];out=RF(1)
    for j in range(len(b)):
        k=next((k for k in range(j,len(b)) if b[k][j]),None)
        if k is None:return RF()
        if k!=j:b[k],b[j]=b[j],b[k];out=-out
        a=b[j][j];out=out*a
        for k in range(j+1,len(b)):
            c=b[k][j]/a
            if c:
                for l in range(j+1,len(b)):b[k][l]=b[k][l]-c*b[j][l]
    return out
def submatrix(m,rr,cc):return [[m[i][j] for j in cc] for i in rr]
def annihilate(m,k):return all(not sum((a*b for a,b in zip(r,c)),RF()) for r in m for c in transpose(k))
PAIRS=list(combinations(range(4),2))
def pf(p):return add(mul(p[0],p[5]),scale(mul(p[1],p[4]),-1),mul(p[2],p[3]))
def cross(p,q):return add(pf([add(a,b) for a,b in zip(p,q)]),scale(pf(p),-1),scale(pf(q),-1))
def divergence(p):
    a={}
    for ij,v in zip(PAIRS,p):a[ij]=v;a[ij[::-1]]=scale(v,-1)
    return [add(*(der(a.get((i,j),{}),i) for i in range(4))) for j in range(4)]
def xycoeff(p):
    out={}
    for e,c in p.items():
        key=e[:4];f=(0,)*4+e[4:];out.setdefault(key,{})[f]=c
    return out
def linearvec(p):return [der(v,j).get(ZERO,RF()) for v in p for j in range(4)]
def roots_supported(poly,allowed):
    found=set();p=poly
    for r in allowed:
        fac=(-r,F(1))
        while len(p)>1:
            q,rem=ud(p,fac)
            if rem:break
            p=q;found.add(r)
    assert len(p)<=1,('uncovered denominator/pivot root',p)
    return found

def main():
    global U
    resource.setrlimit(resource.RLIMIT_CPU,(110,111));resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    started=time.monotonic();root=Path(__file__).resolve().parents[1]
    path=root/'reports/CLASSIFICATION_UNIMODULAR_LOGARITHMIC_COMPLETE_20261107.json'
    data=json.loads(path.read_text());records=[];coverage={}
    for c in data['cases']:
        begin=time.monotonic();U=RF((F(0),F(1))) if c['residue']=='generic u' else RF(F(c['residue']))
        q=parse(c['quadric']);x,y=atom(1,0),atom(1,1);dq=[der(q,i) for i in range(4)]
        gamma=[add(scale(q,(U-1)*int(i==0 and j==1)),scale(mul(y,add(scale(dq[j],int(i==0)),scale(dq[i],-int(j==0)))),U),mul(x,add(scale(dq[j],int(i==1)),scale(dq[i],-int(j==1))))) for i,j in PAIRS]
        Q=[gamma[5],scale(gamma[4],-1),gamma[3],gamma[2],scale(gamma[1],-1),gamma[0]]
        assert Q==list(map(parse,c['highest_bivector'])) and not pf(Q) and not any(divergence(Q))
        columns=[];mon3=set()
        for i in range(6):
            for j in range(4):
                L=[atom(1,j) if k==i else {} for k in range(6)];cr=cross(Q,L);mon3.update(cr)
                columns.append(([v.get(ZERO,RF()) for v in divergence(L)],cr))
        mon3=sorted(mon3,reverse=True)
        M=transpose([d+[cr.get(e,RF()) for e in mon3] for d,cr in columns])
        assert M==matrix(c['full_linear_matrix'])
        K=matrix(c['full_kernel']);r=c['rank'];assert annihilate(M,K) and rank(K)==24-r
        pivot=det(submatrix(M,c['pivot_rows'],c['pivot_columns']));assert pivot and pivot==scalar(c['pivot_determinant']) and len(c['pivot_rows'])==r
        T=transpose([linearvec([der(p,j) for p in Q]) for j in range(4)])
        assert T==matrix(c['translation_basis']) and annihilate(M,T)
        tr=c['translation_rank'];assert rank(T)==tr
        trrows=c['translation_cut_rows'];trcols=[]
        for j in range(4):
            if rank(submatrix(T,trrows,trcols+[j]))>len(trcols):trcols.append(j)
        trminor=det(submatrix(T,trrows,trcols));assert trminor==scalar(c['translation_minor']) and trminor
        Ln=list(map(parse,c['normal_linear_bivector']));nd=c['normal_dimension']
        NK=transpose([linearvec([der(p,4+i) for p in Ln]) for i in range(nd)])
        assert nd+tr==24-r and annihilate(M,NK) and rank(NK)==nd
        assert all(not v for i in trrows for v in NK[i])
        assert rank([a+b for a,b in zip(T,NK)])==24-r
        mon2=[]
        for i,j in combinations_with_replacement(range(4),2):
            e=[0]*N;e[i]+=1;e[j]+=1;mon2.append(tuple(e))
        G=transpose([[cross(Q,[atom() if k==i else {} for k in range(6)]).get(e,RF()) for e in mon2] for i in range(6)])
        assert G==matrix(c['constant_equation_matrix'])
        gr,gc=c['constant_pivot_rows'],c['constant_pivot_columns'];gm=det(submatrix(G,gr,gc))
        assert gm and gm==scalar(c['constant_minor']) and rank(G)==len(gc)
        C=list(map(parse,c['constant_source']));free=[i for i in range(6) if i not in gc]
        assert all(C[i]==atom(1,9+j) for j,i in enumerate(free))
        Pf=pf([add(a,b,d) for a,b,d in zip(Q,Ln,C)]);coeff=xycoeff(Pf)
        assert all(mon2[i][:4] not in coeff for i in gr)
        actual={tuple(sorted(p.items())) for e,p in coeff.items() if sum(e)>0}
        expected={tuple(sorted(parse(p).items())) for p in c['remaining_full_Pfaffian_equations']}
        assert actual==expected and coeff.get((0,)*4,{})==parse(c['constant_Pfaffian'])
        if c['name']!='rank0_full_mixed':assert not parse(c['constant_Pfaffian'])
        else:
            # Under the independently read support implication n0=n1=0,
            # both fiber columns have the single 34 image, including translations.
            L0=[{e:v for e,v in p.items() if not e[4] and not e[5]} for p in Ln]
            assert all(not der(p,j) for p in L0[:5] for j in (2,3))
            assert all(not der(der(p,k),j) for p in Q[:5] for k in range(4) for j in (2,3))
        if c['residue']=='generic u':
            allowed={F(0),F(1)}|{F(d['residue']) for d in data['cases'] if d['name']==c['name'] and d['residue']!='generic u'}
            found=set()
            for v in [pivot,trminor,gm]:found |= roots_supported(v.p,allowed)|roots_supported(v.q,allowed)
            for row in K+NK:
                for v in row:found |= roots_supported(v.q,allowed)
            for p in C:
                for v in p.values():found |= roots_supported(v.q,allowed)
            assert found-{F(0),F(1)}==allowed-{F(0),F(1)}
            coverage[c['name']]=list(map(str,sorted(found-{F(0),F(1)})))
        records.append({'name':c['name'],'residue':c['residue'],'full_linear_rank':r,'translation_rank':tr,'normal_dimension':nd,'constant_rank':len(gc),'all_original_coefficient_matrices_and_complete_kernels':True,'true_translation_slice_and_inverse_constant_solve':True,'all_positive_degree_Pfaffian_and_constant_coefficients':True,'whole_mixed_primitive_rank_bound_after_support_reduction':c['name']=='rank0_full_mixed','all_passed':True,'seconds':time.monotonic()-begin})
        print(c['name'],c['residue'],'PASS',round(time.monotonic()-begin,3),flush=True)
    assert len(records)==16 and len(coverage)==5
    result={'task_id':'G','all_passed':True,'method':'Independent standard-library Fraction(u) and sparse multivariate polynomial arithmetic; no producer or CAS imported. Whole generic identities and every exceptional residue, not point interpolation.','cases':records,'exceptional_coverage':coverage,'input_sha256':{str(path.relative_to(root)):hashlib.sha256(path.read_bytes()).hexdigest()},'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'All sixteen coefficient certificates. Geometric reduction and support implication n0=n1=0 are separately read proofs; parameter-boundary exhaustion is not claimed here.','seconds':time.monotonic()-started}
    (root/'reports/CLASSIFICATION_UNIMODULAR_LOGARITHMIC_REVIEW_20261107.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'cases':16,'seconds':result['seconds']}))
if __name__=='__main__':main()

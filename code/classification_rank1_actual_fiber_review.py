#!/usr/bin/env python3
"""Independent Fraction/Laurent reconstruction of the universal same-jet curve."""
import ast
from fractions import Fraction as F
from itertools import combinations
import hashlib
import json
from pathlib import Path
import resource
import time

resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
resource.setrlimit(resource.RLIMIT_CPU,(55,56))
started=time.monotonic()
N=12
names=('x','y','z','w','a','b','c','k','A','B','C','lam')
def const(a):return {(0,)*N:F(a)} if a else {}
def add(*ps):
    out={}
    for p in ps:
        for m,a in p.items():out[m]=out.get(m,F(0))+a
    return {m:a for m,a in out.items() if a}
def scale(p,a):return {m:a*b for m,b in p.items() if a*b}
def mul(p,q):
    out={}
    for m,a in p.items():
        for n,b in q.items():
            e=tuple(x+y for x,y in zip(m,n));out[e]=out.get(e,F(0))+a*b
    return {m:a for m,a in out.items() if a}
def power(p,k):
    if k<0:
        assert len(p)==1
        m,a=next(iter(p.items()));return {tuple(e*k for e in m):a**k}
    out=const(1)
    for _ in range(k):out=mul(out,p)
    return out
def var(i):
    m=[0]*N;m[i]=1;return {tuple(m):F(1)}
def der(p,i):
    out={}
    for m,a in p.items():
        if m[i]:
            e=list(m);e[i]-=1;out[tuple(e)]=a*m[i]
    return out
def sub(p,mapping):
    out={}
    for m,a in p.items():
        term=const(a)
        for i,k in enumerate(m):
            if k:term=mul(term,power(mapping.get(i,var(i)),k))
        out=add(out,term)
    return out
def parse(text):
    def walk(node):
        if isinstance(node,ast.Constant) and isinstance(node.value,int):return const(node.value)
        if isinstance(node,ast.Name) and node.id in names:return var(names.index(node.id))
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return scale(walk(node.operand),-1)
        if isinstance(node,ast.BinOp):
            a=walk(node.left)
            if isinstance(node.op,ast.Pow):
                assert isinstance(node.right,ast.Constant) and isinstance(node.right.value,int)
                return power(a,node.right.value)
            b=walk(node.right)
            if isinstance(node.op,ast.Add):return add(a,b)
            if isinstance(node.op,ast.Sub):return add(a,scale(b,-1))
            if isinstance(node.op,ast.Mult):return mul(a,b)
            if isinstance(node.op,ast.Div):return mul(a,power(b,-1))
        raise ValueError(ast.dump(node))
    return walk(ast.parse(text,mode='eval').body)

x,y,z,w,a,b,c,k,A,B,C,lam=[var(i) for i in range(N)]
one=const(1);half=F(1,2)
zeta=add(z,const(-1));ell=add(mul(add(b,c),zeta),scale(mul(a,w),-1))
K=add(k,mul(A,power(zeta,2)),mul(B,mul(zeta,w)),mul(C,power(w,2)))
KL=add(mul(lam,K),scale(mul(mul(lam,add(lam,const(-1))),power(ell,2)),half))
aL=mul(lam,a);cL=add(mul(lam,add(b,c)),scale(b,-1))
U_w=add(b,mul(cL,z));H=mul(z,KL)
Fp=add(mul(mul(aL,z),der(KL,2)),mul(U_w,der(KL,3)),mul(add(aL,one),KL))
r=add(y,scale(power(x,2),-half))
entries=[scale(add(Fp,mul(add(aL,one),r)),-1),mul(aL,z),U_w,
         add(mul(mul(aL,x),z),der(H,3)),add(mul(x,U_w),scale(r,-1),scale(der(H,2),-1)),z]
pairs=list(combinations(range(4),2));P={}
for ij,value in zip(pairs,entries):P[ij]=value;P[ij[::-1]]=scale(value,-1)
frame=[[one if i==j else {} for j in range(4)] for i in range(4)]
frame[0][2]=mul(add(one,scale(lam,-1)),add(b,c))
frame[0][3]=scale(mul(add(one,scale(lam,-1)),a),-1)
frame[1][1]=power(lam,-1)
inverse={0:add(x,scale(mul(add(one,scale(lam,-1)),ell),-1)),1:mul(lam,y)}
Q={}
for i,j in pairs:
    value=add(*(mul(mul(frame[i][u],frame[j][v]),P.get((u,v),{})) for u in range(4) for v in range(4)))
    Q[i,j]=sub(value,inverse);Q[j,i]=scale(Q[i,j],-1)
packet_path=Path('reports/CLASSIFICATION_RANK1_LONG_ACTUAL_FIBER_20261107.json')
packet=json.loads(packet_path.read_text())
for ij in pairs:assert Q[ij]==parse(packet['full_curve_brackets'][','.join(map(str,ij))])
for i,j,k0 in combinations(range(4),3):
    jac=add(*(mul(Q.get((u,l),{}),der(Q.get((v,t),{}),l))
              for u,v,t in [(i,j,k0),(j,k0,i),(k0,i,j)] for l in range(4)))
    assert not jac
base={0:{},1:{},2:one,3:{}}
for ij in pairs:
    difference=add(Q[ij],scale(sub(Q[ij],{11:one}),-1))
    assert not sub(difference,base)
    for v in range(4):assert not sub(der(difference,v),base)
    assert all(sum(m[:4])<=2 for m in Q[ij])
derivative=sub(der(Q[1,3],11),{11:one})
target=[0]*N;target[0]=2
coeff={tuple([0]*4+list(m[4:])):v for m,v in derivative.items() if m[:4]==(2,0,0,0)}
assert coeff==const(-half)
source=Path(__file__).read_text()
out={'status':'PASS_INDEPENDENT_UNIVERSAL_SAME_JET_CURVE','arithmetic':'Fraction sparse Laurent polynomial; no producer/SymPy import','input_sha256':hashlib.sha256(packet_path.read_bytes()).hexdigest(),
     'full_affine_source_reconstructed_from_formula':True,'exact_affine_pushforward_matches_all_six_report_brackets':True,
     'all_four_full_parameter_Jacobi_zero':True,'all_thirty_zero_firstjet_differences_zero':True,
     'all_original_coordinate_degrees_at_most_two':True,'universal_nonconstant_derivative':'-1/2',
     'scope':'Exact universal curve arithmetic. Whole-source normalization and dimension argument separately accepted by two geometry readers; no numerical point rank is used as an upper bound.',
     'seconds':time.monotonic()-started,'producer_sha256':hashlib.sha256(source.encode()).hexdigest(),'producer_source':source}
Path('reports/CLASSIFICATION_RANK1_LONG_ACTUAL_FIBER_REVIEW_20261107.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='producer_source'},indent=2))

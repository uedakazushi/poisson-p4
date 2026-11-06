#!/usr/bin/env python3
"""Independent Fraction reconstruction of the unique singular Kronecker pencil linear system."""
import ast
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
import hashlib
import json
from pathlib import Path
import resource
import time

resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
resource.setrlimit(resource.RLIMIT_CPU,(55,56))
started=time.monotonic()
N=8
names=('x1','x2','x3','x4','r0','r1','r2','r3')
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


ROOT=Path(__file__).resolve().parents[1];PAIRS=list(combinations(range(4),2));TRIPLES=list(combinations(range(4),3))
def wedge(a,b):return add(*(scale(mul(a[i],b[j]),sgn) for i,j,sgn in [(0,5,1),(5,0,1),(1,4,-1),(4,1,-1),(2,3,1),(3,2,1)]))
def rank(A):
    if not A:return 0
    B=[[F(v) for v in row] for row in A];r=0
    for j in range(len(B[0])):
        k=next((k for k in range(r,len(B)) if B[k][j]),None)
        if k is None:continue
        B[r],B[k]=B[k],B[r];d=B[r][j];B[r]=[v/d for v in B[r]]
        for k in range(r+1,len(B)):
            d=B[k][j]
            if d:B[k]=[u-d*v for u,v in zip(B[k],B[r])]
        r+=1
    return r

def coefficient_vector(p,degree):
    ex=[]
    for inds in combinations_with_replacement(range(4),degree):
        e=[0]*N
        for i in inds:e[i]+=1
        ex.append(tuple(e))
    assert all(e in ex for e in p)
    return [p.get(e,F(0)) for e in ex]
path=ROOT/'reports/CLASSIFICATION_UNIMODULAR_SINGULAR_PENCIL_20261107.json';saved=json.loads(path.read_text());results=[]
for case in saved['cases']:
    f=parse(case['nondegenerate_quadric']);g=parse(case['second_quadric'])
    df=[der(f,i) for i in range(4)];dg=[der(g,i) for i in range(4)]
    gamma=[add(mul(df[i],dg[j]),scale(mul(df[j],dg[i]),-1)) for i,j in PAIRS]
    assert gamma==list(map(parse,case['highest_six_coefficients']))
    def cond(beta):
        b={ij:v for ij,v in zip(PAIRS,beta)}
        c=[add(der(b[j,k],i),scale(der(b[i,k],j),-1),der(b[i,j],k)) for i,j,k in TRIPLES]
        return [p.get((0,)*N,F(0)) for p in c]+coefficient_vector(wedge(beta,gamma),3)
    cols=[]
    for at in range(6):
        for k in range(4):
            beta=[{} for _ in PAIRS];beta[at]=var(k);cols.append(cond(beta))
    D=list(map(list,zip(*cols)));assert rank(D)==16
    B=[]
    for at in range(8):
        u=[const(int(at==i)) for i in range(4)];v=[const(int(at==4+i)) for i in range(4)]
        beta=[add(mul(df[i],v[j]),scale(mul(df[j],v[i]),-1),mul(u[i],dg[j]),scale(mul(u[j],dg[i]),-1)) for i,j in PAIRS]
        assert cond(beta)==[0]*24
        B.append([der(p,k).get((0,)*N,F(0)) for p in beta for k in range(4)])
    assert rank(B)==8
    Ccols=[]
    for at in range(6):
        om=[const(int(at==i)) for i in range(6)]
        Ccols.append(coefficient_vector(wedge(om,gamma),2))
    C=list(map(list,zip(*Ccols)));K=[[F(v) for v in row] for row in case['constant_orthogonal_kernel']]
    assert rank(C)+len(K)==6 and rank(K)==len(K)
    assert all(sum(a*b for a,b in zip(row,k))==0 for row in C for k in K)
    om=[add(*(scale(var(4+k),K[k][i]) for k in range(len(K)))) for i in range(6)]
    assert scale(wedge(om,om),F(1,2))==parse(case['constant_orthogonal_Pfaffian'])
    assert bool(wedge(om,om))==case['constant_nondegenerate_possible']
    results.append({'name':case['name'],'full24_coefficient_rank':16,'whole_eight_parameter_kernel_verified':True,
                    'full_six_constant_coefficient_rank':rank(C),'constant_kernel_dimension':len(K),
                    'entire_constant_Pfaffian_verified':True,'constant_nondegenerate_possible':case['constant_nondegenerate_possible']})
assert [x['name'] for x in results if x['constant_nondegenerate_possible']]==[]
result={'task_id':'G','all_passed':True,'method':'Independent standard-library Fraction reconstruction, no producer or SymPy arithmetic imported',
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'cases':results,
        'scope':'Full eight-dimensional linear kernel and one-dimensional constant kernel dx1 wedge dx2. The constant-Pfaffian condition forces u3=u4=0, hence the constant form J is degenerate and the symplectic source is empty. Common-kernel pencils are handled by the separate Hilbert-Burch proof.',
        'seconds':time.monotonic()-started}
(ROOT/'reports/CLASSIFICATION_UNIMODULAR_SINGULAR_PENCIL_REVIEW_20261107.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'all_passed':True,'types':len(results),'seconds':result['seconds']}))

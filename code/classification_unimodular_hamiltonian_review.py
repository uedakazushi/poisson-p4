#!/usr/bin/env python3
"""Independent Fraction reconstruction of all six full affine-Hamiltonian sources."""
import ast
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, permutations
import hashlib
import json
from pathlib import Path
import resource
import time

resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
resource.setrlimit(resource.RLIMIT_CPU,(55,56))
started=time.monotonic()
N=14
names=('x','y','z','w')+tuple('n'+str(i) for i in range(10))
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


ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(combinations(range(4),2))
M=[const(1)]+[var(i) for i in range(4)]+[mul(var(i),var(j)) for i,j in combinations_with_replacement(range(4),2)]
MEX=[next(iter(p)) for p in M]
def anti(d):return {**d,**{(j,i):scale(v,-1) for (i,j),v in d.items()}}
def pf(P):return add(mul(P[0,1],P[2,3]),scale(mul(P[0,2],P[1,3]),-1),mul(P[0,3],P[1,2]))
def conditions(P):
    out=[add(*(der(P.get((i,j),{}),i) for i in range(4))) for j in range(4)]
    out.append({e:c for e,c in pf(P).items() if any(e[:4])})
    return out

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

def pdet(A):
    out={}
    for perm in permutations(range(len(A))):
        term=const((-1)**sum(perm[i]>perm[j] for i in range(len(A)) for j in range(i+1,len(A))))
        for i,j in enumerate(perm):term=mul(term,A[i][j])
        out=add(out,term)
    return out

path=ROOT/'reports/CLASSIFICATION_UNIMODULAR_AFFINE_HAMILTONIAN_20261107.json'
saved=json.loads(path.read_text())
cases=[('nil31_long',['y','z','1','0'],3),('nil31_short',['y','z','0','1'],2),
       ('nil22_top',['y','1','w','0'],3),('nil211_long',['y','1','0','0'],3),
       ('nil211_short_chain_covector',['y','0','1','0'],1),('nil211_short_other_covector',['y','0','1','0'],3)]
results=[]
for (name,vs,ell),old in zip(cases,saved['cases']):
    assert old['name']==name
    V=list(map(parse,vs));assert not V[ell]
    fixed={ij:V[ij[1]] if ij[0]==ell else scale(V[ij[0]],-1) if ij[1]==ell else {} for ij in PAIRS}
    P0=anti(fixed);cols=[(ij,e) for ij in PAIRS if ell not in ij for e in MEX];assert len(cols)==45
    C0=conditions(P0);Dcols=[]
    for ij,e in cols:
        p=dict(fixed);p[ij]={e:F(1)}
        Dcols.append([add(a,scale(b,-1)) for a,b in zip(conditions(anti(p)),C0)])
    rows=sorted({(i,e) for con in [C0]+Dcols for i,p in enumerate(con) for e in p})
    D=[[con[i].get(e,F(0)) for con in Dcols] for i,e in rows]
    rhs=[-C0[i].get(e,F(0)) for i,e in rows]
    r=rank(D);ar=rank([row+[b] for row,b in zip(D,rhs)])
    if name=='nil31_short':
        assert r==39 and ar==40 and old['status']=='EMPTY_LINEAR_SOURCE'
        results.append({'name':name,'coefficient_rank':r,'augmented_rank':ar,'all45_source_empty':True});continue
    P=anti({tuple(map(int,k.split(','))):parse(v) for k,v in old['full_six_brackets'].items()})
    assert all(P.get((ell,j),{})==V[j] for j in range(4))
    assert all(not p for p in conditions(P))
    assert pf(P)==parse(old['constant_Pfaffian']) and pf(P)
    for i,j,k in combinations(range(4),3):
        J=add(*(mul(P.get((a,l),{}),der(P.get((b,c),{}),l)) for a,b,c in [(i,j,k),(j,k,i),(k,i,j)] for l in range(4)))
        assert not J
    d=old['source_dimension'];assert r==ar==45-d
    dirs=[]
    for k in range(d):
        v=[]
        for ij,e in cols:
            ex=list(e);ex[4+k]=1;v.append(P[ij].get(tuple(ex),F(0)))
        dirs.append(v)
    assert rank(dirs)==d
    assert all(sum(e[4:])<=1 for ij in PAIRS for e in P[ij])
    L=[[sub(der(P[ij],k),{i:{} for i in range(4)}) for k in range(4)] for ij in PAIRS]
    assert L==[[parse(v) for v in row] for row in old['linear_evaluation_matrix']]
    assert all(not pdet([L[i] for i in rr]) for rr in combinations(range(6),4))
    results.append({'name':name,'full_coefficient_rank':r,'augmented_rank':ar,'source_parameter_rank':d,
                    'entire_affine_linear_solution_verified':True,'full_divergence_and_Pfaffian_conditions_zero':True,
                    'all_four_Jacobi_identically_zero':True,'Hamiltonian_identity':True,'all15_primitive_evaluation_minors_identically_zero':True})
result={'task_id':'G','all_passed':True,'method':'Independent standard-library Fraction sparse polynomial reconstruction; no SymPy or producer arithmetic import',
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'all_six_types':results,
        'scope':'All45 unknown coefficients retained. Exact coefficient rank plus universal solution and independent parameter directions prove completeness. Geometric exhaustion of six pairs and semisimple/constant cases independently read in the focused note.',
        'seconds':time.monotonic()-started}
(ROOT/'reports/CLASSIFICATION_UNIMODULAR_AFFINE_HAMILTONIAN_REVIEW_20261107.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'all_passed':True,'six_types':len(results),'seconds':result['seconds']}))

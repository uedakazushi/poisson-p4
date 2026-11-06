#!/usr/bin/env python3
"""Independent Fraction/Laurent verification of the full q4-nil22 orbit and smooth-quadric deformation."""
import ast
from fractions import Fraction as F
from itertools import combinations, permutations
import hashlib
import json
from pathlib import Path
import resource
import time

resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
resource.setrlimit(resource.RLIMIT_CPU,(55,56))
started=time.monotonic()
N=10
names=('h','x1','x2','x3','x4','t','lam','n2','n3','n6')
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
        if isinstance(node,ast.Name) and node.id=='n0':return power(var(6),-2)
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


ROOT=Path(__file__).resolve().parents[1];PAIRS=list(combinations(range(5),2))
def anti(d):return {**d,**{(j,i):scale(p,-1) for (i,j),p in d.items()}}
def decode(d):return anti({tuple(map(int,k.split(','))):parse(v) for k,v in d.items()})
def jac(P):
    J={ijk:add(*(mul(P.get((a,l),{}),der(P.get((b,c),{}),l)) for a,b,c in [ijk,(ijk[1],ijk[2],ijk[0]),(ijk[2],ijk[0],ijk[1])] for l in range(5))) for ijk in combinations(range(5),3)}
    return {q:add(*(scale(mul(var(q[k]),J[q[:k]+q[k+1:]]),(-1)**k) for k in range(4))) for q in combinations(range(5),4)}
def pf(P):
    out={}
    for m in range(5):
        i,j,k,l=[v for v in range(5) if v!=m]
        p=add(mul(P.get((i,j),{}),P.get((k,l),{})),scale(mul(P.get((i,k),{}),P.get((j,l),{})),-1),mul(P.get((i,l),{}),P.get((j,k),{})))
        out=add(out,scale(mul(var(m),p),F((-1)**m,5)))
    return out

def pdet(A):
    out={}
    for perm in permutations(range(len(A))):
        v=const((-1)**sum(perm[i]>perm[j] for i in range(len(A)) for j in range(i+1,len(A))))
        for i,j in enumerate(perm):v=mul(v,A[i][j])
        out=add(out,v)
    return out

paths=['reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL22_CURVE_H_20261107.json',
       'reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL22_CANONICAL_20261107.json',
       'reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL_COMPONENTS_20261107.json']
curve,can,components=[json.loads((ROOT/p).read_text()) for p in paths]
P=decode(curve['family']);P0=decode(can['canonical_brackets'])
assert all(sub(P.get(ij,{}),{5:{}})==P0.get(ij,{}) for ij in PAIRS)
assert all(not p for p in jac(P).values())
D=pf(P);assert D==parse(curve['Euler_contracted_Pfaffian'])
q=parse(curve['quadratic_factor']);H=[[der(der(q,i),j) for j in range(5)] for i in range(5)]
assert pdet(H)==parse(curve['quadratic_hessian_determinant'])==parse('1394713760400000000*t**5')
# The published factorization is a product of two parenthesized factors divided
# by a nonzero rational scalar; obtain and check the cubic without CAS factoring.
node=ast.parse(curve['Euler_contracted_Pfaffian'],mode='eval').body
assert isinstance(node,ast.BinOp) and isinstance(node.op,ast.Div)
assert isinstance(node.left,ast.BinOp) and isinstance(node.left.op,ast.Mult)
factors=[parse(ast.unparse(node.left.left)),parse(ast.unparse(node.left.right))]
assert q in factors;c=factors[1] if factors[0]==q else factors[0]
assert all(sum(ex[:5])==3 for ex in c)
c0=sub(c,{5:{}})
assert c0==parse('787320000*h**2*x3-787320000*h*x1*x2+157464000*x2**3')
assert der(c0,3)==parse('787320000*h**2')
assert sub(c0,{0:{}})==parse('157464000*x2**3')
# The entire original open four-parameter source is the orbit of P0.
G=[[parse(v) for v in row] for row in can['coordinate_change']]
h,x1,x2,x3,x4,t,lam,n2,n3,n6=[var(i) for i in range(N)];k=scale(n3,F(3,5))
inverse={0:h,1:x1,2:mul(add(x2,mul(k,h)),power(lam,-1)),3:mul(add(x3,mul(k,x1)),power(lam,-1)),
         4:add(x4,scale(mul(add(n2,mul(k,n6)),h),F(-3,10)),scale(mul(n6,x2),F(-3,5)))}
for i in range(5):
    assert add(*(mul(G[i][j],inverse[j]) for j in range(5)))==var(i)
assert pdet(G)==power(lam,2)
case=next(c for c in components['cases'] if c['source']=='nil22_nonzero')
for i,j in PAIRS:
    pushed=add(*(mul(mul(G[i][a],G[j][b]),P0.get((a,b),{})) for a in range(5) for b in range(5)))
    pushed=sub(pushed,inverse)
    expected=scale(parse(case['universal_rational_family'].get(str((i,j)),'0').replace('x0','h')),F(1,125))
    assert pushed==expected
result={'task_id':'E_reduced_nonnormal','all_passed':True,'method':'Independent standard-library Fraction sparse Laurent arithmetic; no SymPy/producer imports',
        'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'all350_curve_Jacobi_identically_zero':True,'actual_Euler_Pfaffian_full_factorization_verified':True,
        'quadratic_hessian_determinant':'1394713760400000000*t^5','cubic_at_zero_absolutely_irreducible_by_primitive_x3_linearity':True,
        'curve_endpoint_is_source_canonical_point':True,'entire_four_parameter_orbit_identity_in_all10_brackets':True,
        'coordinate_change_determinant':'lam^2','open_source_surjectivity':'n0=lam^-2 covers all n0!=0 over C',
        'scope':'The whole q4 nil22-nonzero source lies in MPS by the established smooth-quadric/irreducible-cubic full-support theorem and closedness. Full finite-slice ideal is not a dependency.',
        'seconds':time.monotonic()-started}
(ROOT/'reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL22_REVIEW_20261107.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'all_passed':True,'seconds':result['seconds'],'whole_source_MPS':True}))

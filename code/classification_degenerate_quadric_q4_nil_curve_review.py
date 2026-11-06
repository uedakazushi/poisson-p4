#!/usr/bin/env python3
"""Independent Fraction verification of the literal nil31-II Poisson curve."""
import ast
from fractions import Fraction as F
from itertools import combinations, permutations
import hashlib
import json
from pathlib import Path
import sys
import time

START=time.monotonic();N=7
def add(*polys):
    out={}
    for p in polys:
        for m,a in p.items():out[m]=out.get(m,F(0))+a
    return {m:a for m,a in out.items() if a}
def mul(p,q):
    out={}
    for m,a in p.items():
        for n,b in q.items():
            e=tuple(i+j for i,j in zip(m,n));out[e]=out.get(e,F(0))+a*b
    return {m:a for m,a in out.items() if a}
def scale(p,a):return {m:a*b for m,b in p.items() if a*b}
def term(a,i=None):
    e=[0]*N
    if i is not None:e[i]=1
    return {tuple(e):F(a)} if a else {}
def der(p,i):
    out={}
    for m,a in p.items():
        if m[i]:e=list(m);e[i]-=1;out[tuple(e)]=a*m[i]
    return out
def parse(text):
    names={f'x{i}':i for i in range(5)}|{'h':0,'t':5,'L':6}
    def walk(n):
        if isinstance(n,ast.Constant) and isinstance(n.value,int):return term(n.value)
        if isinstance(n,ast.Name):return term(1,names[n.id])
        if isinstance(n,ast.UnaryOp):return scale(walk(n.operand),-1 if isinstance(n.op,ast.USub) else 1)
        if isinstance(n,ast.BinOp):
            if isinstance(n.op,ast.Add):return add(walk(n.left),walk(n.right))
            if isinstance(n.op,ast.Sub):return add(walk(n.left),scale(walk(n.right),-1))
            if isinstance(n.op,ast.Mult):return mul(walk(n.left),walk(n.right))
            if isinstance(n.op,ast.Div) and isinstance(n.right,ast.Constant):return scale(walk(n.left),F(1,n.right.value))
            if isinstance(n.op,ast.Pow) and isinstance(n.right,ast.Constant):
                out=term(1)
                for _ in range(n.right.value):out=mul(out,walk(n.left))
                return out
        raise ValueError(ast.dump(n))
    return walk(ast.parse(text,mode='eval').body)
def tensor(data):
    P={}
    for key,value in data.items():
        i,j=ast.literal_eval(key);P[i,j]=parse(value);P[j,i]=scale(P[i,j],-1)
    return P
def entry(P,i,j):return P.get((i,j),{})
def divergence(P):return [add(*(der(entry(P,i,j),i) for i in range(5))) for j in range(5)]

nil22=len(sys.argv)>2 and sys.argv[2]=='nil22'
paths=['reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL_CURVE_20261107.json',
       'reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL_CURVE_RECOGNITION_20261107.json']
if nil22:paths=['reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL22_CURVE_H_20261107.json']
curve=json.loads(Path(paths[0]).read_text());P=tensor(curve['family'])
J={}
for i,j,k in combinations(range(5),3):
    J[i,j,k]=add(*(mul(entry(P,a,l),der(entry(P,b,c),l))
                   for a,b,c in [(i,j,k),(j,k,i),(k,i,j)] for l in range(5)))
for inds in combinations(range(5),4):
    assert not add(*(scale(mul(term(1,inds[k]),J[inds[:k]+inds[k+1:]]),(-1)**k) for k in range(4)))
if nil22:
    q=parse(curve['quadratic_factor']);H=[[der(der(q,i),j) for j in range(5)]for i in range(5)]
    determinant={}
    for perm in permutations(range(5)):
        value=term((-1)**sum(perm[i]>perm[j] for i in range(5) for j in range(i+1,5)))
        for i,j in enumerate(perm):value=mul(value,H[i][j])
        determinant=add(determinant,value)
    assert determinant==parse(curve['quadratic_hessian_determinant'])
    D={}
    for i in range(5):
        a,b,c,d=[j for j in range(5) if j!=i]
        pf=add(mul(entry(P,a,b),entry(P,c,d)),scale(mul(entry(P,a,c),entry(P,b,d)),-1),mul(entry(P,a,d),entry(P,b,c)))
        D=add(D,scale(mul(term(1,i),pf),(-1)**i))
    assert scale(D,F(1,5))==parse(curve['Euler_contracted_Pfaffian'])
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,field='QQ',
        method='Independent sparse Fraction and AST parser; no SymPy or producer imported.',
        all_parameter_full_projective_Jacobi_zero=True,actual_Euler_Pfaffian_factorization_checked=True,
        quadratic_Hessian_determinant_checked_by_120_term_determinant=True,
        quadratic_Hessian_determinant=curve['quadratic_hessian_determinant'],
        source_sha256={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths},
        independent_checker_source=Path(__file__).read_text(),seconds=time.monotonic()-START,
        scope='Literal universal curve and smooth-quadric recognition only; whole-source orbit reduction and MPS theorem remain separate geometric inputs.')
    with Path(sys.argv[1]).open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({'all_passed':True,'seconds':result['seconds']}));sys.exit(0)
recognition=json.loads(Path(paths[1]).read_text());U=tensor(recognition['unimodular_brackets'])
div=divergence(P)
for i,j in combinations(range(5),2):
    canonical=add(entry(P,i,j),scale(add(mul(term(1,i),div[j]),scale(mul(term(1,j),div[i]),-1)),F(-1,5)))
    assert canonical==entry(U,i,j)
assert all(not f for f in divergence(U))
q=parse('-t*x2**2+10*x0*x2+5*x1**2-5*x3**2')
X=list(map(parse,recognition['canonical_cofactor']))
for i in range(5):
    assert add(*(mul(entry(U,j,i),der(q,j)) for j in range(5)))==mul(q,X[i])
matrix=[[parse(v) for v in row] for row in recognition['canonical_cofactor_matrix']]
for i in range(5):assert add(*(mul(matrix[i][j],term(1,j)) for j in range(5)))==X[i]
char={}
for perm in permutations(range(5)):
    sign=(-1)**sum(perm[i]>perm[j] for i in range(5) for j in range(i+1,5))
    value=term(sign)
    for i,j in enumerate(perm):value=mul(value,add(term(1,6) if i==j else {},scale(matrix[i][j],-1)))
    char=add(char,value)
assert char==parse(recognition['characteristic_polynomial'])
D={}
for i in range(5):
    a,b,c,d=[j for j in range(5) if j!=i]
    pf=add(mul(entry(P,a,b),entry(P,c,d)),scale(mul(entry(P,a,c),entry(P,b,d)),-1),mul(entry(P,a,d),entry(P,b,c)))
    D=add(D,scale(mul(term(1,i),pf),(-1)**i))
assert scale(D,F(1,5))==parse(curve['Euler_contracted_Pfaffian'])
checker=Path(__file__).read_text()
result=dict(task_id='E_reduced_nonnormal',all_passed=True,field='QQ',
    method='Independent sparse Fraction and AST parser; no SymPy or producer imported.',
    all_parameter_full_projective_Jacobi_zero=True,unimodular_gauge_and_divergence_checked=True,
    all_five_q_cofactor_identities_checked=True,characteristic_polynomial_checked_by_120_term_determinant=True,
    actual_Euler_Pfaffian_factorization_checked=True,
    source_sha256={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths},
    independent_checker_source=checker,seconds=time.monotonic()-START,
    scope='Literal universal curve and regular-semisimple recognition identities only; whole-source orbit reduction and MPS theorem are separate geometric inputs.')
with Path(sys.argv[1]).open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')
print(json.dumps({'all_passed':True,'seconds':result['seconds']}))

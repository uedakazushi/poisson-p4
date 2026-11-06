#!/usr/bin/env python3
"""Independent exact QQ(i) review of rank2 semisimple and mixed smooth-quadric sources."""
import argparse
import ast
from fractions import Fraction as Rational
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

class F:
    """Exact Gaussian rational; two arguments denote a real rational fraction."""
    __slots__=('a','b')
    def __init__(self,value=0,denominator=None):
        if denominator is not None:self.a=Rational(value,denominator);self.b=Rational(0)
        elif isinstance(value,F):self.a=value.a;self.b=value.b
        elif isinstance(value,tuple):self.a=Rational(value[0]);self.b=Rational(value[1])
        else:self.a=Rational(value);self.b=Rational(0)
    def __add__(self,other):
        other=F(other);return F((self.a+other.a,self.b+other.b))
    __radd__=__add__
    def __neg__(self):return F((-self.a,-self.b))
    def __sub__(self,other):return self+-F(other)
    def __rsub__(self,other):return F(other)+-self
    def __mul__(self,other):
        other=F(other);return F((self.a*other.a-self.b*other.b,self.a*other.b+self.b*other.a))
    __rmul__=__mul__
    def __truediv__(self,other):
        other=F(other);d=other.a**2+other.b**2
        return F(((self.a*other.a+self.b*other.b)/d,(self.b*other.a-self.a*other.b)/d))
    def __rtruediv__(self,other):return F(other)/self
    def __bool__(self):return bool(self.a or self.b)
    def __eq__(self,other):
        other=F(other);return self.a==other.a and self.b==other.b
    def __str__(self):return str(self.a) if not self.b else '('+str(self.a)+'+('+str(self.b)+')*I)'

ROOT = Path(__file__).resolve().parents[1]
ZERO = (0,)*16
PAIRS = list(itertools.combinations(range(5), 2))
PARAMETERS = ['l0','l1','l2']+['m'+str(i) for i in range(8)]


def scalar(c):
    return {ZERO: F(c)} if c else {}


def variable(j):
    return {tuple(int(i == j) for i in range(16)): F(1)}


def add(*polys):
    out = {}
    for poly in polys:
        for mon, c in poly.items():
            out[mon] = out.get(mon, F(0))+c
    return {m: c for m, c in out.items() if c}


def scale(poly, c):
    return {m: c*v for m, v in poly.items() if c*v}


def multiply(a, b):
    out = {}
    for m, c in a.items():
        for n, d in b.items():
            word = tuple(x+y for x, y in zip(m, n))
            out[word] = out.get(word, F(0))+c*d
    return {m: c for m, c in out.items() if c}


def power(poly, n):
    out = scalar(1)
    for _ in range(n):
        out = multiply(out, poly)
    return out


def derivative(poly, j):
    return {tuple(v-int(i == j) for i, v in enumerate(m)): c*m[j]
            for m, c in poly.items() if m[j]}


def parse(text, extra=None):
    env = {f'x{i}': variable(i) for i in range(5)}
    env.update({name: variable(i+5) for i, name in enumerate(PARAMETERS)})
    env['I']=scalar(F((0,1)))
    env.update(extra or {})
    def read(node):
        if isinstance(node, ast.Constant):
            return scalar(node.value)
        if isinstance(node, ast.Name):
            return env[node.id]
        if isinstance(node, ast.UnaryOp):
            assert isinstance(node.op, (ast.USub, ast.UAdd))
            return scale(read(node.operand), -1 if isinstance(node.op, ast.USub) else 1)
        assert isinstance(node, ast.BinOp)
        if isinstance(node.op, ast.Pow):
            return power(read(node.left), ast.literal_eval(node.right))
        a, b = read(node.left), read(node.right)
        if isinstance(node.op, ast.Add): return add(a, b)
        if isinstance(node.op, ast.Sub): return add(a, scale(b, -1))
        if isinstance(node.op, ast.Mult): return multiply(a, b)
        assert isinstance(node.op, ast.Div) and set(b) == {ZERO}
        return scale(a, 1/b[ZERO])
    return read(ast.parse(str(text), mode='eval').body)


def entry(P, i, j):
    return P.get((i, j), {}) if i < j else scale(P.get((j, i), {}), -1) if i > j else {}


def jacobi(P):
    out = {}
    for i, j, k in itertools.combinations(range(5), 3):
        terms = []
        for a, b, c in ((i, j, k), (j, k, i), (k, i, j)):
            terms += [multiply(entry(P, a, r), derivative(entry(P, b, c), r)) for r in range(5)]
        out[(i, j, k)] = add(*terms)
    return out


def rank(rows):
    rows = [list(map(F, row)) for row in rows]
    if not rows: return 0
    at = 0
    for column in range(len(rows[0])):
        pivot = next((i for i in range(at, len(rows)) if rows[i][column]), None)
        if pivot is None: continue
        rows[at], rows[pivot] = rows[pivot], rows[at]
        rows[at] = [v/rows[at][column] for v in rows[at]]
        for i in range(at+1, len(rows)):
            if rows[i][column]:
                c = rows[i][column]
                rows[i] = [a-c*b for a, b in zip(rows[i], rows[at])]
        at += 1
        if at == len(rows): break
    return at


def substitute(poly, values):
    out = {}
    for mon, c in poly.items():
        term = scalar(c)
        for i, exponent in enumerate(mon):
            term = multiply(term, power(values.get(i, variable(i)), exponent))
        out = add(out, term)
    return out


def pfaffians(P):
    out = []
    for missing in range(5):
        a, b, c, d = [i for i in range(5) if i != missing]
        out.append(scale(add(multiply(entry(P,a,b), entry(P,c,d)),
                            scale(multiply(entry(P,a,c),entry(P,b,d)), -1),
                            multiply(entry(P,a,d),entry(P,b,c))), (-1)**missing))
    return out

def decode(data):return {tuple(map(int,key)):parse(value) for key,value in data.items()}

def coefficient_polynomials(polys):
    output=[]
    for p in polys:
        xmons=sorted(set(m[:5] for m in p))
        for xm in xmons:
            output.append({(0,)*5+m[5:]:c for m,c in p.items() if m[:5]==xm})
    return output

def polynomial_rank(polys):
    mons=sorted(set().union(*(set(p) for p in polys)))
    return rank([[p.get(m,F(0)) for m in mons] for p in polys])

def affine_system(polys,indices):
    rows=[];constants=[]
    for p in polys:
        assert all(sum(m)<=1 and all(not power or i in indices for i,power in enumerate(m)) for m in p)
        rows.append([p.get(tuple(int(j==i) for j in range(16)),F(0)) for i in indices])
        constants.append(p.get(ZERO,F(0)))
    return rank(rows),rank([r+[c] for r,c in zip(rows,constants)])

def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic()
    paths=['reports/CLASSIFICATION_SMOOTH_QUADRIC_RANK2_'+suffix+'_20261107.json'
           for suffix in ['LINEAR','SUPPORT','MIXED','MIXED_SUPPORT']]
    inputs=[json.loads((ROOT/path).read_text()) for path in paths]
    for saved in inputs:assert hashlib.sha256(saved['producer_source'].encode()).hexdigest()==saved['producer_sha256']
    linear,semisimple,mixed0,mixed=inputs
    assert linear['brackets']==semisimple['brackets']
    assert mixed0['brackets']==mixed['brackets']
    weights=(0,0,0,1,-1);q=parse('x0**2+x1**2+x2**2+x3*x4')
    xmons={n:sorted((e for e in itertools.product(range(n+1),repeat=5) if sum(e)==n),reverse=True) for n in [1,2,3]}
    invariant=[(ij,m) for ij in PAIRS for m in xmons[2] if sum(weights[k]*m[k] for k in range(5))==weights[ij[0]]+weights[ij[1]]]
    assert len(invariant)==46
    def linear_map(P):
        div=[add(*(derivative(entry(P,i,j),i) for i in range(5))) for j in range(5)]
        normal=[add(*(multiply(derivative(q,i),entry(P,i,j)) for i in range(5))) for j in range(5)]
        return div+normal
    rows=[(j,m+(0,)*11) for j in range(10) for m in xmons[1 if j<5 else 3]]
    cols=[]
    for ij,m in invariant:
        vals=linear_map({ij:{m+(0,)*11:F(1)}})
        cols.append([vals[j].get(e,F(0)) for j,e in rows])
    matrix=list(map(list,zip(*cols)))
    rhs_polys=[{} for _ in range(5)]+[scale(multiply(q,variable(j)),weights[j]) for j in range(5)]
    rhs=[rhs_polys[j].get(e,F(0)) for j,e in rows]
    assert rank(matrix)==rank([r+[b] for r,b in zip(matrix,rhs)])==35
    P=decode(semisimple['brackets'])
    assert all((ij,m[:5]) in invariant for ij,p in P.items() for m in p)
    assert linear_map(P)==rhs_polys
    def source_direction_rank(A,indices):
        columns=[]
        for i in indices:
            d={ij:substitute(derivative(p,i),{j:{} for j in range(5,16)}) for ij,p in A.items()}
            columns.append([d[ij].get(m+(0,)*11,F(0)) for ij in PAIRS for m in xmons[2]])
        return rank(columns)
    assert source_direction_rank(P,range(5,16))==11
    raw=coefficient_polynomials(jacobi(P).values())
    for p in raw:
        assert all(sum(m[5:8])==1 and sum(m[8:])<=1 for m in p)
    saved_conditions=list(map(parse,semisimple['residual_jacobi_coefficients']))
    assert polynomial_rank(raw)==polynomial_rank(saved_conditions)==polynomial_rank(raw+saved_conditions)
    orbit_checks=[]
    for name,values in [('zero',(0,0,0)),('nonisotropic',(1,0,0)),('isotropic',(1,F((0,1)),0))]:
        sub={5+i:scalar(v) for i,v in enumerate(values)}
        restricted=[substitute(p,sub) for p in raw]
        ranks=affine_system(restricted,list(range(8,16)))
        if name=='zero':
            assert ranks==(0,0)
            pp={ij:substitute(p,sub) for ij,p in P.items()}
            assert all(not p for p in pfaffians(pp))
        else:assert ranks==(7,8)
        orbit_checks.append({'ell_case':name,'coefficient_rank':ranks[0],'augmented_rank':ranks[1],
            'rank_at_most_two' if name=='zero' else 'empty':True})
    # Construct the mixed shift independently from the displayed constant nilpotent matrix.
    I=F((0,1));N=[[F(0) for _ in range(5)] for _ in range(5)]
    N[0][2]=-I;N[1][2]=F(1);N[2][0]=I;N[2][1]=F(-1)
    V=[add(*(scale(variable(j),N[i][j]) for j in range(5))) for i in range(5)]
    shift={ij:add(scale(q,-F(5,6)*N[ij[0]][ij[1]]),
        scale(add(multiply(variable(ij[0]),V[ij[1]]),scale(multiply(variable(ij[1]),V[ij[0]]),-1)),-F(1,3))) for ij in PAIRS}
    shifted={ij:add(P[ij],shift[ij]) for ij in PAIRS}
    expected=[{} for _ in range(5)]+[multiply(q,add(scale(variable(j),weights[j]),V[j])) for j in range(5)]
    assert linear_map(shifted)==expected
    def lie(A):
        return {ij:add(*(multiply(V[k],derivative(A[ij],k)) for k in range(5)),
            *(scale(entry(A,k,ij[1]),-N[ij[0]][k]) for k in range(5)),
            *(scale(entry(A,ij[0],k),-N[ij[1]][k]) for k in range(5))) for ij in PAIRS}
    assert all(not p for p in lie(shift).values())
    invariance=coefficient_polynomials(lie(shifted).values())
    assert affine_system(invariance,list(range(5,16)))==(8,8)
    sub={5+PARAMETERS.index(k):parse(v) for k,v in mixed['mixed_invariance']['substitution'].items()}
    P3={ij:substitute(p,sub) for ij,p in shifted.items()}
    assert P3==decode(mixed['brackets'])
    assert all(not p for p in lie(P3).values())
    assert source_direction_rank(P3,[6,12,15])==3
    raw3=coefficient_polynomials(jacobi(P3).values())
    conditions=[parse('3*l1+5*I'),parse('l1*(3*m4+1)')]
    assert polynomial_rank(raw3)==polynomial_rank(conditions)==polynomial_rank(raw3+conditions)==2
    final_sub={6:parse('-5*I/3'),12:parse('-1/3')}
    P1={ij:substitute(p,final_sub) for ij,p in P3.items()}
    assert all(not p for p in jacobi(P1).values())
    pfs=pfaffians(P1)
    actual=scale(add(*(multiply(variable(i),pfs[i]) for i in range(5))),F(1,5))
    expected_D=parse('5*(x0+I*x1)*(x0**2+x1**2+x2**2)*(x0**2+x1**2+x2**2+x3*x4)/6')
    assert actual==expected_D==parse(mixed['actual_divisor'])
    assert all(derivative(actual,i)==pfs[i] for i in range(5))
    checker=Path(__file__).read_text()
    result={'task_id':'E_reduced_nonnormal','all_passed':True,'field':'QQ(i)',
        'method':'Independent standard-library Gaussian Fraction sparse arithmetic; no producer or SymPy import.',
        'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        'source_sha256':hashlib.sha256(checker.encode()).hexdigest(),'independent_checker_source':checker,
        'invariant_coefficients':46,'full_linear_matrix_shape':[len(matrix),46],
        'full_linear_rank':35,'full_augmented_rank':35,'displayed_source_rank':11,
        'semisimple_Jacobi_linear_in_ell_and_affine_in_A':True,'all_ell_orbit_checks':orbit_checks,
        'semisimple_rank_four_support_empty':True,'canonical_mixed_shift_exact':True,
        'mixed_invariance_rank':8,'mixed_augmented_rank':8,'mixed_full_source_dimension':3,
        'mixed_Jacobi_coefficient_rowspace':['3*l1+5*I','l1*(3*m4+1)'],
        'mixed_full_Poisson_source':{'l1':'-5*I/3','m4':'-1/3','m7':'arbitrary'},
        'mixed_actual_divisor':'5*(x0+I*x1)*(x0^2+x1^2+x2^2)*(x0^2+x1^2+x2^2+x3*x4)/6',
        'all_five_Pfaffian_minors_equal_divisor_gradient':True,
        'mixed_irreducible_cubic_factor_impossible':True,
        'scope':'Excludes the two rank2-semisimple normal-Hamiltonian types from a smooth quadric times an irreducible cubic. Does not identify the mixed H+Q3+Q5 family component.',
        'global_exhaustion_proved':False,'elapsed_seconds':time.monotonic()-begun}
    with output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:result[k] for k in ['all_passed','full_linear_rank','mixed_invariance_rank','semisimple_rank_four_support_empty','mixed_irreducible_cubic_factor_impossible','elapsed_seconds']},indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output)



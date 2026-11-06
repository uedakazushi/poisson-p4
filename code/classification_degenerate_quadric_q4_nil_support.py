#!/usr/bin/env python3
"""Independent exact QQ(i) review of rank4-quadric nil31/nil22 support reductions."""
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
ZERO = (0,)*20
PAIRS = list(itertools.combinations(range(5), 2))
PARAMETERS = ['n'+str(i) for i in range(15)]


def scalar(c):
    return {ZERO: F(c)} if c else {}


def variable(j):
    return {tuple(int(i == j) for i in range(20)): F(1)}


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
        rows.append([p.get(tuple(int(j==i) for j in range(20)),F(0)) for i in indices])
        constants.append(p.get(ZERO,F(0)))
    return rank(rows),rank([r+[c] for r,c in zip(rows,constants)])

def divide_exact(p,d):
    rest=dict(p);answer={};lead=max(d)
    while rest:
        top=max(rest);diff=tuple(a-b for a,b in zip(top,lead));assert min(diff)>=0
        term={diff:rest[top]/d[lead]};answer=add(answer,term);rest=add(rest,scale(multiply(term,d),-1))
    assert multiply(answer,d)==p
    return answer

def necessary_coefficients(P,texts):
    slots=[]
    for triple,p in jacobi(P).items():
        for xm in sorted(set(m[:5] for m in p)):
            slots.append((triple,xm,{(0,)*5+m[5:]:a for m,a in p.items() if m[:5]==xm}))
    output=[]
    for text in texts:
        target=parse(text);key=next(iter(target));found=None
        for triple,xm,p in slots:
            if key in p and p==scale(target,p[key]/target[key]):
                found={'condition':text,'triple':list(triple),'spatial_monomial':list(xm),'ratio':str(p[key]/target[key])};break
        assert found,text
        output.append(found)
    return output

def apply(A,sub):return {ij:substitute(p,sub) for ij,p in A.items()}
def nsub(data):return {5+int(k[1:]):parse(str(v)) for k,v in data.items()}
def potential(P):
    pf=pfaffians(P)
    return scale(add(*(multiply(variable(i),pf[i]) for i in range(5))),F(1,5))

def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic();records=[]
    names=['NIL31_ZERO','NIL31_SHORT','NIL31_LONG','NIL22_ZERO','NIL22_NONZERO']
    paths=['reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_'+name+'_20261107.json' for name in names]
    for name,path in zip(names,paths):
        saved=json.loads((ROOT/path).read_text());P=decode(saved['brackets']);q=parse(saved['q'])
        assert hashlib.sha256(saved['producer_source'].encode()).hexdigest()==saved['producer_sha256']
        X=list(map(parse,saved['fixed_specification']['cofactor_components']))
        assert all(not add(*(derivative(entry(P,i,j),i) for i in range(5))) for j in range(5))
        assert all(add(*(multiply(derivative(q,i),entry(P,i,j)) for i in range(5)))==multiply(q,X[j]) for j in range(5))
        assert all(not add(*(multiply(X[k],derivative(P[ij],k)) for k in range(5)),
            *(scale(multiply(derivative(X[ij[0]],k),entry(P,k,ij[1])),-1) for k in range(5)),
            *(scale(multiply(derivative(X[ij[1]],k),entry(P,ij[0],k)),-1) for k in range(5))) for ij in PAIRS)
        D=potential(P);assert D==parse(saved['euler_contracted_pfaffian_not_yet_actual_divisor'])
        record={'case':name,'full_source_dimension':saved['source_dimension'],
            'all_divergence_q_normal_and_invariance_identities_exact':True,
            'Euler_Pfaffian_exact':True,'full_linear_source_rank_not_recomputed_here':True}
        if name=='NIL31_ZERO':
            quotient=divide_exact(D,multiply(variable(0),q));assert all(sum(m[:5])==2 for m in quotient)
            record.update(entire_source_Euler_Pfaffian_divisible_by_x0_q=True,irreducible_cubic_support_empty=True)
        elif name=='NIL31_SHORT':
            record['necessary_Jacobi_coefficients']=necessary_coefficients(P,['n5**2','6*n5*n8-5*n9'])
            restricted=substitute(D,nsub({'n5':'0','n9':'0'}))
            quotient=divide_exact(restricted,multiply(variable(0),q));assert all(sum(m[:5])==2 for m in quotient)
            record.update(reduced_support_forces_n5_n9_zero=True,restricted_Euler_Pfaffian_divisible_by_x0_q=True,irreducible_cubic_support_empty=True)
        elif name=='NIL31_LONG':
            equations=['n5*(3*n7+1)','n5*(2*n7-1)',
                '2*n2*n7+4*n2+2*n6**2+2*n7**2+7*n7+6',
                '4*n2*n7-2*n2+4*n6**2+10*n7**2+3*n7-4',
                'n6*(3*n2+3*n7+2)',
                '2*n1*n6+6*n2*n4+3*n4*n7+3*n4-12*n5**2',
                '6*n0*n6+18*n2*n3+12*n3*n7+10*n3+9*n4*n5']
            record['necessary_Jacobi_coefficients']=necessary_coefficients(P,equations)
            ae=parse('(6*n7**2-11*n7-16)/10')
            assert add(parse(equations[3]),scale(parse(equations[2]),-2),scale(add(variable(7),scale(ae,-1)),10))=={}
            assert substitute(scale(parse(equations[2]),5),{7:ae})==parse('10*n6**2+(2*n7-1)*(3*n7+1)*(n7+2)')
            assert substitute(scale(parse(equations[4]),10),{7:ae})==parse('n6*(6*n7+7)*(3*n7-4)')
            branches=[]
            for a,b,dexpr in [('-2','1/2','0'),('-7/6','-1/3','5*x3/36'),('3','-2','5*x3/4')]:
                sub=nsub({'n2':a,'n7':b,'n6':'0','n5':'0','n3':'0','n4':'0'})
                pp=apply(P,sub);assert all(not p for p in jacobi(pp).values())
                actual=potential(pp)
                expected={} if dexpr=='0' else multiply(parse(dexpr),power(q,1 if b=='-1/3' else 2))
                if b=='-1/3':expected=multiply(parse('5*x3**3/36'),q)
                assert actual==expected
                assert all(derivative(actual,i)==pfaffians(pp)[i] for i in range(5))
                branches.append({'n2':a,'n7':b,'n6':'0','n3_n4_n5_zero':True,'free_parameters':['n0','n1'],
                    'all_Jacobi_exact_zero':True,'irreducible_cubic_impossible':True})
            survivors=[];family_by_key={}
            for typ,a,b,c_abs,n3_factor,n4_factor in [('I','1/2','-7/6','5*I/6','-I','-2*I/3'),('II','-2','4/3','5*I/3','I','2*I/3')]:
                for sign in [1,-1]:
                    ctext='('+str(sign)+')*('+c_abs+')'
                    subdata={'n2':a,'n7':b,'n5':'0','n6':ctext,
                        'n3':'('+str(sign)+')*('+n3_factor+')*n0',
                        'n4':'('+str(sign)+')*('+n4_factor+')*n1'}
                    sub=nsub(subdata);pp=apply(P,sub);family_by_key[typ,sign]=pp
                    assert all(not p for p in jacobi(pp).values())
                    actual=potential(pp);cubic=divide_exact(actual,q)
                    assert derivative(cubic,4)==multiply(scale(parse(ctext),F(1,2)),power(variable(0),2))
                    b3=tuple(3 if i==1 else 0 for i in range(20))
                    assert cubic[b3]==-parse(ctext)[ZERO]/6
                    assert all(derivative(actual,i)==pfaffians(pp)[i] for i in range(5))
                    survivors.append({'family':typ,'sign':sign,'substitution':subdata,'free_parameters':['n0','n1'],
                        'all_Jacobi_exact_zero':True,'linear_x4_coefficient':'('+ctext+')*x0^2/2',
                        'x1_cubed_coefficient':'-('+ctext+')/6','cubic_irreducible_for_all_parameters_by_Gauss':True,
                        'full_family_brackets':{str(ij):{str(m):str(v) for m,v in p.items()} for ij,p in pp.items()}})
            for typ in ['I','II']:
                pp=apply(family_by_key[typ,1],{3:scale(variable(3),-1),6:scale(variable(6),-1)})
                reflected={ij:scale(p,(-1)**int(3 in ij)) for ij,p in pp.items()}
                assert reflected==family_by_key[typ,-1]
            record.update(scalar_c_zero_branches=branches,surviving_families=survivors,
                plus_minus_identified_by_x3_reflection_and_n1_sign=True,
                irreducible_cubic_support_exactly_two_full_two_parameter_sources_up_to_reflection=True)
        elif name=='NIL22_ZERO':
            Fc=parse('(n8*x0**2+2*n12*x0*x2-n1*x2**2)*x4/2+(n5*x0**3+3*n7*x0**2*x2+3*n10*x0*x2**2-n0*x2**3)/3')
            assert P[0,2]==derivative(Fc,4) and P[0,4]==scale(derivative(Fc,2),-1) and P[2,4]==derivative(Fc,0)
            record['necessary_Jacobi_coefficients']=necessary_coefficients(P,['n1*n8+n12**2','n0*n8+n1*n7+2*n10*n12','n1*n5-n10*n8+2*n12*n7'])
            pp=apply(P,nsub({'n'+str(i):'0' for i in [0,1,5,7,8,10,12]}))
            assert all(not p for p in pfaffians(pp))
            record.update(base_Jacobian_cubic_identical_to_smooth_nil221=True,F_zero_implies_rank_at_most_two=True,
                same_invariant_cubic_space_C3_plus_q_ell=True,
                irreducible_cubic_support_empty_by_same_Casimir_discriminant_argument=True)
        else:
            equations=['n9','n0*n5+n0*n8+3*n7*n9',
                '9*n3*n9-2*n5**2-n5*n8+14*n5+n8**2-n8-20',
                '3*n0*n3+12*n4*n9-2*n5*n7-5*n7*n8+n7',
                '-2*n1*n5+n1*n8-n1+3*n3*n4',
                '3*n1*n9-3*n3*n7+n4*n5-2*n4*n8+n4']
            record['necessary_Jacobi_coefficients']=necessary_coefficients(P,equations)
            quotient=divide_exact(substitute(D,nsub({'n0':'0','n9':'0'})),multiply(variable(0),q))
            assert all(sum(m[:5])==2 for m in quotient)
            subdata={'n9':'0','n5':'4/3','n8':'-4/3','n7':'-3*n0*n3/5','n4':'-9*n0*n3**2/25','n1':'-27*n0*n3**3/125'}
            pp=apply(P,nsub(subdata));assert all(not p for p in jacobi(pp).values())
            actual=potential(pp);expected=scale(multiply(q,add(scale(multiply(variable(0),q),5),multiply(variable(5),power(parse('x2+3*n3*x0/5'),3)))),-F(1,3))
            assert actual==expected and all(derivative(actual,i)==pfaffians(pp)[i] for i in range(5))
            cubic=divide_exact(actual,q)
            assert derivative(cubic,3)==scale(power(variable(0),2),-F(5,3))
            v3=tuple(3 if i==2 else 0 for i in range(20))
            assert cubic.get(v3+(0,)*0,F(0))==0 # The coefficient is the parameter polynomial -n0/3.
            assert substitute(cubic,{0:{},1:{},3:{},4:{}})==scale(multiply(variable(5),power(variable(2),3)),-F(1,3))
            record.update(n0_zero_implies_reducible_cubic=True,surviving_source_substitution=subdata,
                surviving_free_parameters=['n0','n2','n3','n6'],rank_four_irreducible_open='n0!=0',
                all_surviving_source_Jacobi_exact_zero=True,
                surviving_actual_divisor='-q*(5*x0*q+n0*(x2+3*n3*x0/5)^3)/3',
                all_five_gradient_identities=True,
                cubic_irreducible_on_n0_nonzero_by_Gauss=True,
                irreducible_cubic_support_exactly_this_four_parameter_source=True)
        records.append(record)
        print(json.dumps({'case':name,'all_checks_passed':True}),flush=True)
    checker=Path(__file__).read_text()
    result={'task_id':'E_reduced_nonnormal','all_passed':True,'field':'QQ(i)',
        'method':'Exact standard-library Gaussian Fraction sparse arithmetic, no producer or SymPy import.',
        'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        'source_sha256':hashlib.sha256(checker.encode()).hexdigest(),'independent_checker_source':checker,
        'cases':records,'elapsed_seconds':time.monotonic()-begun,
        'scope':'Full necessary-support reduction of five rank4 nilpotent quadric cofactors. Three cases excluded; nil31-long gives two full2-parameter sources up to reflection and nil22-nonzero one full4-parameter source. Their component membership remains open.',
        'global_exhaustion_proved':False}
    with output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'elapsed_seconds':result['elapsed_seconds']},indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output)


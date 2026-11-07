#!/usr/bin/env python3
"""Independent all-ninety-coefficient affine cofactor source verifier.

Standard-library Fraction Gaussian elimination, integer exponent bookkeeping,
and an AST polynomial parser.  No producer code or CAS is imported.
"""
import argparse
import ast
from fractions import Fraction as Q
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time


PAIRS = list(it.combinations(range(4),2))
COORDS = ('x','y','z','w')
MONS = []
for k,l in it.combinations_with_replacement(range(5),2):
    e = [0]*4
    if k:
        e[k-1] += 1
    if l:
        e[l-1] += 1
    MONS.append(tuple(e))
LINEARS = MONS[:5]


def linear_system(A,c):
    """Full 20 divergence and 90 Lie-derivative rows, including zero rows."""
    rows = [{} for _ in range(110)]
    rhs = [0]*110

    def insert(row,col,value):
        rows[row][col] = rows[row].get(col,0)+value
        if not rows[row][col]:
            del rows[row][col]

    for i in range(4):
        rhs[5*i] = c[i]
        for j in range(4):
            rhs[5*i+j+1] = A[i][j]
    for pair,(a,b) in enumerate(PAIRS):
        for m,e in enumerate(MONS):
            col = pair*15+m
            for coord,output,sign in [(a,b,1),(b,a,-1)]:
                if e[coord]:
                    d = list(e)
                    d[coord] -= 1
                    insert(5*output+LINEARS.index(tuple(d)),col,sign*e[coord])
            for d in range(4):
                if not e[d]:
                    continue
                f = list(e)
                f[d] -= 1
                if c[d]:
                    insert(20+pair*15+MONS.index(tuple(f)),col,c[d]*e[d])
                for t in range(4):
                    if A[d][t]:
                        g = f.copy()
                        g[t] += 1
                        insert(20+pair*15+MONS.index(tuple(g)),col,A[d][t]*e[d])
            for out,(u,v) in enumerate(PAIRS):
                change = 0
                for k in range(4):
                    p_kv = int((k,v)==(a,b))-int((k,v)==(b,a))
                    p_uk = int((u,k)==(a,b))-int((u,k)==(b,a))
                    change -= A[u][k]*p_kv+A[v][k]*p_uk
                if change:
                    insert(20+out*15+m,col,change)
    return rows,rhs


def complete_solution(rows,rhs):
    a = [{k:Q(v) for k,v in row.items()} | ({90:Q(b)} if b else {})
         for row,b in zip(rows,rhs)]
    pivots = []
    for col in range(90):
        r = len(pivots)
        k = next((k for k in range(r,len(a)) if a[k].get(col)),None)
        if k is None:
            continue
        a[r],a[k] = a[k],a[r]
        div = a[r][col]
        a[r] = {j:v/div for j,v in a[r].items()}
        for k in range(r+1,len(a)):
            q = a[k].get(col)
            if q:
                for j,v in a[r].items():
                    value = a[k].get(j,Q(0))-q*v
                    if value:
                        a[k][j] = value
                    else:
                        a[k].pop(j,None)
        pivots.append(col)
    assert all(not row for row in a[len(pivots):]), 'inconsistent affine system'
    free = [j for j in range(90) if j not in pivots]
    n = len(free)
    solution = [[Q(0)]*(n+1) for _ in range(90)]
    for i,j in enumerate(free):
        solution[j][i+1] = Q(1)
    for r,j in reversed(list(enumerate(pivots))):
        solution[j][0] = a[r].get(90,Q(0))
        for k,value in a[r].items():
            if k!=j and k!=90:
                for t,s in enumerate(solution[k]):
                    solution[j][t] -= value*s
    # Substitute every reconstructed coefficient into all original rows.
    for row,b in zip(rows,rhs):
        actual = [sum(v*solution[j][t] for j,v in row.items()) for t in range(n+1)]
        assert actual==[Q(b)]+[Q(0)]*n
    return free,pivots,solution


class Ring:
    def __init__(self,n):
        self.names = tuple(f'n{i}' for i in range(n))+COORDS
        self.n = n
        self.zero = (0,)*len(self.names)

    def constant(self,c):
        return {self.zero:Q(c)} if c else {}

    def variable(self,name):
        e = list(self.zero)
        e[self.names.index(name)] = 1
        return {tuple(e):Q(1)}

    def add(self,*args):
        out = {}
        for a in args:
            for e,c in a.items():
                out[e] = out.get(e,Q(0))+c
        return {e:c for e,c in out.items() if c}

    def scale(self,a,c):
        return {e:v*c for e,v in a.items() if v*c}

    def mul(self,*args):
        out = self.constant(1)
        for a in args:
            tmp = {}
            for e,c in out.items():
                for f,d in a.items():
                    g = tuple(x+y for x,y in zip(e,f))
                    tmp[g] = tmp.get(g,Q(0))+c*d
            out = {e:c for e,c in tmp.items() if c}
        return out

    def derivative(self,a,i):
        pos = self.n+i
        out = {}
        for e,c in a.items():
            if e[pos]:
                f = list(e)
                f[pos] -= 1
                out[tuple(f)] = c*e[pos]
        return out

    def parse(self,expression):
        def walk(node):
            if isinstance(node,ast.Constant) and type(node.value) is int:
                return self.constant(node.value)
            if isinstance(node,ast.Name) and node.id in self.names:
                return self.variable(node.id)
            if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):
                return self.scale(walk(node.operand),-1)
            if isinstance(node,ast.BinOp):
                a,b = walk(node.left),walk(node.right)
                if isinstance(node.op,ast.Add):
                    return self.add(a,b)
                if isinstance(node.op,ast.Sub):
                    return self.add(a,self.scale(b,-1))
                if isinstance(node.op,ast.Mult):
                    return self.mul(a,b)
                if isinstance(node.op,ast.Div):
                    assert set(b)=={self.zero}
                    return self.scale(a,1/b[self.zero])
                if isinstance(node.op,ast.Pow):
                    assert set(b)=={self.zero} or not b
                    exponent = b.get(self.zero,Q(0))
                    assert exponent.denominator==1 and exponent>=0
                    out = self.constant(1)
                    for _ in range(int(exponent)):
                        out = self.mul(out,a)
                    return out
            raise ValueError(ast.dump(node))
        return walk(ast.parse(expression,mode='eval').body)

    def coefficient(self,values):
        out = self.constant(values[0])
        for t,c in enumerate(values[1:]):
            out = self.add(out,self.scale(self.variable(f'n{t}'),c))
        return out

    def encode(self,a):
        return [[list(e),str(c)] for e,c in sorted(a.items())]


def reconstruct_normal_form(lam,variant):
    if isinstance(variant['jordan_blocks'],list):
        assert lam==[0]*4
        variant = dict(jordan_blocks={'0':variant['jordan_blocks']},
                       translation_top_length=variant['translation_top_length'])
    A = [[0]*4 for _ in range(4)]
    c = [0]*4
    index = 0
    done = False
    for value in sorted(set(lam)):
        for length in variant['jordan_blocks'][str(value)]:
            for t in range(length):
                A[index+t][index+t] = value
                if t:
                    A[index+t-1][index+t] = 1
            if value==0 and length==variant['translation_top_length'] and not done:
                c[index+length-1] = 1
                done = True
            index += length
    assert index==4 and [A[i][i] for i in range(4)]==lam
    assert done==(variant['translation_top_length']!=0)
    return A,c


def verify_source(row,index):
    eigenvalues = row.get('eigenvalues',[0]*4)
    A,c = reconstruct_normal_form(eigenvalues,row['variant'])
    assert A==row['field_linear_matrix'] and c==row['field_constant']
    rows,rhs = linear_system(A,c)
    free,pivots,values = complete_solution(rows,rhs)
    assert free==row['free_coefficient_indices']
    assert len(free)==row['free_parameter_count'] and len(pivots)==row['rank']
    R = Ring(len(free))
    source = [R.coefficient(v) for v in values]
    assert source==[R.parse(v) for v in row['all_ninety_coefficients']]
    entries = []
    for k in range(6):
        entry = {}
        for t,m in enumerate(MONS):
            mon = {(0,)*len(free)+m:Q(1)}
            entry = R.add(entry,R.mul(source[k*15+t],mon))
        entries.append(entry)
    assert entries==[R.parse(v) for v in row['brackets']]
    P = [[{} for j in range(4)] for i in range(4)]
    for (i,j),f in zip(PAIRS,entries):
        P[i][j],P[j][i] = f,R.scale(f,-1)
    jac = [R.add(*(R.mul(P[a][t],R.derivative(P[b][d],t))
                   for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4)))
           for i,j,k in it.combinations(range(4),3)]
    assert jac==[R.parse(v) for v in row['all_four_jacobi']]
    coefficients = set()
    for j in jac:
        grouped = {}
        for e,cvalue in j.items():
            grouped.setdefault(e[-4:],{})[e[:-4]+(0,)*4] = cvalue
        coefficients.update(tuple(sorted(v.items())) for v in grouped.values())
    assert coefficients=={tuple(sorted(R.parse(v).items())) for v in row['distinct_nonzero_jacobi_coefficients']}
    pf = R.add(R.mul(entries[0],entries[5]),R.scale(R.mul(entries[1],entries[4]),-1),R.mul(entries[2],entries[3]))
    assert pf==R.parse(row['actual_affine_pfaffian'])
    if row['factorized_pfaffian'] is not None:
        assert pf==R.parse(row['factorized_pfaffian'])
    top = {e:cvalue for e,cvalue in pf.items() if sum(e[-4:])==4}
    assert top==R.parse(row['highest_quartic_pfaffian'])
    divisors = [COORDS[i] for i in range(4) if pf and all(e[len(free)+i]>0 for e in pf)]
    matrix_bytes = json.dumps([[[k,v] for k,v in sorted(r.items())] for r in rows],separators=(',',':')).encode()
    return dict(type=index,eigenvalues=eigenvalues,variant=row['variant'],
                rank=len(pivots),free_indices=free,all_90_coefficients_reconstructed=True,
                all_110_original_linear_rows_satisfied=True,pivot_columns=pivots,
                original_integer_matrix_sha256=hashlib.sha256(matrix_bytes).hexdigest(),
                all_four_full_Jacobi_checked=True,all_distinct_Jacobi_coefficients_checked=True,
                actual_and_highest_pfaffians_checked=True,
                optional_factorized_pfaffian_checked=row['factorized_pfaffian'] is not None,
                unconditional_coordinate_divisors=divisors,pfaffian_identically_zero=not pf,
                pfaffian_term_count=len(pf),jacobi_term_counts=list(map(len,jac)))


def main():
    resource.setrlimit(resource.RLIMIT_CPU,(170,175))
    resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    started = time.monotonic()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--producer-report',type=Path,required=True)
    ap.add_argument('--entrance-report',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--nil-sources',action='store_true',help='All ten pure nilpotent affine cofactor types')
    args = ap.parse_args()
    raw = args.producer_report.read_bytes()
    report = json.loads(raw)
    assert hashlib.sha256(report['source_bytes'].encode()).hexdigest()==report['source_sha256']
    entrance_raw = args.entrance_report.read_bytes()
    entrance = json.loads(entrance_raw)
    if args.nil_sources:
        expected_labels = entrance['resonance_table']['pure_nilpotent_affine_types']
        rows = report['full_nilpotent_sources']
        assert len(rows)==10 and expected_labels==[r['variant'] for r in rows]
    else:
        expected_labels = [(r['eigenvalues'],v) for r in entrance['resonance_table']['semisimple_spectra']
                           for v in r['affine_jordan_variants']]
        rows = report['full_semisimple_sources']
        assert expected_labels==[(r['eigenvalues'],r['variant']) for r in rows]
    out = []
    for i,row in enumerate(rows,1):
        result = verify_source(row,i)
        out.append(result)
        print(json.dumps(dict(type=i,rank=result['rank'],free=len(result['free_indices']),
                              unconditional_divisors=result['unconditional_coordinate_divisors'])),flush=True)
    unconditional = [r['type'] for r in out if r['unconditional_coordinate_divisors'] or r['pfaffian_identically_zero']]
    result = dict(task='E_reduced_nonnormal',no_producer_import=True,no_CAS_import=True,
                  input_sha256={str(args.producer_report):hashlib.sha256(raw).hexdigest(),
                                str(args.entrance_report):hashlib.sha256(entrance_raw).hexdigest()},
                  reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  full_source_count=len(out),all_sources=out,
                  source_table='pure_nilpotent' if args.nil_sources else 'semisimple',
                  unconditionally_reducible_pfaffian_types=unconditional,
                  remaining_types=[r['type'] for r in out if r['type'] not in unconditional],
                  all_passed=True,
                  seconds=time.monotonic()-started,
                  maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  scope='Complete affine linear sources and their full Jacobi/Pfaffians verified. No residual Jacobi support or component coverage inferred for remaining types.')
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(all_passed=True,unconditional_types=unconditional,remaining_types=result['remaining_types'])))


if __name__=='__main__':
    main()

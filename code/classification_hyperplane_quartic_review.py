#!/usr/bin/env python3
"""Independent exact entrance check: Fraction RREF and sparse integer polys.

No producer/CAS imports.  This is a finite cofactor reduction, not coverage.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
import hashlib
import itertools as it
import json
import math
from pathlib import Path


def primitive(values):
    den = math.lcm(*(x.denominator for x in map(Q, values)))
    a = [int(Q(x)*den) for x in values]
    g = math.gcd(*a)
    assert g
    a = tuple(x//g for x in a)
    return a if next(x for x in a if x) > 0 else tuple(-x for x in a)


def wall4(values):
    """Orthogonal projection to trace-zero space, with no coordinates dropped."""
    return primitive([4*x-sum(values) for x in values])


def null_line(rows):
    """Unique projective nullspace, solved by full Fraction row reduction."""
    a = [list(map(Q, r)) for r in rows]
    pivot = []
    for col in range(4):
        k = next((k for k in range(len(pivot), len(a)) if a[k][col]), None)
        if k is None:
            continue
        r = len(pivot)
        a[r], a[k] = a[k], a[r]
        z = a[r][col]
        a[r] = [v/z for v in a[r]]
        for k in range(len(a)):
            if k != r and a[k][col]:
                z = a[k][col]
                a[k] = [v-z*w for v,w in zip(a[k], a[r])]
        pivot.append(col)
    if len(pivot) != 3:
        return None
    f = next(c for c in range(4) if c not in pivot)
    v = [Q(0)]*4
    v[f] = Q(1)
    for row, col in enumerate(pivot):
        v[col] = -a[row][f]
    return primitive(v)


def canonical(v):
    w = primitive(v)
    return min(tuple(sorted(w)), tuple(sorted(-x for x in w)))


def all_partitions(n):
    # Independent finite enumeration: enumerate nonincreasing tuples of
    # positive integers having sum n, rather than recursive producer logic.
    out = set()
    for length in range(1, n+1):
        for p in it.product(range(1,n+1), repeat=length):
            if sum(p)==n and list(p)==sorted(p, reverse=True):
                out.add(p)
    return sorted(out)


def reconstruct_table(report):
    walls = {}
    low = set()
    for i in range(1,5):
        for j in range(i+1,5):
            for k in range(5):
                for l in range(k,5):
                    v = [0]*4
                    for m,sign in [(i,1),(j,1),(k,-1),(l,-1)]:
                        if m:
                            v[m-1] += sign
                    if len(set(v))==1:
                        continue
                    key = wall4(v)
                    walls.setdefault(key, []).append([i,j,k,l])
                    if k==0:
                        low.add(key)
    claimed = report['resonance_table']
    # Compare complete labels, not just counts.
    assert len(claimed['walls'])==len(walls)
    for row in claimed['walls']:
        w = wall4(row['normal']+[0])
        assert sorted(row['monomial_resonances'])==sorted(walls[w])
        assert row['low_degree']==(w in low)
    spectra = set()
    for a in low:
        for b in walls:
            line = null_line([(1,1,1,1), a, b])
            if line is not None:
                spectra.add(canonical(line))
    assert {tuple(row['eigenvalues']) for row in claimed['semisimple_spectra']}==spectra
    variant_count = 0
    for row in claimed['semisimple_spectra']:
        lam = tuple(row['eigenvalues'])
        assert canonical(lam)==lam
        witness = row['intersecting_walls']
        assert len(witness)==2
        w4s = [wall4(w+[0]) for w in witness]
        assert w4s[0] in low and w4s[1] in walls
        assert null_line([(1,1,1,1),*w4s]) is not None
        assert all(sum(x*y for x,y in zip(w,lam))==0 for w in w4s)
        expected = set()
        eigenspaces = sorted(Counter(lam).items())
        for blocks in it.product(*(all_partitions(n) for _,n in eigenspaces)):
            for top in {0,*next((p for (v,_),p in zip(eigenspaces,blocks) if v==0), ())}:
                expected.add((tuple((v,p) for (v,_),p in zip(eigenspaces,blocks)),top))
        actual = {(tuple(sorted((int(k),tuple(v)) for k,v in t['jordan_blocks'].items())),
                   t['translation_top_length']) for t in row['affine_jordan_variants']}
        assert actual==expected and len(actual)==len(row['affine_jordan_variants'])
        variant_count += len(actual)
    expected_nil = {(p,t) for p in all_partitions(4) if p!=(1,1,1,1) for t in {0,*p}}
    actual_nil = {(tuple(r['jordan_blocks']),r['translation_top_length'])
                  for r in claimed['pure_nilpotent_affine_types']}
    assert expected_nil==actual_nil and len(actual_nil)==len(claimed['pure_nilpotent_affine_types'])
    low_orbits = {min(tuple(sorted(w)),tuple(sorted(-x for x in w))) for w in low}
    generic_orbits = {min(tuple(sorted(wall4(w))), tuple(sorted(-x for x in wall4(w))))
                      for w in ([1,1,0,0],[1,0,0,0],[1,1,-1,0])}
    assert low_orbits==generic_orbits
    return dict(walls=len(walls), low_walls=len(low), low_wall_S4_orbits=len(low_orbits),
                spectra=[list(v) for v in sorted(spectra)],
                affine_semisimple_types=variant_count, pure_nilpotent_affine_types=len(actual_nil),
                all_monomial_labels_and_canonical_witnesses_checked=True)


NAMES = ('C','A','B','D','E','F','p','q','r','x','y','z','w')
N = len(NAMES)
ZERO = (0,)*N


def const(n):
    return {ZERO:n} if n else {}


def variable(name):
    v = [0]*N
    v[NAMES.index(name)] = 1
    return {tuple(v):1}


def add(*args):
    out = {}
    for a in args:
        for m,c in a.items():
            out[m] = out.get(m,0)+c
    return {m:c for m,c in out.items() if c}


def scale(a,n):
    return {m:c*n for m,c in a.items() if c*n}


def mul(*args):
    out = const(1)
    for a in args:
        b = {}
        for m,c in out.items():
            for v,d in a.items():
                k = tuple(x+y for x,y in zip(m,v))
                b[k] = b.get(k,0)+c*d
        out = {m:c for m,c in b.items() if c}
    return out


def derivative(a,name):
    idx = NAMES.index(name)
    out = {}
    for m,c in a.items():
        if m[idx]:
            v = list(m)
            v[idx] -= 1
            out[tuple(v)] = c*m[idx]
    return out


def substitute(a,name,value):
    idx = NAMES.index(name)
    out = {}
    for m,c in a.items():
        power = m[idx]
        assert power>=0
        term = const(1)
        for _ in range(power):
            term = mul(term,value)
        v = list(m)
        v[idx] = 0
        out = add(out,mul({tuple(v):c},term))
    return out


def reciprocal_var(name):
    v = [0]*N
    v[NAMES.index(name)] = -1
    return {tuple(v):1}


def eq(a,b):
    assert not add(a,scale(b,-1))


def normalized(a):
    if not a:
        return ()
    g = math.gcd(*a.values())
    sign = 1 if a[min(a)]>0 else -1
    return tuple(sorted((m,c//g*sign) for m,c in a.items()))


def coefficients_in_coords(a):
    out = {}
    for mon,c in a.items():
        key = mon[-4:]
        u = mon[:-4]+(0,)*4
        out.setdefault(key,{})[u] = c
    return list(out.values())


def encode(a):
    return [[list(m),c] for m,c in sorted(a.items())]


def pair_zero_identities():
    C,A,B,D,E,F,p,q,r,x,y,z,w = map(variable,NAMES)
    v1,v2 = add(p,q),add(p,r)
    entries = [add(C,mul(A,x,y),mul(B,z,w)),mul(p,x,z),mul(q,x,w),
               mul(r,y,z),mul(scale(add(p,q,r),-1),y,w),add(D,mul(E,x,y),mul(F,z,w))]
    matrix = [[{} for j in range(4)] for i in range(4)]
    for (i,j),f in zip(it.combinations(range(4),2),entries):
        matrix[i][j],matrix[j][i] = f,scale(f,-1)
    jac = []
    for i,j,k in it.combinations(range(4),3):
        jac.append(add(*(mul(matrix[a][l],derivative(matrix[b][c],name))
                         for a,b,c in [(i,j,k),(j,k,i),(k,i,j)]
                         for l,name in enumerate(('x','y','z','w')))))
    conditions = [mul(B,E),add(mul(v2,C),mul(B,D)),add(mul(E,C),scale(mul(v1,D),-1)),
                  mul(B,add(v2,F)),mul(E,add(A,scale(v1,-1)))]
    actual = {normalized(c) for j in jac for c in coefficients_in_coords(j) if c}
    desired = {normalized(c) for c in conditions}
    assert actual==desired
    pf = add(mul(entries[0],entries[5]),scale(mul(entries[1],entries[4]),-1),mul(entries[2],entries[3]))
    eq(pf,add(mul(add(C,mul(A,x,y),mul(B,z,w)),add(D,mul(E,x,y),mul(F,z,w))),mul(v1,v2,x,y,z,w)))
    # B=E=0: either C=D=0, or v1*v2=0 is imposed by Jacobi.
    base = substitute(substitute(pf,'B',{}),'E',{})
    eq(substitute(substitute(base,'C',{}),'D',{}),mul(add(mul(A,F),mul(v1,v2)),x,y,z,w))
    eq(base,add(mul(add(C,mul(A,x,y)),add(D,mul(F,z,w))),mul(v1,v2,x,y,z,w)))
    # B-open: all five equations imply E=0,D=-v2*C/B,F=-v2,
    # and the remaining equation v1*v2*C=0.
    bopen = substitute(substitute(substitute(pf,'E',{}),'D',scale(mul(v2,C,reciprocal_var('B')),-1)),'F',scale(v2,-1))
    eq(substitute(bopen,'C',{}),mul(v2,z,w,add(mul(add(v1,scale(A,-1)),x,y),scale(mul(B,z,w),-1))))
    eq(substitute(bopen,'q',scale(p,-1)),scale(mul(v2,reciprocal_var('B'),add(C,mul(B,z,w)),add(C,mul(A,x,y),mul(B,z,w))),-1))
    eq(substitute(bopen,'r',scale(p,-1)),{})
    # E-open checked directly, rather than relying on a coordinate symmetry.
    eopen = substitute(substitute(substitute(pf,'B',{}),'C',mul(v1,D,reciprocal_var('E'))),'A',v1)
    eq(substitute(eopen,'D',{}),mul(v1,x,y,add(mul(E,x,y),mul(add(F,v2),z,w))))
    eq(substitute(eopen,'q',scale(p,-1)),{})
    eq(substitute(eopen,'r',scale(p,-1)),mul(v1,reciprocal_var('E'),add(D,mul(E,x,y)),add(D,mul(E,x,y),mul(F,z,w))))
    return dict(coefficient_ring_variables=NAMES, full_jacobi=[encode(j) for j in jac],
                necessary_and_sufficient_conditions=[encode(c) for c in conditions],
                jacobi_coefficients_exactly_signed_five_conditions=True,
                pfaffian=encode(pf), B_and_E_open_identities=True,
                all_closed_complements_retained=True)


def generic_allowed(report):
    # Weights are formal vectors in Z*a+Z*b, not sampled numbers.
    weights = {'pair_zero':[(0,0),(1,0),(-1,0),(0,1),(0,-1)],
               'zero_weight':[(0,0),(0,0),(1,0),(0,1),(-1,-1)],
               'triple_sum':[(0,0),(1,0),(0,1),(1,1),(-2,-2)]}
    answer = {}
    coords = [const(1)]+[variable(c) for c in ('x','y','z','w')]
    for typ,ww in weights.items():
        found = {}
        entries = []
        coeff = iter(NAMES[:9])
        for i,j in it.combinations(range(1,5),2):
            pairs = [[k,l] for k,l in it.combinations_with_replacement(range(5),2)
                     if all(ww[i][t]+ww[j][t]==ww[k][t]+ww[l][t] for t in range(2))]
            found[str(i)+str(j)] = pairs
            if typ!='pair_zero':
                entries.append(add(*(mul(variable(next(coeff)),coords[k],coords[l]) for k,l in pairs)))
        assert found==report['generic_walls']['allowed_terms'][typ]
        if typ!='pair_zero':
            pf = add(mul(entries[0],entries[5]),scale(mul(entries[1],entries[4]),-1),mul(entries[2],entries[3]))
            divisors = ('y','z','w') if typ=='zero_weight' else ('z','w')
            assert all(all(m[NAMES.index(c)]>=1 for c in divisors) for m in pf)
        answer[typ] = found
    return dict(all_formal_weight_terms=answer, zero_weight_pf_divisible_by='y*z*w',
                triple_sum_pf_divisible_by='z*w')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--producer-report',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args = ap.parse_args()
    raw = args.producer_report.read_bytes()
    report = json.loads(raw)
    # Embedded producer bytes are provenance only and are never evaluated.
    assert hashlib.sha256(report['source_bytes'].encode()).hexdigest()==report['source_sha256']
    result = dict(task='E_reduced_nonnormal',
                  producer_report=dict(path=str(args.producer_report),sha256=hashlib.sha256(raw).hexdigest()),
                  reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  no_producer_import=True,no_CAS_import=True,
                  table=reconstruct_table(report),generic_allowed=generic_allowed(report),
                  generic_pair_zero=pair_zero_identities(),all_passed=True,
                  scope='Complete finite entrance only. No residual-source Jacobi solution, componenthood or global coverage claim.')
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(all_passed=True,**{k:result['table'][k] for k in
                     ('walls','low_walls','low_wall_S4_orbits','affine_semisimple_types','pure_nilpotent_affine_types')})))


if __name__=='__main__':
    main()

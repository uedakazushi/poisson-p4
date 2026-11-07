#!/usr/bin/env python3
"""Exact modular resonance entrance for a simple hyperplane plus quartic.

Enumerate the complete weight arrangement, not a bounded eigenvalue search.
The output is a finite reduction, not a classification of its residual sources.
"""
import argparse
import hashlib
import itertools as it
import json
import math
import time
from collections import Counter
from pathlib import Path

import sympy as s


def primitive(v):
    d = math.gcd(*v)
    assert d
    w = tuple(x // d for x in v)
    return w if next(x for x in w if x) > 0 else tuple(-x for x in w)


def canonical_spectrum(v):
    d = math.gcd(*v)
    return min(tuple(sorted(x // d for x in v)), tuple(sorted(-x // d for x in v)))


def partitions(n, maximum=None):
    if not n:
        yield ()
    else:
        for k in range(min(n, maximum or n), 0, -1):
            for tail in partitions(n-k, k):
                yield (k,) + tail


def resonance_table():
    walls = {}
    low = set()
    for i, j in it.combinations(range(1, 5), 2):
        for k, l in it.combinations_with_replacement(range(5), 2):
            v = [0]*4
            for index, sign in ((i, 1), (j, 1), (k, -1), (l, -1)):
                if index:
                    v[index-1] += sign
            w = tuple(v[t]-v[3] for t in range(3))
            if not any(w):
                continue
            w = primitive(w)
            walls.setdefault(w, []).append([i, j, k, l])
            if k == 0:
                low.add(w)
    spectra = {}
    for a in sorted(low):
        for b in sorted(walls):
            cross = tuple(a[(t+1)%3]*b[(t+2)%3]-a[(t+2)%3]*b[(t+1)%3] for t in range(3))
            if not any(cross):
                continue
            lam = canonical_spectrum(cross + (-sum(cross),))
            spectra.setdefault(lam, [list(a), list(b)])
    representatives = []
    for lam in sorted(spectra):
        # Sorting/scalar normalization changes the coordinate labels. Select
        # witnesses in the canonical labels, rather than copying the raw pair.
        incident = [b for b in sorted(walls) if sum(b[t]*lam[t] for t in range(3)) == 0]
        witness = next([list(a), list(b)] for a in incident if a in low for b in incident
                       if any(a[(t+1)%3]*b[(t+2)%3]-a[(t+2)%3]*b[(t+1)%3] for t in range(3)))
        multiplicities = sorted(Counter(lam).items())
        variants = []
        for blocks in it.product(*(list(partitions(n)) for _, n in multiplicities)):
            jordan = {str(v): list(p) for (v, _), p in zip(multiplicities, blocks)}
            # Affine translation is reduced modulo im N. On the zero
            # eigenspace, the centralizer orbits are indexed by the largest
            # block length with nonzero coker component (plus zero).
            for top in [0]+sorted(set(jordan.get('0', []))):
                variants.append(dict(jordan_blocks=jordan, translation_top_length=top))
        representatives.append(dict(eigenvalues=list(lam), intersecting_walls=witness,
                                    affine_jordan_variants=variants))
    pure_nil = []
    for p in partitions(4):
        if p == (1, 1, 1, 1):
            continue
        for top in [0]+sorted(set(p)):
            pure_nil.append(dict(jordan_blocks=list(p), translation_top_length=top))
    return dict(walls=[dict(normal=list(k), monomial_resonances=v, low_degree=k in low)
                       for k, v in sorted(walls.items())],
                wall_count=len(walls), low_degree_wall_count=len(low),
                semisimple_spectra=representatives, semisimple_spectrum_count=len(representatives),
                affine_semisimple_type_count=sum(len(x['affine_jordan_variants']) for x in representatives),
                pure_nilpotent_affine_types=pure_nil)


def allowed_terms(weights):
    ans = {}
    for i, j in it.combinations(range(1, 5), 2):
        ans[f'{i}{j}'] = [[k, l] for k, l in it.combinations_with_replacement(range(5), 2)
                         if s.expand(weights[i]+weights[j]-weights[k]-weights[l]) == 0]
    return ans


def generic_walls():
    a, b = s.symbols('a b')
    allowed = {
        'pair_zero': allowed_terms([0, a, -a, b, -b]),
        'zero_weight': allowed_terms([0, 0, a, b, -a-b]),
        'triple_sum': allowed_terms([0, a, b, a+b, -2*a-2*b]),
    }
    expected = {
        'pair_zero': {'12': [[0,0],[1,2],[3,4]], '13': [[1,3]], '14': [[1,4]],
                      '23': [[2,3]], '24': [[2,4]], '34': [[0,0],[1,2],[3,4]]},
        'zero_weight': {'12': [[0,2],[1,2]], '13': [[0,3],[1,3]], '14': [[0,4],[1,4]],
                        '23': [[2,3]], '24': [[2,4]], '34': [[3,4]]},
        'triple_sum': {'12': [[0,3],[1,2]], '13': [[1,3]], '14': [[1,4]],
                       '23': [[2,3]], '24': [[2,4]], '34': [[3,4]]},
    }
    assert allowed == expected
    x, y, z, w = coords = s.symbols('x y z w')
    C,A,B,D,E,F,p,q,r = s.symbols('C A B D E F p q r')
    t = -p-q-r
    P = s.zeros(4)
    for (i,j), value in zip(it.combinations(range(4),2),
                           [C+A*x*y+B*z*w,p*x*z,q*x*w,r*y*z,t*y*w,D+E*x*y+F*z*w]):
        P[i,j], P[j,i] = value, -value
    jac = [s.expand(sum(P[i,l]*s.diff(P[j,k],coords[l])+
                        P[j,l]*s.diff(P[k,i],coords[l])+
                        P[k,l]*s.diff(P[i,j],coords[l]) for l in range(4)))
           for i,j,k in it.combinations(range(4),3)]
    v1,v2 = p+q,p+r
    conditions = [B*E,v2*C+B*D,E*C-v1*D,B*(v2+F),E*(A-v1)]
    coeffs = [c for f in jac for c in s.Poly(f,*coords).coeffs()]
    variables = (C,A,B,D,E,F,p,q,r)
    basis = s.groebner(conditions,*variables)
    reverse = s.groebner(coeffs,*variables)
    assert all(basis.reduce(f)[1]==0 for f in coeffs)
    assert all(reverse.reduce(f)[1]==0 for f in conditions)
    pf = s.expand(P[0,1]*P[2,3]-P[0,2]*P[1,3]+P[0,3]*P[1,2])
    claimed = (C+A*x*y+B*z*w)*(D+E*x*y+F*z*w)+v1*v2*x*y*z*w
    assert s.expand(pf-claimed)==0
    # Check each factorization argument in the written proof, preserving all
    # closed complements of the B or E localizations.
    generic_factor = s.factor(pf.subs({B:0,E:0,C:0,D:0}))
    assert s.expand(generic_factor - (A*F+v1*v2)*x*y*z*w)==0
    # B!=0: D=-v2*C/B and F=-v2. The remaining equation is v1*v2*C=0.
    bopen = s.factor(pf.subs({E:0,D:-v2*C/B,F:-v2}))
    assert s.expand(bopen.subs(C,0)-v2*z*w*((v1-A)*x*y-B*z*w))==0
    assert s.expand(bopen.subs(q,-p)-(-v2/B)*(C+B*z*w)*(C+A*x*y+B*z*w))==0
    assert s.expand(bopen.subs(r,-p))==0
    return dict(allowed_terms=allowed, pair_zero_jacobi=[str(f) for f in jac],
                pair_zero_conditions=[str(f) for f in conditions],
                jacobi_ideal_equal=True, pfaffian=str(pf), pfaffian_identity=True,
                b_open_factorization=str(bopen),
                factorization_complements=['B=E=0,C=D=0','B=E=0,C or D nonzero',
                                          'B!=0,C=0','B!=0,C!=0,v1=0 or v2=0',
                                          'E!=0: exchange (x,y) with (z,w)'])


def modular_normal_form(spectrum, variant):
    """Keep translations in each zero-eigenvalue coker-N orbit."""
    A = s.zeros(4)
    c = s.zeros(4, 1)
    index = 0
    translation_done = False
    for eigenvalue, blocks in sorted((int(v), p) for v, p in variant['jordan_blocks'].items()):
        for length in blocks:
            for j in range(length):
                A[index+j,index+j] = eigenvalue
                if j+1 < length:
                    A[index+j,index+j+1] = 1
            if eigenvalue == 0 and variant['translation_top_length'] == length and not translation_done:
                c[index+length-1] = 1
                translation_done = True
            index += length
    assert index == 4 and tuple(A[i,i] for i in range(4)) == tuple(spectrum)
    assert variant['translation_top_length'] == 0 or translation_done
    return A, c


def full_source(A, c, factorize=True):
    """Solve div(P)=Z and L_Z(P)=0 in all ninety affine coefficients."""
    begun = time.monotonic()
    coords = s.symbols('x y z w')
    homogeneous = (s.Integer(1),) + coords
    monomials = [homogeneous[k]*homogeneous[l] for k,l in it.combinations_with_replacement(range(5),2)]
    pairs = tuple(it.combinations(range(4),2))
    unknowns = s.symbols('a0:90')
    P = s.zeros(4)
    for k,(i,j) in enumerate(pairs):
        P[i,j] = sum(unknowns[15*k+t]*mon for t,mon in enumerate(monomials))
        P[j,i] = -P[i,j]
    Z = A*s.Matrix(coords)+c
    div = s.Matrix([sum(s.diff(P[i,j],coords[i]) for i in range(4)) for j in range(4)])
    lie = [s.expand(sum(Z[k]*s.diff(P[i,j],coords[k])-A[i,k]*P[k,j]-A[j,k]*P[i,k]
                         for k in range(4))) for i,j in pairs]
    equations = [v for f in list(div-Z)+lie for v in s.Poly(s.expand(f),*coords).coeffs() if v != 0]
    matrix, rhs = s.linear_eq_to_matrix(equations, unknowns)
    rr, pivots = matrix.row_join(rhs).rref()
    if 90 in pivots:
        return dict(linear_system_inconsistent=True, linear_system_shape=list(matrix.shape),
                    augmented_rank=len(pivots), seconds=time.monotonic()-begun)
    free_indices = [j for j in range(90) if j not in pivots]
    parameters = s.symbols(f'n0:{len(free_indices)}')
    free = dict(zip(free_indices,parameters))
    solution = [free.get(j,0) for j in range(90)]
    for row,j in enumerate(pivots):
        solution[j] = s.expand(rr[row,90]-sum(rr[row,k]*free[k] for k in free_indices))
    substitution = dict(zip(unknowns,solution))
    assert matrix*s.Matrix(solution) == rhs
    Pn = P.subs(substitution, simultaneous=True).applyfunc(s.expand)
    jacobi = [s.expand(sum(Pn[i,l]*s.diff(Pn[j,k],coords[l])+
                          Pn[j,l]*s.diff(Pn[k,i],coords[l])+
                          Pn[k,l]*s.diff(Pn[i,j],coords[l]) for l in range(4)))
              for i,j,k in it.combinations(range(4),3)]
    residual = list(dict.fromkeys(s.factor(value) for f in jacobi
                                 for value in s.Poly(f,*coords).coeffs() if value != 0))
    pf = s.expand(Pn[0,1]*Pn[2,3]-Pn[0,2]*Pn[1,3]+Pn[0,3]*Pn[1,2])
    top = sum(mon*coefficient for mon,coefficient in
              ((s.prod(v**e for v,e in zip(coords,exponent)),coefficient)
               for exponent,coefficient in s.Poly(pf,*coords).terms() if sum(exponent)==4))
    return dict(linear_system_shape=list(matrix.shape), rank=len(pivots),
                free_coefficient_indices=free_indices, free_parameter_count=len(parameters),
                parameters=list(map(str,parameters)), all_ninety_coefficients=list(map(str,solution)),
                brackets=[str(Pn[i,j]) for i,j in pairs], all_four_jacobi=list(map(str,jacobi)),
                distinct_nonzero_jacobi_coefficients=list(map(str,residual)),
                actual_affine_pfaffian=str(pf), factorized_pfaffian=str(s.factor(pf)) if factorize else None,
                highest_quartic_pfaffian=str(s.factor(top)) if factorize else str(s.expand(top)),
                full_linear_equations_satisfied=True, seconds=time.monotonic()-begun)


def semisimple_sources():
    result = []
    for row in resonance_table()['semisimple_spectra']:
        for variant in row['affine_jordan_variants']:
            A,c = modular_normal_form(row['eigenvalues'], variant)
            source = full_source(A,c)
            source.update(eigenvalues=row['eigenvalues'], variant=variant,
                          field_linear_matrix=[list(A.row(i)) for i in range(4)],
                          field_constant=list(c))
            # Cast SymPy integers at the JSON boundary, preserving exact strings
            # for all polynomial coefficients and equations.
            source['field_linear_matrix'] = [[int(v) for v in row] for row in source['field_linear_matrix']]
            source['field_constant'] = [int(v) for v in source['field_constant']]
            result.append(source)
            print(json.dumps(dict(type=len(result), eigenvalues=row['eigenvalues'],
                                  dimension=source.get('free_parameter_count'),
                                  seconds=round(source['seconds'],3))), flush=True)
    return result


def nilpotent_sources(checkpoint=None):
    result = []
    for variant in resonance_table()['pure_nilpotent_affine_types']:
        normal = dict(jordan_blocks={'0':variant['jordan_blocks']},
                      translation_top_length=variant['translation_top_length'])
        A,c = modular_normal_form([0,0,0,0],normal)
        # Large purely nilpotent coefficient spaces need only the raw exact
        # Pfaffian. Factoring their unconstrained Pfaffian adds no support
        # information before Jacobi and can exhaust the finite time budget.
        source = full_source(A,c,factorize=False)
        source.update(variant=variant,
                      field_linear_matrix=[[int(A[i,j]) for j in range(4)] for i in range(4)],
                      field_constant=[int(v) for v in c])
        result.append(source)
        if checkpoint:
            raw = Path(__file__).read_bytes()
            checkpoint.write_text(json.dumps(dict(task='E_reduced_nonnormal',
                source_sha256=hashlib.sha256(raw).hexdigest(),source_bytes=raw.decode(),
                full_nilpotent_sources=result,completed_types=len(result),
                all_passed=len(result)==10,scope='Full linear sources only; Jacobi support remains open.'),indent=2)+'\n')
        print(json.dumps(dict(type=len(result),variant=variant,
                              dimension=source.get('free_parameter_count'),
                              seconds=round(source['seconds'],3))),flush=True)
    return result


def semisimple_support(source_path):
    """Small exact factor certificates; branch completeness is in the note."""
    source = json.loads(source_path.read_text())
    x,y,z,w = coords = s.symbols('x y z w')
    n = s.symbols('n0:19')
    names = {str(v):v for v in coords+n}
    rows = source['full_semisimple_sources']
    pf = [s.sympify(row['actual_affine_pfaffian'],locals=names) for row in rows]
    universal = []
    for index,f in enumerate(pf,1):
        terms = s.Poly(f,*coords).terms()
        minimum = tuple(min(exponent[j] for exponent,_ in terms) for j in range(4))
        if any(minimum):
            factor = s.prod(v**e for v,e in zip(coords,minimum))
            quotient = s.cancel(f/factor)
            assert quotient.is_polynomial(*coords,*n)
            universal.append(dict(type=index,coordinate_factor=str(factor),quotient=str(quotient)))
    assert [v['type'] for v in universal] == [1,2,3,5,6,8,9,10,11,13,14,15,16,18,19]
    certificates = []
    def check(index,label,substitution,factor1,factor2,denominators=()):
        image = s.factor(pf[index-1].subs(substitution,simultaneous=True))
        assert s.cancel(image-factor1*factor2)==0,(index,label,s.factor(image-factor1*factor2))
        certificates.append(dict(type=index,branch=label,
                                 substitutions={str(k):str(v) for k,v in substitution.items()},
                                 factor1=str(factor1),factor2=str(factor2),
                                 inverted_elements=list(map(str,denominators)),identity=True))
    v = w*x
    # Type4: n0(n4-1)=0 and n3(n2-n5)=0; the final linear relation is
    # n0*n3=n1*(n2-1). These three equations cover all closed complements.
    f = s.factor(pf[3].subs(n[0],0))
    check(4,'n0=0',{n[0]:0},z,s.cancel(f/z))
    f = s.factor(pf[3].subs(n[3],0))
    check(4,'n3=0',{n[3]:0},y,s.cancel(f/y))
    check(4,'n0*n3!=0',{n[4]:1,n[5]:n[2],n[0]:n[1]*(n[2]-1)/n[3]},
          n[1]*y*y-2*n[3]*x*z,z-(n[2]-1)*w*y/n[3],(n[3],))
    # Type17: true reduced support, not a claim about the unreduced ideal.
    sub = {n[5]:s.Rational(1,2),n[3]:0}
    f = s.factor(pf[16].subs(sub,simultaneous=True))
    check(17,'(2n5-1)^2=0 implies n5=1/2; n3(3n5-1)=0 implies n3=0',sub,w*y,s.cancel(f/(w*y)))
    # Type21. The three equations below imply n4^3=0 on its n5-open branch.
    a,b,c = s.symbols('a b c')
    radical_basis = s.groebner([a*b,a*(3*c+1),b*(2*c+1)+a*a],a,b,c)
    assert radical_basis.reduce(a**3)[1]==0
    f = n[0]+n[1]*z+n[2]*z*z
    g = n[3]+n[4]*z+n[6]*z*z
    check(21,'n5=0,n9=0',{n[5]:0,n[9]:0},f+v,g)
    check(21,'n5=0,n9!=0',{n[5]:0,n[3]:0,n[4]:0},z*z,n[6]*f+(n[6]+n[9])*v,(n[9],))
    sub = {n[9]:s.Rational(1,4),n[2]:0,n[4]:0,n[1]:0,n[0]:0,n[3]:0}
    image = s.factor(pf[20].subs(sub,simultaneous=True))
    check(21,'n5!=0,n3=0; n4^3=0',sub,v,s.cancel(image/v),(n[5],))
    check(21,'n5*n3!=0; n4^3=0',
          {n[9]:s.Rational(1,4),n[2]:0,n[4]:0,n[1]:0,n[0]:n[3]/(2*n[5]),n[6]:-s.Rational(1,2)},
          n[3]/(2*n[5])+v/2,n[3]-z*z/2+n[5]*v,(n[5],n[3]))
    # Type22: n3=0, n4*n5=0, n1*(4*n4+1)=0,
    # (2*n4+1)*n0+2*n1*n2=0.
    t = z*z-2*y
    check(22,'n4!=0',{n[3]:0,n[5]:0},n[0]+n[1]*t+v,n[2]+n[4]*t,(n[4],))
    check(22,'n4=0',{n[3]:0,n[4]:0,n[1]:0,n[0]:0},v,n[2]-n[5])
    # Type24: all n5-open complements, then n5=0 with n6 nonzero/zero.
    sub = {n[10]:s.Rational(1,2),n[2]:0,n[6]:0,n[1]:0,n[0]:0}
    image = s.factor(pf[23].subs(sub,simultaneous=True))
    check(24,'n5!=0,n0=0; n2^2=0',sub,v,s.cancel(image/v),(n[5],))
    check(24,'n5*n0!=0; n2^2=0',
          {n[10]:s.Rational(1,2),n[2]:0,n[6]:0,n[1]:0,n[4]:-1},
          (v+2*n[0])/2,n[5]*v+n[3]-z,(n[5],n[0]))
    check(24,'n5=0,n6!=0',
          {n[5]:0,n[10]:0,n[9]:0,n[2]:0,n[1]:0,n[0]:0},v,g,(n[6],))
    check(24,'n5=n6=0,n4!=0',
          {n[5]:0,n[6]:0,n[9]:n[10]*n[3]/n[4]},
          n[3]+n[4]*z,f+(n[4]+n[10])*v/n[4],(n[4],))
    check(24,'n5=n6=n4=0,n3!=0',
          {n[5]:0,n[6]:0,n[4]:0,n[10]:0},
          s.Integer(1),n[3]*f+(n[3]+n[9])*v,(n[3],))
    # Here degree<=2, so homogenization to degree4 has an h^2 factor.
    certificates[-1]['highest_quartic_zero'] = True
    check(24,'n5=n6=n4=n3=0',
          {n[5]:0,n[6]:0,n[4]:0,n[3]:0},v,n[9]+n[10]*z)
    exclusions = sorted(set(v['type'] for v in universal+certificates))
    # The remaining regular paired spectra7 and12 admit true SNC curves.
    # Every Jacobi coefficient is homogeneous for each indicated rescaling,
    # so a whole source point, not only a tangent vector, is retained.
    curves = []
    for index,zero_index,weights in [
            (7,4,{0:1,7:1,2:2,3:1,1:3}),
            (7,2,{0:1,7:1,4:1,1:1,3:2}),
            (12,4,{0:1,6:1,2:1,3:1,1:2}),
            (12,2,{0:1,6:1,4:1,1:2,3:3})]:
        for expression in rows[index-1]['distinct_nonzero_jacobi_coefficients']:
            f = s.sympify(expression,locals=names).subs(n[zero_index],0)
            monomial_weights = {sum(exponent[j]*weights.get(j,0) for j in range(19))
                                for exponent,coefficient in s.Poly(f,*n).terms() if coefficient}
            assert len(monomial_weights)<=1
        limit = {n[zero_index]:0,**{n[j]:0 for j in weights}}
        kappa = 3*n[5]-n[6]-n[8] if index==7 else 2*n[5]-n[7]-n[8]
        assert s.expand(pf[index-1].subs(limit)-kappa*x*y*z*w)==0
        curves.append(dict(type=index,zero_parameter=str(n[zero_index]),
                           positive_rescaling_weights={str(n[j]):weight for j,weight in weights.items()},
                           all_jacobi_weight_homogeneous=True,SNC_pfaffian=str(kappa*x*y*z*w),
                           inverted_element=str(kappa),complement='kappa=0 retained below'))
    # Type7 kappa=0: put a=n0*n7 and E=a-n2*(n5-1/2).
    # This exact identity proves a^2 belongs to (Jacobi,kappa).
    a = n[0]*n[7]
    E = a-n[2]*(n[5]-s.Rational(1,2))
    kappa = 3*n[5]-n[6]-n[8]
    rhs = a*E+n[2]/2*(a*kappa+n[7]*n[0]*(n[6]-1)-n[0]*n[7]*(n[5]-n[8]))
    assert s.expand(a*a-rhs)==0
    kappa_zero_start = len(certificates)
    for branch,substitution in [
            ('n3!=0,n0=0',{n[4]:0,n[5]:0,n[8]:-n[6],n[1]:0,n[2]:0,n[0]:0}),
            ('n3!=0,n7=0',{n[4]:0,n[5]:0,n[8]:-n[6],n[1]:0,n[2]:0,n[7]:0}),
            ('n3=0,n4!=0,n0=0',{n[3]:0,n[5]:s.Rational(1,2),n[8]:s.Rational(3,2)-n[6],n[1]:0,n[2]:0,n[0]:0}),
            ('n3=0,n4!=0,n7=0',{n[3]:0,n[5]:s.Rational(1,2),n[8]:s.Rational(3,2)-n[6],n[1]:0,n[2]:0,n[7]:0}),
            ('n3=n4=0,n0=0',{n[3]:0,n[4]:0,n[8]:3*n[5]-n[6],n[0]:0}),
            ('n3=n4=0,n7=0',{n[3]:0,n[4]:0,n[8]:3*n[5]-n[6],n[7]:0})]:
        image = s.factor(pf[6].subs(substitution,simultaneous=True))
        factor = x if branch in ['n3!=0,n0=0','n3=0,n4!=0,n0=0'] else (
                 w if branch in ['n3!=0,n7=0','n3=0,n4!=0,n7=0'] else
                 z if branch=='n3=n4=0,n0=0' else y)
        check(7,'kappa=0,'+branch,substitution,factor,s.cancel(image/factor))
    # Type12 kappa=0: if n4 or n2 is nonzero, n5=1/2 and
    # n7+n8=1. Then n0*n6!=0 would force n7=1,n8=1/2,
    # a contradiction. The closed n2=n4=0 complement has K=0.
    for branch,substitution,factor in [
            ('n4!=0,n0=0',{n[5]:s.Rational(1,2),n[8]:1-n[7],n[1]:0,n[2]:0,n[3]:0,n[0]:0},x),
            ('n4!=0,n6=0',{n[5]:s.Rational(1,2),n[8]:1-n[7],n[1]:0,n[2]:0,n[3]:0,n[6]:0},w),
            ('n4=0,n2!=0,n0=0',{n[5]:s.Rational(1,2),n[8]:1-n[7],n[1]:0,n[4]:0,n[3]:0,n[0]:0},z),
            ('n4=0,n2!=0,n6=0',{n[5]:s.Rational(1,2),n[8]:1-n[7],n[1]:0,n[4]:0,n[3]:0,n[6]:0},y)]:
        image = s.factor(pf[11].subs(substitution,simultaneous=True))
        check(12,'kappa=0,'+branch,substitution,factor,s.cancel(image/factor))
    top12 = s.sympify(rows[11]['highest_quartic_pfaffian'],locals=names)
    assert s.expand(top12.subs({n[2]:0,n[4]:0,n[8]:2*n[5]-n[7]},simultaneous=True))==0
    example = {n[0]:1,n[2]:0,n[7]:0,n[1]:s.Rational(3,2),n[3]:1,n[4]:1,
               n[5]:1,n[6]:1,n[8]:s.Rational(1,2)}
    assert all(s.sympify(f,locals=names).subs(example)==0 for f in rows[6]['all_four_jacobi'])
    example_pf = s.factor(pf[6].subs(example))
    expected = s.Rational(3,2)*(v+1)*(v+y*z+1)-w*y**3
    assert s.expand(example_pf-expected)==0
    assert s.gcd(s.Rational(3,2)*y*(v+1),s.Rational(3,2)*(v+1)**2-w*y**3)==1
    return dict(input_report=str(source_path),input_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                unconditional_coordinate_exclusions=universal,factorization_certificates=certificates,
                excluded_types=exclusions,covered_by_MPS_types=[7,12],
                SNC_curves=curves,kappa_zero_factorizations=certificates[kappa_zero_start:],
                type12_kappa_zero_n2_n4_zero_highest_pfaffian_zero=True,
                remaining_types=[20,23],
                reduced_radical_identity={'generators':['a*b','a*(3*c+1)','b*(2*c+1)+a^2'],
                                          'verified_ideal_member':'a^3'},
                type7_actual_irreducible_quartic_example=dict(parameters={str(k):str(v) for k,v in example.items()},
                     affine_pfaffian=str(example_pf),all_jacobi_zero=True,
                     irreducibility_reason='Primitive polynomial linear in z; coefficients (3/2)y(wx+1) and (3/2)(wx+1)^2-wy^3 are coprime over C[x,y,w]. Homogenization has no h factor.'),
                branch_completeness_requires_written_jacobi_proof=True)


def nil_rank_one_quartic_geometry():
    """Weak divisor recognition; does not infer a full Poisson fiber orbit."""
    y,t=s.symbols('y t')
    c,u,v=s.symbols('c u v')
    f=c+u*y+v*y*y+3*y*t
    monomials=[y**i*t**j for i in range(5) for j in range(5-i)]
    unknowns=s.symbols('a0:15')
    F=sum(a*m for a,m in zip(unknowns,monomials))
    residual=s.Poly(s.expand(f*s.diff(F,t)-y*y*s.diff(F,y)),y,t)
    M,_=s.linear_eq_to_matrix(residual.coeffs(),unknowns)
    specialized=M.subs({c:1,u:0,v:0})
    _,pivots=specialized.rref()
    kernel=M.nullspace()
    I=y**3*t+c*y*y/2+u*y**3/3+v*y**4/4
    basis=s.Matrix([[s.Poly(g,y,t).coeff_monomial(m) for g in (s.Integer(1),I)] for m in monomials])
    assert M*basis==s.zeros(M.rows,2) and basis.rank()==2
    assert len(pivots)==13 and len(kernel)==2
    rows=specialized[:,pivots].T.rref()[1]
    minor=s.factor(M[list(rows),list(pivots)].det())
    # This minor has no exceptional parameter except the already inverted c.
    assert set(minor.free_symbols)<= {c} and minor!=0
    assert s.factor(minor).subs(c,1)!=0
    return dict(coefficient_ring='QQ[c,u,v,y,t], c!=0',
                full_PDE='(c+u*y+v*y^2+3*y*t)*F_t-y^2*F_y=0',
                all_degree_at_most_four_monomials=list(map(str,monomials)),
                full_matrix=[list(map(str,M.row(i))) for i in range(M.rows)],
                rank=13,complete_kernel_basis=['1',str(I)],
                pivot_rows=list(rows),pivot_columns=list(pivots),pivot_minor=str(minor),
                every_parameter_exception_retained=True,
                geometric_implication='For nil211,c=0, an irreducible quartic Pfaffian has nonzero restriction to y=0. Invariance and the modular identity force this PDE after affine linear changes; its divisor is projectively MINIMAL_M. Full Poisson source membership is not implied.',
                remaining_obligation='Recognize the whole Poisson fiber with this fixed divisor, or give an all-point deformation to a known closed component.')


def nil_rank_one_quartic_support():
    """Full irreducible MINIMAL_M source, with a fixed minimal modular field."""
    x,y,z,w=s.symbols('x y z w')
    a0,c,u,v,A=s.symbols('a0 c u v A')
    magnetic=s.symbols('b0:6')
    B=sum(a*m for a,m in zip(magnetic,[1,y,w,y*y,y*w,w*w]))
    f=c+u*y+v*y*y+3*y*w
    F=a0*c-y**3*w-c*y*y/2-u*y**3/3-v*y**4/4
    a=a0-a0*u*y/c-3*a0*y*w/c+A*y*y
    d=s.cancel((F-a*f)/(y*y)).expand()
    assert s.denom(d).free_symbols <= {c}
    b=s.expand(x*y+z*(-y-s.diff(a,y)-s.diff(d,w))+B)
    entries=[a,b,d,y*y,0,f]
    coords=(x,y,z,w)
    P=s.zeros(4)
    for (i,j),value in zip(it.combinations(range(4),2),entries):
        P[i,j],P[j,i]=value,-value
    div=[s.expand(sum(s.diff(P[i,j],coords[i]) for i in range(4))) for j in range(4)]
    assert div==[y,0,0,0]
    jac=[s.expand(sum(P[i,l]*s.diff(P[j,k],coords[l])+P[j,l]*s.diff(P[k,i],coords[l])+
                      P[k,l]*s.diff(P[i,j],coords[l]) for l in range(4)))
         for i,j,k in it.combinations(range(4),3)]
    assert all(s.cancel(f)==0 for f in jac)
    pf=s.expand(P[0,1]*P[2,3]-P[0,2]*P[1,3]+P[0,3]*P[1,2])
    assert s.cancel(pf-F)==0
    for f0 in entries:
        assert s.Poly(f0,*coords).total_degree()<=2
    # The old MINIMAL_M representative in these coordinates. Its known
    # uniqueness of the ambient component is an external repository input.
    old=s.symbols('X Y Z W')
    X,Y,Z,W=old
    oldentries=[1+3*X*X-6*X*Z,0,-2*X*X,3*X*Z-18*Z*Z,
                6*X*X-2*X*Y-3*X*W+30*Z*W,1+6*X*Z]
    oldP=s.zeros(4)
    for (i,j),value in zip(it.combinations(range(4),2),oldentries):
        oldP[i,j],oldP[j,i]=value,-value
    T=s.Matrix([[0,s.Rational(1,6),0,0],[1,0,0,0],[0,0,0,s.Rational(-1,2)],[0,0,1,0]])
    transformed=(T*oldP*T.T).subs(dict(zip(old,T.inv()*s.Matrix(coords))),simultaneous=True)
    point={a0:s.Rational(-1,6),c:s.Rational(1,2),u:0,v:0,A:s.Rational(-1,2),
           **dict(zip(magnetic,[0,0,0,s.Rational(-1,2),0,0]))}
    assert all(s.expand(transformed[i,j]-P[i,j].subs(point))==0 for i,j in it.combinations(range(4),2))
    return dict(localization=['a0!=0','c!=0'],complete_source_parameter_count=11,
                parameters=list(map(str,[a0,c,u,v,A,*magnetic])),
                brackets=list(map(str,entries)),actual_affine_pfaffian=str(F),
                all_four_jacobi_identically_zero=True,full_divergence=list(map(str,div)),
                all_affine_degrees_at_most_two=True,
                MINIMAL_M_coordinate_matrix=[[str(t) for t in T.row(i)] for i in range(4)],
                MINIMAL_M_source_point={str(k):str(v0) for k,v0 in point.items()},
                MINIMAL_M_all_six_brackets_exactly_matched=True,
                source_completeness='F irreducible forces F|y=0 nonzero. Modular Casimir equations force a|y=0 and f|y=0 to be nonzero constants, l=alpha*y,m=beta*y. Nonzero(alpha,beta) normalized to(1,0). The full degree<=4 PDE and modular identity force F above; Pf=af+y^2*d and divergence then give exactly these a,d,b with arbitrary B2(y,w).',
                classification_implication='This irreducible eleven-parameter source contains the known MINIMAL_M point, at which the full Poisson scheme has a unique component. Hence the whole source and its affine/PGL saturation lie in that closed component. Zero(alpha,beta) has no irreducible quartic; all normalizations have retained complements.',
                repository_proof_dependency='legacy/v07/minimal_thickness.tex Theorem minimal-full and unique minimal prime; apply only unique reduced component, not recomputed nil exponent.')


def nil_regular_support(source_path):
    """All scalar branches of regular nilpotent modular type4.

    The two remaining normalized sources are retained for recognition;
    their actual quartic is not mislabeled a new component.
    """
    rows=json.loads(source_path.read_text())['full_nilpotent_sources']
    n=s.symbols('n0:36');x,y,z,w=coords=s.symbols('x y z w')
    pairs=tuple(it.combinations(range(4),2))
    J=[[s.sympify(f) for f in row['distinct_nonzero_jacobi_coefficients']] for row in rows[:2]]
    pf=[s.sympify(row['actual_affine_pfaffian']) for row in rows[:2]]
    type2=[]
    for value,substitution in [
            (-s.Rational(1,9),{n[0]:0,n[1]:0,n[2]:0,n[4]:0,n[5]:0,n[6]:n[8]}),
            (-s.Rational(1,10),{n[0]:10*n[5]*n[8],n[1]:5*n[8]**2,n[2]:0,
                               n[3]:-5*n[5]-10*n[8]**2,n[4]:-n[8]/2,n[6]:3*n[8]/2})]:
        expressions=[s.expand(f.subs(n[9],value)) for f in J[1]]
        actual=s.groebner(expressions,*n[:9],domain=s.QQ)
        expected=[variable-expression for variable,expression in substitution.items()]
        basis=s.groebner(expected,*n[:9],domain=s.QQ)
        assert all(basis.reduce(f)[1]==0 for f in expressions)
        assert all(actual.reduce(f)[1]==0 for f in expected)
        image=s.factor(pf[1].subs({n[9]:value,**substitution},simultaneous=True))
        if value==-s.Rational(1,9):assert image==0
        else:
            square=-20*n[5]-10*n[8]*w*w+20*n[8]*z-2*w*y+2*x+z*z
            assert s.expand(image+square*square/80)==0
        type2.append(dict(n9=str(value),complete_parameter_solution={str(k):str(v)for k,v in substitution.items()},
                          ideals_equal_exactly=True,actual_pfaffian=str(image),irreducible_quartic_excluded=True))
    radical=[]
    for label,substitution,variables in [
            ('n8=-1/4',{n[8]:-s.Rational(1,4)},[v for v in n[:13]if v!=n[8]]),
            ('n12=5n8+1',{n[12]:5*n[8]+1},list(n[:12]))]:
        basis=s.groebner([s.expand(f.subs(substitution,simultaneous=True))for f in J[0]],*variables,domain=s.QQ)
        assert basis.reduce(n[11]**3)[1]==0
        radical.append(dict(branch=label,verified_ideal_member='n11^3',reduced_consequence='n11=0'))
    # On the second branch away from the first, the exact Jacobi equations
    # force n5=n2=0 after n11=0; highest Pfaffian is then identically zero.
    top=s.sympify(rows[0]['highest_quartic_pfaffian'])
    assert s.expand(top.subs({n[12]:5*n[8]+1,n[11]:0,n[2]:0,n[5]:0},simultaneous=True))==0
    # On n8=-1/4, n11=0, Jacobi gives n7=0. At n12=-1/4,
    # n10=n6=n1=n0=0 and Pf has a w^2 factor.
    overlap={n[8]:-s.Rational(1,4),n[12]:-s.Rational(1,4),n[11]:0,n[7]:0,
             n[10]:0,n[6]:0,n[1]:0,n[0]:0}
    overlap_pf=s.factor(pf[0].subs(overlap,simultaneous=True))
    assert s.cancel(overlap_pf/(w*w)).is_polynomial(*coords,*n)
    # All remaining n5 are removed by an actual linear centralizer of N.
    d=4*n[12]+1
    initial={n[8]:-s.Rational(1,4),n[11]:0,n[7]:0,n[9]:-3*n[5],n[2]:12*n[5]**2/d}
    P=s.zeros(4)
    for (i,j),text in zip(pairs,rows[0]['brackets']):
        value=s.sympify(text).subs(initial,simultaneous=True)
        P[i,j],P[j,i]=value,-value
    N=s.zeros(4)
    for i in range(3):N[i,i+1]=1
    transform=s.eye(4)-6*n[5]/d*N**2
    assert transform*N==N*transform and transform.det()==1
    changed=(transform*P*transform.T).subs(dict(zip(coords,transform.inv()*s.Matrix(coords))),simultaneous=True)
    changed=changed.applyfunc(s.cancel)
    hcoords=(s.Integer(1),)+coords
    mons=[hcoords[i]*hcoords[j] for i,j in it.combinations_with_replacement(range(5),2)]
    changed_coefficients=[s.Poly(changed[i,j],*coords).coeff_monomial(mon)for i,j in pairs for mon in mons]
    new_parameters={n[j]:s.factor(changed_coefficients[index])for j,index in enumerate(rows[0]['free_coefficient_indices'])}
    assert new_parameters[n[5]]==new_parameters[n[9]]==new_parameters[n[2]]==0
    assert new_parameters[n[8]]==-s.Rational(1,4) and new_parameters[n[12]]==n[12]
    for (i,j),text in zip(pairs,rows[0]['brackets']):
        assert s.cancel(changed[i,j]-s.sympify(text).subs(new_parameters,simultaneous=True))==0
    # Full normalized Jacobi is exactly this small four-equation system.
    normalized={n[8]:-s.Rational(1,4),n[11]:0,n[7]:0,n[9]:0,n[5]:0,n[2]:0}
    residual=[s.expand(f.subs(normalized,simultaneous=True))for f in J[0]]
    equations=[n[10]*(2*n[12]+1),n[6]*(1-2*n[12])-6*n[10]**2,
               n[1]*(3*n[12]+1)-3*n[10]*n[4],n[0]*(2*n[12]+1)-2*n[10]*n[3]]
    variables=[n[j]for j in [0,1,3,4,6,10,12]]
    actual=s.groebner(residual,*variables,domain=s.QQ);expected=s.groebner(equations,*variables,domain=s.QQ)
    assert all(expected.reduce(f)[1]==0 for f in residual)
    assert all(actual.reduce(f)[1]==0 for f in equations)
    I2=2*w*y-z*z;I3=3*w*w*x-3*w*y*z+z**3
    checks=[]
    def verify(label,substitution,claimed):
        image=s.factor(pf[0].subs({**normalized,**substitution},simultaneous=True))
        assert s.expand(image-claimed)==0,(label,s.factor(image-claimed))
        checks.append(dict(branch=label,substitution={str(k):str(v)for k,v in substitution.items()},
                           actual_pfaffian=str(image),claimed_pfaffian=str(claimed),identity=True))
    verify('n12 outside {-1/2,-1/3,1/2}',{n[10]:0,n[6]:0,n[1]:0,n[0]:0},-(4*n[12]+1)*I2**2/48)
    verify('n12=1/2',{n[12]:s.Rational(1,2),n[10]:0,n[1]:0,n[0]:0},-(I2+4*n[6])**2/16)
    verify('n12=-1/2,n10=0',{n[12]:-s.Rational(1,2),n[10]:0,n[6]:0,n[1]:0},I2**2/48-n[0]*w*w/2)
    verify('n12=-1/3',{n[12]:-s.Rational(1,3),n[10]:0,n[6]:0,n[0]:0},I2**2/144-n[1]*w**3/3)
    verify('n12=-1/2,n10!=0',{n[12]:-s.Rational(1,2),n[6]:3*n[10]**2,n[1]:-6*n[10]*n[4],n[3]:0},
           I2**2/48-n[10]*I3/3-3*n[10]**2*I2/2-9*n[10]**4-n[0]*w*w/2+2*n[10]*n[4]*w**3)
    return dict(input_report=str(source_path),input_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                type2_all_parameter_branches=type2,type1_reduced_radical_certificates=radical,
                type1_second_branch_complement='n8!=-1/4 gives n11=n5=n2=0 and highest Pfaffian zero',
                type1_overlap_pfaffian=str(overlap_pf),
                centralizer_normalization=dict(inverted_element=str(d),
                    complement='n12=-1/4 excluded separately by w^2 factor',
                    matrix=[[str(transform[i,j])for j in range(4)]for i in range(4)],
                    transformed_parameters={str(k):str(v)for k,v in new_parameters.items()},
                    all_six_brackets_exactly_matched=True),
                normalized_full_jacobi_equations=list(map(str,equations)),normalized_jacobi_ideals_equal=True,
                all_normalized_pfaffian_identities=checks,
                type1_two_remaining_sources=['n12=-1/3,n1!=0','n12=-1/2,n10!=0'],
                scope='Type2 excluded entirely; type1 reduced to two complete normalized sources; their membership remains unproved.')


def nil_regular_component(source_path, second=False, point=None, quadratic=False, twist=False):
    """Full150-variable one-normal support test, through at most cubic order."""
    import importlib.util
    from sympy.polys.matrices import DomainMatrix
    helper=Path('legacy/v03/code/poisson_p4_components_check.py')
    spec=importlib.util.spec_from_file_location('quartic_ambient',helper)
    ambient=importlib.util.module_from_spec(spec);spec.loader.exec_module(ambient)
    row=json.loads(source_path.read_text())['full_nilpotent_sources'][point['source_row'] if point else 0]
    n=s.symbols('n0:36');h,x,y,z,w=coords=s.symbols('h x y z w')
    substitution={v:0 for v in n}
    substitution.update({n[8]:-s.Rational(1,4),n[12]:-s.Rational(1,3),n[1]:1})
    if second:
        substitution.update({n[12]:-s.Rational(1,2),n[1]:0,n[10]:1,n[6]:3})
    if point:
        substitution={v:0 for v in n}
        substitution.update({s.Symbol(k):s.Rational(v)for k,v in point['substitutions'].items()})
    P=s.zeros(5)
    for (i,j),text in zip(it.combinations(range(1,5),2),row['brackets']):
        f=s.Poly(s.sympify(text).subs(substitution),x,y,z,w)
        P[i,j]=sum(c*h**(2-sum(e))*s.prod(v**a for v,a in zip((x,y,z,w),e))for e,c in f.terms())
        P[j,i]=-P[i,j]
    denominator=s.ilcm(*[c.q for i,j in ambient.PAIRS for c in s.Poly(P[i,j],*coords).coeffs()])
    P*=denominator
    def vector(A):
        return s.Matrix([s.Poly(A[i,j],*coords).coeff_monomial(mon)for (i,j),mon in ambient.BIVECTOR_BASIS])
    def unvector(v):
        A=s.zeros(5)
        for coefficient,((i,j),mon) in zip(v,ambient.BIVECTOR_BASIS):
            value=coefficient*s.prod(v**e for v,e in zip(coords,mon))
            A[i,j]+=value;A[j,i]-=value
        return A.applyfunc(s.expand)
    def jac(A):
        J={}
        for i,j,k in it.combinations(range(5),3):
            J[i,j,k]=s.expand(sum(A[i,l]*s.diff(A[j,k],coords[l])+A[j,l]*s.diff(A[k,i],coords[l])+
                                  A[k,l]*s.diff(A[i,j],coords[l])for l in range(5)))
        W={ij:s.expand(sum((-1)**k*coords[ij[k]]*J[ij[:k]+ij[k+1:]]for k in range(4)))
           for ij in it.combinations(range(5),4)}
        return s.Matrix([s.Poly(W[ij],*coords).coeff_monomial(mon)for ij,mon in ambient.FOUR_VECTOR_BASIS])
    def dm(A):return DomainMatrix.from_Matrix(A).convert_to(s.QQ)
    def minor(A):
        columns=list(dm(A).rref()[1]);rows=list(dm(A[:,columns].T).rref()[1]);det=dm(A[rows,columns]).det()
        assert det
        return dict(rank=len(columns),rows=rows,columns=columns,determinant=str(det))
    def encode(v):return {str(i):str(value)for i,value in enumerate(v)if value}
    def digest(A):return hashlib.sha256(json.dumps(A.tolist(),default=str).encode()).hexdigest()
    assert jac(P)==s.zeros(350,1)
    sparse={(i,j):{mon:int(c)for mon,c in s.Poly(P[i,j],*coords).terms()if c}for i,j in ambient.PAIRS if P[i,j]!=0}
    D=s.Matrix(ambient.linearized_jacobi_matrix(sparse).tolist())
    K=dm(D).nullspace().to_Matrix().T
    assert D*K==s.zeros(350,K.cols)
    directions=[]
    for i,j in it.product(range(5),repeat=2):
        A=s.zeros(5);A[i,j]=1;X=A*s.Matrix(coords)
        directions.append(vector(sum([X[k]*P.diff(coords[k])for k in range(5)],s.zeros(5))-A*P-P*A.T))
    M=s.Matrix(ambient.vertical_bivector_matrix().tolist()).row_join(s.Matrix.hstack(*directions,vector(P)))
    assert D*M==s.zeros(350,M.cols)
    if point:
        dcert,mcert,kcert=minor(D),minor(M),minor(K)
        assert dcert['rank']+K.cols==150 and kcert['rank']==K.cols
        result=dict(input_report=str(source_path),input_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
            helper_path=str(helper),helper_sha256=hashlib.sha256(helper.read_bytes()).hexdigest(),
            point=point,full150variables_and350rows=True,integer_scale=str(denominator),
            witness_brackets={f'{i},{j}':str(P[i,j])for i,j in ambient.PAIRS if P[i,j]!=0},
            derivative=dcert,orbit_vertical_scale=mcert,whole_kernel=kcert,
            full_exact_kernel_sparse_columns=[encode(K[:,j])for j in range(K.cols)],
            all_exact_D_kernel_products_zero=True,
            projective_orbit_dimension=mcert['rank']-25,projective_tangent_dimension=125-dcert['rank'],
            smooth_orbit_component_proved=dcert['rank']+mcert['rank']==150,
            matrix_sha256=dict(D=digest(D),orbit_vertical_scale=digest(M),kernel=digest(K)))
        if point.get('family_report'):
            family_path=Path(point['family_report'])
            family=json.loads(family_path.read_text())[point['family_result_key']]
            if 'family_index' in point:family=family['all_three_complete_families'][point['family_index']]
            parameters=[s.Symbol(v)for v in family['parameters']]
            parameter_env={str(v):v for v in parameters}
            values={s.Symbol(k):s.Rational(v)for k,v in point['family_parameter_point'].items()}
            family_matrix=s.zeros(5)
            for (i,j),text in zip(it.combinations(range(1,5),2),family['universal_brackets']):
                polynomial=s.Poly(s.sympify(text,locals=parameter_env),x,y,z,w)
                entry=sum(c*h**(2-sum(e))*s.prod(v**a for v,a in zip((x,y,z,w),e))for e,c in polynomial.terms())
                family_matrix[i,j]=denominator*entry;family_matrix[j,i]=-denominator*entry
            assert (family_matrix.subs(values)-P).applyfunc(s.expand)==s.zeros(5)
            family_columns=s.Matrix.hstack(*[vector(family_matrix.diff(v).subs(values))for v in parameters])
            actual=M.row_join(family_columns)
            assert D*actual==s.zeros(350,actual.cols)
            family_cert=minor(actual)
            result.update(family_report=str(family_path),family_report_sha256=hashlib.sha256(family_path.read_bytes()).hexdigest(),
                actual_full_family=family_cert,all_actual_family_directions_annihilated=True,
                actual_family_projective_rank=family_cert['rank']-25,
                smooth_full_source_component_proved=dcert['rank']+family_cert['rank']==150)
            if quadratic:
                pivots=list(dm(actual.row_join(K)).rref()[1])
                extra=[j-actual.cols for j in pivots if j>=actual.cols]
                normals=[]
                for at in extra:
                    value=K[:,at];common=s.ilcm(*[v.q for v in value]);value*=common
                    divisor=math.gcd(*(int(v)for v in value));normals.append(value/divisor)
                normal_matrix=s.Matrix.hstack(*normals)
                matrices=[unvector(v)for v in normals]
                pair_order=list(it.combinations_with_replacement(range(len(extra)),2))
                Q=s.Matrix.hstack(*[jac(matrices[i])if i==j else jac(matrices[i]+matrices[j])-jac(matrices[i])-jac(matrices[j])for i,j in pair_order])
                rr,cc=dcert['rows'],dcert['columns']
                projected=Q-D[:,cc]*dm(D[rr,cc]).inv().matmul(dm(Q[rr,:])).to_Matrix()
                independent_rows=list(dm(projected.T).rref()[1])
                variables=s.symbols(f'e0:{len(extra)}')
                polynomials=[s.expand(sum(projected[row,j]*variables[i]*variables[k]for j,(i,k)in enumerate(pair_order)))for row in independent_rows]
                basis=s.groebner(polynomials,*variables,domain=s.QQ)
                zero_dimensional=bool(polynomials)and basis.is_zero_dimensional
                result.update(true_normal_count=len(extra),independent_normal_sparse_vectors=[encode(v)for v in normals],
                    quadratic_monomial_pairs=[list(v)for v in pair_order],
                    full_quadratic_sparse_columns=[encode(Q[:,j])for j in range(Q.cols)],
                    quadratic_cokernel_rows=independent_rows,all_quadratic_residual_generators=list(map(str,polynomials)),
                    exact_quadratic_groebner_basis=[str(v.as_expr())for v in basis.polys],
                    quadratic_initial_support_zero_dimensional=zero_dimensional,
                    reduced_full_source_component_proved_by_true_slice=zero_dimensional,
                    point_unique_component_not_claimed=True)
            if twist:
                Z=s.zeros(5,1)
                for i in range(4):Z[i+1]=row['field_constant'][i]*h+sum(row['field_linear_matrix'][i][j]*coords[j+1]for j in range(4))
                N=Z.jacobian(coords);Rcolumns=[];commutators=[]
                for i,j in it.product(range(5),repeat=2):
                    B=s.zeros(5);B[i,j]=1;Y=B*s.Matrix(coords)
                    Rcolumns.append(vector(Z*Y.T-Y*Z.T))
                    commutators.append(s.Matrix(N*B-B*N).reshape(25,1))
                W=s.Matrix.hstack(*Rcolumns)
                linear=s.Matrix.hstack(*commutators).col_join(D*W)
                Ykernel=dm(linear).nullspace().to_Matrix().T
                assert linear*Ykernel==s.zeros(linear.rows,Ykernel.cols)
                allowed=W*Ykernel
                Hcolumns=[i for i in range(Ykernel.cols)if any(Ykernel[j,i]!=0 for j in range(1,5))]
                source_rank=family_cert['rank']
                transverse=[i for i in range(Ykernel.cols)if minor(actual.row_join(allowed[:,i]))['rank']>source_rank]
                info=dict(complete_linear_system_shape=list(linear.shape),commuting_compatible_Y_dimension=Ykernel.cols,
                    complete_Y_kernel_sparse_columns=[encode(Ykernel[:,i])for i in range(Ykernel.cols)],
                    H_removing_columns=Hcolumns,source_transverse_columns=transverse,
                    H_removing_transverse_twist_exists=bool(Hcolumns and transverse))
                if Hcolumns and transverse:
                    at=Hcolumns[0];value=Ykernel[:,at];Rvec=allowed[:,at]
                    if minor(actual.row_join(Rvec))['rank']==source_rank:
                        other=transverse[0]
                        multiplier=1
                        if not any((value+Ykernel[:,other])[j]!=0 for j in range(1,5)):multiplier=2
                        value+=multiplier*Ykernel[:,other];Rvec=W*value
                    assert any(value[j]!=0 for j in range(1,5))
                    R=unvector(Rvec)
                    assert jac(R)==s.zeros(350,1) and D*Rvec==s.zeros(350,1)
                    info.update(Y_matrix=[[str(value[5*i+j])for j in range(5)]for i in range(5)],
                        compatible_R_brackets={f'{i},{j}':str(R[i,j])for i,j in ambient.PAIRS if R[i,j]!=0},
                        all_parameter_projective_Poisson_pencil=True,
                        actual_source_plus_R=minor(actual.row_join(Rvec)))
                result['commuting_twist']=info
        return result
    pivots=list(dm(M.row_join(K)).rref()[1]);extra=[j-M.cols for j in pivots if j>=M.cols]
    if second:
        assert len(extra)==3 and K.cols==50
        normals=[]
        for at in extra:
            value=K[:,at]
            common=s.ilcm(*[v.q for v in value]);value*=common
            divisor=math.gcd(*(int(v)for v in value));value/=divisor
            normals.append(value)
        normal_matrix=s.Matrix.hstack(*normals)
        matrices=[unvector(v)for v in normals]
        pair_order=list(it.combinations_with_replacement(range(3),2))
        quadratic_columns=[jac(matrices[i]) if i==j else jac(matrices[i]+matrices[j])-jac(matrices[i])-jac(matrices[j])
                           for i,j in pair_order]
        Q=s.Matrix.hstack(*quadratic_columns)
        dcert,mcert,kcert=minor(D),minor(M),minor(M.row_join(normal_matrix))
        assert dcert['rank']==100 and mcert['rank']==47 and kcert['rank']==50
        rr,cc=dcert['rows'],dcert['columns']
        projected=Q-D[:,cc]*dm(D[rr,cc]).inv().matmul(dm(Q[rr,:])).to_Matrix()
        independent_rows=list(dm(projected.T).rref()[1])
        variables=s.symbols('u v w')
        polynomials=[s.expand(sum(projected[row,j]*variables[i]*variables[k]
                                 for j,(i,k) in enumerate(pair_order)))for row in independent_rows]
        basis=s.groebner(polynomials,*variables,domain=s.QQ)
        zero_dimensional=bool(polynomials) and basis.is_zero_dimensional
        return dict(input_report=str(source_path),input_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
            helper_path=str(helper),helper_sha256=hashlib.sha256(helper.read_bytes()).hexdigest(),
            full150variables_and350rows=True,point_substitutions={str(k):str(v)for k,v in substitution.items()if v},
            integer_scale=str(denominator),witness_brackets={f'{i},{j}':str(P[i,j])for i,j in ambient.PAIRS if P[i,j]!=0},
            derivative=dcert,orbit_vertical_scale=mcert,whole_kernel=kcert,
            independent_normal_sparse_vectors=[encode(v)for v in normals],all_exact_D_kernel_products_zero=True,
            quadratic_monomial_pair_order=[list(p)for p in pair_order],
            full_quadratic_columns=[encode(Q[:,j])for j in range(Q.cols)],
            full_D_projection_rows=rr,full_D_projection_columns=cc,
            independent_cokernel_rows=independent_rows,all_quadratic_residual_generators=list(map(str,polynomials)),
            exact_quadratic_groebner_basis=[str(p.as_expr())for p in basis.polys],
            quadratic_initial_ideal_zero_dimensional=zero_dimensional,
            projective_orbit_dimension=22,projective_tangent_dimension=25,
            support_component_proved=zero_dimensional,
            next_step=None if zero_dimensional else 'Use genuine source geometry in the surviving normal directions; do not infer support from only these quadrics or repeat higher jets.',
            matrix_sha256=dict(D=digest(D),orbit_vertical_scale=digest(M),normals=digest(normal_matrix),quadratic=digest(Q),projected_quadratic=digest(projected)))
    assert len(extra)==1 and K.cols==48
    Rvec=K[:,extra[0]]
    common=s.ilcm(*[v.q for v in Rvec]);Rvec*=common
    divisor=math.gcd(*(int(v)for v in Rvec));Rvec/=divisor
    R=unvector(Rvec)
    assert D*Rvec==s.zeros(350,1)
    Q=jac(R)
    dcert,mcert,kcert,qcert=minor(D),minor(M),minor(M.row_join(Rvec)),minor(D.row_join(Q))
    assert dcert['rank']==102 and mcert['rank']==47 and kcert['rank']==48
    result=dict(input_report=str(source_path),input_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
        helper_path=str(helper),helper_sha256=hashlib.sha256(helper.read_bytes()).hexdigest(),
        full150variables_and350rows=True,point_substitutions={str(k):str(v)for k,v in substitution.items()if v},
        integer_scale=str(denominator),witness_brackets={f'{i},{j}':str(P[i,j])for i,j in ambient.PAIRS if P[i,j]!=0},
        derivative=dcert,orbit_vertical_scale=mcert,whole_kernel=kcert,
        independent_normal_sparse_vector=encode(Rvec),all_exact_D_kernel_products_zero=True,
        quadratic_obstruction=qcert,quadratic_coefficient_sparse=encode(Q),
        matrix_sha256=dict(D=digest(D),orbit_vertical_scale=digest(M),normal=digest(Rvec),quadratic=digest(Q)),
        projective_orbit_dimension=22,projective_tangent_dimension=23,
        support_component_proved=qcert['rank']>dcert['rank'])
    if result['support_component_proved']:
        result['first_nonzero_obstruction_order']=2
        return result
    # One actual quadratic correction; do not continue beyond the requested
    # cubic support test or claim an exponent from a vanishing obstruction.
    Tvec=s.zeros(150,1)
    rr,cc=dcert['rows'],dcert['columns']
    solution=dm(D[rr,cc]).inv().matmul(dm(-Q[rr,0])).to_Matrix()
    for i,value in zip(cc,solution):Tvec[i]=value
    assert D*Tvec+Q==s.zeros(350,1)
    T=unvector(Tvec)
    cubic=jac(R+T)-jac(R)-jac(T)
    mixed=s.Matrix.hstack(*[jac(R+unvector(K[:,j]))-jac(R)-jac(unvector(K[:,j]))for j in range(K.cols)])
    mixed_cert=minor(D.row_join(mixed));cubic_cert=minor(D.row_join(mixed).row_join(cubic))
    result.update(quadratic_correction_sparse=encode(Tvec),exact_full_quadratic_correction=True,
                  all_kernel_correction_mixed_rank=mixed_cert,cubic_obstruction=cubic_cert,
                  cubic_coefficient_sparse=encode(cubic),
                  support_component_proved=cubic_cert['rank']>mixed_cert['rank'])
    if result['support_component_proved']:
        result['first_nonzero_obstruction_order']=3
    else:result['next_step']='Quadratic and cubic obstruction routes did not close support; change to finite weighted-slice or global geometric route.'
    return result


def nil_regular_graph(source_path,point=None):
    point_spec=point
    """Exact whole T39/T40 graph test at the second regular-nil witness."""
    from s61_model import S61,X,pushforward
    rows=json.loads(source_path.read_text())['full_nilpotent_sources']
    n=s.symbols('n0:36');coords=s.symbols('h x y z w')
    substitution={v:0 for v in n}
    substitution.update({n[8]:-s.Rational(1,4),n[12]:-s.Rational(1,2),n[10]:1,n[6]:3})
    P=s.zeros(5)
    if point:
        substitution={v:0 for v in n}
        substitution.update({s.Symbol(k):s.Rational(v)for k,v in point['substitutions'].items()})
    row=rows[point['source_row'] if point else 0]
    for (i,j),text in zip(it.combinations(range(1,5),2),row['brackets']):
        f=s.Poly(s.sympify(text).subs(substitution),*coords[1:])
        P[i,j]=sum(c*coords[0]**(2-sum(e))*s.prod(v**a for v,a in zip(coords[1:],e))for e,c in f.terms())
        P[j,i]=-P[i,j]
    P=P.subs(dict(zip(coords,X)),simultaneous=True)
    translation=s.eye(5)
    if point:
        for i,value in enumerate(point.get('evaluation',[0,0,0,0]),1):translation[i,0]=-s.Rational(value)
        P=pushforward(P,translation)
    K=P.subs(dict(zip(X,(1,0,0,0,0))))[1:5,1:5]
    assert K.det()!=0
    unit=[s.eye(4)[:,i]for i in range(4)]
    i,j=next((i,j)for i,j in it.combinations(range(4),2)if K[i,j])
    r0=unit[i];r1=unit[j]/K[i,j];remaining=[]
    for r in unit:
        q=r-(r.T*K*r1)[0]*r0+(r.T*K*r0)[0]*r1
        if q!=s.zeros(4,1) and (not remaining or s.Matrix.hstack(*remaining,q).rank()>len(remaining)):
            remaining.append(q)
    assert len(remaining)==2
    r2,r3=remaining;r3=r3/(r2.T*K*r3)[0]
    T=s.Matrix.vstack(*(r.T for r in(r0,r1,r2,r3)))
    J=s.zeros(4);J[0,1]=J[2,3]=1;J[1,0]=J[3,2]=-1
    assert T*K*T.T==J
    frame=s.eye(5);frame[1:5,1:5]=T
    model=S61();normalized=model.normalize(pushforward(P,frame))
    model.verify_point(normalized,check_equations=True)
    point=normalized['parameters'];us=s.symbols('u1:62');graphs={}
    paths=[source_path,Path('code/s61_model.py')]
    for label,path in [('T39',Path('reports/Q51_UNIVERSAL_SECTION_20261101.json')),
                       ('T40',Path('reports/Q52_UNIVERSAL_SECTION_20261101.json'))]:
        graph=json.loads(path.read_text());paths.append(path)
        residual=[s.cancel(point[k]-s.sympify(q).subs(dict(zip(us,point)),simultaneous=True))
                  for k,q in zip(graph['weight2_u_indices_zero_based'],graph['quadratic_section_q'])]
        graphs[label]=dict(all45_residuals=list(map(str,residual)),whole_graph_membership=all(v==0 for v in residual))
    return dict(input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
                source_point=point_spec,translation=[[str(v)for v in translation.row(i)]for i in range(5)],
                nondegenerate_constant_matrix=[[str(v)for v in K.row(i)]for i in range(4)],
                Darboux_frame=[[str(v)for v in frame.row(i)]for i in range(5)],
                perspective=list(map(str,normalized['perspective'])),all61_parameters=list(map(str,point)),
                full_quotient_reconstruction_and_all134_equations_checked=True,graph_membership=graphs,
                implication='If a known25component graph matches, combine with the separate fullambient tangent bound25 to prove local regularity/uniqueness, then use irreducibility for the complete source.')


def nil_short_translation_geometry(source_path):
    """Entire28-parameter nil211,translation-length1 potential reduction."""
    row=json.loads(source_path.read_text())['full_nilpotent_sources'][8]
    n=s.symbols('n0:36');x,y,z,w=s.symbols('x y z w');u,t=s.symbols('u t')
    old=(x,y,z,w);new=(u,y,w,t);pairs=list(it.combinations(range(4),2));P=s.zeros(4)
    for (i,j),text in zip(pairs,row['brackets']):P[i,j]=s.sympify(text);P[j,i]=-P[i,j]
    transformation=s.Matrix([[1,-z,-y,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
    assert transformation.det()==-1
    inverse={x:u+y*t,z:t}
    changed=(transformation*P*transformation.T).subs(inverse,simultaneous=True).applyfunc(s.expand)
    assert all(s.diff(changed[i,j],t)==0 for i,j in pairs)
    c=n[18]+n[19]*y+n[20]*w-n[25]*y*y/3
    Phi=u*c+n[0]*w+n[1]*w*w/2+n[2]*w**3/3-n[11]*y*w+(n[23]-n[13])*y*y*w/2-n[14]*y*w*w\
        +n[26]*y**3*w/3+n[27]*y*y*w*w/2-n[9]*y-(n[10]-n[21])*y*y/2-(n[12]-n[22])*y**3/3+n[24]*y**4/4
    yy=(u,y,w)
    for ij,derivative in [((0,1),s.diff(Phi,w)),((0,2),-s.diff(Phi,y)),((1,2),s.diff(Phi,u))]:
        assert s.expand(changed[ij]-derivative)==0
    V=s.Matrix([changed[3,j]for j in range(3)])
    assert s.expand(sum(s.diff(V[i],yy[i])for i in range(3))+1)==0
    image=s.expand(sum(V[i]*s.diff(Phi,yy[i])for i in range(3)))
    affine_pf=s.sympify(row['actual_affine_pfaffian']).subs(inverse,simultaneous=True)
    assert s.expand(affine_pf-image)==0
    Q=s.expand(image+Phi)
    def strict(i,j,k):
        return s.expand(sum(changed[i,l]*s.diff(changed[j,k],new[l])+changed[j,l]*s.diff(changed[k,i],new[l])+
                            changed[k,l]*s.diff(changed[i,j],new[l])for l in range(4)))
    assert strict(0,1,2)==0
    for ij,derivative in [((0,1),s.diff(Q,w)),((0,2),-s.diff(Q,y)),((1,2),s.diff(Q,u))]:
        assert s.expand(strict(*ij,3)-derivative)==0
    for exponent,coefficient in s.Poly(Phi,u,y,w).terms():assert 2*exponent[0]+exponent[1]+exponent[2]<=4
    point={v:0 for v in n}
    point.update({n[0]:1,n[13]:-1,n[16]:1,n[18]:-s.Rational(3,2),n[23]:-1,n[25]:-3})
    assert all(s.expand(s.diff(Q,v).subs(point))==0 for v in yy)
    point_pf=s.factor(affine_pf.subs(point))
    assert s.expand(point_pf+(y*y-s.Rational(3,2))*u+w)==0
    return dict(input_report=str(source_path),input_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                entire28_parameter_source_preserved=True,straightened_coordinates=['u=x-y*z','y','w','t=z'],
                inverse_coordinates=['x=u+y*t','y=y','z=t','w=w'],transformation_determinant=-1,
                transformed_brackets=[str(changed[i,j])for i,j in pairs],all_brackets_independent_of_t=True,
                potential=str(Phi),vector_field=list(map(str,V)),vector_field_divergence=-1,
                weighted_degree_weights=[2,1,1],potential_weighted_degree_at_most4=True,
                actual_affine_pfaffian_equals_V_potential=True,
                all_jacobi_equivalent_to_three_derivatives_of=str(Q),
                support_equation='V(Phi)+Phi=constant, equivalently originalPf=-Phi+constant',
                retained_point_substitutions={str(k):str(v)for k,v in point.items()if v},retained_point_pfaffian=str(point_pf),
                all_point_jacobi_zero=True,
                scope='Complete potential reduction only; full eigenpotential support and component membership remain to classify.')


def nil22_long_dimension(source_path,core_only=False,fixed_n14=None,a_rank_one=False):
    """Whole type7 dimension bound, certified by lifted leading monomials."""
    import os
    import subprocess
    import tempfile
    from functools import lru_cache
    row=json.loads(source_path.read_text())['full_nilpotent_sources'][6]
    n=s.symbols('n0:20')
    generators=[];core_indices=[2,8,11,12,14,17,18,19]
    for text in row['distinct_nonzero_jacobi_coefficients']:
        p=s.Poly(s.sympify(text),*n,domain=s.QQ).clear_denoms()[1].primitive()[1]
        if core_only and not p.as_expr().free_symbols<={n[i]for i in core_indices}:continue
        if p.LC()<0:p=-p
        if p not in generators:generators.append(p)
    if fixed_n14 is not None:generators.append(s.Poly(n[14]-fixed_n14,*n,domain=s.QQ))
    if a_rank_one:
        J=list(map(s.sympify,row['distinct_nonzero_jacobi_coefficients']))
        Amatrix,_=s.linear_eq_to_matrix([J[i]for i in [2,8]],[n[i]for i in [1,7,10,13,16]])
        for columns in it.combinations(range(5),2):
            p=s.Poly(Amatrix[:,list(columns)].det(),*n,domain=s.QQ)
            if p!=0:generators.append(p)
    source='ring R=0,('+','.join(map(str,n))+'),dp;\nshort=0;\n'
    source+='ideal I='+','.join(str(p.as_expr()).replace('**','^')for p in generators)+';\n'
    source+='matrix T; ideal G=liftstd(I,T);\n'
    source+='print("DIM:"+string(dim(G)));\n'
    source+='int i;int j;for(i=1;i<=size(G);i++){print("G:"+string(i)+":"+string(G[i]));}\n'
    source+='for(i=1;i<=nrows(T);i++){for(j=1;j<=ncols(T);j++){if(T[i,j]!=0){print("L:"+string(i)+":"+string(j)+":"+string(T[i,j]));}}}\nquit;\n'
    prefix=Path('work/nonsquare-v62/singular-runtime/prefix').resolve()
    env=dict(os.environ)
    env['LD_LIBRARY_PATH']='/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu:'+str(prefix/'usr/lib/x86_64-linux-gnu')+':'+env.get('LD_LIBRARY_PATH','')
    env['SINGULAR_PROCS_DIR']=str(prefix/'usr/lib/x86_64-linux-gnu/singular/MOD')
    with tempfile.TemporaryDirectory(prefix='quartic-nil22-dim-')as directory:
        path=Path(directory)/'dimension.sing';path.write_text(source)
        started=time.monotonic()
        run=subprocess.run([str(prefix/'usr/bin/Singular'),'-q',str(path)],env=env,
                           capture_output=True,text=True,timeout=80,check=False)
        assert run.returncode==0,run.stderr
        elapsed=time.monotonic()-started
    rows=run.stdout.splitlines();G={};lifts={};singular_dimension=None
    for line in rows:
        if line.startswith('G:'):
            _,j,text=line.split(':',2);G[int(j)]=s.Poly(s.sympify(text.replace('^','**')),*n,domain=s.QQ)
        if line.startswith('L:'):
            _,i,j,text=line.split(':',3);lifts[int(i),int(j)]=s.Poly(s.sympify(text.replace('^','**')),*n,domain=s.QQ)
        if line.startswith('DIM:'):singular_dimension=int(line.split(':')[1])
    assert G and singular_dimension is not None and '?'not in run.stdout
    zero=s.Poly(0,*n,domain=s.QQ)
    for j,g in G.items():
        rhs=sum((generators[i-1]*coefficient for(i,k),coefficient in lifts.items()if k==j),zero)
        assert rhs==g
    leading=[g.monoms(order='grevlex')[0]for g in G.values()]
    supports=tuple(sorted(set(sum(1<<i for i,e in enumerate(m)if e)for m in leading)))
    @lru_cache(None)
    def minimum_cover(edges):
        if not edges:return 0
        edge=min(edges,key=int.bit_count)
        assert edge
        return 1+min(minimum_cover(tuple(e for e in edges if not e&(1<<v)))for v in range(20)if edge&(1<<v))
    height=21 if 0 in supports else minimum_cover(supports)
    bound=20-height
    assert bound==singular_dimension
    A=s.Matrix(row['field_linear_matrix']);c=s.Matrix(row['field_constant'])
    B=s.Matrix(4,4,s.symbols('b0:16'));translation=s.Matrix(s.symbols('t0:4'))
    affine_orbit=s.Matrix(list(B*A-A*B)+list(B*c-A*translation)).jacobian(list(B)+list(translation))
    linear_rank=s.Matrix(list(B*A-A*B)).jacobian(list(B)).rank()
    orbit_rank=affine_orbit.rank()
    assert linear_rank==8 and orbit_rank==12
    result=dict(input_report=str(source_path),input_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
        entire20_parameter_nil22_top2_source_preserved=True,localizations=[],all_complements_retained=True,
        primitive_full_Jacobi_generators=[str(p.as_expr())for p in generators],
        exact_Singular_input=source,exact_Singular_stdout=run.stdout,Singular_stderr=run.stderr,Singular_seconds=elapsed,
        lifted_polynomials=[str(G[j].as_expr())for j in sorted(G)],
        exact_ideal_membership_lifts=[dict(generator=i-1,polynomial=j-1,coefficient=str(v.as_expr()))for(i,j),v in sorted(lifts.items())],
        all_lift_identities_independently_checked_in_QQ=True,
        leading_monomial_exponents=[list(m)for m in leading],leading_monomial_minimum_vertex_cover=height,
        minimum_cover_recursion_states=minimum_cover.cache_info().currsize,
        fixed_source_dimension_upper=bound,affine_cofactor_orbit_dimension=orbit_rank,
        exact_affine_orbit_tangent_matrix=[[str(v)for v in affine_orbit.row(i)]for i in range(20)],
        hyperplane_choice_dimension=4,projectivization_subtracts=1,
        saturated_projective_dimension_upper=bound+orbit_rank+4-1,
        excludes_unknown_G25_to27=(bound+orbit_rank+4-1<25),exceptional_membership_not_proved=True,
        scope='Whole source dimension bound from exact ideal-membership identities and leading-monomial height; no restriction to irreducible-quartic points, sampling, or generic-rank inference.')
    if core_only:
        result.update(core_variable_indices=core_indices,core_dimension_upper=bound-12 if bound>=0 else -1,
                      fixed_n14_value=str(fixed_n14)if fixed_n14 is not None else None,
                      imposed_A_rank_at_most_one=a_rank_one,
                      full_source_dimension_not_yet_proved=True,
                      scope='Complete nonlinear eightvariable core only; all12fiber variables and their remaining equations retained as unresolved. No generic non-dominance conclusion yet.')
    return result


def nil22_long_support(source_path):
    """Uniform triangular fibers over the complete two-dimensional type7 core."""
    row=json.loads(source_path.read_text())['full_nilpotent_sources'][6]
    n=s.symbols('n0:20');J=list(map(s.sympify,row['distinct_nonzero_jacobi_coefficients']))
    a_variables=[n[i]for i in [1,7,10,13,16]]
    b_variables=[n[i]for i in [0,3,4,5,6,9,15]]
    Amatrix,Arhs=s.linear_eq_to_matrix([J[i]for i in [2,8]],a_variables)
    Bmatrix,Brhs=s.linear_eq_to_matrix([J[i]for i in [13,23,9,14,15]],b_variables)
    assert Arhs==s.zeros(2,1)
    # This constant identity shows rankA>=1 for every core, without localization.
    assert s.expand(Amatrix[1,1]+6*Amatrix[0,4]-2*Amatrix[0,3]-2*Amatrix[1,2])==-4
    Cmatrix=Bmatrix.extract([1,2],[0,4,5,6])
    assert Bmatrix.extract([1,2],[3,2,1])==s.zeros(2,3)
    c_unknowns=[n[i]for i in [17,18,11,14,19,12]]
    forced_matrix,forced_rhs=s.linear_eq_to_matrix(list(Cmatrix),c_unknowns)
    forced=s.Matrix([0,0,2,1,s.Rational(5,2),-3*n[8]])
    assert forced_matrix.rank()==6 and forced_matrix*forced==forced_rhs
    forced_rows=list(forced_matrix.T.rref()[1])
    forced_minor=s.factor(forced_matrix.extract(forced_rows,range(6)).det())
    assert forced_minor!=0 and not forced_minor.free_symbols
    contradiction=s.expand(J[25].subs(dict(zip(c_unknowns,forced)),simultaneous=True))
    assert contradiction==3
    regular_minor=s.factor(Bmatrix.extract([0,3,4],[3,2,1]).det())
    assert s.expand(regular_minor+(n[14]+1)*(2*n[14]+1)*(3*n[14]+1))==0
    critical_minors=[]
    for value,rows,columns,expected in [(-s.Integer(1),[0,3],[3,2],s.Integer(2)),
            (-s.Rational(1,2),[0,4],[3,1],-s.Rational(1,4)),
            (-s.Rational(1,3),[3,4],[2,1],s.Rational(2,9))]:
        minor=s.factor(Bmatrix.extract(rows,columns).det().subs(n[14],value))
        assert minor==expected
        critical_minors.append(dict(n14=str(value),rows=rows,columns=columns,value=str(minor)))
    core_path=Path('reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_CORE_20261107.json')
    critical_path=Path('reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_CRITICAL_CORES_20261107.json')
    low_path=Path('reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_LOW_RANK_CORE_20261107.json')
    core=json.loads(core_path.read_text())['nil22_long_core']
    critical=json.loads(critical_path.read_text())['nil22_long_critical_cores']
    low=json.loads(low_path.read_text())['nil22_long_low_rank_core']
    input_hash=hashlib.sha256(source_path.read_bytes()).hexdigest()
    assert all(c['input_sha256']==input_hash for c in [core,*critical,low])
    assert core['core_dimension_upper']==2
    assert [(c['fixed_n14_value'],c['core_dimension_upper'])for c in critical]==[('-1',1),('-1/2',2),('-1/3',1)]
    assert low['fixed_n14_value']=='-1/2'and low['imposed_A_rank_at_most_one']and low['core_dimension_upper']==1
    strata=[dict(name='n14 outside {-1,-1/2,-1/3}',core_dimension=2,rank_A_at_least=1,rank_B_at_least=4),
            dict(name='n14=-1',core_dimension=1,rank_A_at_least=1,rank_B_at_least=3),
            dict(name='n14=-1/3',core_dimension=1,rank_A_at_least=1,rank_B_at_least=3),
            dict(name='n14=-1/2,rankA=2',core_dimension=2,rank_A_at_least=2,rank_B_at_least=3),
            dict(name='n14=-1/2,rankA=1',core_dimension=1,rank_A_at_least=1,rank_B_at_least=3)]
    for case in strata:
        case['first_fiber_dimension_upper']=5-case['rank_A_at_least']
        case['second_fiber_dimension_upper']=7-case['rank_B_at_least']
        case['total_source_dimension_upper']=case['core_dimension']+case['first_fiber_dimension_upper']+case['second_fiber_dimension_upper']
        assert case['total_source_dimension_upper']==9
    A=s.Matrix(row['field_linear_matrix']);c=s.Matrix(row['field_constant'])
    B=s.Matrix(4,4,s.symbols('b0:16'));translation=s.Matrix(s.symbols('t0:4'))
    orbit=s.Matrix(list(B*A-A*B)+list(B*c-A*translation)).jacobian(list(B)+list(translation))
    assert orbit.rank()==12
    return dict(input_report=str(source_path),input_sha256=input_hash,
        dimension_certificate_dependencies=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())for p in [core_path,critical_path,low_path]],
        entire20_parameter_nil22_top2_source_preserved=True,all_boundary_strata_retained=True,
        first_fiber_variables=list(map(str,a_variables)),second_fiber_variables=list(map(str,b_variables)),
        A_matrix=[[str(v)for v in Amatrix.row(i)]for i in range(2)],
        A_rank_at_least_one_constant_identity='A[1,1]+6*A[0,4]-2*A[0,3]-2*A[1,2]=-4',
        B_matrix=[[str(v)for v in Bmatrix.row(i)]for i in range(5)],B_right_hand_side=list(map(str,Brhs)),
        C_matrix=[[str(v)for v in Cmatrix.row(i)]for i in range(2)],
        C_zero_forced_substitution={str(k):str(v)for k,v in zip(c_unknowns,forced)},
        C_zero_forcing_minor=dict(rows=forced_rows,value=str(forced_minor)),
        C_zero_contradiction_core_coefficient=25,C_zero_contradiction_value='3',
        B_regular_three_pivot_minor=str(regular_minor),B_critical_two_pivot_minors=critical_minors,
        complete_dimension_strata=strata,fixed_source_dimension_upper=9,
        affine_cofactor_orbit_dimension=12,exact_affine_orbit_tangent_matrix=[[str(v)for v in orbit.row(i)]for i in range(20)],
        hyperplane_choice_dimension=4,projectivization_subtracts=1,saturated_projective_dimension_upper=24,
        excludes_unknown_G25_to27=True,exceptional_membership_not_proved=True,
        scope='Whole20parameter source including every reducible/zero Pfaffian boundary has projective saturateddimension<=24. No component membership or newcomponent assertion.')


def nil31_long_support(source_path):
    """Complete nil31longtranslation support: empty complement plus two full5sources."""
    row=json.loads(source_path.read_text())['full_nilpotent_sources'][4]
    n=s.symbols('n0:36');x,y,z,w=s.symbols('x y z w');a=s.Symbol('a')
    J=[s.sympify(v)for v in row['distinct_nonzero_jacobi_coefficients']]
    Pf=s.sympify(row['actual_affine_pfaffian'])
    K=s.sympify(row['highest_quartic_pfaffian'])
    assert s.expand(K.subs({n[14]:0,n[11]:0,n[9]:0,n[7]:0},simultaneous=True))==0
    # On n14=0, n11!=0 would force n13=-1, then n9=n7=0,
    # and the fourth Jacobi coefficient is -2*n11^2/3. Hence n11=0.
    impossible={n[14]:0,n[13]:-1,n[9]:0,n[7]:0}
    assert s.expand(J[3].subs(impossible,simultaneous=True)+2*n[11]**2/3)==0
    # n9!=0 then forces n13=-3/8. Every other forced coefficient is
    # triangular; its actual Pfaffian has a quadratic square factor.
    boundary={n[14]:0,n[13]:-s.Rational(3,8),n[11]:0,n[7]:0,n[10]:0,
              n[2]:-24*n[9]**2,n[6]:s.Rational(16,3)*n[12]*n[9],
              n[1]:s.Rational(448,3)*n[12]*n[9]**2+8*n[8]*n[9],
              n[5]:-s.Rational(64,9)*n[12]**2*n[9],
              n[0]:-s.Rational(2048,9)*n[12]**2*n[9]**2-s.Rational(64,3)*n[12]*n[8]*n[9]}
    assert all(s.expand(v.subs(boundary,simultaneous=True))==0 for v in J)
    boundary_pf=s.factor(Pf.subs(boundary,simultaneous=True))
    assert s.expand(boundary_pf+n[9]*(-8*n[12]+3*w)**2*(64*n[12]*n[9]+8*n[8]-16*n[9]*w+2*y-z*z)/18)==0
    boundary_vars=[n[2],n[6],n[1],n[5],n[0]]
    boundary_base={n[14]:0,n[13]:-s.Rational(3,8),n[11]:0,n[7]:0,n[10]:0}
    boundary_matrix,boundary_rhs=s.linear_eq_to_matrix([v.subs(boundary_base,simultaneous=True)for v in J],boundary_vars)
    boundary_rows=[3,7,8,9,18]
    boundary_minor=s.factor(boundary_matrix.extract(boundary_rows,range(5)).det())
    assert boundary_minor==-s.Rational(45,8192)
    assert (boundary_matrix*s.Matrix([boundary[v]for v in boundary_vars])-boundary_rhs).applyfunc(s.expand)==s.zeros(len(J),1)
    empty={n[14]:0,n[13]:s.Rational(1,4),n[10]:0,n[11]:0,n[9]:0,n[2]:0,n[6]:0,n[1]:0,n[0]:0,
           n[5]:-4*n[7]*n[8]}
    empty_pf=s.factor(Pf.subs(empty,simultaneous=True))
    I=2*y-z*z
    assert s.expand(empty_pf+n[7]*(I+4*n[8])**2/4)==0
    assert all(s.expand(v.subs(empty,simultaneous=True))==0 for v in J)
    empty_vars=[n[10],n[2],n[6],n[1],n[0],n[5]]
    empty_base={n[14]:0,n[13]:s.Rational(1,4),n[11]:0,n[9]:0}
    empty_matrix,empty_rhs=s.linear_eq_to_matrix([v.subs(empty_base,simultaneous=True)for v in J],empty_vars)
    empty_rows=[3,7,8,9,14,18]
    empty_minor=s.factor(empty_matrix.extract(empty_rows,range(6)).det())
    assert empty_minor==-s.Rational(315,1024)
    assert (empty_matrix*s.Matrix([empty[v]for v in empty_vars])-empty_rhs).applyfunc(s.expand)==s.zeros(len(J),1)
    # n14!=0 forces n13=-1/6. Four fullJacobi coefficients then force
    # n11,n9,n7,n10 exactly, leaving one scalar product for n2.
    first={n[14]:a,n[13]:-s.Rational(1,6),n[11]:-9*a*n[4],n[9]:6*a*n[2],
           n[7]:-6*a*(n[6]+2*n[8]),n[10]:6*a*(3*n[4]*n[12]-n[3])}
    initial_vars=[n[11],n[9],n[7],n[10]]
    initial_matrix,initial_rhs=s.linear_eq_to_matrix([J[i].subs({n[14]:a,n[13]:-s.Rational(1,6)},simultaneous=True)for i in [0,2,13,14]],initial_vars)
    assert initial_matrix.det()==s.Rational(2,27)
    assert (initial_matrix*s.Matrix([first[v]for v in initial_vars])-initial_rhs).applyfunc(s.expand)==s.zeros(4,1)
    scalar=s.expand((3*n[13]*n[2]+n[2]-3*n[9]**2).subs(first,simultaneous=True))
    assert s.expand(scalar-n[2]*(1-216*a*a*n[2])/2)==0
    cases=[]
    for label,at2,linear_solution,last_solution,free in [
        ('n2_zero',s.Integer(0),
         {n[6]:0,n[1]:-27*a*n[4]**2},
         {n[0]:18*a*n[4]*(3*n[4]*n[12]-n[3]),n[5]:36*a*n[8]**2},[a,n[3],n[4],n[8],n[12]]),
        ('n2_nonzero',1/(216*a*a),
         {n[8]:-n[12]/(6*a),n[1]:-27*a*n[4]**2-n[6]/(12*a)-n[12]/(36*a*a)},
         {n[0]:54*a*n[4]**2*n[12]-18*a*n[3]*n[4]+n[6]*n[12]/(2*a),
          n[5]:n[12]**2/a-6*n[6]*n[12]},[a,n[3],n[4],n[6],n[12]])]:
        fixed={**first,n[2]:at2}
        # Sequential substitution is essential: the four forced expressions
        # contain the branch n2 and the later linear unknowns.
        residual=[s.cancel(v.subs(first,simultaneous=True).subs(n[2],at2))for v in J]
        linear_equations=[v for v in residual if not v.has(n[0],n[5])]
        matrix,rhs=s.linear_eq_to_matrix(linear_equations,[n[1],n[6],n[8]])
        assert matrix.rank()==2
        first_columns=[0,1]if at2==0 else[0,2]
        first_minor=s.factor(matrix.extract([3,8],first_columns).det())
        assert first_minor==(s.Rational(1,3)if at2==0 else-s.Rational(1,3))
        solution=s.Matrix([linear_solution.get(v,v)for v in [n[1],n[6],n[8]]])
        assert (matrix*solution-rhs).applyfunc(s.cancel)==s.zeros(matrix.rows,1)
        second=[s.cancel(v.subs(linear_solution,simultaneous=True))for v in residual]
        matrix2,rhs2=s.linear_eq_to_matrix(second,[n[0],n[5]])
        assert matrix2.rank()==2
        second_minor=s.factor(matrix2.extract([9,18],[0,1]).det())
        assert second_minor==(-s.Rational(5,9)if at2==0 else-s.Rational(1,2))
        assert (matrix2*s.Matrix([last_solution[n[0]],last_solution[n[5]]])-rhs2).applyfunc(s.cancel)==s.zeros(matrix2.rows,1)
        # Resolve all expressions, retaining exactly the five named free parameters.
        replacements={n[2]:at2,**linear_solution,**last_solution}
        complete={}
        for j in range(15):
            value=first.get(n[j],n[j])
            value=s.cancel(value.subs(replacements,simultaneous=True))
            if n[j]in replacements:value=replacements[n[j]]
            complete[n[j]]=s.cancel(value)
        for variable in list(complete):
            complete[variable]=s.cancel(complete[variable].subs({v:complete[v]for v in n[:15]if v not in free},simultaneous=True))
        assert set().union(*(v.free_symbols for v in complete.values()))==set(free)
        assert all(s.cancel(v.subs(complete,simultaneous=True))==0 for v in J)
        affine_pf=s.factor(Pf.subs(complete,simultaneous=True))
        cases.append(dict(name=label,localization=['a=n14!=0'],parameters=list(map(str,free)),
            complete_source_substitutions={str(k):str(v)for k,v in complete.items()},
            universal_brackets=[str(s.factor(s.sympify(text).subs(complete,simultaneous=True)))for text in row['brackets']],
            actual_pfaffian=str(affine_pf),all_four_Jacobi_identically_zero=True,
            first_complete_linear_rank=2,second_complete_linear_rank=2,full_source_parameter_count=5))
        cases[-1].update(first_uniform_rank_minor=dict(rows=[3,8],columns=first_columns,value=str(first_minor)),
                         second_uniform_rank_minor=dict(rows=[9,18],columns=[0,1],value=str(second_minor)))
    A=s.Matrix(row['field_linear_matrix']);c=s.Matrix(row['field_constant'])
    B=s.Matrix(4,4,s.symbols('b0:16'));translation=s.Matrix(s.symbols('t0:4'))
    linear_orbit=s.Matrix(list(B*A-A*B)).jacobian(list(B))
    affine_orbit=s.Matrix(list(B*A-A*B)+list(B*c-A*translation)).jacobian(list(B)+list(translation))
    assert linear_orbit.rank()==10 and affine_orbit.rank()==14
    return dict(input_report=str(source_path),input_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
        entire15_parameter_nil31_long_top3_source_preserved=True,
        n14_zero_complement=dict(n11_nonzero_contradiction='n13=-1,n9=n7=0 force -2*n11^2/3=0; hence n11=0.',
                                n9_nonzero_forced_substitution={str(k):str(v)for k,v in boundary.items()},
                                n9_nonzero_exact_pfaffian=str(boundary_pf),
                                n9_nonzero_complete_linear_minor=dict(rows=boundary_rows,value=str(boundary_minor)),
                                n9_n7_zero='HighestPf identically zero after n11=0.',n9_zero_n7_nonzero_forced_substitution={str(k):str(v)for k,v in empty.items()},
                                n9_zero_n7_nonzero_complete_linear_minor=dict(rows=empty_rows,value=str(empty_minor)),
                                exact_pfaffian=str(empty_pf),irreducible_quartic_excluded=True),
        nonzero_n14_scalar_equation='n13=-1/6, n2*(1-216*a^2*n2)=0',
        nonzero_n14_initial_linear_minor=dict(rows=[0,2,13,14],variables=list(map(str,initial_vars)),value='2/27'),
        all_two_complete_families=cases,
        dimension_bound=dict(fixed_source_dimension_at_most=5,nil31_linear_orbit_dimension=10,
            long_translation_space_dimension=4,affine_cofactor_orbit_dimension=14,hyperplane_choice_dimension=4,
            projectivization_subtracts=1,total_saturated_projective_dimension_at_most=22,
            exact_affine_orbit_tangent_matrix=[[str(v)for v in affine_orbit.row(i)]for i in range(affine_orbit.rows)],
            reason='For a fixed Jordan31 linear map, longtranslation is the open c_top3!=0 in affine4space, one centralizer/translation orbit. Its linear orbit has dimension10. Both complete source parameter spaces have dimension5.'),
        generic_classification_implication='Unknown beta-dominant fullspace components have dimensions25,26,27 by the established dominant-component geometry. A union of these full22dim source closures cannot be their general support.',
        exceptional_membership_not_proved=True,
        scope='Whole H+irreduciblequartic nil31longtop3 source is two explicit5parameter families, with all empty complements retained. E membership remains separate; dimension22 excludes its generic G dominance.')


def source_magnetic_normalform(point_path):
    """Actual affine magnetic normalization and finite monomial-pencil test."""
    import importlib.util
    from sympy.polys.matrices import DomainMatrix
    point=json.loads(point_path.read_text());source_path=Path(point['source_report'])
    row=json.loads(source_path.read_text())['full_nilpotent_sources'][point['source_row']]
    h,x,y,z,w=coords=s.symbols('h x y z w');n=s.symbols('n0:36');mag=s.symbols('m0:6');a=s.symbols('a0:6')
    values={s.Symbol(k):s.Rational(v)for k,v in point['substitutions'].items()}
    indices=[3,4,5,6,7,8];oldmag=s.Matrix([values[n[i]]for i in indices])
    universal={**values,**dict(zip([n[i]for i in indices],mag))}
    P=s.zeros(5)
    for (i,j),text in zip(it.combinations(range(1,5),2),row['brackets']):
        f=s.Poly(s.sympify(text).subs(universal,simultaneous=True),x,y,z,w)
        P[i,j]=sum(c*h**(2-sum(e))*s.prod(v**k for v,k in zip((x,y,z,w),e))for e,c in f.terms());P[j,i]=-P[i,j]
    G=s.eye(5);G[1,0],G[1,2],G[1,4]=a[:3];G[3,0],G[3,2],G[3,4]=a[3:]
    def pushed(A,T):return (T*A*T.T).subs(dict(zip(coords,T.inv()*s.Matrix(coords))),simultaneous=True).applyfunc(s.expand)
    changed=pushed(P,G);delta=s.Poly(s.expand(changed[1,3]-P[1,3]),h,y,w)
    mons=[h*h,h*y,h*w,y*y,y*w,w*w]
    dv=s.Matrix([delta.coeff_monomial(m)for m in mons]);matrix=dv.jacobian(a)
    assert dv==matrix*s.Matrix(a) and matrix.rank()==5
    left=matrix.T.nullspace();assert len(left)==1
    invariant=left[0];index=next(i for i,v in enumerate(invariant)if v)
    coefficient=s.cancel((invariant.T*oldmag)[0]/invariant[index]);assert coefficient!=0
    target=s.zeros(6,1);target[index]=coefficient
    answers=next(iter(s.linsolve((matrix,target-oldmag),a)))
    free=set().union(*(v.free_symbols for v in answers));answers=[v.subs({f:0 for f in free})for v in answers]
    Gnum=G.subs(dict(zip(a,answers)));scale=1/coefficient
    diagonal=s.diag(1,scale,1,scale,1);total=diagonal*Gnum
    original=P.subs(dict(zip(mag,oldmag)))
    simplemag=[int(i==index)for i in range(6)];simple=P.subs(dict(zip(mag,simplemag)))
    assert pushed(original,total)==scale*simple
    simple_values={**values,**dict(zip([n[i]for i in indices],simplemag))}
    helper=Path('legacy/v03/code/poisson_p4_components_check.py')
    spec=importlib.util.spec_from_file_location('quartic_magnetic_ambient',helper)
    ambient=importlib.util.module_from_spec(spec);spec.loader.exec_module(ambient)
    denominator=s.ilcm(*[c.q for i,j in ambient.PAIRS for c in s.Poly(simple[i,j],*coords).coeffs()]);simple*=denominator
    sparse={(i,j):{e:int(c)for e,c in s.Poly(simple[i,j],*coords).terms()if c}for i,j in ambient.PAIRS if simple[i,j]!=0}
    D=s.Matrix(ambient.linearized_jacobi_matrix(sparse).tolist())
    def vec(A):return s.Matrix([s.Poly(A[i,j],*coords).coeff_monomial(e)for(i,j),e in ambient.BIVECTOR_BASIS])
    directions=[]
    for i,j in it.product(range(5),repeat=2):
        B=s.zeros(5);B[i,j]=1;X=B*s.Matrix(coords)
        directions.append(vec(sum([X[k]*simple.diff(coords[k])for k in range(5)],s.zeros(5))-B*simple-simple*B.T))
    M=s.Matrix(ambient.vertical_bivector_matrix().tolist()).row_join(s.Matrix.hstack(*directions,vec(simple)))
    def dm(A):return DomainMatrix.from_Matrix(A).convert_to(s.QQ)
    rank=len(dm(M).rref()[1]);assert rank==48 and D*M==s.zeros(350,M.cols)
    zero_columns=[j for j in range(150)if D[:,j]==s.zeros(350,1)]
    transverse=[]
    for j in zero_columns:
        v=s.zeros(150,1);v[j]=1
        if len(dm(M.row_join(v)).rref()[1])>rank:
            pair,exponent=ambient.BIVECTOR_BASIS[j]
            transverse.append(dict(column=j,pair=list(pair),monomial_exponents=list(exponent)))
    return dict(point_specification=str(point_path),point_sha256=hashlib.sha256(point_path.read_bytes()).hexdigest(),
                source_report=str(source_path),source_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                affine_magnetic_translation_matrix=[[str(Gnum[i,j])for j in range(5)]for i in range(5)],
                magnetic_translation_rank=5,magnetic_invariant_row=list(map(str,invariant)),
                invariant_value=str(coefficient),canonical_magnetic_index=index,canonical_magnetic_values=simplemag,
                complete_projective_transform=[[str(total[i,j])for j in range(5)]for i in range(5)],
                tensor_scale=str(scale),all_six_brackets_conjugacy_checked=True,
                canonical_source_substitutions={str(k):str(v)for k,v in simple_values.items()},
                canonical_integer_scale=str(denominator),canonical_brackets={f'{i},{j}':str(simple[i,j])for i,j in ambient.PAIRS if simple[i,j]!=0},
                exact_actual_orbit_rank=rank,zero_D_columns=zero_columns,
                source_transverse_monomial_Poisson_pencils=transverse,
                scope='Exact sourcepoint PGLnormalform and complete150monomial pencil test. Empty test excludes this sparse pencil route, not other deformations.')


def nil_short_translation_zero_c(source_path):
    """Entire c=0 support: two empty top types and three full11parameter sources."""
    n,(u,y,w),Phi,V,Q=nil_short_translation_potential()
    row=json.loads(source_path.read_text())['full_nilpotent_sources'][8]
    c0={n[18]:0,n[19]:0,n[20]:0,n[25]:0}
    # Normalize n27!=0 to n27=1,n26=0 by affine y,w scalings/shears.
    fixed={**c0,n[27]:1,n[26]:0,n[2]:0,n[17]:-n[14],
           n[13]:2*n[16]-n[23]-1,n[12]:n[24]*n[14]}
    image=s.Poly(s.expand(Q.subs(fixed,simultaneous=True)),u,y,w)
    assert s.expand(image.coeff_monomial(w**3)-n[1]+n[14]**2)==0
    assert s.expand(image.coeff_monomial(y*w*w)+n[14]*n[16]+n[15])==0
    fixed.update({n[1]:n[14]**2,n[15]:-n[14]*n[16]})
    phi=s.expand(Phi.subs(fixed,simultaneous=True));restriction=s.Poly(phi.subs(y,n[14]),w)
    assert restriction.degree()<=1
    root_phi=restriction.as_expr();linear=restriction.coeff_monomial(w)
    root_Q=s.Poly(s.expand(Q.subs(fixed,simultaneous=True).subs(y,n[14])),w)
    assert s.expand(root_Q.coeff_monomial(w*w)-linear)==0
    assert s.expand(V[1].subs(fixed,simultaneous=True).subs(y,n[14]))==0
    empty27=dict(normalization='n27=1,n26=0',forced_parameters={str(k):str(v)for k,v in fixed.items()},
                 potential_at_y_n14=str(root_phi),root_Q_w2_equals_root_Phi_w=True,
                 conclusion='Qconstant forces linearwcoefficient zero and k=Phi(n14,w); F=-Phi+k has factor y-n14.')
    # Normalize n27=0,n26!=0 to n26=3,n24=0.
    fixed={**c0,n[27]:0,n[26]:3,n[24]:0,n[2]:0,n[14]:-3*n[17]/5,
           n[13]:3*n[23]-6*n[16]+2,n[12]:0}
    image=s.Poly(s.expand(Q.subs(fixed,simultaneous=True)),u,y,w)
    assert s.expand(image.coeff_monomial(w**3)+3*n[17]**2/5)==0
    fixed.update({n[17]:0,n[14]:0})
    image=s.Poly(s.expand(Q.subs(fixed,simultaneous=True)),u,y,w)
    assert image.coeff_monomial(y*w*w)==3*n[1]
    fixed[n[1]]=0
    phi=s.Poly(s.expand(Phi.subs(fixed,simultaneous=True)),w)
    A=phi.coeff_monomial(w);D=phi.coeff_monomial(1);L=y*y+n[16]*y+n[15]
    assert s.Poly(A,y).degree()==3
    image=s.Poly(s.expand(Q.subs(fixed,simultaneous=True)),w)
    assert s.expand(image.coeff_monomial(w)+L*s.diff(A,y)-(3*y+n[23]+1)*A)==0
    assert s.expand(image.coeff_monomial(1)+L*s.diff(D,y)-(n[22]*y+n[21])*A-D)==0
    empty26=dict(normalization='n27=0,n26=3,n24=0',forced_parameters={str(k):str(v)for k,v in fixed.items()},
                 cubic_A=str(A),D=str(D),quadratic_L=str(L),
                 root_lemma='Qw=0 gives A|(L*Aprime). At every complex root rho of cubicA, valuation forces L(rho)=0; Qconstant then gives k=D(rho). Hence F=-A*w-D+k has factor y-rho.')
    # The only remaining quartic top has n26=n27=0,r=n24!=0.
    base={**c0,n[26]:0,n[27]:0}
    top=s.Poly(s.expand(Q.subs(base,simultaneous=True)),u,y,w)
    assert top.coeff_monomial(y*y*w*w)==n[24]*n[2]
    assert s.expand(top.coeff_monomial(y**3*w)+n[24]*(2*n[14]+n[17]))==0
    assert s.expand(top.coeff_monomial(y**4)-n[24]*(s.Rational(1,4)+(n[23]-n[13])/2-n[16]))==0
    high={**base,n[2]:0,n[13]:n[23]-2*n[16]+s.Rational(1,2),n[17]:-2*n[14]}
    residual=s.Poly(s.expand(Q.subs(high,simultaneous=True)),u,y,w)
    assert residual.coeff_monomial(w**3)==-2*n[14]**2
    high.update({n[14]:0,n[17]:0})
    residual=s.Poly(s.expand(Q.subs(high,simultaneous=True)),u,y,w)
    assert s.expand(residual.coeff_monomial(w*w)-n[1]*(n[23]+s.Rational(1,2)))==0
    assert s.expand(residual.coeff_monomial(y*y*w)-n[24]*n[1]-(n[16]-s.Rational(1,4))*(n[23]-2*n[16]+1))==0
    r,C,E,H,alpha=s.symbols('r C E H alpha');mag=s.symbols('m0:6')
    families=[]
    for power in [1,3]:
        b=s.Rational(1,power);d=-s.Rational(1,2);B=b-s.Rational(1,4)
        N=2*B*B/r;root=-alpha/b;A=B*y*y+E*y+H;R=C*(y-root)**power
        k=s.expand(R.subs(y,0)-H*H/(2*N));D=s.expand(k+A*A/(2*N)-R)
        assert s.Poly(D,y).coeff_monomial(1)==0 and s.Poly(D,y).coeff_monomial(y**4)==r/4
        q=-((s.Rational(1,2)-b)*E+2*b*B*root)/N
        p=-(H/2+b*E*root)/N
        full={**c0,n[26]:0,n[27]:0,n[2]:0,n[14]:0,n[17]:0,n[24]:r,n[16]:b,n[23]:d,
              n[13]:-2*b,n[15]:alpha,n[21]:p,n[22]:q,n[0]:H,n[1]:N,n[11]:-E,
              n[9]:-s.Poly(D,y).coeff_monomial(y),n[10]:p-2*s.Poly(D,y).coeff_monomial(y*y),
              n[12]:q-3*s.Poly(D,y).coeff_monomial(y**3)}
        full.update(dict(zip([n[i]for i in [3,4,5,6,7,8]],mag)))
        assert len(full)==28
        assert s.cancel(Q.subs(full,simultaneous=True)-k)==0
        assert all(s.cancel(s.sympify(v).subs(full,simultaneous=True))==0 for v in row['all_four_jacobi'])
        target=s.expand(R-N*(w+A/N)**2/2)
        actual=s.sympify(row['actual_affine_pfaffian']).subs(full,simultaneous=True)
        assert s.cancel(actual-target)==0
        families.append(dict(name='quadratic_w_power'+str(power),localization=['r!=0','C!=0'],
            parameters=list(map(str,[r,C,E,H,alpha,*mag])),
            complete_source_substitutions={str(k):str(v)for k,v in full.items()},
            universal_brackets=[str(s.factor(s.sympify(v).subs(full,simultaneous=True)))for v in row['brackets']],
            actual_pfaffian=str(target),constant_k=str(k),all_four_Jacobi_identically_zero=True,
            irreducibility='Quadratic in w with discriminant a nonzero constant times (y-root)^odd_power, nonsquare over C(y).'))
    r,L,alpha,p,q=s.symbols('r L alpha p q')
    D3=4*alpha*r;D2=24*alpha*alpha*r-2*r*L;D1=8*alpha*D2/3-4*q*L/3
    D=r*y**4/4+D3*y**3+D2*y*y+D1*y;k=-alpha*D1+p*L
    full={**c0,n[26]:0,n[27]:0,n[2]:0,n[14]:0,n[17]:0,n[24]:r,n[16]:s.Rational(1,4),n[23]:-1,
          n[13]:-1,n[15]:alpha,n[21]:p,n[22]:q,n[0]:L,n[1]:0,n[11]:0,
          n[9]:-D1,n[10]:p-2*D2,n[12]:q-3*D3}
    full.update(dict(zip([n[i]for i in [3,4,5,6,7,8]],mag)))
    assert len(full)==28 and s.cancel(Q.subs(full,simultaneous=True)-k)==0
    assert all(s.cancel(s.sympify(v).subs(full,simultaneous=True))==0 for v in row['all_four_jacobi'])
    target=s.expand(k-L*w-D)
    assert s.cancel(s.sympify(row['actual_affine_pfaffian']).subs(full,simultaneous=True)-target)==0
    families.append(dict(name='linear_w',localization=['r!=0','L!=0'],parameters=list(map(str,[r,L,alpha,p,q,*mag])),
        complete_source_substitutions={str(k):str(v)for k,v in full.items()},
        universal_brackets=[str(s.factor(s.sympify(v).subs(full,simultaneous=True)))for v in row['brackets']],
        actual_pfaffian=str(target),constant_k=str(k),all_four_Jacobi_identically_zero=True,
        irreducibility='Primitive linear polynomial in w; homogeneous leading y^4 is nonzero.'))
    return dict(input_report=str(source_path),input_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
        whole_c_zero_source_retained=True,nonzero_n27_empty=empty27,nonzero_n26_empty=empty26,
        remaining_top='n26=n27=0,n24!=0; n24=0 has highestPf=0',
        quadratic_w_completeness='Complete square F=R(y)-n1*W^2/2. Jacobi forces (alpha+b*y)*Rprime=R. Irreducibility permits only nonzero odd degree1or3, hence b=1or1/3; all even/zero cases factor.',
        linear_w_completeness='Qw=(−alpha−b*y)*Aprime+(d+1)A=0. Any nonconstantA has a root shared by −alpha−b*y and forces an Ffactor. Remaining Aconstant nonzero forces b=1/4,d=-1 and displayed fullsource.',
        all_three_complete_families=families,all_universal_Jacobi_and_Pfaffians_verified=True,
        scope='All c=0 irreduciblequartic support reduced to three complete elevenparameter sources; their ambient membership remains open.')


def source_family_generic(point_path):
    """One fullambient family-rank comparison over QQ(b), no point sampling."""
    import importlib.util
    from sympy.polys.matrices import DomainMatrix
    point=json.loads(point_path.read_text());family_path=Path(point['family_report'])
    family=json.loads(family_path.read_text())[point['family_result_key']]
    helper=Path('legacy/v03/code/poisson_p4_components_check.py')
    spec=importlib.util.spec_from_file_location('quartic_generic_ambient',helper)
    ambient=importlib.util.module_from_spec(spec);spec.loader.exec_module(ambient)
    b=s.Symbol('b');h,x,y,z,w=coords=s.symbols('h x y z w')
    parameters=[s.Symbol(v)for v in family['parameters']]
    parameter_env={str(v):v for v in parameters}
    values={s.Symbol(k):s.Rational(v)for k,v in point['family_parameter_point'].items()}
    values[b]=b
    universal=s.zeros(5)
    for (i,j),text in zip(it.combinations(range(1,5),2),family['universal_brackets']):
        polynomial=s.Poly(s.sympify(text,locals=parameter_env),x,y,z,w)
        entry=sum(c*h**(2-sum(e))*s.prod(v**a for v,a in zip((x,y,z,w),e))for e,c in polynomial.terms())
        universal[i,j]=entry;universal[j,i]=-entry
    P=universal.subs(values).applyfunc(s.expand)
    denominator=s.ilcm(*[c.q for i,j in ambient.PAIRS for c in s.Poly(P[i,j],*coords,b).coeffs()])
    P*=denominator;universal*=denominator
    maximum=max(s.Poly(P[i,j],b).degree()for i,j in ambient.PAIRS if P[i,j]!=0)
    D=s.zeros(350,150);coefficient_hashes=[]
    for degree in range(maximum+1):
        sparse={}
        for i,j in ambient.PAIRS:
            entry=s.Poly(P[i,j],b).coeff_monomial(b**degree)
            if entry!=0:sparse[i,j]={mon:int(c)for mon,c in s.Poly(entry,*coords).terms()if c}
        matrix=s.Matrix(ambient.linearized_jacobi_matrix(sparse).tolist())
        D+=b**degree*matrix
        coefficient_hashes.append(hashlib.sha256(json.dumps(matrix.tolist(),default=str).encode()).hexdigest())
    def vector(A):return s.Matrix([s.Poly(A[i,j],*coords).coeff_monomial(mon)for(i,j),mon in ambient.BIVECTOR_BASIS])
    orbit=[]
    for i,j in it.product(range(5),repeat=2):
        B=s.zeros(5);B[i,j]=1;X=B*s.Matrix(coords)
        orbit.append(vector(sum([X[k]*P.diff(coords[k])for k in range(5)],s.zeros(5))-B*P-P*B.T))
    M=s.Matrix(ambient.vertical_bivector_matrix().tolist()).row_join(s.Matrix.hstack(*orbit,vector(P)))
    actual=M.row_join(s.Matrix.hstack(*[vector(universal.diff(v).subs(values).applyfunc(s.expand))for v in parameters]))
    assert (D*actual).applyfunc(s.expand)==s.zeros(350,actual.cols)
    field=s.QQ.frac_field(b)
    def dm(A):return DomainMatrix.from_Matrix(A).convert_to(field)
    def certificate(A):
        cols=list(dm(A).rref()[1]);rows=list(dm(A[:,cols].T).rref()[1]);det=dm(A[rows,cols]).det()
        assert det
        return dict(rank=len(cols),rows=rows,columns=cols,determinant=str(det))
    dcert,mcert,fcert=certificate(D),certificate(M),certificate(actual)
    return dict(point_specification=str(point_path),point_sha256=hashlib.sha256(point_path.read_bytes()).hexdigest(),
        family_report=str(family_path),family_sha256=hashlib.sha256(family_path.read_bytes()).hexdigest(),
        helper_path=str(helper),helper_sha256=hashlib.sha256(helper.read_bytes()).hexdigest(),coefficient_field='QQ(b)',
        generic_parameter='b',other_parameter_values={str(k):str(v)for k,v in values.items()if k!=b},
        integer_scale=str(denominator),full_derivative_shape=[350,150],
        derivative_polynomial_degree=int(maximum),derivative_coefficient_matrix_sha256=coefficient_hashes,
        derivative=dcert,orbit_vertical_scale=mcert,actual_full_family=fcert,
        all_actual_family_directions_annihilated_exactly=True,
        projective_tangent_dimension=125-dcert['rank'],actual_family_projective_rank=fcert['rank']-25,
        smooth_full_source_component_proved=dcert['rank']+fcert['rank']==150,
        exceptional_parameters='Zeros/poles of displayed minors. If component criterion passes, the entire irreducible source closure retains all exceptional points; otherwise no coverage is claimed.')


def nil_short_translation_constant_c(source_path):
    """Complete c=1 irreducible-quartic source as one irreducible10parameter family."""
    n,(u,y,w),Phi,V,Q=nil_short_translation_potential()
    row=json.loads(source_path.read_text())['full_nilpotent_sources'][8]
    r,b,alpha,p,q,A,D,G,H,k=s.symbols('r b alpha p q A D G H k')
    base={n[18]:1,n[19]:0,n[20]:0,n[25]:0,n[26]:0,n[27]:0,n[23]:n[16]}
    top=s.Poly(s.expand(Q.subs(base,simultaneous=True)),u,y,w)
    assert s.diff(top.as_expr(),u)==0
    assert s.expand(top.coeff_monomial(y**4)-n[24]*(s.Rational(1,4)-(n[16]+n[13])/2))==0
    assert s.expand(top.coeff_monomial(y**3*w)+n[24]*(n[17]+2*n[14]))==0
    assert top.coeff_monomial(y*y*w*w)==n[24]*n[2]
    high={**base,n[2]:0,n[13]:s.Rational(1,2)-n[16],n[17]:-2*n[14]}
    residual=s.Poly(s.expand(Q.subs(high,simultaneous=True)),u,y,w)
    assert residual.coeff_monomial(w**3)==-2*n[14]**2
    F=-(1-b)*(b-s.Rational(1,4))/(2*r)
    E=alpha-((1-3*b)*A+q*(b-s.Rational(1,4)))/r
    phi=r*y**4/4+A*y**3+(b-s.Rational(1,4))*y*y*w+D*y*y+E*y*w+F*w*w+G*y+H*w
    planar=s.Matrix([-alpha-b*y,p+q*y+b*w+r*y*y])
    A2=s.expand(phi+planar[0]*s.diff(phi,y)+planar[1]*s.diff(phi,w)-k)
    assert s.Poly(A2,y,w).total_degree()<=2
    coefficients=s.Poly(A2,y,w)
    full={n[18]:1,n[19]:0,n[20]:0,n[25]:0,n[26]:0,n[27]:0,n[23]:b,n[16]:b,
          n[2]:0,n[13]:s.Rational(1,2)-b,n[17]:0,n[14]:0,n[24]:r,n[15]:alpha,
          n[21]:p,n[22]:q,n[12]:q-3*A,n[10]:p-2*D,n[9]:-G,n[0]:H,n[1]:2*F,n[11]:-E}
    for index,monomial in zip([3,4,5,6,7,8],[1,y,w,y*y,y*w,w*w]):
        full[n[index]]=coefficients.coeff_monomial(monomial)
    assert len(full)==28
    assert s.cancel(Phi.subs(full,simultaneous=True)-u-phi)==0
    assert s.cancel(Q.subs(full,simultaneous=True)-k)==0
    assert all(s.cancel(s.sympify(text).subs(full,simultaneous=True))==0 for text in row['all_four_jacobi'])
    brackets=[str(s.factor(s.sympify(text).subs(full,simultaneous=True)))for text in row['brackets']]
    return dict(input_report=str(source_path),input_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
        c_normalization='Nonzero constant c scaled to1 by an invertible affine w scaling.',
        localization=['r=n24!=0'],complement='r=0 gives degree<=3 potential and highestPf=0 after fullu equations.',
        parameters=list(map(str,[r,b,alpha,p,q,A,D,G,H,k])),
        reduced_necessity='Quarticcoefficients force n2=0,n13=1/2-n16,n17=-2n14; cubicw3=-2n14^2 forces n14=n17=0 on reducedsupport. Remaining two cubiccoefficients uniquely solve n1,n11; all lower equations determine A2 up to constant k.',
        potential_phi=str(phi),potential_full=str(u+phi),planar_vector_field=list(map(str,planar)),
        complete_source_substitutions={str(a):str(v)for a,v in full.items()},universal_brackets=brackets,
        all_four_Jacobi_identically_zero=True,actual_pfaffian=str(k-u-phi),
        full_source_parameter_count=10,
        scope='Complete c=1 source retained; ambient component and existingfamily recognition require separate fullambient checks.')


def nil_short_translation_potential():
    n=s.symbols('n0:36');u,y,w=s.symbols('u y w')
    Phi=u*(n[18]+n[19]*y+n[20]*w-n[25]*y*y/3)+n[0]*w+n[1]*w*w/2+n[2]*w**3/3-n[11]*y*w\
        +(n[23]-n[13])*y*y*w/2-n[14]*y*w*w+n[26]*y**3*w/3+n[27]*y*y*w*w/2\
        -n[9]*y-(n[10]-n[21])*y*y/2-(n[12]-n[22])*y**3/3+n[24]*y**4/4
    V=s.Matrix([(n[16]-n[23]-1-n[26]*y/3-n[27]*w)*u-(n[3]+n[4]*y+n[5]*w+n[6]*y*y+n[7]*w*y+n[8]*w*w),
                -n[15]-n[16]*y-n[17]*w-n[26]*y*y/3-n[27]*w*y,
                n[21]+n[22]*y+n[23]*w+n[24]*y*y+n[26]*w*y+n[27]*w*w-n[25]*u])
    Q=s.expand(sum(V[i]*s.diff(Phi,v)for i,v in enumerate((u,y,w)))+Phi)
    return n,(u,y,w),Phi,V,Q


def nil_short_translation_empty_c(source_path):
    n,(u,y,w),Phi,V,Q=nil_short_translation_potential()
    cases=[]
    # c=w. The u coefficient forces these five coefficients, then restriction
    # w=0 is an elementary polynomial translation ODE with only zero solution.
    fixed={n[20]:1,n[18]:0,n[19]:0,n[25]:0,n[26]:0,n[21]:0,n[22]:0,n[24]:0,n[16]:0}
    image=s.expand(Q.subs(fixed,simultaneous=True));assert s.diff(image,u)==0
    D=-n[9]*y-n[10]*y*y/2-n[12]*y**3/3
    assert s.expand(image.subs(w,0)-(D-n[15]*s.diff(D,y)))==0
    forced={**fixed,n[9]:0,n[10]:0,n[12]:0}
    phi=s.expand(Phi.subs(forced,simultaneous=True));assert s.cancel(phi/w).is_polynomial(u,y,w,*n)
    cases.append(dict(c='w',forced_parameters={str(k):str(v)for k,v in forced.items()},
                      restriction_ODE=str(D-n[15]*s.diff(D,y)),polynomial_ODE_forces_D_zero=True,
                      potential=str(phi),constant_k=0,Pfaffian_factor='w'))
    # c=y. Restriction y=0 is the same invertible triangular ODE in w.
    fixed={n[19]:1,n[18]:0,n[20]:0,n[25]:0,n[15]:0,n[23]:0,n[17]:0,n[26]:0,n[27]:0}
    image=s.expand(Q.subs(fixed,simultaneous=True));assert s.diff(image,u)==0
    A=n[0]*w+n[1]*w*w/2+n[2]*w**3/3
    assert s.expand(image.subs(y,0)-(A+n[21]*s.diff(A,w)))==0
    forced={**fixed,n[0]:0,n[1]:0,n[2]:0}
    phi=s.expand(Phi.subs(forced,simultaneous=True));assert s.cancel(phi/y).is_polynomial(u,y,w,*n)
    cases.append(dict(c='y',forced_parameters={str(k):str(v)for k,v in forced.items()},
                      restriction_ODE=str(A+n[21]*s.diff(A,w)),polynomial_ODE_forces_A_zero=True,
                      potential=str(phi),constant_k=0,Pfaffian_factor='y'))
    # c=y^2. The restriction at y=0 first forces n17^2=0, then n15^2=0.
    fixed={n[25]:-3,n[18]:0,n[19]:0,n[20]:0,n[2]:0,n[0]:0,n[1]:0,
           n[13]:(-2*n[16]+n[23])/3,n[14]:-n[17]/3,n[11]:-2*n[15]/3}
    image=s.expand(Q.subs(fixed,simultaneous=True));assert s.diff(image,u)==0
    restriction=s.expand(image.subs(y,0))
    expected=-(n[15]+n[17]*w)*(n[17]*w*w/3+2*n[15]*w/3-n[9])
    assert s.expand(restriction-expected)==0
    assert s.Poly(restriction,w).coeff_monomial(w**3)==-n[17]**2/3
    assert s.Poly(restriction.subs(n[17],0),w).coeff_monomial(w)==-2*n[15]**2/3
    forced={**fixed,n[17]:0,n[15]:0,n[14]:0,n[11]:0}
    phi=s.expand(Phi.subs(forced,simultaneous=True));assert s.cancel(phi/y).is_polynomial(u,y,w,*n)
    cases.append(dict(c='y^2',forced_parameters={str(k):str(v)for k,v in forced.items()},
                      restriction=str(restriction),reduced_square_consequences=['n17=0','n15=0'],
                      potential=str(phi),constant_k=0,Pfaffian_factor='y'))
    return dict(input_report=str(source_path),input_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                all_three_complete_cases=cases,all_cases_exclude_irreducible_quartic=True,
                normalizations='Affine centralizers preserving y*d_x+d_z normalize n20!=0 to c=w; n25=n20=0,n19!=0 to c=y; n25!=0,double root to c=y^2.',
                complements=['Distinct-root quadratic c: separate complete five-parameter source','Nonzero constant c: open','c=0: open'],
                scope='Entire c=w,c=y,c=y^2 supports excluded by fullJacobi and actualPfaffian, including all coefficient-zero complements.')


def nil_short_translation_distinct_c(source_path):
    """Complete n25!=0, discriminant(c)!=0 source, not its ambient component."""
    row=json.loads(source_path.read_text())['full_nilpotent_sources'][8]
    n=s.symbols('n0:36');x,y,z,w=s.symbols('x y z w');u=s.Symbol('u')
    pairs=list(it.combinations(range(4),2))
    Phi=u*(n[18]+n[19]*y+n[20]*w-n[25]*y*y/3)+n[0]*w+n[1]*w*w/2+n[2]*w**3/3-n[11]*y*w\
        +(n[23]-n[13])*y*y*w/2-n[14]*y*w*w+n[26]*y**3*w/3+n[27]*y*y*w*w/2\
        -n[9]*y-(n[10]-n[21])*y*y/2-(n[12]-n[22])*y**3/3+n[24]*y**4/4
    V=s.Matrix([(n[16]-n[23]-1-n[26]*y/3-n[27]*w)*u-(n[3]+n[4]*y+n[5]*w+n[6]*y*y+n[7]*w*y+n[8]*w*w),
                -n[15]-n[16]*y-n[17]*w-n[26]*y*y/3-n[27]*w*y,
                n[21]+n[22]*y+n[23]*w+n[24]*y*y+n[26]*w*y+n[27]*w*w-n[25]*u])
    Q=s.expand(sum(V[i]*s.diff(Phi,v)for i,v in enumerate((u,y,w)))+Phi)
    initial={n[25]:-3,n[18]:-1,n[19]:0,n[20]:0,n[2]:0,
             n[13]:(-2*n[16]+n[23])/3,n[14]:-n[17]/3,n[1]:-n[27]/3,
             n[0]:(n[16]-n[23])/3,n[11]:n[26]/9-2*n[15]/3}
    normal_Q=s.expand(Q.subs(initial,simultaneous=True))
    assert s.diff(normal_Q,u)==0
    cubic_w=s.Poly(normal_Q,w).coeff_monomial(w**3)
    for root in [-1,1]:assert s.expand(cubic_w.subs(y,root)+(n[17]+root*n[27])**2/3)==0
    reduced=s.expand(normal_Q.subs({n[17]:0,n[27]:0}))
    linear_w=s.Poly(reduced,w).coeff_monomial(w)
    for root in [-1,1]:
        b=n[15]+root*n[16]+n[26]/3
        assert s.expand(linear_w.subs(y,root)-2*b*(root-b)/3)==0
    transformations=[]
    for positive,negative in it.product([0,1],repeat=2):
        a=s.Rational(positive+negative,2);b=s.Rational(positive-negative,2)
        flag={n[17]:0,n[27]:0,n[16]:a,n[15]:b-n[26]/3}
        P=s.zeros(4)
        for (i,j),text in zip(pairs,row['brackets']):
            P[i,j]=s.sympify(text).subs(initial,simultaneous=True).subs(flag,simultaneous=True);P[j,i]=-P[i,j]
        G=s.eye(4);G[0,3]=(a+n[23])/3;G[2,3]=-n[26]/3
        changed=(G*P*G.T).subs(dict(zip((x,y,z,w),G.inv()*s.Matrix((x,y,z,w)))),simultaneous=True).applyfunc(s.expand)
        hc=(s.Integer(1),x,y,z,w)
        mons=[hc[i]*hc[j]for i,j in it.combinations_with_replacement(range(5),2)]
        coefficients=[s.Poly(changed[i,j],x,y,z,w).coeff_monomial(mon)for i,j in pairs for mon in mons]
        new_parameters={n[j]:s.factor(coefficients[index])for j,index in enumerate(row['free_coefficient_indices'])}
        assert new_parameters[n[26]]==0 and new_parameters[n[23]]==-a
        assert new_parameters[n[16]]==a and new_parameters[n[15]]==b
        for (i,j),text in zip(pairs,row['brackets']):
            assert s.expand(changed[i,j]-s.sympify(text).subs(new_parameters,simultaneous=True))==0
        canonical={n[25]:-3,n[18]:-1,n[19]:0,n[20]:0,n[2]:0,n[17]:0,n[27]:0,
                   n[26]:0,n[23]:-a,n[16]:a,n[15]:b,n[0]:2*a/3,n[1]:0,n[14]:0,
                   n[13]:-a,n[11]:-2*b/3}
        A=(s.Integer(positive-negative)*y+positive+negative)/3
        transformed_phi=s.expand(Phi.subs(canonical,simultaneous=True))
        D=s.expand(transformed_phi-(y*y-1)*u-A*w)
        assert not D.has(u,w)
        Qcan=s.expand(Q.subs(canonical,simultaneous=True))
        inactive=[]
        for root,active in [(1,positive),(-1,negative)]:
            if not active:
                assert s.expand(Qcan.subs(y,root)-D.subs(y,root))==0
                inactive.append(root)
        transformations.append(dict(flags=[positive,negative],matrix=[[str(G[i,j])for j in range(4)]for i in range(4)],
            all_six_brackets_exactly_matched=True,canonical_parameters={str(k):str(v)for k,v in canonical.items()},
            potential_w_coefficient=str(A),inactive_roots_forcing_Pf_factor=inactive,
            transformed_parameters={str(k):str(v)for k,v in new_parameters.items()}))
    canonical={n[25]:-3,n[18]:-1,n[19]:0,n[20]:0,n[2]:0,n[17]:0,n[27]:0,
               n[26]:0,n[23]:-1,n[16]:1,n[15]:0,n[0]:s.Rational(2,3),n[1]:0,n[14]:0,n[13]:-1,n[11]:0}
    unknowns=[n[j]for j in [3,4,5,6,7,8,9,10,12,21,22,24]]
    Qcan=s.Poly(s.expand(Q.subs(canonical,simultaneous=True)),u,y,w)
    equations=[coefficient for exponent,coefficient in Qcan.terms()if any(exponent)]
    matrix,rhs=s.linear_eq_to_matrix(equations,unknowns)
    assert matrix.rank()==7
    D1,D2,D3,D4,d0=s.symbols('D1 D2 D3 D4 d0')
    solution={n[3]:-D2-D4/3,n[4]:-2*D3,n[5]:0,n[6]:-3*D4,n[7]:0,n[8]:0,n[9]:-D1,
              n[10]:d0-2*D2,n[12]:0,n[21]:d0,n[22]:3*D3,n[24]:4*D4}
    assert matrix*s.Matrix([solution[v]for v in unknowns])==rhs
    assert s.Matrix([solution[v]for v in unknowns]).jacobian([D1,D2,D3,D4,d0]).rank()==5
    full={**canonical,**solution}
    assert all(s.expand(s.sympify(f).subs(full,simultaneous=True))==0 for f in row['all_four_jacobi'])
    phi=s.expand(Phi.subs(full,simultaneous=True));k=2*d0/3-D2-D4/3
    affine_pf=s.expand(s.sympify(row['actual_affine_pfaffian']).subs(full,simultaneous=True).subs({x:u+y*z},simultaneous=True))
    assert s.expand(affine_pf+phi-k)==0
    return dict(input_report=str(source_path),input_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                hypotheses=['n25!=0','discriminant of c(y) nonzero'],
                c_normalization='Affine y shift/scaling, x transformed by x_new=lambda*x-lambda*beta*z, and w scaling normalize c=y^2-1 while preserving Z=y*d_x+d_z. All maps invertible over C.',
                complements=['n25=0','double root c(y) retained separately'],
                c_normalized_forced_parameters={str(k):str(v)for k,v in initial.items()},
                all_u_coefficients_zero=True,root_w_cubic_squares_force_n17_n27_zero=True,
                root_w_linear_equation='(2/3)b(r)*(r-b(r))=0 for r=+/-1',
                all_four_root_flags=transformations,
                three_inactive_flag_cases_have_reducible_Pfaffian=True,
                surviving_flags=[1,1],canonical_five_parameter_source={str(k):str(v)for k,v in full.items()},
                complete_linear_Jacobi_rank=7,free_source_parameters=list(map(str,[D1,D2,D3,D4,d0])),
                all_four_Jacobi_identically_zero=True,potential=str(phi),constant_k=str(k),actual_pfaffian=str(affine_pf),
                full_source_dimension_before_normalization=10,
                scope='Whole n25-nonzero distinct-root source reduced to one irreducible five-parameter canonical source and invertible affine normalizers. Ambient membership remains open.')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--nil-rank1-geometry', action='store_true', help='Recognize the quartic divisor, not the full Poisson fiber')
    ap.add_argument('--nil-rank1-support', action='store_true', help='Whole nil211,c0 irreducible-quartic source to MINIMAL_M')
    ap.add_argument('--nil-regular-support',type=Path,help='Whole regular-nilpotent scalar support reduction from frozen10source report')
    ap.add_argument('--nil-regular-component',type=Path,help='One-normal fullambient support test for nil4 quartic source')
    ap.add_argument('--nil-regular-second',type=Path,help='Fullthree-normal quadratic support test for the other nil4quartic source')
    ap.add_argument('--nil-regular-graph',type=Path,help='Exact T39/T40 wholegraph recognition at the second nil4quartic source')
    ap.add_argument('--nil-short-translation-geometry',type=Path,help='Entire nil211translation1 potential reduction')
    ap.add_argument('--source-point',type=Path,help='Fullambient exact D/orbit/kernel at a specified classified-source point')
    ap.add_argument('--source-graph',type=Path,help='Exact T39/T40graph at a structurally specified sourcepoint and rank4evaluation')
    ap.add_argument('--nil-short-translation-distinct-c',type=Path,help='Complete distinct-root quadratic-c branch of nil211translation1')
    ap.add_argument('--nil-short-translation-empty-c',type=Path,help='Exclude entire w,y,y^2 cofactor-potential branches')
    ap.add_argument('--nil-short-translation-constant-c',type=Path,help='Complete10parameter c=1 irreducible-quartic source')
    ap.add_argument('--source-family-generic',type=Path,help='Fullambient actualfamily comparison over QQ(b)')
    ap.add_argument('--source-point-quadratic',type=Path,help='Fulltrue normal quadratic support transverse to actualPGL and complete family')
    ap.add_argument('--nil-short-translation-zero-c',type=Path,help='Complete c=0 support to three11parameter sources')
    ap.add_argument('--source-point-twist',type=Path,help='Complete linear search for commuting cofactor twists and actualPoisson pencils')
    ap.add_argument('--source-magnetic-normalform',type=Path,help='Actual6magnetic affine normalform and complete150monomial compatiblepencil test')
    ap.add_argument('--nil31-long-support',type=Path,help='Complete nil31longtop3 quartic source and saturateddimension bound')
    ap.add_argument('--nil22-long-dimension',type=Path,help='Whole nil22top2 source dimension via lifted QQ leading monomials')
    ap.add_argument('--nil22-long-core',type=Path,help='Complete eightvariable nonlinear core of nil22top2 source')
    ap.add_argument('--nil22-long-critical-cores',type=Path,help='Complete nil22top2 core strata n14=-1,-1/2,-1/3')
    ap.add_argument('--nil22-long-low-rank-core',type=Path,help='Complete nil22top2 n14=-1/2 low-rank-A core')
    ap.add_argument('--nil22-long-support',type=Path,help='Uniform whole nil22top2 source bound and generic non-dominance')
    ap.add_argument('--sources', action='store_true', help='All24 full semisimple/mixed modular sources')
    ap.add_argument('--nil-sources', action='store_true', help='All10 full purely nilpotent modular sources')
    ap.add_argument('--support',type=Path,help='Verify finite support factorizations from a frozen24source report')
    args = ap.parse_args()
    source = Path(__file__).read_bytes()
    report = dict(task='E_reduced_nonnormal', coefficient_ring='QQ, exact symbolic polynomials',
                  source_sha256=hashlib.sha256(source).hexdigest(), source_bytes=source.decode(),
                  scope='Finite modular cofactor entrance only; residual full Jacobi sources and coverage remain open.')
    if args.nil22_long_support:
        result=nil22_long_support(args.nil22_long_support)
        report.update(nil22_long_support=result,all_passed=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,fixed_source_dimension_upper=9,
              saturated_projective_dimension_upper=24,excludes_unknown_G25_to27=True)))
        return
    if args.nil22_long_low_rank_core:
        result=nil22_long_dimension(args.nil22_long_low_rank_core,core_only=True,fixed_n14=-s.Rational(1,2),a_rank_one=True)
        report.update(nil22_long_low_rank_core=result,all_passed=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,low_rank_core_dimension=result['core_dimension_upper'])))
        return
    if args.nil22_long_critical_cores:
        result=[nil22_long_dimension(args.nil22_long_critical_cores,core_only=True,fixed_n14=v)
                for v in [-s.Integer(1),-s.Rational(1,2),-s.Rational(1,3)]]
        report.update(nil22_long_critical_cores=result,all_passed=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,critical_core_dimensions=[r['core_dimension_upper']for r in result])))
        return
    if args.nil22_long_core:
        result=nil22_long_dimension(args.nil22_long_core,core_only=True)
        report.update(nil22_long_core=result,all_passed=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,core_dimension_upper=result['core_dimension_upper'],
              full_source_dimension_not_yet_proved=True)))
        return
    if args.nil22_long_dimension:
        result=nil22_long_dimension(args.nil22_long_dimension)
        report.update(nil22_long_dimension=result,all_passed=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,fixed_dimension_upper=result['fixed_source_dimension_upper'],
              saturated_dimension_upper=result['saturated_projective_dimension_upper'],
              excludes_unknown_G25_to27=result['excludes_unknown_G25_to27'])))
        return
    if args.nil31_long_support:
        result=nil31_long_support(args.nil31_long_support)
        report.update(nil31_long_support=result,all_passed=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,complete_sources=2,fixed_dimension=5,saturated_dimension_upper=22)))
        return
    if args.source_magnetic_normalform:
        result=source_magnetic_normalform(args.source_magnetic_normalform)
        report.update(source_magnetic_normalform=result,all_passed=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,magnetic_index=result['canonical_magnetic_index'],
              actual_orbit_rank=result['exact_actual_orbit_rank'],transverse_monomial_pencils=len(result['source_transverse_monomial_Poisson_pencils']))))
        return
    if args.source_point_twist:
        point=json.loads(args.source_point_twist.read_text())
        result=nil_regular_component(Path(point['source_report']),point=point,twist=True)
        report.update(source_point_twist=result,all_passed=True,
                      scope='Complete commuting-linear-field twist search; no component recognition from a negative search.')
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        info=result['commuting_twist']
        print(json.dumps(dict(all_passed=True,Y_dimension=info['commuting_compatible_Y_dimension'],
             H_removing_transverse_twist=info['H_removing_transverse_twist_exists'])))
        return
    if args.nil_short_translation_zero_c:
        result=nil_short_translation_zero_c(args.nil_short_translation_zero_c)
        report.update(nil_short_translation_zero_c=result,all_passed=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,full_sources=[v['name']for v in result['all_three_complete_families']])))
        return
    if args.source_point_quadratic:
        point=json.loads(args.source_point_quadratic.read_text())
        result=nil_regular_component(Path(point['source_report']),point=point,quadratic=True)
        report.update(source_point_quadratic=result,all_passed=True,
                      scope='Entire true-normal quadratic initial support; only a zero-dimensional result proves source componenthood.')
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,true_normal_count=result['true_normal_count'],
             zero_dimensional=result['quadratic_initial_support_zero_dimensional'],
             quadratic_generators=result['all_quadratic_residual_generators'])))
        return
    if args.source_family_generic:
        result=source_family_generic(args.source_family_generic)
        report.update(source_family_generic=result,all_passed=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,tangent_dimension=result['projective_tangent_dimension'],
                             actual_family_rank=result['actual_family_projective_rank'],
                             smooth_component=result['smooth_full_source_component_proved'])))
        return
    if args.nil_short_translation_constant_c:
        result=nil_short_translation_constant_c(args.nil_short_translation_constant_c)
        report.update(nil_short_translation_constant_c=result,all_passed=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,whole_source_parameters=10)))
        return
    if args.nil_short_translation_empty_c:
        result=nil_short_translation_empty_c(args.nil_short_translation_empty_c)
        report.update(nil_short_translation_empty_c=result,all_passed=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,entire_empty_branches=['w','y','y^2'])))
        return
    if args.nil_short_translation_distinct_c:
        result=nil_short_translation_distinct_c(args.nil_short_translation_distinct_c)
        report.update(nil_short_translation_distinct_c=result,all_passed=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,canonical_source_parameters=5,ambient_membership='open')))
        return
    if args.source_graph:
        point=json.loads(args.source_graph.read_text())
        result=nil_regular_graph(Path(point['source_report']),point=point)
        report.update(source_graph=result,all_passed=True,
                      scope='Exact pointgraph recognition; whole-source membership still requires the explicit source-containment argument.')
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,membership={k:v['whole_graph_membership']for k,v in result['graph_membership'].items()})))
        return
    if args.source_point:
        point=json.loads(args.source_point.read_text())
        result=nil_regular_component(Path(point['source_report']),point=point)
        report.update(source_point=result,all_passed=True,
                      scope='Exact fullambient point D/orbit/kernel; allsource classification and distinctness are separate.')
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,tangent_dimension=result['projective_tangent_dimension'],
                             orbit_dimension=result['projective_orbit_dimension'],
                             smooth_orbit_component=result['smooth_orbit_component_proved'],
                             actual_family_rank=result.get('actual_family_projective_rank'),
                             smooth_full_source_component=result.get('smooth_full_source_component_proved'))))
        return
    if args.nil_short_translation_geometry:
        result=nil_short_translation_geometry(args.nil_short_translation_geometry)
        report.update(nil_short_translation_geometry=result,all_passed=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,whole_source_parameters=28)))
        return
    if args.nil_regular_graph:
        result=nil_regular_graph(args.nil_regular_graph)
        report.update(nil_regular_graph=result,all_passed=True,
                      scope='Exact pointgraph recognition; whole-source coverage additionally uses tangent25 and irreducibility.')
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,membership={k:v['whole_graph_membership']for k,v in result['graph_membership'].items()})))
        return
    if args.nil_regular_second:
        result=nil_regular_component(args.nil_regular_second,second=True)
        report.update(nil_regular_second_component=result,all_passed=True,
                      scope='Exact fullthree-normal quadratic initial support; remaining directions are retained if not zero dimensional.')
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,support_component_proved=result['support_component_proved'],
                             quadratic_generators=result['all_quadratic_residual_generators'])))
        return
    if args.nil_regular_component:
        result=nil_regular_component(args.nil_regular_component)
        report.update(nil_regular_component=result,all_passed=True,
                      scope='Exact local support test at one fullsource point; distinctness and global source recognition are separate.')
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,support_component_proved=result['support_component_proved'],
                              first_nonzero_obstruction_order=result.get('first_nonzero_obstruction_order'))))
        return
    if args.nil_regular_support:
        report.update(nil_regular_support=nil_regular_support(args.nil_regular_support),all_passed=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,type2_excluded=True,type1_remaining_sources=2)))
        return
    if args.nil_rank1_support:
        report.update(nil_rank_one_support=nil_rank_one_quartic_support(),all_passed=True,
                      scope='Whole nil211,c0 irreducible-quartic source coverage, conditional on known MINIMAL_M unique ambient component.')
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,whole_source_parameters=11)))
        return
    if args.nil_rank1_geometry:
        report.update(nil_rank_one_divisor_geometry=nil_rank_one_quartic_geometry(),all_passed=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(dict(all_passed=True,rank=13,scope='Divisor recognition only')))
        return
    if args.sources:
        report.update(full_semisimple_sources=semisimple_sources(), all_passed=True)
        args.output.write_text(json.dumps(report, indent=2)+'\n')
        return
    if args.nil_sources:
        report.update(full_nilpotent_sources=nilpotent_sources(args.output),all_passed=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        return
    if args.support:
        report.update(semisimple_support=semisimple_support(args.support),all_passed=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        return
    report.update(generic_walls=generic_walls(), resonance_table=resonance_table(), all_passed=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:report['resonance_table'][k] for k in
                      ('wall_count','low_degree_wall_count','semisimple_spectrum_count','affine_semisimple_type_count')}))


if __name__ == '__main__':
    main()

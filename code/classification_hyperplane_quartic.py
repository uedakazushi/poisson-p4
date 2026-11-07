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


def nil_regular_component(source_path, second=False, point=None):
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
        return dict(input_report=str(source_path),input_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
            helper_path=str(helper),helper_sha256=hashlib.sha256(helper.read_bytes()).hexdigest(),
            point=point,full150variables_and350rows=True,integer_scale=str(denominator),
            witness_brackets={f'{i},{j}':str(P[i,j])for i,j in ambient.PAIRS if P[i,j]!=0},
            derivative=dcert,orbit_vertical_scale=mcert,whole_kernel=kcert,
            full_exact_kernel_sparse_columns=[encode(K[:,j])for j in range(K.cols)],
            all_exact_D_kernel_products_zero=True,
            projective_orbit_dimension=mcert['rank']-25,projective_tangent_dimension=125-dcert['rank'],
            smooth_orbit_component_proved=dcert['rank']+mcert['rank']==150,
            matrix_sha256=dict(D=digest(D),orbit_vertical_scale=digest(M),kernel=digest(K)))
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
    ap.add_argument('--sources', action='store_true', help='All24 full semisimple/mixed modular sources')
    ap.add_argument('--nil-sources', action='store_true', help='All10 full purely nilpotent modular sources')
    ap.add_argument('--support',type=Path,help='Verify finite support factorizations from a frozen24source report')
    args = ap.parse_args()
    source = Path(__file__).read_bytes()
    report = dict(task='E_reduced_nonnormal', coefficient_ring='QQ, exact symbolic polynomials',
                  source_sha256=hashlib.sha256(source).hexdigest(), source_bytes=source.decode(),
                  scope='Finite modular cofactor entrance only; residual full Jacobi sources and coverage remain open.')
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
                             smooth_orbit_component=result['smooth_orbit_component_proved'])))
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

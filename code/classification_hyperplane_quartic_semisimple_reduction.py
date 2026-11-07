#!/usr/bin/env python3
"""Exact structural reductions of the two large purely semisimple sources.

Reuses independent reviewer arithmetic, not a producer/CAS implementation.
The type-20 constant-rank normalization is by GL2 x GL2 preserving its field.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools as it
import json
from pathlib import Path

from classification_quadratic_cofactor_review import Ring, PAIRS


def substitute(R,f,name,value):
    pos = R.names.index(name)
    out = {}
    for e,c in f.items():
        mon = list(e)
        power = mon[pos]
        mon[pos] = 0
        term = {tuple(mon):c}
        for _ in range(power):
            term = R.mul(term,value)
        out = R.add(out,term)
    return out


def eq(R,a,b):
    assert not R.add(a,R.scale(b,-1))


def make_jac(R,entries):
    P = [[{} for j in range(4)] for i in range(4)]
    for (i,j),f in zip(PAIRS,entries):
        P[i][j],P[j][i] = f,R.scale(f,-1)
    return [R.add(*(R.mul(P[a][t],R.derivative(P[b][d],t))
                    for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4)))
            for i,j,k in it.combinations(range(4),3)]


def pfaffian(R,p):
    return R.add(R.mul(p[0],p[5]),R.scale(R.mul(p[1],p[4]),-1),R.mul(p[2],p[3]))


def coefficients(R,f,degree=None):
    out = {}
    for e,c in f.items():
        if degree is None or sum(e[-4:])==degree:
            out.setdefault(e[-4:],{})[e[:-4]+(0,)*4] = c
    return list(out.values())


def normalize(f):
    if not f:
        return ()
    scale = f[min(f)]
    return tuple(sorted((e,c/scale) for e,c in f.items()))


def type23(row):
    R = Ring(19)
    x,y,z,w = [R.variable(t) for t in ('x','y','z','w')]
    n = [R.variable(f'n{i}') for i in range(19)]
    # The complete divergence-normal source, verified separately in all90.
    p = [R.parse(s) for s in row['brackets']]
    u = R.mul(x,w)
    g = R.add(p[3],R.scale(R.mul(n[9],u),-1))
    delta = R.add(n[14],n[18])
    k = R.add(R.constant(1),R.scale(delta,-1))
    f = R.add(p[2],R.scale(R.mul(k,u),-1))
    # Divide x from the first two brackets by exponent bookkeeping.
    def divide_x(poly):
        ans = {}
        for e,c in poly.items():
            assert e[19]>0
            v = list(e)
            v[19] -= 1
            ans[tuple(v)] = c
        return ans
    a,b = divide_x(p[0]),divide_x(p[1])
    Vg = R.add(R.mul(a,R.derivative(g,1)),R.mul(b,R.derivative(g,2)))
    I = R.add(Vg,R.scale(R.mul(delta,g),-1),R.mul(n[9],f))
    II = R.mul(n[9],R.add(R.constant(1),R.scale(delta,-2)))
    eq(R,p[4],R.mul(w,R.add(a,R.scale(R.derivative(g,2),-1))))
    eq(R,p[5],R.mul(w,R.add(b,R.derivative(g,1))))
    fg = R.mul(f,g)
    predicted = [R.mul(x,R.add(I,R.mul(II,u))),
                 R.scale(R.add(R.derivative(fg,2),R.mul(u,R.derivative(I,2))),-1),
                 R.add(R.derivative(fg,1),R.mul(u,R.derivative(I,1))),
                 R.scale(R.mul(w,R.add(I,R.mul(II,u))),-1)]
    actual = make_jac(R,p)
    assert actual==predicted==[R.parse(s) for s in row['all_four_jacobi']]
    desired_coeffs = [v for poly in [I,II,R.derivative(fg,1),R.derivative(fg,2)] for v in coefficients(R,poly)]
    actual_coeffs = [v for poly in actual for v in coefficients(R,poly)]
    assert {normalize(v) for v in desired_coeffs}=={normalize(v) for v in actual_coeffs}
    pf = pfaffian(R,p)
    eq(R,pf,R.add(fg,R.mul(u,g),R.mul(k,n[9],u,u),R.mul(u,I)))
    return dict(type=23,all_four_Jacobi_structural_identities=True,
                exact_Jacobi_coefficient_ideal_equality=True,
                f=R.encode(f),g=R.encode(g),a=R.encode(a),b=R.encode(b),
                delta=R.encode(delta),I=R.encode(I),II=R.encode(II),
                pfaffian_identity='Pf = f*g + u*g + (1-delta)*l*u^2 + u*I',
                branches=['g nonconstant => f=0 => u factor',
                          'g=0 => Pf=(1-delta)*l*u^2',
                          'g=gamma!=0,l=0 => delta=0,Pf=gamma*(phi+u); homogeneous quartic has h^2 factor',
                          'g=gamma!=0,l!=0 => delta=1/2,l*phi=gamma/2,Pf=l/2*(u+gamma/l)^2'])


def solve_linear(n,rows,rhs):
    a = [{k:Q(v) for k,v in row.items()} | ({n:Q(b)} if b else {}) for row,b in zip(rows,rhs)]
    pivots = []
    for j in range(n):
        r = len(pivots)
        k = next((k for k in range(r,len(a)) if a[k].get(j)),None)
        if k is None:
            continue
        a[r],a[k] = a[k],a[r]
        scale = a[r][j]
        a[r] = {k:v/scale for k,v in a[r].items()}
        for k in range(r+1,len(a)):
            scale = a[k].get(j)
            if scale:
                for t,v in a[r].items():
                    value = a[k].get(t,Q(0))-scale*v
                    if value:
                        a[k][t] = value
                    else:
                        a[k].pop(t,None)
        pivots.append(j)
    assert all(not r for r in a[len(pivots):])
    free = [j for j in range(n) if j not in pivots]
    values = [[Q(0)]*(n+1) for _ in range(n)]
    for j in free:
        values[j][j+1] = Q(1)
    for r,j in reversed(list(enumerate(pivots))):
        values[j][0] = a[r].get(n,Q(0))
        for k,v in a[r].items():
            if k not in (j,n):
                values[j] = [x-v*y for x,y in zip(values[j],values[k])]
    return free,values


def type20(row):
    R = Ring(19)
    p = [R.parse(s) for s in row['brackets']]
    jac = make_jac(R,p)
    # Rank of constant crossing matrix C, normalized by GL2 x GL2.
    branches = []
    for rank,fixed in [(1,{0:1,3:0,6:0,11:0}),(2,{0:1,3:0,6:0,11:1})]:
        pj = jac
        for idx,value in fixed.items():
            pj = [substitute(R,f,f'n{idx}',R.constant(value)) for f in pj]
        equations = [v for f in pj for v in coefficients(R,f,degree=1)]
        rows,rhs = [],[]
        for eqn in equations:
            linear = {}
            const = Q(0)
            for exp,value in eqn.items():
                assert sum(exp)<=1
                if sum(exp)==0:
                    const += value
                else:
                    idx = exp.index(1)
                    assert idx<19
                    linear[idx] = value
            rows.append(linear)
            rhs.append(-const)
        for idx,value in fixed.items():
            rows.append({idx:1})
            rhs.append(value)
        free,values = solve_linear(19,rows,rhs)
        substitutions = [R.coefficient(v) for v in values]
        source = p
        for idx,value in enumerate(substitutions):
            if idx not in free:
                source = [substitute(R,f,f'n{idx}',value) for f in source]
        residual = make_jac(R,source)
        assert all(not coefficients(R,f,degree=1) for f in residual)
        pf = pfaffian(R,source)
        residual_coeffs = {normalize(v) for f in residual for v in coefficients(R,f)}
        branches.append(dict(constant_matrix_rank=rank,normalized_constant_entries=fixed,
                             free_indices=free,substitutions=[R.encode(v) for v in substitutions],
                             brackets=[R.encode(f) for f in source],
                             residual_Jacobi_coefficients=[[[list(e),str(c)] for e,c in v] for v in sorted(residual_coeffs)],
                             pfaffian=R.encode(pf)))
    return dict(type=20,constant_rank_zero='Pure quadratic; reduced irreducible quartic excluded by closed rational form lemma.',
                nonzero_constant_rank_branches=branches)


def followup_type20(reduction):
    R = Ring(19)
    def decode(poly):
        return {tuple(e):Q(c) for e,c in poly}
    rank1,rank2 = reduction['nonzero_constant_rank_branches']
    p = [decode(f) for f in rank1['brackets']]
    factor_witnesses = []
    for substitutions,divisor in [({16:0,4:0},'w'),({7:0,8:0},'y')]:
        source = p
        for idx,value in substitutions.items():
            source = [substitute(R,f,f'n{idx}',R.constant(value)) for f in source]
        pf = pfaffian(R,source)
        assert all(e[R.names.index(divisor)]>0 for e in pf)
        factor_witnesses.append(dict(substitutions=substitutions,divisor=divisor,all_coefficients_checked=True))
    # A whole-source check of one irreducible actual quartic, illustrating why
    # rank2 cannot be discarded.  The source has SNC top after X=x-y.
    source = [decode(f) for f in rank2['brackets']]
    example_values = {idx:0 for idx in rank2['free_indices']}
    example_values.update({5:2,15:1,17:1})
    for idx,value in example_values.items():
        source = [substitute(R,f,f'n{idx}',R.constant(value)) for f in source]
    assert not any(make_jac(R,source))
    x,y,z,w = [R.variable(c) for c in ('x','y','z','w')]
    expected = [R.add(R.mul(x,y),R.scale(R.mul(y,y),-1)),R.add(R.constant(1),R.mul(x,z)),
                R.add(R.scale(R.mul(x,w),-1),R.scale(R.mul(y,w),2)),
                R.mul(y,z),R.add(R.constant(1),R.mul(y,w)),R.mul(z,w)]
    assert source==expected
    pf = pfaffian(R,source)
    eq(R,pf,R.scale(R.add(R.mul(R.add(R.constant(1),R.mul(x,z)),R.add(R.constant(1),R.mul(y,w))),
                           R.scale(R.mul(y,y,z,w),-1)),-1))
    changed = [source[0],R.add(source[1],R.scale(source[3],-1)),
               R.add(source[2],R.scale(source[4],-1)),source[3],source[4],source[5]]
    changed = [substitute(R,f,'x',R.add(x,y)) for f in changed]
    tops = [{e:c for e,c in f.items() if sum(e[-4:])==2} for f in changed]
    expected_tops = [R.mul(x,y),R.mul(x,z),R.scale(R.mul(x,w),-1),R.mul(y,z),R.mul(y,w),R.mul(z,w)]
    assert tops==expected_tops
    assert not any(make_jac(R,changed))
    eq(R,pfaffian(R,tops),R.scale(R.mul(x,y,z,w),-1))
    return dict(rank1_factor_witnesses=factor_witnesses,
                rank1_branch_exhaustion='n7*n16=0; n7!=0=>n16=n4=0; n16!=0=>n7=n8=0; n7=n16=0=>n4*n8=0.',
                rank2_verified_example=dict(free_parameter_values=example_values,
                                           brackets=[R.encode(f) for f in source],all_four_Jacobi_zero=True,
                                           actual_pfaffian=R.encode(pf),
                                           irreducibility_argument='Homogeneous F=h^4+h^2(xz+yw)+y(x-y)zw is primitive linear in z: neither h nor h^2+yw divides its z coefficient.',
                                           coordinate_change='X=x-y,Y=y,Z=z,W=w',
                                           transformed_full_brackets=[R.encode(f) for f in changed],
                                           transformed_top_six_logcanonical_coefficients=[1,1,-1,1,1,1],
                                           transformed_top_pfaffian='-X*Y*Z*W',
                                           containment='E_poly: complete SNC-top quadratic plus all Poisson lower terms.'),
                all_rank1_points_excluded=True,rank2_full_support_open=True)


def complete_type20(reduction):
    R = Ring(19)
    p = [{tuple(e):Q(c) for e,c in f} for f in reduction['nonzero_constant_rank_branches'][1]['brackets']]
    a,b,c,d,e,f,g,h,t = [R.variable(f'n{i}') for i in (4,5,7,8,12,13,14,15,17)]
    one = R.constant(1)
    def plus(*v):
        return R.add(*v)
    def neg(v):
        return R.scale(v,-1)
    def times(*v):
        return R.mul(*v)
    def scalar(v,q):
        return R.scale(v,Q(q))
    def put(poly,assignments):
        for idx,value in assignments.items():
            poly = substitute(R,poly,f'n{idx}',R.constant(value) if isinstance(value,(int,Q)) else value)
        return poly
    def changed_x(source,coefficient):
        out = [source[0],plus(source[1],times(coefficient,source[3])),
               plus(source[2],times(coefficient,source[4])),source[3],source[4],source[5]]
        return [substitute(R,v,'x',plus(R.variable('x'),neg(times(coefficient,R.variable('y'))))) for v in out]
    def coefficient_ideal(source):
        return {normalize(v) for j in make_jac(R,source) for v in coefficients(R,j)}
    # All nine equations on the complete r=t*z*w source.
    semip = [put(v,{16:0,18:0}) for v in p]
    E = [plus(times(plus(one,neg(h)),plus(one,scalar(h,-2))),scalar(times(b,f),2)),
         plus(times(f,plus(one,t,neg(h))),neg(times(b,d))),
         plus(neg(times(f,g)),times(e,plus(one,neg(h))),neg(times(b,c)),neg(times(a,d))),
         plus(times(g,plus(t,neg(e))),neg(times(a,c))),
         plus(times(g,plus(R.constant(3),neg(t),scalar(h,-4))),times(b,plus(one,neg(t),scalar(e,-2))),scalar(times(a,f),-2)),
         plus(scalar(times(g,g),2),times(a,plus(R.constant(-1),scalar(t,2),scalar(e,2)))),
         plus(times(f,plus(one,neg(t),scalar(e,-4))),scalar(times(d,g),2),times(c,plus(R.constant(-1),t,scalar(h,2)))),
         plus(scalar(times(f,f),2),times(d,plus(one,scalar(t,-2),scalar(h,-2)))),
         plus(times(e,plus(scalar(e,2),R.constant(-1))),scalar(times(c,g),-2))]
    assert coefficient_ideal(semip)=={normalize(v) for v in E}
    eq(R,plus(times(c,E[5]),times(plus(R.constant(-1),scalar(t,2),scalar(e,2)),E[3]),times(g,E[8])),
       times(g,t,plus(scalar(t,2),R.constant(-1))))
    eq(R,plus(times(f,E[0]),neg(times(b,E[7])),neg(times(plus(one,scalar(t,-2),scalar(h,-2)),E[1]))),
       times(f,t,plus(one,scalar(t,2))))

    # Every actual point in these seven full families has an SNC initial
    # tensor. The listed scalar equations are retained, never inverted away.
    families = [
        ('generic_e0_h1_c0',{16:0,18:0,13:0,14:0,4:0,8:0,12:0,15:1,7:0},[0,1,0,0],None,'t!=0,+/-1/2; (t-1)*b=0'),
        ('generic_e0_h1_b0',{16:0,18:0,13:0,14:0,4:0,8:0,12:0,15:1,5:0},[0,-1,0,0],None,'t!=0,+/-1/2; (t+1)*c=0'),
        ('generic_e0_hhalf',{16:0,18:0,13:0,14:0,4:0,8:0,12:0,15:Q(1,2),7:0},[0,1,0,0],None,'t!=0,+/-1/2; (t-1)*b=0'),
        ('generic_ehalf_h1',{16:0,18:0,13:0,14:0,4:0,8:0,12:Q(1,2),15:1,5:0},[0,-1,0,0],None,'t!=0,+/-1/2; (t+1)*c=0'),
        ('special_e0_h1',{16:0,18:0,17:Q(1,2),13:0,14:0,7:0,8:0,5:0,12:0,15:1},[0,1,0,0],None,'t=1/2; a arbitrary'),
        ('special_e0_hhalf',{16:0,18:0,17:Q(1,2),13:0,14:0,7:0,8:0,5:0,12:0,15:Q(1,2)},[0,1,0,0],None,'t=1/2; a arbitrary'),
        ('special_ehalf_h1',{16:0,18:0,17:Q(1,2),13:0,7:0,8:0,12:Q(1,2),15:1,
                             4:scalar(times(g,g),-2),5:scalar(g,-3)},[0,-1,0,0],scalar(g,2),'t=1/2; g arbitrary')]
    family_reports = []
    for label,assignments,weights,xchange,hypotheses in families:
        source = [put(v,assignments) for v in p]
        if xchange is not None:
            source = changed_x(source,xchange)
        top = [{ex:cv for ex,cv in v.items() if sum(ex[-4:])==2} for v in source]
        initial = []
        positive = []
        for (i,j),v in zip(PAIRS,top):
            ini = {}
            for ex,cv in v.items():
                weight = sum(power*w for power,w in zip(ex[-4:],weights))-weights[i]-weights[j]
                assert weight>=0
                if weight==0:
                    expected = [0]*4
                    expected[i] += 1
                    expected[j] += 1
                    assert tuple(expected)==ex[-4:]
                    ini[ex] = cv
                else:
                    positive.append(weight)
            initial.append(ini)
        ipf = pfaffian(R,initial)
        assert len(ipf)==1
        exponent,value = next(iter(ipf.items()))
        assert exponent==(0,)*19+(1,1,1,1) and value!=0
        # M=4 makes every constant exponent strictly positive for these weights.
        assert all(4-weights[i]-weights[j]>0 for i,j in PAIRS)
        family_reports.append(dict(label=label,hypotheses=hypotheses,
                                   substitutions={str(i):R.encode(R.constant(v) if isinstance(v,(int,Q)) else v) for i,v in assignments.items()},
                                   coordinate_change='X=x+2*g*y' if xchange is not None else 'identity',
                                   four_coordinate_weights=weights,constant_scaling_power=4,
                                   all_quadratic_weights_nonnegative=True,positive_off_log_weights=sorted(set(positive)),
                                   initial_six_brackets=[R.encode(v) for v in initial],initial_pfaffian_coefficient=str(value)))

    # The other special value is carried to t=+1/2 by (x,y,z,w)->(y,x,w,z).
    # Nonzero square r is impossible.  Check its complete nine equations and
    # the exact finite contradiction on g=0, plus the unipotent g->0 formula.
    nilp = [put(v,{16:0,17:0,18:1}) for v in p]
    F = [plus(times(plus(one,neg(h)),plus(one,scalar(h,-2))),g,scalar(b,3),scalar(times(b,f),2)),
         plus(h,f,neg(times(f,h)),neg(e),neg(times(b,d))),
         plus(scalar(g,2),neg(times(f,g)),times(e,plus(one,neg(h))),neg(times(b,c)),neg(times(a,d))),
         plus(neg(times(e,g)),neg(times(a,c))),
         plus(times(g,plus(R.constant(3),scalar(h,-4))),times(b,plus(one,scalar(e,-2))),scalar(a,-4),scalar(times(a,f),-2)),
         plus(scalar(times(g,g),2),times(a,plus(R.constant(-1),scalar(e,2)))),
         plus(times(f,plus(one,scalar(e,-4))),scalar(times(d,g),2),times(c,plus(R.constant(-1),scalar(h,2)))),
         plus(f,scalar(times(f,f),2),times(d,plus(one,scalar(h,-2))),c),
         plus(times(e,plus(scalar(e,2),R.constant(-1))),scalar(times(c,g),-2))]
    assert coefficient_ideal(nilp)=={normalize(v) for v in F}
    # e=0,g=0 => a=b=0; h=1 gives1=0; h=1/2 gives two incompatible linear equations.
    eq(R,put(F[1],{12:0,14:0,4:0,5:0,15:1}),one)
    eq(R,put(plus(scalar(F[1],2),neg(F[6])),{12:0,14:0,4:0,5:0,15:Q(1,2)}),one)
    # e=1/2,g=0: a!=0 forces c=f=0, where F5=-4a.
    eq(R,put(F[4],{12:Q(1,2),14:0,7:0,13:0}),scalar(a,-4))
    # On a=0, the following exact identity forces h=1, without dividing c.
    eq(R,put(plus(times(c,F[0]),times(plus(R.constant(3),scalar(f,2)),F[2]),times(plus(one,neg(h)),F[6])),
             {12:Q(1,2),14:0,4:0}),scalar(plus(one,neg(h)),Q(3,2)))
    eq(R,put(plus(F[1],times(scalar(d,Q(1,3)),plus(F[0],scalar(F[2],2),scalar(times(b,F[6]),2)))),
             {12:Q(1,2),14:0,4:0,15:1}),R.constant(Q(1,2)))
    # Preserve C=I and r=w^2 under y'=y-k*x,z'=z+k*w.
    k = R.variable('n0')
    transformed = [nilp[0],plus(nilp[1],times(k,nilp[2])),nilp[2],
                   plus(nilp[3],times(k,nilp[4]),neg(times(k,nilp[1])),neg(times(k,k,nilp[2]))),
                   plus(nilp[4],neg(times(k,nilp[2]))),nilp[5]]
    transformed = [substitute(R,substitute(R,v,'y',plus(R.variable('y'),times(k,R.variable('x')))),
                              'z',plus(R.variable('z'),neg(times(k,R.variable('w'))))) for v in transformed]
    yz = [0]*4
    yz[1]=yz[2]=1
    yz_coefficient = {ex[:-4]+(0,)*4:cv for ex,cv in transformed[4].items() if ex[-4:]==tuple(yz)}
    eq(R,yz_coefficient,plus(g,neg(times(k,a))))
    eq(R,transformed[5],times(R.variable('w'),R.variable('w')))

    # If both same-sign brackets vanish, all Jacobi are the two-by-two minors
    # of a symmetric rank-one matrix.  Its nonzero trace form gives a
    # semisimple J and a Cartan projection on Mat2, hence a product Pfaffian.
    zero_p = [put(v,{16:0,17:0,18:0,7:neg(f),5:g,15:plus(one,neg(e))}) for v in p]
    M = [[plus(R.constant(Q(1,2)),neg(e)),neg(g),neg(f)],
         [neg(g),a,e],[neg(f),e,d]]
    minors = [plus(times(M[i][j],M[k][l]),neg(times(M[i][l],M[k][j])))
              for i,k in it.combinations(range(3),2) for j,l in it.combinations(range(3),2)]
    assert coefficient_ideal(zero_p)=={normalize(v) for v in minors if v}
    eq(R,plus(M[0][0],M[1][2]),R.constant(Q(1,2)))
    return dict(all_type20_points_covered=True,
                rank_zero_and_rank_one='Excluded by earlier certificates.',
                rank_two_nonzero_square_binary_bracket_impossible=True,
                nilpotent_nine_equations=[R.encode(v) for v in F],
                nilpotent_g_zero_contradiction_identities_checked=True,
                nilpotent_unipotent_g_transform_checked=True,
                nilpotent_g_nonzero_chart='F6 implies a!=0; k=g/a sends g to0 and preserves C=I,r=w^2.',
                zero_zero_Jacobi_exactly_rank_one_symmetric_matrix=True,
                zero_zero_matrix=[[R.encode(v) for v in row] for row in M],
                zero_zero_semisimple_projection_argument='M=vv^T, v1^2+v2*v3=1/2; J^2=I/2. GL2 conjugates Cartan projection to diagonal, giving product Pfaffian.',
                semisimple_nine_equations=[R.encode(v) for v in E],
                two_generic_resonance_syzygies_checked=True,
                all_seven_weighted_families=family_reports,
                negative_special_value='Permutation (x,y,z,w)->(y,x,w,z) preserves C=I,S and changes r=tzw to -tzw.',
                r_zero_q_nonzero='Exchange the two eigenspaces and multiply P by -1; renormalize invertible C. Covered by nonzero r cases.',
                SNC_curve='P_s=D_s^{-1}(Q(D_s x)+s^4 C)D_s^{-T}; all coefficients polynomial in s, P_1=P, P_0 nondegenerate logcanonical with five independent hyperplanes.',
                containment='Every irreducible actual curve meets SNC; MPS Theorem1.4 places each containing full-space component in the established MPS union.')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-report',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--rank-two-followup',action='store_true')
    ap.add_argument('--complete-type20',action='store_true')
    args = ap.parse_args()
    raw = args.source_report.read_bytes()
    rows = json.loads(raw)['full_semisimple_sources']
    result = dict(task='E_reduced_nonnormal',
                  source_report=dict(path=str(args.source_report),sha256=hashlib.sha256(raw).hexdigest()),
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  source_bytes=Path(__file__).read_text(),
                  type23=type23(rows[22]),type20=type20(rows[19]),
                  type23_reducible=True,type20_support_open=True,
                  scope='Whole type23 excluded. Type20 normalized full rank1/rank2 constant sources only; residual supports retained.',
                  all_passed=True)
    if args.rank_two_followup:
        result['type20_followup'] = followup_type20(result['type20'])
    if args.complete_type20:
        result['type20_complete_support'] = complete_type20(result['type20'])
        result.update(type20_support_open=False,type20_all_points_covered=True,
                      scope='Whole type23 excluded. Whole type20 is reducible or covered by actual SNC curves and MPS Theorem1.4.')
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(type23_reducible=True,type20_rank_branch_free_indices=
                         [r['free_indices'] for r in result['type20']['nonzero_constant_rank_branches']])) )


if __name__=='__main__':
    main()

#!/usr/bin/env python3
"""Exact R0 restrictions of the unsaturated S61 infinity equations.

Certifies a binary-quadratic locus and a dominant nonbinary incidence family.
It does not identify the saturated projective closure or exhaust infinity.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_known_fiber_audit import equation_jacobian
from s61_model import S61, X, J as J5
from weight1_structure import COORDS, J, infinitesimal_action, minor_witness, sp_basis

ROOT = Path(__file__).resolve().parents[1]


def sm(matrix):
    return [[str(value) for value in row] for row in matrix.tolist()]


def build():
    model = S61()
    weights = model.reduced['parameter_weights']
    one = [i for i, a in enumerate(weights) if a == 1]
    two = [i for i, a in enumerate(weights) if a == 2]
    allu = s.symbols('u1:62')
    u = [allu[i] for i in one]
    g = s.symbols('g0:9')
    x, y, z, w = X[1:]
    monomials = [s.prod(v) for v in combinations_with_replacement([x, z, w], 3)][1:]
    G = sum(a*b for a, b in zip(g, monomials, strict=True))
    bivector = X[0]**2 * J5
    bivector[2, 3] += s.diff(G, w)
    bivector[3, 2] -= s.diff(G, w)
    bivector[2, 4] -= s.diff(G, z)
    bivector[4, 2] += s.diff(G, z)
    normalized = model.normalize(bivector)['parameters']
    assert all(normalized[i] == 0 for i in one)
    values = [allu[i] if i in one else normalized[i] for i in range(61)]
    extra_path = ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    extra = json.loads(extra_path.read_text())['universal_new_equation']
    encoded = model.reduced['equations'] + [extra]
    top = [[term for term in f if sum(i in two for i in term[0]) == 2] for f in encoded]
    equations = [s.expand(sum(s.Rational(n, d)*s.prod(values[i] for i in mon)
                              for mon, n, d in f)) for f in top]
    assert equations[:70] == [0]*70
    F, G6 = equations[70:134], equations[134]

    # A fixed Lagrangian plane <d_y,d_z>, with a completely arbitrary binary F.
    A, B, C = s.symbols('A B C')
    binary_sub = dict(zip(g, [0, A, 0, 0, B/2, 0, 0, 0, C/3], strict=True))
    binary = [s.expand(f.subs(binary_sub)) for f in F]
    r, t = allu[11], allu[12]
    h1, h2 = 2*A*r-B*t, B*r-2*C*t
    four = [-A*h1, -C*h2, -B*h1-A*h2, -C*h1-B*h2]
    assert set(binary)-{s.S.Zero} == set(map(s.expand, four))
    assert s.expand(G6.subs(binary_sub)) == 0
    assert all(s.expand(f.subs({r: 0, t: 0})) == 0 for f in binary)
    a, b = s.symbols('a b')
    assert s.expand(h1.subs({A:a*a, B:2*a*b, C:b*b})-2*a*(a*r-b*t)) == 0
    assert s.expand(h2.subs({A:a*a, B:2*a*b, C:b*b})-2*b*(a*r-b*t)) == 0

    # Recognition away from the vertex condition u4=u12=0. Constant row
    # combinations are stored, not an unverified Groebner calculation.
    open_F = [s.expand(f.subs({allu[3]:1, allu[11]:0})) for f in F]
    recognition = [g[6]**2-3*g[5]*g[7], g[6]*g[7]-9*g[5]*g[8],
        g[7]**2-3*g[6]*g[8], g[0]*g[6]-3*g[1]*g[5],
        2*g[2]*g[6]-3*g[3]*g[5], g[3]*g[6]-6*g[4]*g[5],
        2*g[1]*g[2]-g[0]*g[3],
        -2*g[0]*g[7]+2*g[1]*g[6]+4*g[2]*g[4]-g[3]**2,
        6*g[2]*g[8]-3*g[3]*g[7]+4*g[4]*g[6]]
    polys = [s.Poly(f, *g, *u) for f in open_F+recognition]
    mon = sorted(set().union(*(set(p.monoms()) for p in polys)))
    matrix = s.Matrix([[p.coeff_monomial(m) for p in polys[:64]] for m in mon])
    target = s.Matrix([[p.coeff_monomial(m) for p in polys[64:]] for m in mon])
    combinations, free = matrix.gauss_jordan_solve(target)
    combinations = combinations.subs(dict.fromkeys(free, 0))
    assert matrix*combinations == target
    case1 = {g[i]:0 for i in [2,3,5,6,7]}
    assert s.expand(open_F[5].subs(case1)-18*g[8]**2*allu[12]) == 0
    assert s.expand(open_F[3].subs(case1)-3*g[8]*(-g[0]+6*g[4]*allu[12])) == 0
    case2 = {g[i]:0 for i in [2,3,5,6,7,8]}
    assert s.expand(open_F[1].subs(case2)-g[4]*(-3*g[0]+4*g[4]*allu[12])) == 0
    assert s.expand(open_F[14].subs(case2)-g[0]*(2*g[0]+g[4]*allu[12])) == 0

    # A nonbinary solution over the vertex stratum; six beta variables stay free.
    k = s.symbols('k0:5')
    beta_sub = {allu[3]:0, allu[11]:0, allu[4]:0, allu[12]:0,
                allu[13]:k[0], allu[48]:2*k[0]+1, allu[25]:k[1],
                allu[34]:k[2], allu[49]:k[3], allu[50]:k[4]}
    AA = -(5*k[2]+k[3])/(3*k[0]-1)
    BB = (-5*k[1]+k[4])/(3*k[0]+4)
    nonbinary_sub = dict(zip(g, [AA, BB, 0, 1, 0, 0, 0, 0, 0], strict=True))
    assert all(s.cancel(f.subs(beta_sub).subs(nonbinary_sub)) == 0 for f in equations)
    generic_shape = dict(zip(g, [A, B, 0, 1, 0, 0, 0, 0, 0], strict=True))
    assert s.expand(G6.subs(generic_shape)) == 0
    gradient = s.Matrix([[g[0],2*g[2],g[3],3*g[5],2*g[6],g[7]],
                         [g[1],g[3],2*g[4],g[6],2*g[7],3*g[8]]])
    # The xz,xw columns are already independent, irrespective of A and B.
    assert gradient[:,[1,2]].subs(generic_shape).det() == -1
    sample_k = dict(zip(k, range(1,6), strict=True))
    sample_beta = dict.fromkeys(u, s.S.Zero)
    sample_beta.update({key:s.sympify(value).subs(sample_k) for key,value in beta_sub.items()})
    sample_g = {key:s.cancel(value.subs(sample_k)) if hasattr(value,'subs') else value
                for key,value in nonbinary_sub.items()}
    fixed_jacobian = s.Matrix(F+[G6]).jacobian(g).subs(sample_beta).subs(sample_g)
    fixed_witness = minor_witness(fixed_jacobian)
    assert fixed_witness['rank'] == 8

    # Direct coefficient proof that beta -> tau_beta in H0(N(2)) is injective.
    zz, uu, beta = primitive_beta(model)
    assert uu == u
    v = s.Matrix(s.symbols('v0:4'))
    at_v = beta.subs(dict(zip(zz,v,strict=True)))
    lam = at_v.T*v
    omega = J.T*v
    cubic = [s.expand(lam[i]*omega[j]-lam[j]*omega[i])
             for i in range(4) for j in range(i+1,4)]
    cubic_polys = [s.Poly(f,*v) for f in cubic]
    vmons = sorted(set().union(*(set(p.monoms()) for p in cubic_polys)))
    coefficient = s.Matrix([p.coeff_monomial(mm) for p in cubic_polys for mm in vmons]).jacobian(u)
    tau_witness = minor_witness(coefficient)
    assert tau_witness['rank'] == 16
    vv = s.Matrix([0,1,0,0])
    beta_v = beta.subs(dict(zip(zz,vv,strict=True)))
    assert beta_v[1,2] == -allu[11] and beta_v[1,3] == allu[3]

    # Dominance of Sp4 times the eleven-dimensional normalized beta family.
    vec = s.Matrix([s.diff(beta[i,j],zz[h]) for (i,j),h in COORDS])
    basis = vec.jacobian(u)
    basis_rows = list(basis.T.rref()[1])
    left = basis[basis_rows,:].inv()
    beta_point = vec.subs(sample_beta)
    action = s.Matrix.hstack(*[left*infinitesimal_action(beta_point,aa)[basis_rows,:]
                               for aa in sp_basis()])
    fixed = {allu[i] for i in [3,11,4,12,48]}
    family, family_variables = [], []
    for value in u:
        if value in fixed:
            continue
        column = s.zeros(16,1)
        column[u.index(value)] = 1
        if value == allu[13]:
            column[u.index(allu[48])] = 2
        family.append(column)
        family_variables.append(str(value))
    dominance = s.Matrix.hstack(*family,action)
    dominance_witness = minor_witness(dominance)
    assert len(family) == 11 and dominance_witness['rank'] == 16
    assert s.Rational(dominance_witness['determinant']) == 48
    # Pulling back beta by exp(tA) transports its distinguished line by exp(-tA).
    line_action = s.Matrix.hstack(*[(-aa*vv)[[0,2,3],:] for aa in sp_basis()])
    incidence = dominance.col_join(s.zeros(3,11).row_join(line_action))
    incidence_witness = minor_witness(incidence)
    assert incidence_witness['rank'] == 17
    assert s.Rational(incidence_witness['determinant']) == 20
    point = [s.cancel(value.subs(sample_beta).subs(sample_g)) if hasattr(value,'subs') else value
             for value in values]
    full_jacobian = equation_jacobian(top,point)[:,two]
    full_witness = minor_witness(full_jacobian)
    assert full_witness['rank'] == 42

    inputs = [Path(__file__), ROOT/'code/s61_model.py',
              ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
              ROOT/'code/weight1_structure.py', ROOT/'code/s61_known_fiber_audit.py',
              model.source/'symplectic_slice_reduced_equations.json', extra_path]
    return {'id':'S61_R0_INFINITY_20261101','field':'QQ','all_passed':True,
        'g_variables':list(map(str,g)),'beta_variables':list(map(str,u)),
        'G':str(G),'q_coordinates':list(map(str,[normalized[i] for i in two])),
        'weight_two_u_indices_zero_based':two,
        'leading_equations_70_through_133':list(map(str,F)), 'leading_G6':str(G6),
        'first_seventy_restrictions_zero':True,
        'binary_plane':{'plane':'span(d_y,d_z)','quadratic':'A*x^2+B*x*w+C*w^2',
            'four_equations':list(map(str,map(s.expand,four))),
            'h1':str(h1),'h2':str(h2),'G6_identically_zero':True},
        'open_recognition':{'normalization':'u4=1,u12=0',
            'nine_consequences':list(map(str,recognition)),
            'constant_64_row_combinations':sm(combinations.T),
            'exceptional_case_rows_zero_based':[5,3,1,14],
            'gradient_matrix':sm(gradient)},
        'nonbinary':{'beta_substitution':{str(a):str(b) for a,b in beta_sub.items()},
            'G':'A*x^2*z+B*x^2*w+x*z*w','A':str(AA),'B':str(BB),
            'all_135_identities':True,'gradient_minor':'-1',
            'sample_beta':{str(a):str(b) for a,b in sample_beta.items()},
            'sample_g':{str(a):str(b) for a,b in sample_g.items()},
            'sample_S61_coordinates':list(map(str,point)),
            'fixed_v_jacobian_minor':fixed_witness,
            'full_infinity_jacobian_minor':full_witness,
            'full_projective_tangent_dimension':2,
            'family_variables':family_variables,'dominance_matrix':sm(dominance),
            'dominance_minor':dominance_witness,'beta_line_incidence_matrix':sm(incidence),
            'beta_line_incidence_minor':incidence_witness},
        'null_correlation':{'tau_at_dy':['-u12','u4'],
            'tau_zero_projective_cubics':list(map(str,cubic)),
            'injectivity_coefficient_monomials':list(map(list,vmons)),
            'injectivity_coefficient_matrix':sm(coefficient),'injectivity_minor':tau_witness},
        'scope':'Unsaturated leading equations only. No affine extension, saturation, '
                'full infinity exhaustion, or full S61 generic-fiber upper bound.',
        'proof_note':'research/notes/S61_R0_INFINITY_20261101.md',
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                        for path in inputs}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=build()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps({'all_passed':True,'dominance_rank':16,'fixed_v_rank':8,
                      'full_infinity_rank':42,'output':str(args.output)}))


if __name__=='__main__':
    main()

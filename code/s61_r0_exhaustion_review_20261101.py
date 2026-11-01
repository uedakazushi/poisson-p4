#!/usr/bin/env python3
"""Independent source/normal-form audit for the reduced generic R0 infinity.

Rebuilds the 135 leading restrictions from S61, uses a different Groebner
order, and retains the two nuisance coordinates suppressed in the main code.
The scope is R0 in unsaturated infinity, not the entire S61 fiber.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
from pathlib import Path
import time
import sympy as s
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61, X, J as J5
from weight1_structure import COORDS, J

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    return json.loads((ROOT/'reports'/name).read_text())


def verify():
    started = time.monotonic()
    old = load('S61_R0_INFINITY_20261101.json')
    uniform = load('S61_R0_UNIFORM_VERTEX_20261101.json')
    cases = load('S61_R0_GENERIC_EXHAUSTION_20261101.json')
    extra = load('S61_UNIVERSAL_WEIGHT6_20261031.json')['universal_new_equation']
    model = S61()
    weights = model.reduced['parameter_weights']
    one = [i for i, w in enumerate(weights) if w == 1]
    two = [i for i, w in enumerate(weights) if w == 2]
    au = s.symbols('u1:62')
    g = s.symbols('g0:9')
    x, y, z, w = X[1:]
    cubics = [x*x*z, x*x*w, x*z*z, x*z*w, x*w*w, z**3, z*z*w, z*w*w, w**3]
    G = sum(a*b for a, b in zip(g, cubics, strict=True))
    P = X[0]**2*J5
    P[2, 3] += s.diff(G, w); P[3, 2] -= s.diff(G, w)
    P[2, 4] -= s.diff(G, z); P[4, 2] += s.diff(G, z)
    q = model.normalize(P)['parameters']
    assert all(q[i] == 0 for i in one)
    values = [au[i] if i in one else q[i] for i in range(61)]
    raw = []
    for equation in model.reduced['equations']+[extra]:
        raw.append(s.expand(sum(s.Rational(n, d)*s.prod(values[i] for i in mon)
            for mon, n, d in equation if sum(i in two for i in mon) == 2)))
    assert raw[:70] == [0]*70
    F = raw[70:]
    assert all(s.expand(a-s.sympify(b)) == 0 for a, b in zip(F,
        old['leading_equations_70_through_133']+[old['leading_G6']], strict=True))
    assert [str(q[i]) for i in two] == old['q_coordinates']

    # Recheck the nine universal consequences and the exceptional rows used
    # in the off-vertex recognition proof; no Groebner inference is used here.
    opened = s.Matrix([f.subs({au[3]: 1, au[11]: 0}) for f in F[:64]])
    rec = old['open_recognition']
    combination = s.Matrix([[s.Rational(v) for v in row]
                            for row in rec['constant_64_row_combinations']])
    desired = s.Matrix(list(map(s.sympify, rec['nine_consequences'])))
    assert (combination*opened-desired).applyfunc(s.expand) == s.zeros(9, 1)
    sub1 = dict.fromkeys([g[i] for i in [2, 3, 5, 6, 7]], 0)
    assert s.expand(opened[5].subs(sub1)-18*g[8]**2*au[12]) == 0
    assert s.expand(opened[3].subs(sub1)-3*g[8]*(-g[0]+6*g[4]*au[12])) == 0
    sub2 = dict.fromkeys([g[i] for i in [2, 3, 5, 6, 7, 8]], 0)
    assert s.expand(opened[1].subs(sub2)-g[4]*(-3*g[0]+4*g[4]*au[12])) == 0
    assert s.expand(opened[14].subs(sub2)-g[0]*(2*g[0]+g[4]*au[12])) == 0

    h = s.symbols('h0:9'); k = s.Symbol('k')
    fixed = {au[3]: 0, au[11]: 0, au[4]: 0, au[12]: 0,
             au[13]: k, au[48]: 2*k+1, au[25]: 0, au[34]: 0, au[49]: 0, au[50]: 0}
    UF = [s.expand(f.subs(fixed).subs(dict(zip(g, h, strict=True)))) for f in F]
    assert all(s.expand(a-s.sympify(b)) == 0 for a, b in zip(UF, uniform['original_equations'], strict=True))
    # Main certificate uses grevlex; this review recomputes grlex from the
    # newly reconstructed original generators.
    UG = s.groebner(UF, *h, k, order='grlex', domain=s.QQ)
    assert len(UG.polys) == 153
    generators = list(map(s.sympify, uniform['radical_generators']))
    for v, power in zip(generators, uniform['generator_powers_in_original_ideal'], strict=True):
        assert UG.reduce(v**power)[1] == 0
    coefficients = []
    edges = []
    for v in generators:
        poly = s.Poly(v, *h)
        assert len(poly.terms()) == 1
        ex, coefficient = poly.terms()[0]
        assert sum(ex) == 2 and max(ex) == 1
        coefficients.append(coefficient)
        edges.append({i for i, e in enumerate(ex) if e})
    for f in UF:
        for ex, co in s.Poly(f, *h).terms():
            if co:
                edge = {i for i, e in enumerate(ex) if e}
                assert max(ex) == 1 and edge in edges
                assert s.rem(co, coefficients[edges.index(edge)], k) == 0
    generic_faces = list(map(set, uniform['generic_coordinate_faces']))
    vertical = {s.Rational(a): list(map(set, fs))
                for a, fs in uniform['extra_vertical_prime_faces'].items()}
    support_count = 0
    for mask in range(512):
        support = {i for i in range(9) if mask & (1 << i)}
        coeff = s.S.Zero
        for edge, c in zip(edges, coefficients, strict=True):
            if edge <= support: coeff = s.gcd(coeff, c)
        if coeff: coeff = s.Poly(coeff, k).monic().as_expr()
        expect = (s.S.Zero if any(support <= f for f in generic_faces) else
            s.prod(k-a for a, fs in vertical.items() if any(support <= f for f in fs)))
        assert s.expand(coeff-expect) == 0
        support_count += 1

    case_sizes = {}
    for case in cases['fixed_normal_form_radicals']:
        sub = {s.Symbol(a): s.sympify(b) for a, b in case['beta_substitution'].items()}
        eq = [s.expand(f.subs(sub)) for f in F]
        assert all(s.expand(a-s.sympify(b)) == 0 for a, b in zip(eq, case['equations'], strict=True))
        gb = s.groebner(eq, *g, order='grlex', domain=s.QQ)
        gs = list(map(s.sympify, case['radical_generators']))
        for v, power in zip(gs, case['generator_powers_in_original_ideal'], strict=True):
            assert gb.reduce(v**power)[1] == 0
        mon = [s.Poly(v, *g).monoms()[0] for v in gs]
        assert all(max(v) == 1 for v in mon)
        for f in eq:
            assert all(any(all(a <= b for a, b in zip(v, e, strict=True)) for v in mon)
                       for e, co in s.Poly(f, *g).terms() if co)
        case_sizes[case['label']] = len(gb.polys)

    # Check all normalizing transformations with BOTH nuisance coefficients
    # u26,u35 retained, rather than only the representatives with them zero.
    zz, u, beta = primitive_beta(model)
    coefficient_vector = s.Matrix([beta[i, j].diff(zz[a]) for (i, j), a in COORDS])
    basis = coefficient_vector.jacobian(u)
    rows = list(basis.T.rref()[1]); inverse_basis = basis[rows, :].inv()
    K, L, a1, a2, A, B, t = s.symbols('K L a1 a2 A B t', nonzero=True)
    T = s.Matrix([[1, 0, 0, 0], [0, 1, A, -B], [B, 0, 1, 0], [A, 0, 0, 1]])
    assert T.T*J*T == J
    for rank_one in [False, True]:
        sub = {au[3]: 0, au[11]: 0, au[4]: -int(rank_one), au[12]: 0,
            au[13]: k, au[48]: 2*k+int(not rank_one), au[25]: a1, au[34]: a2,
            au[49]: K-5*a2, au[50]: 5*a1-L}
        assert all(not (f.subs(sub).expand().free_symbols & {a1, a2}) for f in F)
        def pull(M):
            bb = beta.subs(sub)
            pb = M.T*bb.subs(dict(zip(zz, M*s.Matrix(zz), strict=True)), simultaneous=True)*M
            vec = s.Matrix([pb[i, j].diff(zz[a]) for (i, j), a in COORDS])
            return (inverse_basis*vec[rows, :]).applyfunc(s.cancel)
        nu = pull(T.inv())
        assert all(nu[u.index(au[i])] == sub[au[i]] for i in [3, 11, 4, 12, 13, 48])
        expectedK = 3*k*A+K if rank_one else (3*k-1)*A+K
        expectedL = -5*A+3*k*B+L if rank_one else (3*k+4)*B+L
        assert s.expand(5*nu[u.index(au[34])]+nu[u.index(au[49])]-expectedK) == 0
        assert s.expand(5*nu[u.index(au[25])]-nu[u.index(au[50])]-expectedL) == 0
        torus = s.diag(t*t, t**-2, t**-1, t) if rank_one else s.diag(1, 1, t**-1, t)
        assert torus.T*J*torus == J
        nu = pull(torus)
        assert nu[u.index(au[13])] == (k/t**2 if rank_one else k)
        assert nu[u.index(au[4])] == -int(rank_one)
        assert s.cancel(5*nu[u.index(au[34])]+nu[u.index(au[49])]-K/t) == 0
        assert s.cancel(5*nu[u.index(au[25])]-nu[u.index(au[50])]-L*t) == 0

    v = s.Matrix([0, 1, 0, 0]); a, b, c = s.symbols('a b c')
    beta_v = beta.subs(dict(zip(zz, v, strict=True)))
    vv = s.Matrix([0, a, b, c])
    qv = s.expand((vv.T*beta.subs(dict(zip(zz, vv, strict=True)))*v)[0])
    qp = s.Poly(qv, a, b, c)
    evaluation = s.Matrix([beta_v[1, 2], beta_v[1, 3],
        qp.coeff_monomial(b*b), qp.coeff_monomial(b*c), qp.coeff_monomial(c*c),
        beta_v[0, 1]]).jacobian(u)
    pivots = list(evaluation.rref()[1]); assert len(pivots) == 6
    determinant = evaluation[:, pivots].det(); assert determinant != 0

    paths = [Path(__file__), ROOT/'code/s61_model.py',
        ROOT/'code/primitive_jet_elliptic_curve_20261101.py', ROOT/'code/weight1_structure.py',
        model.source/'symplectic_slice_reduced_equations.json']
    paths += [ROOT/'reports'/name for name in ['S61_R0_INFINITY_20261101.json',
        'S61_R0_UNIFORM_VERTEX_20261101.json', 'S61_R0_GENERIC_EXHAUSTION_20261101.json',
        'S61_UNIVERSAL_WEIGHT6_20261031.json']]
    return {'id': 'S61_R0_EXHAUSTION_REVIEW_20261101', 'all_passed': True,
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'independently_rebuilt_leading_restrictions': 135,
        'off_vertex_nine_consequences_and_four_exceptional_rows_checked': True,
        'alternate_groebner_order': 'grlex', 'uniform_groebner_size': len(UG.polys),
        'uniform_radical_power_memberships_checked': len(generators),
        'relative_prime_intersection_supports_checked': support_count,
        'fixed_case_groebner_sizes': case_sizes,
        'normalizing_maps_checked_with_u26_u35_retained': True,
        'evaluation_rank': 6, 'evaluation_minor_columns': pivots,
        'evaluation_minor': str(determinant),
        'excluded_incidence_dimension_bounds': [14, 15],
        'scope': 'Reduced R0 part of unsaturated infinity for sufficiently general beta only; no full Q45 or affine-fiber exhaustion.',
        'proof_note': 'research/notes/S61_R0_EXHAUSTION_REVIEW_20261101.md',
        'elapsed_seconds': round(time.monotonic()-started, 3),
        'input_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); result = verify()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({'all_passed': True, 'seconds': result['elapsed_seconds'],
                      'alternate_groebner_sizes': result['fixed_case_groebner_sizes']}))


if __name__ == '__main__':
    main()

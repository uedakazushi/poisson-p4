#!/usr/bin/env python3
"""Audit the full fixed-trace Pfaffian reductions directly from archived forms."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def verify():
    p05 = ROOT/'legacy/v05/data/minimal_nilpotent_certificate.json'
    p06 = ROOT/'legacy/v06/data/nilpotent22_certificate.json'
    src = ROOT/'reports/TORIC_ZERO_JET_NONZERO_TRACE_20261101.json'
    old5, old6, saved = (json.loads(p.read_text()) for p in [p05, p06, src])
    x, y, z, w = s.symbols('x y z w')
    a, b, c, d, u, v, A, B, p3, p4, p5 = s.symbols('a b c d u v A B p3 p4 p5')
    coords = dict(zip(s.symbols('x1:5'), [x, y, z, w]))
    t = s.symbols('t1:16')
    sub = coords | {t[9]: c, t[10]: a, t[11]: b, t[12]: 0, t[13]: 0, t[14]: 0}
    q5 = {k: s.sympify(f).subs(sub, simultaneous=True) for k, f in old5['normal_form_brackets'].items()}
    eq5 = [s.expand(s.sympify(f).subs(sub, simultaneous=True)) for f in old5['original_Jacobi_generators']]
    # Verify the archived nilpotence identities before passing to reduced support.
    full = [s.sympify(f) for f in old5['original_Jacobi_generators']]
    assert s.expand(full[13]-t[12]**2) == 0
    assert s.expand(full[9]-t[14]**2) == 0
    assert s.expand(2*t[13]*full[7]-2*t[14]*full[6]-t[13]**3) == 0
    e = [s.sympify(f) for f in saved['minimal']['jacobi_generators']]
    chosen = [0, 1, 2, 10, 11, 14]
    factors = [1, 3, s.Rational(3, 2), 1, 3, -3]
    assert all(s.expand(ei-factor*eq5[j]) == 0 for ei, j, factor in zip(e, chosen, factors))
    assert all(not f for j, f in enumerate(eq5) if j not in chosen)
    def pf(q):
        return s.expand(q['12']*q['34']-q['13']*q['24']+q['14']*q['23'])
    assert s.expand(pf(q5)-s.sympify(saved['minimal']['reduced_pfaffian'])
                    -s.sympify(saved['minimal']['exact_ideal_correction'])) == 0
    assert s.expand(q5['12']+q5['34']-x*x/2) == 0
    assert s.factor(e[-1].subs({a: 0, b: 0})) == 3*c*(4*c-1)
    params = [a, b, p3, p4, p5, c, v, d, A, B, u]
    sub = coords | dict(zip(t[:11], params))
    q6 = {k: s.sympify(f).subs(sub, simultaneous=True) for k, f in old6['normal_form_brackets'].items()}
    eq6 = [s.expand(s.sympify(f).subs(sub, simultaneous=True)) for f in old6['Jacobi_generators']]
    P = s.Poly(pf(q6), x, y, z, w)
    assert s.expand(P.coeff_monomial(x**3*y)-A-eq6[0]) == 0
    assert s.expand(P.coeff_monomial(z**3*w)+B-eq6[3]) == 0
    assert s.expand(q6['12']+q6['34']-x*z) == 0
    core = [s.sympify(f) for f in saved['nil22']['core_generators']]
    assert all(s.expand(f-factor*eq6[i].subs({A: 0, B: 0})) == 0
               for f, i, factor in zip(core, [4, 5, 6, 7, 8], [1, 1, 2, 1, 3]))
    assert s.expand(P.as_expr().subs({A: 0, B: 0})
                    -s.sympify(saved['nil22']['reduced_pfaffian'])
                    -s.sympify(saved['nil22']['exact_ideal_correction'])) == 0
    r, q, ell = s.symbols('r q ell')
    H = r*x*x+q*x*z+ell*z*z
    dets = []
    for sign in [-1, 1]:
        Hess = s.hessian(x*z+sign*2*H, (x, z))/2
        dets.append(s.expand(Hess.det()))
    assert s.expand(dets[0]-dets[1]) == 2*q
    assert dets[0].subs(q, 0) == 4*r*ell-s.Rational(1, 4)
    ab = {A: 0, B: 0, u: 0, v: 0, c: 0, d: s.Rational(1, 2), b: -1/(16*a)}
    assert all(s.cancel(f.subs(ab, simultaneous=True)) == 0 for f in eq6)
    for key, f in q6.items():
        assert s.cancel(f.subs(ab, simultaneous=True)-s.sympify(saved['nil22']['final_AB_bivector'][key])) == 0
    paths = [Path(__file__), p05, p06, src, ROOT/'reports/TORIC_ZERO_JET_MINIMAL_REVIEW_20261101.json']
    return dict(id='TORIC_ZERO_JET_NONZERO_TRACE_REVIEW_20261101', all_passed=True,
        minimal_nilpotence_identities_verified=3,
        minimal_remaining_generators_verified=6, nil22_original_generators_verified=9,
        exact_Pfaffian_identities_verified=True,
        accepted_conclusion='All nonzero-trace points of actual reduced F0 lie in closed AB.',
        scope='Historical exhaustive linear fixed-trace reductions retained as dependencies; no zero-modular conclusion.',
        proof_note='research/notes/TORIC_ZERO_JET_NONZERO_TRACE_REVIEW_20261101.md',
        input_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); result = verify()
    with args.output.open('x') as f:
        json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps(dict(all_passed=True)))


if __name__ == '__main__':
    main()

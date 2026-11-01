#!/usr/bin/env python3
"""All nonzero-trace zero-jet candidates reduce to AB by Pfaffian identities."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    x, y, z, w = s.symbols('x y z w')
    a, b, c = s.symbols('a b c')
    t = s.symbols('t1:10')
    q12 = (s.Rational(1, 2)-c)*x*x-a*x*z-b*x*w
    q13, q14 = b*x*x/3, -a*x*x/3
    q23 = t[0]*x*x+(-2*c-t[5])*x*z-(a+t[7])*z*z/2+b*x*y/3-(b+2*t[8])*z*w+t[1]*x*w+t[2]*w*w
    q24 = -a*x*y/3+t[3]*x*x+t[4]*x*z+t[5]*x*w+t[6]*z*z+t[7]*z*w+t[8]*w*w
    q34 = c*x*x+a*x*z+b*x*w
    e = [-12*c*b-2*a*t[1]-2*b*t[5]+b,
         5*a*b-2*a*t[8]+b*t[7], a*t[2]+3*b*b+b*t[8],
         -8*c*a+2*a*t[5]+a-2*b*t[4],
         5*a*a-a*t[7]+2*b*t[6], 4*a*t[0]+12*c*c-3*c+4*b*t[3]]
    pf = s.expand(q12*q34-q13*q24+q14*q23)
    target = c*x**4/4+(a*z+b*w)*x**3/3
    correction = -e[5]*x**4/12+e[3]*x**3*z/6+e[0]*x**3*w/6-e[4]*x*x*z*z/6-e[1]*x*x*z*w/3-e[2]*x*x*w*w/3
    assert s.expand(pf-target-correction) == 0
    assert s.expand(e[5].subs({a: 0, b: 0})-3*c*(4*c-1)) == 0

    historical = ROOT/'legacy/v06/data/nilpotent22_certificate.json'
    source = json.loads(historical.read_text())
    p3, p4, p5, v, d, A, B, u = s.symbols('p3 p4 p5 v d A B u')
    params = [a, b, p3, p4, p5, c, v, d, A, B, u]
    sub = {s.Symbol('t'+str(i+1)): value for i, value in enumerate(params)}
    sub.update(dict(zip(s.symbols('x1:5'), [x, y, z, w])))
    q = {key: s.sympify(value).subs(sub, simultaneous=True) for key, value in source['normal_form_brackets'].items()}
    pf22 = s.Poly(s.expand(q['12']*q['34']-q['13']*q['24']+q['14']*q['23']), x, y, z, w)
    E0 = 6*B*b-c*A+5*d*A-2*A
    E3 = 6*a*A+B*c+3*B*d-2*B
    assert s.expand(pf22.coeff_monomial(x**3*y)-A+E0/2) == 0
    assert s.expand(pf22.coeff_monomial(z**3*w)+B+E3/2) == 0
    core = [2*u*b+v*(c-2*d+1), -4*a*b-56*u*v+(c-d)*(c-d+1),
            2*a*v-u*c, a*(1-2*c-2*d)-4*u*u, b*(2*c-6*d+3)-4*v*v]
    closed_pf = s.expand(pf22.as_expr().subs({A: 0, B: 0}))
    closed_target = b*x**4/4-v*x**3*z+(d-c)*x*x*z*z/4+u*x*z**3-a*z**4/4
    closed_correction = core[4]*x**4/4-2*core[0]*x**3*z-core[1]*x*x*z*z/4-2*core[2]*x*z**3+core[3]*z**4/4
    assert s.expand(closed_pf-closed_target-closed_correction) == 0
    r, s0, s1 = s.symbols('r s0 s1')
    square_root = r*x*x+s0*x*z+s1*z*z
    determinants = []
    for sign in (1, -1):
        form = s.Poly(x*z+sign*2*square_root, x, z)
        determinants.append(s.expand(form.coeff_monomial(x*x)*form.coeff_monomial(z*z)-form.coeff_monomial(x*z)**2/4))
    assert s.expand(determinants[0]-determinants[1]) == -2*s0
    assert determinants[0].subs(s0, 0) == 4*r*s1-s.Rational(1, 4)
    AB = {A: 0, B: 0, u: 0, v: 0, c: 0, d: s.Rational(1, 2), b: -1/(16*a)}
    assert all(s.cancel(f.subs(AB, simultaneous=True)) == 0 for f in core)
    abq = {key: s.factor(value.subs(AB, simultaneous=True)) for key, value in q.items()}
    assert abq['12'] == x*z/2 and abq['34'] == x*z/2 and abq['13'] == 0
    assert s.expand(abq['14']-a*z*z-x*x/4) == 0
    assert s.expand(abq['23']+x*x/(16*a)+z*z/4) == 0
    assert s.expand(abq['24']-p3*x*x-p4*x*z-p5*z*z) == 0
    paths = [Path(__file__), historical, ROOT/'legacy/v05/proofs/minimal_nilpotent_theorem.txt',
             ROOT/'research/notes/TORIC_ZERO_JET_MINIMAL_20261101.md',
             ROOT/'research/notes/TORIC_ZERO_JET_FIBER_20261101.md']
    report = dict(id='TORIC_ZERO_JET_NONZERO_TRACE_20261101', all_passed=True, field='QQ',
                  minimal=dict(jacobi_generators=list(map(str, e)), pfaffian=str(pf),
                               reduced_pfaffian=str(target), exact_ideal_correction=str(correction)),
                  nil22=dict(outside_support_coefficients=[str(pf22.coeff_monomial(x**3*y)), str(pf22.coeff_monomial(z**3*w))],
                             support_corrections=[str(E0), str(E3)], core_generators=list(map(str, core)),
                             reduced_pfaffian=str(closed_target), exact_ideal_correction=str(closed_correction),
                             square_form_determinants=list(map(str, determinants)),
                             final_AB_bivector={key: str(value) for key, value in abq.items()}),
                  conclusion='Every point of the reduced actual T01 zero-first-jet fiber with nonzero trace belongs to closed AB.',
                  limitations=['Uses the historical exhaustive reduced fixed-trace normal forms, not only generic orbit representatives.',
                               'Zero-trace C0/L0 boundary points remain separate; no full zero-fiber exhaustion.'],
                  proof_note='research/notes/TORIC_ZERO_JET_NONZERO_TRACE_20261101.md',
                  source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps(dict(all_passed=True, nonzero_trace_intersection_contained_in_AB=True)))


if __name__ == '__main__':
    main()

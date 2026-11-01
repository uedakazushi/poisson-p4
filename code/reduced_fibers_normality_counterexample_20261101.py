#!/usr/bin/env python3
"""Reduced connected fibers of a proper smooth-source map need not imply normality."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    u, v, up, vp, x, y, z, r, q = s.symbols('u v up vp x y z r q')
    chart0 = {x: u, y: u*v*v, z: u*v}
    chart1 = {x: up*vp*vp, y: up, z: up*vp}
    assert s.expand((x*y-z*z).subs(chart0)) == 0
    assert s.expand((x*y-z*z).subs(chart1)) == 0
    transition = {up: u*v*v, vp: 1/v}
    assert all(s.cancel(chart1[a].subs(transition, simultaneous=True)-chart0[a]) == 0 for a in (x, y, z))
    generators = [x, y, x*z, y*z]
    pull0 = [s.expand(a.subs(chart0)) for a in generators]
    pull1 = [s.expand(a.subs(chart1)) for a in generators]
    assert list(s.groebner(pull0, u, v).polys)[0].as_expr() == u
    assert len(s.groebner(pull0, u, v).polys) == 1
    assert list(s.groebner(pull1, up, vp).polys)[0].as_expr() == up
    assert len(s.groebner(pull1, up, vp).polys) == 1
    # Eliminate the normalization coordinate z to display the nonnormal affine ring.
    graph = s.groebner([x*y-z*z, r-x*z, q-y*z], z, r, q, x, y, order='lex')
    relations = [g.as_expr() for g in graph.polys if not g.as_expr().has(z)]
    assert relations and all(s.rem(g.subs({r: x*z, q: y*z}), x*y-z*z, z) == 0 for g in relations)
    report = dict(id='REDUCED_FIBERS_NORMALITY_COUNTEREXAMPLE_20261101', all_passed=True, field='QQ',
                  normal_ring='A=C[x,y,z]/(xy-z^2)', nonnormal_ring='B=C[x,y,xz,yz] subset A',
                  B_relations_in_r_q_x_y=list(map(str, relations)),
                  normalization_fiber='A/m_B A = C[z]/(z^2)',
                  resolution='Tot(O_P1(-2)) -> Spec(A) -> Spec(B), projective and birational',
                  chart0={str(k): str(a) for k, a in chart0.items()}, chart1={str(k): str(a) for k, a in chart1.items()},
                  transition={str(k): str(a) for k, a in transition.items()},
                  pullback_origin_ideal_chart0=list(map(str, pull0)), pullback_origin_ideal_chart1=list(map(str, pull1)),
                  reduced_origin_fiber='P1; chart ideals (u) and (up)',
                  other_fibers='Reduced points; B=A after localizing x or y.',
                  scope='Counterexample to deducing AB normality from the reduced connected fibers of its normal proper birational incidence.',
                  proof_note='research/notes/REDUCED_FIBERS_NORMALITY_COUNTEREXAMPLE_20261101.md',
                  source_sha256={str(Path(__file__).relative_to(ROOT)): hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps(dict(all_passed=True, shortcut_refuted=True)))


if __name__ == '__main__':
    main()

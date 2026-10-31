#!/usr/bin/env python3
"""Bounded exact one-variable viability probe for the elliptic five-variable model.

Fix lambda=2, p=(t,2,3,4), target at t=1, and eliminate the original scalar.
Two necessary invariant equations have numerator gcd t-1 for BOTH families.
This proves the restricted open line fiber is one reduced point. It does not
compute the five-variable fiber length or either elliptic component degree.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import signal
import time

import sympy as s
from sympy.polys.domains import QQ

from elliptic_six_variable_fiber import normalized_jet, covariant_frame, affine_family, COORDS, J

ROOT = Path(__file__).resolve().parents[1]


def coefficients(polynomial, variable):
    return [str(v) for v in s.Poly(polynomial.as_expr(), variable).all_coeffs()]


def probe(kind):
    signal.alarm(60)
    started = time.monotonic()
    t = s.Symbol('t')
    field = QQ.frac_field(t)
    value = normalized_jet(kind, 2, [t, 2, 3, 4], 1).applyfunc(s.cancel)
    tensor = [[[field.zero for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for coefficient, ((b, c), a) in zip(value, COORDS, strict=True):
        tensor[a][b][c] = field.from_sympy(coefficient)
        tensor[a][c][b] = -tensor[a][b][c]
    jet_seconds = time.monotonic()-started
    nz = [(i, j, int(J[i, j])) for i in range(4) for j in range(4) if J[i, j]]
    moment = [[sum((tensor[i][j][k]*tensor[a][b][c]*u*v
                    for j, b, u in nz for k, c, v in nz), field.zero)
               for a in range(4)] for i in range(4)]
    operator = [[sum((int(J[i, k])*moment[k][j] for k in range(4)), field.zero)
                 for j in range(4)] for i in range(4)]
    raised = [[sum((operator[i][k]*int(J[j, k]) for k in range(4)), field.zero)
               for j in range(4)] for i in range(4)]
    covector = [sum((tensor[a][b][i]*raised[a][b]
                    for a in range(4) for b in range(4)), field.zero) for i in range(4)]
    vector = [sum((int(J[i, j])*covector[j] for j in range(4)), field.zero)
              for i in range(4)]
    moved = [sum((operator[i][j]*vector[j] for j in range(4)), field.zero)
             for i in range(4)]
    i4 = sum((operator[i][j]*operator[j][i]
              for i in range(4) for j in range(4)), field.zero)
    i14 = sum((moved[a]*vector[b]*moved[c]*tensor[a][b][c]
               for a in range(4) for b in range(4) for c in range(4)), field.zero)
    h = i14/i4**3
    g8 = sum((vector[i]*int(J[i, j])*moved[j]
              for i in range(4) for j in range(4)), field.zero)
    at_one = lambda f: field.from_sympy(f.as_expr().subs(t, 1))
    recovered_scalar = at_one(h)/h
    necessary_i4 = recovered_scalar**2*i4-at_one(i4)
    necessary_g8 = recovered_scalar**4*g8-at_one(g8)
    gcd = necessary_i4.numer.gcd(necessary_g8.numer)
    assert gcd.degree() == 1 and gcd.as_expr().subs(t, 1) == 0
    assert at_one(i4) != 0 and at_one(i14) != 0
    assert necessary_i4.denom.as_expr().subs(t, 1) != 0
    assert necessary_g8.denom.as_expr().subs(t, 1) != 0
    numeric = normalized_jet(kind, 2, [1, 2, 3, 4], 1)
    assert value.subs(t, 1) == numeric
    assert covariant_frame(numeric)[-1].det() != 0
    assert recovered_scalar.as_expr().subs(t, 1) == 1
    lam, z, affine, _ = affine_family(kind)
    affine_line = affine.subs(dict(zip(z, [t, 2, 3, 4], strict=True)) | {lam: 2})
    pfaffian = field.from_sympy(s.expand(
        affine_line[0, 1]*affine_line[2, 3]-affine_line[0, 2]*affine_line[1, 3]
        + affine_line[0, 3]*affine_line[1, 2]))
    denominator_ratios = {}
    for name, invariant, exponent in [('I4', i4, 1), ('I14', i14, 3), ('G8', g8, 2)]:
        ratio = field.from_sympy(invariant.denom.as_expr())/pfaffian**exponent
        assert ratio.numer.degree() == ratio.denom.degree() == 0
        denominator_ratios[name] = {'pfaffian_exponent': exponent,
                                    'nonzero_constant_ratio': str(ratio.as_expr())}
    invariants = {}
    for label, f in [('I4', i4), ('I14', i14), ('h', h), ('G8', g8)]:
        invariants[label] = {
            'numerator_degree': f.numer.degree(),
            'denominator_degree': f.denom.degree(),
            'numerator_coefficients_descending': coefficients(f.numer, t),
            'denominator_coefficients_descending': coefficients(f.denom, t)}
    result = {
        'kind': kind, 'component_id': 'T39' if kind == 1 else 'T40',
        'line': 'lambda=2, p=(t,2,3,4), source c=1 before scalar recovery',
        'target': 'normalized_jet(kind,2,[1,2,3,4],1)',
        'invariants': invariants,
        'affine_pfaffian': str(pfaffian.as_expr()),
        'reduced_denominator_relations_on_this_line_only': denominator_ratios,
        'necessary_equations': {
            'I4': 'c(t)^2*I4(t)-I4(1)', 'G8': 'c(t)^4*G8(t)-G8(1)',
            'recovered_scalar': 'h(1)/h(t)',
            'numerator_degrees': [necessary_i4.numer.degree(), necessary_g8.numer.degree()],
            'I4_numerator_coefficients_descending': coefficients(necessary_i4.numer, t),
            'G8_numerator_coefficients_descending': coefficients(necessary_g8.numer, t),
            'monic_numerator_gcd': str(gcd.monic().as_expr())},
        'target_satisfies_all_open_conditions': True,
        'restricted_line_fiber_one_reduced_point': True,
        'five_variable_fiber_degree_computed': False,
        'timings_seconds': {'jet': jet_seconds, 'total': time.monotonic()-started}}
    signal.alarm(0)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = {
        'id': 'ELLIPTIC_FIBER_LINE_PROBE_20261031',
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'proof_note': 'research/notes/ELLIPTIC_FIBER_LINE_PROBE_20261031.md',
        'input_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in [Path(__file__).resolve(),
                                      ROOT/'code/elliptic_six_variable_fiber.py',
                                      ROOT/'code/elliptic_weight1_dominance.py']},
        'all_passed': True, 'families': [probe(1), probe(2)]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print('PASS: both restricted line fibers are one reduced point; '
          'no five-variable or component degree inferred.')


if __name__ == '__main__':
    main()

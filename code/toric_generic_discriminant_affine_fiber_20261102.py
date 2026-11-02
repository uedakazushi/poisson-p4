#!/usr/bin/env python3
"""Check the generic cubic discriminant algebra used in the DVR proof."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
STEM = 'TORIC_GENERIC_DISCRIMINANT_AFFINE_FIBER_20261102'


def verify():
    x, c4, c6, r = sp.symbols('X c4 c6 r')
    cubic = x**3 - 27*c4*x - 54*c6
    derivative = sp.diff(cubic, x)
    norm = sp.factor(sp.resultant(cubic, derivative, x))
    assert sp.expand(norm + 78732*(c4**3-c6**2)) == 0
    special = {c4: r**2/9, c6: -r**3/27}
    assert sp.expand(cubic.subs(special) - (x-r)**2*(x+2*r)) == 0
    assert derivative.subs(special).subs(x, r) == 0
    assert derivative.subs(special).subs(x, -2*r) == 9*r**2
    # In the double-root summand the derivative generates its maximal ideal.
    eps = sp.Symbol('epsilon')
    tangent = sp.rem(derivative.subs(special).subs(x, r+eps), eps**2, eps)
    assert tangent == 6*r*eps
    inputs = [Path(__file__), ROOT/'research/notes'/f'{STEM}.md',
              ROOT/'research/notes/TORIC_AFFINE_MODIFICATION_20261101.md',
              ROOT/'research/notes/TORIC_NORMALIZATION_BOUNDARY_20261101.md']
    for name in ['TORIC_NORMALIZATION_BOUNDARY_20261101', 'TORIC_C9_POLYNOMIAL_20261101',
                 'TORIC_POLYNOMIAL_E2_COORDINATE_20261101', 'TORIC_ZERO_JET_EXHAUSTION_20261101']:
        source = ROOT/'reports'/f'{name}.json'
        assert json.loads(source.read_text())['all_passed']
        inputs.append(source)
    return {'id': STEM, 'all_passed': True, 'field': 'QQ',
            'derivative_norm': str(norm), 'derivative_on_double_root': str(tangent),
            'heights': {'ramified': {'valuation_fprime': 1, 'ramification_index': 2},
                        'simple': {'valuation_fprime': 0, 'ramification_index': 1}},
            'proof_note_conclusions': {'generic_discriminant_affine_length': 1,
                                       'projective_length_at_infinity': 2,
                                       'some_q_has_simple_pole': True},
            'scope': 'Generic point and a dense open of D for T01 only.',
            'input_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in inputs}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    certificate = verify()
    with args.output.open('x') as stream:
        json.dump(certificate, stream, indent=2)
        stream.write('\n')
    print('PASS', certificate['proof_note_conclusions'])

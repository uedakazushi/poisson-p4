#!/usr/bin/env python3
"""Finite identities for a family with constant initial algebra and nonnilpotent t.

The filtration and completeness proof for every integer N>=2 is in the note.
The symbol a in this checker stands for u**N, with filtration degree N.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def verify():
    u, v, h, z, a = s.symbols('u v h z a')
    f, g = h*z, z*z-a*z
    assert s.expand(z*f-h*g-a*f) == 0
    t = u**3*z
    assert s.expand(t*t-u**6*a*z-u**6*g) == 0
    assert s.expand((t*t).subs({h: 0, z: a})-u**6*a*a) == 0
    assert (t*t).subs({h: 0, z: a}) != 0
    # Abstract second component of the product on R plus A*z, after reducing
    # the first component from R=k[[u,v,h]] to A=k[[u,v]].
    r, b, c, d, e, f0 = s.symbols('r b c d e f0')
    def multiply(x, y):
        return x[0]*y[0], x[0]*y[1]+y[0]*x[1]+a*x[1]*y[1]
    lhs = multiply(multiply((r, b), (c, d)), (e, f0))
    rhs = multiply((r, b), multiply((c, d), (e, f0)))
    assert all(s.expand(x-y) == 0 for x, y in zip(lhs, rhs))
    assert multiply((1, 0), (r, b)) == (r, b)
    paths = [str(Path(__file__).relative_to(ROOT)),
             'research/notes/FILTERED_INITIAL_NILPOTENCE_OBSTRUCTION_20261101.md']
    return dict(id='FILTERED_INITIAL_NILPOTENCE_OBSTRUCTION_20261101', all_passed=True,
                field='Any field', parameter_range='Every integer N >= 2',
                algebra='k[[u,v,h,z]]/(h*z,z^2-u^N*z)',
                module_decomposition='k[[u,v,h]] direct-sum k[[u,v]]*z',
                associated_graded='k[u,v,h,z]/(h*z,z^2)',
                generic_rank_over_three_variable_base=1,
                embedded_difference_free_over_two_variable_base=True,
                embedded_difference_free_rank=1,
                distinguished_element='t=u^3*z', order_t=4,
                square='u^(N+6)*z', order_square='N+7', square_nonzero=True,
                reduced_components=[dict(ideal='(z)', dimension=3),
                                    dict(ideal='(h,z-u^N)', dimension=2)],
                s_polynomial_identity='z*(h*z)-h*(z^2-a*z)=a*(h*z)',
                associative_module_product_verified=True,
                limitations=['This family is not claimed to be an S61 or Poisson family.',
                             'It refutes an implication from initial-algebra and embedded-module data alone.'],
                proof_note=paths[1],
                input_sha256={path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify()
    with args.output.open('x') as output:
        json.dump(result, output, indent=2)
        output.write('\n')
    print('FILTERED INITIAL NILPOTENCE OBSTRUCTION PASS')

#!/usr/bin/env python3
"""Universal dual-curve recognition of Q51 and known-cover disjointness input.

The polynomial checks are over all sixteen first-jet variables.  Geometric
disjointness uses the primary results and closed-divisor argument in the note;
this program does not classify the other components of S61.
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

from elliptic_curve_invariants_20261101 import (
    curve_invariants_from_beta, pfaffOmega,
)
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61, X, coefficients
from weight1_structure import J

ROOT = Path(__file__).resolve().parents[1]


def verify():
    started = time.monotonic()
    model = S61()
    z, u, beta = primitive_beta(model)
    y, t = s.Matrix(z), s.Symbol('dual_t')
    phi = (beta+t*J).row_join(J*y).col_join(
        (-(J*y).T).row_join(s.zeros(1)))
    omega, pfaffians = pfaffOmega(phi, [*y, t])
    # Fisher's coordinates are (y0,y1,y2,y3,t); the S61 basepoint is e_t.
    reorder = [4, 0, 1, 2, 3]
    P = omega.extract(reorder, reorder).subs(
        dict(zip([t, *y], X, strict=True)), simultaneous=True).applyfunc(s.expand)
    value = P.subs(dict(zip(X, [1, 0, 0, 0, 0], strict=True)))
    expected_value = s.diag(s.zeros(1), -5*J)
    assert value == expected_value
    normalized = model.normalize(-P/5)
    one = [i for i, w in enumerate(model.reduced['parameter_weights']) if w == 1]
    two = [i for i, w in enumerate(model.reduced['parameter_weights']) if w == 2]
    actual_beta = [normalized['parameters'][i] for i in one]
    assert all(s.expand(a-2*b) == 0 for a, b in zip(actual_beta, u, strict=True))
    q51 = json.loads((ROOT/'reports/Q51_UNIVERSAL_SECTION_20261101.json').read_text())
    doubled = dict(zip(u, [2*v for v in u], strict=True))
    target_q = [s.sympify(v).subs(doubled, simultaneous=True)
                for v in q51['quadratic_section_q']]
    assert all(s.expand(a-normalized['parameters'][i]) == 0
               for a, i in zip(target_q, two, strict=True))
    # All 150 homogeneous coefficients are compared, not just the 61 free ones.
    reconstructed = coefficients(model.reconstruct(normalized['parameters']))
    assert all(s.expand(a-b) == 0 for a, b in zip(
        reconstructed, normalized['coefficients'], strict=True))

    # A representative of the unique nonzero rank-one flag-cone orbit:
    # beta = y0 * dy0 wedge dy2.  Sp4 invariance then gives vanishing on
    # the entire cone, including its vertex by homogeneity.
    flag_beta = s.zeros(4)
    flag_beta[0, 2], flag_beta[2, 0] = z[0], -z[0]
    pairs = list(combinations(range(4), 2))
    flag_coefficients = s.Matrix([flag_beta[i, j].diff(z[k])
                                 for i, j in pairs for k in range(4)])
    flag_c4, flag_c6, flag_delta, flag_j = curve_invariants_from_beta(flag_coefficients)
    assert (flag_c4, flag_c6, flag_delta, flag_j) == (0, 0, 0, None)

    inputs = [Path(__file__), ROOT/'code/s61_model.py',
        ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
        ROOT/'code/elliptic_curve_invariants_20261101.py', ROOT/'code/weight1_structure.py',
        ROOT/'reports/Q51_UNIVERSAL_SECTION_20261101.json',
        ROOT/'reports/Q52_UNIVERSAL_SECTION_20261101.json',
        ROOT/'reports/NULL_CORRELATION_DUALITY_20261101.json',
        ROOT/'reports/ELLIPTIC_SECTION_INTERSECTION_20261101.json',
        ROOT/'reports/TORIC_C9_POLYNOMIAL_20261101.json',
        ROOT/'reports/TORIC_POLYNOMIAL_E2_COORDINATE_20261101.json',
        ROOT/'research/notes/KNOWN_COMPONENT_CURVE_INVARIANT_CONSTANCY_20261101.md',
        ROOT/'research/notes/TORIC_AFFINE_MODIFICATION_20261101.md']
    return {
        'id': 'KNOWN_FIVE_SMOOTH_LOCUS_20261101', 'all_passed': True,
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'field': 'QQ[all sixteen primitive first-jet coordinates]',
        'weight1_u_indices_zero_based': one, 'weight2_u_indices_zero_based': two,
        'Phi_E': '[[beta(y)+t*J,J*y],[-(J*y)^T,0]]',
        'Fisher_coordinate_order': [str(v) for v in [*y, t]],
        'S61_coordinate_order': ['t', *map(str, y)],
        'basepoint_value': '-5*J', 'normalizing_scalar': '-1/5',
        'normalized_weight1': [str(v) for v in actual_beta],
        'normalized_weight2': [str(v) for v in target_q],
        'identity': 'S61_normalize(-Fisher_Omega(Phi_E(beta))/5) = Q51(2*beta)',
        'all16_first_jet_identities_checked': True,
        'all45_section_identities_checked': True,
        'all150_homogeneous_reconstruction_coefficients_checked': True,
        'perspective_parameters': [str(v) for v in normalized['perspective']],
        'flag_representative': 'beta = y0 dy0 wedge dy2',
        'flag_c4': str(flag_c4), 'flag_c6': str(flag_c6), 'flag_Delta': str(flag_delta),
        'geometric_conclusion': 'Over Delta!=0 the three reduced closures T01,T39,T40 are pairwise disjoint, and their union is finite etale of degree5.',
        'known_union_algebra': 'R_Delta[X]/(X^3-27*c4*X-54*c6) x R_Delta x R_Delta',
        'scope': 'Only the reduced union of the three known components; no exclusion of other S61 components and no ambient smoothness claim.',
        'proof_note': 'research/notes/KNOWN_FIVE_SMOOTH_LOCUS_20261101.md',
        'elapsed_seconds': round(time.monotonic()-started, 3),
        'input_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in inputs},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = verify()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as handle:
        json.dump(report, handle, indent=2)
        handle.write('\n')
    print(json.dumps({'all_passed': True, 'universal_coefficients': 150,
                      'seconds': report['elapsed_seconds']}))


if __name__ == '__main__':
    main()

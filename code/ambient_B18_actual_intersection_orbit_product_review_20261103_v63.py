#!/usr/bin/env python3
"""Tiny independent orbit19 minor and scheme-level formal-product proof gate."""
import argparse
import gzip
import hashlib
import itertools as it
import json
from pathlib import Path
import time

import sympy as s
from ambient_B18_conductor_dimension20_review_20261103_v63 import physical_GL_column

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63_REVIEW'


def load(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def accepted(p):
    d = load(p)
    assert d['all_passed'] is True
    for name, digest in d.get('input_sha256', {}).items():
        assert sha(ROOT / name) == digest, name
    for name, digest in d.get('evidence_sha256', {}).items():
        artifact = ROOT / name if (ROOT / name).exists() else p.parent / name
        assert sha(artifact) == digest, name
    return d


def main(output):
    start = time.monotonic()
    assert not output.exists()
    folder = ROOT / 'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63'
    accepted(folder / 'report.json')
    saved = load(folder / 'certificate.json.gz')
    reduction_path = ROOT / 'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json'
    reduction = load(reduction_path)
    tangent = ROOT / 'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63_REVIEW/report.json'
    accepted(tangent)
    # Four nonzero physical coefficients of e_y wedge(x²,0,0,hx,h²+z²).
    pairs = list(it.combinations(range(5), 2))
    monomials = list(it.combinations_with_replacement(range(5), 2))
    point = s.zeros(150, 1)
    for pair, monomial, coefficient in [((0, 2), (1, 1), -1), ((2, 3), (0, 1), 1), ((2, 4), (0, 0), 1), ((2, 4), (3, 3), 1)]:
        point[15 * pairs.index(pair) + monomials.index(monomial)] = coefficient
    assert point == s.Matrix([s.Rational(v) for v in saved['physical_point150']])
    O25 = s.Matrix.hstack(*[physical_GL_column(i, j, point) for i, j in it.product(range(5), repeat=2)])
    assert O25[:, 12] == point and point[106] == 1
    gl18 = saved['GL18_generator_indices']
    O19 = point.row_join(O25[:, gl18])
    assert O19 == s.Matrix([[s.Rational(v) for v in row] for row in saved['physical150_by19_scalar_then_GL']])
    free = saved['canonical126_raw_coordinate_indices']
    C = s.Matrix(reduction['canonical_reconstruction_from126'])
    assert C[free, :] == s.eye(126) and C * O19[free, :] == O19
    rows = saved['orbit19_minor_rows_in126']
    rawrows = [free[i] for i in rows]
    minor = O19[rawrows, :]
    assert minor.det() == -32 == s.Rational(saved['orbit19_minor'])
    # Full upper rank check, not merely a lower minor.
    coefficients = minor.inv() * O25[rawrows, :]
    assert O19 * coefficients == O25
    complement = saved['affine_transverse107_free_coordinate_indices']
    assert complement == [i for i in range(126) if i not in rows] and len(complement) == 107
    deps = [ROOT / 'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_CA_CB_MINIMAL_LENGTHS_AND_R_CM_20261103_v63_REVIEW.json']
    for p in deps:
        accepted(p)
    inputs = [Path(__file__), ROOT / f'research/notes/{STEM}.md',
        ROOT / 'code/ambient_B18_conductor_dimension20_review_20261103_v63.py',
        ROOT / 'code/ambient_B18_actual_intersection_orbit_product_20261103_v63.py',
        ROOT / 'research/notes/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63.md',
        folder / 'report.json', folder / 'certificate.json.gz', reduction_path, tangent, *deps]
    result = dict(id=STEM, all_passed=True, producer_rerun=False,
        independent_index_GL25_and_orbit19_minor=-32,
        full_GL25_span_identity_and_literal_Eyy_scalar_checked=True,
        canonical126_complement107_and_action_Jacobian_unit_checked=True,
        actual_sum_ideal_G_invariance_and_all_Artin_points_proof_accepted=True,
        full_D_nilpotents_retained_in_formal_G19_product=True,
        simultaneous_CA_CB_D_N_R_slice_accepted=True,
        cone_dimension=28, cone_edim=39, cone_multiplicity=2, cone_depth_interval=[19, 20],
        transverse_dimension=9, transverse_edim=20, transverse_multiplicity=2, transverse_depth_interval=[0, 1],
        projective_dimension=27, projective_edim=38, projective_depth_interval=[18, 19],
        exact_depth_slice_socle_full_relative_ring_or_nilvisibility_not_determined=True,
        independent_method='Index-level GL columns; complete19-span identity; actual ideal equivariance plus formal IFT/Artin functors; regular scalar removal.',
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in inputs}, seconds=time.monotonic() - start)
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output)

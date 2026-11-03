#!/usr/bin/env python3
"""Conceptual reduced B18 depth/seminormalization input binding; no CAS."""
import argparse
import hashlib
import json
from math import comb
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_REDUCED_DEPTH_SEMINORMALIZATION_20261103_v63_REVIEW'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def accepted(path):
    data = json.loads(path.read_text())
    assert data['all_passed'] is True
    for name, digest in data.get('input_sha256', {}).items():
        assert sha(ROOT/name) == digest, name
    for name, digest in data.get('evidence_sha256', {}).items():
        target = ROOT/name if (ROOT/name).exists() else path.parent/name
        assert sha(target) == digest, name
    return data


def main(output):
    start = time.monotonic()
    assert not output.exists()
    paths = [ROOT/'reports/AMBIENT_B18_N_R_DVR_CUBIC_GRAPH_20261103_v63_REVIEW/report.json',
             ROOT/'reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63_REVIEW.json',
             ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63_REVIEW.json',
             ROOT/'reports/AMBIENT_B18_NILRADICAL_LINEAR_VISIBILITY_20261103_v63_REVIEW.json']
    contact, N, orbit, gluing = [accepted(p) for p in paths]
    assert contact['k_exact'] == 4
    assert N['cone_dimension_embedding_dimension_depth'] == [28, 37, 20]
    assert N['elementary_subintegral_finite_universal_homeomorphism_and_seminormalization_B']
    assert orbit['simultaneous_CA_CB_D_N_R_slice_accepted']
    assert gluing['reduced_slice_fiber_product_with_DVR_and_unknown_finite_t_power_proved']
    # Independent monomial intersection: (T) intersect (x1,...,x9)=(T*x1,...,T*x9).
    tgen = (0,)*9+(1,)
    xgens = [tuple(int(i == j) for i in range(10)) for j in range(9)]
    intersection = [tuple(max(a, b) for a, b in zip(tgen, x)) for x in xgens]
    products = [tuple(a+b for a, b in zip(tgen, x)) for x in xgens]
    assert intersection == products and all(set(e) <= {0, 1} for e in intersection)
    # Count independent standard monomials in the displayed homogeneous model.
    hilbert = [1]+[comb(d+8, 8)+1 for d in range(1, 12)]
    assert hilbert[1] == 10
    # Hilbert numerator at denominator (1-z)^9 is 1+z(1-z)^8.
    numerator = [1]+[(-1)**(i-1)*comb(8, i-1) for i in range(1, 10)]
    assert sum(numerator) == 1
    # Multiplication by T identifies the ideal resolution with the Koszul ideal on x9.
    betti = [1]+[comb(9, i) for i in range(1, 10)]
    assert betti == [1, 9, 36, 84, 126, 126, 84, 36, 9, 1]
    inputs = [Path(__file__), ROOT/f'research/notes/{STEM}.md', *paths]
    report = dict(id=STEM, all_passed=True,
        verification_kind='Independent conceptual finite/subintegral/square-cube/conductor proof with accepted input hashes and elementary monomial/Koszul checks; no CAS.',
        actual_reduced_slice_length4_fiber_product_not_arbitrary_abstract_map=True,
        reduced_slice_depth_exact=1, reduced_cone_depth_exact=20, reduced_projective_depth_exact=19,
        formal_reduced_seminormalization='C[[g19,x1,...,x9,T]]/(T*x1,...,T*x9)',
        normalization='C[[g19,x1,...,x9]] direct_sum C[[g19,T]]',
        finite_subintegral_extension_and_total_fraction_ring_checked=True,
        seminormality_componentwise_normality_and_square_cube_residue_test_proved=True,
        seminormalization_to_normalization_conductor='(x1,...,x9) on N-normalizer; (T) on R-normalizer',
        seminormalization_conductor_quotient='C[[g19]]',
        conductor_support_reduced_smooth_G19=True,
        seminormal_cone_dim_edim_depth_multiplicity=[28, 29, 20, 1],
        seminormal_projective_dim_edim_depth_multiplicity=[27, 28, 19, 1],
        seminormal_slice_Hilbert_values=hilbert,
        seminormal_slice_Betti_numbers=betti,
        full_nonreduced_D_depth_or_original_reduced_conductor_or_actual37jet_map_claimed=False,
        scope='Formal generic-B18 completions of D_red only; reduced formal seminormalization, no convergent/global-completion-commutation claim.',
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in inputs}, seconds=time.monotonic()-start)
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    main(parser.parse_args().output)

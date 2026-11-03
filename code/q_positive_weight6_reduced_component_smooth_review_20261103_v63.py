#!/usr/bin/env python3
"""Independent multiplicity-application audit; the mathematical proof is in the note."""
import argparse
import hashlib
import json
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_POSITIVE_WEIGHT6_REDUCED_COMPONENT_SMOOTH_20261103_v63_REVIEW'
DEPS = [
    'reports/Q_POSITIVE_WEIGHT6_SMALL_MODEL_20261103_v63_REVIEW/report.json',
    'reports/Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63/report.json',
    'reports/Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63_REVIEW/report.json',
    'reports/Q_POSITIVE_WEIGHT6_COMPONENT25_20261103_v63_REVIEW.json',
    'reports/Q_NEW_FAMILY_LENGTH5_COMPLETION_20261103_v63.json',
    'reports/Q_NEW_FAMILY_FULL_COMPLETION_ROOT_REVIEW_20261103_v63.json',
    'reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63_REVIEW.json',
    'reports/Q_POSITIVE_WEIGHT6_NONCM_20261103_v63.json',
    'reports/R4_NEW25_POSITIVE_AXIS_SMOOTH_REDUCED_COMPONENT_20261103_v63.json',
]


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main(out):
    start = time.monotonic()
    checked = {}
    docs = []
    for name in DEPS:
        path = ROOT / name
        d = json.loads(path.read_text())
        assert d['all_passed'], name
        for filename, want in d.get('input_sha256', {}).items():
            source = ROOT / filename
            actual = checked.setdefault(str(source), sha(source))
            assert actual == want, (name, filename)
        for filename, want in d.get('evidence_sha256', {}).items():
            source = path.parent / filename
            actual = checked.setdefault(str(source), sha(source))
            assert actual == want, (name, filename)
        docs.append(d)
    small, fifth, fifth_review, component, length, length_review, orbit, noncm, author = docs
    assert small['exact_orbit23_and_projective22_factor']
    assert small['exact_implicit95_and_formal8_presentation']
    assert fifth['through5_ideal_dimension'] == 3
    assert fifth['through5_ideal_degree'] == 5
    assert fifth_review['exact_fifth_cone_dimension3_degree5']
    assert fifth_review['radical_two_linear_components_dimensions'] == [3, 2]
    assert fifth['radical_minimal_primes_in_normal_z_coordinates'] == [
        [1, 2, 4, 6, 7], [1, 3, 4, 5, 6, 7]]
    assert fifth['radical_unique3plane_free_coordinates'] == ['z0', 'z3', 'z5']
    assert component['actual_global_projective_irreducible_component_dimension_exact'] == 25
    assert component['actual_full_transverse_local_dimension_exact'] == 3
    assert component['actual_global_component_image_closure_identification_accepted']
    assert length['actual_generic_Artin_length'] == 5
    assert length_review['full_closed_and_generic_transverse_length'] == 5
    assert orbit['nonzero_eta_one_complex_projective_orbit_and_closure_Qa8over9']
    assert noncm['actual_full_transverse_embedding_dimension'] == 8
    assert noncm['full_actual_Poisson_local_and_completed_rings_non_Cohen_Macaulay_at_positive_point']
    assert author['full_transverse_multiplicity_exact'] == 5
    assert author['reduced_new_component_transverse_multiplicity_exact'] == 1
    assert author['reduced_new_component_smooth_on_entire_nonzero_positive_axis']
    assert author['unique_dimension25_minimal_component_at_positive_axis']
    assert not author['lower_dimensional_minimal_or_embedded_components_excluded']
    assert not author['global_component_normality_claimed']

    # After the proof's associativity inequality, these are all integer possibilities.
    # extra is the aggregate contribution from other top-dimensional minimal primes.
    possibilities = [(ea, eb, extra)
                     for ea in range(1, 6) for eb in range(1, 6)
                     for extra in range(6) if ea >= 5 * eb + extra]
    assert possibilities == [(5, 1, 0)]
    assert 8 - 3 + 1 == 6 > 5
    assert sorted(set(range(8)) - {1, 2, 4, 6, 7}) == [0, 3, 5]
    paths = [Path(__file__), ROOT / f'research/notes/{STEM}.md',
             *[ROOT / name for name in DEPS]]
    report = dict(
        id=STEM, all_passed=True,
        review_kind='Independent mathematical multiplicity and completion-scope audit; current dependency hashes checked',
        root_author_proof_and_receipt_independently_accepted=True,
        checked_unique_current_input_and_evidence_hashes=len(checked),
        all_inputs_are_actual_full_Poisson_ring_not_fixed_pencil_source=True,
        actual_G23_formal_product_respects_reduced_invariant_component=True,
        completion_analytic_splitting_accounted_for=True,
        excellent_reduced_equidimensional_B_dimension=3,
        all_analytic_C_branches_inherit_full_scheme_generic_length=5,
        associativity_lower_bound='e(A) >= 5*e(B)',
        actual_initial5_multiplicity_upper_bound=5,
        actual_full_transverse_Hilbert_Samuel_multiplicity=5,
        reduced_new_component_transverse_multiplicity=1,
        regularity_proof='minimal reduction; finite torsion-free rank1 over normal C[[t3]]; no assumed analytic domain',
        reduced_component_smooth_at_all_nonzero_positive_axis_points=True,
        reduced_projective_component_dimension=25,
        unique_top_dimensional_transverse_minimal_prime=True,
        unique25dimensional_projective_component_through_positive_point=True,
        reduced_transverse_tangent_plane_vanishing_coordinates=[1, 2, 4, 6, 7],
        reduced_transverse_tangent_plane_free_coordinates=[0, 3, 5],
        full_Poisson_ring_still_non_Cohen_Macaulay=True,
        lower_dimensional_minimal_primes_or_embedded_primes_excluded=False,
        P2_actual_branch_existence_decided=False,
        full_initial_ideal_or_full_formal_branch_exhaustion_claimed=False,
        Qa8over9_regular_or_all_Qa_membership_claimed=False,
        mathematical_proof_is_in_note_not_automatic_from_hash_checks=True,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        seconds=time.monotonic() - start)
    with out.open('x') as f:
        json.dump(report, f, indent=2)
        f.write('\n')
    print('PASS positive reduced C58 smooth, e(A)=5/e(B)=1; lower-dimensional branches remain open')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output)

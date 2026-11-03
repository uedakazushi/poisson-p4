#!/usr/bin/env python3
"""Input binding for the conceptual DVR-union/linear-nil existence proof; no CAS."""
import argparse
import hashlib
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_NILRADICAL_LINEAR_VISIBILITY_20261103_v63_REVIEW'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def accepted(path):
    data = json.loads(path.read_text())
    assert data['all_passed'] is True
    for name, digest in data.get('input_sha256', {}).items():
        assert sha(ROOT / name) == digest, name
    for name, digest in data.get('evidence_sha256', {}).items():
        artifact = ROOT / name if (ROOT / name).exists() else path.parent / name
        assert sha(artifact) == digest, name
    return data


def main(output):
    start = time.monotonic()
    assert not output.exists()
    paths = [ROOT / 'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63_REVIEW.json']
    orbit, tangent, N, R, union = [accepted(path) for path in paths]
    assert orbit['full_D_nilpotents_retained_in_formal_G19_product']
    assert orbit['simultaneous_CA_CB_D_N_R_slice_accepted']
    assert orbit['transverse_edim'] == 20
    assert tangent['actual_D_cone_dim28_edim39_e2_depth_atmost20_accepted']
    assert N['cone_dimension_embedding_dimension_depth'][:2] == [28, 37]
    assert N['actual_N_domain_and_no_embedded_associated_primes']
    assert R['entire_tzero_nilpotent_Jordanflag_cover_and_reduced_R_intersection_N_equals_irreducible_B18']
    assert R['proper_finite_Nakayama_target_smoothness_at_generic_boundary_and_FULL_R_R1_accepted']
    assert 37 - 19 == 18 and 20 - (18 + 1) == 1
    inputs = [Path(__file__), ROOT / f'research/notes/{STEM}.md', *paths]
    report = dict(id=STEM, all_passed=True,
        verification_kind='Independent conceptual proof with accepted input hashes; no CAS or source computation.',
        abstract_DVR_fiber_product_principal_kernel_lemma_proved=True,
        both_fiber_product_projections_surjective_used=True,
        maximal_ideal_generated_by_lifted18_plus_one_kernel_element=True,
        actual_R_regular_on_open_B18_orbit_from_accepted_generic_smoothness=True,
        reduced_slice_fiber_product_with_DVR_and_unknown_finite_t_power_proved=True,
        reduced_transverse_embedding_dimension_upper_bound=19,
        actual_transverse_embedding_dimension=20,
        linear_nil_cotangent_kernel_dimension_lower_bound=1,
        actual_nilradical_not_contained_in_maximal_ideal_squared=True,
        completed_reduced_D_embedding_dimension_upper_bound=38,
        excellence_reduction_completion_and_actual_local_descent_checked=True,
        exact_Dred_embedding_dimension_or_explicit_nil_generator_or_exact_depth_claim=False,
        input_sha256={str(path.relative_to(ROOT)): sha(path) for path in inputs}, seconds=time.monotonic() - start)
    output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)

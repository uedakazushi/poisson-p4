#!/usr/bin/env python3
"""Gate-only application: target/source equality and exact spectral degrees.

No invariant evaluation, CAS, new point, Krylov replay or source expansion.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_NEW_FAMILY_TARGET_FIELD_AND_SPECTRAL_DEGREES_20261104_v63'
REVIEW = 'Q_NEW_FAMILY_ZN_PRIMITIVE_KRYLOV_20261104_v63_REVIEW'
GATES = {
    'source_constants_and_known_finite_involution':
        'Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW',
    'whole_fixed_source_coordinates':
        'Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63_REVIEW',
    'rational_q_B_inverse_and_actual_normalized_maps':
        'Q_NEW_FAMILY_TWO_SHAPE_RECOVERY_20261104_v63_REVIEW',
    'spectral_base_independence':
        'Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63_REVIEW',
    'shape_degree10_and_content_maximality':
        'Q_NEW_FAMILY_TWO_SHAPE_CONTENT_PRIMITIVITY_20261104_v63_REVIEW',
    'shape_geometric_birationality_by_coprime_degrees3_5':
        'Q_NEW_FAMILY_TWO_SHAPE_COPRIME_DEGREES_20261104_v63_REVIEW',
    'all_degree_source_homogeneity':
        'Q_NEW_FAMILY_GENERAL_SOURCE_SCALE_20261104_v63_REVIEW',
    'equation_degree_at_least26':
        'Q_NEW_FAMILY_FULL_DEGREE24_SOURCE_RELATION_20261104_v63_REVIEW',
    'ZN_primitive_independent_gate': REVIEW,
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(output, allow_pending):
    start = time.monotonic()
    status, hashes, pending, reports = {}, {}, [], {}
    for key, stem in GATES.items():
        path = ROOT / 'reports' / stem / 'report.json'
        if not path.exists():
            pending.append(str(path.relative_to(ROOT)))
            status[key] = False
            continue
        data = json.loads(path.read_text())
        assert data['all_passed'] is True, path
        status[key] = True
        reports[key] = data
        hashes[str(path.relative_to(ROOT))] = sha(path)
    if pending and not allow_pending:
        raise RuntimeError('Independent prerequisite still pending: ' + ', '.join(pending))
    if 'shape_degree10_and_content_maximality' in reports:
        data = reports['shape_degree10_and_content_maximality']
        assert data['exact_source_SHAPE_extension_degree'] == 10
        assert data['primitive_G_irreducible_over_QQ_R_Z0_Z1_Z2_in_r']
    if 'all_degree_source_homogeneity' in reports:
        data = reports['all_degree_source_homogeneity']
        assert data['all_degrees_invariant_source_homogeneity_c_power_D_over2_proved']
        assert data['degree14_form_xi7_times_shape_and_nonvanishing_by_lower26_accepted']
    if 'spectral_base_independence' in reports:
        assert reports['spectral_base_independence'][
            'NEW25_dominates_A5_pencil_and_full_irreducible_characteristic_base']
    if 'ZN_primitive_independent_gate' in reports:
        data = reports['ZN_primitive_independent_gate']
        assert data['ZN_primitive_degree10_independently_proved']
        assert data['independent_Bareiss_determinant'] == 21349
    author = ROOT / 'reports/Q_NEW_FAMILY_ZN_PRIMITIVE_KRYLOV_20261104_v63/report.json'
    data = json.loads(author.read_text())
    assert data['all_passed'] is True and data['ZN_primitive_degree10_proved'] is True
    assert data['Krylov_determinant_mod32003'] == 21349
    hashes[str(author.relative_to(ROOT))] = sha(author)
    note = ROOT / 'research/notes' / (STEM + '.md')
    hashes[str(note.relative_to(ROOT))] = sha(note)
    hashes[str(Path(__file__).relative_to(ROOT))] = sha(Path(__file__))
    complete = not pending
    output.mkdir(parents=True, exist_ok=False)
    report = {
        'id': STEM,
        'all_passed': complete,
        'status': 'proved_by_accepted_prerequisites' if complete else 'conditional_pending_peer_gate',
        'gates': status,
        'pending': pending,
        'method': 'exact field-inclusion and valuation proof; prerequisite/hash checks only',
        'target': 'primitive first-jet hypersurface B_NEW of dimension15 in W16; not the entire Poisson component25',
        'coefficient_fields': ['QQ', 'C'],
        'actual_target_field_equals_known_tau_fixed_source_field': complete,
        'actual_target_field_is_pure_rational_trdeg5': complete,
        'QQ_rationality_is_for_B_NEW_not_ambient_W16_or_entire_component': True,
        'known_tau_fixed_source_to_actual_target_index': 1 if complete else None,
        'projective_target_to_normalized_spectral_base_degree': 10 if complete else None,
        'affine_target_to_F4_and_four_ratios_degree': 20 if complete else None,
        'affine_target_to_five_octic_pencil_coefficients_degree': 40 if complete else None,
        'proof_obligations': {
            'E_in_L_is_actual_dominant_source_pullback_not_blanket_invariant_fraction_claim': True,
            'Kshape_in_E_from_actual_polynomial_invariant_ratios_and_ZN_primitivity': True,
            'J14_nonzero_on_NEW_from_prime_equation_degree_lower26': True,
            'xi_recovery_uses_J14_over_F4_cubed_and_nonzero_shape_coefficients': True,
            'complex_base_extension_uses_geometric_curve_and_preserved_Bezout_not_QQ_irreducibility_alone': True,
            'degree2_xi_square_extension_by_F4_adic_odd_valuation': True,
            'degree2_F4_over_octic_base_by_linear_polynomial_odd_valuation': True,
            'projective_field_is_scalar_fixed_field_xi_weight2': True,
            'QQ_source_and_actual_tau_identity_and_Lie_invariant_kernel_descent': True,
            'geometric_firstjet_image_defined_over_QQ_and_J14_nonzero_after_complexification': True,
        },
        'generic_function_field_degrees_not_global_finite_flat_morphisms': True,
        'curve_genus_actual_global_component_rationality_or_defining_equation_claimed': False,
        'no_new_CAS_point_prime_invariant_expansion_source26_or_Krylov_rerun': True,
        'input_sha256': hashes,
        'seconds': time.monotonic() - start,
    }
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'all_passed': complete, 'pending': pending,
                      'seconds': report['seconds']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--allow-pending', action='store_true')
    args = parser.parse_args()
    main(args.output, args.allow_pending)

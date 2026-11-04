#!/usr/bin/env python3
"""Finite coefficient schema for exact B18 depth; no coefficient/CAS rerun."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_FOUR_MODULE_DEPTH_CRITERION_20261104_v63'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def main(output):
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    signal.alarm(60)
    relative = ROOT / 'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63'
    data, cert = read(relative / 'report.json'), read(relative / 'certificate.json')
    assert data['all_frozen_source_weight_and_sparsity_checks_passed']
    assert sha(relative / 'certificate.json') == data['certificate_sha256']
    for name, digest in data['input_sha256'].items():
        assert sha(ROOT / name) == digest, name
    peer = ROOT / 'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63_REVIEW/report.json'
    assert read(peer)['all_passed'] is True
    orbit = ROOT / 'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63_REVIEW.json'
    assert read(orbit)['all_passed'] is True
    gauge = ROOT / 'reports/AMBIENT_B18_CORE_GAUGE_PREFLIGHT_20261104_v63'
    assert read(gauge / 'report.json')['rank2_normalizer_algebraic_power_series_recipe_proved']
    implicit = read(gauge / 'implicit-system.json')
    base = cert['remaining_actual_even_base10_raw_coordinates']
    positive = cert['positive_even_base_raw_coordinate']
    zero = [i for i in base if i != positive]
    assert positive == 140
    assert zero == [27, 108, 109, 112, 113, 117, 118, 119, 61]
    assert implicit['target_parameter_names'] == [f'x{i}' for i in zero]
    assert cert['whole_D_bound_positive_degree'] == 3
    assert cert['whole_D_number_of_generators_upper_bound'] == 285
    assert cert['graph_discrepancies_retained_all87']
    assert cert['eliminated87_ambient_graph_weights'] == {'0': 4, '1': 53, '2': 24, '3': 6}
    cubics = [rec['raw'] for rec in cert['possible_cubic_discrepancy_raw_coordinates']]
    assert cubics == [9, 39, 54, 84, 99, 144]
    for family in ('CA', 'CB'):
        weights = cert['component_records'][family]['actual_odd_minimal10_weights']
        assert weights.count(0) == 8 and weights.count(1) == 2
        assert cert['component_records'][family]['tracefree_rank2_generator_weight'] == 0
    # Exact monomial example: C0 has exponent semigroup {0,2,3,...}.
    # Multiplying its generator ideal (tau) by C0 omits exactly exponent2.
    def C0_exponent(e):
        return e == 0 or e >= 2
    def Q3_exponent(e):
        return e >= 1 and C0_exponent(e - 1)
    def Q4_exponent(e):
        return e >= 1
    assert C0_exponent(2) and not Q3_exponent(2) and Q4_exponent(2)
    schema = {
        'regular_zero_base9_actual_rows': zero,
        'positive_even_base_raw_coordinate': positive,
        'algebras': {'C0': 'R0+I0*tau', 'H': 'R0+(I0+J0)*tau',
                     'normalizer_relation': 'tau^2=nu',
                     'pair_product': '(a,b)*(c,d)=(a*c+nu*b*d,a*d+b*c)'},
        'generator_normalization': 'q_(d,alpha)=w^d*c_(d,alpha), c0 in C0, c_(d>=1) in H',
        'L_d': 'C0-span of all normalized degree-d generators; L_d=0 for d>3',
        'Jinf': 'sum_(d=0..3) H*L_d inside H subset R0^2',
        'Qn': {'0': 'L0', '1': 'H*L0+L1', '2': 'H*L0+H*L1+L2',
               '3': 'H*L0+H*L1+H*L2+L3', 'n>=4': 'Jinf'},
        'Tn': {'0': '(C0 intersect Jinf)/L0', '1': 'Jinf/Q[1]',
               '2': 'Jinf/Q[2]', '3': 'Jinf/Q[3]', 'n>=4': '0'},
        'depth0_iff': 'one of the four finite R0-modules Tn has nonzero m_R0 socle',
        'socle0': '((L0:m0) intersect C0 intersect Jinf)/L0',
        'socle_n1to3': '((Q[n]:m0) intersect Jinf)/Q[n]',
        'top_T3_at_most12_generators': {
            'raw_cubic_graph_rows': cubics,
            'positive_odd_coefficient_generators': ['j1', 'j2'],
            'pairs': ['(nu*j_e*b_alpha,j_e*a_alpha) for e=1,2; alpha=1..6'],
            'quotient': '(H*L0+H*L1+H*L2+H*L3)/(H*L0+H*L1+H*L2+L3)',
        },
        'actual_coefficients_needed': ['nu', 'eight I0 generators', 'two J0 generators',
                                       'all normalized actual Q generator pairs, including all87 graphs'],
        'actual_coefficients_computed_here': False,
        'sharpness_example': {'R0': 'k[[nu]]', 'tau_square': 'nu', 'I0': '(nu)',
                              'J0': 'R0', 'Q': '(w^3*tau)',
                              'w3_tau2_nonzero': True, 'w4_tau2_zero': True},
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / 'criterion-schema.json').write_text(json.dumps(schema, indent=2) + '\n')
    paths = [Path(__file__), ROOT / 'research/notes' / (STEM + '.md'),
             relative / 'report.json', relative / 'certificate.json', peer, orbit,
             gauge / 'report.json', gauge / 'implicit-system.json']
    report = {
        'id': STEM, 'all_passed': True,
        'exact_structural_reduction_not_actual_coefficient_solution': True,
        'actual_CA_component_embedding_and_all87_relative_graphs_retained': True,
        'homogeneous_ideal_generator_degree_upper_bound': 3,
        'positive_components_stable_from_degree': 4,
        'w_torsion_nonzero_possible_degrees': [0, 1, 2, 3],
        'uniform_w_torsion_annihilating_power': 4,
        'finite_zero_regular_base_dimension': 9,
        'number_of_finite_module_socle_tests': 4,
        'all_module_subobjects_embed_in_R0_squared': True,
        'top_weight3_torsion_generator_bound': 12,
        'cone_depth19_or20_decided': False,
        'full_nilradical_computed_or_w_torsion_identified_with_nil': False,
        'no_new_CAS_GB_ungauged_core_or_source_rerun_points_or_coefficients': True,
        'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in paths},
        'evidence_sha256': {'criterion-schema.json': sha(output / 'criterion-schema.json')},
        'seconds': time.monotonic() - start,
        'max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ['all_passed', 'seconds', 'max_rss_kib']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)

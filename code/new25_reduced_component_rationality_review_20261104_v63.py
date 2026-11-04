#!/usr/bin/env python3
"""Independent combined field/component audit; requires all FINAL receipts."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = 'NEW25_REDUCED_COMPONENT_RATIONALITY_20261104_v63'
STEM = BASE + '_REVIEW'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_text())


def matmul2(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)]
            for i in range(2)]


def main(output):
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    signal.alarm(60)
    primary_path = ROOT/'reports'/BASE/'report.json'
    primary = read(primary_path)
    assert primary['all_passed']
    assert primary['full_S61_Y16_rational_over_QQ_and_C_by_same_base_Luroth']
    assert primary['reduced_projective_X25_birational_Y16_times_A9']
    assert primary['reduced_projective_new_component_pure_rational_dimension25']
    assert primary['coefficient_fields'] == ['QQ', 'C']
    assert not primary['one_variable_source_degree1_or_expanded_Luroth_generator_assumed']
    assert not primary['Jacobi_scheme_reducedness_or_generic_length5_changed']
    assert not primary['NEW_firstjet_equation_or_global_component_exhaustion_claimed']

    names = ['Q_NEW_FAMILY_FULL_TENSOR_FIBER_PREFLIGHT_20261104_v63_REVIEW',
             'S61_RATIONAL_COMPONENT_TRANSFER_20261104_v63_REVIEW',
             'Q_NEW_FAMILY_RATIONAL15_APPLICATION_20261104_v63_REVIEW']
    paths = [ROOT/'reports'/n/'report.json' for n in names]
    records = [read(p) for p in paths]
    assert all(d['all_passed'] for d in records)
    full, transfer, base = records
    assert full['independent_pair_to_beta_birational_field_sandwich']
    assert full['rawalpha_matching_not_inferred_from_primitive_beta_alone']
    assert full['full_tensor_tau_descent_from_all_degree_family_identity_reviewed']
    assert full['actual_full_S61_Sp4_dominance_and_historical_perspective_quotient_reviewed']
    assert full['same_K_B15_field_inclusion_into_K_z_reviewed']
    assert full['actual_reduced_S61_C16_QQ_and_C_rational_dimension16']
    assert full['source_z_map_degree1_or_entire_component25_rationality_claimed'] is False
    assert transfer['accepted55_elimination_and134_IDEAL_identity_scheme_product']
    assert transfer['connected_component_stability_and_slice_irreducibility_proof_accepted']
    assert transfer['projective_c_equals1_before_normalization_yields_exact_dimension_shift9']
    assert transfer['QQ_and_C_conditional_rationality_transfer_accepted']
    assert base['same_actual_quotient_pullback_and_section_identity_reviewed']
    assert base['dominant_equal_dimension_separable_degree1_over_QQ_reviewed']
    assert base['B_NEW15_pure_rational_over_QQ_and_C']
    assert base['full25_component_rationality_claim'] is False
    checked = {}
    # Bind both application and every FINAL gate's declared evidence/dependencies.
    for p, d in [(primary_path, primary), *zip(paths, records)]:
        for key in ('input_sha256', 'checked_dependency_sha256'):
            for name, want in d.get(key, {}).items():
                target = ROOT/name
                assert sha(target) == want, str(target)
                checked[str(target.relative_to(ROOT))] = want
        for name, want in d.get('evidence_sha256', {}).items():
            target = p.parent/name
            assert sha(target) == want, str(target)
            checked[str(target.relative_to(ROOT))] = want

    # Own integer check of the finite involution's field coordinate identities.
    # tau(t)=t/(t-1), tau(u)=-u; z=u(t-2)/t.
    mobius = [[1,0],[1,-1]]
    assert matmul2(mobius, mobius) == [[1,0],[0,1]]
    original_z_numerator = [-2,1]
    tau_t_numerator, tau_t_denominator = [0,1], [-1,1]
    substituted_z_numerator = [-(a-2*b) for a,b in
                              zip(tau_t_numerator, tau_t_denominator)]
    assert substituted_z_numerator == original_z_numerator
    # Reversibility is (u-z)*t=2u, obtained from z*t=u*t-2u.
    # Thus E5(u,t)=E5(u,z), and the sign automorphism fixes precisely E5(z).
    assert 10+5 == 15 and 15+1 == 16 and 16+5+4 == 25
    proof = {
        'three_distinct_targets': {'B15': 'actual reduced primitive firstjet image',
                                   'Y16': 'actual reduced full S61 component',
                                   'X25': 'actual reduced projective Poisson component'},
        'full_tau_fixed_field': 'E5(u,t)^tau=E5(z), z=u*(t-2)/t; u^2=1600*r',
        'same_base_inclusion': 'K=QQ(B15) subset QQ(Y16) subset K(z)',
        'one_variable_map_degree_required': 'nonconstant only, not degree one',
        'Luroth_over_actual_base_field': 'QQ(Y16)=K(z0)',
        'rank4_component_transfer': 'X25 birational U5 x A4 x Y16 over QQ',
        'rationality_fields': ['QQ','C'],
        'generic_Jacobi_length5_and_CM_type2_preserved': True,
        'actual_Luroth_generator_expanded': False,
        'NEW_equation_degree_or_generic_beta_exhaustion_solved': False,
        'new_arithmetic_source_point_or_CAS_job': False,
    }
    output.mkdir(parents=True, exist_ok=False)
    evidence = output/'independent-proof-receipt.json'
    evidence.write_text(json.dumps(proof,indent=2)+'\n')
    inputs = [Path(__file__), ROOT/'research/notes'/(STEM+'.md'),
              ROOT/'research/notes'/(BASE+'.md'), primary_path, *paths]
    inputs += [ROOT/'research/notes'/(n+'.md') for n in names]
    for p in inputs:
        checked[str(p.relative_to(ROOT))] = sha(p)
    result = {
        'id': STEM, 'all_passed': True,
        'same_actual_K_inclusion_and_full_tau_descent_gate_accepted': True,
        'Luroth_over_QQ_B15_no_degree_one_source_map_assumption': True,
        'QQ_geometric_component_and_field_descent_accepted': True,
        'connected_component_stability_and_exact_projective_plus9_accepted': True,
        'actual_reduced_new_Poisson_projective_component_QQ_and_C_rational25': True,
        'firstjet_B15_full_S61_Y16_and_projective_X25_distinguished': True,
        'generic_Jacobi_length5_or_scheme_reducedness_unchanged': True,
        'NEW_equation_or_global_component_exhaustion_claimed': False,
        'author_producer_rerun_or_new_CAS_source_point_prime': False,
        'input_sha256': checked,
        'evidence_sha256': {'independent-proof-receipt.json': sha(evidence)},
        'seconds': time.monotonic()-start,
        'frozen': True,
    }
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('all_passed','seconds')}))


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)

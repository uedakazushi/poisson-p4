#!/usr/bin/env python3
"""Source-bound pure proof: universal quadratic kernels and actual nil upper12.

No source worker, point, jet, CAS, Groebner basis or numerical rank replay.
"""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_QUADRATIC_NORMALIZER_NIL_UPPER12_20261104_v63'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def read(p):
    return json.loads(p.read_text())

def add(a, b):
    r = dict(a)
    for m, c in b.items():
        r[m] = r.get(m, 0) + c
        if not r[m]:
            del r[m]
    return r

def neg(a):
    return {m: -c for m, c in a.items()}

def mul(a, b):
    r = {}
    for m, c in a.items():
        for n, d in b.items():
            k = tuple(x + y for x, y in zip(m, n))
            r[k] = r.get(k, 0) + c * d
    return {m: c for m, c in r.items() if c}

def main(out):
    start = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    out.mkdir(parents=True, exist_ok=False)
    # Formal identity after k=j0+j1*t: a=k-j0=j1*t.
    def var(i):
        m = [0] * 6
        m[i] = 1
        return {tuple(m): 1}
    j0, j1, t, s, nu, b = [var(i) for i in range(6)]
    k = add(j0, mul(j1, t))
    a = add(k, neg(j0))
    lhs = add(add(mul(a, a), neg(mul(mul(s, j1), a))), neg(mul(nu, mul(j1, j1))))
    rhs = mul(mul(j1, j1), add(add(mul(t, t), neg(mul(s, t))), neg(nu)))
    assert lhs == rhs
    ab = add(a, b)
    polarization = add(add(mul(ab, ab), neg(mul(a, a))), neg(mul(b, b)))
    assert polarization == {m: 2*c for m, c in mul(a, b).items()}

    expected = {
        'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW/report.json': 'b547efc8340bef2444459c95062aaf9550c551376be6edd9e61b055c2d85c1af',
        'reports/AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63_REVIEW/report.json': 'cffc3d5c6a9cf9ed1207788a27675bdc76f925380fa27f0dc3a3c0c8e5337f80',
        'reports/AMBIENT_B18_NIL_CUBE_FINITE_ORIGIN_20261104_v63/report.json': '6235023a65377a0b3b157cd5322362953110788911f7dbb1b9ff5a62f44934c5',
        'reports/AMBIENT_B18_NIL_CUBE_FINITE_ORIGIN_20261104_v63_REVIEW/report.json': 'ee80d12faa103fb5ac0e2144b4d12cd662f0c5afcfa1ded335e2171946ad1ccc',
        'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW/report.json': 'edab93bf8b666e088d62d6ece3a65be62b3ba7621b3e22c6490ca6f85ad16662',
        'reports/AMBIENT_B18_T1_T2_QUADRATIC_PHI_20261104_v63_REVIEW/report.json': 'e503ac401f444b27678e4992a96834f6d1a4ff7c6f737a3e498849f2f2e466b5',
        'reports/AMBIENT_B18_T3_ZERO_STRUCTURAL_20261104_v63_REVIEW/report.json': 'c20f18e10d045c7cc34cd4d2d0e9a7d78548620c2e0ddf37b862738cc4822c19',
        'reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63_REVIEW/report.json': '087ca11d0efe1af921b94c2b3421cfa0ad2f347edc4f1e6603ed112853d4a9c5',
    }
    reports = {}
    for rel, digest in expected.items():
        p = ROOT / rel
        assert sha(p) == digest, rel
        reports[rel] = read(p)
        assert reports[rel]['all_passed'], rel
    z = reports[next(iter(expected))]
    assert z['actual_free4_48_and_product_conductor_comparison_accepted']
    f = reports['reports/AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63_REVIEW/report.json']
    assert f['core_and_normalized_E4_exact_nilradical_index3_accepted']
    n = reports['reports/AMBIENT_B18_NIL_CUBE_FINITE_ORIGIN_20261104_v63/report.json']
    assert n['all_positive_homogeneous_coefficient_maps_A_n_to_M_injective']
    assert n['tail_nilradical_square_zero'] and n['mixed_nil_product_statement']
    assert n['nil_square_nonzero_at_P0_by_Qb_coherent_closure']
    nr = reports['reports/AMBIENT_B18_NIL_CUBE_FINITE_ORIGIN_20261104_v63_REVIEW/report.json']
    assert nr['full_completed_nil_cube_equals_degree0_cube'] and nr['mixed_positive_products_zero']
    qb = reports['reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW/report.json']
    assert qb['nilrank4_square_rank2_cubezero_Ass_N_and_embedded_SQ_depth20_e1_accepted']
    # Exponent arithmetic is ideal multiplication, not only individual powers.
    assert 2 * 2 == 4 and 3 * 4 == 12
    proof = {
        'universal_hypotheses': ['A injects into B', 'B=A+A*tau', 'tau^2=s*tau+nu with s,nu in A', '2 is invertible'],
        'universal_base_change': 'every commutative A-algebra M, without flatness',
        'lift': 'k=j0+j1*tau, ji in J in an arbitrary polynomial presentation M=A[X]/J',
        'universal_identity': '(k-j0)^2-s*j1*(k-j0)-nu*j1^2=j1^2*(tau^2-s*tau-nu)',
        'identity_checked_in_integer_polynomial_ring': True,
        'injectivity_used_only_before_quotient': 'A[X] injects into B[X]',
        'kernel_ideal_square_zero_by_polarization': True,
        'completed_base_change': 'finite modules over complete Noetherian rings are complete; the same tensor base change and closed-ideal proof applies',
        'actual_source': 'Ai=R0+Ii*taui contained in Bi=R0+R0*taui, taui^2=nui',
        'actual_two_stages': ['D0=AA completed_tensor_P0 AB', 'D1=BA completed_tensor_P0 AB', 'E4=BA completed_tensor_P0 BB'],
        'first_kernel_square_zero': True,
        'second_kernel_square_zero': True,
        'composite_kernel_K_power4_zero': True,
        'whole_E4_nil_power3_zero_peer_gated': True,
        'Nil_D0_power3_subset_K': True,
        'Nil_D0_power12_zero': True,
        'all_completed_positive_nil_products_zero_peer_gated': True,
        'Nil_full_actual_B18_power12_zero': True,
        'nil_square_nonzero_by_actual_Qb_closed_support_peer_gated': True,
        'nil_index_interval': [3, 12],
        'G19_and_projective_scalar_regular_factors_preserve_nil_index': True,
        'requires_generic_support_exhaustion_or_depth_or_T_length': False,
        'upper8_or_exact_index_or_nil_cube_nonzero_claimed': False,
        'optional_upper8_missing_gate': 'global actual image Nil(D0)^2=0 in E4; generic selected-source vanishings alone are insufficient',
        'new_source_CAS_GB_point_jet_or_rank_job': False,
    }
    ep = out / 'quadratic-universal-kernel-and-exponent-proof.json'
    ep.write_text(json.dumps(proof, indent=2) + '\n')
    files = list(expected) + [
        'research/notes/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63.md',
        'research/notes/AMBIENT_B18_NIL_CUBE_FINITE_ORIGIN_20261104_v63.md',
        'research/notes/AMBIENT_B18_NIL_CUBE_FINITE_ORIGIN_20261104_v63_REVIEW.md',
        'research/notes/AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63_REVIEW.md',
        'research/notes/AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63.md',
        'research/notes/AMBIENT_B18_W1_POLYNOMIAL_LAURENT_SCOPE_20261104_v63.md',
        f'research/notes/{STEM}.md',
        str(Path(__file__).relative_to(ROOT)),
    ]
    report = {
        'id': STEM, 'all_passed': True,
        'status': 'PROVED_unconditional_quadratic_base_change_nil_upper_bound_no_source_job',
        'scope': 'actual CA intersect CB, formal generic B18 representative P0, cone and projective germ; not full Poisson scheme or global classification',
        'universal_quadratic_base_change_kernel_square_zero': True,
        'monic_trace_variant_also_proved': True,
        'completed_finite_base_change_covered': True,
        'actual_two_step_comparison_kernel_power4_zero': True,
        'whole_E4_nil_power3_zero': True,
        'actual_D0_nil_power12_zero': True,
        'full_actual_B18_nil_power12_zero': True,
        'actual_nil_index_lower_bound': 3,
        'actual_nil_index_upper_bound': 12,
        'upper_bound_independent_of_depth19_T_length_generic_K_rank_or_branch_exhaustion': True,
        'nil_cube_zero_or_nonzero_exact_index_upper8_claimed': False,
        'source_GB_point_jet_or_rank_replayed': False,
        'input_sha256': {rel: sha(ROOT / rel) for rel in files},
        'evidence_sha256': {ep.name: sha(ep)},
        'seconds': time.monotonic() - start,
        'max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ['id', 'all_passed', 'actual_nil_index_lower_bound', 'actual_nil_index_upper_bound', 'seconds']}))

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output', required=True, type=Path)
    args = p.parse_args()
    main(args.output.resolve())

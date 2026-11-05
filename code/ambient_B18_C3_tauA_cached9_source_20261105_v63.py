#!/usr/bin/env python3
"""Source-only nine tauA cubic adjoints; never pairs cubic cache leaves."""
import hashlib
import json
import math
import signal
import time
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_C3_TAUA_CACHED9_SOURCE_20261105_v63'
PATHS = {
    'C2': 'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
    'raw33': 'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
    'ell': 'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63/all164-cached-pure-ell-readbacks-and-physical12-second-values.json',
    'pair': 'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
    'frame': 'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63/actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json',
    'frame_peer': 'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63_REVIEW/report.json',
    'C2_peer': 'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
    'raw_peer': 'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63_REVIEW/report.json',
    'ell_peer': 'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63_REVIEW/report.json',
}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def mul(a, b):
    return [[sum((x * y for x, y in zip(row, col)), Q(0))
             for col in zip(*b)] for row in a]


def enc(a):
    return [[str(v) for v in r] for r in a]


def main():
    started = time.monotonic()
    signal.alarm(60)
    out = ROOT / 'reports' / STEM
    out.mkdir(exist_ok=False)
    pins = {p: sha(ROOT / p) for p in PATHS.values()}
    d = {k: json.loads((ROOT / p).read_text()) for k, p in PATHS.items()}
    assert all(v['all_passed'] for k, v in d.items() if k.endswith('_peer'))
    basis = d['raw33']['actual33_basis']
    i3 = basis.index([0, 3, 0])
    i21 = basis.index([0, 2, 1])
    assert [i3, i21] == [29, 30]
    # Every old ambient normal form preserves these literal pure-plane Taylor
    # functionals, not merely total degree.  No cubic scalar is paired here.
    readback_receipts = []
    for m, nf in zip(d['raw33']['basis74'], d['raw33']['all74_normal_forms']):
        a = dict(nf)
        got = [Q(a.get(i3, '0')), Q(a.get(i21, '0'))]
        want = [Q(m == [0, 3, 0]), Q(m == [0, 2, 1])]
        assert got == want
        readback_receipts.append([m, [str(v) for v in got]])
    # The beta2 and all mixed correction terms are retained in the proof;
    # only the eA projection has no such carry.
    frame = d['frame']['all12_verified_source_records']
    rows = [r['physical_row'] for r in frame]
    b = [[Q(r['bB_second_plane'][j]) for r in frame] for j in range(3)]
    piv = [0, 4, 5]
    p = [[b[i][j] for j in piv] for i in range(3)]
    det = (p[0][0] * (p[1][1] * p[2][2] - p[1][2] * p[2][1])
           - p[0][1] * (p[1][0] * p[2][2] - p[1][2] * p[2][0])
           + p[0][2] * (p[1][0] * p[2][1] - p[1][1] * p[2][0]))
    assert det == Q(5, 768)
    # The cyclic two-by-two cross-product formula gives adj(P) directly.
    adj = [[(p[(j + 1) % 3][(i + 1) % 3] * p[(j + 2) % 3][(i + 2) % 3]
             - p[(j + 1) % 3][(i + 2) % 3] * p[(j + 2) % 3][(i + 1) % 3])
            for j in range(3)] for i in range(3)]
    inv = [[v / det for v in row] for row in adj]
    eye = [[Q(i == j) for j in range(3)] for i in range(3)]
    assert mul(p, inv) == eye == mul(inv, p)
    free = [i for i in range(12) if i not in piv]
    k = [[Q(0) for _ in free] for _ in range(12)]
    for c, f in enumerate(free):
        k[f][c] = Q(1)
        for i, pi in enumerate(piv):
            k[pi][c] = -sum((inv[i][j] * b[j][f] for j in range(3)), Q(0))
    assert mul(b, k) == [[Q(0)] * 9 for _ in range(3)]
    assert [[k[i][j] for j in range(9)] for i in free] == [
        [Q(i == j) for j in range(9)] for i in range(9)]
    # Compose the already reviewed cubic gA(r,s) with s=r+ell/5.
    ga = [Q(0)] * 4
    for mon, coefficient in d['pair']['gA_plane3']:
        nr, ns = mon
        assert nr + ns == 3
        for ell_degree in range(ns + 1):
            ga[ell_degree] += Q(coefficient) * math.comb(ns, ell_degree) / (5 ** ell_degree)
    assert ga == [Q(125, 24), Q(5, 3), Q(1, 6), Q(2, 375)]
    ratio = ga[1] / ga[0]
    assert ratio == Q(8, 25)
    folded = []
    for j in range(9):
        terms = []
        for i in range(12):
            if k[i][j]:
                terms += [[i, 'delta', i21, str(k[i][j])],
                          [i, 'delta', i3, str(-ratio * k[i][j])]]
        folded.append(terms)
    coefficients = [Q(t[3]) for terms in folded for t in terms]
    max_weight_height = max(max(abs(q.numerator).bit_length(), q.denominator.bit_length())
                            for q in coefficients)
    # SOURCE-only height bound, based on already literal cached leaf strings.
    # The leaves are inspected for numerator/denominator sizes, never summed.
    leaves = []
    for i, r in enumerate(d['C2']['physical12_coefficients']):
        assert r['physical_row'] == rows[i]
        values = dict(r['delta'])
        for idx in [i3, i21]:
            q = Q(values.get(idx, '0'))
            leaves.append([i, idx, max(abs(q.numerator).bit_length(), q.denominator.bit_length())])
    leaf_height = max(x[2] for x in leaves)
    max_terms = max(map(len, folded))
    # Each product height <= weight+leaf; summing <=8 terms using a common
    # denominator costs at most n*productheight+ceil(log2 n)+1.
    product_bits = max_weight_height + leaf_height
    dot_bits = max_terms * product_bits + math.ceil(math.log2(max_terms)) + 1
    assert dot_bits <= 2048
    data = {
        'comparison_scope': 'Original C=coker(D0→E4) over the SAME physical R0; C3=C/ell^3 C. These are modules, not rings.',
        'first_frame_source': 'Accepted RA_i=zi eA, RB_i=(zi+beta_i)eB, RH_i=F_i tauB. All first delta=0; RH2_plane=0.',
        'physical_row_order': rows,
        'bB_plane2_rows3_columns12': enc(b),
        'fixed_pivot_columns': piv,
        'fixed_minor': str(det),
        'fixed_inverse': enc(inv),
        'both_inverse_products_identity': True,
        'free_columns9': free,
        'complete_constant_kernel12_by9': enc(k),
        'all27_kernel_equations_and_free_I9': True,
        'raw33_literal_pureplane_cubic_dual_receipts': readback_receipts,
        'pureplane_cubic_cache_indices': {'r^3': i3, 'r^2*ell': i21},
        'gA_cubic_r_ell_order': [str(q) for q in ga],
        'gA_reduction_ratio': str(ratio),
        'tauA_first_cancellation': 'q_l=F_l tauA+sum_i (bB_l)_zi RH_i. Its first is zero, its eA component is exactly delta_l, and its plane2 component is -bB_l,2 eH.',
        'nine_source_cycles_through_order2': 'For each k in ker(bB2), q_k has no plane2 component. Clear every zz and z*r/z*ell term by the ACTUAL first21 frame. This gives an actual order>=3 relation with no arbitrary unavailable coefficients set zero.',
        'full_plane_cubic_formula': 'eA=sum_l k_l delta_l,3. eB=sum_l,i k_l h_li delta_i,3 - sum_i (q_k,2,eB,zi*r*r+q_k,2,eB,zi*ell*ell) beta_i,2. eH is retained but dies modulo the known plane2*eH initial ideal. All zz corrections retain a z multiplier and have no plane3 contribution.',
        'eA_carry_zero_source_reason': 'Only RA corrects eA components, and RA_i=zi eA exactly. RB/RH and every beta2 or delta_i3 correction live in eB/eH; none can contribute to eA.',
        'nine_cached_eA_cubic_adjoints': folded,
        'cubic_values_not_evaluated': True,
        'what_the_adjoints_read': 'Coefficient of r^2*ell in eA after reduction modulo (gA3,ell^3): delta_(2,1)-(8/25)delta_(3,0), then the constant kernel12x9. Unknown delta_(1,2), canonical gauge3 and all other physical thirds have coefficient zero in this functional, not assigned zero.',
        'nonzero_value_consequence': 'One nonzero dot yields an actual cubic initial relation whose eA component has nonzero r^2*ell. Multiplication by ell gives a nonzero eA quartic r^2*ell^2 in the known CI(3,3) model. Hence two independent graded relations lower known length21 to <=19. This would force at least one of the two remaining C torsion exponents to equal2, but would not prove both2.',
        'zero_value_consequence': 'All nine dots zero is only a negative result for these nine projections; other cubic directions and quartic corrections remain open.',
        'height_profile': {'max_weight_bits': max_weight_height, 'max_cached_leaf_bits': leaf_height,
                           'max_dot_terms': max_terms, 'product_bits_bound': product_bits,
                           'source_dot_bits_bound': dot_bits, 'hard_bits': 2048},
        'cache_leaf_size_receipts_only_no_dots': leaves,
        'worker_profile': {'static_owned_numeric_envelope': 18000, 'hard_owned_numeric': 30000,
                           'hard_seconds': 60, 'hard_memory_bytes': 2147483648,
                           'fresh_RAW_input_plus_output_bytes_bound': 400000,
                           'hard_RAW_bytes': 33554432,
                           'new_physical_source_gauge_or_jet_coefficients': 0,
                           'cache_inputs': 24, 'kernel_sources': 9,
                           'max_folded_products': sum(map(len, folded))},
        'no_new_source_value_rank_C3_exponent_or_nil_claim': True,
    }
    evidence = out / 'nine-complete-tauA-source-kernel-and-unpaired-cubic-cache-adjoints.json'
    evidence.write_text(json.dumps(data, indent=2) + '\n')
    report = {
        'id': STEM, 'all_passed': True, 'status': 'SOURCE_ONLY_nine_cached_tauA_cubic_adjoints_NOT_PAIRED',
        'cubic_values_evaluated': False, 'input_sha256': pins,
        'code_sha256': sha(Path(__file__)),
        'note_sha256': sha(ROOT / 'research/notes' / (STEM + '.md')),
        'evidence_sha256': {evidence.name: sha(evidence)},
        'seconds': time.monotonic() - started,
    }
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))
    signal.alarm(0)


if __name__ == '__main__':
    main()

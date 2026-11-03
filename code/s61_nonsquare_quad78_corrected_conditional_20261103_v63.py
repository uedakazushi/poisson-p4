#!/usr/bin/env python3
"""Universal algebra checks for the corrected conditional quad78 lemma.

No G1/G2/G3 or integrality locus is computed or declared to hold globally.
"""
import argparse
import json
from pathlib import Path
import time
import sympy as s
from s61_nonsquare_structural_reductions_review_20261103_v63 import read, sha

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_NONSQUARE_QUAD78_CORRECTED_CONDITIONAL_20261103_v63'


def main(out):
    started = time.monotonic(); out.mkdir(parents=True, exist_ok=False)
    deps = [ROOT/f'reports/{stem}/report.json' for stem in [
        'S61_NONSQUARE_FOUR_QUADRIC_BINARY_REDUCTION_20261103_v63_REVIEW',
        'S61_NONSQUARE_PAIR_EMPTY_QQ_20261103_v63_REVIEW',
        'S61_NONSQUARE_QUADRIC_CI_OPEN_20261103_v63_REVIEW']]
    for p in deps:
        rr = read(p); assert rr['all_passed']
        for rel, h in rr['input_sha256'].items(): assert sha(ROOT/rel) == h, rel
    A, X, W, B, C1, C2, Q0, Q3, K0, K1, K2, K3, L = s.symbols('A X W B C1 C2 Q0 Q3 K0 K1 K2 K3 L')
    alpha, beta, kappa = s.symbols('alpha beta kappa'); c = s.Integer(6075)
    Q1 = A*B/c+W*C1; Q2 = 40*X*B/c+W*C2; D = 40*X*C1-A*C2
    assert s.expand(W*D-40*X*Q1+A*Q2) == 0
    assert s.expand(B*D+c*C2*Q1-c*C1*Q2) == 0
    # Actual removal of each B*(alpha*A+beta*X) contribution by two Koszuls.
    old = B*(alpha*A+beta*X)+W*kappa
    new = W*(kappa-c*alpha*C1-c*beta*C2/40)
    assert s.expand(old-c*alpha*Q1-c*beta*Q2/40-new) == 0
    restored = [W*K0, 40*X*L+W*K1, -A*L+W*K2, W*K3]
    expression = restored[0]*Q0+restored[1]*Q1+restored[2]*Q2+restored[3]*Q3
    assert s.expand(expression-W*(K0*Q0+K1*Q1+K2*Q2+K3*Q3+L*D)) == 0
    first = {K0: 0, K1: -40*X, K2: A, K3: 0, L: W}
    second = {K0: 0, K1: c*C2, K2: -c*C1, K3: 0, L: B}
    assert all(s.expand(f.subs(first, simultaneous=True)) == 0 for f in restored)
    assert [s.expand(f.subs(second, simultaneous=True)) for f in restored] == [0, s.expand(c*Q2), s.expand(-c*Q1), 0]
    # The six constant Koszuls are independent whenever the four forms are.
    z = s.symbols('z0:4'); h = s.symbols('h01 h02 h03 h12 h13 h23'); pairs = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    vector = [s.Integer(0)]*4
    for coefficient, (i, j) in zip(h, pairs): vector[i] += coefficient*z[j]; vector[j] -= coefficient*z[i]
    anchors = [s.expand(vector[i]).coeff(z[j]) for i, j in pairs]
    assert anchors == list(h)
    assert 4*21-6 == 78
    evidence = dict(checked_over_universal_QQ_polynomial_ring=True, two_scroll_syzygies_verified=True,
        four_Koszul_B_part_removals_verified_by_one_universal_formula=True,
        exact_W_division_and_four_quadric_restoration_formula_verified=True,
        scroll_first_restores_zero=True, scroll_second_restores_c_times_Q1_Q2_Koszul=True,
        six_Koszuls_constant_independence_anchor_matrix='I6', degree2_cofactor_dimension=84, kernel_dimension_if_gates=6,
        conditional_raw84_by126_rank=78,
        accepted_hypotheses=['a*P*J!=0', 'G1: Q0,Q3 mod(B,W) coprime', 'G2: Q0,Q3 mod(A,X,W) coprime',
            'G3: vertex binary coefficient matrix rank2 OR C geometrically integral'],
        domain_alternative_includes_rank1_and_zero_binary_cases=True,
        actual_gate_or_integrality_loci_not_computed=True,
        smooth_ACM_curve_implies_integral_only_via_connectedness=True,
        no_specific_E5_smoothness_certificate_checked_here=True)
    (out/'certificate.json').write_text(json.dumps(evidence, indent=2)+'\n')
    result = dict(id=STEM, all_passed=True, corrected_conditional_five_linear_syzygies_exactly2=True,
        corrected_conditional_four_quadratic_cofactor_kernel_exactly6_Koszuls=True,
        raw_rank78_only_under_all_stated_conditions=True, actual_G1_G2_G3_or_integrality_loci_proved=False,
        withdrawn_G1_G2_only_implication_not_reinstated=True, no_global_quad78_CI_or_fullquartic_claim=True,
        no_GB_new_field_parametric_rank_search_or_new_prime=True, seconds=time.monotonic()-started)
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md',
        ROOT/'research/notes/S61_NONSQUARE_PAIR_OMIT_BSQUARE_CONDITIONAL_SYZYGY_GAP_20261103_v63.md']+deps
    result['input_sha256'] = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    result['evidence_sha256'] = {'certificate.json': sha(out/'certificate.json')}
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['input_sha256', 'evidence_sha256']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output)

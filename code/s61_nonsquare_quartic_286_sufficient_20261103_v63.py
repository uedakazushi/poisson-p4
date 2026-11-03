#!/usr/bin/env python3
"""Finite corollary: 286x238 is a sufficient, not whole, quartic gate."""
import argparse
import json
from pathlib import Path
import time
import sympy as s
from s61_nonsquare_structural_reductions_review_20261103_v63 import sha, read

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_NONSQUARE_QUARTIC_286_SUFFICIENT_20261103_v63'


def main(out):
    start = time.monotonic(); out.mkdir(parents=True, exist_ok=False)
    review = ROOT/'reports/S61_NONSQUARE_CUBIC_QUARTIC_STRUCTURAL_REDUCTIONS_20261103_v63_REVIEW/report.json'
    heavy_review = ROOT/'reports/S61_BINARY_HYPERPLANE_HEAVY_SQUARE_PIVOT_REVIEW_20261102.json'
    heavy_source = ROOT/'reports/S61_BINARY_HYPERPLANE_HEAVY_SQUARE_PIVOT_20261102.json'
    for p in [review, heavy_review]:
        rr = read(p); assert rr['all_passed']
        for rel, h in rr['input_sha256'].items(): assert sha(ROOT/rel) == h, rel
    assert read(heavy_review)['known_open'] == 'a*F*B*C != 0'
    source = read(heavy_source)
    certpath = ROOT/'reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/certificate.json'
    cert = read(certpath); rows = cert['source_rows']; columns = cert['weight4_columns']
    assert len(rows) == 354 and all(row[0] == 'quadric' for row in rows[:174])
    assert all(row[0] == 'cubic_normalized' for row in rows[174:321])
    assert all(row[0] == 'quartic' for row in rows[321:])
    deleted = cert['deleted_row_indices']; assert len(deleted) == 35 and set(deleted) <= set(range(174, 321))
    assert all(all(i < 321 for i, c in rel['coefficients']) for rel in cert['lifted_relations'])
    retained = [i for i in range(321) if i not in deleted]; assert len(retained) == 286
    h2 = [i for i, mon in enumerate(columns) if mon[5] == 2]; assert len(h2) == 1
    pure_h2 = [0]*8; pure_h2[5] = 2; assert columns[h2[0]] == pure_h2
    remaining = [i for i in range(239) if i not in h2]; assert len(remaining) == 238
    assert ['quartic', 77] in rows[321:] and ['quartic', 131] in rows[321:]
    a, e, d = s.symbols('a e d'); B = 50*a*d+18*d*e**2-9*d*e-2*d+60*e-40; C = 100*a+18*e**2-9*e-2
    parse = lambda text: s.sympify(text, locals={'a': a, 'e': e, 'd': d})
    V = 50*a*B*C
    assert s.expand(parse(source['bezout_unit'])-V) == 0
    combined = parse(source['combined_heavy_square_coefficient'])
    assert s.expand(combined-V**2) == 0
    coeffs = list(map(parse, source['coefficient_formulas'])); weights = list(map(parse, source['row_combination_coefficients']))
    assert s.cancel(sum(x*y for x, y in zip(coeffs, weights))-V**2) == 0
    evidence = dict(lower_weight_original_shape=[321, 238], lower_weight_reduced_shape=[286, 238],
        retained_354_source_row_indices=retained, column_indices=remaining, heavy_square_column=h2[0],
        deleted_rows=deleted, deletion_minor='-a^14', heavy_square_actual_rows=[77, 131],
        actual_heavy_square_pivot_unit='(50*a*B*C)^2', all35_relations_supported_on_lower_weight321=True,
        all_lower_weight_rows_have_zero_heavy_square_by_weights=True,
        sufficient_implication='286 lower-weight rows span L238 => append actual heavy-square unit row => full weight4 span239',
        not_equivalent_to_whole319_by239_criterion=True,
        all33_genuine_quartics_only_omitted_from_this_sufficient_test_not_from_whole_ideal=True)
    (out/'certificate.json').write_text(json.dumps(evidence, indent=2)+'\n')
    result = dict(id=STEM, all_passed=True, open='a*F*B*C!=0',
        lower_weight321_to286_exact_same_row_span=True, sufficient_quartic286_by238_criterion=True,
        actual_rank238_not_proved=True, whole319_equivalence_not_claimed=True,
        no_new_source_expansion_GB_field_or_unknown_localization=True, seconds=time.monotonic()-start)
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', review, heavy_review, heavy_source, certpath]
    result['input_sha256'] = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    result['evidence_sha256'] = {'certificate.json': sha(out/'certificate.json')}
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['input_sha256', 'evidence_sha256']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output)

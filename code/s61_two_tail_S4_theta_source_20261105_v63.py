#!/usr/bin/env python3
"""Compare only already-saved exact S4/Theta/p coefficient tables."""
import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main(output):
    output.mkdir(parents=True, exist_ok=False)
    pins = {}

    def load(path):
        file = ROOT / path
        pins[path] = sha(file)
        return json.loads(file.read_text())

    factors = load('reports/S61_TWO_TAIL_H_EXACT_QQ_ENTRY_20261105_v63/canonical-factor-candidates.json')
    peer = load('reports/S61_TWO_TAIL_H_EXACT_QQ_ROOT_DENSE_REVIEW_20261105_v63/report.json')
    assert peer['all_passed'] and peer['S4_geometrically_irreducible_proved']
    old = load('reports/S61_NONSQUARE_RAW_H_P_BRANCH_PROFILES_20261104_v63/certificate.json')
    old_report = load('reports/S61_NONSQUARE_RAW_H_P_BRANCH_PROFILES_20261104_v63/report.json')
    assert pins['reports/S61_NONSQUARE_RAW_H_P_BRANCH_PROFILES_20261104_v63/certificate.json'] == old_report['certificate_sha256']
    coverage = load('reports/S61_NONSQUARE_PREPOWER_FIFTH_PURE_PIVOT_20261104_v63/coverage.json')
    load('reports/S61_NONSQUARE_PREPOWER_FIFTH_PURE_PIVOT_20261104_v63_REVIEW/report.json')
    s4 = next(record for record in factors['factors'] if record['total_degree'] == 4)
    assert s4['multiplicity'] == 1 and len(s4['terms']) == 11

    def decode(records):
        return {tuple(m): Fraction(c) for m, c in records if Fraction(c)}

    theta = decode(old['Theta'])
    factor = decode(s4['terms'])
    assert theta == {m: 20 * c for m, c in factor.items()}
    g1 = decode(old['g1'])
    product = {}
    for m, c in g1.items():
        for n, d in factor.items():
            k = tuple(x + y for x, y in zip(m, n))
            product[k] = product.get(k, Fraction()) + c * d
    product = {m: c for m, c in product.items() if c}
    p = decode(old['source_p'])
    assert product == {m: 37500 * c for m, c in p.items()}
    coefficient_receipt = [dict(exponent=list(m), S4=str(factor[m]), Theta=str(theta[m]), ratio=20)
                           for m in sorted(theta)]
    certificate = dict(
        source_only=True, variables=['a', 'e', 'd'], exact_coefficient_receipt=coefficient_receipt,
        all11_Theta_equals20S4=True, old_source_p_equals_g1S4_over37500=True,
        old_factor_product_pairs=33, source_p_terms=len(p),
        accepted_H_factorization_rewritten='H=R7*Theta/20',
        geometric_component_identified='V(S4)=V(Theta)',
        S4_absolute_irreducibility_from_accepted_peer=True,
        scope='Exact saved polynomial tables; no source evaluation or new rank/minor/factor computation.',
        accepted_Theta_heavy_scope=coverage,
        all_Theta_full57_NOT_claimed=True, actual_bad_Fitting_equality_NOT_claimed=True,
        selected_tail_criterion_vanishes_on_Theta=True,
        remaining_other_selected_component='V(R7), absolute irreducibility pending R7 ABS14 value and peer')
    evidence = output / 'coefficient-identity.json'
    evidence.write_text(json.dumps(certificate, indent=2) + '\n')
    for path in [str(Path(__file__).relative_to(ROOT)),
                 'research/notes/S61_TWO_TAIL_S4_THETA_SOURCE_20261105_v63.md']:
        pins[path] = sha(ROOT / path)
    report = dict(id=output.name, all_passed=True, frozen=True, source_only=True,
                  exact_S4_Theta_recognition=True, no_new_physical_query=True,
                  no_new_factor_GCD_rank_minor_or_point=True,
                  no_global_full57_or_exhaustion_claim=True,
                  input_sha256=pins, evidence_sha256={evidence.name: sha(evidence)})
    file = output / 'report.json'
    file.write_text(json.dumps(report, indent=2) + '\n')
    raw = sum(p.stat().st_size for p in output.iterdir())
    assert raw < 1024 * 1024
    print(json.dumps(dict(report_sha256=sha(file), raw_new_bytes=raw,
                          identity='Theta=20*S4; 37500*p=g1*S4'), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output.resolve())

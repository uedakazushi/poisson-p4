#!/usr/bin/env python3
"""ONE gated existing-cache nine-dot job; no physical or implicit evaluation."""
import argparse
import hashlib
import json
import re
import resource
import signal
import time
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_C3_TAUA_CACHED9_VALUES_20261105_v63'
SOURCE = 'reports/AMBIENT_B18_C3_TAUA_CACHED9_SOURCE_20261105_v63/report.json'
SOURCE_SHA = '8234e80934cba6d761eb72666d3d4337fcfa07a1c5d12ca9b0776f85d9339bd0'
SOURCE_EVIDENCE = 'reports/AMBIENT_B18_C3_TAUA_CACHED9_SOURCE_20261105_v63/nine-complete-tauA-source-kernel-and-unpaired-cubic-cache-adjoints.json'
C2 = 'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json'
LIMIT = 2048


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def numeric_count(x):
    if isinstance(x, (int, float)):
        return 1
    if isinstance(x, str):
        return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?', x)))
    if isinstance(x, list):
        return sum(map(numeric_count, x))
    if isinstance(x, dict):
        return sum(numeric_count(v) for v in x.values())
    return 0


def bits(x):
    return max(abs(x.numerator).bit_length(), x.denominator.bit_length())


def parse(s):
    p = str(s).split('/')
    assert max(abs(int(v)).bit_length() for v in p) <= LIMIT
    q = Q(s)
    assert bits(q) <= LIMIT
    return q


def product(a, b):
    if not a or not b:
        return Q(0)
    assert max(abs(a.numerator).bit_length() + abs(b.numerator).bit_length(),
               a.denominator.bit_length() + b.denominator.bit_length()) <= LIMIT
    return a * b


def plus(a, b):
    if not a:
        return b
    if not b:
        return a
    assert max(max(abs(a.numerator).bit_length() + b.denominator.bit_length(),
                   abs(b.numerator).bit_length() + a.denominator.bit_length()) + 1,
               a.denominator.bit_length() + b.denominator.bit_length()) <= LIMIT
    return a + b


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-peer-report', required=True)
    args = parser.parse_args()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2147483648, 2147483648))
    started = time.monotonic()
    assert sha(ROOT / SOURCE) == SOURCE_SHA
    source = json.loads((ROOT / SOURCE).read_text())
    peer_path = Path(args.source_peer_report)
    if not peer_path.is_absolute():
        peer_path = ROOT / peer_path
    peer = json.loads(peer_path.read_text())
    assert peer['all_passed']
    # Review must explicitly bind this immutable PRIMARY SOURCE report.
    assert SOURCE_SHA in set(peer.get('input_sha256', {}).values())
    assert sha(ROOT / SOURCE_EVIDENCE) == source['evidence_sha256'][Path(SOURCE_EVIDENCE).name]
    assert sha(ROOT / C2) == source['input_sha256'][C2]
    evidence = json.loads((ROOT / SOURCE_EVIDENCE).read_text())
    c2 = json.loads((ROOT / C2).read_text())
    input_count = sum(map(numeric_count, [source, peer, evidence, c2]))
    # 6000 covers every parsed leaf copy, term record, saved product/result,
    # output metadata and source pin including simultaneous serialization.
    owned_bound = input_count + 6000
    assert owned_bound <= 18000 < 30000
    inputs_raw = sum(p.stat().st_size for p in [ROOT / SOURCE, peer_path,
                                               ROOT / SOURCE_EVIDENCE, ROOT / C2])
    assert inputs_raw + 100000 <= 400000
    out = ROOT / 'reports' / STEM
    out.mkdir(exist_ok=False)
    leaves = [{i: s for i, s in row['delta']} for row in c2['physical12_coefficients']]
    actual = []
    terms_total = 0
    for j, terms in enumerate(evidence['nine_cached_eA_cubic_adjoints']):
        values = []
        value = Q(0)
        for physical_index, name, coefficient_index, weight in terms:
            assert name == 'delta'
            x = parse(leaves[physical_index].get(coefficient_index, '0'))
            w = parse(weight)
            p = product(w, x)
            value = plus(value, p)
            values.append([physical_index, coefficient_index, str(w), str(x), str(p)])
        actual.append({'source_index': j,
                       'free_physical_row': evidence['physical_row_order'][evidence['free_columns9'][j]],
                       'all_term_products': values, 'value': str(value)})
        terms_total += len(terms)
    assert terms_total == evidence['worker_profile']['max_folded_products']
    nonzero = [x['source_index'] for x in actual if parse(x['value'])]
    result = {
        'all9_actual_cached_eA_cubic_dots': actual,
        'nonzero_source_indices': nonzero,
        'allzero_is_only_nine_projection_negative': not nonzero,
        'actual_C3_length_upper19_from_nonzero_cubic_and_ell_multiple': bool(nonzero),
        'actual_C3_length18_or_full_torsion1122_determined': False,
        'new_physical_implicit_gauge_or_jet_coefficients': 0,
        'owned_numeric_bound_including_output_and_metadata': owned_bound,
        'input_RAW_bytes': inputs_raw,
        'source_derived_dot_bits_bound': evidence['height_profile']['source_dot_bits_bound'],
        'every_QQ_operation_prechecked': True,
    }
    ev = out / 'all9-cached-cubic-dots-and-termwise-source-products.json'
    ev.write_text(json.dumps(result, indent=2) + '\n')
    report = {
        'id': STEM, 'all_passed': True,
        'status': 'POSITIVE_C3_UPPER19' if nonzero else 'ALLZERO_NINE_eA_PROJECTION_SUBSET_ONLY',
        'nonzero_count': len(nonzero),
        'actual_C3_length_upper19_proved': bool(nonzero),
        'actual_C3_length_or_all_torsion_exponents_determined': False,
        'input_sha256': {SOURCE: SOURCE_SHA, SOURCE_EVIDENCE: sha(ROOT / SOURCE_EVIDENCE),
                         C2: sha(ROOT / C2), str(peer_path.relative_to(ROOT)): sha(peer_path)},
        'code_sha256': sha(Path(__file__)),
        'note_sha256': sha(ROOT / 'research/notes' / (STEM + '.md')),
        'evidence_sha256': {ev.name: sha(ev)},
        'seconds': time.monotonic() - started,
    }
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    assert inputs_raw + sum(p.stat().st_size for p in out.iterdir()) <= 400000
    print(json.dumps(report))
    signal.alarm(0)


if __name__ == '__main__':
    main()

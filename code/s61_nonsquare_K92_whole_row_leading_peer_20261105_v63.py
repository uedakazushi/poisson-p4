#!/usr/bin/env python3
"""Independent direct bordered5 determinant, source bound before one QQ value."""
import argparse
import hashlib
import json
import resource
import signal
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_NONSQUARE_K92_WHOLE_ROW_LEADING_PEER_20261105_v63'
CAP = 6663
SOURCE = ['reports/S61_NONSQUARE_PREPOWER_K92_EXACT_ABC72_20261105_v63/ABC72.jsonl',
          'reports/S61_NONSQUARE_PREPOWER_K92_EXACT_D49_20261105_v63/D49.jsonl']
PLAN = 'reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_ROW_SOURCE_20261105_v63/whole-row-source-and-operator.json'
VALUE = 'reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_ROW_ENTRY_20261105_v63/whole-row-column0.json'

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def dump(p, data):
    p.write_text(json.dumps(data, indent=2) + '\n')

def inputs():
    # Streaming parser retains only the selected literal entries, below CAP.
    plan = json.loads((ROOT / PLAN).read_text())
    den = int(plan['uniform_integer_denominator'])
    literals = {}
    for name in SOURCE:
        with (ROOT / name).open() as stream:
            for line in stream:
                e = json.loads(line)
                block = e.get('block', 'D')
                i = e['block_row'] if block != 'D' else e['residual_row']
                j = e['block_column'] if block != 'D' else e['residual_column']
                wanted = block == 'A' or (block == 'B' and j == 0) or (block == 'C' and i in (0, 2)) or (block == 'D' and i in (0, 2) and j == 0)
                if not wanted:
                    continue
                content = Fraction(e['QQ_content'])
                assert den % content.denominator == 0
                scalar = content.numerator * (den // content.denominator)
                literals[block, i, j] = {int(k): int(v) * scalar for k, v in e['sparse_integer_polynomial'] if int(k) <= CAP}
    matrix = [[{} for _ in range(5)] for _ in range(5)]
    for i in range(4):
        for j in range(4):
            matrix[i][j] = literals['A', i, j]
        matrix[i][4] = literals['B', i, 0]
    for j in range(5):
        block, column = ('C', j) if j < 4 else ('D', 0)
        p = dict(literals[block, 2, column])
        for k, v in literals[block, 0, column].items():
            p[k] = p.get(k, 0) - 59375 * v
            if not p[k]:
                del p[k]
        matrix[4][j] = p
    return den, matrix

def determinant(matrix, absolute, source_gate=None):
    # Full subset determinant with permutation signs; no 2+2 fused cleanup.
    dp = {0: {0: 1}}
    peaks = []
    products = 0
    maxbits = 0
    for row in range(5):
        nxt = {}
        for mask, p in dp.items():
            for col in range(5):
                if mask & (1 << col):
                    continue
                q = matrix[row][col]
                if not q:
                    continue
                sign = 1 if absolute else (-1) ** ((mask >> (col + 1)).bit_count())
                result = nxt.setdefault(mask | (1 << col), {})
                for g, a in p.items():
                    for h, b in q.items():
                        if g + h > CAP:
                            continue
                        value = abs(a * b) if absolute else sign * a * b
                        result[g + h] = result.get(g + h, 0) + value
                        if not result[g + h]:
                            del result[g + h]
                        products += 1
                        if products % 4096 == 0:
                            owned = 47864 + 15000 + 2 * (sum(map(len, dp.values())) + sum(map(len, nxt.values())))
                            assert owned < 150000, ('owned', owned)
        slots = sum(map(len, dp.values())) + sum(map(len, nxt.values()))
        peaks.append(slots)
        maxbits = max(maxbits, max((abs(v).bit_length() for p in nxt.values() for v in p.values()), default=0))
        if source_gate:
            assert slots <= source_gate['simultaneous_DP_slots'][row]
            assert maxbits <= source_gate['absolute_coefficient_bits']
        dp = nxt
    return dp.get(31, {}).get(CAP, 0), peaks, products, maxbits

def main(out, mode):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    out.mkdir(parents=True, exist_ok=False)
    den, matrix = inputs()
    pins = {p: sha(ROOT / p) for p in SOURCE + [PLAN]}
    gatepath = ROOT / 'reports' / (STEM + '_SOURCE') / 'resource-gate.json'
    if mode == 'source':
        bound, peaks, products, bits = determinant(matrix, True)
        owned = 47864 + 15000 + 2 * max(peaks)
        data = dict(target_degree=CAP, absolute_target_upper=str(bound), simultaneous_DP_slots=peaks,
                    absolute_coefficient_bits=bits, absolute_product_pairs=products,
                    source_and_parser_allowance=47864, numeric_metadata_allowance=15000,
                    whole_owned_upper=owned, input_RAW=sum((ROOT / p).stat().st_size for p in SOURCE + [PLAN]),
                    future_new_RAW_upper=2000000, source_only_no_signed_target=True)
        assert owned < 150000 and data['input_RAW'] + 2000000 < 33554432
        dump(out / 'resource-gate.json', data)
    else:
        gate = json.loads(gatepath.read_text())
        pins[str(gatepath.relative_to(ROOT))] = sha(gatepath)
        n, peaks, products, bits = determinant(matrix, False, gate)
        expected = json.loads((ROOT / VALUE).read_text())
        pins[VALUE] = sha(ROOT / VALUE)
        actual = dict((int(k), int(v)) for k, v in expected['integer_polynomial'])
        assert n == actual[CAP] and n != 0
        coefficient = Fraction(n, den**5)
        assert coefficient == Fraction(17149069629000000000, 7)
        data = dict(direct5_determinant_coefficient_matches=True, degree=CAP,
                    QQ_coefficient=str(coefficient), coefficient_NONZERO=True,
                    all4131_coefficients_checked=False, no_full_det11_or_global57_claim=True,
                    actual_product_pairs=products, simultaneous_DP_slots=peaks, coefficient_bits=bits,
                    whole_owned_upper=gate['whole_owned_upper'])
        dump(out / 'direct5-coefficient-peer.json', data)
    pins[str(Path(__file__).relative_to(ROOT))] = sha(__file__)
    dump(out / 'report.json', dict(id=STEM + ('_SOURCE' if mode == 'source' else ''), all_passed=True,
         source_only=mode == 'source', input_sha256=pins,
         evidence_sha256={p.name: sha(p) for p in out.iterdir()}, seconds=time.monotonic()-started))
    print(json.dumps({k: v for k, v in data.items() if k != 'absolute_target_upper'}))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['source', 'value'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    main(args.output, args.mode)

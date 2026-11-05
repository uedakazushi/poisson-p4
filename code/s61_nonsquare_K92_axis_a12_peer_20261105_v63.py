#!/usr/bin/env python3
"""Independent binary-source restriction and subset determinant coefficient."""
import argparse
import hashlib
import json
import math
import resource
import signal
import struct
import time
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ID = 'S61_NONSQUARE_K92_AXIS_A12_PEER_20261105_v63'
STATE = 'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin'
BIND = 'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json'
AXIS = 'reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_SOURCE_V2_20261105_v63/whole-axis-literal-source.json'
FACE = 'reports/S61_NONSQUARE_PREPOWER_K92_ORIGINAL_TERMINAL_FACE_SOURCE_20261105_v63/terminal-face-and-normalization.json'
PRIMARY = 'reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_BAREISS_ENTRY_20261105_v63/axis-rank-and-determinant-normalization.json'
CAP = 7

def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def write(path, obj):
    path.write_text(json.dumps(obj, indent=2) + '\n')

def original_axis():
    # Implement the saved K92SP001 format here; import no producer routines.
    binding = json.loads((ROOT / BIND).read_text())
    roots = binding['actual_pure_matrix_roots']
    assert len(roots) == 11 and all(len(r) == 8 for r in roots)
    needed = {s for row in roots for s in row} | set(range(575, 593))
    by_source, tail = {}, {}
    kept = 0
    assert sha(ROOT / STATE) == 'cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02'
    with (ROOT / STATE).open('rb') as f:
        def get(fmt):
            size = struct.calcsize('<' + fmt)
            data = f.read(size)
            assert len(data) == size
            return struct.unpack('<' + fmt, data)[0]
        def integer(bits, retain):
            size, negative = get('H'), get('B')
            assert negative in (0, 1) and size <= (bits + 7) // 8
            if not retain:
                f.seek(size, 1)
                return None
            data = f.read(size)
            assert len(data) == size and (not size or data[0])
            value = int.from_bytes(data, 'big')
            assert value.bit_length() <= bits
            return -value if negative else value
        assert f.read(8) == b'K92SP001'
        size = get('I')
        assert size < 65536
        header = json.loads(f.read(size))
        assert header['source_completed'] == 572
        count = get('I')
        assert count <= 1024
        for _ in range(count):
            pid, role, src = get('I'), get('B'), get('i')
            term, row, mask, done = get('i'), get('i'), get('i'), get('B')
            retain = src in needed or (role == 2 and term == -20)
            numerator, denominator = integer(4096, retain), integer(4096, retain)
            nt = get('I')
            axis = {}
            if retain:
                assert done
                content = Q(numerator, denominator)
                assert content > 0
                kept += nt
            for _ in range(nt):
                a, e, d = get('i'), get('i'), get('i')
                if retain:
                    assert min(a, e, d) >= 0
                value = integer(1024, retain and e == 0 and d == 0)
                if retain and e == 0 and d == 0:
                    assert a not in axis and value
                    axis[a] = content * value
            if retain:
                if src >= 0:
                    by_source[src] = axis
                else:
                    assert row == 1 and mask in (0, 1, 2)
                    tail[mask] = axis
        if get('B'):
            for _ in range(3): get('i')
            for _ in range(2): get('Q')
            for _ in range(3): get('q')
            integer(16384, False)
            get('B')
            integer(16384, False)
        assert not f.read(1)
    assert len(by_source) == 106 and len(tail) == 3 and kept == 21690
    matrix = [[by_source[s] for s in row] + [{}, {}, {}] for row in roots]
    for i in range(6):
        for j in range(3):
            matrix[5+i][8+j] = tail[j] if i == 1 else by_source[575+3*i+j]
    for j in range(11):
        p = dict(matrix[1][j])
        for degree, value in matrix[2][j].items():
            p[degree] = p.get(degree, 0) + 2 * value
            if not p[degree]: del p[degree]
        matrix[1][j] = p
    return matrix

def matrix_and_source_check():
    original = original_axis()
    axis = json.loads((ROOT / AXIS).read_text())
    face = json.loads((ROOT / FACE).read_text())
    rp, cp = face['row_dual'], face['column_dual']
    assert sum(rp) + sum(cp) == -5
    lcd = axis['axis_row_denominators']
    matrix = [[{} for _ in range(11)] for _ in range(11)]
    for entry in axis['axis_coefficients']:
        i, j = entry['row'], entry['column']
        shift = rp[i] + cp[j]
        assert shift == entry['a_rowcol_shift']
        assert lcd[i] == entry['axis_row_integer_clear']
        p = {}
        for degree, q in original[i][j].items():
            value = q * lcd[i]
            assert degree + shift >= 0 and value.denominator == 1
            p[degree + shift] = value.numerator
        assert p == dict((int(k), int(v)) for k, v in entry['axis_integer_polynomial'])
        matrix[i][j] = p
    assert sum(len(p) for row in matrix for p in row) == 579
    return matrix, math.prod(lcd)

def coefficient(matrix, absolute, gate=None):
    dp, peaks, products, bits = {0: {0: 1}}, [], 0, 0
    for row in range(11):
        nxt = {}
        for mask, p in dp.items():
            for col, q in enumerate(matrix[row]):
                if mask & (1 << col) or not q: continue
                result = nxt.setdefault(mask | (1 << col), {})
                sign = 1 if absolute else (-1) ** ((mask >> (col+1)).bit_count())
                for g, a in p.items():
                    for h, b in q.items():
                        if g+h > CAP: continue
                        v = abs(a*b) if absolute else sign*a*b
                        result[g+h] = result.get(g+h, 0) + v
                        if not result[g+h]: del result[g+h]
                        products += 1
        slots = sum(map(len, dp.values())) + sum(map(len, nxt.values()))
        peaks.append(slots)
        bits = max(bits, max((abs(v).bit_length() for p in nxt.values() for v in p.values()), default=0))
        assert 30000 + 2*slots < 150000
        if gate:
            assert slots <= gate['DP_simultaneous_slots'][row]
            assert bits <= gate['absolute_bit_bound']
        dp = nxt
    return dp.get(2047, {}).get(CAP, 0), peaks, products, bits

def main(out, mode):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    out.mkdir(parents=True, exist_ok=False)
    matrix, denominator = matrix_and_source_check()
    pins = {s: sha(ROOT/s) for s in (STATE, BIND, AXIS, FACE)}
    if mode == 'source':
        bound, peaks, products, bits = coefficient(matrix, True)
        data = dict(all121_generic_axis_restrictions_verified=True, actual579_integer_coefficients_verified=True,
                    no_producer_parser_or_Bareiss_reused=True, target_Caxis_degree=7, target_original_degree=12,
                    source_parser_and_all_numeric_metadata_allowance=30000,
                    DP_simultaneous_slots=peaks, whole_owned_upper=30000+2*max(peaks),
                    absolute_target_bound=str(bound), absolute_bit_bound=bits,
                    product_pairs_upper=products, input_RAW=sum((ROOT/s).stat().st_size for s in pins),
                    future_new_RAW_upper=2000000, no_signed_target_coefficient_computed=True)
        assert data['input_RAW'] + 2000000 < 33554432
        write(out/'resource-gate.json', data)
    else:
        gp = ROOT/'reports'/(ID+'_SOURCE')/'resource-gate.json'
        gate = json.loads(gp.read_text())
        pins[str(gp.relative_to(ROOT))] = sha(gp)
        n, peaks, products, bits = coefficient(matrix, False, gate)
        value = Q(n, denominator)
        expected = dict((int(k), Q(v)) for k, v in json.loads((ROOT/PRIMARY).read_text())['axis_detM_e0d0_QQ_coefficients'])
        pins[PRIMARY] = sha(ROOT/PRIMARY)
        assert value == expected[12] == Q(488042709959022324503443014699328840826617856,17514302730560302734375)
        assert value != 0
        data = dict(all121_generic_axis_restrictions_verified=True, source_nonnegative_e_d_verified=True,
                    original_a12_coefficient=str(value), coefficient_NONZERO=True,
                    direct11_subset_determinant_matches=True, primary_Bareiss_divisions_not_reused=True,
                    all24_primary_coefficients_checked=False, actual_products=products,
                    actual_DP_slots=peaks, actual_bits=bits,
                    generic_original_M11_determinant_NONZERO=True,
                    no_fullU57_allpoint_or_badFitting_claim=True)
        write(out/'independent-source-and-a12-coefficient.json', data)
    pins[str(Path(__file__).relative_to(ROOT))] = sha(__file__)
    write(out/'report.json', dict(id=ID+('_SOURCE' if mode=='source' else '_REVIEW'), all_passed=True,
          source_only=mode=='source', input_sha256=pins,
          evidence_sha256={p.name: sha(p) for p in out.iterdir()}, seconds=time.monotonic()-started))
    print(json.dumps({k:v for k,v in data.items() if k!='absolute_target_bound'}))

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--mode', choices=['source','value'], required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    main(a.output, a.mode)

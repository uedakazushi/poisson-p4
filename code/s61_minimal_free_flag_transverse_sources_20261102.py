#!/usr/bin/env python3
"""Exact source-preserving transverse reduction; all assertions are over F_32003."""
import argparse
import hashlib
import inspect
import json
import subprocess
import tarfile
import tempfile
import time
from itertools import combinations, combinations_with_replacement
from pathlib import Path

import numpy as np
import s61_minimal_free_flag_initial_probe_20261101 as raw
from truncated_local_ideal_io_20261101 import write_input

ROOT = Path(__file__).resolve().parents[1]
P = 32003


def echelon(matrix):
    matrix = matrix.copy() % P
    row = 0
    pivots = []
    for col in range(matrix.shape[1]):
        choices = np.flatnonzero(matrix[row:, col])
        if not len(choices):
            continue
        other = row + int(choices[0])
        matrix[[row, other]] = matrix[[other, row]]
        matrix[row] = matrix[row] * pow(int(matrix[row, col]), -1, P) % P
        ids = np.flatnonzero(matrix[row + 1:, col]) + row + 1
        if len(ids):
            matrix[ids] = (matrix[ids] - matrix[ids, col, None] * matrix[row]) % P
        pivots.append(col)
        row += 1
        if row == matrix.shape[0]:
            break
    return matrix[:row], pivots


def run():
    start = time.time()
    archive = ROOT / 'research/wip/S61_MINIMAL_FREE_FLAG_R_ADIC_ACTIONS_20261102.tar.gz'
    with tarfile.open(archive) as tf:
        candidate = json.load(tf.extractfile('candidate.json'))
    old = json.loads((ROOT / 'reports/S61_MINIMAL_FREE_FLAG_R_ADIC_ACTIONS_20261102.json').read_text())
    border = old['minimal_border']
    mons = [m for d in range(9) for m in combinations_with_replacement(range(8), d)]
    def divides(a, b):
        return all(b.count(j) >= a.count(j) for j in set(a))
    basis = [m for m in mons if not any(divides(b, m) for b in border)]
    assert len(basis) == 141 and max(map(len, basis)) == 7
    matrices = []
    for action in candidate['actions'][:8]:
        matrix = np.zeros((141, 141), dtype=np.int64)
        for col, sparse in enumerate(action['columns']):
            for row, coefficient in sparse:
                if row < 141:
                    matrix[row, col] = coefficient
        matrices.append(matrix)
    d1 = np.hstack(matrices)
    d2 = np.zeros((1128, 3948), dtype=np.int64)
    for k, (i, j) in enumerate(combinations(range(8), 2)):
        d2[141*i:141*(i+1), 141*k:141*(k+1)] = -matrices[j] % P
        d2[141*j:141*(j+1), 141*k:141*(k+1)] = matrices[i]
    assert not np.any(d1 @ d2 % P)
    assert len(echelon(d1)[1]) == 140
    image, image_pivots = echelon(d2.T)
    assert len(image_pivots) == 949
    normal_forms = {(): np.eye(141, dtype=np.int64)[:, 0]}
    for degree in range(1, 8):
        for mon in combinations_with_replacement(range(8), degree):
            normal_forms[mon] = matrices[mon[-1]] @ normal_forms[mon[:-1]] % P
    for i, mon in enumerate(basis):
        assert np.array_equal(normal_forms[mon], np.eye(141, dtype=np.int64)[:, i])
    # The original 135 equations and the frozen frame, with only base variables set to zero.
    source = inspect.getsource(raw.build_jets).replace('n=T.cols', 'T=T[:,:8];n=T.cols')
    source = source.replace('def build_jets(', 'def transverse_jets(')
    namespace = vars(raw).copy()
    exec(source, namespace)
    n, jets, _ = namespace['transverse_jets'](8)
    assert n == 8
    with tempfile.TemporaryDirectory(prefix='poisson-source-reduction-') as temp:
        td = Path(temp)
        inp = td / 'raw8.txt'
        write_input(inp, P, 8, jets, None, maximum_degree=8)
        lines = inp.read_text().splitlines()
        inputs = []
        for line in lines[10:145]:
            z = list(map(int, line.split()))
            inputs.append(dict(zip(z[1::2], z[2::2])))
        cycles = np.zeros((135, 1128), dtype=np.int64)
        for i, poly in enumerate(inputs):
            for k, coefficient in poly.items():
                mon = mons[k]
                j = mon[0]
                cycles[i, 141*j:141*(j+1)] = (cycles[i, 141*j:141*(j+1)] + coefficient * normal_forms[mon[1:]]) % P
        assert not np.any(d1 @ cycles.T % P)
        for row, pivot in zip(image, image_pivots):
            ids = np.flatnonzero(cycles[:, pivot])
            cycles[ids] = (cycles[ids] - cycles[ids, pivot, None] * row) % P
        _, selected = echelon(cycles.T)
        assert len(selected) == 39 == 1128 - 140 - 949
        small = td / 'selected39.txt'
        header = lines[0].split()
        header[3] = '39'
        small.write_text('\n'.join([' '.join(header)] + lines[1:10] + [lines[10+i] for i in selected] + [lines[-1]]) + '\n')
        native = ['free_m0_transverse_source_trace_20261102', 'free_m0_transverse_source_replay_20261102']
        for stem in native:
            subprocess.run(['clang++', '-O3', '-std=c++17', str(ROOT / 'code' / (stem + '.cpp')), '-o', str(td / stem)], check=True, capture_output=True)
        trace, full, constant = [td / name for name in ['trace.txt', 'full.txt', 'constant.txt']]
        produced = subprocess.run([str(td / native[0]), str(small), str(td / 'initial.txt'), str(trace), str(full)], check=True, capture_output=True, text=True, timeout=60)
        replay = subprocess.run([str(td / native[1]), str(small), str(full), str(trace), str(constant)], check=True, capture_output=True, text=True, timeout=60)
        full_lines = full.read_text().splitlines()
        rows = []
        for line in full_lines[1:]:
            z = list(map(int, line.split()))
            rows.append(dict(zip(z[1::2], z[2::2])))
        assert len(rows) == 167
        leads = [mons[min(row)] for row in rows]
        minimal = [i for i, mon in enumerate(leads) if not any(k != i and (len(other) < len(mon) or (len(other) == len(mon) and k < i)) and divides(other, mon) for k, other in enumerate(leads))]
        assert len(minimal) == 89
        actual_basis = [mon for mon in mons if not any(divides(leads[i], mon) for i in minimal)]
        assert actual_basis == basis
        text = constant.read_text().splitlines()
        coefficient_matrix = np.array([list(map(int, line.split())) for line in text[1:]], dtype=np.int64)
        assert coefficient_matrix.shape == (167, 39)
        _, minor_columns = echelon(coefficient_matrix[minimal].T)
        assert len(minor_columns) == 39
        minor_rows = [minimal[i] for i in minor_columns]
        matrix = coefficient_matrix[minor_rows].copy()
        determinant = 1
        for i in range(39):
            other = i + int(np.flatnonzero(matrix[i:, i])[0])
            if other != i:
                matrix[[i, other]] = matrix[[other, i]]
                determinant = -determinant
            determinant = determinant * int(matrix[i, i]) % P
            matrix[i] = matrix[i] * pow(int(matrix[i, i]), -1, P) % P
            matrix[i+1:] = (matrix[i+1:] - matrix[i+1:, i, None] * matrix[i]) % P
        assert determinant == 9252
        result = dict(id='S61_MINIMAL_FREE_FLAG_TRANSVERSE_SOURCES_20261102', all_passed=True, prime=P,
                      transverse_length=141, transverse_top_degree=7, Koszul_d1_rank=140, Koszul_d2_rank=949,
                      minimal_transverse_generators=39, original135_indices=selected,
                      trace_nodes=167, minimal_monic_rewrites=89, source_constant_rank=39,
                      minimal_rows=minimal, source_minor_rows=minor_rows, source_minor_determinant=determinant,
                      source_constant_matrix=coefficient_matrix[minimal].tolist(),
                      fresh_raw8_sha256=hashlib.sha256(inp.read_bytes()).hexdigest(),
                      trace_sha256=hashlib.sha256(trace.read_bytes()).hexdigest(),
                      full_rows_sha256=hashlib.sha256(full.read_bytes()).hexdigest(),
                      source_replay_stdout=replay.stdout.strip(), native_stdout=produced.stdout.strip(),
                      refined_all_r_cutoff='m^(r+7) subset a^r for all r>=1',
                      limitations=['Only the specified finite-field point.', 'The 39 equations generate the transverse ideal; the other 96 must be retained in the relative ideal.', 'No full relative finite presentation, t-nilpotence, QQ lift, or exhaustion follows.'])
    paths = [archive, ROOT / 'reports/S61_MINIMAL_FREE_FLAG_R_ADIC_ACTIONS_20261102.json', ROOT / 'reports/S61_MINIMAL_FREE_FLAG_R_ADIC_CUTOFF_20261102.json', ROOT / 'reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json', Path(__file__), ROOT / 'code/s61_minimal_free_flag_initial_probe_20261101.py', ROOT / 'code/truncated_local_ideal_io_20261101.py']
    paths += [ROOT / 'code' / (stem + '.cpp') for stem in native]
    paths += [ROOT / 'research/notes/S61_MINIMAL_FREE_FLAG_TRANSVERSE_SOURCES_20261102.md']
    result['input_sha256'] = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    result['seconds'] = time.time() - start
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = run()
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print('PASS', report['seconds'])

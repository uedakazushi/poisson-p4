#!/usr/bin/env python3
"""Actual F_32003 S/a^4, via source rewrites and a 141-generator module."""
import argparse
import hashlib
import inspect
import json
import subprocess
import tarfile
import tempfile
import time
from itertools import combinations_with_replacement as cwr
from pathlib import Path

import numpy as np
import s61_minimal_free_flag_initial_probe_20261101 as raw
from truncated_local_ideal_io_20261101 import write_input
from s61_minimal_free_flag_transverse_sources_20261102 import echelon

ROOT = Path(__file__).resolve().parents[1]
P = 32003
STEM = 'S61_MINIMAL_FREE_FLAG_THIRD_BASE_QUOTIENT_20261102'
NATIVE = ['free_m0_transverse_source_trace_20261102', 'free_m0_transverse_source_replay_20261102',
          'free_m0_relative_source_lift_20261102', 'free_m0_relative_source_canonical_20261102',
          'free_m0_relative_module_actions_20261102', 'free_m0_relative_original_residuals_20261102',
          'free_m0_relative_module_closure_20261102']


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def generate_raw(path):
    source = inspect.getsource(raw.build_jets).replace('def build_jets(', 'def weighted_jets(')
    source = source.replace('list(combinations_with_replacement(range(n),d))', '[m for m in combinations_with_replacement(range(n),d) if sum(j>=8 for j in m)<4]')
    old = 'indices=np.array([idx[degree][tuple(sorted(ma+mb))]for ma in ms[da]for mb in ms[db]])'
    new = 'indices=np.array([idx[degree].get(tuple(sorted(ma+mb)),-1)for ma in ms[da]for mb in ms[db]]);mask=indices>=0'
    assert old in source
    source = source.replace(old, new).replace('np.bincount(indices,weights=cross.ravel(),', 'np.bincount(indices[mask],weights=cross.ravel()[mask],')
    source = source.replace('return n,jet,ts', 'return n,jet,ts,ms')
    namespace = vars(raw).copy()
    exec(source, namespace)
    n, jets, t, monomials = namespace['weighted_jets'](10)
    assert n == 11
    for degree in range(2, 11):
        full = list(cwr(range(11), degree))
        index = {m: i for i, m in enumerate(full)}
        selected = [index[m] for m in monomials[degree]]
        matrix = np.zeros((150, len(full)), dtype=np.int64)
        matrix[:135, selected] = jets[degree]
        observable = np.zeros(len(full), dtype=np.int64)
        observable[selected] = t[degree]
        if degree == 4:
            for i, mon in enumerate(cwr([8, 9, 10], 4)):
                matrix[135+i, index[mon]] = 1
        jets[degree], t[degree] = matrix, observable
    write_input(path, P, 11, jets, t, maximum_degree=10)


def transverse_input(raw_path, td, certificate):
    normal = [m for d in range(9) for m in cwr(range(8), d)]
    index = {m: i for i, m in enumerate(normal)}
    all_monomials = [m for d in range(11) for m in cwr(range(11), d)]
    projection = {i: index[m] for i, m in enumerate(all_monomials) if m in index}
    rows = []
    with raw_path.open() as handle:
        assert list(map(int, handle.readline().split())) == [P, 11, len(all_monomials), 150]
        for _ in range(12):
            handle.readline()
        for row in range(150):
            z = list(map(int, handle.readline().split()))
            assert len(z) == 2*z[0]+1
            if row < 135:
                rows.append({projection[i]: c for i, c in zip(z[1::2], z[2::2]) if i in projection})
    header = [f'{P} 8 {len(normal)} 135', ' '.join(str(len(m)) for m in normal)]
    header += [' '.join(str(index.get(tuple(sorted(m+(j,))), -1)) for m in normal) for j in range(8)]
    body = [str(len(row)) + ' ' + ' '.join(f'{i} {c}' for i, c in row.items()) for row in rows]
    ending = ' '.join('0' for _ in normal)
    fresh = td / 'transverse8.txt'
    fresh.write_text('\n'.join(header + body + [ending]) + '\n')
    assert sha(fresh) == certificate['fresh_raw8_sha256']
    header[0] = f'{P} 8 {len(normal)} 39'
    selected = td / 'selected39.txt'
    selected.write_text('\n'.join(header + [body[i] for i in certificate['original135_indices']] + [ending]) + '\n')
    selection = td / 'selection.txt'
    selection.write_text('39\n' + ' '.join(map(str, certificate['original135_indices'])) + '\n89\n' + ' '.join(map(str, certificate['minimal_rows'])) + '\n')
    return selected, selection


def run(archive, cached_raw=None):
    started = time.time()
    source_certificate_path = ROOT / 'reports/S61_MINIMAL_FREE_FLAG_TRANSVERSE_SOURCES_20261102.json'
    certificate = json.loads(source_certificate_path.read_text())
    times, logs = {}, {}
    with tempfile.TemporaryDirectory(prefix='poisson-free-m0-A4-') as temp:
        td = Path(temp)
        def stage(name, command, timeout):
            tick = time.time()
            done = subprocess.run(command, check=True, capture_output=True, text=True, timeout=timeout)
            times[name] = time.time() - tick
            logs[name] = done.stdout + done.stderr
            print(name, times[name], flush=True)
        for stem in NATIVE:
            stage('compile_' + stem, ['clang++', '-O3', '-std=c++17', str(ROOT / 'code' / (stem + '.cpp')), '-o', str(td / stem)], 60)
        raw_path = cached_raw or td / 'raw10.txt'
        if cached_raw is None:
            tick = time.time()
            generate_raw(raw_path)
            times['fresh_original135_raw10'] = time.time() - tick
            print('fresh_original135_raw10', times['fresh_original135_raw10'], flush=True)
        else:
            times['fresh_original135_raw10'] = None
        raw_digest = sha(raw_path)
        selected, selection = transverse_input(raw_path, td, certificate)
        trace, full, constant = [td / name for name in ['trace.txt', 'full.txt', 'constant.txt']]
        stage('transverse_source_trace', [str(td / NATIVE[0]), str(selected), str(td / 'initial.txt'), str(trace), str(full)], 60)
        stage('independent_source_replay', [str(td / NATIVE[1]), str(selected), str(full), str(trace), str(constant)], 60)
        assert sha(trace) == certificate['trace_sha256'] and sha(full) == certificate['full_rows_sha256']
        lifted, metadata = td / 'lifted.bin', td / 'lifted-metadata.txt'
        stage('relative_source_lift', [str(td / NATIVE[2]), str(raw_path), str(trace), str(selection), str(full), str(lifted), str(metadata)], 180)
        canonical, actions, residuals, closure = [td / name for name in ['canonical.txt', 'actions.txt', 'residuals.txt', 'closure.txt']]
        stage('canonical_source_rewrites', [str(td / NATIVE[3]), str(lifted), str(selection), str(canonical)], 180)
        text = canonical.read_text().splitlines()
        assert list(map(int, text[0].split())) == [P, 20, 89, 141]
        constants = np.array([list(map(int, line.split())) for line in text[3+20*89:]], dtype=np.int64)
        assert constants.shape == (89, 89)
        # The division matrix is triangular in leading-monomial order, with unit diagonal.
        lead = list(map(int, text[2].split()))
        order = np.argsort(lead)
        triangular = constants[np.ix_(order, order)]
        assert np.array_equal(np.diag(triangular), np.ones(89))
        assert not np.any(np.tril(triangular, -1))
        stage('eight_normal_actions', [str(td / NATIVE[4]), str(canonical), str(actions)], 120)
        stage('all135_and_observable_remainders', [str(td / NATIVE[5]), str(canonical), str(raw_path), str(residuals)], 120)
        stage('actual_module_closure_and_cutoff', [str(td / NATIVE[6]), str(actions), str(residuals), str(canonical), str(closure)], 120)
        log = logs['actual_module_closure_and_cutoff']
        assert 'CUTOFF checked 223158 added 0 final rank 294' in log
        assert 'CLOSED rank 294 dim 2526' in log and 'T2 support 0' in log and 'LOW A3 dim 1303' in log
        rows = closure.read_text().splitlines()
        assert list(map(int, rows[0].split())) == [P, 2820, 294]
        relations = np.zeros((294, 2820), dtype=np.int64)
        for i, line in enumerate(rows[1:295]):
            z = list(map(int, line.split()))
            relations[i, z[1::2]] = z[2::2]
        powers = [m for d in range(4) for m in cwr(range(3), d)]
        pi = {m: i for i, m in enumerate(powers)}
        shifts = np.zeros((882, 2820), dtype=np.int64)
        for j in range(3):
            for a, mon in enumerate(powers):
                if len(mon) < 3:
                    target = pi[tuple(sorted(mon+(j,)))]
                    shifts[294*j:294*(j+1), 141*target:141*(target+1)] = relations[:, 141*a:141*(a+1)]
        shift_basis, shift_pivots = echelon(shifts)
        assert len(shift_pivots) == 244
        classes = relations.copy()
        for row, pivot in zip(shift_basis, shift_pivots):
            ids = np.flatnonzero(classes[:, pivot])
            classes[ids] = (classes[ids] - classes[ids, pivot, None] * row) % P
        relation_orders = [len(powers[int(np.flatnonzero(row)[0])//141]) for row in relations]
        order = sorted(range(294), key=lambda i: -relation_orders[i])
        _, selected_minimal = echelon(classes[order].T)
        minimal = [order[i] for i in selected_minimal]
        assert len(minimal) == 50
        orders = [len(powers[int(np.flatnonzero(relations[i])[0])//141]) for i in minimal]
        counts = {str(d): orders.count(d) for d in [1, 2, 3]}
        assert counts == {'1': 22, '2': 25, '3': 3}
        data = dict(prime=P, coefficient_powers=powers, minimal_relation_rows=minimal, minimal_relation_orders=counts)
        (td / 'module_metadata.json').write_text(json.dumps(data, indent=2) + '\n')
        members = [canonical, actions, residuals, closure, constant, metadata, td / 'module_metadata.json']
        manifest = {path.name: sha(path) for path in members}
        (td / 'MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
        archive.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive, 'w:gz') as tf:
            for path in [td / 'MANIFEST.json'] + members:
                tf.add(path, arcname=path.name)
        generated = dict(raw10_sha256=raw_digest, source_trace_sha256=sha(trace), lifted_source_sha256=sha(lifted), compact_data_sha256=manifest)
    paths = [Path(__file__), source_certificate_path, ROOT / 'code/s61_minimal_free_flag_transverse_sources_20261102.py', ROOT / 'code/s61_minimal_free_flag_initial_probe_20261101.py', ROOT / 'code/truncated_local_ideal_io_20261101.py', ROOT / 'reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json', ROOT / 'reports/S61_MINIMAL_FREE_FLAG_R_ADIC_ACTIONS_20261102.json', ROOT / 'reports/S61_MINIMAL_FREE_FLAG_R_ADIC_CUTOFF_20261102.json', ROOT / 'research/notes' / (STEM + '.md'), archive]
    paths += [ROOT / 'code' / (stem + '.cpp') for stem in NATIVE]
    return dict(id=STEM, all_passed=True, prime=P, coefficient_ring='F_32003[u0,u1,u2]/(u)^4', generators=141,
                coefficient_columns=2820, relation_vector_space_rank=294, dimension=2526, previous_dimension=1303,
                third_base_layer_dimension=1223, minimal_module_relations=50, minimal_relation_orders=counts,
                base_times_kernel_rank=244, normal_actions=8, observable_actions=1,
                total_degree11_monomials_checked=223158, additional_cutoff_relations=0, maximal_ideal_power11_zero=True,
                t_nonzero=True, t_squared_zero=True, canonical_source_change_determinant=1,
                generated_data=generated, fresh_original_raw10=cached_raw is None, stage_seconds=times, logs=logs,
                limitations=['Only the specified F_32003 point and the finite quotient S/a^4.', 'No QQ lift, full-S t-nilpotence, relative ambient removal, or exhaustion.'],
                input_sha256={str(path.relative_to(ROOT)): sha(path) for path in paths}, seconds=time.time()-started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--cached-raw', type=Path)
    args = parser.parse_args()
    report = run(args.archive.resolve(), args.cached_raw)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print('PASS', report['seconds'])

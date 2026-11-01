#!/usr/bin/env python3
"""Validate a frozen WIP archive; do not rerun its native research experiments."""
import argparse
import hashlib
import json
import tarfile
import time
from itertools import combinations_with_replacement as cwr
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = 'research/wip/S61_MINIMAL_FREE_FLAG_R_ADIC_PREPARATION_20261102.tar.gz'
PRIME = 32003


def sparse_rows(data):
    lines = data.decode().splitlines()
    assert int(lines[0]) == len(lines) - 1
    result = []
    for line in lines[1:]:
        row = list(map(int, line.split()))
        assert len(row) == 2 * row[0] + 1
        value = dict(zip(row[1::2], row[2::2]))
        assert len(value) == row[0] and all(0 < x < PRIME for x in value.values())
        result.append(value)
    return result


def modular_rank(matrix):
    matrix = matrix.copy() % PRIME
    row = 0
    for column in range(matrix.shape[1]):
        candidates = np.flatnonzero(matrix[row:, column])
        if not len(candidates):
            continue
        k = row + int(candidates[0])
        matrix[[row, k]] = matrix[[k, row]]
        matrix[row] = matrix[row] * pow(int(matrix[row, column]), -1, PRIME) % PRIME
        targets = np.flatnonzero(matrix[row + 1:, column]) + row + 1
        if len(targets):
            matrix[targets] = (matrix[targets] - matrix[targets, column, None] * matrix[row]) % PRIME
        row += 1
        if row == matrix.shape[0]:
            break
    return row


def staircase(data, n, degree):
    monomials = [m for d in range(degree + 1) for m in cwr(range(n), d)]
    rows = sparse_rows(data)
    lead = [monomials[min(row)] for row in rows]
    assert all(len({len(monomials[i]) for i in row}) == 1 for row in rows)
    exponents = [tuple(m.count(i) for i in range(n)) for m in lead]
    minimal = []
    for exponent in sorted(set(exponents), key=lambda x: (sum(x), x)):
        if not any(all(a <= b for a, b in zip(old, exponent)) for old in minimal):
            minimal.append(exponent)
    levels = [[()]]
    for d in range(1, degree + 2):
        candidates = {tuple(sorted(m + (i,))) for m in levels[-1] for i in range(n)}
        standard = []
        for m in sorted(candidates):
            exponent = tuple(m.count(i) for i in range(n))
            if not any(all(a <= b for a, b in zip(old, exponent)) for old in minimal):
                standard.append(m)
        levels.append(standard)
        if not standard:
            break
    assert not levels[-1], 'The saved leading ideal has not established an Artin cutoff.'
    return dict(hilbert=[len(level) for level in levels[:-1]],
                length=sum(map(len, levels)), nilpotence_degree=len(levels) - 1,
                initial_rows=len(rows)), [m for level in levels for m in level]


def columns(rows):
    support = sorted(set().union(*(row.keys() for row in rows)))
    index = {value: i for i, value in enumerate(support)}
    matrix = np.zeros((len(support), len(rows)), dtype=np.int64)
    for j, row in enumerate(rows):
        for i, value in row.items():
            matrix[index[i], j] = value
    return matrix


def relation_matrix(rows, height):
    matrix = np.zeros((height, len(rows)), dtype=np.int64)
    for j, row in enumerate(rows):
        for i, value in row:
            matrix[i, j] = value
    return matrix


def verify():
    started = time.time()
    with tarfile.open(ROOT / ARCHIVE, 'r:gz') as archive:
        data = {}
        for member in archive.getmembers():
            assert member.isfile() and not member.name.startswith('/')
            assert '..' not in Path(member.name).parts and member.name not in data
            data[member.name] = archive.extractfile(member).read()
    manifest = json.loads(data.pop('MANIFEST.json'))
    assert set(data) == set(manifest)
    for name, content in data.items():
        assert len(content) == manifest[name]['bytes']
        assert hashlib.sha256(content).hexdigest() == manifest[name]['sha256']

    a0, basis = staircase(data['tmp/minimal-free-flag-transverse10-initial.txt'], 8, 10)
    a1, _ = staircase(data['tmp/minimal-free-flag-base-linear-initial.txt'], 11, 8)
    a2, _ = staircase(data['tmp/minimal-free-flag-base-quadratic-initial.txt'], 11, 9)
    assert [a0['length'], a1['length'], a2['length']] == [141, 542, 1303]
    assert [a0['nilpotence_degree'], a1['nilpotence_degree'], a2['nilpotence_degree']] == [8, 9, 10]
    saved_basis = json.loads(data['work/minimal-free-flag-module-query-basis.json'])['base_monomials']
    assert basis == list(map(tuple, saved_basis))

    first = json.loads(data['work/minimal-free-flag-first-base-relations.json'])
    second = json.loads(data['work/minimal-free-flag-second-base-relations.json'])
    answers1 = sparse_rows(data['tmp/minimal-free-flag-module-answers.txt'])
    answers2 = sparse_rows(data['tmp/minimal-free-flag-module-quadratic-answers.txt'])
    assert len(answers1) == 564 and len(answers2) == 1410
    matrix1 = columns(answers1[141:])
    matrix2 = columns(answers2[564:])
    n1 = relation_matrix(first['linear_relations'], 423)
    n2 = relation_matrix(second['quadratic_relations'], 846)
    assert modular_rank(matrix1) == first['first_u_rank'] == 401
    assert modular_rank(matrix2) == second['second_u_rank'] == 761
    assert n1.shape[1] == 22 and n2.shape[1] == 85
    assert modular_rank(n1) == 22 and modular_rank(n2) == 85
    assert not np.any(matrix1 @ n1 % PRIME)
    assert not np.any(matrix2 @ n2 % PRIME)
    # Every first-layer relation multiplied by one base variable appears
    # among second-layer relations. There are sixty independent multiples.
    powers = list(cwr(range(3), 2))
    multiples = np.zeros((846, 66), dtype=np.int64)
    for j in range(22):
        for k in range(3):
            for ell in range(3):
                position = powers.index(tuple(sorted((k, ell))))
                multiples[141 * position:141 * (position + 1), 3 * j + k] += n1[141 * ell:141 * (ell + 1), j]
    multiples %= PRIME
    assert not np.any(matrix2 @ multiples % PRIME)
    assert modular_rank(multiples) == 60
    assert '(MODULE, 141, 22, 125, 3, 125)' in first['m2']
    assert '(MODULE, 141, 107, 108, 3, 108)' in second['m2']
    t_record = json.loads(data['work/minimal-free-flag-parameter-t-action.json'])
    assert t_record['t_terms'] == 1053 and t_record['t2_terms'] == t_record['t3_terms'] == 0

    paths = [ARCHIVE, str(Path(__file__).relative_to(ROOT)),
             'research/notes/S61_MINIMAL_FREE_FLAG_R_ADIC_PREPARATION_20261102.md',
             'code/truncated_local_module_basis_20261102.cpp',
             'code/truncated_local_module_normal_form_20261102.cpp',
             'code/s61_minimal_free_flag_initial_probe_20261101.py',
             'code/truncated_local_standard_basis_20261101.cpp',
             'code/truncated_local_ideal_io_20261101.py',
             'code/s61_model.py', 'code/s61_flag_transverse_obstruction_20261101.py',
             'reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json',
             'reports/S61_MINIMAL_FREE_FLAG_FLAT4_PROBE_20261101.json',
             'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_MINIMAL_FREE_FLAG_R_ADIC_PREPARATION_20261102', all_passed=True,
                status='work in progress; validator of saved artifacts, not a new ambient theorem',
                prime=PRIME, base_tangent_columns=[8, 9, 10],
                saved_initial_quotients=[a0, a1, a2], parameter_adic_layer_dimensions=[141, 401, 761],
                linear_relations=22, quadratic_relations=85, independent_linear_multiples=60,
                new_quadratic_relations=25, recorded_polynomial_module_ranks=[125, 108],
                recorded_t_square_mod_base_cubed_zero=True,
                archive_members=len(manifest) + 1, archive_bytes=(ROOT / ARCHIVE).stat().st_size,
                bounded_failures=[dict(task='third parameter-adic layer, original ordering', seconds=240),
                                  dict(task='third parameter-adic layer, parameter generators first', seconds=180),
                                  dict(task='Artin local Mora with base squared', seconds=90)],
                verified_here=['archive SHA256', 'saved leading-ideal Artin staircases',
                               'normal-form matrix ranks', 'all saved relation kernels and sixty multiples'],
                not_reexecuted_here=['raw135 jets', 'native standard basis', 'M2 function-field module ranks',
                                     't-square normal-form reduction'],
                limitations=['No characteristic-zero or relative conclusion.',
                             'The module with relations through degree two still has recorded rank108, not the target70.',
                             'No full initial-ideal equality, actual torsion recognition, t-nilpotence or saturation removal.'],
                input_sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
                seconds=time.time() - started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print('FREE M0 R-ADIC WIP SNAPSHOT PASS', result['seconds'])

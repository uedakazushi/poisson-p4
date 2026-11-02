#!/usr/bin/env python3
"""Verify compressed 141-generator module and nine actions over Fp[u]/(u)^3."""
import argparse
import hashlib
import json
import tarfile
import time
from itertools import combinations_with_replacement as cwr
from pathlib import Path
import numpy as np
from s61_minimal_free_flag_r_adic_preparation_20261102 import modular_rank, sparse_rows

ROOT = Path(__file__).resolve().parents[1]
OLD = 'research/wip/S61_MINIMAL_FREE_FLAG_R_ADIC_PREPARATION_20261102.tar.gz'
DATA = 'research/wip/S61_MINIMAL_FREE_FLAG_R_ADIC_ACTIONS_20261102.tar.gz'
P = 32003


def read_archive(path):
    with tarfile.open(ROOT / path, 'r:gz') as archive:
        data = {m.name: archive.extractfile(m).read() for m in archive.getmembers() if m.isfile()}
    manifest = json.loads(data.pop('MANIFEST.json'))
    assert set(data) == set(manifest)
    for name, value in data.items():
        assert len(value) == manifest[name]['bytes']
        assert hashlib.sha256(value).hexdigest() == manifest[name]['sha256']
    return data


def verify():
    started = time.time()
    old, data = read_archive(OLD), read_archive(DATA)
    saved = json.loads(data['candidate.json'])
    C = sparse_rows(old['tmp/minimal-free-flag-module-quadratic-answers.txt'])
    answers = sparse_rows(data['answers.txt'])
    queries = sparse_rows(data['queries.txt'])
    B = list(map(tuple, json.loads(old['work/minimal-free-flag-module-query-basis.json'])['base_monomials']))
    assert len(C) == 1410 and len(B) == 141 and len(answers) == len(queries) == 1269
    powers = [()] + list(cwr(range(3), 1)) + list(cwr(range(3), 2))
    index = {power: i for i, power in enumerate(powers)}
    assert [list(v) for v in powers] == saved['coefficient_powers']

    def axpy(v, w, a):
        for k, c in w.items():
            value = (v.get(k, 0) + a*c) % P
            if value:
                v[k] = value
            else:
                v.pop(k, None)
        return v

    def encode(v):
        result = {}
        for k, c in v.items():
            axpy(result, C[k], c)
        return result

    def shift(v, power):
        result = {}
        for k, c in v.items():
            product = powers[k // 141] + power
            if len(product) <= 2:
                result[141*index[tuple(sorted(product))] + k % 141] = c
        return result

    def action(M, v):
        result = {}
        for k, c in v.items():
            axpy(result, shift(M[k % 141], powers[k // 141]), c)
        return result

    relations = [dict(r) for r in saved['lifted_linear_relations'] + saved['new_quadratic_relations']]
    assert len(relations) == 47
    assert all(all(k >= 141 for k in r) for r in relations[:22])
    assert all(all(k >= 564 for k in r) for r in relations[22:])
    old_linear = json.loads(old['work/minimal-free-flag-first-base-relations.json'])['linear_relations']
    for previous, lifted in zip(old_linear, relations[:22]):
        assert {k: c for k, c in lifted.items() if k < 564} == {i + 141: c for i, c in previous}
    closure = relations + [shift(r, (j,)) for r in relations[:22] for j in range(3)]
    matrix = np.zeros((1410, len(closure)), dtype=np.int64)
    for j, relation in enumerate(closure):
        assert not encode(relation)
        for i, c in relation.items():
            matrix[i, j] = c
    assert modular_rank(matrix) == 107
    # The previously verified normal-form image has dimension1303. Thus
    # these stable A3 relations are the whole kernel, not just a subspace.
    assert 1410 - 107 == saved['dimension'] == 1303
    A = [[dict(column) for column in row['columns']] for row in saved['actions']]
    assert len(A) == 9 and all(len(M) == 141 for M in A)
    for i, M in enumerate(A):
        for j, column in enumerate(M):
            assert encode(column) == answers[141*i + j]

    monomials = [m for d in range(10) for m in cwr(range(11), d)]
    mi = {m: i for i, m in enumerate(monomials)}
    for i in range(8):
        for j, b in enumerate(B):
            assert queries[141*i + j] == {mi[tuple(sorted(b + (i,)))]: 1}
    old_log = json.loads(old['work/minimal-free-flag-base-quadratic-module.json'])['log']
    observable = {int(a): int(b) for a, b in (word.split(':') for word in old_log.split('\nt ')[1].split())}
    for j, b in enumerate(B):
        expected = {mi[tuple(sorted(monomials[k] + b))]: c for k, c in observable.items()
                    if len(monomials[k]) + len(b) <= 9}
        assert queries[1128 + j] == expected

    stable = 0
    for M in A:
        for relation in relations:
            assert not encode(action(M, relation))
            stable += 1
    for j in range(141):
        assert not encode(action(A[8], A[8][j]))
    commutators = 0
    for i in range(9):
        for j in range(i):
            for k in range(141):
                value = action(A[i], A[j][k])
                axpy(value, action(A[j], A[i][k]), -1)
                assert not encode(value)
                commutators += 1
        print('ACTION COMMUTATORS', i, time.time() - started, flush=True)
    bindex = {m: i for i, m in enumerate(B)}
    for j, mon in enumerate(B):
        if mon:
            previous = bindex[mon[:-1]]
            assert encode(A[mon[-1]][previous]) == C[j]
    border = {tuple(sorted(mon + (i,))) for mon in B for i in range(8)} - set(B)
    minimal = [mon for mon in border if all(mon[:i] + mon[i+1:] in bindex for i in range(len(mon)))]
    assert len(border) == 537 and len(minimal) == 89
    paths = [OLD, DATA, str(Path(__file__).relative_to(ROOT)),
             'research/notes/S61_MINIMAL_FREE_FLAG_R_ADIC_ACTIONS_20261102.md',
             'code/s61_minimal_free_flag_r_adic_preparation_20261102.py',
             'code/truncated_local_module_normal_form_20261102.cpp',
             'reports/S61_MINIMAL_FREE_FLAG_R_ADIC_PREPARATION_20261102.json',
             'reports/S61_MINIMAL_FREE_FLAG_R_ADIC_CUTOFF_20261102.json']
    return dict(id='S61_MINIMAL_FREE_FLAG_R_ADIC_ACTIONS_20261102', all_passed=True,
                prime=P, coefficient_ring='Fp[u0,u1,u2]/(u0,u1,u2)^3', generators=141,
                linear_relations_with_quadratic_tails=22, new_quadratic_relations=25,
                coefficient_vector_dimension=1410, relation_kernel_dimension=107, dimension=1303,
                normal_form_unit_columns=1284, residual_inverse_size=19,
                quadratic_image_residual_inverse_size=24,
                normal_variable_actions=8, observable_actions=1,
                relation_stability_checks=stable, commutator_checks=commutators,
                t_squared_checks=141, t_squared_zero_in_this_quotient=True,
                distinct_border_monomials=537, minimal_border_monomials=89,
                next_base_coefficient_unknowns_per_monomial=89*141,
                minimal_border=[list(m) for m in sorted(minimal, key=lambda m: (len(m), m))],
                data_archive_bytes=(ROOT / DATA).stat().st_size,
                reduction_source=json.loads(data['reduction_source.json']),
                limitations=['The saved finite-field full local basis and its normal forms are the source-identification dependency.',
                             'This verifier checks module kernels, all nine actions and all commutators; it does not rebuild raw135 jets or the full basis.',
                             'Base-degree-three lifting still needs tracked membership in the original135 equations.',
                             'No full-S t-nilpotence, characteristic-zero conclusion or relative saturation removal.'],
                input_sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
                seconds=time.time() - started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = verify()
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print('FREE M0 R-ADIC ACTIONS PASS', report['seconds'])

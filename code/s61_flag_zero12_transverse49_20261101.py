#!/usr/bin/env python3
"""Uniform actual ZERO12 transverse length49 via polynomial filtered row operations."""
import argparse
import hashlib
import json
import time
from itertools import combinations_with_replacement
from pathlib import Path

import sympy as s
from sympy.polys.matrices import DomainMatrix

from s61_flag_zero12_normal_six_jets_20261101 import normal_jets

ROOT = Path(__file__).resolve().parents[1]


def filtered_certificate(cache):
    start = time.time()
    parameters = list(map(s.Symbol, cache['parameters']))
    K = s.QQ.frac_field(*parameters)
    P, PR = K.get_ring(), K.field.ring
    n, maximum_degree = 6, 6
    zero = (0,) * len(parameters)
    mons = {d: list(combinations_with_replacement(range(n), d)) for d in range(2, 7)}
    counts = {d: len(mons[d]) for d in mons}
    dok, offset = {}, 0
    for degree in range(2, 7):
        for i, row in enumerate(cache['A' + str(degree)]):
            for j, polynomial in enumerate(row):
                if polynomial:
                    dok[i, offset + j] = PR.from_dict({tuple(e): s.QQ(c) for e, c in polynomial})
        offset += counts[degree]
    current = DomainMatrix.from_dok(dok, (135, offset), P)
    records = []
    for degree in range(2, 7):
        count = counts[degree]
        lead = current.extract(range(current.shape[0]), range(count))
        if degree == 6:
            prefix = records[-1]['rank'] * n
            lower = lead.extract(range(prefix), range(count))
            assert all(not any(any(mon) for mon in value) for value in lower.to_dok().values())
            numeric = DomainMatrix.from_dok({ij: value[zero] for ij, value in lower.to_dok().items()},
                                           lower.shape, s.QQ)
            reduced, pivots = numeric.rref()
            assert len(pivots) == 455
            nonpivots = [j for j in range(count) if j not in pivots]
            projected = (lead.extract(range(lead.shape[0]), nonpivots)
                         - lead.extract(range(lead.shape[0]), pivots)
                         * reduced.extract(range(len(pivots)), nonpivots).convert_to(P))
            _, extra_columns = projected.to_field().rref()
            _, extra_rows = projected.to_field().transpose().rref()
            assert len(extra_columns) == len(extra_rows) == 7
            block = projected.extract(extra_rows, extra_columns)
            determinant = block.det()
            assert set(determinant) == {zero} and determinant[zero] != 0
            inverse = block.to_field().inv().convert_to(P)
            assert inverse * block == DomainMatrix.eye(7, P)
            records.append(dict(degree=6, rank=462, lower_rank=455, extra_rank=7,
                                quotient_dimension=0, extra_source_rows=list(extra_rows),
                                extra_columns=list(extra_columns), quotient_columns=nonpivots,
                                constant_extra_minor=str(determinant[zero]),
                                polynomial_inverse_verified=True))
            print(json.dumps(dict(stage='final_constant_minor', seconds=time.time() - start)), flush=True)
            break

        reduced, pivots = lead.to_field().rref()
        rank = len(pivots)
        assert rank == {2: 8, 3: 41, 4: 118, 5: 246}[degree]
        head = reduced.extract(range(rank), range(count))
        assert all(not any(any(m) for m in value.numer)
                   and not any(any(m) for m in value.denom) for value in head.to_dok().values())
        _, source_rows = lead.to_field().transpose().rref()
        block = lead.extract(source_rows, pivots)
        inverse = block.to_field().inv().convert_to(P)
        assert inverse * block == DomainMatrix.eye(rank, P)
        basis = inverse * current.extract(source_rows, range(current.shape[1]))
        assert basis.extract(range(rank), range(count)) == head.convert_to(P)
        coefficients = lead.extract(range(lead.shape[0]), pivots)
        difference = current - coefficients * basis
        assert difference.extract(range(difference.shape[0]), range(count)).is_zero_matrix
        tail = difference.extract(range(difference.shape[0]), range(count, difference.shape[1]))
        active_rows = sorted({i for i, _ in tail.to_dok()})
        tail = tail.extract(active_rows, range(tail.shape[1]))
        source_mons = sum([mons[d] for d in range(degree, 7)], [])
        target_mons = sum([mons[d] for d in range(degree + 1, 7)], [])
        index = {mon: i for i, mon in enumerate(target_mons)}
        multiples = {}
        for (i, j), coefficient in basis.to_dok().items():
            mon = source_mons[j]
            if len(mon) == maximum_degree:
                continue
            for variable in range(n):
                multiples[i * n + variable, index[tuple(sorted(mon + (variable,)))]] = coefficient
        multiplied = DomainMatrix.from_dok(multiples, (rank * n, len(target_mons)), P)
        constant_head = [[i, j, str(value.numer[zero] / value.denom[zero])]
                         for (i, j), value in sorted(head.to_dok().items())]
        records.append(dict(degree=degree, rank=rank, source_row_count=current.shape[0],
                            quotient_dimension=count - rank, source_rows=list(source_rows),
                            pivot_columns=list(pivots), constant_initial_basis_sparse=constant_head,
                            polynomial_inverse_verified=True, full_filtered_row_identity_verified=True,
                            retained_tail_rows=active_rows))
        current = multiplied.vstack(tail)
        print(json.dumps(dict(stage='filtered_degree', degree=degree, rank=rank,
                              next_shape=current.shape, seconds=time.time() - start)), flush=True)
    return records


def run():
    start = time.time()
    low = json.loads((ROOT / 'reports/S61_FLAG_ZERO12_LOW_JETS_20261101.json').read_text())
    assert low['all_passed']
    for path, expected in low['input_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == expected, path
    cache = normal_jets()
    print(json.dumps(dict(stage='raw135_normal_six_jets', seconds=time.time() - start)), flush=True)
    records = filtered_certificate(cache)
    paths = [Path(__file__), ROOT / 'code/s61_flag_zero12_normal_six_jets_20261101.py',
             ROOT / 'code/parametric_implicit_jets_20261101.cpp', ROOT / 'code/s61_model.py',
             ROOT / 'code/s61_flag_transverse_obstruction_20261101.py',
             ROOT / 'reports/S61_FLAG_ZERO12_LOW_JETS_20261101.json',
             ROOT / 'reports/S61_FLAG_UNIFORM_OPEN_20261101.json',
             ROOT / 'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',
             ROOT / 'reports/S61_FLAG_ZERO12_RANK_FIVE_20261101.json']
    return dict(id='S61_FLAG_ZERO12_TRANSVERSE49_20261101', all_passed=True,
                parameters=cache['parameters'], normal_columns=cache['normal_columns'],
                parameter_columns=cache['parameter_columns'], implicit_jacobian_minor=cache['inverse_minor'],
                independent_low_jet_restriction_verified=True, maximum_degree=6,
                residual_coefficient_term_counts={str(d): sum(len(p) for row in cache['A' + str(d)] for p in row)
                                                   for d in range(2, 7)},
                filtered_certificates=records, transverse_hilbert_function=[1, 6, 13, 15, 8, 6],
                actual_transverse_length=49, transverse_maximal_ideal_power6_zero=True,
                all_seven_parameter_values=True, transverse_family_finite_free_rank=49,
                ambient_finite_free_rank49_claimed=False, ambient_t_nilpotence_proved=False,
                limitations=['The quotient sets the three actual support coordinates v0,v1,v3 to zero.',
                             'No full local-ring rank49 or ambient t-nilpotence follows.',
                             'No whole-boundary saturation removal or full affine exhaustion is claimed.'],
                proof_note='research/notes/S61_FLAG_ZERO12_TRANSVERSE49_20261101.md',
                native_log=cache['native_log'],
                input_sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in paths}, seconds=time.time() - start)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = run()
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print(json.dumps(dict(all_passed=True, seconds=report['seconds'])))

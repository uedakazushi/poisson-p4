#!/usr/bin/env python3
"""Independent Cauchy-character Euler calculation and AB Hilbert review."""
import argparse
import hashlib
import json
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement
from math import comb
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]


def symmetric_character(n, weights):
    result = Counter()
    if n < 0:
        return result
    for mon in combinations_with_replacement(range(len(weights)), n):
        result[tuple(sum(weights[i][j] for i in mon) for j in range(2))] += 1
    return result


def core_character(n):
    # Rank-one matrices = Segre diagonal. Divide by the trace equation,
    # of GL2 weight (1,1), before taking the degree-n determinant twist.
    weights = [(2, 0), (1, 1), (0, 2)]
    result = Counter()
    for degree, sign in [(n, 1), (n - 1, -1)]:
        part = symmetric_character(degree, weights)
        for a, ca in part.items():
            for b, cb in part.items():
                result[n + a[0] - b[0], n + a[1] - b[1]] += sign * ca * cb
    assert min(result.values()) >= 0
    assert sum(result.values()) == (n + 1) ** 3
    return Counter({a: c for a, c in result.items() if c})


def tensor_character(left, right):
    result = Counter()
    for a, ca in left.items():
        for b, cb in right.items():
            result[a[0] + b[0], a[1] + b[1]] += ca * cb
    return result


def euler_character(character):
    # GL2 highest-weight multiplicities are successive weight differences.
    # The signed Sp4 Weyl dimension polynomial already gives the Euler
    # characteristic; no main-producer Bott sorting routine is used.
    result = Fraction(0)
    rows = []
    reconstructed = Counter()
    for (a, b), coefficient in sorted(character.items(), reverse=True):
        if a < b:
            continue
        multiplicity = coefficient - character[a + 1, b - 1]
        assert multiplicity >= 0
        if not multiplicity:
            continue
        value = Fraction((a - b + 1) * (b + 1) * (a + 2) * (a + b + 3), 6)
        assert value.denominator == 1
        result += multiplicity * value
        rows.append(dict(weight=[a, b], multiplicity=multiplicity, Euler=int(value)))
        for j in range(a - b + 1):
            reconstructed[a - j, b + j] += multiplicity
    assert reconstructed == character
    assert result.denominator == 1
    return int(result), rows


def verify():
    started = perf_counter()
    primary_path = 'reports/AB_RATIONAL_SINGULARITIES_HILBERT_20261102.json'
    primary = json.loads((ROOT / primary_path).read_text())
    assert primary['all_passed']
    for path, digest in primary['input_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest, path

    values = []
    details = []
    for n in range(9):
        total = 0
        pieces = []
        for vertical in range(n + 1):
            character = tensor_character(core_character(n - vertical),
                                         symmetric_character(vertical, [(3, 1), (2, 2), (1, 3)]))
            value, decomposition = euler_character(character)
            total += value
            pieces.append(dict(vertical_degree=vertical, Euler=value, decomposition=decomposition))
        values.append(total)
        details.append(pieces)
    assert values == primary['Hilbert_values']
    assert values[5] == 166818

    # Enumerate all local Tor shifts independently; the exact Y resolution
    # is a frozen input, not recomputed by this review.
    tor = Counter()
    for i, d, rank in [(0, 0, 1), (1, 2, 9), (2, 3, 16), (3, 4, 9), (4, 6, 1)]:
        for koszul in range(35):
            tor[i + koszul, d + koszul] += rank * comb(34, koszul)
    assert max(d for i, d in tor) == 40
    assert max(d - i for i, d in tor) == 2
    for p in range(39):
        allowed = [d for i, d in tor if 0 <= i - p <= 3]
        if allowed:
            assert max(allowed) == primary['global_resolution_shift_upper'][str(p)]
            assert max(allowed) <= p + 5

    # The numerator bound 40 - (45-10)=5 precedes this finite inversion.
    numerator = [sum((-1) ** j * comb(10, j) * values[n - j]
                     for j in range(n + 1)) for n in range(9)]
    assert numerator == [1, 35, 356, 986, 715, 91, 0, 0, 0]
    h = numerator[:6]
    assert h == primary['Hilbert_numerator']
    assert sum(h) == 2184
    for n in range(9):
        assert values[n] == sum(h[i] * comb(n - i + 9, 9) for i in range(min(n, 5) + 1))
        assert comb(n + 2, 2) ** 2 - comb(n + 1, 2) ** 2 == (n + 1) ** 3
    assert primary['regularity'] == 5
    assert primary['a_invariant'] == -5
    assert primary['projective_dimension'] == 35
    assert primary['ideal_generator_degree_upper'] == 6
    assert comb(46, 2) - values[2] == primary['quadratic_relations'] == 274

    paths = sorted(set(primary['input_sha256']) | {
        primary_path,
        str(Path(__file__).relative_to(ROOT)),
        'research/notes/AB_RATIONAL_SINGULARITIES_HILBERT_REVIEW_20261102.md',
        'research/notes/AB_ORIGIN_NORMALITY_REVIEW_20261102.md',
        'code/ab_origin_degree_four_cohomology_review_20261101.py',
    })
    return dict(id='AB_RATIONAL_SINGULARITIES_HILBERT_REVIEW_20261102', all_passed=True,
                independent_Cauchy_Euler_values=values, independent_Cauchy_details=details,
                independently_checked_degree_five=values[5], maximum_Tor_internal_shift=40,
                maximum_Tor_internal_minus_homological=2,
                Hilbert_numerator_degree_bound=5, independent_Hilbert_numerator=h,
                degree=2184, regularity=5, a_invariant=-5, projective_dimension=35,
                ideal_generator_degree_upper=6, geometric_proof_accepted=True,
                limitations=['The exact fiber resolution and global normality/higher-image vanishing are frozen dependencies.',
                             'The blowup and derived spectral-sequence arguments are proof audits, not program-only conclusions.',
                             'This concerns reduced AB, not the nonreduced T01 zero fiber or full S61 exhaustion.'],
                input_sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
                seconds=perf_counter() - started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print('AB RATIONAL/HILBERT REVIEW PASS', result['seconds'])

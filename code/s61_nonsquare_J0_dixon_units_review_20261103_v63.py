#!/usr/bin/env python3
"""Independent compact Fraction audit of the two frozen Dixon J0 units.

No solver, native pivot, or large normalization is rerun. Large matrices are
stream-hashed; only the two complete target rows are parsed as polynomials.
"""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from math import comb
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_NONSQUARE_J0_DIXON_UNITS_20261103_v63'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)


def polynomial_line(line):
    words = line.split()
    count = int(words[0])
    assert len(words) == 1 + 4 * count
    result = {}
    for i in range(count):
        coefficient = Fraction(words[1 + 4 * i])
        power = tuple(map(int, words[2 + 4 * i:5 + 4 * i]))
        assert coefficient and power not in result and power[2] == 0
        assert min(power) >= 0
        result[power] = coefficient
    return result


def decode(rows):
    result = {}
    for power, coefficient in rows:
        power = tuple(power)
        assert len(power) == 3 and power[2] == 0 and min(power) >= 0
        assert power not in result and Fraction(coefficient)
        result[power] = Fraction(coefficient)
    return result


def image(row, witness):
    result = {}
    for entry, multiplier in zip(row, witness):
        for monomial, coefficient in entry.items():
            for other, scalar in multiplier.items():
                key = tuple(x + y for x, y in zip(monomial, other))
                result[key] = result.get(key, Fraction(0)) + coefficient * scalar
    return {key: value for key, value in result.items() if value}


def bound_row(source, shape, target):
    # Scan every source line, retaining only the specified complete row.
    result = []
    with source.open() as f:
        assert list(map(int, f.readline().split())) == [*shape, 1]
        assert len(f.readline().split()) == 1
        labels = list(map(int, f.readline().split()))
        assert len(labels) == shape[0] and len(set(labels)) == shape[0]
        for i in range(shape[0]):
            for _ in range(shape[1]):
                line = f.readline()
                assert line
                if i == target:
                    result.append(polynomial_line(line))
        assert not f.read().strip()
    return result


def main(output):
    started = time.monotonic()
    folder = ROOT / 'reports' / STEM
    report = read(folder / 'report.json')
    payloads = read(folder / 'identities.json.gz')
    assert report['all_passed'] and len(payloads) == 2
    inputs = {}
    for name, expected in report['input_sha256'].items():
        path = ROOT / name
        actual = digest(path)
        assert actual == expected, name
        inputs[name] = actual
    previous = read(ROOT / 'reports/S61_NONSQUARE_J0_DEGREE1_UNIT_20261103_v63/report.json')
    assert previous['all_passed'] and previous['total_safe_pivots_from57basis'] == 37
    base = ROOT / 'work/nonsquare-J0-pid-v63'
    sources = [base / 'degree1-row4-v1/scalar-primitive/residual.txt',
               base / 'degree1-dixon-v1/residual.txt']
    rowsaved = [base / 'degree1-chain-v1/stage-1/row1.txt',
                base / 'dixon-chain-v1/stage-1/row0.txt']
    witnesses = [folder / 'unit38-row1.witness', folder / 'unit39-row0.witness']
    results = []
    for number, payload in enumerate(payloads):
        shape, target = ([20, 54], 1) if number == 0 else ([19, 53], 0)
        assert payload['source_shape'] == shape and payload['target_row'] == target
        assert digest(sources[number]) == payload['source_sha256']
        assert (number or payload['source_sha256'] == previous['normalized_residual_sha256'])
        row = list(map(decode, payload['row_polynomials']))
        multipliers = list(map(decode, payload['cofactor_polynomials']))
        assert len(row) == len(multipliers) == shape[1]
        assert bound_row(sources[number], shape, target) == row
        with rowsaved[number].open() as f:
            assert list(map(int, f.readline().split())) == [1, shape[1], 0]
            assert [polynomial_line(f.readline()) for _ in range(shape[1])] == row
            assert not f.read().strip()
        with witnesses[number].open() as f:
            assert int(f.readline()) == shape[1]
            assert [polynomial_line(f.readline()) for _ in range(shape[1])] == multipliers
            assert not f.read().strip()
        # Construct the unit by a binomial formula, independently of factor division.
        apower, epower = (0, 3) if number == 0 else (1, 6)
        expected = {(apower, k, 0): Fraction(comb(epower, k) * 3**k * 2**(epower-k))
                    for k in range(epower + 1)}
        assert image(row, multipliers) == expected == decode(payload['target_polynomial'])
        replacement = 7 if number == 0 else 4
        assert set(multipliers[replacement]) == {(0, 0, 0)}
        determinant = multipliers[replacement][(0, 0, 0)]
        assert determinant
        log = (folder / f'unit{number+38}-pivot.log').read_text()
        assert f'EXACT_WITNESS_REPLACE {replacement} TARGET {target}' in log and 'PASS 1' in log
        active = sum(bool(q) for q in multipliers)
        terms = sum(map(len, multipliers))
        assert (active, terms) == ((26, 34) if number == 0 else (44, 84))
        results.append(dict(unit=number+38, source_shape=shape, target_row=target,
            full_source_stream_hash_and_complete_target_row_binding=True,
            archived_full_row_and_witness_match=True,
            independent_Fraction_full_polynomial_identity=True,
            target='(3e+2)^3' if number == 0 else 'a*(3e+2)^6',
            active_cofactors=active, cofactor_terms=terms,
            invertible_column_replacement_index=replacement,
            replacement_determinant_nonzero_QQ=str(determinant),
            output_shape=[shape[0]-1, shape[1]-1]))
    assert (folder / 'units.txt').read_text().splitlines() == [
        'a', '500*a+90*e^2+3*e-24', '300*a*e+54*e^3+9*e^2-24*e-4',
        '100*a+18*e^2-9*e-2', '3*e+2']
    # Confirm that the archived QQ normalization is only a nonzero rational
    # diagonal change. Its full 447025-term replay is the pinned producer receipt.
    with (folder / 'unit38-QQ-factors.txt').open() as f:
        assert list(map(int, f.readline().split())) == [19, 53]
        for kind, count in [('ROW', 19), ('COLUMN', 53)]:
            for i in range(count):
                tag, index, value = f.readline().split()
                assert (tag, int(index)) == (kind, i) and Fraction(value)
        assert not f.read().strip()
    assert report['unit38_all_QQ_diagonal_terms_independently_verified'] == 447025
    assert report['unit39_optional_normalization_failed_and_excluded'] is True
    residual = ROOT / report['last_confirmed_state']
    assert residual.stat().st_size == report['new_residual_bytes'] == 1637754878
    assert digest(residual) == report['new_residual_sha256']
    for path in [Path(__file__), ROOT/'research/notes'/f'{STEM}_REVIEW.md',
                 folder/'report.json', folder/'identities.json.gz',
                 ROOT/'research/notes'/f'{STEM}.md', folder/'units.txt',
                 folder/'unit38-QQ-factors.txt', *witnesses,
                 folder/'unit38-pivot.log', folder/'unit39-pivot.log', residual]:
        inputs[str(path.relative_to(ROOT))] = digest(path)
    result = dict(id=STEM+'_REVIEW', all_passed=True,
        base_ring='QQ[a,e,1/(a*F*Bs*C*(3e+2))]',
        independent_compact_identity_and_source_binding_results=results,
        both_constant_cofactor_column_changes_unimodular_over_QQ_polynomials=True,
        both_pivots_units_using_only_previously_declared_open_factors=True,
        J_not_inverted=True, unknown_factor_not_inverted=True,
        previous37_safe_pivots_reused=True, safe_pivots_now39=39,
        cokernel_preserved_by_invertible_operations=True, residual_shape=[18,52],
        final_large_residual_stream_hash_checked=True,
        large_native_pivot_and_QQ_normalization_not_rerun=True,
        pinned_producer_full447025_term_QQ_diagonal_check_reused=True,
        no_finite_defect_full_coker_or_quartic_coverage_conclusion=True,
        input_sha256=inputs, seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','safe_pivots_now39','residual_shape','seconds']},indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)

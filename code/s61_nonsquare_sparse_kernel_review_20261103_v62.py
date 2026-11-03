#!/usr/bin/env python3
"""Independent bounded review of frozen V4 sparse kernels and row identities."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import tarfile
import time

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

ROOT = Path(__file__).resolve().parents[1]
P = (1 << 61) - 1
R, aa, ee, dd = ring('a,e,d', QQ)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def parse_poly(line):
    tokens = line.split()
    count = int(tokens[0])
    assert len(tokens) == 1 + 4 * count
    answer = {}
    for i in range(count):
        q = Fraction(tokens[1 + 4 * i])
        exponent = tuple(map(int, tokens[2 + 4 * i:5 + 4 * i]))
        assert len(exponent) == 3 and min(exponent) >= 0
        answer[exponent] = answer.get(exponent, Fraction()) + q
    return {exponent: q for exponent, q in answer.items() if q}


def parse_matrix(data):
    lines = data.decode().splitlines()
    n, m, target = map(int, lines[0].split())
    assert len(lines) == 1 + n * m
    polys = [parse_poly(line) for line in lines[1:]]
    return [polys[i * m:(i + 1) * m] for i in range(n)], target


def parse_kernel(data):
    lines = data.decode().splitlines()
    n, m, _ = map(int, lines[0].split())
    assert len(lines) == 3 + n * m
    polys = [parse_poly(line) for line in lines[3:]]
    return [polys[i * m:(i + 1) * m] for i in range(n)]


def multiply(x, y):
    answer = {}
    for u, c in x.items():
        for v, d in y.items():
            monomial = tuple(u[i] + v[i] for i in range(3))
            answer[monomial] = answer.get(monomial, Fraction()) + c * d
    return {u: c for u, c in answer.items() if c}


def sum_products(coefficients, matrix, column):
    answer = {}
    for i, polynomial in coefficients.items():
        for monomial, q in multiply(polynomial, matrix[i][column]).items():
            answer[monomial] = answer.get(monomial, Fraction()) + q
    return {u: c for u, c in answer.items() if c}


def to_ring(polynomial):
    return R.from_dict({u: QQ(q.numerator, q.denominator) for u, q in polynomial.items()})


def from_ring(polynomial):
    return {u: Fraction(str(q)) for u, q in polynomial.items()}


def known_unit(polynomial, units):
    value = to_ring(polynomial)
    assert value
    for unit in units:
        assert unit and set(unit) != {(0, 0, 0)}
        while True:
            quotient, remainder = value.div(unit)
            if remainder:
                break
            value = quotient
    assert set(value) == {(0, 0, 0)} and value[(0, 0, 0)]


def monomials(degree):
    return [(a, e, total-a-e) for total in range(degree+1)
            for a in range(total+1) for e in range(total-a+1)]


def check_forced_zero(matrix, degree, sequence):
    mus = monomials(degree)
    ids, keys = {}, []
    for j in range(len(matrix[0])):
        for row in matrix:
            for exponent in sorted(row[j], reverse=True):
                for mu in mus:
                    key = (j, *(exponent[k] + mu[k] for k in range(3)))
                    if key not in ids:
                        ids[key] = len(keys)
                        keys.append(key)
    zeros = set()
    for unknown, equation in sequence:
        assert unknown not in zeros and 0 <= equation < len(keys)
        column, *exponent = keys[equation]
        live = []
        for i, row in enumerate(matrix):
            for j, mu in enumerate(mus):
                u = i * len(mus) + j
                if u in zeros:
                    continue
                difference = tuple(exponent[k] - mu[k] for k in range(3))
                if min(difference) >= 0 and row[column].get(difference):
                    live.append(u)
        assert live == [unknown]
        zeros.add(unknown)
    return len(keys), len(zeros)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    handoff = ROOT / 'research/handoff/20261103_v62'
    frozen = handoff / 'v4_sparse_kernel'
    manifest = json.loads((frozen / 'manifest.json').read_text())
    archive = frozen / manifest['archive']
    assert digest(archive.read_bytes()) == manifest['archive_sha256']
    with tarfile.open(archive) as tar:
        assert len(tar.getmembers()) == len(manifest['files']) == 101
        contents = {}
        for member in tar.getmembers():
            assert member.isfile() and member.name in manifest['files']
            data = tar.extractfile(member).read()
            expected = manifest['files'][member.name]
            assert len(data) == expected['bytes'] and digest(data) == expected['sha256']
            contents[member.name] = data
    v3 = json.loads((handoff / 'nonsquare_core_wip_v3_manifest.json').read_text())
    v3archive = handoff / v3['archive']
    assert digest(v3archive.read_bytes()) == v3['archive_sha256']

    def get(name):
        if name in contents:
            return contents[name]
        expected = v3['files'][name]
        path = ROOT / name
        if path.exists() and digest(path.read_bytes()) == expected['sha256']:
            return path.read_bytes()
        with tarfile.open(v3archive) as tar:
            data = tar.extractfile(name).read()
        assert len(data) == expected['bytes'] and digest(data) == expected['sha256']
        return data

    prefix = 'work/nonsquare-v62/kernel-v4/'
    public_dir = ROOT / 'reports/S61_BINARY_NONSQUARE_CUBIC_REDUNDANCIES_20261103_v62'
    public = json.loads((public_dir / 'report.json').read_text())
    assert public['all_passed'] and not public['all_point_cubic_rank57_claimed']
    for name, value in public['input_sha256'].items():
        assert digest((ROOT / name).read_bytes()) == value
    for name, value in public['generated_sha256'].items():
        assert digest((public_dir / name).read_bytes()) == value
    certificate = json.loads((ROOT / 'code/data/s61_binary_nonsquare_cubic_redundancies_20261103_v62.json').read_text())
    full_bytes = get('work/nonsquare-v62/fresh-cubic/full')
    assert digest(full_bytes) == certificate['original_full_matrix_sha256']
    full, target = parse_matrix(full_bytes)
    normalized_bytes = get('work/nonsquare-v62/fresh-cubic/syzygies-normalized-degree1.normalized-input')
    A, _ = parse_matrix(normalized_bytes)
    assert len(A) == 63 and len(A[0]) == 91 and target == 57
    normalization_lines = (public_dir / 'normalization.txt').read_text().splitlines()
    assert normalization_lines[0] == '63' and len(normalization_lines) == 64
    multipliers = [parse_poly(line) for line in normalization_lines[1:]]
    units = [R.from_expr(sp.sympify(value, locals={'a': sp.Symbol('a'), 'e': sp.Symbol('e'), 'd': sp.Symbol('d')}))
             for value in certificate['known_units']]
    for i in range(63):
        known_unit(multipliers[i], units)
        for j in range(91):
            assert multiply(A[i][j], multipliers[i]) == full[i][j]
    exact_identities = []
    deleted = []
    expected_coefficients = [{(0, 0, 0): Fraction(1)}, {(1, 0, 0): Fraction(12, 5)}, {(0, 0, 0): Fraction(1)}]
    for r, relation in enumerate(certificate['relations']):
        coefficients = {item['original_row']: parse_poly(str(len(item['terms'])) + ''.join(' ' + ' '.join(map(str, term)) for term in item['terms']))
                        for item in relation['coefficients']}
        drop, degree = relation['deleted_original_row'], relation['coefficient_degree_bound']
        assert drop == [48, 14, 55][r] and coefficients[drop] == expected_coefficients[r]
        assert all(not coefficients.get(i) for i in deleted)
        assert all(sum(exponent) <= degree for polynomial in coefficients.values() for exponent in polynomial)
        known_unit(coefficients[drop], units)
        for column in range(91):
            assert not sum_products(coefficients, A, column)
        exact_identities.append(dict(deleted_original_row=drop, coefficient_degree=degree,
                                     all_QQ_columns_checked=91, unit_coefficient=str(to_ring(coefficients[drop]))))
        deleted.append(drop)
    remaining_bytes = gzip.decompress((public_dir / 'remaining.txt.gz').read_bytes())
    remaining, _ = parse_matrix(remaining_bytes)
    assert remaining == [row for i, row in enumerate(A) if i not in deleted]
    assert digest(remaining_bytes) == manifest['omitted_regenerable_large_inputs'][prefix+'full-minus-14-48-55']['sha256']
    twice_removed = [row for i, row in enumerate(A) if i not in [48, 14]]
    assert twice_removed == parse_matrix(get('work/nonsquare-v62/fresh-cubic/full-minus-14-48'))[0]
    assert all(q.denominator % 32003 and q.denominator % P for row in A for polynomial in row for q in polynomial.values())

    forced = []
    for degree, count, active in [(4, 90, 2045), (5, 118, 3298)]:
        stem = f'certified-degree{degree}'
        sequence = [tuple(map(int, line.split())) for line in get(prefix+stem+'.forced-zero').decode().splitlines()]
        equations, zero_count = check_forced_zero(twice_removed, degree, sequence)
        assert zero_count == count
        assert 61 * math.comb(degree+3, 3) - count == active
        selected = list(map(int, get(prefix+stem+'.selected-equations').split()))
        assert len(selected) == len(set(selected)) == active and max(selected) < equations
        log = get(prefix+stem+'.log').decode()
        assert f'MODP_COMPLETE_EQUATION_RANK {active} NULLITY 0' in log
        assert not parse_kernel(get(prefix+stem))
        forced.append(dict(degree=degree, exact_singleton_steps=count,
                           coefficient_equations=equations, remaining_rank_from_audited_native=active))

    # Regenerate the chart exactly, preserving the input term order for its hash.
    chart = json.loads(get(prefix+'P-chart.json'))
    shifts = chart['row_multiplier_powers_of_a']
    chart_rows, lines = [], ['60 91 57']
    for i, row in enumerate(remaining):
        converted = []
        for polynomial in row:
            terms = [((exponent[0]-exponent[2]+shifts[i], exponent[1], exponent[2]), q)
                     for exponent, q in polynomial.items()]
            assert all(min(exponent) >= 0 for exponent, q in terms)
            converted.append(dict(terms))
            lines.append(str(len(terms)) + ''.join(' '+str(q)+' '+' '.join(map(str, exponent)) for exponent, q in terms))
        chart_rows.append(converted)
    chart_bytes = ('\n'.join(lines)+'\n').encode()
    assert digest(chart_bytes) == chart['output_sha256']
    assert chart['output_sha256'] == manifest['omitted_regenerable_large_inputs'][prefix+'P-chart-full']['sha256']
    # Check that precisely the old known units are being used in this new chart.
    p_units = [R.from_expr(sp.sympify(value, locals={'a':sp.Symbol('a'),'e':sp.Symbol('e'),'d':sp.Symbol('d')}))
               for value in get(prefix+'P-chart-units').decode().splitlines()]
    assert p_units == [aa, 500*aa+90*ee**2+3*ee-24,
                       (50*aa+18*ee**2-9*ee-2)*dd+(60*ee-40)*aa,
                       100*aa+18*ee**2-9*ee-2, 5*dd-6*ee-4, dd]

    work = ROOT / 'work/nonsquare-v62/connectivity-v4-review'
    work.mkdir(parents=True, exist_ok=True)
    flags = subprocess.check_output(['pkg-config','--cflags','--libs','flint'], text=True).split()
    binaries = {}
    for stem in ['sparse-early-certified-kernel', 'certified-modular-low-syzygies']:
        source = work / (stem+'.cpp')
        source.write_bytes(get(prefix+stem+'.cpp'))
        executable = work / stem
        result = subprocess.run(['c++','-O2','-std=c++17',str(source),*flags,'-o',str(executable)],
                                capture_output=True,text=True,check=True,timeout=60)
        (args.output / (stem+'-compile.log')).write_text(result.stdout+result.stderr)
        binaries[stem] = executable

    def native_run(stem, data, degree, engine='sparse-early-certified-kernel', extra=None, expected=0):
        inp, out = work/(stem+'.in'), work/stem
        inp.write_bytes(data)
        command = [str(binaries[engine]),str(inp),str(degree),str(out)]
        if extra is not None:
            command.append(str(extra))
        result = subprocess.run(command,capture_output=True,text=True,timeout=60)
        assert result.returncode == expected
        (args.output/(stem+'.log')).write_text(result.stdout+result.stderr)
        if expected:
            assert 'FINAL_QQ_KERNEL_RECONSTRUCTION_NOT_CERTIFIED' in result.stderr
            return None, out
        assert 'PASS' in result.stderr
        return parse_kernel(out.read_bytes()), out

    # Multi-dimensional kernel: independent QQ nullity and free-vector independence.
    fixture = b'3 1 0\n1 1 1 0 0\n1 1 0 1 0\n2 1 1 0 0 1 0 1 0\n'
    vectors, _ = native_run('five-dimensional-fixture', fixture, 1)
    assert len(vectors) == 5
    mus = monomials(1)
    flat = sp.Matrix([[sp.Rational(polynomial.get(mu,Fraction())) for polynomial in vector for mu in mus]
                      for vector in vectors])
    assert flat.rank() == 5
    fixture_matrix, _ = parse_matrix(fixture)
    for vector in vectors:
        assert not sum_products(dict(enumerate(vector)),fixture_matrix,0)
    # A bad prime cannot produce an accepted QQ kernel.
    bad = f'1 1 0\n1 {P} 1 0 0\n'.encode()
    native_run('bad-prime-refusal',bad,0,expected=4)
    forced_fixture = b'2 2 0\n1 1 1 0 0\n1 1 0 0 0\n1 1 0 1 0\n0\n'
    v, out = native_run('forced-zero-fixture',forced_fixture,1,'certified-modular-low-syzygies')
    assert not v
    seq = [tuple(map(int,line.split())) for line in out.with_name(out.name+'.forced-zero').read_text().splitlines()]
    assert check_forced_zero(parse_matrix(forced_fixture)[0],1,seq)[1] == 8
    # Fresh bounded real-input replay, not the expensive degree6/7 discoveries.
    v, out = native_run('original-degree2-fresh',normalized_bytes,2)
    assert len(v) == 1 and v[0][48] == {(0,0,0):Fraction(1)}
    assert out.read_bytes() == get(prefix+'early-degree2')

    units_path = work / 'P-units'
    units_path.write_bytes(get(prefix+'P-chart-units'))
    _, out = native_run('P-normalization-degree0-fresh',chart_bytes,0,extra=units_path)
    normalized_p_bytes = out.with_name(out.name+'.normalized-input').read_bytes()
    assert digest(normalized_p_bytes) == manifest['omitted_regenerable_large_inputs'][prefix+'P-degree3.normalized-input']['sha256']
    p_rows, _ = parse_matrix(normalized_p_bytes)
    for raw, normalized in zip(chart_rows, p_rows):
        first = next(j for j in range(91) if normalized[j])
        multiplier = to_ring(raw[first]).exquo(to_ring(normalized[first]))
        polynomial = from_ring(multiplier)
        known_unit(polynomial,p_units)
        assert all(multiply(normalized[j],polynomial)==raw[j] for j in range(91))

    # Accept saved large finite ranks only after auditing their exact algorithm;
    # no repeated 80-90 second computations are hidden in this reviewer.
    status = json.loads((frozen/'STATUS.json').read_text())
    saved_negative = []
    for item in status['negative_results']:
        if item['stem'].startswith('certified-'):
            continue
        stem = item['stem']
        log = get(prefix+stem+'.log').decode()
        assert f'SPARSE_COMPLETE_RANK {item["rank"]} NULLITY 0' in log
        assert not parse_kernel(get(prefix+stem))
        selected = list(map(int,get(prefix+stem+'.selected-equations').split()))
        assert len(selected)==len(set(selected))==item['rank']
        saved_negative.append(item)
    receipts = json.loads(get(prefix+'execution-receipts.json'))
    assert all(item['status']=='fulfilled' and item['value']['receipt']['exit_code']==0 for item in receipts)
    assert sp.isprime(P)
    examples = [0,1,2,P//2,P-2,P-1]
    for x in examples:
        for y in examples:
            product=x*y
            reduced=(product&P)+(product>>61)
            if reduced>=P:
                reduced-=P
            assert reduced==product%P
    inputs = [Path(__file__),frozen/'manifest.json',frozen/'STATUS.json',frozen/'README.md',archive,
              handoff/'nonsquare_core_wip_v3_manifest.json',v3archive,public_dir/'report.json',
              ROOT/'research/notes/S61_BINARY_NONSQUARE_SPARSE_KERNEL_20261103_v62_REVIEW.md']
    report = dict(id='S61_BINARY_NONSQUARE_SPARSE_KERNEL_20261103_v62_REVIEW',all_passed=True,
                  frozen_v4_members_verified=101,original_full_matrix_shape=[63,91],
                  pure_Fraction_all_column_QQ_identities=exact_identities,
                  exact_known_unit_row_reconstruction=True,exact_forced_zero_sequences=forced,
                  rational_reconstruction_rank_sandwich_accepted=True,
                  early_stop_full_identity_rank_sandwich_accepted=True,
                  fresh_multidimensional_kernel_and_unlucky_prime_tests_passed=True,
                  fresh_original_degree2_replay_passed=True,
                  P_chart_all_entry_transform_and_unit_normalization_verified=True,
                  saved_negative_finite_degree_results=saved_negative,
                  expensive_degree6_degree7_ranks_independently_rerun=False,
                  all_point_rank57_or_product41_proved_by_this_packet=False,
                  full_quartic_or_ambient_removal_proved=False,
                  scope='Three localized row redundancies and exact finite coefficient-degree kernels only.',
                  input_sha256={str(p.relative_to(ROOT)):digest(p.read_bytes()) for p in inputs},
                  evidence_sha256={p.name:digest(p.read_bytes()) for p in args.output.iterdir() if p.is_file()},
                  seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS V4 independent sparse-kernel review; row48/14/55; exact finite scopes;',report['seconds'])


if __name__=='__main__':
    main()

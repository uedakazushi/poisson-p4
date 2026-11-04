#!/usr/bin/env python3
"""Independent Python integer coefficient words from actual source5 to seven invariants."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import gzip
import hashlib
import json
import math
from pathlib import Path
import resource
import signal
import time

import sympy as s
from sympy.polys.rings import ring

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63_REVIEW'
PWR = [1 << (8*i) for i in range(4)]
MAX_TERMS = 0
MULS = 0


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for a in iter(lambda: f.read(1 << 20), b''):
            h.update(a)
    return h.hexdigest()


def guard(a):
    global MAX_TERMS
    MAX_TERMS = max(MAX_TERMS, len(a))
    assert len(a) <= 150000
    return a


def add(a, b, c=1):
    r = dict(a)
    for m, v in b.items():
        r[m] = r.get(m, 0)+c*v
        if not r[m]:
            del r[m]
    return guard(r)


def mul(a, b):
    global MULS
    MULS += 1
    r = defaultdict(int)
    for m, x in a.items():
        for n, y in b.items():
            r[m+n] += x*y
    return guard({m: x for m, x in r.items() if x})


def total(values):
    z = {}
    for a in values:
        z = add(z, a)
    return z


def scale(a, c):
    assert int(c) == c
    return {m: int(c)*x for m, x in a.items() if c*x}


def mm(A, B):
    return [[total(mul(A[i][k], B[k][j]) for k in range(4)) for j in range(4)] for i in range(4)]


def tr(A):
    return total(A[i][i] for i in range(4))


def pair(A, B):
    return total(mul(A[i][j], B[j][i]) for i in range(4) for j in range(4))


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    out.mkdir(parents=True, exist_ok=False)
    folder = ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63'
    ar = read(folder/'report.json')
    assert ar['all_passed']
    source_path = ROOT/'reports/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63/source-five-circuit.json.gz'
    spin_path = ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    gate = ROOT/'reports/Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW/report.json'
    spin_gate = ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63_REVIEW/report.json'
    assert read(gate)['all_passed'] and read(spin_gate)['all_passed']
    data, spin = map(read, [source_path, spin_path])
    assert len(data['nodes']) == 1032
    # Source5 has already been proved dominant for source constants; derive the
    # sixteen numerators from its actual arithmetic circuit, not the native input.
    R, a, X, W, B = ring('a,X,W,B', s.QQ)
    K = R.to_field()
    replacement = [K(a*X), K.one, K(X), K(W), K(B)]
    cache = {}
    def value(i):
        if i in cache:
            return cache[i]
        op, arg = data['nodes'][i]
        if op == 'q':
            z = K(s.Rational(arg))
        elif op == 'v':
            z = replacement[arg]
        elif op == 'inv':
            z = 1/value(arg)
        elif op == '+':
            z = sum((value(j) for j in arg), K.zero)
        else:
            assert op == '*'
            z = K.one
            for j in arg:
                z *= value(j)
        cache[i] = z
        return z
    D = 147456*X**2-(3072*W-1)**2
    clear = K(1080*B*D)
    bp = []
    for i in data['outputs']['beta16']:
        p = value(i)*clear
        assert p.denom == 1 and all(s.denom(s.Rational(x)) == 1 for x in p.numer.values())
        assert all(all(e <= 4 for e in m) for m in p.numer)
        bp.append({sum(e*q for e, q in zip(m, PWR)): int(c) for m, c in p.numer.items()})
    # All words below have beta degree<=12, so each exponent<=48. Packed
    # exponent base256 cannot carry between coordinates at any multiplication.
    with gzip.open(folder/'native-input.txt.gz', 'rt') as f:
        tokens = iter(f.read().split())
    assert [int(next(tokens)) for _ in range(4)] == [4, 16, 136, 1080]
    for expected in bp:
        n = int(next(tokens))
        own = {}
        for _ in range(n):
            c, eu, ex, ew, eb = [int(next(tokens)) for _ in range(5)]
            key = eu*PWR[0]+(eu+ex)*PWR[1]+ew*PWR[2]+eb*PWR[3]
            own[key] = own.get(key, 0)+c
        assert {m: c for m, c in own.items() if c} == expected
    monomials = [[int(next(tokens)), int(next(tokens))] for _ in range(136)]
    assert monomials == spin['quadratic_monomials']
    coefficient_rows = [[int(next(tokens)) for _ in range(136)] for _ in range(32)]
    assert next(tokens, None) is None
    rows = spin['spin_adjoint_C_coefficients16x136']+spin['vector_adjoint_D_coefficients16x136']
    assert coefficient_rows == [[int(4*F(a)) for a in row] for row in rows]
    assert all((4*F(a)).denominator == 1 for row in rows for a in row)
    products = [mul(bp[i], bp[j]) for i, j in monomials]
    cov = [total(scale(products[i], c) for i, c in enumerate(row) if c) for row in coefficient_rows]
    C = [cov[4*i:4*i+4] for i in range(4)]
    E = [cov[16+4*i:16+4*i+4] for i in range(4)]
    C2, CE, E2 = mm(C, C), mm(C, E), mm(E, E)
    qC, qE, qCE = tr(C2), tr(E2), tr(CE)
    assert add(qC, qE, -2) == {} and add(qCE, qE) == {}
    fword = add(add(qC, qCE, -4), qE, 4)
    # Different determinant engine: Newton for traceless4x4 matrices. The
    # author used a direct sum over all24 permutations.
    EC = mm(E, C)
    cross = [[add(CE[i][j], EC[i][j]) for j in range(4)] for i in range(4)]
    mats = [C2, cross, E2]
    traces = [tr(A) for A in mats]
    pencil = []
    for d in range(5):
        aa, bb = {}, {}
        for i in range(3):
            j = d-i
            if 0 <= j < 3:
                aa = add(aa, mul(traces[i], traces[j]))
                bb = add(bb, pair(mats[i], mats[j]))
        pencil.append(add(aa, bb, -2))
    # Compute the degree12 trace in the opposite cyclic order E²C E²C.
    H = mm(E2, C)
    nword = add(add(scale(pair(H, H), 8), mul(qE, pair(C2, E2)), -4),
        add(mul(mul(qE, qE), qC), mul(pair(E2, E2), qC), -2))
    words = [fword, *pencil, nword]
    names = ['F4', 'p0', 'p1', 'p2', 'p3', 'p4', 'N12']
    degrees = [4, 8, 8, 8, 8, 8, 12]
    scales = [16, *([2048]*5), 32768]
    scalar = [2, *([4]*5), 6]
    counts = []
    for name, raw, d, z, cw in zip(names, words, degrees, scales, scalar):
        saved = read(folder/(name+'.json.gz'))
        assert saved['name'] == name and saved['source_c_power'] == cw
        assert saved['denominator_D_power'] == d
        den, bpow = int(saved['denominator_integer']), saved['denominator_B_power']
        assert den > 0 and 0 <= bpow <= d
        numerator = {}
        for coefficient, ea, ey, ew, eb in saved['terms']:
            key = ea*PWR[0]+2*ey*PWR[1]+ew*PWR[2]+eb*PWR[3]
            assert key not in numerator and coefficient != '0'
            numerator[key] = int(coefficient)
        # D^degree is deliberately retained on BOTH sides. Exact coefficient
        # equality after B/integer clearing proves every serialized coefficient.
        lhs = {m+bpow*PWR[3]: den*c for m, c in raw.items()}
        rhs = {m+d*PWR[3]: z*1080**d*c for m, c in numerator.items()}
        assert lhs == rhs
        assert all(((m//PWR[1]) % 256) % 2 == 0 for m in raw)
        assert math.gcd(den, *numerator.values()) == 1
        assert min((m//PWR[3]) % 256 for m in numerator) == 0
        counts.append(len(numerator))
    assert counts == [96, 336, 1040, 1061, 973, 893, 1700]
    # p0,p4 here are Newton numerators8det(C+tD). Clear that known8.
    assert add(add(scale(mul(fword, fword), 8), pencil[0], -400), pencil[4], -1600) == {}
    certificate = dict(source_input_names=data['input_names'], source_nodes=1032,
        beta_numerator_terms=[len(p) for p in bp], independent_output_term_counts=counts,
        polynomial_engine='Python arbitrary precision integer sparse packed-exponent dictionaries, no native worker',
        determinant_engine='Newton trace identity; author permutation determinant not reused',
        N12_trace_order='tr(E^2 C E^2 C), opposite cyclic order to author',
        shape_map='u=aX, Y=X^2; all invariant numerators even in X',
        word_denominator_scales=dict(F4='16*(1080BD)^4', pencil_Newton='2048*(1080BD)^8', N12='32768*(1080BD)^12'),
        max_intermediate_terms=MAX_TERMS, independent_polynomial_multiplications=MULS,
        source_c_powers=scalar, all7_complete_polynomials_verified=True)
    (out/'independent-certificate.json').write_text(json.dumps(certificate, indent=2)+'\n')
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', source_path, spin_path, gate, spin_gate,
        folder/'report.json', folder/'native-report.json', folder/'native-input.txt.gz',
        *[folder/(name+'.json.gz') for name in names],
        ROOT/'code/q_new_family_low_invariant_expand_20261104_v63.py',
        ROOT/'code/q_new_family_low_invariant_expand_20261104_v63.cpp']
    result = dict(id=STEM, all_passed=True,
        actual_reviewed1032_source_circuit_to_all16_cleared_beta_polynomials=True,
        all_native_input_coefficients_and32_covariant_rows_bound=True,
        independent_integer_word_arithmetic_all7_complete_coefficient_maps=True,
        exact_Newton_pencil_and_cyclic_N12_definition_match=True,
        shape_parity_integer_content_B_denominators_and_restored_scalar_weights=True,
        output_term_counts=counts,
        no_author_native_replay_new_point_or_degree24_relation_search=True,
        no_all_Inv24_exclusion_or_NEW_equation_claim=True,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'independent-certificate.json': sha(out/'independent-certificate.json')},
        seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'output_term_counts', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())

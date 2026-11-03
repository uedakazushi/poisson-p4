#!/usr/bin/env python3
"""Independent direct Laurent tensor audit of the complete CA self-flag 4jet."""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_CA_ORIGIN_SELF_FLAG_4JET_20261103_v63_REVIEW'
Z = (0, 0, 0, 0, 0)  # epsilon, t, y0, y1, y2
PAIRS = list(it.combinations(range(5), 2))
MONS = list(it.combinations_with_replacement(range(5), 2))
PI = {p: i for i, p in enumerate(PAIRS)}
MI = {m: i for i, m in enumerate(MONS)}


class Ring(dict):
    """QQ[t,t^-1,y0,y1,y2][epsilon]/epsilon^5, independently implemented."""
    def __init__(self, x=0):
        if isinstance(x, dict):
            super().__init__((tuple(k), Q(v)) for k, v in x.items() if v and k[0] <= 4)
        else:
            super().__init__({Z: Q(x)} if x else {})

    def __add__(self, x):
        x = Ring(x)
        out = dict(self)
        for m, c in x.items():
            out[m] = out.get(m, Q(0)) + c
            if not out[m]:
                del out[m]
        return Ring(out)

    __radd__ = __add__

    def __neg__(self):
        return Ring({m: -c for m, c in self.items()})

    def __sub__(self, x):
        return self + -Ring(x)

    def __rsub__(self, x):
        return Ring(x) + -self

    def __mul__(self, x):
        x = Ring(x)
        out = {}
        for m, c in self.items():
            for n, d in x.items():
                if m[0] + n[0] > 4:
                    continue
                k = tuple(a + b for a, b in zip(m, n))
                out[k] = out.get(k, Q(0)) + c * d
        return Ring(out)

    __rmul__ = __mul__

    def __truediv__(self, n):
        return self * (Q(1) / Q(n))

    def __pow__(self, n):
        assert n >= 0
        result = Ring(1)
        for _ in range(n):
            result = result * self
        return result

    def shift(self, n):
        return Ring({(m[0] + n, *m[1:]): c for m, c in self.items() if m[0] + n <= 4})

    def coefficient(self, n):
        return Ring({(0, *m[1:]): c for m, c in self.items() if m[0] == n})

    def evaluate(self, values):
        out = {}
        for m, c in self.items():
            k = (m[0], 0, 0, 0, 0)
            out[k] = out.get(k, Q(0)) + c * Q(values[0]) ** m[1] * Q(values[1]) ** m[2] * Q(values[2]) ** m[3] * Q(values[3]) ** m[4]
        return Ring(out)

    def save(self):
        return [[list(m), str(c)] for m, c in sorted(self.items())]


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def matrix(x):
    return s.Matrix([[s.Rational(v) for v in row] for row in x])


def fraction(x):
    return Q(int(s.numer(x)), int(s.denom(x)))


def constmul(A, v):
    return [sum((fraction(A[i, j]) * v[j] for j in range(A.cols) if A[i, j]), Ring()) for i in range(A.rows)]


def addv(a, b):
    return [x + y for x, y in zip(a, b)]


def negv(a):
    return [-x for x in a]


def matmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(5)), Ring()) for j in range(5)] for i in range(5)]


def matadd(A, B):
    return [[A[i][j] + B[i][j] for j in range(5)] for i in range(5)]


def direct_push(g, source):
    """Literal inverse substitution followed by exterior g, through epsilon^4."""
    identity = [[Ring(i == j) for j in range(5)] for i in range(5)]
    nil = [[g[i][j] - identity[i][j] for j in range(5)] for i in range(5)]
    inverse = identity
    power = identity
    for n in range(1, 5):
        power = matmul(power, nil)
        inverse = matadd(inverse, [[((-1) ** n) * x for x in row] for row in power])
    assert matmul(g, inverse) == identity == matmul(inverse, g)
    wedge = [[g[i][k] * g[j][l] - g[i][l] * g[j][k] for k, l in PAIRS] for i, j in PAIRS]
    substitute = []
    for k, l in MONS:
        substitute.append([inverse[k][i] * inverse[l][j] + (inverse[k][j] * inverse[l][i] if i != j else Ring()) for i, j in MONS])
    # First substitute coordinates in all ten coefficients; then apply exterior g.
    middle = [[sum((source[15 * k + m] * substitute[m][n] for m in range(15) if source[15 * k + m] and substitute[m][n]), Ring()) for n in range(15)] for k in range(10)]
    return [sum((wedge[i][k] * middle[k][n] for k in range(10) if wedge[i][k] and middle[k][n]), Ring()) for i in range(10) for n in range(15)]


def first_derivative(A, v):
    """Differentiate the literal pushforward, without producer's Lie operators."""
    out = [Q(0)] * 150
    for p, (k, l) in enumerate(PAIRS):
        for m, (r, t) in enumerate(MONS):
            c = fraction(v[15 * p + m])
            if not c:
                continue
            for a in range(5):
                if a != l and A[a, k]:
                    out[15 * PI[tuple(sorted((a, l)))] + m] += (1 if a < l else -1) * fraction(A[a, k]) * c
                if k != a and A[a, l]:
                    out[15 * PI[tuple(sorted((k, a)))] + m] += (1 if k < a else -1) * fraction(A[a, l]) * c
            for old, other in ((r, t), (t, r)):
                for a in range(5):
                    if A[old, a]:
                        out[15 * p + MI[tuple(sorted((a, other)))]] -= fraction(A[old, a]) * c
    return s.Matrix([s.Rational(c.numerator, c.denominator) for c in out])


def reconstruct_CA():
    h, x, y, z, w = X = s.symbols('h x y z w')
    cubic = [s.prod((x, h, z, w)[i] for i in inds) for inds in it.combinations_with_replacement(range(4), 3) if any(i >= 2 for i in inds)]
    assert len(cubic) == 16
    columns = []
    for index in range(20):
        f = (x * x, h * x, h * h)[index] if index < 3 else 0
        g = int(index == 3)
        H = cubic[index - 4] if index >= 4 else s.Integer(0)
        P = s.zeros(5)
        P[0, 2] = f
        P[3, 4] = g * x * x
        P[2, 3] = -s.diff(H, w)
        P[2, 4] = s.diff(H, z)
        P = P - P.T
        divergence = s.Matrix([sum(s.diff(P[i, j], X[i]) for i in range(5)) for j in range(5)])
        P = P - (s.Matrix(X) * divergence.T - divergence * s.Matrix(X).T) / 5
        columns.append(s.Matrix([s.Poly(s.expand(P[i, j]), *X).coeff_monomial(X[a] * X[b]) for i, j in PAIRS for a, b in MONS]))
    return s.Matrix.hstack(*columns)


def saved_lp(data):
    return Ring({(0, *m): Q(c) for m, c in data})


def main(output):
    started = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    signal.alarm(120)
    output.mkdir(parents=True, exist_ok=False)
    base = ROOT / 'reports/AMBIENT_CA_ORIGIN_SELF_FLAG_3JET_20261103_v63'
    producer = ROOT / 'reports/AMBIENT_CA_ORIGIN_SELF_FLAG_4JET_20261103_v63'
    physical = ROOT / 'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
    data = read(physical)
    old = read(base / 'certificate.json.gz')
    setup = read(base / 'setup.json')
    current = read(producer / 'certificate.json.gz')
    producer_report = read(producer / 'report.json')
    assert producer_report['all_passed']
    for name, digest in producer_report['input_sha256'].items():
        assert sha(ROOT / name) == digest, name
    for name, digest in producer_report['evidence_sha256'].items():
        assert sha(producer / name) == digest, name

    C = reconstruct_CA()
    assert C == matrix(data['fixed_fiber_linear_basis150_by20'])
    selected = list(C.T.rref()[1])
    assert len(selected) == 20
    left = C[selected, :].inv()
    projection = s.eye(150) - C * left * s.eye(150)[selected, :]
    assert projection * C == s.zeros(150, 20)
    free = data['affine_slice_free_indices']
    assert free == [1, 2, 3, 9, 10, 14, 15, 17, 18, 19]
    full_flags = list(map(tuple, data['flag_positions']))
    assert len(full_flags) == 9 and (1, 0) in full_flags
    flags = [p for p in full_flags if p != (1, 0)]
    assert flags == list(map(tuple, setup['second_and_third_order_flag8']))
    M = s.zeros(5)
    M[1, 0], M[0, 3] = 1, -2
    assert M == matrix(setup['tangent_M'])
    P0 = C * s.Matrix([s.Rational(v) for v in data['center_fixed20']])
    F = C[:, free]
    directions = [first_derivative(M, F[:, j]) for j in range(10)]
    for i, j in flags:
        A = s.zeros(5)
        A[i, j] = 1
        directions.append(first_derivative(A, P0))
    D = projection * s.Matrix.hstack(*directions)
    assert D == matrix(setup['all150_outside_fiber_order2_linear_matrix'])
    pc, pr = setup['pivot_columns'], setup['pivot_rows']
    inverse = matrix(setup['inverse_pivot'])
    assert D[pr, pc] * inverse == s.eye(15) == inverse * D[pr, pc]
    assert D == D[:, pc] * inverse * D[pr, :]
    T = matrix(old['order2_affine_embedding'])
    assert T.shape == (18, 3) and D * T == s.zeros(150, 3)
    assert T[old['order2_free_original_columns'], :] == s.eye(3)
    d = [Ring(Q(v)) for v in old['order2_affine_origin']]
    t, ti = Ring({(0, 1, 0, 0, 0): 1}), Ring({(0, -1, 0, 0, 0): 1})
    ys = [Ring({(0, 0, *(int(j == i) for j in range(3))): 1}) for i in range(3)]
    low = addv(d, constmul(T, [-7 * t / 4 - 1 + 3 * ti / 4, t - 1, Ring()]))
    assert low == list(map(saved_lp, current['full_order2_solution']))

    def group(B_values, C_values=None, fourth=None):
        g = [[Ring(i == j) + Ring(fraction(M[i, j])).shift(1) for j in range(5)] for i in range(5)]
        for values, degree in ((B_values, 2), (C_values, 3), (fourth, 4)):
            if values is not None:
                for k, (i, j) in enumerate(flags):
                    g[i][j] = g[i][j] + values[k].shift(degree)
        return g

    def source(first, second=None, third=None):
        v = [Ring(fraction(c)) for c in P0]
        for values, degree in ((first, 1), (second, 2), (third, 3)):
            if values is not None:
                v = addv(v, [c.shift(degree) for c in constmul(F, values)])
        return v

    # Derive intrinsic3 directly; do not use the saved Lie/log/quotient matrix.
    early = direct_push(group(low[10:]), source(low[:10]))
    intrinsic3 = constmul(projection, [v.coefficient(3) for v in early])
    high = constmul(T, ys)
    correction = constmul(inverse, [intrinsic3[j] for j in pr])
    for k, col in enumerate(pc):
        high[col] = high[col] - correction[k]
    assert high == list(map(saved_lp, current['full_order3_solution_with_all3_higher_freedoms']))
    assert not any(addv(constmul(D, high), intrinsic3))

    actual = direct_push(group(low[10:], high[10:]), source(low[:10], high[:10]))
    projected = [constmul(projection, [v.coefficient(degree) for v in actual]) for degree in range(5)]
    assert all(not any(projected[k]) for k in range(4))
    fourth = projected[4]
    assert fourth == list(map(saved_lp, current['all150_order4_before_fresh18_elimination']))
    obstruction = addv(fourth, negv(constmul(D[:, pc], constmul(inverse, [fourth[j] for j in pr]))))
    assert obstruction == list(map(saved_lp, current['all150_order4_coker']))
    monomials = sorted({m for v in obstruction for m in v})
    coefficients = s.Matrix([[s.Rational(v.get(m, Q(0)).numerator, v.get(m, Q(0)).denominator) for m in monomials] for v in obstruction])
    eqrows = list(coefficients.T.rref()[1])
    assert eqrows == [46, 63, 72, 95] == current['independent_coker_rows']
    pivots = list(coefficients[eqrows, :].rref()[1])
    assert len(pivots) == 4
    selected_inverse = coefficients[eqrows, pivots].inv()
    assert coefficients == coefficients[:, pivots] * selected_inverse * coefficients[eqrows, :]
    expected = [-ys[2] * (t * t - 3) * ti / 5,
                3 * (t * t - 1) ** 2 * ti / 2,
                -(4 * t * t * ys[0] + (7 * t * t + 3) * ys[1]) * ti / 2,
                -3 * ys[2] * (3 * t * t + 1) * ti / 20]
    q1, q2, q3, q4 = [obstruction[j] for j in eqrows]
    assert expected == [q1, q2, q3, q4]
    generators = [ys[2], (t * t - 1) ** 2, 4 * t * t * ys[0] + (7 * t * t + 3) * ys[1]]
    assert generators == [3 * t * q1 / 2 - 2 * t * q4 / 3, 2 * t * q2 / 3, -2 * t * q3]
    # Also verify reverse containment by the explicit multiples.
    assert [q1, q2, q3, q4] == [-(t * t - 3) * ti * generators[0] / 5, 3 * ti * generators[1] / 2, -ti * generators[2] / 2, -3 * (3 * t * t + 1) * ti * generators[0] / 20]

    # Scheme structure, including the nonzero square-zero parameter and CRT.
    u = s.symbols('u')
    f = (u * u - 1) ** 2
    plus = (2 + 3 * u - u ** 3) / 4
    assert s.rem(plus * plus - plus, f, u) == 0
    assert s.rem(plus - 1, (u - 1) ** 2, u) == 0 and s.rem(plus, (u + 1) ** 2, u) == 0
    assert s.rem((u * u - 1) ** 2, f, u) == 0 and s.rem(u * u - 1, f, u) != 0
    assert s.gcd(u, f) == 1 and s.degree(f, u) == 4
    for sign in (1, -1):
        e = s.symbols('e')
        local = s.rem((7 * (sign + e) ** 2 + 3) * (1 - 2 * sign * e) / 4, e ** 2, e)
        assert s.expand(local - (s.Rational(5, 2) - s.Rational(3, 2) * sign * e)) == 0

    witness_evidence = []
    for saved in current['finite4jet_witnesses']:
        sign = int(saved['t'])
        values = [sign, 0, 0, 0]
        solved = constmul(inverse, [fourth[j] for j in pr])
        fresh = [Ring() for _ in range(18)]
        for k, col in enumerate(pc):
            fresh[col] = -solved[k].evaluate(values)
        assert fresh == [Ring(Q(v)) for v in saved['fresh_order4_source10_and_flag8']]
        if sign == 1:
            assert fresh[:10] == [Ring(v) for v in (0, 0, 0, 0, 0, 0, 2, 0, 0, Q(1, 3))]
        else:
            assert not any(fresh)
        assert not any(fresh[10:])
        G = [[v.evaluate(values) for v in row] for row in group(low[10:], high[10:], fresh[10:])]
        P = [v.evaluate(values) for v in source(low[:10], high[:10], fresh[:10])]
        transformed = direct_push(G, P)
        for degree in range(5):
            assert not any(constmul(projection, [v.coefficient(degree) for v in transformed]))
        def rational_coeff(v, degree):
            return str(v.get((degree, 0, 0, 0, 0), Q(0)))
        src = [[rational_coeff(v, degree) for degree in range(5)] for v in P]
        target = [[rational_coeff(v, degree) for degree in range(5)] for v in transformed]
        g_saved = [[[rational_coeff(G[i][j], degree) for j in range(5)] for i in range(5)] for degree in range(5)]
        assert src == saved['source_physical150_by5']
        assert target == saved['second_source_physical150_by5']
        assert g_saved == saved['group5_by5_by5']
        witness_evidence.append(dict(t=sign, source_sha256=hashlib.sha256(json.dumps(src).encode()).hexdigest(), second_source_sha256=hashlib.sha256(json.dumps(target).encode()).hexdigest(), all_750_coefficients_direct_checked=True))

    evidence = dict(CA_basis_150_by20_reconstructed=True, all9_flag_directions_retained=True,
                    correction_matrix_shape=[150, 18], correction_rank=15, complete_kernel_dimension=3,
                    direct_Laurent_inverse_substitution_all150_orders0_to4=True,
                    earlier_free3_retained=True, coker_independent_rows=eqrows,
                    exact_Laurent_ideal=['y2', '(t^2-1)^2', '4*t^2*y0+(7*t^2+3)*y1'],
                    both_ideal_containments=True, nonzero_nilpotent_square_zero=True,
                    two_doubled_A1_plus_fresh_A3=True, normalized_scheme_two_doubled_A4=True,
                    finite_witnesses=witness_evidence,
                    formal_or_algebraic_branch_claim=False, whole_raw_jet_moduli_claim=False)
    (output / 'independent-evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
    paths = [Path(__file__), ROOT / f'research/notes/{STEM}.md',
             ROOT / 'code/ambient_CA_origin_self_flag_4jet_20261103_v63.py',
             ROOT / 'research/notes/AMBIENT_CA_ORIGIN_SELF_FLAG_4JET_20261103_v63.md',
             physical, base / 'setup.json', base / 'certificate.json.gz', base / 'report.json',
             producer / 'certificate.json.gz', producer / 'report.json',
             ROOT / 'reports/AMBIENT_CA_ORIGIN_SELF_FLAG_3JET_20261103_v63_REVIEW/report.json']
    result = dict(id=STEM, all_passed=True, accepted_scope='Complete normalized finite ordinary self-flag 4jet; not a formal or algebraic branch.',
                  producer_rerun=False, independent_method='Raw CA construction; direct Laurent gP(g^-1X)g^T, all150 coefficients and complete cokernel; explicit ideal and CRT identities.',
                  input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
                  evidence_sha256={'independent-evidence.json': sha(output / 'independent-evidence.json')},
                  seconds=time.monotonic() - started)
    (output / 'report.json').write_text(json.dumps(result, indent=2) + '\n')
    signal.alarm(0)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)

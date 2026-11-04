#!/usr/bin/env python3
"""One exact universal centered-triangle interpolation/cofactor calculation.

No actual K5 field evaluation, normal Gram, source point, or Groebner solver.
All polynomial arithmetic is independent Fraction dictionary arithmetic.
"""
import argparse
import gzip
import hashlib
import json
import resource
import signal
import time
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ID = 'S61_UNIVERSAL_FIVE_PLANE_INTERPOLATION_20261104_v63'
ZERO = (0,) * 6
ONE = {ZERO: Q(1)}
MAX_TERMS = 150000
METRICS = {'multiplications': 0, 'max_intermediate_terms': 1}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add(a, b, scalar=Q(1)):
    c = dict(a)
    for m, v in b.items():
        c[m] = c.get(m, Q(0)) + scalar * v
        if not c[m]:
            del c[m]
    assert len(c) <= MAX_TERMS
    METRICS['max_intermediate_terms'] = max(METRICS['max_intermediate_terms'], len(c))
    return c


def scale(a, scalar):
    return {m: v * scalar for m, v in a.items() if v * scalar}


def mul(a, b):
    METRICS['multiplications'] += 1
    c = {}
    for ma, va in a.items():
        for mb, vb in b.items():
            m = tuple(x + y for x, y in zip(ma, mb))
            c[m] = c.get(m, Q(0)) + va * vb
    c = {m: v for m, v in c.items() if v}
    assert len(c) <= MAX_TERMS
    METRICS['max_intermediate_terms'] = max(METRICS['max_intermediate_terms'], len(c))
    return c


def product(polys):
    p = ONE
    for a in polys:
        p = mul(p, a)
    return p


def variable(i):
    e = list(ZERO)
    e[i] = 1
    return {tuple(e): Q(1)}


def encode(poly):
    return [[list(m), str(v)] for m, v in sorted(poly.items())]


def invert(a):
    n = len(a)
    b = [[Q(a[i][j]) for j in range(n)] + [Q(i == j) for j in range(n)] for i in range(n)]
    for j in range(n):
        i = next(i for i in range(j, n) if b[i][j])
        b[i], b[j] = b[j], b[i]
        c = b[j][j]
        b[j] = [v / c for v in b[j]]
        for i in range(n):
            if i != j:
                c = b[i][j]
                b[i] = [v - c * w for v, w in zip(b[i], b[j])]
    assert all(b[i][j] == (i == j) for i in range(n) for j in range(n))
    return [row[n:] for row in b]


def substitute_h(p, value):
    assert all(not any(m[4:]) for m in p)
    powers = [ONE]
    for k in range(max((m[3] for m in p), default=0)):
        powers.append(mul(powers[-1], value))
    ans = {}
    for m, v in p.items():
        coeff = {m[:3] + (0, 0, 0): v}
        ans = add(ans, mul(coeff, powers[m[3]]))
    return ans


def verify_identity(target, cofactors, fs):
    rhs = {}
    for c, f in zip(cofactors, fs):
        rhs = add(rhs, mul(c, f))
    assert rhs == target


def make_rules(fs, a, b, sigma):
    matrix = [[-1, 1, -1], [-1, 2, 0], [0, 2, -1]]
    mi = invert(matrix)
    assert all(sum(Q(matrix[i][k]) * mi[k][j] for k in range(3)) == (i == j)
               for i in range(3) for j in range(3))
    rhs = [add(add(ONE, sigma), mul(sigma, sigma), Q(-2)), a, b]
    borders = [mul(a, a), mul(a, b), mul(b, b), mul(a, sigma), mul(b, sigma), product([sigma] * 3)]
    rules = []
    for i in range(3):
        normal = {}
        for j in range(3):
            normal = add(normal, rhs[j], -mi[i][j])
        cof = [scale(ONE, mi[i][j]) for j in range(3)] + [{}, {}]
        rules.append((next(iter(borders[i]))[3:], normal, cof))
    rules.append(((1, 0, 1), {}, [{}, {}, {}, ONE, {}]))
    rules.append(((0, 1, 1), {}, [{}, {}, {}, {}, ONE]))
    normal = scale(add(sigma, mul(sigma, sigma)), Q(1, 2))
    cof = [scale(sigma, Q(-1, 2)), {}, {}, scale(add(b, a, Q(-1)), Q(1, 2)), scale(b, Q(-1, 2))]
    rules.append(((0, 0, 3), normal, cof))
    for border, normal, cof in rules:
        target = add({(0, 0, 0) + border: Q(1)}, normal, Q(-1))
        verify_identity(target, cof, fs)
    return matrix, mi, rules


def reduce_with_cofactors(p, rules):
    """Fixed graded-lex triangular reduction, linear over QQ[alpha,beta,gamma]."""
    work, normal = dict(p), {}
    cofs = [{} for _ in range(5)]
    steps = 0
    while work:
        plane = max((m[3:] for m in work), key=lambda x: (sum(x),) + x)
        coeff = {m[:3] + (0, 0, 0): v for m, v in work.items() if m[3:] == plane}
        work = {m: v for m, v in work.items() if m[3:] != plane}
        rule = next((r for r in rules if all(x >= y for x, y in zip(plane, r[0]))), None)
        if rule is None:
            normal = add(normal, mul(coeff, {(0, 0, 0) + plane: Q(1)}))
            continue
        border, remainder, c = rule
        quotient = tuple(x - y for x, y in zip(plane, border))
        multiplier = mul(coeff, {(0, 0, 0) + quotient: Q(1)})
        work = add(work, mul(multiplier, remainder))
        cofs = [add(u, mul(multiplier, v)) for u, v in zip(cofs, c)]
        steps += 1
        assert steps <= 56
    return normal, cofs, steps


def main(out):
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    signal.alarm(60)
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    alpha, beta, gamma, a, b, sigma = [variable(i) for i in range(6)]
    h_symbol = a
    h = add(add(mul(alpha, a), mul(beta, b)), mul(gamma, sigma))
    theta = [alpha, beta, scale(add(alpha, beta), Q(-1)), scale(gamma, Q(-1, 2)), gamma]
    points = [[1, 0, 0], [0, 1, 0], [-1, -1, 0], [0, 0, Q(-1, 2)], [0, 0, 1]]
    d = product([add(theta[i], theta[j], Q(-1)) for i in range(5) for j in range(i + 1, 5)])
    assert d
    p_cleared = [{} for _ in range(3)]
    for i in range(5):
        others = [j for j in range(5) if j != i]
        d_over_di = scale(product([add(theta[k], theta[l], Q(-1)) for k in others for l in others if k < l]), Q((-1) ** i))
        di = product([add(theta[i], theta[j], Q(-1)) for j in others])
        assert mul(d_over_di, di) == d
        li = mul(d_over_di, product([add(h_symbol, theta[j], Q(-1)) for j in others]))
        for k in range(3):
            p_cleared[k] = add(p_cleared[k], li, Q(points[i][k]))
    g = product([add(h_symbol, t, Q(-1)) for t in theta])
    assert {m: v for m, v in g.items() if m[3] == 5} == {(0, 0, 0, 5, 0, 0): Q(1)}
    point_checks = []
    for i in range(5):
        for k in range(3):
            assert substitute_h(p_cleared[k], theta[i]) == scale(d, Q(points[i][k]))
        assert not substitute_h(g, theta[i])
        point_checks.append({'point': [str(v) for v in points[i]], 'theta': encode(theta[i]), 'p_values': [str(v) for v in points[i]], 'g_value': '0'})
    identity_h = add(add(mul(alpha, p_cleared[0]), mul(beta, p_cleared[1])), mul(gamma, p_cleared[2]))
    assert identity_h == mul(d, h_symbol)
    aa, ab, bb, ss = mul(a, a), mul(a, b), mul(b, b), mul(sigma, sigma)
    fs = [add(add(add(add(add(ONE, sigma), ss, Q(-2)), aa, Q(-1)), ab), bb, Q(-1)),
          add(add(a, aa, Q(-1)), ab, Q(2)), add(add(b, ab, Q(2)), bb, Q(-1)), mul(a, sigma), mul(b, sigma)]
    for point in points:
        for f in fs:
            assert sum(v * Q(point[0]) ** m[3] * Q(point[1]) ** m[4] * Q(point[2]) ** m[5] for m, v in f.items()) == 0
    matrix, mi, rules = make_rules(fs, a, b, sigma)
    targets = [add(mul(d, z), substitute_h(p, h), Q(-1)) for z, p in zip([a, b, sigma], p_cleared)] + [substitute_h(g, h)]
    all_cofs, reduction_steps = [], []
    for target in targets:
        normal, cof, steps = reduce_with_cofactors(target, rules)
        assert not normal
        verify_identity(target, cof, fs)
        all_cofs.append(cof)
        reduction_steps.append(steps)
    assert not add(add(mul(alpha, targets[0]), mul(beta, targets[1])), mul(gamma, targets[2]))
    for i in range(5):
        assert not add(add(mul(alpha, all_cofs[0][i]), mul(beta, all_cofs[1][i])), mul(gamma, all_cofs[2][i]))
    assert all(sum(m[:3]) == 10 and sum(m[3:]) <= 4 for p in targets[:3] for m in p)
    assert all(sum(m[:3]) == 5 and sum(m[3:]) <= 5 for m in targets[3])
    assert all(sum(m[3:]) <= (2 if k < 3 else 3) for k, cofs in enumerate(all_cofs) for p in cofs for m in p)
    saved_polys = [d, g] + p_cleared + fs + targets + [p for cofs in all_cofs for p in cofs] + [p for _, rem, cof in rules for p in [rem] + cof]
    total_terms = sum(map(len, saved_polys))
    max_bits = max(max(abs(v.numerator).bit_length(), v.denominator.bit_length()) for p in saved_polys for v in p.values())
    assert total_terms <= 30000 and max_bits <= 512
    cert = dict(id=ID, variables=['alpha', 'beta', 'gamma', 'A_or_formal_h', 'B', 'sigma'],
                formal_h_convention='In g5_h and p_cleared_h ONLY, exponentposition3 is formal h; substitute h=alpha*A+beta*B+gamma*sigma.',
                f5=[encode(p) for p in fs], theta=[encode(t) for t in theta], points=point_checks,
                Vandermonde=encode(d), g5_h=encode(g), p_cleared_h=[encode(p) for p in p_cleared],
                targets=[encode(p) for p in targets], cofactors_4x5=[[encode(p) for p in cofs] for cofs in all_cofs],
                identities='targets[k]=sum_i cofactors[k][i]*f5[i]; target0..2=D*(coordinate-p_coordinate(h)); target3=g5(h).',
                coherence='alpha*p_A+beta*p_B+gamma*p_sigma=D*h AND alpha*C_Ai+beta*C_Bi+gamma*C_sigma_i=0 for ALL5i.',
                M3=[[str(v) for v in row] for row in matrix], M3_inverse=[[str(v) for v in row] for row in mi],
                borders=[{'leading': list(border), 'normal': encode(rem), 'cofactors': [encode(p) for p in cof]} for border, rem, cof in rules],
                reduction_steps=reduction_steps,
                actual_K5_specialization_gate='NOTcheckedhere. Need accepted exact actualplane->standard triangle/ell GL3 map before alpha=x+tau*z,beta=y+tau*w,gamma=tau² specialization.',
                no_actual_K5_coefficients_or_normal_expansion=True)
    raw = (json.dumps(cert, indent=2) + '\n').encode()
    assert len(raw) <= 12 * 1024 * 1024
    with gzip.GzipFile(filename=str(out / 'interpolation-certificate.json.gz'), mode='wb', mtime=0) as z:
        z.write(raw)
    inputs = [Path(__file__), ROOT / f'research/notes/{ID}.md',
              ROOT / 'reports/S61_SPARSE_GENERIC5_FIVE_MODULE_PREFLIGHT_20261104_v63/report.json',
              ROOT / 'reports/S61_SPARSE_GENERIC5_FIVE_MODULE_PREFLIGHT_20261104_v63/planned-universal-plane-interpolation.json']
    report = dict(id=ID, all_passed=True, exact_QQ_identities=True, all_five_point_values=True,
                  four_full_f5_cofactor_identities=True, five_cofactor_coherence_identities=True,
                  fixed_linear_border_reduction=True, canonical_only_actual_adapter_not_evaluated=True,
                  no_K5_Gram_D_or_membership_computation=True,
                  term_counts=dict(total_saved_polynomial_occurrences=total_terms, Vandermonde=len(d),
                                  p_cleared_h=list(map(len, p_cleared)), g5_h=len(g), targets=list(map(len, targets)),
                                  cofactor_rows=[[len(p) for p in cofs] for cofs in all_cofs]),
                  max_coefficient_bits=max_bits, metrics=METRICS, reduction_steps=reduction_steps,
                  artifact_uncompressed_bytes=len(raw), artifact_gzip_bytes=(out / 'interpolation-certificate.json.gz').stat().st_size,
                  input_sha256={str(p.relative_to(ROOT)): sha(p) for p in inputs},
                  evidence_sha256={'interpolation-certificate.json.gz': sha(out / 'interpolation-certificate.json.gz')},
                  seconds=time.monotonic() - start)
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)

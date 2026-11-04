#!/usr/bin/env python3
"""Independent fixed integer-convolution checks; no producer or solver imports."""
import argparse
import gzip
import hashlib
import json
import math
import resource
import signal
import time
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ID = 'S61_UNIVERSAL_FIVE_AND_MODULE_20261104_v63_REVIEW'
Z = (0,) * 6
ONE = {Z: Q(1)}
INPUTS = {}
CHECKS = {}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(path):
    p = ROOT / path
    INPUTS[str(p.relative_to(ROOT))] = sha(p)
    return json.loads(gzip.decompress(p.read_bytes()) if p.suffix == '.gz' else p.read_bytes())


def gate(path):
    d = read(path)
    assert d['all_passed'], path
    for p, expected in d.get('input_sha256', {}).items():
        p = Path(p)
        p = p if p.is_absolute() else ROOT / p
        assert sha(p) == expected, str(p)
        INPUTS[str(p.relative_to(ROOT))] = expected
    for p, expected in d.get('evidence_sha256', {}).items():
        p = Path(p)
        p = p if p.is_absolute() else (ROOT / path).parent / p
        assert sha(p) == expected, str(p)
        INPUTS[str(p.relative_to(ROOT))] = expected
    return d


def decode(p):
    d = {tuple(m): Q(v) for m, v in p}
    assert len(d) == len(p) and all(v and len(m) == 6 and min(m) >= 0 for m, v in d.items())
    return d


def add(*args):
    d = {}
    for p in args:
        for m, v in p.items():
            d[m] = d.get(m, 0) + v
    return {m: v for m, v in d.items() if v}


def scaled(p, c):
    c = Q(c)
    return {m: v * c for m, v in p.items() if c}


def var(k):
    m = list(Z)
    m[k] = 1
    return {tuple(m): Q(1)}


def mul(p, q):
    """Clear both operands once, convolve integers, restore exact QQ coefficients."""
    a = math.lcm(*(v.denominator for v in p.values())) if p else 1
    b = math.lcm(*(v.denominator for v in q.values())) if q else 1
    d = {}
    for m, x in p.items():
        xx = x.numerator * (a // x.denominator)
        for n, y in q.items():
            e = tuple(m[i] + n[i] for i in range(6))
            d[e] = d.get(e, 0) + xx * y.numerator * (b // y.denominator)
    assert len(d) <= 30000
    return {e: Q(v, a * b) for e, v in d.items() if v}


def product(ps):
    ans = ONE
    for p in ps:
        ans = mul(ans, p)
    return ans


def formal_h(p, h):
    assert all(m[4:] == (0, 0) for m in p)
    powers = [ONE]
    for k in range(max((m[3] for m in p), default=0)):
        powers.append(mul(powers[-1], h))
    return add(*(mul({m[:3] + (0, 0, 0): v}, powers[m[3]]) for m, v in p.items()))


def point_value(p, point):
    d = {}
    for m, v in p.items():
        e = m[:3] + (0, 0, 0)
        d[e] = d.get(e, 0) + v * math.prod(point[i] ** m[i + 3] for i in range(3))
    return {m: v for m, v in d.items() if v}


def check_interpolation(d, report):
    al, be, ga, a, b, s = [var(i) for i in range(6)]
    fs = [add(ONE, s, scaled(mul(s, s), -2), scaled(mul(a, a), -1), mul(a, b), scaled(mul(b, b), -1)),
          add(a, scaled(mul(a, a), -1), scaled(mul(a, b), 2)),
          add(b, scaled(mul(a, b), 2), scaled(mul(b, b), -1)), mul(a, s), mul(b, s)]
    assert fs == [decode(p) for p in d['f5']]
    theta = [al, be, scaled(add(al, be), -1), scaled(ga, Q(-1, 2)), ga]
    assert theta == [decode(p) for p in d['theta']]
    vand = product([add(theta[i], scaled(theta[j], -1)) for i in range(5) for j in range(i + 1, 5)])
    assert vand and vand == decode(d['Vandermonde'])
    g = product([add(a, scaled(t, -1)) for t in theta])
    assert g == decode(d['g5_h'])
    ps = [decode(p) for p in d['p_cleared_h']]
    points = [[Q(1), Q(0), Q(0)], [Q(0), Q(1), Q(0)], [-Q(1), -Q(1), Q(0)],
              [Q(0), Q(0), -Q(1, 2)], [Q(0), Q(0), Q(1)]]
    for i, point in enumerate(points):
        assert [Q(x) for x in d['points'][i]['point']] == point
        assert all(not point_value(f, point) for f in fs)
        assert all(formal_h(ps[k], theta[i]) == scaled(vand, point[k]) for k in range(3))
        assert not formal_h(g, theta[i])
    # Degrees <=4 make these fifteen interpolation values a uniqueness proof.
    assert all(max((m[3] for m in p), default=0) <= 4 for p in ps)
    assert add(mul(al, ps[0]), mul(be, ps[1]), mul(ga, ps[2])) == mul(vand, a)
    h = add(mul(al, a), mul(be, b), mul(ga, s))
    targets = [add(mul(vand, v), scaled(formal_h(p, h), -1)) for v, p in zip([a, b, s], ps)] + [formal_h(g, h)]
    assert targets == [decode(p) for p in d['targets']]
    cs = [[decode(p) for p in row] for row in d['cofactors_4x5']]
    assert len(cs) == 4 and all(len(row) == 5 for row in cs)
    for k, row in enumerate(cs):
        assert add(*(mul(c, f) for c, f in zip(row, fs))) == targets[k]
        assert all(sum(m[3:]) <= (2 if k < 3 else 3) for p in row for m in p)
    for i in range(5):
        assert not add(mul(al, cs[0][i]), mul(be, cs[1][i]), mul(ga, cs[2][i]))
    assert all(sum(m[:3]) == 10 and sum(m[3:]) <= 4 for p in targets[:3] for m in p)
    assert all(sum(m[:3]) == 5 and sum(m[3:]) <= 5 for m in targets[3])
    expected = [(2, 0, 0), (1, 1, 0), (0, 2, 0), (1, 0, 1), (0, 1, 1), (0, 0, 3)]
    for border, leading in zip(d['borders'], expected):
        assert tuple(border['leading']) == leading
        rem = decode(border['normal'])
        assert add({(0, 0, 0) + leading: Q(1)}, scaled(rem, -1)) == add(*(mul(decode(c), f) for c, f in zip(border['cofactors'], fs)))
        assert all(m[3:] in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 2)] for m in rem)
    matrix = [[Q(x) for x in row] for row in d['M3']]
    inverse = [[Q(x) for x in row] for row in d['M3_inverse']]
    assert matrix == [[-1, 1, -1], [-1, 2, 0], [0, 2, -1]]
    for x, y in [(matrix, inverse), (inverse, matrix)]:
        assert [[sum(x[i][k] * y[k][j] for k in range(3)) for j in range(3)] for i in range(3)] == [[Q(i == j) for j in range(3)] for i in range(3)]
    det = sum(matrix[0][j] * ((-1) ** j) * (matrix[1][(j + 1) % 3] * matrix[2][(j + 2) % 3] - matrix[1][(j + 2) % 3] * matrix[2][(j + 1) % 3]) * ((-1) ** j) for j in range(3))
    assert det == 3
    assert d['no_actual_K5_coefficients_or_normal_expansion'] and report['canonical_only_actual_adapter_not_evaluated']
    CHECKS['canonical_QQ_arithmetic'] = dict(five_fibers=True, fifteen_interpolation_values=True,
        four_full_cofactor_identities=True, five_cofactor_coherences=True, six_border_identities=True,
        M3_both_inverse_products=True, determinant=str(det), no_producer_reduction_rerun=True)
    return fs, points


def check_recipe(fs, points):
    base = 'reports/S61_SPARSE_GENERIC5_FIVE_MODULE_PREFLIGHT_20261104_v63/'
    plan = read(base + 'source-bound-next-plane-recipe.json')
    counts = read(base + 'support-and-dimension-preflight.json')
    target = read(base + 'finite-target-theorem.json')
    toy = read(base + 'toy-countermodel.json')
    r = read(plan['source_files']['residual129'])
    b = read(plan['source_files']['full134_blocks'])
    f = read(plan['source_files']['full134_direct'])
    op = read(plan['source_files']['normal_operator'])
    selected = [50, 10, 52, 51, 15]
    assert plan['selected_five_rows'] == r['selected_original134_rows'] == op['selected_original134_rows'] == selected
    assert plan['selected_plane10_roots'] == [b['rows134'][i]['plane_z_degree_le2'] for i in selected]
    assert plan['adjB5_roots'] == r['adjugate5_nodes'] and plan['inverseDelta5_root'] == r['inverseDelta5_node']
    assert plan['pivot_q_rows'] == f['pivot_q_rows'] == r['pivot_q_rows'] == [3, 6, 13]
    assert len(set(f['normal_q_rows42'])) == 42 and len(f['F134']) == len(b['rows134']) == 134
    assert op['residual130_row_order'] == r['residual_original134_rows'] + ['G6']
    assert f['all_inverses_parameter_only'] and r['all_inverses_parameter_only'] and plan['NO_normal_dependent_inverse']
    assert not plan['actual_new_field_polynomials_or_coefficients_computed']
    assert counts['global_D_shape'] == [5, 675] and counts['first_normaljet_relation_matrix_slots'] == 141750
    assert counts['all_normal_monomials_degree5'] == math.comb(46, 5) == 1370754
    assert counts['generic_dense_normal_monomials_upto6'] == math.comb(48, 6)
    assert counts['generic_dense_normal_monomials_upto10'] == math.comb(52, 10)
    # Bigrading, rather than independent maxima, proves the two degree bounds.
    residual_terms = [(4, 1), (0, 5)]
    selected_terms = [(8, 0), (4, 4), (0, 8), (4, 1), (0, 5)]
    bound = lambda terms: max(n + 5 * ((h + k) // 5) for h, n in terms for k in range(5))
    assert bound(residual_terms) == 6 and bound(selected_terms) == 10
    # Verify both idempotent fibers of the exact toy, without a GB calculation.
    polynomials = [{tuple(e): Q(v) for e, v in p} for p in toy['generators']]
    def n_value(poly, x):
        result = {}
        for e, v in poly.items():
            result[e[:3]] = result.get(e[:3], 0) + v * x ** e[3]
        return {e: v for e, v in result.items() if v}
    assert [n_value(p, Q(0)) for p in polynomials[:5]] == [{m[3:]: v for m, v in p.items()} for p in fs]
    assert all(not n_value(p, Q(0)) for p in polynomials[5:])
    p1 = [n_value(p, Q(1)) for p in polynomials]
    assert p1[5:8] == [{(1, 0, 0): Q(1)}, {(0, 1, 0): Q(1)}, {(0, 0, 1): Q(1)}]
    assert all((0, 0, 0) not in p for p in p1)
    assert sorted(p[0] + 4 * p[2] for p in points) == [-2, -1, 0, 1, 4]
    assert target['encoded_tau_degree'] == 41 and target['no_membership_annihilator_or_exhaustion_computed']
    CHECKS['source_recipe_and_finite_theorem'] = dict(selected_five_roots_and_adjugate_bound=True,
        all135_rows_kept=True, faithful_transcendental_primitive_element=True,
        parameter_only_inverse_recipe=True, coherent_h_reverse_map=True,
        cross_word_reverse_order=True, residual_bigrading_degree6=True,
        selected_bigrading_degree10=True, five_by675_recipe_only=True,
        finite42_targets_and_single_tau_target_accepted=True, finite_colon5_accepted=True,
        exact_toy_nonempty_residual=True, actual_K5_adapter_and_coefficients_still_unconstructed=True)


def main(out):
    start = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    interpolation = 'reports/S61_UNIVERSAL_FIVE_PLANE_INTERPOLATION_20261104_v63/'
    rep = gate(interpolation + 'report.json')
    gate('reports/S61_SPARSE_GENERIC5_FIVE_MODULE_PREFLIGHT_20261104_v63/report.json')
    for stem in ['S61_SPARSE_GENERIC5_NORMAL42_EXACT_PRESENTATION', 'S61_SPARSE_GENERIC5_RESIDUAL129_CIRCUIT', 'S61_SPARSE_GENERIC5_G6_RESIDUAL']:
        gate(f'reports/{stem}_20261104_v63_REVIEW/report.json')
    d = read(interpolation + 'interpolation-certificate.json.gz')
    fs, points = check_interpolation(d, rep)
    check_recipe(fs, points)
    for p in [Path(__file__), ROOT / f'research/notes/{ID}.md']:
        INPUTS[str(p.relative_to(ROOT))] = sha(p)
    report = dict(id=ID, all_passed=True, input_sha256=INPUTS, evidence_sha256={}, checks=CHECKS,
        independent_fixed_QQ_integer_convolution=True, proof_scope='Canonical universal coefficients and actual-source-bound construction/target theorem only; no actual K5 interpolation, normal coefficients, 675-column matrix or membership witness constructed.',
        no_points_primes_GB_producer_solver_or_unknown_dependent_inverse=True,
        seconds=time.monotonic() - start, max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    p = out / 'report.json'
    p.write_text(json.dumps(report, indent=2) + '\n')
    assert p.stat().st_size < 32 * 1024 * 1024
    print(json.dumps(dict(all_passed=True, report=str(p), sha256=sha(p), seconds=report['seconds'], max_rss_kib=report['max_rss_kib'])))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output)

#!/usr/bin/env python3
"""Independent cubic bordered-chart and quartic unit-deletion audit.

No author producer is executed.  QQ point linear algebra uses Fraction
Gauss elimination.  Parametric checks are coefficient cross-products of
the 35 lifted identities, not expanded parametric determinants.
"""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_NONSQUARE_CUBIC_QUARTIC_STRUCTURAL_REDUCTIONS_20261103_v63_REVIEW'


def read(p): return json.loads(p.read_text())


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()


def select(m, rows, cols): return [[m[i][j] for j in cols] for i in rows]


def mm(a, b):
    bt = list(zip(*b))
    return [[sum((x*y for x, y in zip(row, col) if x and y), Q(0)) for col in bt] for row in a]


def gauss(a, inverse=False):
    n = len(a)
    assert n and all(len(row) == n for row in a)
    m = [list(row)+([Q(i == j) for j in range(n)] if inverse else []) for i, row in enumerate(a)]
    det = Q(1)
    for j in range(n):
        k = next((i for i in range(j, n) if m[i][j]), None)
        assert k is not None
        if k != j: m[k], m[j] = m[j], m[k]; det = -det
        pivot = m[j][j]; det *= pivot
        m[j] = [x/pivot for x in m[j]]
        for i in range(n):
            if i != j and m[i][j]:
                c = m[i][j]
                m[i] = [x-c*y for x, y in zip(m[i], m[j])]
    return (det, [row[n:] for row in m]) if inverse else det


def add(p, q):
    r = p.copy()
    for ex, c in q.items(): r[ex] = r.get(ex, Q(0))+c
    return {ex: c for ex, c in r.items() if c}


def mul(p, q):
    r = {}
    for ex, c in p.items():
        for ey, d in q.items():
            z = tuple(x+y for x, y in zip(ex, ey)); r[z] = r.get(z, Q(0))+c*d
    return {ex: c for ex, c in r.items() if c}


ONE = {(0, 0, 0): Q(1)}


def main(out):
    started = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000)); signal.alarm(120)
    out.mkdir(parents=True, exist_ok=False)
    cubic = ROOT/'reports/S61_NONSQUARE_CUBIC_NINE_ROW_CRITERION_20261103_v63'
    quartic = ROOT/'reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63'
    basis = ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62'
    review = ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW/report.json'
    source_path = ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json'
    br, accepted, source = read(basis/'report.json'), read(review), read(source_path)
    assert accepted['all_passed'] and accepted['exact_all6_QQ_relations_all91_columns']
    assert accepted['exact_all63_known_unit_reconstructions'] and accepted['normalization_uses_only_a_F_B_C']
    assert [accepted[k] for k in ['global_free_row_module_rank', 'quadric_product_row_submodule_free_rank', 'abstract_row_module_quotient_free_rank']] == [57, 41, 16]
    reports = [cubic/'report.json', quartic/'report.json', basis/'report.json', review,
        ROOT/'reports/S61_NONSQUARE_CUBIC_NINE_ROW_COUNT_CLARIFICATION_20261103_v63.json']
    for rp in reports:
        rr = read(rp); assert rr['all_passed']
        for rel, h in rr['input_sha256'].items(): assert sha(ROOT/rel) == h, rel
        if 'evidence_sha256' in rr:
            for rel, h in rr['evidence_sha256'].items(): assert sha(rp.parent/rel) == h, rel
    circuit, qc = read(cubic/'circuit.json'), read(quartic/'certificate.json')
    mons, labels, original_ids = source['normal_monomials'], source['labels'], br['remaining_original_row_indices']
    light = [0, 1, 2, 3, 4, 6, 7]
    heavy = [i for i, mon in enumerate(mons) if mon[5]]
    pure = [i for i, mon in enumerate(mons) if not mon[5]]
    assert len(heavy) == 7 and len(pure) == 84 and len(mons) == 91
    assert all(mons[j][5] == 1 and sum(mons[j]) == 2 for j in heavy)
    assert all(sum(mons[j]) == 3 for j in pure)
    assert all(labels[i][0] == 'product' for i in original_ids[:41])
    assert all(labels[i][0] == 'cubic' for i in original_ids[41:])
    assert circuit['source_row_original_ids'] == original_ids and circuit['heavy_columns'] == heavy
    point = tuple(circuit['point']); assert point == (1, 2, 3)
    point_matrix = []
    count = 0; maxdegree = 0; digest = hashlib.sha256()
    with gzip.open(basis/'remaining.txt.gz', 'rb') as f:
        header = f.readline(); digest.update(header); assert header.strip() == b'57 91 57'
        for i in range(57):
            row = []
            for j in range(91):
                line = f.readline(); assert line; digest.update(line)
                ts = line.split(); n = int(ts[0]); assert len(ts) == 1+4*n
                assert not (i < 41 and j in heavy and n)
                value = Q(0)
                for k in range(n):
                    c = Q(ts[1+4*k].decode()); ex = tuple(map(int, ts[2+4*k:5+4*k])); assert min(ex) >= 0
                    value += c*point[0]**ex[0]*point[1]**ex[1]*point[2]**ex[2]
                    maxdegree = max(maxdegree, sum(ex))
                count += n; row.append(value)
            point_matrix.append(row)
        assert not f.read()
    assert count == 164333 and maxdegree == 16
    assert digest.hexdigest() == br['remaining_matrix_sha256'] == circuit['source_uncompressed_sha256']
    lines = (cubic/'point-matrix.txt').read_text().splitlines()
    assert lines[0] == '57 91' and list(map(int, lines[1].split())) == heavy
    assert [[Q(x) for x in line.split()] for line in lines[2:]] == point_matrix
    ip, gr, rr, cc = [circuit[k] for k in ['product_pivot_columns', 'heavy_pivot_rows_basis_positions', 'remaining_rows_basis_positions', 'remaining_columns']]
    assert len(ip) == len(set(ip)) == 41 and set(ip) <= set(pure)
    assert len(gr) == len(set(gr)) == 7 and set(gr) <= set(range(41, 57))
    assert rr == [i for i in range(41, 57) if i not in gr]
    assert cc == [j for j in pure if j not in ip] and len(cc) == 43
    rows48, cols48 = list(range(41))+gr, ip+heavy
    assert circuit['pivot48_rows'] == rows48 and circuit['pivot48_columns'] == cols48
    assert circuit['product_rows_basis_positions'] == list(range(41))
    D = select(point_matrix, rows48, cols48)
    dp = gauss(select(point_matrix, list(range(41)), ip))
    dh = gauss(select(point_matrix, gr, heavy))
    determinant, inverse = gauss(D, True)
    assert determinant == dp*dh
    assert mm(D, inverse) == [[Q(i == j) for j in range(48)] for i in range(48)]
    correction = mm(mm(select(point_matrix, rr, cols48), inverse), select(point_matrix, rows48, cc))
    schur = [[point_matrix[i][j]-correction[ii][jj] for jj, j in enumerate(cc)] for ii, i in enumerate(rr)]
    saved = (cubic/'point-schur.txt').read_text().splitlines()
    assert saved[0] == 'EXACT_QQ_NINE_ROW_CRITERION 41 7 48 9 57'
    assert [list(map(int, row.split())) for row in saved[1:5]] == [ip, gr, rr, cc]
    local = circuit['point_nonzero_9minor_local_columns']
    assert list(map(int, saved[5].split())) == local and len(set(local)) == len(local) == 9
    ds = gauss(select(schur, range(9), local))
    assert list(map(Q, saved[6:10])) == [dp, dh, determinant, ds]
    assert list(map(Q, saved[10:])) == [c for row in schur for c in row]
    assert Q(circuit['point_DeltaP']) == dp and Q(circuit['point_DeltaH']) == dh
    assert Q(circuit['point_Delta48']) == determinant and Q(circuit['point_Schur9_determinant']) == ds
    assert Q(circuit['point_border9_determinant']) == ds*determinant**9
    assert circuit['entry_count'] == 387 and circuit['entry_determinant_size'] == 49
    assert circuit['bordered_entry_total_degree_upper'] == 49*maxdegree
    # A direct bordered determinant independently checks sign and conventions.
    border = select(point_matrix, rows48+[rr[0]], cols48+[cc[0]])
    assert gauss(border) == determinant*schur[0][0]
    ca = read(ROOT/'reports/S61_NONSQUARE_CUBIC_NINE_ROW_COUNT_CLARIFICATION_20261103_v63.json')
    assert ca['correct_matrix_entry_count'] == 57*91 == 5187 and ca['previous_prose_entry_count'] == 5190

    # All35 lifted rational relations, independently constructed by monomial maps.
    a, e, d = s.symbols('a e d'); vs = (a, e, d)
    def pd(expr):
        return {tuple(ex): Q(str(c)) for ex, c in s.Poly(expr, *vs, domain=s.QQ).terms() if c}
    def rational(text):
        num, den = s.sympify(text, locals={'a': a, 'e': e, 'd': d}).as_numer_denom()
        return pd(num), pd(den)
    known = [s.Poly(f, *vs, domain=s.QQ) for f in [a, 500*a+90*e**2+3*e-24,
        (50*a+18*e**2-9*e-2)*d+60*e-40, 100*a+18*e**2-9*e-2]]
    def isunit(p):
        q = s.Poly.from_dict({ex: s.Rational(c.numerator, c.denominator) for ex, c in p.items()}, vs, domain=s.QQ)
        assert not q.is_zero
        for f in known:
            while True:
                div, rem = q.div(f)
                if not rem.is_zero: break
                q = div
        assert q.total_degree() == 0
    def parse(line):
        ts = line.split(); n = int(ts[0]); assert len(ts) == 1+4*n
        return {tuple(map(int, ts[2+4*k:5+4*k])): Q(ts[1+4*k]) for k in range(n)}
    norms = (basis/'normalization.txt').read_text().splitlines(); assert norms[0] == '63' and len(norms) == 64
    factors = [parse(line) for line in norms[1:]]
    for factor in factors: isunit(factor)
    exponent = lambda indices: tuple(indices.count(i) for i in range(8))
    w2 = sorted([exponent(list(z)) for z in itertools.combinations_with_replacement(light, 2)]+[exponent([5])])
    w4 = sorted([exponent(list(z)) for z in itertools.combinations_with_replacement(light, 4)]+
        [exponent([5]+list(z)) for z in itertools.combinations_with_replacement(light, 2)]+[exponent([5, 5])])
    assert len(w2) == 29 and len(w4) == 239 and qc['weight4_columns'] == [list(z) for z in w4]
    chart_path = ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json'; chart = read(chart_path)
    q4 = [i for i, weight in enumerate(chart['equation_weights']) if weight == 4 and i not in chart['selected34_rows']]
    cubic_original = [i for i, weight in enumerate(chart['equation_weights']) if weight == 3 and i not in chart['selected34_rows']]
    assert [label[1] for label in labels[42:]] == cubic_original
    assert len(q4) == 33 and qc['genuine_quartic_original_rows_retained'] == q4
    rows = [['quadric', i, list(mon)] for i in range(6) for mon in w2]
    rows += [['cubic_normalized', i, labels[i][1], k] for i in range(42, 63) for k in light]
    rows += [['quartic', i] for i in q4]
    assert len(rows) == 354 and qc['source_rows'] == rows
    lookup = {(i, k): 174+(i-42)*7+light.index(k) for i in range(42, 63) for k in light}
    source_cert = read(ROOT/'code/data/s61_binary_nonsquare_cubic_row_basis_20261103_v62.json')
    deleted, triangular = [], ONE
    relation_number = 0; rational_count = 0
    for relation, original, pivot in zip(source_cert['relations'][:5], [48, 55, 59, 57, 61], [1, 1, 1, a, -a]):
        assert relation['deleted_original_row'] == original
        coeffs = {term['original_row']: {tuple(map(int, z[1:])): Q(z[0]) for z in term['terms']} for term in relation['coefficients']}
        assert coeffs[original] == pd(pivot)
        for k in light:
            expected = {}
            for i, co in coeffs.items():
                if i < 42:
                    kind, qi, j = labels[i]; assert kind == 'product'
                    index = qi*29+w2.index(exponent([j, k])); den = factors[i]
                else: index = lookup[i, k]; den = ONE
                if index in expected:
                    num0, den0 = expected[index]
                    expected[index] = add(mul(num0, den), mul(co, den0)), mul(den0, den)
                else: expected[index] = co, den
            record = qc['lifted_relations'][relation_number]; relation_number += 1
            drop = lookup[original, k]
            assert record['deleted_row'] == drop and record['deleted_source'] == rows[drop]
            assert pd(s.sympify(record['pivot'], locals={'a': a})) == pd(pivot)
            supplied = {i: rational(c) for i, c in record['coefficients']}
            assert len(supplied) == len(record['coefficients']) and set(supplied) == set(expected)
            for i, (num, den) in supplied.items():
                en, ed = expected[i]
                assert mul(num, ed) == mul(en, den)
                isunit(den); rational_count += 1
            assert supplied[drop][0] == mul(pd(pivot), supplied[drop][1])
            assert all(i not in supplied for i in deleted)
            triangular = mul(triangular, pd(pivot)); deleted.append(drop)
    assert relation_number == 35 and rational_count == 728 and len(set(deleted)) == 35
    assert deleted == qc['deleted_row_indices']
    retained = [i for i in range(354) if i not in deleted]
    assert len(retained) == 319 and retained == qc['retained_row_indices']
    assert triangular == {(14, 0, 0): Q(-1)} and qc['triangular35minor'] == '-a^14'
    # Every239 coordinate of each relation is a monomial translate of an
    # accepted all91-column cubic identity; no new quartic expansion is needed.
    maps = []
    for k in light:
        image = [w4.index(tuple(x+int(j == k) for j, x in enumerate(mon))) for mon in mons]
        assert len(image) == len(set(image)) == 91
        maps.append(image)
    evidence = dict(cubic_actual_shape=[57, 91], matrix_entry_count=5187, actual_coefficient_terms=count,
        maximum_entry_degree=maxdegree, independent_Fraction_QQ_point=[1, 2, 3], all387_schur_coefficients_match=True,
        full48_inverse_and_DeltaP_times_DeltaH_checked=True, nonzero9minor_and_border_scale_checked=True,
        direct49_border_sign_checked=True, extra_chart_units=['DeltaP', 'DeltaH'],
        quartic_original_rows=354, retained_rows=319, columns=239, all35_rational_relations_reconstructed=True,
        all728_nonzero_rational_coefficients_cross_multiplied=True, all63_normalizers_and728_denominators_aFBC_units=True,
        seven_injective_91_to_239_monomial_column_maps=maps, deleted_rows=deleted,
        triangular_minor='-a^14', all33_genuine_quartics_retained=True,
        source_dependency='Accepted independent exact all6 cubic identities/all91 columns and all63 source reconstructions',
        no_d_J_product_lower41_heavy7_or_full57_assumption_for_quartic_row_span_equality=True,
        no_full57_finiteness_quartic_surjectivity_length65_or_ambient_claim=True)
    (out/'independent-certificate.json').write_text(json.dumps(evidence, indent=2)+'\n')
    result = dict(id=STEM, all_passed=True, cubic_9_by_43_bordered_criterion_on_extra_chart_accepted=True,
        exact_QQ_point_only_proves_chart_nonempty=True, full57_all_points_or_defect_finite_claim=False,
        quartic354_to319_row_span_equal_on_aFBC=True, all35_lifted_rational_row_deletions_independently_checked=True,
        deleted_minor='-a^14', quartic_surjectivity_or_ambient_claim=False,
        no_GB_new_field_or_parametric_minor_expansion=True, seconds=time.monotonic()-started)
    paths = reports+[Path(__file__), ROOT/f'research/notes/{STEM}.md', cubic/'circuit.json', cubic/'point-schur.txt',
        cubic/'point-matrix.txt', cubic/'module-countermodel.json', quartic/'certificate.json', source_path,
        basis/'remaining.txt.gz', basis/'normalization.txt', chart_path,
        ROOT/'code/data/s61_binary_nonsquare_cubic_row_basis_20261103_v62.json']
    result['input_sha256'] = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    result['evidence_sha256'] = {'independent-certificate.json': sha(out/'independent-certificate.json')}
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n'); signal.alarm(0)
    print(json.dumps({k: v for k, v in result.items() if k not in ['input_sha256', 'evidence_sha256']}, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)

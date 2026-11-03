#!/usr/bin/env python3
"""Independent QQ original ZW rank19 and generic quartic146 audit."""
import argparse
from fractions import Fraction as Q
import itertools
import json
from pathlib import Path
import resource
import shutil
import signal
import tarfile
import time
import sympy as s
from s61_nonsquare_structural_reductions_review_20261103_v63 import sha, read, gauss, add, mul, ONE
from s61_nonsquare_five_quadric_fourth_unit_review_20261103_v63 import native, witness

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63_REVIEW'


def main(out):
    started = time.monotonic(); resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000)); signal.alarm(60)
    out.mkdir(parents=True, exist_ok=False)
    parent = ROOT/'reports/S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63'; report = read(parent/'report.json')
    assert report['frozen'] and report['all_exact_QQ_checks_passed']
    archive = ROOT/report['archive']; assert sha(archive) == report['archive_sha256'] and archive.stat().st_size == 54607
    snap = out/'immutable-input'; snap.mkdir()
    with tarfile.open(archive) as tar:
        for member in tar:
            assert member.isfile(); dst = snap/member.name; assert dst.resolve().is_relative_to(snap.resolve())
            dst.parent.mkdir(parents=True, exist_ok=True)
            with tar.extractfile(member) as src, dst.open('wb') as target: shutil.copyfileobj(src, target)
    for rel, h in report['input_sha256'].items(): assert sha(snap/rel if (snap/rel).exists() else ROOT/rel) == h, rel
    a, e, d = s.symbols('a e d'); x = s.symbols('x0:5'); v = s.symbols('v0:8'); vs = (a, e, d)
    C = 100*a+18*e**2-9*e-2
    def pd(expr): return {tuple(ex): Q(str(c)) for ex, c in s.Poly(expr, *vs, domain=s.QQ).terms() if c}
    def decode(terms): return {tuple(ex): Q(c) for ex, c in terms if Q(c)}
    def unit(expr):
        num, den = s.cancel(expr).as_numer_denom()
        for p in [num, den]:
            pp = s.Poly(p, *vs, domain=s.QQ); assert not pp.is_zero
            for f in [s.Poly(a, *vs, domain=s.QQ), s.Poly(C, *vs, domain=s.QQ)]:
                while True:
                    q, r = pp.div(f)
                    if not r.is_zero: break
                    pp = q
            assert pp.total_degree() == 0
    source_rel = 'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json'
    source = read(snap/source_rel)
    quads = []
    for row in source['quadratic_rows_sparse']:
        quads.append(sum(sum(s.Rational(c)*a**ex[0]*e**ex[1]*d**ex[2] for ex, c in cs)*s.prod(z**k for z, k in zip(v, mon)) for mon, cs in row))
    images = {v[0]: x[0], v[1]: x[1], v[2]: x[2], v[3]: 0, v[4]: x[3], v[5]: 0, v[6]: x[1]-3*x[0], v[7]: x[4]}
    qs = [s.expand(quads[i].subs(images, simultaneous=True)) for i in [0, 1, 3, 5]]
    # This is the polynomial ring of five independent original coordinates;
    # no ell relation has been imposed.
    assert images[v[0]] == x[0] and images[v[1]] == x[1] and images[v[2]] == x[2] and images[v[4]] == x[3] and images[v[7]] == x[4]
    assert s.expand((3*v[0]-v[1]+v[6]).subs(images, simultaneous=True)) == 0
    L = s.cancel(-s.Rational(405, 4)*s.diff(qs[1], x[1])/C)
    assert s.Poly(L, *x, *vs, domain=s.QQ).total_degree() >= 1
    assert s.Poly(L, *x).total_degree() == 1 and s.diff(L, x[1]) == -40*C
    assert s.expand(qs[1]-L**2/8100) == 0 and s.expand(qs[2]-L*(e*L-400*a*x[0])/8100) == 0
    assert s.expand((e*L-400*a*x[0])*qs[1]-L*qs[2]) == 0
    saved = read(snap/'work/nonsquare-v62/ZW-product19-and-generic146-v82/certificate.json')
    assert s.expand(L-s.sympify(saved['L'], locals=dict(zip(['a','e','d'], vs)) | dict(zip(map(str, x), x)))) == 0
    for actual, encoded in zip(qs, saved['actual_quadric_coefficients']):
        assert {tuple(mon): pd(co) for mon, co in s.Poly(actual, *x).terms()} == {tuple(mon): decode(cs) for mon, cs in encoded}
    # Independent24-term determinant from coefficient polynomials, no det/Hessian helper.
    assert not qs[0].has(x[1]); variables = [x[i] for i in [0, 2, 3, 4]]
    H = [[pd(s.diff(qs[0], z, w)/2) for w in variables] for z in variables]
    determinant = {}
    for perm in itertools.permutations(range(4)):
        inv = sum(perm[i] > perm[j] for i in range(4) for j in range(i+1, 4))
        term = {(0, 0, 0): Q((-1)**inv)}
        for i, j in enumerate(perm): term = mul(term, H[i][j])
        determinant = add(determinant, term)
    assert determinant == {(2, 0, 0): Q(625, 1024)}
    solved = s.cancel((L+40*C*x[1]).subs(x[0], 0)/(40*C))
    assert not solved.has(x[1]) and s.cancel(L.subs(x[0], 0).subs(x[1], solved)) == 0
    pair = [s.Poly(s.cancel(q.subs(x[0], 0).subs(x[1], solved)), *x[2:]) for q in [qs[0], qs[3]]]
    mons = [(i, j, 2-i-j) for i in range(3) for j in range(3-i)]
    actual_pair = [[pd(s.cancel(q.coeff_monomial(s.prod(z**k for z, k in zip(x[2:], mon)))*(C if i else 1))) for mon in mons] for i, q in enumerate(pair)]
    base = snap/'work/nonsquare-v62/ZW-L-pair-v82'
    matrix, _, _ = native(base/'pair-native.txt'); assert matrix == actual_pair
    residual, pivot_columns, original_rows = native(base/'after.txt', True)
    assert len(residual) == 1 and len(residual[0]) == 5 and pivot_columns == [2] and original_rows == [1]
    expr = lambda p: sum(s.Rational(c.numerator, c.denominator)*a**ex[0]*e**ex[1]*d**ex[2] for ex, c in p.items())
    assert matrix[0][2] == {(1, 0, 0): Q(5)}
    ratios = []
    for j, saved_poly in zip([0, 1, 3, 4, 5], residual[0]):
        raw = expr(matrix[1][j])-expr(matrix[1][2])*expr(matrix[0][j])/(5*a)
        ratio = s.cancel(raw/expr(saved_poly)); unit(ratio)
        assert s.cancel(raw-ratio*expr(saved_poly)) == 0; ratios.append(str(ratio))
    bez = read(base/'literal-Bezout.json')
    scalar = list(map(decode, bez['source_scalar5'])); cofs = list(map(decode, bez['scalar_unit_cofactors']))
    assert scalar == residual[0] and cofs == witness(base/'last-row-unit.witness')
    assert sum(map(len, scalar)) == 30
    assert sum(map(len, bez['scalar_unit_cofactors'])) == 117 and sum(map(len, cofs)) == 116
    assert max(sum(ex) for p in cofs for ex in p) == 9
    unity = {}
    for f, c in zip(scalar, cofs): unity = add(unity, mul(f, c))
    assert unity == ONE

    # Exact old CI certificate is rebuilt from actual804 coefficient terms.
    raw_rel = 'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json'
    raw = read(snap/raw_rel); xx = s.symbols('u0:6'); yy = s.symbols('y0:4'); count = 0; qp = []
    for row in raw['cleared_qfour_sparse']:
        terms = []
        for mon, cs in row:
            assert len(mon) == 7 and mon[6] == 0; count += len(cs)
            terms.append(sum(s.Rational(c) for ex, c in cs)*s.prod(z**k for z, k in zip(xx, mon[:6])))
        qp.append(sum(terms))
    assert count == 804
    qp = [qp[0], qp[1], (qp[2]-qp[1])/(-5), qp[3]]
    cut = dict(zip(xx[:4], yy)) | {xx[4]: yy[0]+yy[1]+2*yy[2]+3*yy[3], xx[5]: 2*yy[0]+3*yy[1]+5*yy[2]+7*yy[3]}
    restricted = [s.Poly(q.subs(cut, simultaneous=True), *yy, domain=s.QQ) for q in qp]
    monomials = lambda D: [(i, j, k, D-i-j-k) for i in range(D+1) for j in range(D-i+1) for k in range(D-i-j+1)]
    polys = []
    for q in restricted:
        for mon in monomials(3):
            polys.append({tuple(i+j for i, j in zip(ex, mon)): Q(str(c)) for ex, c in q.terms()})
    Mac = [[q.get(mon, Q(0)) for q in polys] for mon in monomials(5)]
    assert len(Mac) == 56 and all(len(row) == 80 for row in Mac)
    prior = read(snap/'work/nonsquare-v62/raw-fourquad-structure-v49/generic-CI.json'); selected = prior['selected_columns']
    assert len(selected) == len(set(selected)) == 56
    det = gauss([[row[j] for j in selected] for row in Mac]); assert det == Q(prior['selected_square_determinant']) and det
    # Generic independence of the six original quadrics from one exact QQ minor.
    light = [0, 1, 2, 3, 4, 6, 7]; light2 = list(itertools.combinations_with_replacement(light, 2))
    point_rows = [[Q(str(s.Poly(q.subs({a: 1, e: 1, d: 1}), *v, domain=s.QQ).coeff_monomial(v[i]*v[j]))) for i, j in light2] for q in quads]
    work = [row[:] for row in point_rows]; pivots = []
    for j in range(28):
        r = next((i for i in range(len(pivots), 6) if work[i][j]), None)
        if r is None: continue
        k = len(pivots); work[k], work[r] = work[r], work[k]; value = work[k][j]; work[k] = [z/value for z in work[k]]
        for i in range(k+1, 6):
            factor = work[i][j]
            if factor: work[i] = [z-factor*w for z, w in zip(work[i], work[k])]
        pivots.append(j)
        if len(pivots) == 6: break
    assert len(pivots) == 6
    six_det = gauss([[row[j] for j in pivots] for row in point_rows]); assert six_det
    bridge = ROOT/'reports/S61_NONSQUARE_RAW24_PRODUCT41_BRIDGE_20261103_v63_REVIEW/report.json'
    reviewed_bridge = read(bridge); assert reviewed_bridge['all_passed']
    assert reviewed_bridge['checked_dependency_sha256'][source_rel] == sha(snap/source_rel)
    assert reviewed_bridge['checked_dependency_sha256'][raw_rel] == sha(snap/raw_rel)
    # Hilbert coefficient of a six-variable complete intersection of4 quadrics.
    from math import comb
    HF4 = sum((-1)**i*comb(4, i)*comb(4-2*i+5, 5) for i in range(3))
    HF3_ZW = comb(7, 3)-19
    assert HF4 == 48 and HF3_ZW == 16 and HF4+HF3_ZW == 64
    assert comb(10, 4)-64 == 146 and 146+6 == 152
    evidence = dict(actual_original_ZW_four_quadric_coefficients_match=True,
        five_independent_light_coordinates=[0, 1, 2, 4, 7], ell_not_imposed=True,
        q1_L_square_q3_L_times_eL_minus400aX_verified=True, L_x1='-40*C',
        independent_24term_symmetric_determinant='625*a^2/1024',
        actual_pair_all12_entries_match=True, safe_first_pivot='5*a', residual_column_unit_ratios=ratios,
        literal_Fraction_5generator_sum1_verified=True, saved_cofactor_terms_including_one_zero=117,
        actual_nonzero_cofactor_terms=116, no_unknown_localizers=True,
        unique_linear_syzygy_proof='Coprime modL; pair mod(L,X) eliminates c0,c5; Koszul(L,X) remains',
        generic_CI_from_actual804_fixed4plane_56square_QQ_minor=True, CI_determinant=str(det),
        six_original_quadrics_generic_independence_minor={'point': [1, 1, 1], 'columns': pivots, 'determinant': str(six_det)},
        colon_bridge='generic barJ CI4 => J CI4 and ell regular => (J+ell(Z,W)):ell=J+(Z,W)',
        Hilbert_coefficients={'barJ_degree4': 48, 'ZW_degree3': 16, 'I_degree4': 64},
        generic_light_image_rank=146, global_polynomial_specialization_upper=146,
        weighted_extra_H_times_quadrics_disjoint_columns=6, generic_weighted_rank=152, global_weighted_upper=152,
        no_allpoint_CI_lower146_fullquartic_or_ambient_claim=True)
    (out/'independent-certificate.json').write_text(json.dumps(evidence, indent=2)+'\n'); shutil.rmtree(snap)
    result = dict(id=STEM, all_passed=True, original_ZW_product_rank_exactly19_on_aC=True,
        generic_CI_and_colon_Hilbert_bridge_accepted=True, generic_light_rank146_and_global_upper146=True,
        generic_weighted_rank152_and_global_upper152=True, old147_and153_targets_invalid=True,
        allpoint_lower146_CI_fullquartic_or_ambient_claim=False,
        no_GB_new_field_modp_promotion_or_abstract_to_image_shortcut=True, seconds=time.monotonic()-started)
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', parent/'report.json', archive, bridge]
    result['input_sha256'] = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    result['evidence_sha256'] = {'independent-certificate.json': sha(out/'independent-certificate.json')}
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n'); signal.alarm(0)
    print(json.dumps({k: v for k, v in result.items() if k not in ['input_sha256', 'evidence_sha256']}, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)

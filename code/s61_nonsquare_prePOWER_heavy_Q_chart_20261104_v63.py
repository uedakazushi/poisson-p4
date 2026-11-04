#!/usr/bin/env python3
"""Typed ONE actual H5 chart and nine prePOWER Q9 charts; no parameter expansion."""
import argparse
import hashlib
import importlib.util
import itertools
import json
import resource
import signal
import time
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ID = 'S61_NONSQUARE_PREPOWER_HEAVY_Q_CHART_20261104_v63'


def read(p):
    return json.loads(p.read_text())


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def main(out):
    began = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    signal.alarm(60)
    out.mkdir(parents=True, exist_ok=True)
    assert not (out/'report.json').exists(), 'Immutable completed output exists'
    phase = read(out/'phase.json') if (out/'phase.json').exists() else {}
    assert phase.get('phase') == 'source_recipe_persisted_before_builder'
    paths = {}

    def source(stem, filename='certificate.json'):
        p = ROOT/'reports'/stem/filename
        paths[str(p.relative_to(ROOT))] = p
        return read(p), p

    def bind(p):
        paths[str(p.relative_to(ROOT))] = p

    def ref(p, **kw):
        return dict(path=str(p.relative_to(ROOT)), **kw)

    q, qp = source('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63')
    pair, pp = source('S61_NONSQUARE_CUBIC_ZW25_PAIR_UNIT_COMPLETION_20261104_v63')
    prefix, xp = source('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63', 'checkpoint.json')
    hc, hp = source('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63')
    rank, rp = source('S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63')
    h9, h9p = source('S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63', 'prePOWER-h9-source.json')
    cover, cp = source('S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63', 'H5-source-and-chart-profile.json')
    ab, ap = source('S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63', 'input.json')
    unit, up = source('S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63')
    heavy, dp = source('S61_NONSQUARE_CUBIC_PROJECTED15_HEAVY_PREFLIGHT_20261104_v63')
    for stem in [
        'S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63',
        'S61_NONSQUARE_CUBIC_ZW25_PAIR_UNIT_COMPLETION_20261104_v63',
        'S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63',
        'S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63',
        'S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63',
        'S61_NONSQUARE_CUBIC_LAST_Q9_POWER_BEZOUT_20261104_v63',
        'S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63',
    ]:
        review, _ = source(stem+'_REVIEW', 'report.json')
        assert review['all_passed']
    for stem in ['S61_NONSQUARE_CUBIC_PURE_HEAVY_KERNEL_SCOPE_20261104_v63',
                 'S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63']:
        source(stem, 'report.json')
    bind(Path(__file__))
    bind(ROOT/f'research/notes/{ID}.md')
    helper = ROOT/'code/s61_nonsquare_pure_kernel9_charts_20261104_v63.py'
    bind(helper)
    spec = importlib.util.spec_from_file_location('frozen_kernel9_slp', helper)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    o = module.SLP()
    external = []

    def atom(p, key, root):
        if root in (0, 1):
            return root
        n = o.put(('accepted_scalar', str(p.relative_to(ROOT)), key, root))
        external.append(ref(p, key=key, root=root, scalar_node=n))
        return n

    def poly(p, file, key, **more):
        return o.poly(p, ref(file, key=key, **more))

    def summ(xs):
        z = 0
        for x in xs:
            z = o.add(z, x)
        return z

    def dot(x, y):
        assert len(x) == len(y)
        return summ(o.mul(a, b) for a, b in zip(x, y))

    def rational(nd, file, key, **more):
        return o.mul(poly(nd[0], file, [key, 'numerator'], **more),
                     o.inv(poly(nd[1], file, [key, 'denominator'], **more)))

    named_polys = {
        'a': [[[1, 0, 0], '1']], 'd': [[[0, 0, 1], '1']],
        'F': [[[1, 0, 0], '500'], [[0, 2, 0], '90'], [[0, 1, 0], '3'], [[0, 0, 0], '-24']],
        'B': [[[1, 0, 1], '50'], [[0, 2, 1], '18'], [[0, 1, 1], '-9'], [[0, 0, 1], '-2'], [[0, 1, 0], '60'], [[0, 0, 0], '-40']],
        'C': [[[1, 0, 0], '100'], [[0, 2, 0], '18'], [[0, 1, 0], '-9'], [[0, 0, 0], '-2']],
        'J': [[[1, 0, 1], '5'], [[0, 1, 0], '-6'], [[0, 0, 0], '-4']],
    }
    named = {k: o.poly(p, dict(definition='FULLU_named_factor', name=k)) for k, p in named_polys.items()}

    def known_unit(c):
        n = o.const(c['constant'])
        for name, power in zip(['a', 'd', 'F', 'B', 'C', 'J'], c['a_d_F_B_C_J_powers']):
            assert power >= 0
            for _ in range(power):
                n = o.mul(n, named[name])
        return n

    # Cost gate is made before constructing parameter operators.
    prefix_terms = sum(sum(map(len, st['pivot'])) +
                       sum(sum(map(len, nd)) for _, nd in st['eliminations'])
                       for st in prefix['all_completed_unit_steps'])
    preflight = dict(source_shape=[15, 22], pure_shape=[15, 17], heavy_shape=[14, 5],
        prefix_parameter_terms_counting_repeats=prefix_terms,
        new_explicit_H5_terms=sum(len(p) for row in hc['raw_compact14x5'] for p in row),
        raw_HZ_HW_terms=sum(len(p) for row in ab['raw16_HZ_HW_coefficients'] for p in row),
        existing_Delta_terms=len(rank['nonzero_minor_of14x5']), h9_terms=sum(map(len, h9['source_h9'])),
        upper_explicit_terms_with_repeats=25000, hard_explicit_terms=150000,
        hard_seconds=60, hard_address_space_bytes=2000000000, hard_output_bytes=33554432,
        parameter_products_expanded=False, source_pure_entries='external frozen scalar roots',
        one_H5_chart_only=True, no_new_scalar_evaluation_minor_rank_GB_point_prime=True)
    assert prefix_terms == 12752
    assert sum(preflight[k] for k in ['prefix_parameter_terms_counting_repeats', 'new_explicit_H5_terms',
               'raw_HZ_HW_terms', 'existing_Delta_terms', 'h9_terms']) < 25000
    (out/'preflight.json').write_text(json.dumps(preflight, indent=2)+'\n')
    phase.update(phase='preflight_saved_building_typed_operators', no_new_scalar_arithmetic=True)
    (out/'phase.json').write_text(json.dumps(phase, indent=2)+'\n')

    rows15 = q['new15_source_original25_labels']
    rows14 = hc['final14_source_labels']
    cols22 = q['new22_target_original32_columns']
    assert rows15 == [9]+rows14 and len(cols22) == 22 and cols22[17:] == [27, 28, 29, 30, 31]
    pure_cols = cols22[:17]
    pure = [[atom(qp, 'scalar_DAG.nodes', n) for n in row[:17]]
            for row in q['remaining15x22_scalar_roots']]
    assert q['remaining15x22_scalar_roots'][0][17:] == [0]*5
    qrow = pure[0][:]
    h_nodes, aliases = [], []
    for k, (label, pol) in enumerate(zip(h9['original_target_labels'], h9['source_h9'])):
        j = pure_cols.index(label)
        assert q['remaining15x22_scalar_roots'][0][j] == h9['source_prePOWER_DAG_roots'][k]
        n = poly(pol, h9p, 'source_h9', index=k)
        aliases.append(dict(accepted_Q9_root=qrow[j], actual_small_h_root=n,
                            source_DAG_root=h9['source_prePOWER_DAG_roots'][k], pure_position=j))
        qrow[j] = n
        h_nodes.append(n)
    pure[0] = qrow
    rho = [known_unit(x) for x in hc['known_row_units_rho14']]
    gamma = known_unit(hc['known_column_units'][4])
    P = [[o.mul(o.inv(rho[i]), x) for x in pure[i+1]] for i in range(14)]
    N = [[poly(p, hp, 'raw_compact14x5', index=[i, j]) for j, p in enumerate(row)]
         for i, row in enumerate(hc['raw_compact14x5'])]
    I, pi = [1, 3, 6, 12, 9], [0, 4, 1, 3, 2]
    assert [rank['input14_remaining_rows'][i] for i in I] == rank['nonzero_minor_of14x5_rows']
    assert [rank['input5_columns'][j] for j in pi] == rank['nonzero_minor_of14x5_cols']
    Delta = poly(rank['nonzero_minor_of14x5'], rp, 'nonzero_minor_of14x5')
    assert len(rank['nonzero_minor_of14x5']) == 1752 and Delta
    A = [[N[i][j] for j in pi] for i in I]
    detcache = {}

    def det(M):
        key = tuple(tuple(row) for row in M)
        if key not in detcache:
            if not M:
                ans = 1
            elif len(M) == 1:
                ans = M[0][0]
            else:
                ans = summ(o.mul(o.const((-1)**j), o.mul(M[0][j], det([r[:j]+r[j+1:] for r in M[1:]])))
                           for j in range(len(M)))
            detcache[key] = ans
        return detcache[key]

    adj = [[o.mul(o.const((-1)**(i+j)), det([row[:i]+row[i+1:] for row in A[:j]+A[j+1:]]))
            for j in range(5)] for i in range(5)]
    delta_expression = det(A)
    rest = [i for i in range(14) if i not in I]
    K, Ksource, lambdas = [], [], []
    for u in rest:
        lam = [dot([N[u][k] for k in pi], [adj[k][j] for k in range(5)]) for j in range(5)]
        v = [0]*15
        v[u+1] = o.mul(Delta, o.inv(rho[u]))
        for j, i in enumerate(I):
            v[i+1] = o.neg(o.mul(lam[j], o.inv(rho[i])))
        Ksource.append(v)
        K.append([o.sub(o.mul(Delta, P[u][j]), dot(lam, [P[i][j] for i in I])) for j in range(17)])
        lambdas.append(lam)

    # Exact source15 -> original25 operator. This is covector transport, not
    # transposed source matrices. Accepted source matrices act on the LEFT.
    S16 = [[atom(qp, 'scalar_DAG.nodes', n) for n in row] for row in q['source16_roots']]
    S17 = [[atom(pp, 'scalar_DAG.nodes', n) for n in row] for row in pair['source17_roots']]
    old16, old17 = q['old16_source_original25_rows'], pair['old17_source_original25_rows']
    assert rows15 == old16[1:] and old16 == old17[1:]
    assert prefix['remaining_rows'] == old17
    for i in range(16):
        assert q['source16_roots'][i][2:] == [int(j == i) for j in range(2, 16)]
    for i in range(17):
        assert pair['source17_roots'][i][3:] == [int(j == i) for j in range(3, 17)]
    steps = []
    for k, st in enumerate(prefix['all_completed_unit_steps']):
        pivot = known_unit(st['known_unit'])
        cs = [(t, rational(nd, xp, 'all_completed_unit_steps.eliminations', step=k, target_row=t))
              for t, nd in st['eliminations']]
        steps.append(dict(pivot_row=st['pivot_row'], pivot_column=st['pivot_column'], pivot=pivot,
                          inverse_pivot=o.inv(pivot), eliminations=[[t, c] for t, c in cs],
                          accepted_source_pivot_fraction=ref(xp, key='all_completed_unit_steps', index=k)))
    T15 = []
    for z in range(15):
        y16 = [0]+[int(i == z) for i in range(15)]
        y16 = [dot(y16, [S16[i][j] for i in range(16)]) for j in range(16)]
        y17 = [0]+y16
        y17 = [dot(y17, [S17[i][j] for i in range(17)]) for j in range(17)]
        y = [0]*25
        for lab, value in zip(old17, y17):
            y[lab] = value
        for st in reversed(steps):
            p = st['pivot_row']
            for t, c in reversed(st['eliminations']):
                y[p] = o.sub(y[p], o.mul(c, y[t]))
            y[p] = o.mul(y[p], st['inverse_pivot'])
        T15.append(y)
    # The only genuine source change before the two Q completions is the
    # accepted last HL pivot. Earlier pivot rows are Q and have zero H.
    last = steps[-1]
    lastc = {t: c for t, c in last['eliminations']}
    map_g = dict(zip([10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24],
                    [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15]))
    assert last['pivot_row'] == 11 and last['pivot_column'] == 26
    for j, label in enumerate(rows15):
        expected = [0]*15
        if label != 9:
            expected[label-10] = 1
            if label in lastc:
                expected[1] = o.neg(o.mul(lastc[label], last['inverse_pivot']))
        assert T15[j][10:] == expected, ('genuine source transport', j)

    mu = [poly(p, ap, 'genuine16_source_clearing_units', index=i)
          for i, p in enumerate(ab['genuine16_source_clearing_units'])]
    HZW = [[poly(p, ap, 'raw16_HZ_HW_coefficients', index=[i, j]) for j, p in enumerate(row)]
           for i, row in enumerate(ab['raw16_HZ_HW_coefficients'])]
    cz = poly(unit['Bezout_cZ'], up, 'Bezout_cZ')
    cw = poly(unit['Bezout_cW'], up, 'Bezout_cW')
    eta = poly(unit['known_unit_eta'], up, 'known_unit_eta')
    zd = poly(heavy['delta_scaled_HZ_HW_coefficients']['Z'], dp, ['delta_scaled_HZ_HW_coefficients', 'Z'])
    wd = poly(heavy['delta_scaled_HZ_HW_coefficients']['W'], dp, ['delta_scaled_HZ_HW_coefficients', 'W'])
    charts = []
    for hi, (h, label) in enumerate(zip(h_nodes, h9['original_target_labels'])):
        pivotcol = pure_cols.index(label)
        keepcols = [j for j in range(17) if j != pivotcol]
        result = []
        for k, u in enumerate(rest):
            v15 = [o.mul(h, x) for x in Ksource[k]]
            v15[0] = o.neg(K[k][pivotcol])
            image = [o.sub(o.mul(h, K[k][j]), o.mul(K[k][pivotcol], qrow[j])) for j in keepcols]
            pivot_zero = o.sub(o.mul(h, K[k][pivotcol]), o.mul(K[k][pivotcol], h))
            assert pivot_zero == 0
            y25 = [dot(v15, [T15[i][j] for i in range(15)]) for j in range(25)]
            y26 = y25[:10]+[0]*16
            for oldlab, gi in map_g.items():
                y26[10+gi] = y25[oldlab]
            hz = summ(o.mul(o.mul(y26[10+i], mu[i]), HZW[i][0]) for i in range(16))
            hw = summ(o.mul(o.mul(y26[10+i], mu[i]), HZW[i][1]) for i in range(16))
            correction = o.mul(o.add(o.mul(cz, hz), o.mul(cw, hw)), o.inv(eta))
            y26[16] = o.sub(y26[16], correction)
            y26[11] = o.add(y26[11], correction)
            raw_genuine = [o.mul(y26[10+i], mu[i]) for i in range(16)]
            result.append(dict(unselected_H5_row14=u, source15_coefficients=v15,
                original26_coefficients_after_Cdelta=y26, original57_genuine_coefficients=raw_genuine,
                before_Cdelta_HZ=hz, before_Cdelta_HW=hw, Cdelta_correction=correction,
                pure_pivot_zero_root=pivot_zero, pure16_roots=image))
        charts.append(dict(h_index=hi, h_original32_target=label, h_root=h, inverse_h=o.inv(h),
            pure16_target_labels=[pure_cols[j] for j in keepcols], source9=result,
            basis_minor_in_source15='h^9 * Delta^9 * product(rho_unselected^-1)',
            projected_rank9_NOT_computed=True))

    # Universal adjugate identities in 25 independent formal entries. These
    # are NOT expanded parameter polynomials and do not evaluate a minor.
    def fadd(a, b, sign=1):
        z = dict(a)
        for m, c in b.items():
            z[m] = z.get(m, 0)+sign*c
            if not z[m]:
                del z[m]
        return z

    def fmul(a, b):
        z = {}
        for m, c in a.items():
            for n, d in b.items():
                mn = tuple(sorted(m+n))
                z[mn] = z.get(mn, 0)+c*d
        return {m: c for m, c in z.items() if c}

    def fdet(M):
        z = {}
        for perm in itertools.permutations(range(len(M))):
            sign = (-1)**sum(perm[i] > perm[j] for i in range(len(M)) for j in range(i+1, len(M)))
            p = {(): sign}
            for i, j in enumerate(perm):
                p = fmul(p, M[i][j])
            z = fadd(z, p)
        return z

    FA = [[{(5*i+j,): 1} for j in range(5)] for i in range(5)]
    FD = fdet(FA)
    Fadj = [[{m: (-1)**(i+j)*c for m, c in fdet([row[:i]+row[i+1:] for row in FA[:j]+FA[j+1:]]).items()}
             for j in range(5)] for i in range(5)]
    for i in range(5):
        for j in range(5):
            for left in (True, False):
                val = {}
                for k in range(5):
                    val = fadd(val, fmul(FA[i][k], Fadj[k][j]) if left else fmul(Fadj[i][k], FA[k][j]))
                assert val == (FD if i == j else {})
    inverse_nodes = []
    for i, n in enumerate(o.nodes):
        if n[0] == 'inv':
            inverse_nodes.append(i)
    assert o.explicit_terms < 25000 and len(o.nodes) < 150000
    known_localizers = dict(FULLU=['a', 'd', 'F', 'B', 'C', 'J'], rho14=rho,
        gamma=gamma, accepted_8_prefix_pivots=[st['pivot'] for st in steps], eta=eta,
        external_source_Q_completions='Only a and nonzero QQ denominators, accepted unit-completion receipts',
        no_old_normalizer_or_POWER207_cofactor=True)
    chart_localizers = dict(actual_nonzero_generic_Delta=Delta, inverse_Delta=o.inv(Delta),
        actual_h9=h_nodes, all_inverse_h=[ch['inverse_h'] for ch in charts],
        nine_h_charts_cover_FULLU='Reused accepted sum(H_i*h_i)=a^29; degree207 H_i not evaluated',
        selected_H5_chart_only=True, no_H5_chart_complement_removed=True,
        other35_Dp_charts=cover['Dp_35_H5_row_charts'],
        Dp_coverage=cover['Dp_coverage'], p0_warning=cover['p0_warning'],
        support_potential_H5_charts=cover['all_H5_support_potential_charts'],
        support_potential_is_NOT_nonzero_determinant_certificate=True)
    proof = dict(
        accepted_source_binding='All 15x22 entries and both source matrices are immutable accepted finite scalar roots. All eight Q9 tails retained.',
        precondition='H_actual=diag(rho)*N*diag(1,1,1,1,gamma); source G_i/rho_i and target heavy5/gamma give N.',
        adjugate_both_products_50_universal_polynomial_identities=True,
        actual_Delta_alias=dict(formal_determinant_root=delta_expression, accepted_polynomial_root=Delta,
                               ordered_rows14=I, ordered_columns5=pi, old_raw_rows=rank['nonzero_minor_of14x5_rows'],
                               old_raw_columns=rank['nonzero_minor_of14x5_cols']),
        H5_kernel='For u outside I, lambda=N[u,pi]*adj(A); Ksource=Delta*e_u-sum(lambda_i*e_Ii). All5 H entries vanish by both adjugate identities.',
        H5_kernel_source_minor='Delta^9 in normalized G14 source basis; kernel10 includes unchanged Q9.',
        Q_kernel='At h=Q9[q], use h*K_u-K_u[q]*Q9; the pure pivot coefficient is literally zero in every81 row.',
        all81_pivot_zero_roots_are_literal_zero=True,
        source15_to25='Insert deleted Q5=0; right-multiply accepted source16; insert deleted Q4=0; right-multiply source17; embed in25; reverse8 ELIMINATE coefficient operations then divide each pivot coefficient by its known-unit pivot.',
        all15_genuine_subblocks_match_exact_last_HL_correction=True,
        last_HL='Gi_mu - ci/p_HL*G1_mu, p_HL=2*F*B^2*C^2/15625. Prior7 Q pivots have zero genuine and zero H coefficients.',
        Cdelta='Lift from quotient H5 to raw H7 by subtracting (cZ*HZ+cW*HW)/eta*(G6_mu-G1_mu).',
        raw_Hzero_identities=['psi=-wdelta*HZ+zdelta*HW=0 follows from H5zero',
          'eta*HZ-(cZ*HZ+cW*HW)*zdelta=-cW*psi',
          'eta*HW-(cZ*HZ+cW*HW)*wdelta=cZ*psi'],
        Cdelta_pi_and_psi_zero='Accepted normalized duplicate G6-G1; correction does not change any projected pure entry.',
        exact_projected_rank='On FULLU & D(Delta*h): rank(aug25x32)=16+rank(N9x16); rank(old17x24)=8+rank(N9x16).',
        full_cubic_sufficiency='If any new9x16 has row rank9 on its chart, the actual9 rawH7-kernel vectors are independent modulo the pure product rows, hence full57. No rank9/minor calculated here.',
        raw_full_cubic_failure_not_equivalent='The16 projected pure target coordinates do not give the full43 pure quotient. A projected rank drop is not a proof of fullcubic defect.',
        chart_scope='Only actualDelta andactualhi are extra chart localizers; other35/p0 complement retained. No global FULLU57 or quartic239 claim.',
        no_expanded_parameter_product_new_minor_rank_point_prime_GB=True)
    data = dict(scalar_nodes=o.nodes, scalar_type='QQ(a,e,d), typed scalar SLP with accepted external roots',
        external_scalar_semantics='accepted_scalar(path,key,root) is that immutable accepted scalar node; its complete finite parent DAG remains in the pinned input.',
        Q9_h_aliases=aliases, source15_labels=rows15, pure17_target_labels=pure_cols,
        source15_original25_operator=T15, source8_reverse_covector_steps=list(reversed(steps)),
        normalized_N14x5=N, rho14=rho, gamma=gamma, H5_ordered_rows14=I, H5_ordered_columns5=pi,
        H5_matrix5x5=A, H5_adjugate5x5=adj, H5_Delta=Delta, H5_formal_Delta_root=delta_expression,
        unselected_H5_rows14=rest, H5_lambda9x5=lambdas, source_K9x15=Ksource,
        projected_10x17_roots=[qrow]+K, nine_h_charts=charts,
        raw_Cdelta=dict(cZ=cz,cW=cw,eta=eta,zdelta=zd,wdelta=wd,mu16=mu,raw_HZHW16x2=HZW),
        known_FULLU_localizers=known_localizers, extra_chart_localizers=chart_localizers)
    for name, value in [('circuit.json', data), ('symbolic-proof.json', proof)]:
        b = json.dumps(value, separators=(',', ':')).encode()+b'\n'
        assert len(b) < 30*1024*1024
        (out/name).write_bytes(b)
    phase.update(phase='typed_construction_complete', scalar_nodes=len(o.nodes),
                 active_unique_explicit_parameter_terms=o.explicit_terms)
    (out/'phase.json').write_text(json.dumps(phase, indent=2)+'\n')
    total = sum(p.stat().st_size for p in out.iterdir() if p.is_file())
    assert total < 32*1024*1024
    report = dict(id=ID, all_passed=True, frozen=True, independent_review_pending=True,
        one_actual_H5_chart_and_nine_h_Q9_charts_constructed=True,
        source_shapes=dict(prePOWER=[15,22], H5kernel=[10,17], Q9kernel=[9,16], number_Q9_charts=9),
        all81_actual_raw_H7_lifts_saved=True, all8_Q9_tails_retained=True,
        old_normalizer_or_POWER207_factors_used=False, scalar_nodes=len(o.nodes),
        active_unique_explicit_parameter_terms=o.explicit_terms, output_bytes_before_report=total,
        no_new_scalar_coeff_evaluation_minor_rank_GB_point_prime=True,
        no_complement_or_global_full57_exhaustion_claim=True,
        input_sha256={name:sha(p) for name,p in paths.items()},
        evidence_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()},
        seconds=time.monotonic()-began)
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir() if p.is_file()) < 32*1024*1024
    signal.alarm(0)
    print(json.dumps({k:v for k,v in report.items() if k not in ['input_sha256','evidence_sha256']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    a = p.parse_args()
    try:
        main(a.output.resolve())
    except Exception as e:
        if a.output.exists():
            (a.output/'failure.json').write_text(json.dumps(dict(all_passed=False,error=repr(e),no_auto_retry=True),indent=2)+'\n')
        raise

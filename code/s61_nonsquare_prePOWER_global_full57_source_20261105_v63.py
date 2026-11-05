#!/usr/bin/env python3
"""Audit the GLOBAL unit-ring old25/current15/full57 bridge; no new matrix values."""
import argparse
import collections
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

ROOT = Path(__file__).resolve().parents[1]
ID = 'S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63'
Z = (0, 0, 0)


def sha(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def poly(a):
    return {tuple(m): Q(c) for m, c in a if Q(c)}


def enc(a):
    return [[list(m), str(c)] for m, c in sorted(a.items()) if c]


def add(a, b):
    r = a.copy()
    for m, c in b.items():
        r[m] = r.get(m, 0) + c
        if not r[m]:
            del r[m]
    return r


def mul(a, b):
    r = {}
    for m, c in a.items():
        for n, d in b.items():
            k = tuple(x+y for x, y in zip(m, n))
            r[k] = r.get(k, 0) + c*d
            if not r[k]:
                del r[k]
    assert len(r) < 150000
    return r


def power(a, n):
    r = {Z: Q(1)}
    for _ in range(n):
        r = mul(r, a)
    return r


def main(out):
    began = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    signal.alarm(60)
    out.mkdir(parents=True, exist_ok=False)
    pins = {}

    def load(st, name='certificate.json'):
        p = ROOT/'reports'/st/name
        pins[str(p.relative_to(ROOT))] = sha(p)
        return json.loads(p.read_text())

    def accepted(st, flags=()):
        r = load(st, 'report.json')
        assert r['all_passed']
        for k in flags:
            assert r[k] is True, (st, k)
        return r

    def literal_pin(review, rel):
        assert review['input_sha256'][rel] == sha(ROOT/rel), rel

    def dump(name, value):
        (out/name).write_text(json.dumps(value, indent=2)+'\n')

    f = [poly(x) for x in [
        [[[1,0,0], '1']], [[[0,0,1], '1']],
        [[[1,0,0], '500'], [[0,2,0], '90'], [[0,1,0], '3'], [Z, '-24']],
        [[[1,0,1], '50'], [[0,2,1], '18'], [[0,1,1], '-9'], [[0,0,1], '-2'], [[0,1,0], '60'], [Z, '-40']],
        [[[1,0,0], '100'], [[0,2,0], '18'], [[0,1,0], '-9'], [Z, '-2']],
        [[[1,0,1], '5'], [[0,1,0], '-6'], [Z, '-4']]]]

    def unit(c, exps):
        r = {Z: Q(c)}
        assert len(exps) == 6 and all(n >= 0 for n in exps)
        for g, n in zip(f, exps):
            r = mul(r, power(g, n))
        return r

    # Use the actual all-point product theorem, never a generic rank label.
    product = accepted('S61_NONSQUARE_PRODUCT41_FULL_OPEN_20261103_v63', [
        'actual_original42_by84_product_rank_exact41_on_whole_FULLU'])
    pr = accepted('S61_NONSQUARE_PAIR_EMPTY_QQ_20261103_v63_REVIEW', [
        'raw24_module_surjective_on_FULLU_via_accepted_gcd_binary_finiteSigma_reviews'])
    br = accepted('S61_NONSQUARE_RAW24_PRODUCT41_BRIDGE_20261103_v63_REVIEW', [
        'pointwise_RAW24_rank24_implies_original_product_rank41'])
    for st in ['S61_NONSQUARE_PAIR_EMPTY_QQ_20261103_v63_REVIEW',
               'S61_NONSQUARE_RAW24_PRODUCT41_BRIDGE_20261103_v63_REVIEW']:
        literal_pin(product, f'reports/{st}/report.json')
    basis = accepted('S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62', [
        'each_relation_checked_on_all_91_columns', 'each_deleted_coefficient_is_a_known_unit',
        'previous_deleted_coefficients_zero_in_later_relations'])
    accepted('S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW')
    assert basis['deleted_original_zero_based_rows'] == [48,55,59,57,61,14]
    assert basis['primary_triangular_minor'] == '-12*a^3/5'
    assert len([i for i in basis['remaining_original_row_indices'] if i < 42]) == 41

    # Pin the initial physical projection and all coefficient normalizations.
    lf = load('S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63')
    lfn = load('S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63', 'native.json')
    lr = accepted('S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63_REVIEW', [
        'only_known_a_C_unit_normalizations',
        'exactly9_monomial_unit_pivots_and_all26_plus5_heavy_columns_retained'])
    ur = accepted('S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63_REVIEW', [
        'exact_Fraction_polynomial_Bezout_eta_known_FULLU_unit_verified'])
    up = load('S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63')
    hp = load('S61_NONSQUARE_CUBIC_PROJECTED15_HEAVY_PREFLIGHT_20261104_v63')
    fsr = accepted('S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63_REVIEW', [
        'actual_full57_original63_raw91_source_lifts_and_eta_cleanup_accepted'])
    for name in ['operator.json', 'proof.json']:
        literal_pin(fsr, f'reports/S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63/{name}')
    fullop = load('S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63', 'operator.json')
    fullproof = load('S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63', 'proof.json')
    assert fullop['nodes'][169]['shape'] == [91,32]
    assert fullop['nodes'][169]['args'] == [164,168]
    assert fullop['exports']['lift25_to_full57'] == 49
    assert poly(up['known_unit_eta']) == unit('192/125', [1,0,1,2,3,0])
    assert poly(hp['known_unit_delta']) == unit('-64000', [0,0,0,0,2,0])
    del fullop

    aug = load('S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63', 'input.json')
    assert aug['removed_genuine_relative_row'] == 6 and aug['kept_genuine_relative_row'] == 1
    assert aug['duplicate_zero_in_both_pi_and_psi'] is True
    assert len(aug['source25_labels']) == 25
    mu = aug['genuine16_source_clearing_units']
    cancels = lfn['genuine_C_factors_cancelled']
    assert len(cancels) == 16 and lfn['q5_C_factors_cancelled'] == 2
    for i in range(16):
        assert poly(mu[i]) == unit('-64000', [0,0,0,0,3-cancels[i],0])
    assert mu[1] == mu[6]
    aug_meta = {k: aug[k] for k in ['source25_labels', 'target31_recipe',
        'rawHZ_HW_indices_in_frozen91', 'removed_genuine_relative_row',
        'kept_genuine_relative_row', 'full57_sufficient_condition']}
    del aug

    # Exactly eight source operations. Coefficients stay external; only the
    # eight already saved pivots are remultiplied against named unit factors.
    eight = 'S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63'
    er = accepted(eight+'_REVIEW', ['all8_known_unit_fixed_operations_and_all_checkpoint_entries_verified'])
    literal_pin(er, f'reports/{eight}/checkpoint.json')
    cp = load(eight, 'checkpoint.json')
    assert cp['completed_pivots'] == 8
    expected = [(8,20,'-5/4',[0]*6),(1,10,'-5/4',[0]*6),
        (2,16,'5',[1,0,0,0,0,0]),(3,12,'-5/4',[0]*6),
        (0,0,'25/4',[1,0,0,0,0,0]),(6,23,'-1/6',[1,0,0,0,0,0]),
        (7,24,'-1/6',[1,0,0,0,0,0]),(11,26,'2/15625',[0,0,1,2,2,0])]
    step_inventory = []
    for k, (s, e) in enumerate(zip(cp['all_completed_unit_steps'], expected)):
        i,j,c,exps = e
        assert (s['pivot_row'],s['pivot_column']) == (i,j)
        assert Q(s['known_unit']['constant']) == Q(c)
        assert s['known_unit']['a_d_F_B_C_J_powers'] == exps
        assert poly(s['pivot'][0]) == mul(poly(s['pivot'][1]), unit(c,exps))
        assert all(t != i for t,v in s['eliminations'])
        step_inventory.append(dict(step=k,row=i,column=j,unit_constant=c,unit_powers=exps,
            inverse_constant=str(1/Q(c)),inverse_powers=[-n for n in exps],
            pre_elimination_coefficients_pointer=f'/all_completed_unit_steps/{k}/eliminations',
            elimination_count=len(s['eliminations']),
            forward='first R_p <- R_p/u; then R_i <- R_i-c_i*R_p for exactly the saved rows',
            inverse='first R_i <- R_i+c_i*R_p for those rows; then R_p <- u*R_p',
            no_target_or_tail_dropped_until_unit_block_split=True))
    assert cp['selected_minor_unit_determinant'] == {'constant':'-1/4608','a_d_F_B_C_J_powers':[4,0,1,2,2,0]}
    eight_meta = dict(remaining_rows=cp['remaining_rows'],remaining_columns=cp['remaining_columns'],
        source_terms=sum(sum(len(p) for p in v) for row in cp['matrix25x32'] for v in row)
          +sum(sum(len(p) for p in v) for s in cp['all_completed_unit_steps'] for _,v in s['eliminations']))
    assert eight_meta['source_terms'] < 150000
    del cp

    nine = 'S61_NONSQUARE_CUBIC_ZW25_PAIR_UNIT_COMPLETION_20261104_v63'
    ten = 'S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63'
    r9 = accepted(nine+'_REVIEW', ['actual_one_new_unit_and_same_cokernel_verified'])
    r10 = accepted(ten+'_REVIEW', ['actual_GL3_unit_completion_and_same_cokernel_verified',
                                'all3914_DAG_nodes_and330_remaining_entries_verified'])
    literal_pin(r9, f'reports/{nine}/certificate.json')
    literal_pin(r10, f'reports/{ten}/certificate.json')
    n9, n10 = load(nine), load(ten)
    assert n10['scalar_DAG']['nodes'][:2138] == n9['scalar_DAG']['nodes']
    assert n9['old17_source_original25_rows'] == eight_meta['remaining_rows']
    assert n9['old24_target_original32_columns'] == eight_meta['remaining_columns']
    assert n10['old16_source_original25_rows'] == n9['new16_source_original25_rows']
    assert n10['old23_target_original32_columns'] == n9['new23_target_original32_columns']
    assert n9['tail_clear_original_columns'] == [1,2,3]
    assert n10['tail_clear_original32_columns'] == [1,2,3,4,6,7,8,9,11,13,14,15,17,18,19]
    assert n9['all_5_remaining_heavy_and_psi_columns_unchanged'] == n10['all5_remaining_heavy_quotient_columns_unchanged'] == [27,28,29,30,31]
    allowed = {'z','one','a','m','s','checkpoint','accepted_pair_second','new_unit_atom'}
    counts = collections.Counter(n[0] for n in n10['scalar_DAG']['nodes'])
    assert set(counts) <= allowed and not any('inv' in x for x in counts)
    for i,n in enumerate(n10['scalar_DAG']['nodes']):
        if n[0] in {'a','m','s'}:
            assert all(isinstance(k,int) and 0 <= k < i for k in n[1:])
    dual = n9['actual_six_dot_and_binding']['transported_pair_second_column']
    atoms = n10['scalar_DAG']['new_unit_atoms_num_den']
    dual_den_inventory = []
    for group, values in [('pair6',dual),('Q2new6',atoms)]:
        for i, (num, den) in enumerate(values):
            assert len(den) == 1 and tuple(den[0][0])[1:] == (0,0)
            assert Q(den[0][1]) and den[0][0][0] in [0,1,2]
            dual_den_inventory.append(dict(group=group,index=i,denominator=den,
                                          inverse_factors_only_QQ_and_a=True))
    completions = []
    for st,c,skey,tkey,tikey,pivrow,pivcol in [(nine,n9,'source17_roots','target24_roots','target24_inverse_roots',4,5),
         (ten,n10,'source16_roots','target23_roots','target23_inverse_roots',5,21)]:
        S = c[skey]; size=len(S)
        assert all(S[i][i] == 1 for i in range(size))
        assert all(S[i][j] == 0 for i in range(size) for j in range(size) if i != j and j != 0)
        assert all(S[0][j] == int(j == 0) for j in range(size))
        completions.append(dict(source=st,source_size=size,source_forward_pointer='/'+skey,
            source_inverse='I-(S-I); (S-I)^2=0 since only column0 off-diagonal is used',
            target_forward_pointer='/'+tkey,target_inverse_pointer='/'+tikey,
            source_det=c['source_change_determinant'],target_det=c['target_change_determinant'],
            pivot_old25=pivrow,pivot_original32=pivcol,
            tail_columns=c.get('tail_clear_original_columns',c.get('tail_clear_original32_columns')),
            tails_kept_as_actual_source_DAG=True))
    assert len(n10['remaining15x22_scalar_roots']) == 15
    assert all(len(row)==22 for row in n10['remaining15x22_scalar_roots'])
    q9roots = n10['remaining15x22_scalar_roots'][0]
    tail_labels = [j for j in n10['new22_target_original32_columns'] if j not in [22,11,13,25,15,19,14,18,17,27,31,28,30,29]]
    assert len(tail_labels) == 8

    # The additional DIRECT determinant argument is kept source-only here.
    dr = accepted('S61_NONSQUARE_PREPOWER_DIRECT_M12_FULL15_SOURCE_20261105_v63_REVIEW')
    literal_pin(dr,'research/notes/S61_NONSQUARE_PREPOWER_DIRECT_M12_FULL15_SOURCE_20261105_v63.md')
    recipe = load('S61_NONSQUARE_PREPOWER_M12_AXIS_WORKER_GATE_20261105_v63','recipe.json')
    assert all(Q(c) for c in recipe['constants'].values())
    assert recipe['no_kappa3_denominator'] is True
    assert recipe['current22_original_labels'] == n10['new22_target_original32_columns']
    assert recipe['old_row_operation_also_on_newcolumn'] == 'row1+2row2'

    dump('unit-inventory.json',dict(ring='QQ[a,e,d,1/(a*d*F*B*C*J)]',factor_order=['a','d','F','B','C','J'],
        named_factors=[enc(p) for p in f],eight_actual_pivots=step_inventory,
        initial_projection_units=dict(L_coordinate_inverse='(-40*C)^-1',
          quad_generator_change='L^2=8100*q1; L*x0=(8100*e*q1-8100*q3)/(400*a)',
          genuine16_clearing=[dict(row=i,constant='-64000',C_power=3-cancels[i]) for i in range(16)],
          q5_clearing='(-40*C)^2/C^2=1600',delta='-64000*C^2',eta='192*a*F*B^2*C^3/125'),
        completion_leaf_denominators=dual_den_inventory,
        direct_M12_units=dict(L='a*F*C^2',gamma='-4*F*B*C^2/625',
          constants=recipe['constants'],row_extra_units=['a^-1','a^-1','a'],
          H2='L^-1,(L*gamma)^-1',H3='L^-1 on 28,30,29',pure_tail2='B',
          rho13_over_L='8000000/16119',row_precondition='row1+2row2; determinant1'),
        NOT_inverted=['Delta407','Saux','N6','N7bar','P_K92','h22','p','R434','kappa3','A_anchor'],
        later_corollary_localizers_only=['T','D12']))
    dump('global-operator-recipe.json',dict(
        external_old25_source='reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63/input.json#/matrix25x32',
        old25_source_meta=aug_meta,
        physical57_lift='reports/S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63/operator.json#/nodes/49',
        physical91_to_augmented32='same operator /nodes/169 (31-column pi plus fullHZ/HW psi)',
        eight_unit_forward_inverse=step_inventory,
        eight_block_split='after saved pivot row/column permutations: [[U8,V],[0,M17]], diag(U8)=1, U8 upper triangular; U8^-1=sum((-strictU8)^k,k=0..7); target shear uses only this finite sum',
        completions=completions,actual_scalar_DAG_source=ten+'/certificate.json#/scalar_DAG',
        DAG_node_counts=dict(counts),new15_labels=n10['new15_source_original25_labels'],
        new22_labels=n10['new22_target_original32_columns'],all330_actual_matrix_roots=n10['remaining15x22_scalar_roots'],
        remaining_Q9_all22_roots=q9roots,other8_Q9_tail_columns=tail_labels,
        raw7_cleanup='for any row with oldH5=psi=0: subtract ((cZ*hZ+cW*hW)/eta)*Cdelta; all pure tails change exactly by this original-source row',
        rank_identity_at_every_FULLU_point='rank(M25)=10+rank(M15)',
        no_H5_adjugate_or_normalization_by_unknown_factors=True))
    dump('source-proof.json',dict(
        all_point_product41=dict(actual_all42_rank=41,
          independent_raw24_surjection=pr['id'],pointwise_17_plus24_bridge=br['id'],
          unit_row_change_det='80*a',ell_coordinate_unit='-16*C*J/243',
          four_modZW_pivots=['5*a','a*C/240','16*C^2/81','160*a*C/81'],
          nonzero_syzygy_coefficient='10*a',
          chosen41_span_all42='the deleted product14 has coefficient12*a/5 in the last triangular relation, with no deleted genuine rows present'),
        phi_product_factorization='ZW restriction; quotient by (L^2,L*x0) in degree3 removes9 monomials; q2,q4 are linear combinations of q1,q3 modZW by the80a identities; all product images lie in Q10; psi(product)=0',
        Cdelta_phi_zero=True,Cdelta_nonzero_mod_products='cZ*zdelta+cW*wdelta=eta unit; products have all rawH7 coefficients zero',
        genuine_scalings_all_units=True,kept_genuine15=[i for i in range(16) if i != 6],
        proof_at_every_geometric_point=[
          'rank(M25)=25 makes its Q10 rows independent and retainedG15 independent modulo their span',
          'any relation among retainedG15 modulo full products maps by phi to a relation among the25 independent rows, hence all15 coefficients vanish',
          'add Cdelta, whose eta-unit rawH pair prevents its lying in the pure product span; get16 independent classes',
          'allpoint product rank41 gives fullrank>=41+16=57; the six global source relations give<=57'],
        direct_identity='A*D15_normalized=-T*D12 in the named FULLU localization, before any A inversion',
        conditional_new_chart='on FULLU D(T*D12), RHS is a unit, hence A and D15 are units; current15rank15 => old25rank25 => actualfull57',
        independent_global_bridge_review_pending=True,
        old_Omega_and_all_frozen_conditional_files_unchanged=True,
        no_global_FULLU57_quartic239_exhaustion_or_bad_Fitting_ideal_equality=True))
    note = ROOT/f'research/notes/{ID}.md'
    pins[str(note.relative_to(ROOT))] = sha(note)
    pins[str(Path(__file__).relative_to(ROOT))] = sha(Path(__file__))
    total = sum(p.stat().st_size for p in out.iterdir())
    assert total < 32*1024*1024
    report = dict(id=ID,all_passed=True,frozen=True,source_only=True,
        global_rank_old25_equals10_plus_current15=True,
        global_rank_old25_25_implies_actual_full57=True,
        explicit_forward_inverse_named_FULLU_unit_inventory=True,
        all_H3_Q9_and_other8_tails_retained=True,
        conditional_chart_AFTER_independent_global_source_review='FULLU D(T*D12)',
        no_new_scalar_matrix_value_rank_or_determinant=True,
        no_chart_reduction_adopted_before_independent_review=True,
        independent_review_pending=True,largest_loaded_explicit_parameter_terms=eight_meta['source_terms'],
        output_raw_bytes=total,input_sha256=pins,
        evidence_sha256={p.name:sha(p) for p in out.iterdir()},seconds=time.monotonic()-began)
    dump('report.json',report)
    signal.alarm(0)
    print(json.dumps({k:v for k,v in report.items() if 'sha256' not in k},indent=2))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output.resolve())

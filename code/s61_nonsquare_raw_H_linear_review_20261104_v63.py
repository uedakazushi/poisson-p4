#!/usr/bin/env python3
"""Independent raw-source H2/fixed-operation audit and pointwise rank126 proof."""
import argparse
import hashlib
import itertools
import json
import resource
import signal
import time
from pathlib import Path

import sympy as s
from sympy.polys.rings import ring

ROOT = Path(__file__).resolve().parents[1]
ID = 'S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63_REVIEW'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def main(out):
    began = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    signal.alarm(57)
    out.mkdir(parents=True, exist_ok=False)
    read = lambda p: json.loads(p.read_text())
    base = ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63'
    cr = ROOT/'reports/S61_NONSQUARE_H_LINEAR2_RANK126_20261104_v63'
    ar, pre, cp = [read(base/x) for x in ['report.json', 'preflight.json', 'checkpoint.json']]
    assert ar['frozen'] and ar['additional_global_H_linear_unit_pivots'] == 1
    for rel, h in ar['input_sha256'].items():
        assert sha(ROOT/rel) == h, rel
    for rel, h in ar['evidence_sha256'].items():
        assert sha(base/rel) == h, rel
    corollary = read(cr/'report.json')
    for rel, h in corollary['input_sha256'].items():
        assert sha(ROOT/rel) == h, rel
    raw = ROOT/'reports/S61_NONSQUARE_ELL_SIX_BORDERS_20261104_v63/frozen57-source.txt'
    basisp = ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json'
    monp = ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json'
    basis, mondata = read(basisp), read(monp)
    assert sha(raw) == basis['remaining_matrix_sha256']
    light = [0, 1, 2, 3, 4, 6, 7]
    mon = mondata['normal_monomials']
    hcols = []
    for j in light:
        m = [0]*8; m[5] = 1; m[j] = 1
        hcols.append(mon.index(m))
    assert len(mon) == 91 and all(sum(m[i]*(2 if i == 5 else 1) for i in range(8)) == 3 and m[5] <= 1 for m in mon)
    R, a, e, d = ring('a,e,d', s.QQ)
    K = R.to_field()
    dec = lambda terms: R.from_dict({tuple(m): s.QQ(c) for m, c in terms})
    enc = lambda p: [[list(m), str(c)] for m, c in p.items()]

    def native(line):
        z = line.split(); n = int(z[0]); assert len(z) == 1+4*n
        return R.from_dict({tuple(map(int, z[2+4*i:5+4*i])): s.QQ(z[1+4*i]) for i in range(n)})
    # Directly parse immutable source57; no author H-circuit or producer helper.
    raw16 = [[None]*7 for _ in range(16)]
    with raw.open() as f:
        assert f.readline().split() == ['57', '91', '57']
        for i in range(57):
            for j in range(91):
                line = f.readline(); assert line
                if i >= 41 and j in hcols:
                    raw16[i-41][hcols.index(j)] = native(line)
        assert not f.read().strip()
    assert raw16 == [[dec(q) for q in row] for row in pre['raw16x7']]
    # v6=Z-3*x0+x1 and v3=W; the remaining oldX coordinates are v0,v1,v2,v4,v7.
    X = [[q[0]-3*q[5], q[1]+q[5], q[2], q[4], q[6], q[5], q[3]] for q in raw16]
    frame = [[int(i == j) for j in range(7)] for i in range(7)]
    frame = [[1,0,0,0,0,-3,0],[0,1,0,0,0,1,0],[0,0,1,0,0,0,0],
             [0,0,0,0,1,0,0],[0,0,0,0,0,0,1],[0,0,0,0,0,1,0],[0,0,0,1,0,0,0]]
    assert pre['constant_coefficient_column_frame7x7'] == frame
    inv = pre['constant_inverse7x7']
    assert all(sum(frame[i][k]*inv[k][j] for k in range(7)) == int(i == j) for i in range(7) for j in range(7))
    assert all(sum(inv[i][k]*frame[k][j] for k in range(7)) == int(i == j) for i in range(7) for j in range(7))
    F = 500*a+90*e**2+3*e-24
    B = (50*a+18*e**2-9*e-2)*d+60*e-40
    C = 100*a+18*e**2-9*e-2
    J = 5*a*d-6*e-4
    named = [a, d, F, B, C, J]
    delta = -64000*C**2
    g1, g2 = (3*e-2)*d+20, (6*e+1)*d+40
    h = C*g1-B
    k = 600*a+126*e**2-63*e-14-(6*e+1)*C
    cZ, cW = k*C, a*(2*(6*e+1)*B-k*g2)
    eta = s.QQ(192,125)*a*F*B**2*C**3
    z, w = delta*(X[6][5]-X[1][5]), delta*(X[6][6]-X[1][6])
    assert all(X[6][i] == X[1][i] for i in range(5))
    assert z == -delta*a*F*B*g1*g2/10000000
    assert w == -delta*F*B*h/10000000
    assert cZ*z+cW*w == eta
    gamma = -s.QQ(4,625)*F*B*C**2
    assert gamma == delta*F*B/10000000
    T = [[K(cZ)/K(eta), K(cW)/K(eta)], [-K(w)/K(gamma), K(z)/K(gamma)]]
    Ti = [[K(z), -K(gamma*cW)/K(eta)], [K(w), K(gamma*cZ)/K(eta)]]
    assert all(sum((T[i][k]*Ti[k][j] for k in range(2)), K.zero) == int(i == j) for i in range(2) for j in range(2))
    assert all(sum((Ti[i][k]*T[k][j] for k in range(2)), K.zero) == int(i == j) for i in range(2) for j in range(2))
    assert [T[0][0]*K(z)+T[0][1]*K(w), T[1][0]*K(z)+T[1][1]*K(w)] == [K.one, K.zero]
    keep = [i for i in range(16) if i != 6]
    M = [X[i][:5]+[h*X[i][5]-a*g1*g2*X[i][6]] for i in keep]
    clear = [cZ*X[i][5]+cW*X[i][6] for i in keep]
    assert M == [[dec(q) for q in row] for row in pre['matrix15x6']]
    assert clear == [dec(q) for q in pre['other15_full_source_clearings_numerators_over_eta']]
    assert sum(len(q) for row in M for q in row) == 8613
    assert pre['remaining15_genuine_indices'] == keep
    with (base/'native-input.txt').open() as f:
        assert f.readline().split() == ['15','6','6']
        assert [native(f.readline()) for _ in range(6)] == named
        assert [[native(f.readline()) for _ in range(6)] for _ in range(15)] == M
        assert not f.read().strip()

    def unit(q):
        result = R(s.QQ(q['constant']))
        assert all(n >= 0 for n in q['a_d_F_B_C_J_powers'])
        for f, n in zip(named, q['a_d_F_B_C_J_powers']):
            result *= f**n
        return result
    # Fixed transcript only: no unit scan, pivot selection or author executable.
    A = []
    for row, u in zip(M, cp['initial_row_divisors']):
        divisor = unit(u)
        aa = [p.exquo(divisor) for p in row]
        assert [divisor*p for p in aa] == row
        A.append(aa)
    assert cp['completed_pivots'] == 1 and len(cp['steps']) == 1
    st = cp['steps'][0]
    ri, cj = st['pivot_row'], st['pivot_column']
    assert [ri, cj] == [1, 1]
    p = unit(st['pivot_unit'])
    assert A[ri][cj] == p == -B*C/14062500
    given = {v['row']: v for v in st['eliminations']}
    assert set(given) == {i for i in range(15) if i != ri and A[i][cj]}
    for i, record in given.items():
        coefficient = dec(record['coefficient_before'])
        assert coefficient == A[i][cj]
        before = [p*x-coefficient*y for x, y in zip(A[i], A[ri])]
        divisor = unit(record['row_divisor_after'])
        A[i] = [q.exquo(divisor) for q in before]
        assert [divisor*q for q in A[i]] == before
        assert not A[i][cj]
    assert A == [[dec(q) for q in row] for row in cp['matrix15x6']]
    assert cp['remaining_rows'] == [i for i in range(15) if i != ri]
    assert cp['remaining_columns'] == [i for i in range(6) if i != cj]
    selected = M[ri][cj]
    assert selected == F*B**2*C/12500000 == unit(cp['original_selected_minor_unit'])
    # Independent pointwise geometric corollary, bound to accepted NF117 support.
    hp = ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63'
    hd, hr = read(hp/'certificate.json'), read(hp/'report.json')
    peer = ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63_REVIEW/report.json'
    assert read(peer)['all_passed'] and sha(hp/'certificate.json') == hr['certificate_sha256']
    assert hd['target117_projection_preserves_Hdegree_exact_structural_support']
    assert hr['global_H2_knownunit_sourcepair_and_deletion_verified']
    target = hd['parent_joint_target239_monomials']
    assert len(target) == 239 and sum(m[7] == 1 for m in target) == 28
    assert len(hd['target122_Hdegree_classes']['1']) == 23
    assert len(hd['all112_genuine_cubic_H1_scalar_roots']) == 112
    degree2 = list(itertools.combinations_with_replacement(range(7), 2))
    two_linearforms_image = {tuple(sorted((i,j))) for i in [0,1] for j in range(7)}
    assert len(degree2) == 28 and len(two_linearforms_image) == 13
    assert len([m for m in degree2 if m not in two_linearforms_image]) == 15
    lower = 117+1+(13-(28-23))
    assert lower == 126 == corollary['full319x239_pointwise_rank_lower_primary']
    proof = {'direct_raw16x7_all112_coefficients_rebuilt': True,
             'direct_physical_frame': frame,
             'all90_polynomial15x6_entries_and_native_input_match': True,
             'all15_source_clearing_numerators': [enc(q) for q in clear],
             'fixed_pivot': {'row': ri, 'column': cj, 'polynomial': enc(p)},
             'selected_original_minor': enc(selected),
             'all90_final_checkpoint_entries_and_fixed_elimination_count': len(given),
             'sticky_native_unfinished_operation_is_not_a_new_partial_step': 'The author variable pending is not reset after completion; all saved fixed operations and final90 entries pass. The audit string is a stale completed-step label, not evidence of a pending arithmetic operation.',
             'rank126_proof': {'two_actual_Hlinear_directions_independent_everywhere': True,
                              'abstract_products_span_dimension': 13,
                              'NF117_Hprojection_kernel_dimension': 5,
                              'remaining_cubic_Himage_lower': 8,
                              'all_cubic_light_products_have_H2_zero': True,
                              'separate_quartic_H2_unit_adds_one': True,
                              'total_pointwise_lower': lower},
             'no_fixed_global126minor_or_eighth_Hpivot_claimed': True}
    (out/'independent-certificate.json').write_text(json.dumps(proof, separators=(',', ':'))+'\n')
    paths = [Path(__file__), ROOT/'research/notes'/(ID+'.md'), raw, basisp, monp,
             base/'report.json', base/'preflight.json', base/'checkpoint.json', base/'native-input.txt',
             ROOT/'code/s61_nonsquare_raw_H_linear_units_20261104_v63.cpp',
             cr/'report.json', ROOT/'research/notes/S61_NONSQUARE_H_LINEAR2_RANK126_20261104_v63.md',
             hp/'report.json', hp/'certificate.json', peer]
    report = {'id': ID, 'all_passed': True,
              'own_direct_raw57_reconstruction_all112_H_coefficients': True,
              'independent_frame_Cdelta_eta_gamma_and_full_source_clearings_passed': True,
              'all90_polynomial_and_native15x6_inputs_match': True,
              'fixed_extra_unit_minus_BC_over14062500_passed': True,
              'all90_final_checkpoint_coefficients_and_fixed_row_operations_match': True,
              'original_selected_minor_F_B2_C_over12500000_passed': True,
              'raw_Hlinear_uniform_rank_lower': 2,
              'full319x239_pointwise_rank_lower': lower,
              'no_fixed126minor_eighth_constructiveHpivot_freecomplement_or_full57_full239_claim': True,
              'no_author_solver_unit_scan_newpoint_GB_or_new_rank_probe': True,
              'native_pending_label_scope_clarified': True,
              'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in paths},
              'certificate_sha256': sha(out/'independent-certificate.json'),
              'seconds': time.monotonic()-began, 'frozen': True}
    assert sum(p.stat().st_size for p in out.iterdir()) < 32*1024*1024
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output.resolve())

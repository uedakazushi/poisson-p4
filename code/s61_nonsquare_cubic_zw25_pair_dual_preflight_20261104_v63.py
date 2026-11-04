#!/usr/bin/env python3
"""Fixed six-dot-product pair-dual preflight; no Gaussian or minor search."""
import argparse
import hashlib
import json
import resource
import signal
import time
from pathlib import Path

import sympy as s
from sympy.polys.rings import ring

ROOT = Path(__file__).resolve().parents[1]
ID = 'S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def main(out):
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    signal.alarm(57)
    out.mkdir(parents=True, exist_ok=False)
    paths = {
        'checkpoint': ROOT / 'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json',
        'unit_report': ROOT / 'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/report.json',
        'input': ROOT / 'reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63/input.json',
        'dual': ROOT / 'reports/S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63/global-P22.json',
        'dual_review': ROOT / 'reports/S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63_REVIEW/report.json',
        'Lframe_pair': ROOT / 'reports/S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63/pair2x6.txt',
        'Lframe_review': ROOT / 'reports/S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63_REVIEW/report.json',
        'source': Path(__file__),
    }
    read = lambda p: json.loads(p.read_text())
    st, ur, inp, dual = [read(paths[k]) for k in ['checkpoint', 'unit_report', 'input', 'dual']]
    assert st['completed_pivots'] == 8 and ur['frozen']
    assert not ur['unfinished_step_not_in_checkpoint']
    assert ur['stop_reason'] == 'no_individual_named_FULLU_unit_entry'
    assert sha(paths['checkpoint']) == ur['evidence_sha256']['checkpoint.json']
    assert sha(paths['input']) == st['source_input_sha256']
    assert read(paths['dual_review'])['all_passed']
    assert read(paths['Lframe_review'])['all_passed']
    R, a, e, d = ring('a,e,d', s.QQ)
    K = R.to_field()
    decp = lambda ts: R.from_dict({tuple(m): s.QQ(c) for m, c in ts})
    dec = lambda q: K(decp(q[0])) / K(decp(q[1]))
    encp = lambda p: [[list(m), str(c)] for m, c in p.items()]
    enc = lambda q: [encp(q.numer), encp(q.denom)]
    F = 500*a + 90*e**2 + 3*e - 24
    B = (50*a + 18*e**2 - 9*e - 2)*d + 60*e - 40
    C = 100*a + 18*e**2 - 9*e - 2
    J = 5*a*d - 6*e - 4
    known = [a, d, F, B, C, J]

    def unit(q):
        if not q:
            return None
        exps, constants = [], []
        for p in [q.numer, q.denom]:
            ee = []
            for f in known:
                n = 0
                while not p.rem(f):
                    p = p.exquo(f)
                    n += 1
                ee.append(n)
            if not p.is_ground:
                return None
            exps.append(ee)
            constants.append(p)
        return {'constant': str(K(constants[0])/K(constants[1])),
                'a_d_F_B_C_J_powers': [x-y for x, y in zip(*exps)]}

    # Only ancestors of the twelve accepted pair-dual roots, not all P22.
    cir = dual['circuit']
    roots = dual['pair_dual6x2_scalar_roots']
    need = set()

    def visit(i):
        if i in need:
            return
        need.add(i)
        node = cir['scalar_nodes'][i]
        if node[0] in ('a', 's', 'm'):
            visit(node[1]); visit(node[2])
        else:
            assert node[0] in ('z', 'c')
    for row in roots:
        for i in row:
            visit(i)
    vals = {}
    for i in sorted(need):
        n = cir['scalar_nodes'][i]
        if n[0] == 'z':
            vals[i] = K.zero
        elif n[0] == 'c':
            vals[i] = dec(cir['atoms_num_den'][n[1]])
        elif n[0] == 'a':
            vals[i] = vals[n[1]] + vals[n[2]]
        elif n[0] == 's':
            vals[i] = vals[n[1]] - vals[n[2]]
        else:
            vals[i] = vals[n[1]] * vals[n[2]]
    Dactual = [[vals[i] for i in row] for row in roots]
    # Old recursive-ascending pair order reverses the L-frame descending order.
    # Old Dactual dualizes unscaled q5; the L-frame source uses 1600*q5.
    DL = [[Dactual[5-i][0], Dactual[5-i][1]/1600] for i in range(6)]
    lines = paths['Lframe_pair'].read_text().splitlines()
    assert lines[0] == '2 6'

    def native_poly(line):
        ts = line.split(); n = int(ts[0]); assert len(ts) == 1+4*n
        return K(R.from_dict({tuple(map(int, ts[2+4*k:5+4*k])): s.QQ(ts[1+4*k]) for k in range(n)}))
    pair = [[native_poly(lines[1+6*i+j]) for j in range(6)] for i in range(2)]
    pair_checks = [[sum((pair[i][k]*DL[k][j] for k in range(6)), K.zero)
                    for j in range(2)] for i in range(2)]
    assert pair_checks == [[K.one, K.zero], [K.zero, K.one]]
    # Bind the old pair to the literal initial 25x32 source columns.
    initial = [[K(decp(q)) for q in row] for row in inp['matrix25x32']]
    assert [[initial[i][j] for j in range(20, 26)] for i in [8, 9]] == pair
    assert st['selected_minor_rows'] == [8, 1, 2, 3, 0, 6, 7, 11]
    assert st['selected_minor_columns'] == [20, 10, 16, 12, 0, 23, 24, 26]
    row_ids, col_ids = [4, 5, 9], [21, 22, 25]
    assert all(i in st['remaining_rows'] for i in row_ids)
    assert all(j in st['remaining_columns'] for j in col_ids)
    Qsmall = [[dec(st['matrix25x32'][i][j]) for j in col_ids] for i in row_ids]
    Dsub = [DL[j-20] for j in col_ids]
    six = [[sum((Qsmall[i][k]*Dsub[k][j] for k in range(3)), K.zero)
            for j in range(2)] for i in range(3)]
    outcomes = []
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        det = six[i][0]*six[j][1] - six[i][1]*six[j][0]
        outcomes.append({'residual_row_pair': [row_ids[i], row_ids[j]],
                         'determinant': enc(det), 'named_unit': unit(det)})
    success = [v for v in outcomes if v['named_unit'] is not None]
    # Fixed preflight only: a positive two-minor would authorize a later exact
    # GL3 completion. No additional Gaussian or source/target operation here.
    certificate = {
        'pair_order_reverse': list(reversed(range(6))),
        'Lframe_q5_scale': 1600,
        'evaluated_dual_scalar_ancestors': sorted(need),
        'Lframe_pair2x6': [[enc(q) for q in row] for row in pair],
        'Lframe_dual6x2': [[enc(q) for q in row] for row in DL],
        'all_four_pair_dual_identities': True,
        'initial_pair_source_rows': [8, 9], 'initial_pair_columns': list(range(20, 26)),
        'residual_row_ids': row_ids, 'residual_row_labels': [inp['source25_labels'][i] for i in row_ids],
        'surviving_original_pair_columns': col_ids,
        'current_residual_column_positions': [st['remaining_columns'].index(j) for j in col_ids],
        'Qsmall3x3': [[enc(q) for q in row] for row in Qsmall],
        'Dsub3x2': [[enc(q) for q in row] for row in Dsub],
        'six_dot_products3x2': [[enc(q) for q in row] for row in six],
        'three_fixed_minor_outcomes': outcomes,
        'all_three_minors_zero': all(not dec(v['determinant']) for v in outcomes),
        'known_unit_two_minor_found': bool(success),
        'GL3_completion_criterion_only': 'If a selected F_I=Qsmall_I*Dsub is a unit 2x2 matrix, put L=F_I^-1*Qsmall_I. Then L*Dsub=I2 and [Dsub, cross(L_row0,L_row1)] has determinant one. No such completion is claimed without this gate.',
    }
    (out/'certificate.json').write_text(json.dumps(certificate, separators=(',', ':'))+'\n')
    report = {
        'id': ID, 'all_source_order_scaling_and_four_pair_dual_identities_passed': True,
        'only_six_fixed_dots_and_three_fixed_minors_evaluated': True,
        'known_unit_two_minor_found': bool(success),
        'all_three_fixed_minors_zero': certificate['all_three_minors_zero'],
        'six_dot_zero_pattern': [[not bool(q) for q in row] for row in six],
        'no_new_global_split_or_pivot_promoted': True,
        'no_Gaussian_retry_newpoint_GB_unknown_unit_inversion_or_Delta_expansion': True,
        'whole_FULLU_cubic57_quartic239_and_kernel38_missing42_unchanged': True,
        'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in paths.values()},
        'certificate_sha256': sha(out/'certificate.json'),
        'seconds': time.monotonic()-start, 'frozen': True, 'independent_review_pending': True,
    }
    assert sum(p.stat().st_size for p in out.iterdir()) < 32*1024*1024
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output.resolve())

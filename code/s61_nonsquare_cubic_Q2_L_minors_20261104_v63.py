#!/usr/bin/env python3
"""Exactly three fixed actual Q2 L-block minors, stripping only named units."""
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
ID = 'S61_NONSQUARE_CUBIC_Q2_L_MINORS_20261104_v63'


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def main(out):
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    signal.alarm(57)
    out.mkdir(parents=True, exist_ok=False)
    read = lambda p: json.loads(p.read_text())
    old = ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63'
    new = ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_UNIT_COMPLETION_20261104_v63'
    paths = [Path(__file__), old/'checkpoint.json', old/'report.json',
             new/'certificate.json', new/'report.json']
    cp, oldr, tr, newr = [read(p) for p in paths[1:]]
    assert sha(paths[1]) == oldr['evidence_sha256']['checkpoint.json']
    assert sha(paths[3]) == newr['certificate_sha256']
    assert cp['completed_pivots'] == 8 and newr['actual_new_shape_same_cokernel'] == [16, 23]
    rows, cols = [5, 9], [21, 22, 25]
    positions = [[tr['new16_source_original25_rows'].index(i),
                  [tr['new23_target_original32_columns'].index(j) for j in cols]] for i in rows]
    nodes = tr['scalar_DAG']['nodes']
    for i, (ri, cjs) in zip(rows, positions):
        for j, cj in zip(cols, cjs):
            root = tr['remaining16x23_scalar_roots'][ri][cj]
            assert nodes[root] == ['checkpoint', i, j]
    R, a, e, d = ring('a,e,d', s.QQ)
    K = R.to_field()
    decp = lambda ts: R.from_dict({tuple(m): s.QQ(c) for m, c in ts})
    dec = lambda q: K(decp(q[0]))/K(decp(q[1]))
    encp = lambda p: [[list(m), str(c)] for m, c in p.items()]
    enc = lambda q: [encp(q.numer), encp(q.denom)]
    M = [[dec(cp['matrix25x32'][i][j]) for j in cols] for i in rows]
    profile = lambda q: {'numerator_terms': len(q.numer), 'denominator_terms': len(q.denom),
                         'numerator_degree': max([sum(m) for m in q.numer]+[0]),
                         'denominator_degree': max([sum(m) for m in q.denom]+[0])}
    assert all(q.denom.is_ground for row in M for q in row)
    bounds = []
    pairs = [(0, 1), (0, 2), (1, 2)]
    for i, j in pairs:
        bounds.append(len(M[0][i].numer)*len(M[1][j].numer)+len(M[0][j].numer)*len(M[1][i].numer))
    assert sum(bounds) == 414 and sum(bounds) < 10000
    F = 500*a+90*e**2+3*e-24
    B = (50*a+18*e**2-9*e-2)*d+60*e-40
    C = 100*a+18*e**2-9*e-2
    J = 5*a*d-6*e-4
    known = [a, d, F, B, C, J]
    outcomes = []
    for i, j in pairs:
        det = M[0][i]*M[1][j]-M[0][j]*M[1][i]
        p = det.numer
        powers = []
        if p:
            for f in known:
                n = 0
                while not p.rem(f):
                    p = p.exquo(f); n += 1
                powers.append(n)
            restore = K(p)/K(det.denom)
            for f, n in zip(known, powers):
                restore *= K(f)**n
            assert restore == det
        else:
            powers = [0]*6
        outcomes.append({'original_target_column_pair': [cols[i], cols[j]],
                         'minor': enc(det), 'profile': profile(det),
                         'named_unit_numerator_factor_powers_a_d_F_B_C_J': powers,
                         'untouched_unknown_residual_numerator': encp(p),
                         'nonzero_QQ_denominator': str(det.denom),
                         'named_unit_verdict': bool(p) and p.is_ground,
                         'raw_product_term_bound': bounds[len(outcomes)]})
    certificate = {'source_old25_rows': rows, 'source_new16_positions': [v[0] for v in positions],
                   'target_old32_columns': cols, 'target_new23_positions': positions[0][1],
                   'all_six_new_presentation_roots_literal_old_checkpoint_references': True,
                   'actual2x3': [[enc(q) for q in row] for row in M],
                   'actual2x3_profiles': [[profile(q) for q in row] for row in M],
                   'exact_three_minor_outcomes': outcomes,
                   'only_named_factors_divided': ['a', 'd', 'F', 'B', 'C', 'J'],
                   'all_factor_reconstruction_identities_exact': True,
                   'no_unknown_factor_factored_inverted_or_promoted': True}
    (out/'certificate.json').write_text(json.dumps(certificate, separators=(',', ':'))+'\n')
    report = {'id': ID, 'exactly_three_fixed_actual2x3_minors_verified': True,
              'literal_final16x23_coefficient_binding_passed': True,
              'known_unit_minor_count': sum(v['named_unit_verdict'] for v in outcomes),
              'all_minors_zero': all(not v['minor'][0] for v in outcomes),
              'total_raw_product_term_bound': sum(bounds),
              'minor_profiles': [v['profile'] for v in outcomes],
              'unit_factor_exponents': [v['named_unit_numerator_factor_powers_a_d_F_B_C_J'] for v in outcomes],
              'no_Gaussian_transport_point_GB_Delta_inversion_or_automatic_enlargement': True,
              'full_FULLU57_full239_and_kernel38_gap42_unchanged': True,
              'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in paths},
              'certificate_sha256': sha(out/'certificate.json'),
              'seconds': time.monotonic()-start, 'frozen': True, 'independent_review_pending': True}
    assert sum(p.stat().st_size for p in out.iterdir()) < 32*1024*1024
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output.resolve())

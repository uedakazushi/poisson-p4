#!/usr/bin/env python3
"""Read-only support/cost preflight for actual saved residual Q3; no minors."""
import argparse
import hashlib
import itertools
import json
import resource
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ID = 'S61_NONSQUARE_CUBIC_ZW25_Q3_SUPPORT_PREFLIGHT_20261104_v63'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def main(out):
    began = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    out.mkdir(parents=True, exist_ok=False)
    parent = ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63'
    old = ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63'
    read = lambda p: json.loads(p.read_text())
    cp, rp = parent/'checkpoint.json', parent/'report.json'
    state, prior = read(cp), read(rp)
    assert state['completed_pivots'] == 8 and prior['frozen']
    assert sha(cp) == prior['evidence_sha256']['checkpoint.json']
    pc = read(old/'certificate.json')
    assert sha(old/'certificate.json') == read(old/'report.json')['certificate_sha256']
    rows = [4, 5, 9]
    cols = state['remaining_columns']
    mons = [tuple(t.count(i) for i in range(5))
            for t in itertools.combinations_with_replacement(range(5), 3)]
    keep = list(range(15, 35))+list(range(9, 15))
    labels = [mons[keep[j]] if j < 26 else None for j in cols]

    def profile(q):
        n, d = q
        return {'nonzero': bool(n), 'numerator_terms': len(n),
                'denominator_terms': len(d),
                'numerator_degree': max([sum(m) for m, _ in n]+[0]),
                'denominator_degree': max([sum(m) for m, _ in d]+[0])}

    matrix = [[state['matrix25x32'][i][j] for j in cols] for i in rows]
    pp = [[profile(q) for q in row] for row in matrix]
    assert all(p['denominator_degree'] == 0 for row in pp for p in row)
    assert all(not pp[i][k]['nonzero'] for i in range(3)
               for k, j in enumerate(cols) if j >= 26)
    blocks = {}
    for name, allowed in [('noL_pure', [j for j in cols if j < 20]),
                          ('L_pure', [j for j in cols if 20 <= j < 26]),
                          ('heavy_and_psi', [j for j in cols if j >= 26])]:
        blocks[name] = {'original_columns': allowed,
                        'row_nonzero_counts': [sum(pp[i][cols.index(j)]['nonzero'] for j in allowed) for i in range(3)],
                        'row_numerator_terms': [sum(pp[i][cols.index(j)]['numerator_terms'] for j in allowed) for i in range(3)]}
    active = [j for j in cols if any(pp[i][cols.index(j)]['nonzero'] for i in range(3))]
    candidates = []
    for cc in itertools.combinations(active, 3):
        # Counts only: no coefficients are multiplied and no determinant is evaluated.
        terms = []
        for perm in itertools.permutations(range(3)):
            ps = [pp[i][cols.index(cc[perm[i]])] for i in range(3)]
            if all(p['nonzero'] for p in ps):
                terms.append({'raw_product_term_bound': ps[0]['numerator_terms']*ps[1]['numerator_terms']*ps[2]['numerator_terms'],
                              'degree_bound': sum(p['numerator_degree'] for p in ps)})
        if terms:
            candidates.append({'original_columns': list(cc),
                               'L_column_count': sum(20 <= j < 26 for j in cc),
                               'potential_nonzero_permutation_products': len(terms),
                               'raw_determinant_term_bound': sum(t['raw_product_term_bound'] for t in terms),
                               'degree_bound': max(t['degree_bound'] for t in terms)})
    candidates.sort(key=lambda x: (x['raw_determinant_term_bound'], x['degree_bound'], x['original_columns']))
    # The row4*x0 pair-dual transport is a proposed exact identity, NOT computed here.
    D = pc['Lframe_dual6x2']
    assert D[1][1][0] == []
    eliminations4 = [[v['pivot_row'], v['pivot_column']]
                     for v in state['all_completed_unit_steps']
                     if any(t == 4 for t, _ in v['eliminations'])]
    assert eliminations4 == [[0, 0]]
    cert = {'source_rows': rows, 'source_labels': ['scaled_q5*x0', 'scaled_q5*y1', 'scaled_q5*L'],
            'original_target_columns': cols, 'target_monomials_L_x0_y1_y2_y3_or_null': labels,
            'matrix3x24': matrix, 'entry_profiles': pp, 'support_blocks': blocks,
            'active_columns': active, 'active_column_count': len(active),
            'pattern_viable_three_minor_candidate_count': len(candidates),
            'lowest_cost_12_fixed_candidates_NOT_EVALUATED': candidates[:12],
            'all_candidate_profile_only': {'maximum_term_bound': max(c['raw_determinant_term_bound'] for c in candidates),
                                           'maximum_degree_bound': max(c['degree_bound'] for c in candidates)},
            'row4_modified_only_by_q0_x0_step': eliminations4,
            'proposed_row4_pair_dual_column_target_columns': list(range(4, 10)),
            'proposed_row4_pair_dual_column6': [r[1] for r in D],
            'dual_zero_component_local_index': 1,
            'fixed_six_dot_identity_not_evaluated_in_this_packet': True,
            'GL6_completion_gate': 'If a=(abar,a_k), v=(vbar,0) in the indicated six target coordinates and abar*vbar=1, U0=[[I-vbar*abar,vbar],[abar,0]] is an involution. With E=[[I,0],[-a_k*abar,1]], a*U0*E=(0,1) and det(U0*E)=-1. No unknown factor is inverted. The exact dot and finite operator are not yet verified here.'}
    (out/'certificate.json').write_text(json.dumps(cert, separators=(',', ':'))+'\n')
    sources = [Path(__file__), cp, rp, old/'report.json', old/'certificate.json',
               ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63_REVIEW/report.json']
    report = {'id': ID, 'actual_saved_Q3_support_bound': True,
              'support_blocks': blocks, 'active_columns': active,
              'lowest_cost_fixed_candidate_profiles': candidates[:4],
              'remaining_Q3_rank3_ALL_FULLU_uses_accepted_Q19_and_saved_7_Q_pivots': True,
              'no_determinant_dot_product_Gaussian_newpoint_or_GB_evaluated': True,
              'no_new_global_pivot_dual_completion_or_full57_claim': True,
              'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in sources},
              'certificate_sha256': sha(out/'certificate.json'),
              'seconds': time.monotonic()-began, 'frozen': True}
    assert sum(p.stat().st_size for p in out.iterdir()) < 32*1024*1024
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    main(ap.parse_args().output.resolve())

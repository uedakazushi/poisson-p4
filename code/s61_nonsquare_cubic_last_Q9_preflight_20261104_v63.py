#!/usr/bin/env python3
"""Exactly nine small actual last-Q coefficients; named-factor strips only."""
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
ID = 'S61_NONSQUARE_CUBIC_LAST_Q9_PREFLIGHT_20261104_v63'


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def main(out):
    began = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    signal.alarm(57)
    out.mkdir(parents=True, exist_ok=False)
    read = lambda p: json.loads(p.read_text())
    q = ROOT/'reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63'
    p = ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_UNIT_COMPLETION_20261104_v63'
    cpfile = ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json'
    paths = [Path(__file__), q/'report.json', q/'certificate.json', p/'report.json', p/'certificate.json', cpfile]
    qr, qc, pr, pc, cp = [read(x) for x in paths[1:]]
    assert sha(paths[2]) == qr['certificate_sha256'] and sha(paths[4]) == pr['certificate_sha256']
    assert qc['new15_source_original25_labels'][0] == 9
    assert qr['actual_new_shape_same_cokernel'] == [15,22]
    pure = [11,13,14,15,17,18,19]
    selected = pure+[22,25]
    Lcols = [21,22,25]
    R, a, e, d = ring('a,e,d', s.QQ)
    K = R.to_field()
    decp = lambda ts: R.from_dict({tuple(m): s.QQ(c) for m,c in ts})
    dec = lambda x: K(decp(x[0]))/K(decp(x[1]))
    encp = lambda v: [[list(m),str(c)] for m,c in v.items()]
    enc = lambda v: [encp(v.numer),encp(v.denom)]
    datum = lambda i,j: dec(cp['matrix25x32'][i][j])
    b = [datum(5,j) for j in Lcols]
    sv = [datum(9,j) for j in Lcols]
    v = [K(6*(18*e+7))/K(5*a), K(-18), K(-9)]
    cbound = 2*len(sv[0].numer)+len(sv[1].numer)+len(sv[2].numer)
    assert cbound == 52
    bounds = [len(datum(9,j).numer)+cbound*len(datum(5,j).numer) for j in pure]
    bounds += [len(sv[i].numer)+cbound*len(b[i].numer) for i in [0,1]]
    assert sum(bounds) == 5890 and sum(bounds) < 6000
    c = sum((sv[i]*v[i] for i in range(3)), K.zero)
    expected = [datum(9,j)-c*datum(5,j) for j in pure]+[sv[i]-c*b[i] for i in [0,1]]
    # Evaluate ONLY ancestors of the nine selected actual15x22 roots.
    nodes = qc['scalar_DAG']['nodes']
    roots = [qc['remaining15x22_scalar_roots'][0][qc['new22_target_original32_columns'].index(j)] for j in selected]
    needed = set()
    def visit(n):
        if n in needed:
            return
        needed.add(n)
        node = nodes[n]
        if node[0] in ['a','s','m']:
            visit(node[1]);visit(node[2])
        else:
            assert node[0] in ['z','one','checkpoint','new_unit_atom']
    for n in roots:
        visit(n)
    allowed = {(i,j) for i in [5,9] for j in pure+Lcols}
    values, refs = {}, set()
    for n in sorted(needed):
        node = nodes[n]
        if node[0] == 'z':
            value = K.zero
        elif node[0] == 'one':
            value = K.one
        elif node[0] == 'checkpoint':
            ref = tuple(node[1:]);assert ref in allowed
            refs.add(ref);value = datum(*ref)
        elif node[0] == 'new_unit_atom':
            value = dec(qc['scalar_DAG']['new_unit_atoms_num_den'][node[1]])
        elif node[0] == 'a':
            value = values[node[1]]+values[node[2]]
        elif node[0] == 's':
            value = values[node[1]]-values[node[2]]
        else:
            value = values[node[1]]*values[node[2]]
        values[n] = value
    assert [values[n] for n in roots] == expected
    F = 500*a+90*e**2+3*e-24
    B = (50*a+18*e**2-9*e-2)*d+60*e-40
    C = 100*a+18*e**2-9*e-2
    J = 5*a*d-6*e-4
    known = [a,d,F,B,C,J]
    outcomes = []
    for j, value in zip(selected, expected):
        stripped, exponents = [], []
        for pol in [value.numer,value.denom]:
            ee = []
            for f in known:
                nn = 0
                while pol and not pol.rem(f):
                    pol = pol.exquo(f);nn += 1
                ee.append(nn)
            stripped.append(pol);exponents.append(ee)
        assert stripped[1].is_ground and stripped[1]
        powers = [x-y for x,y in zip(*exponents)]
        restore = K(stripped[0])/K(stripped[1])
        for f,nn in zip(known,powers):
            restore *= K(f)**nn
        assert restore == value
        outcomes.append({'original_target_column':j,'actual_coefficient':enc(value),
                         'named_factor_powers_a_d_F_B_C_J':powers,
                         'unknown_numerator_left_unfactored':encp(stripped[0]),
                         'remaining_QQ_denominator':str(stripped[1]),
                         'named_unit_verdict':bool(stripped[0]) and stripped[0].is_ground,
                         'numerator_terms':len(value.numer),
                         'numerator_degree':max([sum(m) for m in value.numer]+[0])})
    # Full source coefficient formula, left lazy in the two accepted circuits:
    # row9 - c*row5 + (c*d5-d9)*row4 in the eight-step checkpoint source frame.
    row4 = pc['old17_source_original25_rows'].index(4)
    row5 = pc['old17_source_original25_rows'].index(5)
    row9 = pc['old17_source_original25_rows'].index(9)
    negd5 = pc['source17_roots'][row5][row4]
    negd9 = pc['source17_roots'][row9][row4]
    negc = qc['source16_roots'][qc['old16_source_original25_rows'].index(9)][qc['old16_source_original25_rows'].index(5)]
    cert = {'actual_row_original25_label':9,
            'query_only_original32_target_columns':selected,
            'old_checkpoint_b_row5_L3':[enc(x) for x in b],
            'old_checkpoint_s_row9_L3':[enc(x) for x in sv],
            'literal_dual_v':[enc(x) for x in v], 'c_s_dot_v':enc(c),
            'actual_nine_coefficients':outcomes,
            'raw_product_term_bounds':bounds,'raw_total_bound':sum(bounds),
            'all9_actual_scalar_DAG_roots_match_direct_formula':True,
            'queried_roots':roots,'evaluated_ancestor_count':len(needed),
            'evaluated_checkpoint_references':sorted([list(x) for x in refs]),
            'full_source_cofactor_recipe_in_eight_step_checkpoint_frame':{
                'row9_coefficient':'1','row5_coefficient':{'actual_expanded':enc(-c),'Q2_GL3_source16_root':negc},
                'row4_coefficient':{'formula':'negd9+negc*negd5','parent_GL6_negd9_root':negd9,'parent_GL6_negd5_root':negd5},
                'other_checkpoint_source_rows_coefficient':'0',
                'reference':'Full source relation is this linear combination followed by the two saved invertible target changes. Pure/heavy tails are retained. Parent scalar DAGs are pinned; the row4 coefficient is deliberately not expanded.'},
            'no_other_last_row_coefficient_expanded':True,
            'only_named_factor_strips_no_unknown_factorization_or_inversion':True}
    (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
    report = {'id':ID,'all9_actual_DAG_coefficients_match_source_formula':True,
              'raw_product_bound':sum(bounds),'known_unit_coefficient_count':sum(v['named_unit_verdict'] for v in outcomes),
              'actual_total_numerator_terms':sum(v['numerator_terms'] for v in outcomes),
              'actual_maximum_numerator_degree':max(v['numerator_degree'] for v in outcomes),
              'actual_coefficient_profiles':[{k:v[k] for k in ['original_target_column','numerator_terms','numerator_degree','named_factor_powers_a_d_F_B_C_J','named_unit_verdict']} for v in outcomes],
              'no_Gaussian_minor_point_GB_or_automatic_extra_query':True,
              'global_primitive_dual_full57_full239_not_claimed':True,
              'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
              'certificate_sha256':sha(out/'certificate.json'),
              'seconds':time.monotonic()-began,'frozen':True,'independent_review_pending':True}
    assert sum(p.stat().st_size for p in out.iterdir()) < 32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps(report,indent=2),flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output.resolve())

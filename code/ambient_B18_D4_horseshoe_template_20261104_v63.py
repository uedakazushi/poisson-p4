#!/usr/bin/env python3
"""Functorial horseshoe builder for the actual N cone and cyclic nilcurve.

build_horseshoe is an arithmetic API, with actual coefficient nodes supplied by
the separate delta extractor.  The command checks the universal identities;
it does not replace the pending actual coefficient binding.
"""
import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

from ambient_B18_T3_syzygy_bockstein_preflight_20261104_v63 import add, mul, scale, var, ZERO, ONE

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_D4_HORSESHOE_TEMPLATE_20261104_v63'


def decode(col):
    return {int(i): {tuple(mon): int(c) for c, mon in enc} for i, enc in col}


def encode(col):
    return [[i, [[c, list(mon)] for mon, c in sorted(p.items())]] for i, p in sorted(col.items()) if p]


def put(col, i, p):
    col[i] = add(col.get(i, ZERO), p)
    if not col[i]: del col[i]


def compose(left, col):
    out = {}
    for j, c in col.items():
        for i, p in left[j].items(): put(out, i, mul(c, p))
    return out


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()


def linear_coefficients_of_d2(ndata):
    ans = []
    for col in ndata['d2_universal_sparse_polynomial_columns']:
        entries = {}
        for alpha, enc in col:
            coeff = {}
            for c, mon in enc:
                assert len(mon) == 1 and mon[0] in list(range(8))+[9]
                coeff[int(mon[0])] = coeff.get(int(mon[0]), 0)+int(c)
            entries[int(alpha)] = coeff
        ans.append(entries)
    assert len(ans) == 84
    return ans


def build_horseshoe(D, ndata, c_nodes36, psi_nodes8):
    """Append literal d1/d2, given actual c(Z),psi(Z) in the SAME scalar DAG.

    c_alpha is POSITIVE: Nrelation_alpha(lifted target generators)=c_alpha*F.
    psi order is z2,...,z8,v, parameter Z=z1.  All functions must already be
    correctly composed with Z=z1 by a unit input, not numeric substitution.
    The returned obstruction nodes must be zero as functions on that curve.
    """
    assert len(c_nodes36) == 36 and len(psi_nodes8) == 8
    z = ndata['z1_through_z8_actual_odd_coefficient_nodes']
    vv = ndata['actual_v_node']
    P = [D.add(a, D.neg(b)) for a, b in zip(z[1:]+[vv], psi_nodes8)]
    n1 = [{int(i): q for i, q in c.items()} for c in ndata['R0_sparse_d1_8_by36_columns']]
    n2 = [{int(i): q for i, q in c.items()} for c in ndata['R0_sparse_d2_36_by84_columns']]
    linear = linear_coefficients_of_d2(ndata)
    d1 = [{0: p} for p in P]
    for alpha, col in enumerate(n1):
        d1.append({0: D.neg(c_nodes36[alpha]), **{i+1: q for i, q in col.items()}})
    d2 = []
    for i, j in itertools.combinations(range(8), 2):
        d2.append({j: P[i], i: D.neg(P[j])})
    obstruction = []
    offdiag = []
    curved = {0: z[0], **{j: psi_nodes8[j-1] for j in range(1, 8)}, 9: psi_nodes8[7]}
    for beta, col in enumerate(n2):
        h = []
        for j in list(range(1, 8))+[9]:
            h.append(D.total(D.mul(D.const(co.get(j, 0)), c_nodes36[alpha])
                             for alpha, co in linear[beta].items() if co.get(j, 0)))
        obs = D.total(D.mul(c_nodes36[alpha], D.total(D.mul(D.const(a), curved[j]) for j, a in co.items()))
                      for alpha, co in linear[beta].items())
        d2.append({**{i: q for i, q in enumerate(h)}, **{8+alpha: q for alpha, q in col.items()}})
        offdiag.append(h)
        obstruction.append(obs)
    assert len(d1) == 44 and len(d2) == 112
    return {'nil_annihilator8_nodes': P, 'D4_d1_9_by44_columns': d1,
            'D4_d2_44_by112_columns': d2, 'D2_offdiagonal8_by84_columns': offdiag,
            'all84_curve_cocycle_obstruction_nodes': obstruction,
            'required_actual_identity': 'obstructions=0 via actual delta extractor and F free basis'}


def main(out):
    start = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    out.mkdir(parents=True, exist_ok=False)
    np = ROOT/'reports/AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63/actual-N-literal-cone-DAG.json'
    nr = np.parent/'report.json'
    ndata = json.loads(np.read_text())
    assert json.loads(nr.read_text())['all_passed']
    n1 = [decode(c) for c in ndata['d1_universal_sparse_polynomial_columns']]
    n2 = [decode(c) for c in ndata['d2_universal_sparse_polynomial_columns']]
    linear = linear_coefficients_of_d2(ndata)
    z = [var(i) for i in range(8)]
    nu, vv = var(8), var(9)
    cc = [var(10+i) for i in range(36)]
    psi = [var(46+i) for i in range(8)]
    P = [add(a, scale(b, -1)) for a, b in zip(z[1:]+[vv], psi)]
    d1 = [{0: p} for p in P]
    for alpha, col in enumerate(n1):
        d1.append({0: scale(cc[alpha], -1), **{i+1: p for i, p in col.items()}})
    d2 = [{j: P[i], i: scale(P[j], -1)} for i, j in itertools.combinations(range(8), 2)]
    obstructions = []
    curved = {0: z[0], **{j: psi[j-1] for j in range(1, 8)}, 9: psi[7]}
    for beta, col in enumerate(n2):
        h = []
        for j in list(range(1, 8))+[9]:
            h.append(add(*(scale(cc[alpha], co[j]) for alpha, co in linear[beta].items() if j in co)))
        obs = add(*(mul(cc[alpha], add(*(scale(curved[j], a) for j, a in co.items())))
                    for alpha, co in linear[beta].items()))
        d2.append({**{i: p for i, p in enumerate(h) if p}, **{8+alpha: p for alpha, p in col.items()}})
        obstructions.append(obs)
    assert len(d1) == 44 and len(d2) == 112
    for col in d2[:28]: assert not compose(d1, col)
    for beta, col in enumerate(d2[28:]):
        got = compose(d1, col)
        expected = {0: scale(obstructions[beta], -1)} if obstructions[beta] else {}
        assert got == expected

    # Matrix functional calculus is finite: T0^3=q*T0 on the ORIGINAL free9.
    T = [{1: ONE}] + [{0: mul(nu, mul(z[0], zi))} for zi in z]
    T2 = [compose(T, c) for c in T]
    T3 = [compose(T, c) for c in T2]
    q = mul(nu, mul(z[0], z[0]))
    assert T3 == [{i: mul(q, p) for i, p in c.items() if mul(q, p)} for c in T]
    # (v-T)(v^2+vT+T^2-q)= (v^3-qv)*I, coefficientwise.
    D = []
    for j in range(9):
        col = {j: add(mul(vv, vv), scale(q, -1))}
        for i, p in T[j].items(): put(col, i, mul(vv, p))
        for i, p in T2[j].items(): put(col, i, p)
        D.append(col)
    vmT = []
    for j in range(9):
        col = {j: vv}
        for i, p in T[j].items(): put(col, i, scale(p, -1))
        vmT.append(col)
    for j, col in enumerate(D):
        assert compose(vmT, col) == {j: add(mul(vv, mul(vv, vv)), scale(mul(q, vv), -1))}

    ev = {
        'variables': ['z1','z2','z3','z4','z5','z6','z7','z8','nu','v']+
                     [f'c{i}' for i in range(36)]+[f'psi{i}' for i in range(8)],
        'sparse_polynomial_encoding': 'coefficient and sorted variable-index multiset',
        'generic_D4_d1_9_by44': [encode(c) for c in d1],
        'generic_D4_d2_44_by112': [encode(c) for c in d2],
        'all84_cocycle_obstructions_on_curve': [[[c, list(m)] for m, c in sorted(p.items())] for p in obstructions],
        'whole_matrix_identity': 'd1*d2=0 on Nil28; on N84 equals only row0 minus curve obstruction',
        'no_nu_in_d2_N_and_no_curve_Hadamard_needed_for_offdiagonal': True,
        'T0_columns': [encode(c) for c in T],
        'T0_square_columns': [encode(c) for c in T2],
        'T0_cubic_identity': 'T0^3=nu*z1^2*T0',
        'matrix_monic_division_identity': '(v-T0)(v^2+vT0+T0^2-nu*z1^2 I)=(v^3-nu*z1^2 v)I',
        'actual_delta_coefficients_bound_in_this_packet': False,
        'actual2069_old_relations_lifted_in_this_packet': False,
    }
    ep = out/'horseshoe-and-lifting-template.json'
    ep.write_text(json.dumps(ev, separators=(',', ':'))+'\n')
    assert ep.stat().st_size < 32*1024*1024
    paths = [Path(__file__), ROOT/'research/notes'/f'{STEM}.md', np, nr,
             ROOT/'code/ambient_B18_T3_syzygy_bockstein_preflight_20261104_v63.py']
    report = {
        'id': STEM, 'all_passed': True,
        'status': 'functorial_literal_horseshoe_and_finite_cone_homotopy_template_actual_coefficients_pending',
        'horseshoe_free_ranks_through2': [9,44,112],
        'nil_Koszul_d2_columns': 28, 'N_d2_columns': 84,
        'full_universal_d1_d2_identities_checked': True,
        'finite_matrix_monic_division_identities_checked': True,
        'actual_delta36_and_nil_curve8_inputs_supplied': False,
        'old_r_g_restricted_or_redivided': False,
        'actual_T3_presentation_or_depth_claimed': False,
        'no_GB_new_jets_points_primes_or_source_replay': True,
        'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in paths},
        'evidence_sha256': {ep.name: sha(ep)},
        'seconds': time.monotonic()-start,
        'max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'certificate_bytes': ep.stat().st_size,
    }
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ['all_passed','horseshoe_free_ranks_through2','seconds','max_rss_kib','certificate_bytes']}))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output)

#!/usr/bin/env python3
"""A homogeneous initial-submodule certificate, not the whole initial module."""
import argparse, hashlib, itertools, json, resource, signal, time
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_T_KNOWN_INITIAL_MODULE_LENGTH_UPPER24_20261104_v63'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p, x): p.write_text(json.dumps(x, separators=(',', ':')) + '\n')
def enc(p): return [[t, list(m), str(c)] for (t, m), c in sorted(p.items())]
def add(p, q, c=F(1), shift=(0,)*8):
    for (t, m), a in q.items():
        key = t, tuple(x+y for x, y in zip(m, shift))
        p[key] = p.get(key, F(0)) + c*a
        if not p[key]: del p[key]
def order(k): return (-k[0],) + k[1]
def leading(p): return max(p, key=order)
def divides(a, b): return a[0] == b[0] and all(x <= y for x, y in zip(a[1], b[1]))
def exponent(indices): return tuple(indices.count(i) for i in range(8))

def reduce(p, gs):
    p = dict(p); remainder = {}; trace = []
    while p:
        k = leading(p)
        for j, g in enumerate(gs):
            l = leading(g)
            if divides(l, k):
                c = p[k]/g[l]; shift = tuple(x-y for x, y in zip(k[1], l[1]))
                add(p, g, -c, shift); trace.append([j, list(shift), str(c)])
                assert len(trace) < 5000 and len(p) < 150000
                break
        else: remainder[k] = p.pop(k)
    return remainder, trace

def main(out):
    start = time.monotonic(); signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    out.mkdir(parents=True, exist_ok=False)
    paths = [
      'reports/AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63/linear-module-rewrite-certificate.json',
      'reports/AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63/actual-known-quadratic-cubic-closure-certificate.json',
      'reports/AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63_REVIEW/report.json',
      'reports/AMBIENT_B18_CONDUCTOR_SQUARE_CLOSED348_20261104_v63_REVIEW/report.json',
      'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW/report.json',
      'reports/AMBIENT_B18_SPECIALIZED_COMPARISON_T_PREFLIGHT_20261104_v63_REVIEW/report.json']
    docs = [json.loads((ROOT/p).read_text()) for p in paths]
    for d in docs[2:]: assert d['all_passed']
    assert docs[3]['actual_base_m10_zero_r9_nonzero_accepted']
    linear = [{(i//8, exponent([i%8])): F(c) for i,c in p} for p in docs[0]['linear_basis37']]
    mons2 = list(itertools.combinations_with_replacement(range(8), 2))
    quad = [{(i//36, exponent(list(mons2[i%36]))): F(c) for i,c in p} for p in docs[1]['quadratic_basis178']]
    assert len(linear) == 37 and len(quad) == 178
    # Only 4 standard degree-two terms remain after the fixed linear rewrite.
    # All source relations and traces stay tiny: this preflight is independent
    # of source-jet orders and of the 348-dimensional scalar matrix.
    bound = sum(map(len, linear+quad)) + 178*200 + 741*100 + 5000
    assert bound < 150000
    dump(out/'preflight.json', dict(input_terms=sum(map(len,linear+quad)),
        conservative_QQ_and_trace_slot_bound=bound, cap=150000,
        no_new_source_jet_or_348_matrix=True))
    extra = []; selected = []; source_reductions = []
    for j, q in enumerate(quad):
        r, tr = reduce(q, linear)
        r0 = dict(r); r, etr = reduce(r, extra) if r else ({}, [])
        if r:
            k = leading(r); c = r[k]
            extra.append({i:a/c for i,a in r.items()})
            selected.append(dict(source_quad_index=j, linear_trace=tr,
                previous_extra_trace=etr, divide=str(c)))
        source_reductions.append(dict(source_quad_index=j, linear_remainder=enc(r0)))
    assert len(extra) == 2
    expected = [
        {(4, exponent([5,5])): F(1)},
        {(3, exponent([7,7])): F(1), (4, exponent([5,7])): F(11767,3652)}]
    assert extra == expected
    gs = linear + extra; pairs = []; trace_terms = 0
    for i in range(len(gs)):
        for j in range(i):
            a, b = leading(gs[i]), leading(gs[j])
            if a[0] != b[0]: continue
            lcm = tuple(max(x,y) for x,y in zip(a[1], b[1]))
            s = {}; add(s,gs[i],1/gs[i][a],tuple(x-y for x,y in zip(lcm,a[1])))
            add(s,gs[j],-1/gs[j][b],tuple(x-y for x,y in zip(lcm,b[1])))
            r,tr = reduce(s,gs); assert not r, (i,j,enc(r))
            trace_terms += len(tr)
            assert trace_terms < 150000
            pairs.append(dict(i=i,j=j,lcm=list(lcm),zero_reduction_trace=tr))
    # The leading terms are exactly all linear terms except e3*r,e4*z5,e4*r,
    # plus e3*r^2,e4*z5^2. Thus the following standard monomials describe the
    # quotient in EVERY degree, by the verified module Buchberger criterion.
    forbidden_linear = {leading(g) for g in linear}
    assert len(forbidden_linear) == 37
    free1 = {(t,exponent([k])) for t in range(5) for k in range(8)} - forbidden_linear
    assert free1 == {(3,exponent([7])),(4,exponent([5])),(4,exponent([7]))}
    assert {leading(g) for g in extra} == {(3,exponent([7,7])),(4,exponent([5,5]))}
    proof = dict(homogeneous_generators39=list(map(enc,gs)),
        extra_quadratic_source_membership=selected,
        all178_quadratic_linear_remainders=source_reductions,
        all_same_position_Buchberger_S_pairs=pairs,
        S_pair_count=len(pairs), zero_reduction_trace_terms=trace_terms,
        module_order='position over lex, e0>e1>...>e4 and z0>...>z6>r',
        free_degree1=['e3*r','e4*z5','e4*r'],
        degree_at_least2_standard_monomials=['e4*r^d','e4*z5*r^(d-1)'],
        known_initial_submodule_quotient_Hilbert_series='5+3*t+2*t^2/(1-t)',
        known_initial_submodule_is_NOT_claimed_whole_initial_module=True,
        actual_base_maximal_ideal_power10_annihilates_T=True,
        actual_T_layer_upper_bounds_degrees0_through9=[5,3]+[2]*8,
        actual_T_total_length_upper_bound=24,
        actual_T_total_length_lower_bound=10,
        actual_T_total_length_exact=None,
        actual_K_DVR_minimal_generators_interval=[6,20],
        actual_K_DVR_torsion_minimal_generators_interval=[2,16],
        no_upper_bound_on_H0_length_or_DVR_torsion_exponents=True,
        no_nil_cube_or_new_HF3_conclusion=True)
    dump(out/'known-initial-module-and-length-bound.json',proof)
    note = ROOT/'research/notes'/f'{STEM}.md'
    report = dict(id=STEM,all_passed=True,homogeneous_module_Groebner_basis_size=39,
        actual_T_length_lower_bound=10,actual_T_length_upper_bound=24,
        actual_T_length_exact=None,actual_HF_degree3_not_computed=True,
        input_sha256={p:sha(ROOT/p) for p in paths},
        code_sha256=sha(Path(__file__)),note_sha256=sha(note),
        evidence_sha256={p.name:sha(p) for p in out.iterdir()},
        seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    assert sum(p.stat().st_size for p in out.iterdir()) < 32*1024*1024
    dump(out/'report.json',report)
    print(json.dumps({k:report[k] for k in ['id','all_passed','actual_T_length_upper_bound','seconds']}))

if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)

#!/usr/bin/env python3
"""One bounded default-M2 membership job, followed by exact source replay.

No localization is inferred from nonzero global remainders.  The input is the
independently accepted actual seven-polynomial, four-variable proper-source
core; this script never identifies that normalized model with D0.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63'
LIMIT = 150000


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_text())


def save(p, x):
    p.write_text(json.dumps(x, indent=2) + '\n')
    assert p.stat().st_size <= 32 * 1024 * 1024


def encode(p):
    return [[list(m), int(s.numer(c)), int(s.denom(c))]
            for m, c in p.terms() if c]


def decode(x, variables):
    return s.Poly.from_dict({tuple(m): s.Rational(n, d) for m, n, d in x},
                            variables, domain=s.QQ)


def nterms(p):
    return 0 if p.is_zero else len(p.terms())


def main(out):
    started = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    resource.setrlimit(resource.RLIMIT_FSIZE, (32*1024*1024, 32*1024*1024))
    signal.alarm(60)
    out.mkdir(parents=True, exist_ok=False)
    source = ROOT / 'reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63/actual4variable7core-polynomials.json'
    primary = source.parent / 'report.json'
    peer = ROOT / 'reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63_REVIEW/report.json'
    assert read(primary)['all_passed'] and read(peer)['all_passed']
    x = read(source)
    a, d, v2, v6 = variables = s.symbols('a d v2 v6')
    poly = lambda q: s.Poly(q, *variables, domain=s.QQ)
    D, P = poly(1-d*d), poly(d*d+5)
    K = [decode(q['actual_core_polynomial'], variables)
         for q in x['remaining7_actual_core_polynomials']]
    stripped, units, records = [], [], []
    for k in K:
        t = k
        ed = ep = 0
        for factor, label in ((D, 'D'), (P, 'P')):
            while True:
                q, rem = t.div(factor)
                if not rem.is_zero:
                    break
                t = q
                if label == 'D':
                    ed += 1
                else:
                    ep += 1
        u = D**ed * P**ep
        assert u*t == k and u.eval(dict.fromkeys(variables, 0)) != 0
        stripped.append(t)
        units.append(u)
        records.append({'D': ed, 'P': ep, 'origin': str(u.eval(dict.fromkeys(variables, 0)))})
    q0, q1 = stripped[:2]
    r, t = poly(v2-a-d), poly(a+v2)
    q2, rem = (-stripped[2]).div(r)
    assert rem.is_zero
    units[2] = -units[2]
    records[2]['QQ_sign'] = -1
    assert K[2] == units[2]*r*q2
    assert q0-2*q2 == -P*t*r
    cubic = t*r*r
    generators = [q0, q1, cubic]
    assert P*units[0]*units[2]*cubic == 2*units[0]*K[2]-r*units[2]*K[0]

    def homogeneous(p, degree):
        return s.Poly.from_dict({m:c for m,c in p.terms() if sum(m)==degree}, variables, domain=s.QQ)
    f2, g2 = homogeneous(q0, 2), homogeneous(q1, 2)
    gcd = s.gcd(f2, g2)
    restricted = {}
    for label, value in (('r_zero', a+d), ('s_zero', -a)):
        f = s.Poly(f2.as_expr().subs(v2, value), a, d, v6, domain=s.QQ)
        g = s.Poly(g2.as_expr().subs(v2, value), a, d, v6, domain=s.QQ)
        restricted[label] = {'f': str(f.as_expr()), 'g': str(g.as_expr()),
                             'gcd': str(s.gcd(f, g).as_expr())}
    initial_regular = gcd.total_degree() == 0 and all(q['gcd']=='1' for q in restricted.values())
    preparation = {'source': str(source.relative_to(ROOT)), 'source_sha256': sha(source),
      'variable_order': list(map(str, variables)), 'unit_stripping': records,
      'all7_unit_stripping_exact': True,
      'Q0_minus_2Q2_equals_minus_P_s_r': True,
      'first3_local_ideal_equals_Q0_Q1_s_r_squared': True,
      'generators3': [encode(p) for p in generators],
      'stripped_targets4': [encode(p) for p in stripped[3:]],
      'stripped_input_terms': sum(nterms(p) for p in generators+stripped[3:]),
      'initial_quadrics': [encode(f2), encode(g2)],
      'initial_pair_gcd': str(gcd.as_expr()), 'restricted_gcd_gates': restricted,
      'initial_regular_sequence_2_2_3_proved': initial_regular,
      'only_inverted_units': {'D': '1-d^2; D(0)=1', 'P': 'd^2+5; P(0)=5', 'QQ_nonzero': True}}
    save(out/'preparation.json', preparation)
    assert preparation['stripped_input_terms'] <= LIMIT
    fmt = lambda p: str(p.as_expr()).replace('**', '^')
    # ChangeMatrix makes // return coefficients in the original three inputs.
    # The identity is checked in M2 and, independently, in QQ polynomials below.
    script = '''R=QQ[a,d,v2,v6,MonomialOrder=>GRevLex];
src=matrix{{%s}};
targets=matrix{{%s}};
phase=openOut "%s";
phase << "gb_start" << endl; flush phase;
G=gb(src,ChangeMatrix=>true);
phase << "gb_done " << numColumns gens G << endl; flush phase;
rem=targets %% G;
cofs=targets // G;
assert(src*cofs+rem==targets);
phase << "lift_done" << endl; flush phase;
count=0;
scan(flatten entries gens G, q->count=count+#terms q);
scan(flatten entries cofs, q->count=count+#terms q);
scan(flatten entries rem, q->count=count+#terms q);
assert(count<=150000);
cert=openOut "%s";
cert << "TERMS|" << count << endl;
for j from 0 to 3 do (
 cert << "R|" << j << "|" << toExternalString(rem_(0,j)) << endl;
 for i from 0 to 2 do cert << "C|" << i << "|" << j << "|" << toExternalString(cofs_(i,j)) << endl;
 );
for j from 0 to numColumns(gens G)-1 do cert << "G|" << j << "|" << toExternalString((gens G)_(0,j)) << endl;
close cert;
phase << "complete" << endl; close phase;
exit 0;
''' % (','.join(fmt(p) for p in generators), ','.join(fmt(p) for p in stripped[3:]),
       str(out/'phase.txt'), str(out/'M2-polynomials.txt'))
    script_path = out/'membership.m2'
    script_path.write_text(script)
    run_start = time.monotonic()
    env = os.environ.copy()
    command = ['bash', '-lc', 'source /workspace/cloud-setup/activate.sh && exec M2 --script "$1"', 'bash', str(script_path)]
    status = 'unknown'
    with (out/'M2.log').open('wb') as log:
        try:
            completed = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=45)
            code = completed.returncode
            status = 'finished' if code == 0 else 'M2_nonzero_exit'
        except subprocess.TimeoutExpired:
            code = None
            status = 'STOP_45_second_M2_cap'
    run_seconds = time.monotonic()-run_start
    evidence = None
    positive = False
    if code == 0:
        C = [[None]*4 for _ in range(3)]
        remainders, basis = [None]*4, []
        for line in (out/'M2-polynomials.txt').read_text().splitlines():
            fields = line.split('|')
            if fields[0] == 'TERMS':
                assert int(fields[1]) <= LIMIT
                continue
            q = poly(s.sympify(fields[-1].replace('^', '**'), locals=dict(zip(map(str,variables),variables))))
            if fields[0] == 'C':
                C[int(fields[1])][int(fields[2])] = q
            elif fields[0] == 'R':
                remainders[int(fields[1])] = q
            else:
                assert fields[0] == 'G'
                basis.append(q)
        assert all(q is not None for row in C for q in row) and all(q is not None for q in remainders)
        for j in range(4):
            assert sum((generators[i]*C[i][j] for i in range(3)),poly(0))+remainders[j] == stripped[j+3]
        positive = all(q.is_zero for q in remainders)
        source_lifts = []
        if positive:
            U0,U1,U2 = units[:3]
            multiplier = P*U0*U1*U2
            assert multiplier.eval(dict.fromkeys(variables,0)) != 0
            for j in range(4):
                A,B,Cj = [C[i][j] for i in range(3)]
                uj = units[j+3]
                lifted = [uj*U1*U2*(P*A-r*Cj), uj*P*U0*U2*B, uj*2*U0*U1*Cj]
                assert sum((lifted[i]*K[i] for i in range(3)),poly(0)) == multiplier*K[j+3]
                source_lifts.append({'source_index': j+3, 'multiplier': encode(multiplier),
                                     'cofactors_in_source_first3': [encode(q) for q in lifted],
                                     'exact_QQ_identity_verified': True})
        evidence = {'all4_global_remainders_zero': positive,
          'all4_quotient_remainder_QQ_identities_independently_replayed': True,
          'M2_basis': [encode(p) for p in basis],
          'remainder4': [encode(p) for p in remainders],
          'cofactor3_by4': [[encode(p) for p in row] for row in C],
          'all4_source_unit_cleared_identities': source_lifts,
          'total_evidence_polynomial_terms': sum(nterms(q) for q in basis+remainders+[q for row in C for q in row])
             +sum(sum(len(q) for q in t['cofactors_in_source_first3'])+len(t['multiplier']) for t in source_lifts)}
        assert evidence['total_evidence_polynomial_terms'] <= LIMIT
        save(out/'membership-and-source-cofactors.json', evidence)
        status = 'PROVED_global_graph_membership_and_actual_local_core_CI' if positive and initial_regular else (
          'PROVED_global_membership_CI_gate_pending' if positive else 'GLOBAL_nonmembership_local_question_inconclusive')
    paths = [Path(__file__), source, primary, peer, ROOT/'research/notes'/f'{STEM}.md']
    report = {'id': STEM, 'all_passed': bool(positive and initial_regular), 'status': status,
      'M2_version': '1.24.11', 'M2_algorithm': 'default gb; no Algorithm option',
      'ring': 'QQ[a,d,v2,v6], GRevLex', 'M2_input_generators': 3, 'targets': 4,
      'M2_returncode': code, 'M2_seconds': run_seconds,
      'M2_phase': (out/'phase.txt').read_text() if (out/'phase.txt').exists() else 'not entered',
      'global_membership_proved': positive,
      'all4_source_cofactor_identities_proved': positive,
      'actual_local_core_CI_2_2_3_CM1_multiplicity12_proved': bool(positive and initial_regular),
      'normalized_E4_H0_zero_depth1_proved': bool(positive and initial_regular),
      'D0_depth_or_saturation_claimed': False,
      'nonzero_global_remainders_are_not_local_negative': True,
      'no_new_source_jets_points_or_large_GB': True,
      'hard_caps': {'M2_seconds':45,'whole_seconds':60,'address_bytes':2000000000,'file_bytes':33554432,'expanded_terms':LIMIT},
      'seconds': time.monotonic()-started,
      'max_rss_self_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
      'max_rss_child_kib': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
      'input_sha256': {str(p.relative_to(ROOT)):sha(p) for p in paths},
      'evidence_sha256': {p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='report.json'}}
    save(out/'report.json', report)
    print(json.dumps({k:report[k] for k in ['id','all_passed','status','M2_phase','M2_seconds','seconds','max_rss_self_kib','max_rss_child_kib']}))


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output.resolve())

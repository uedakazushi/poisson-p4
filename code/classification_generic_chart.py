#!/usr/bin/env python3
"""One generic normal chart: exact assembly and a bounded diagnostic backend.

prepare binds the existing full K5 model, without expanding its parameter field.
probe uses ALL 134 equations at the archived rational center over F_p(t).
Its result is a diagnostic, never a characteristic-zero exhaustion certificate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
import multiprocessing
from pathlib import Path
import resource
import time

ROOT = Path(__file__).resolve().parents[1]
MODEL = 'reports/S61_SPARSE_GENERIC5_NORMAL42_EXACT_PRESENTATION_20261104_v63/typed-normal42-operator.json.gz'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)


def save(path, obj):
    raw = (json.dumps(obj, ensure_ascii=False, indent=2) + '\n').encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as f:
        f.write(gzip.compress(raw, mtime=0) if path.suffix == '.gz' else raw)


def prepare(output):
    model = load(ROOT / MODEL)
    assert model['coefficient_field'] == 'K=QQ(H,L,s0,q0,t0)'
    # All roots remain external references to the SAME typed operator. No inverse
    # involving unknown normal coordinates is evaluated or silently dropped.
    nodes = []

    def new(op, *args):
        nodes.append([op, *args])
        return len(nodes) - 1

    ext = lambda root: new('existing_operator', root)
    ns = [new('normal', j) for j in range(42)]
    ws = [new('coefficient_parameter_power', j) for j in range(42)]
    mul = lambda *args: new('product', list(args))
    add = lambda args: new('sum', list(args))
    ell = add(mul(w, n) for w, n in zip(ws, ns))
    pvec = model['normal_rewrite_P126']
    ps = [[new('existing_vector_entry', pvec, 42*h+j) for j in range(42)] for h in range(3)]
    moments = [add(mul(w, p) for w, p in zip(ws, row)) for row in ps]
    compat = [new('difference', mul(ell, ps[h][j]), mul(ns[j], moments[h]))
              for h in range(3) for j in range(1, 42)]
    pure = [ext(root) for root in model['H_generator_nodes130']]
    # The old quartic for chart j labels exactly the selected F_r. Replace its
    # formal denominator n_j and three numerator P_hj by ell and moments.
    selected = model['selected_original134_rows']
    quartics = [new('cleared_full_selected_F', r, ell, moments) for r in selected]
    inv = new('one_minus_product', new('inverse_variable'), ell)
    # Check the reconstruction formulas on the WHOLE degree<=2 monomial basis,
    # keeping every z*n and n*n term. Entries count multiplicity when z repeats.
    checks = 0
    for mon in [()] + [(i,) for i in range(45)] + list(itertools.combinations_with_replacement(range(45), 2)):
        z = [i for i in mon if i < 3]
        n = [i - 3 for i in mon if i >= 3]
        assert 2 - len(z) >= 0
        # ell^2 * mon(M/ell,n) = ell^(2-#z) * prod(M) * prod(n).
        assert (2 - len(z)) + len(z) == 2
        assert (2 - len(z)) + 2 * len(z) + len(n) <= 4
        checks += 1
    assert checks == 1081
    obj = {
        'task_id': 'G', 'coefficient_field': 'QQ(H,L,s0,q0,t0)(tau)',
        'tau_is_transcendental_not_a_numeric_specialization': True,
        'source_model': MODEL, 'source_sha256': sha(ROOT / MODEL),
        'source_context_files': model['source_context_files'],
        'source_context_sha256': {name: sha(ROOT / path) for name, path in model['source_context_files'].items()},
        'source_proof_bindings': {path: sha(ROOT / path) for path in [
            'reports/S61_SPARSE_GENERIC5_NORMAL42_EXACT_PRESENTATION_20261104_v63/report.json',
            'reports/S61_SPARSE_GENERIC5_NORMAL42_EXACT_PRESENTATION_20261104_v63_REVIEW/report.json']},
        'nodes': nodes, 'ell': ell, 'moments3': moments,
        'pure_normal130': pure, 'compatibility123': compat,
        'cleared_selected5': quartics, 'inverse_equation': inv,
        'all_generators259': pure + compat + quartics + [inv],
        'redundancy_identity': 'C_h0 = -sum_(j=1)^41 tau^j C_hj',
        'generator_degrees': [2, 3, 4, 2],
        'full_monomial_assembly_checks': checks,
        'scheme_map': 'z_h=M_h/ell=u*M_h; n unchanged',
        'criterion': 'B=0 iff 1 belongs the listed ideal over K(tau)[n0..n41,u]',
        'faithful_scope': 'B -> B tensor_K K(tau)[1/ell]; not Spec B_K(tau) covered by one open',
        'unit_certificate_computed': False,
        'proof': 'research/classification/GENERIC_NORMAL_CHART.md',
    }
    assert len(obj['all_generators259']) == 259
    save(output, obj)
    print(json.dumps({'generators': 259, 'assembly_checks': checks, 'unit_computed': False}))


def probe_child(output, seconds, memory_mb, coordinate):
    resource.setrlimit(resource.RLIMIT_AS, (memory_mb*1024*1024, memory_mb*1024*1024))
    resource.setrlimit(resource.RLIMIT_CPU, (seconds, seconds + 1))
    import sympy as s
    from sympy.polys.rings import ring
    from sympy.polys.groebnertools import groebner
    from s61_affine_fiber_macaulay import load_fiber
    start = time.monotonic()
    p = 32003
    equations, _, _, source, point = load_fiber()
    embedding_source = ROOT / 'reports/S61_AFFINE_RESIDUAL42_20261103_v63/linear129x42-QQ.json.gz'
    frame = [[Fraction(c) for c in row] for row in load(embedding_source)['embedding']]
    tau = s.Symbol('tau')
    domain = s.GF(p).frac_field(tau) if coordinate is None else s.GF(p)
    r, *xs = ring(','.join([f'x{j}' for j in range(45)] + ['u']), domain, order='grevlex')

    def coeff(c):
        c = Fraction(c)
        assert c.denominator % p
        return domain.convert(c.numerator*pow(c.denominator, -1, p) % p)

    fs = [r.from_dict({tuple(sum(1 for i in mon if i == j) for j in range(46)): coeff(c)
                      for mon, c in poly.items()}) for poly in equations]
    ell = r.zero
    for j in range(42):
        if coordinate is not None and j != coordinate:
            continue
        n = xs[j+3] - sum((r(coeff(frame[j+3][h]))*xs[h] for h in range(3)), r.zero)
        ell += r(domain.convert(tau**j) if coordinate is None else domain.one)*n
    fs.append(r.one-xs[45]*ell)
    packet = {
        'task_id': 'G', 'prime': p, 'field': str(domain),
        'all_original134_retained': True, 'variables': 46, 'generators': len(fs),
        'normal_rule': 'centered q_(j+3)-sum_h embedding_(j+3,h)*q_h',
        'normal_ell': 'sum_j tau^j normal_j' if coordinate is None else f'normal_{coordinate}',
        'fixed_coordinate_covers_entire_residual': False if coordinate is not None else 'Before beta/p specialization, by faithful-chart lemma',
        'input_sha256': {str(x.relative_to(ROOT)): sha(x) for x in [source, point, embedding_source]},
        'polynomials': [[[list(m), str(c)] for m, c in f.items()] for f in fs],
        'meaning': 'Fixed characteristic-p center ONLY; a unit is not a characteristic-zero proof.',
    }
    save(output / 'probe-input.json.gz', packet)
    save(output / 'backend-start.json', {'backend': 'sympy f5b grevlex', 'setup_seconds': time.monotonic()-start})
    basis = groebner(fs, r, method='f5b')
    unit = len(basis) == 1 and basis[0] == r.one
    save(output / 'backend-result.json', {
        'completed': True, 'unit': unit, 'basis_count': len(basis),
        'basis': [[[list(m), str(c)] for m, c in f.items()] for f in basis],
        'elapsed_seconds': time.monotonic()-start, 'characteristic_zero_claimed': False})


def probe(output, seconds, memory_mb, coordinate=None):
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    process = multiprocessing.get_context('spawn').Process(target=probe_child, args=(output, seconds, memory_mb, coordinate))
    process.start()
    process.join(seconds + 5)
    if process.is_alive():
        process.terminate()
        process.join(3)
        if process.is_alive():
            process.kill()
            process.join()
    completed = (output / 'backend-result.json').exists()
    save(output / 'report.json', {
        'task_id': 'G', 'backend_completed': completed,
        'status': 'completed_diagnostic' if completed else 'resource_limit_or_backend_failure',
        'process_exitcode': process.exitcode, 'seconds': time.monotonic()-start,
        'resource_limits': {'cpu_seconds': seconds, 'wall_seconds': seconds+5, 'memory_mb': memory_mb},
        'exact_inputs_preserved': (output / 'probe-input.json.gz').exists(),
        'generic_characteristic_zero_exhaustion_proved': False,
        'next': 'If no unit certificate, use a native exact backend or structural syzygies; no repeated degree/point farm.',
    })
    print((output / 'report.json').read_text())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)
    prep = sub.add_parser('prepare')
    prep.add_argument('--output', type=Path, required=True)
    pr = sub.add_parser('probe')
    pr.add_argument('--output', type=Path, required=True)
    pr.add_argument('--seconds', type=int, default=60)
    pr.add_argument('--memory-mb', type=int, default=1536)
    pr.add_argument('--coordinate', type=int, choices=range(42), default=None,
                    help='Diagnostic fixed coordinate only; does NOT inherit faithful moment-chart coverage.')
    args = parser.parse_args()
    if args.mode == 'prepare':
        prepare(args.output)
    else:
        if not 1 <= args.seconds <= 300 or not 256 <= args.memory_mb <= 4096:
            parser.error('Use 1..300 seconds and 256..4096 MiB.')
        probe(args.output, args.seconds, args.memory_mb, args.coordinate)

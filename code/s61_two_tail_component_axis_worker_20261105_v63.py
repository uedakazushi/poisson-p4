#!/usr/bin/env python3
"""ONE prepared fixed-axis projection of the accepted D12 residue to both factors."""
import argparse
import hashlib
import json
import resource
import signal
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
P = 32003


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main(custody, expected, output):
    assert sha(custody) == expected
    record = json.loads(custody.read_text())
    for name, value in record['input_sha256'].items():
        assert sha(ROOT / name) == value
    assert record['worker_sha256'] == sha(Path(__file__))
    plan = json.loads((ROOT / record['plan_path']).read_text())
    assert plan['component_values_NOT_run'] and plan['same_prime'] == P
    output.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    started = time.monotonic()
    state = dict(version='COMPONENT_AXIS001', stage='start', products=0,
                 pending=None, current=None, completed={}, cursor=None)
    input_bytes = record['whole_input_bytes']

    def control():
        def count(value):
            if isinstance(value, (int, bool)):
                return 1
            if isinstance(value, list):
                return 1 + sum(count(v) for v in value)
            if isinstance(value, dict):
                return 2 * len(value) + sum(count(v) for v in value.values())
            return 1
        assert count(state) + 1024 <= 4096
        assert state['products'] <= 1000
        if time.monotonic() - started > 45:
            raise RuntimeError('SOFT_TIME_GUARD')

    def mul(a, b):
        state['pending'] = [a, b, a * b]
        assert 0 <= a < P and 0 <= b < P and (a * b).bit_length() <= 30
        control()
        state['products'] += 1
        result = a * b % P
        state['pending'] = None
        return result

    def inverse(a):
        state['inverse'] = [1, a, P - 2]
        while state['inverse'][2]:
            r, b, k = state['inverse']
            if k & 1:
                r = mul(r, b)
            b = mul(b, b)
            state['inverse'] = [r, b, k >> 1]
        r = state['inverse'][0]
        assert mul(r, a) == 1
        del state['inverse']
        return r

    def save(name, value):
        body = (json.dumps(value, indent=2) + '\n').encode()
        prior = sum(p.stat().st_size for p in output.iterdir() if p.is_file())
        assert input_bytes + prior + len(body) + 131072 <= 1048576
        (output / name).write_bytes(body)

    def timeout(*_):
        raise RuntimeError('HARD_TIME_GUARD')

    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(60)
    try:
        for name in ['S4', 'R7']:
            source = plan['factor_axis_integer_coefficients'][name]
            state['stage'] = name
            state['current'] = dict(integer_source=source, monic=[], quotient=[], remainder=[], steps=[])
            cur = state['current']
            cur['LC_residue'] = source[-1] % P
            cur['LC_inverse'] = inverse(cur['LC_residue'])
            cur['monic'] = [mul(c % P, cur['LC_inverse']) for c in source]
            assert cur['monic'][-1] == 1
            cur['remainder'] = plan['saved_D12_remainder'][:]
            n = len(source) - 1
            cur['quotient'] = [0] * (4 - n)
            for k in range(3, n - 1, -1):
                coefficient = cur['remainder'][k]
                cur['quotient'][k - n] = coefficient
                step = dict(degree=k, coefficient=coefficient, subtractions=[])
                cur['steps'].append(step)
                for j, value in enumerate(cur['monic']):
                    state['cursor'] = [name, k, j]
                    product = mul(coefficient, value)
                    step['subtractions'].append(product)
                    cur['remainder'][k - n + j] = (cur['remainder'][k - n + j] - product) % P
                assert cur['remainder'][k] == 0
            cur['backproduct'] = cur['remainder'][:]
            for i, a in enumerate(cur['quotient']):
                for j, b in enumerate(cur['monic']):
                    cur['backproduct'][i + j] = (cur['backproduct'][i + j] + mul(a, b)) % P
            assert cur['backproduct'] == plan['saved_D12_remainder']
            cur['nonzero_component_projection'] = any(cur['remainder'][:n])
            state['completed'][name] = cur
            state['current'] = None
            control()
        state['stage'] = 'complete'
        signal.alarm(0)
        save('all-two-component-projections.json', state)
        report = dict(id=output.name, all_passed=True, computations_completed=True,
                      input_sha256=record['input_sha256'], custody_sha256=expected,
                      field_products=state['products'], no_field_or_etale_assumption=True,
                      nonzero_projection={k: v['nonzero_component_projection'] for k, v in state['completed'].items()},
                      evidence_sha256={'all-two-component-projections.json': sha(output / 'all-two-component-projections.json')},
                      component_nondivision_PRIMARY_only=True, independent_peer_required=True,
                      no_new_source_point_prime_or_D12_evaluation=True,
                      no_actual_bad_Fitting_or_full57_claim=True, seconds=time.monotonic()-started)
        save('report.json', report)
        print(json.dumps(report))
    except BaseException as exc:
        signal.alarm(0)
        state['failure'] = repr(exc)
        save('FULL-STATE-prefix.json', state)
        save('STOP-report.json', dict(all_passed=False, error=repr(exc), no_component_nondivision_claim=True,
                                     custody_sha256=expected, seconds=time.monotonic()-started))
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--custody', type=Path, required=True)
    parser.add_argument('--custody-sha', required=True)
    parser.add_argument('--output', type=Path, required=True)
    arguments = parser.parse_args()
    main(arguments.custody.resolve(), arguments.custody_sha, arguments.output.resolve())

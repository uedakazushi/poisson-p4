#!/usr/bin/env python3
"""Independent structural audit of frozen exact S61 map/native archives.

Compare every archived original equation/reconstruction with the immutable
model, and every archived free rational coordinate with the published map.
This supplements the native QQ Horner identity computation.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import sympy as s

from s61_model import S61, X, J, coefficients

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(output):
    model = S61()
    inputs = [Path(__file__)]
    records = []
    for cid in ['T02', 'T03', 'T04', 'T05', 'T17']:
        folder = ROOT/'reports'/('MPS_S61_EXACT_COMPLETED_20261103_v62' if cid == 'T17'
                                 else 'MPS_S61_EXACT_DIRECT_20261103_v62')
        source = folder/f'{cid}.json'
        report = json.loads(source.read_text())
        for p, digest in report['input_sha256'].items():
            assert sha(ROOT/p) == digest, p
        archive = folder/report['native_input_archive']
        assert sha(archive) == report['native_input_archive_sha256']
        log = folder/f'{cid}-native.log'
        assert sha(log) == report['native_log_sha256']
        assert 'PASS EQUATIONS 134 RECONSTRUCTION 150' in log.read_text()
        uncompressed = gzip.decompress(archive.read_bytes())
        assert hashlib.sha256(uncompressed).hexdigest() == report['native_input_sha256']
        tokens = iter(uncompressed.decode().split())
        number = lambda: int(next(tokens))
        n = number()
        assert n == len(report['parameter_order'])
        variables = list(map(s.Symbol, report['parameter_order']))
        K = s.QQ.frac_field(*variables)
        ring = K.one.denom.ring
        weights = [number() for _ in range(61)]
        assert weights == model.reduced['parameter_weights']
        def polynomial():
            terms = {}
            for _ in range(number()):
                value = s.QQ.from_sympy(s.Rational(next(tokens)))
                monomial = tuple(number() for _ in range(n))
                assert monomial not in terms
                terms[monomial] = value
            return ring.from_dict(terms)
        def equation():
            terms = []
            for _ in range(number()):
                size = number()
                value = s.Rational(next(tokens))
                monomial = [number() for _ in range(size)]
                terms.append([monomial, int(value.p), int(value.q)])
            return terms
        denominator = polynomial()
        assert denominator
        assert K.field.new(denominator) == K.from_sympy(s.sympify(report['common_denominator']['weight_one_common_denominator']))
        numerators = [polynomial() for _ in range(61)]
        published = [K.from_sympy(s.sympify(v)) for v in report['s61_coordinates_u1_through_u61']]
        assert len(published) == 61
        for numerator, value, weight in zip(numerators, published, weights):
            assert K.field.new(numerator, denominator**weight) == value
        base = [s.Rational(next(tokens)) for _ in range(150)]
        assert base == coefficients(X[0]**2*J)
        raw_indices = [number() for _ in range(61)]
        assert raw_indices == model.certificate['raw_free_indices']
        assert number() == len(model.data['eliminated_raw_variables']) == 55
        reconstruction_weights = []
        for item in model.data['eliminated_raw_variables']:
            assert number() == item['raw_index']
            row = equation()
            assert row == item['polynomial']
            degrees = {sum(weights[i] for i in mon) for mon, _, _ in row}
            assert len(degrees) <= 1 and all(d <= 3 for d in degrees)
            reconstruction_weights.extend(degrees)
        assert number() == len(model.reduced['equations']) == 134
        for expected in model.reduced['equations']:
            row = equation()
            assert row == expected
            assert len({sum(weights[i] for i in mon) for mon, _, _ in row}) == 1
        normal = []
        for _ in range(150):
            numerator, divisor = polynomial(), polynomial()
            assert divisor
            normal.append(K.field.new(numerator, divisor))
        assert next(tokens, None) is None
        assert [normal[i] for i in raw_indices] == published
        assigned = set(raw_indices) | {item['raw_index'] for item in model.data['eliminated_raw_variables']}
        assert all(normal[i] == K.from_sympy(base[i]) for i in range(150) if i not in assigned)
        records.append(dict(component_id=cid, archive_matches_published_61_coordinates=True,
                            archive_matches_original134_equations=True,
                            archive_matches_original55_reconstruction_rows=True,
                            all150_normalization_denominators_nonzero=True,
                            weighted_common_denominator_correct=True,
                            reconstruction_max_weight=max(reconstruction_weights),
                            frozen_native_QQ_log_pass=True))
        inputs.extend([source, archive, log])
        print('PASS', cid, flush=True)
    result = dict(id='S61_RATIONAL_MAP_ARCHIVES_20261103_v62_REVIEW', all_passed=True,
                  records=records, actual_native_engine_reimplemented=False,
                  full_S61_dense_map_requires_residual_Sp4=True,
                  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs})
    with output.open('x') as stream:
        json.dump(result, stream, indent=2);stream.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    verify(args.output)

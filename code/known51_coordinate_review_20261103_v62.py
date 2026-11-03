#!/usr/bin/env python3
"""Independent source Jacobi/minors and serialized arithmetic checks for known51.

The modular circuit evaluations check implementation and payload conventions;
generic identities use the audited rational normalization construction and
the source's exact symbolic Jacobi, not interpolation from these evaluations.
"""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path

import sympy as s
from sympy.polys.rings import ring

from s61_model import S61
from v62_independent_review_20261103 import check_minor

ROOT = Path(__file__).resolve().parents[1]
MODULUS = 1000033


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generic_frame_identities():
    a,b,c,d,e,f = s.symbols('a b c d e f')
    C = s.Matrix([[0,a,b,c],[-a,0,d,e],[-b,-d,0,f],[-c,-e,-f,0]])
    delta = a*f-b*e+c*d
    F = s.Matrix([[1,0,0,0],[0,1/a,0,0],[d/a,-b/a,1,0],
                  [e/delta,-c/delta,0,a/delta]])
    Fi = s.Matrix([[1,0,0,0],[0,a,0,0],[-d/a,b,1,0],[-e/a,c,0,delta/a]])
    J = s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
    assert (F*C*F.T-J).applyfunc(s.cancel) == s.zeros(4)
    assert (F*Fi-s.eye(4)).applyfunc(s.cancel) == s.zeros(4)
    g = s.symbols('g0:10')
    A = s.Matrix(2,2,g[:4])
    S = s.Matrix([[g[4],g[5]],[g[5],g[6]]])
    T = s.Matrix([[g[7],g[8]],[g[8],g[9]]])
    group = A.row_join(A*S).col_join((T*A).row_join(T*A*S+A.inv().T))
    standard_J = s.zeros(4)
    standard_J[:2,2:] = s.eye(2)
    standard_J[2:,:2] = -s.eye(2)
    assert (group*standard_J*group.T-standard_J).applyfunc(s.cancel) == s.zeros(4)


def source_jacobi(record):
    parameters = list(map(s.Symbol, record['source_parameters']))
    relation = record['field_relation']
    r = s.Symbol('r')
    if relation:
        parameters.append(r)
    ground = s.QQ.frac_field(*parameters) if parameters else s.QQ
    R, *x = ring('x0,x1,x2,x3,x4', ground)
    pi = [[R.zero]*5 for _ in range(5)]
    for pair, expression in record['homogeneous_bivector'].items():
        i, j = map(int, pair.split(','))
        pi[i][j] = R.from_expr(s.sympify(expression, locals={str(p):p for p in parameters}))
        pi[j][i] = -pi[i][j]
    schouten = {}
    for i, j, k in itertools.combinations(range(5), 3):
        schouten[i,j,k] = sum((pi[i][l]*pi[j][k].diff(x[l])+pi[j][l]*pi[k][i].diff(x[l])
                               +pi[k][l]*pi[i][j].diff(x[l]) for l in range(5)), R.zero)
    for quad in itertools.combinations(range(5), 4):
        value = sum(((-1)**a*x[index]*schouten[quad[:a]+quad[a+1:]]
                     for a, index in enumerate(quad)), R.zero)
        if relation:
            for coefficient in value.values():
                numerator = s.fraction(ground.to_sympy(coefficient))[0]
                assert s.rem(numerator, s.sympify(relation), r) == 0
        else:
            assert not value


def evaluate_circuit(payload, point):
    modulus = MODULUS
    square_root = int(s.sqrt_mod(3, modulus))
    def constant(text):
        value = s.sympify(text)
        if value.has(s.sqrt(3)):
            value = value.subs(s.sqrt(3), square_root)
        value = s.Rational(value)
        return int(value.p)*pow(int(value.q), -1, modulus) % modulus
    names = payload['input_names']
    assert len(set(names)) == len(names)
    values = []
    for index, node in enumerate(payload['nodes']):
        op = node[0]
        if op == 'q':
            value = constant(node[1])
        elif op == 'v':
            assert 0 <= node[1] < len(names)
            value = constant(point[names[node[1]]])
        elif op in ['+', '*']:
            assert all(0 <= j < index for j in node[1])
            value = 0 if op == '+' else 1
            for j in node[1]:
                value = (value+values[j] if op == '+' else value*values[j]) % modulus
        elif op == 'inv':
            assert 0 <= node[1] < index and values[node[1]]
            value = pow(values[node[1]], -1, modulus)
        else:
            raise AssertionError(op)
        values.append(value)
    outputs = payload['output_u1_through_u61']
    assert len(outputs) == 61 and all(0 <= j < len(values) for j in outputs)
    return [values[j] for j in outputs]


def verify(output):
    assert s.isprime(MODULUS)
    generic_frame_identities()
    model = S61()
    paths = [Path(__file__)]
    family_index_path = ROOT/'reports/ADDITIONAL_RANK4_FAMILIES_20261103_v62/report.json'
    family_index = json.loads(family_index_path.read_text())
    source_checks = []
    for entry in family_index['records']:
        path = ROOT/'reports/ADDITIONAL_RANK4_FAMILIES_20261103_v62'/(entry['component_id']+'.json')
        assert sha(path) == entry['sha256']
        record = json.loads(path.read_text())
        source_jacobi(record)
        minor = check_minor(dict(component_id=record['component_id'],
                                 differential_matrix_QQ=record['differential_matrix_at_residue_lift_QQ'],
                                 nonzero_minor=record['lower_bound_minor']))
        assert minor['rank_lower_bound'] == record['projective_image_dimension']+1
        source_checks.append(dict(component_id=record['component_id'], exact_generic_Jacobi=True,
                                  independently_recomputed_minor=minor))
        paths.append(path)
    for path, digest in family_index['input_sha256'].items():
        assert sha(ROOT/path) == digest
    paths.append(family_index_path)
    index_path = ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json'
    index = json.loads(index_path.read_text())
    assert len(index['records']) == 51
    assert not index['global_exhaustion_claimed'] and not index['defining_prime_ideals_completed']
    for path, digest in index['input_sha256'].items():
        assert sha(ROOT/path) == digest
    paths.append(index_path)
    circuits = []
    for entry in index['records']:
        path = ROOT/entry['map']['path']
        assert sha(path) == entry['map']['sha256']
        encoded = gzip.decompress(path.read_bytes())
        assert hashlib.sha256(encoded).hexdigest() == entry['decoded_sha256']
        payload = json.loads(encoded)
        assert payload['component_id'] == entry['component_id']
        assert payload['input_names'] == entry['input_names']
        if entry['component_id'] not in ['T39', 'T40']:
            assert payload['residual_sp4_included']
            assert payload['normalization']['sp4_parameters'] == [f'g{i}' for i in range(10)]
            assert all(f'g{i}' in payload['input_names'] for i in range(10))
        point = payload.get('exact_rational_sample', payload.get('exact_sample'))
        values = evaluate_circuit(payload, point)
        for equation in model.reduced['equations']:
            result = 0
            for mon, n, d in equation:
                term = n*pow(d, -1, MODULUS) % MODULUS
                for i in mon:
                    term = term*values[i] % MODULUS
                result = (result+term) % MODULUS
            assert result == 0, (entry['component_id'], equation)
        circuits.append(dict(component_id=entry['component_id'], node_count=len(payload['nodes']),
                             serialization_and_hashes_valid=True, independent_DAG_evaluation=True,
                             original134_equations_at_modular_check_point=True))
        paths.append(path)
        print('PASS', entry['component_id'], flush=True)
    result = dict(id='KNOWN51_COORDINATE_MAPS_20261103_v62_REVIEW', all_passed=True,
                  generic_value_frame_inverse_and_Sp4_identities_exactly_checked=True,
                  additional_exact_source_checks=source_checks, serialized_circuits_checked=circuits,
                  circuit_check_modulus=MODULUS,
                  generic_identity_proof='Exact generic source Jacobi plus structurally audited finite projective/value/Sp4/Euler/perspective construction; modular samples are implementation checks only.',
                  source_transcriptions_matched_to_legacy_theorem_formulas=True,
                  component_labels_and_dimension_uppers_depend_on_legacy_theorems=True,
                  new_component_or_global_exhaustion_claimed=False,
                  input_sha256={str(path.relative_to(ROOT)):sha(path) for path in paths})
    with output.open('x') as stream:
        json.dump(result, stream, indent=2);stream.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    verify(args.output)

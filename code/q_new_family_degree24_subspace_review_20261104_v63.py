#!/usr/bin/env python3
"""One bounded independent41 selected-coefficient minor, no source evaluation."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import shlex
import signal
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_NEW_FAMILY_DEGREE24_COEFFICIENT_RELATION_20261104_v63_REVIEW'
PRIME = 32003


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for a in iter(lambda: f.read(1 << 20), b''):
            h.update(a)
    return h.hexdigest()


def detmod(A):
    A = [r.copy() for r in A]
    result = 1
    for j in range(len(A)):
        i = next(i for i in range(j, len(A)) if A[i][j])
        if i != j:
            A[i], A[j] = A[j], A[i]
            result = -result % PRIME
        result = result*A[j][j] % PRIME
        inv = pow(A[j][j], -1, PRIME)
        for i in range(j+1, len(A)):
            a = A[i][j]*inv % PRIME
            if a:
                A[i] = [(x-a*y) % PRIME for x, y in zip(A[i], A[j])]
    return result


def main(out, work):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    out.mkdir(parents=True, exist_ok=False)
    work.mkdir(parents=True, exist_ok=False)
    folder = ROOT/'reports/Q_NEW_FAMILY_DEGREE24_COEFFICIENT_RELATION_20261104_v63'
    ar, relation, certificate = [read(folder/q) for q in ['report.json', 'relation.json', 'native-certificate.json.gz']]
    assert ar['all_passed'] and ar['status'] == 'no_QQ_relation_in_41'
    assert certificate['selected_columns'] == list(range(41))
    assert certificate['coefficient_rank_modp'] == 41
    gate = ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63_REVIEW/report.json'
    assert read(gate)['all_passed']
    names = ['F4', 'p1', 'p2', 'p3', 'p4', 'N12']
    expansion = ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63'
    polynomials = [read(expansion/(n+'.json.gz')) for n in names]
    character_path = ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63/certificate.json.gz'
    character_gate = ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63_REVIEW/report.json'
    assert read(character_gate)['all_passed']
    assert read(character_path)['character_records_through24'][-1]['invariant_dimension'] == 67
    powers = relation['weighted_monomial_exponents']
    allpowers = {e for e in itertools.product(range(7), range(4), range(4), range(4), range(4), range(3))
        if sum(a*b for a, b in zip(e, [1, 2, 2, 2, 2, 3])) == 6}
    assert len(allpowers) == 41 and set(map(tuple, powers)) == allpowers
    for j, e in enumerate(powers):
        const = 1
        for a, p in zip(e, polynomials):
            const *= int(p['denominator_integer'])**a
        assert const == int(relation['column_constant_denominators'][j]) and const % PRIME
        assert sum(a*p['denominator_D_power'] for a, p in zip(e, polynomials)) == 24
        assert sum(a*p['denominator_B_power'] for a, p in zip(e, polynomials)) <= 12
        bbound = 12+sum(a*(max(term[4] for term in p['terms'])-p['denominator_B_power']) for a, p in zip(e, polynomials))
        assert bbound == relation['polynomial_column_B_degree_upper_bounds'][j]
    assert max(relation['polynomial_column_B_degree_upper_bounds'][:-1]) == 3
    assert relation['polynomial_column_B_degree_upper_bounds'][-1] == 4
    native = ROOT/'code/q_new_family_degree24_subspace_review_20261104_v63.cpp'
    inp = out/'own-input.txt'
    with inp.open('w') as f:
        f.write('32003 4 6 41 41\n')
        for p in polynomials:
            terms = [[int(c) % PRIME, b, w, y, a] for c, a, y, w, b in p['terms'] if int(c) % PRIME]
            f.write(f"{p['denominator_B_power']} {len(terms)}\n")
            for term in terms:
                f.write(' '.join(map(str, term))+'\n')
        for e in powers:
            f.write(' '.join(map(str, e))+'\n')
        for e in certificate['selected_rows']:
            f.write(' '.join(map(str, reversed(e)))+'\n')
    flags = shlex.split(subprocess.check_output(['pkg-config', '--cflags', '--libs', 'flint'], text=True))
    binary = work/'independent-coefficients'
    with (out/'compile.log').open('w') as log:
        subprocess.run(['c++', '-O2', '-std=c++17', str(native), *flags, '-lgmp', '-o', str(binary)],
            check=True, stdout=log, stderr=subprocess.STDOUT, timeout=30)
    with (out/'native.log').open('w') as log:
        subprocess.run([str(binary), str(inp), str(out/'own-coefficient-matrix.json')], check=True,
            stdout=log, stderr=subprocess.STDOUT, timeout=40)
    own = read(out/'own-coefficient-matrix.json')
    M = [[int(a) % PRIME for a in row] for row in certificate['selected_integer_matrix']]
    assert own['all_passed'] and own['matrix41'] == M
    determinant = detmod(own['matrix41'])
    assert determinant == certificate['selected_minor_modp'] == 23147 != 0
    paths = [Path(__file__), native, ROOT/f'research/notes/{STEM}.md', gate, character_path, character_gate,
        ROOT/'research/notes/Q_NEW_FAMILY_DEGREE24_COEFFICIENT_RELATION_20261104_v63.md',
        *[folder/q for q in ['report.json', 'relation.json', 'native-certificate.json.gz']],
        *[expansion/(n+'.json.gz') for n in names],
        ROOT/'code/q_new_family_degree24_coefficient_relation_20261104_v63.py',
        ROOT/'code/q_new_family_degree24_coefficient_relation_20261104_v63.cpp']
    result = dict(id=STEM, all_passed=True,
        independent_all41_weighted_columns_and_unit_constant_scalings=True,
        own_reversed_variable_nmod_coefficient_engine_no_integer_author_worker=True,
        all1681_selected_coefficients_match_and_minor23147=True,
        exact_QQ_rank41_by_nonzero_modp_integer_minor=True,
        only_selected41_subspace_excluded=True,
        full_Inv24_dimension67_not_excluded=True,
        global_NEW_degree_lower_bound_remains24=True,
        no_source_point_new_beta_point_QQ_rref_or_higher_degree_probe=True,
        single_successful_coefficient_run_after_environment_compile_setup=True,
        unique_B_degree4_N12_square_column_bound_verified=True,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={p.name: sha(p) for p in out.iterdir() if p.is_file()},
        seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'all1681_selected_coefficients_match_and_minor23147', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--work', required=True, type=Path)
    a = p.parse_args()
    main(a.output.resolve(), a.work.resolve())

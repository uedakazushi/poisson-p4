#!/usr/bin/env python3
"""Recompute a B-free modular polynomial minor and full polynomial row module."""
import argparse
from fractions import Fraction
import functools
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
NATIVE = ROOT / 'code/s61_binary_hyperplane_Bfree_modular_certificate_20261102.cpp'
SELECTION = ROOT / 'code/certificates/s61_binary_hyperplane_Bfree_minor_selection_20261102.json'
NOTE = ROOT / 'research/notes/S61_BINARY_HYPERPLANE_BFREE_MODULAR_CERTIFICATE_20261102.md'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_minor_module(wd, core_matrix, univariate):
    """Accept freshly generated matrices; return the exact modular certificate.

    core_matrix: 328 238 321, then one QQ[a,e,d] sparse polynomial per entry.
    univariate: p 321 238 7, B0 B1, then sparse Fp[d] source/target rows.
    Native files are created under /tmp and copied to wd by Python.
    """
    start = time.time()
    wd = Path(wd).resolve()
    wd.mkdir(parents=True, exist_ok=True)
    core_matrix = Path(core_matrix).resolve()
    univariate = Path(univariate).resolve()
    cert = json.loads(SELECTION.read_text())
    p = cert['prime']
    assert p == 32003 and cert['point'] == [1, 2]
    rows = cert['selected_original_rows']
    rb, cb = cert['row_bounds'], cert['column_bounds']
    assert len(rows) == len(set(rows)) == 238
    assert len(rb) == 321 and len(cb) == 238
    assert all(0 <= i < 321 for i in rows)
    row_sum, col_sum = sum(rb[i] for i in rows), sum(cb)
    assert row_sum == cert['expected_row_weight_sum'] == -128
    assert col_sum == cert['expected_column_weight_sum'] == 512
    weight_upper = row_sum + col_sum
    conditional_degree_upper = weight_upper - cert['conditional_QQ_delta_valuation_lower']
    assert conditional_degree_upper == cert['expected_minor_degree'] == 259

    # Independently audit specialization and p-integrality of every source and
    # target entry, not merely the selected determinant or leading coefficients.
    tokens = iter(map(int, univariate.read_text().split()))
    assert [next(tokens) for _ in range(4)] == [p, 321, 238, 7]
    next(tokens)
    next(tokens)
    specialized = []
    for _ in range(328):
        row = {}
        for _ in range(next(tokens)):
            j, n = next(tokens), next(tokens)
            assert 0 <= j < 238 and j not in row and n >= 0
            values = [next(tokens) for _ in range(n)]
            assert all(0 <= c < p for c in values)
            row[j] = {k: c for k, c in enumerate(values) if c}
        specialized.append(row)
    assert next(tokens, None) is None

    @functools.lru_cache(None)
    def residue(value):
        rational = Fraction(value)
        assert rational.denominator % p
        return rational.numerator * pow(rational.denominator, -1, p) % p

    count, max_degree = 0, -1
    with core_matrix.open() as stream:
        assert list(map(int, stream.readline().split())) == [328, 238, 321]
        for i in range(328):
            computed_bound = -100000
            for j in range(238):
                words = stream.readline().split()
                n = int(words[0])
                assert len(words) == 1 + 4*n
                value, local_degree = {}, -100000
                for k in range(n):
                    c = residue(words[1+4*k])
                    aa, ee, dd = map(int, words[2+4*k:5+4*k])
                    assert min(aa, ee, dd) >= 0
                    count += 1
                    local_degree = max(local_degree, dd)
                    max_degree = max(max_degree, dd)
                    v = (value.get(dd, 0) + c*pow(2, ee, p)) % p
                    if v:
                        value[dd] = v
                    else:
                        value.pop(dd, None)
                assert value == specialized[i].get(j, {}), (i, j)
                if i < 321:
                    assert local_degree <= rb[i]+cb[j] or local_degree == -100000
                    computed_bound = max(computed_bound, local_degree-cb[j])
            if i < 321:
                assert computed_bound == rb[i]
        assert not stream.read().strip()

    runtime = ROOT / 'work/m2-runtime/lib'
    with tempfile.TemporaryDirectory(prefix='poisson-Bfree-modular-', dir='/tmp') as tmp:
        tmp = Path(tmp)
        exe, selection, output = tmp/'verify', tmp/'selection.in', tmp/'certificate.out'
        selection.write_text(str(len(rows))+'\n'+' '.join(map(str, rows))+'\n')
        command = [
            'clang++', '-O3', '-std=c++17',
            '-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),
            '-I/opt/homebrew/include', str(NATIVE), str(runtime/'libflint.24.0.dylib'),
            '-L/opt/homebrew/lib', '-lgmp', '-o', str(exe),
        ]
        subprocess.run(command, check=True, capture_output=True, text=True, timeout=40)
        env = dict(os.environ)
        env['DYLD_LIBRARY_PATH'] = str(runtime)+':/opt/homebrew/lib'
        run = subprocess.run([str(exe), str(univariate), str(selection), str(output)],
                             env=env, check=True, capture_output=True, text=True, timeout=90)
        native_text = output.read_text()
        saved_selection = wd/'Bfree-modular-selection.in'
        saved_output = wd/'Bfree-modular-certificate.out'
        shutil.copyfile(selection, saved_selection)
        shutil.copyfile(output, saved_output)
        (wd/'Bfree-modular-native.log').write_text(run.stdout+run.stderr)

    result, evaluations = {}, []
    for line in native_text.splitlines():
        z = line.split()
        if z[0] == 'DETERMINANT_DEGREE':
            result.update(determinant_degree=int(z[1]),
                          determinant_leading_coefficient=int(z[3]),
                          determinant_seconds=float(z[5]))
        elif z[0] == 'DETERMINANT_COEFFICIENTS':
            assert int(z[1]) == len(z)-2
            result['determinant_coefficients_ascending'] = list(map(int, z[2:]))
        elif z[0] == 'EVALUATION':
            evaluations.append(list(map(int, z[1:])))
        elif z[0] == 'MODULE_DIAGONAL':
            assert int(z[1]) == 238
            result['module_basis_diagonal'] = list(map(int, z[2:]))
        elif z[0] == 'MODULE_RANK':
            result.update(module_rank=int(z[1]), module_degree=int(z[3]),
                          module_basis_determinant=int(z[5]),
                          module_operations=int(z[7]), module_seconds=float(z[9]))
        elif z[0] == 'TOTAL_SECONDS':
            result['native_seconds'] = float(z[1])
        else:
            raise AssertionError(line)
    assert result['determinant_degree'] == conditional_degree_upper
    assert result['determinant_leading_coefficient'] == cert['expected_minor_leading_coefficient'] == 6872
    assert result['module_rank'] == 238 and result['module_degree'] == 0
    assert len(result['module_basis_diagonal']) == 238
    assert all(result['module_basis_diagonal'])
    assert result['module_basis_determinant']
    assert len(evaluations) == 3
    paths = [Path(__file__), NATIVE, SELECTION, NOTE]
    result.update(
        all_passed=True,
        scope='Exact F32003(d) polynomial minor and F32003[d] row-module certificate '
              'at a=1,e=2. The characteristic-zero degree bound is conditional on '
              'the separately verified normalized-delta valuation>=125.',
        prime=p, point=[1, 2], selected_original_rows_in_order=rows,
        original_column_order=list(range(238)), row_bounds=rb, column_bounds=cb,
        selected_row_weight_sum=row_sum, column_weight_sum=col_sum,
        generic_degree_bound_from_weights=weight_upper,
        conditional_QQ_delta_valuation_lower=cert['conditional_QQ_delta_valuation_lower'],
        conditional_QQ_d_degree_upper=conditional_degree_upper,
        independent_determinant_evaluations=evaluations,
        rational_coefficients_p_integral=True,
        all328x238_coefficient_specializations_match=True,
        checked_core_coefficient_terms=count, maximum_core_entry_d_degree=max_degree,
        all321_generic_row_bounds_rechecked=True,
        row_module_uses_only_unit_constant_divisions=True,
        generated_inputs_sha256={'core_matrix': sha(core_matrix), 'univariate': sha(univariate)},
        generated_outputs_sha256={'selection': sha(saved_selection), 'native_certificate': sha(saved_output)},
        input_sha256={str(path.relative_to(ROOT)): sha(path) for path in paths},
        seconds=time.time()-start,
    )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--core-matrix', type=Path, required=True)
    parser.add_argument('--univariate', type=Path, required=True)
    parser.add_argument('--workdir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run_minor_module(args.workdir, args.core_matrix, args.univariate)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({key: result[key] for key in [
        'all_passed', 'determinant_degree', 'determinant_leading_coefficient',
        'module_rank', 'module_degree', 'seconds']}))


if __name__ == '__main__':
    main()

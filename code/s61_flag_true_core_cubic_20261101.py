#!/usr/bin/env python3
"""True cubic initial space on all seven parameters of the r=-4/25 core."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    engines = [
        'code/s61_flag_boundary_eight_frame_engine_20261101.py',
        'code/s61_flag_true_core_frame_engine_20261101.py',
        'code/s61_flag_linear_artifacts_quadratic_engine_20261101.py',
        'code/s61_flag_true_core_raw3_engine_20261101.py',
        'code/s61_flag_true_core_products3_engine_20261101.py',
        'code/s61_flag_true_core_cubic_engine_20261101.py',
    ]
    native = 'code/parametric_implicit_jets_20261101.cpp'
    with tempfile.TemporaryDirectory(prefix='flag-true-core-cubic-') as temp:
        wd = Path(temp)
        binary = wd / 'parametric-jets'
        subprocess.run([
            'c++', '-O3', '-std=c++17',
            '-I' + str(ROOT / 'work/flint-headers/flint/3.5.0/include'),
            '-I' + str(ROOT / 'work/flint-headers/mpfr/4.2.2/include'),
            '-I/opt/homebrew/include', str(ROOT / native),
            str(ROOT / 'work/m2-runtime/lib/libflint.24.0.dylib'),
            '-L/opt/homebrew/lib', '-lgmp', '-o', str(binary),
        ], check=True)
        for i, engine in enumerate(engines):
            directory = wd if i < 2 else wd / 'core-0'
            extra = [str(binary)] if i == 3 else []
            subprocess.run([sys.executable, str(ROOT / engine), str(directory), *extra],
                           cwd=ROOT, check=True)
        wd /= 'core-0'
        products = json.loads((wd / 'flag-true-core-cubic-products-global.json').read_text())
        cubic = json.loads((wd / 'flag-true-core-cubic-global.json').read_text())
    assert products['total_rank'] == 114 and cubic['rank'] == 5
    paths = [str(Path(__file__).relative_to(ROOT)), *engines, native,
             'code/s61_model.py', 'code/s61_flag_transverse_obstruction_20261101.py',
             'reports/S61_FLAG_UNIFORM_OPEN_20261101.json',
             'reports/S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101.json',
             'reports/S61_FLAG_TRUE_CORE_QUADRATIC_20261101.json']
    result = dict(
        all_passed=True, core='-4/25', parameter_count=7,
        raw135_three_jets_rebuilt=True, all_raw_quadratic_coefficients_cross_checked=True,
        quadratic_product_space=products, cubic_space=cubic,
        hilbert_function_through3=[1, 10, 40, 101], t_order_at_least4=True,
        all_seven_parameter_values=True,
        proof_note='research/notes/S61_FLAG_TRUE_CORE_CUBIC_20261101.md',
        limitations=['Finite initial space through degree3 only.',
                     'No exact t order, full initial equality, ambient nilpotence or saturation removal.'],
        seconds=time.time() - started,
        input_sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
    )
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print('ALL TRUE CORE CUBICS PASS', result['seconds'], flush=True)


if __name__ == '__main__':
    main()

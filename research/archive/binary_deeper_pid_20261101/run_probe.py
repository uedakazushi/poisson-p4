#!/usr/bin/env python3
"""Reconstruct the WIP core slice computation. This is not an exhaustion proof."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'code'))
from s61_binary_deeper_cubic_nuisance_20261101 import prepare, CERTIFICATE, LINEAR
from s61_binary_deeper_quartic_nuisance_20261101 import prepare_quartic


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output', required=True)
    p.add_argument('--work-dir', help='An absent /tmp directory; retain large reconstructed caches there.')
    p.add_argument('--exact-inputs', action='store_true', help='Also reconstruct the large QQ matrix and its 55 e-only pivots.')
    args = p.parse_args()
    output = Path(args.output)
    assert not output.exists()
    started = time.time()
    owner = None
    if args.work_dir:
        tmp = Path(args.work_dir)
        assert str(tmp).startswith('/tmp/') and not tmp.exists()
        tmp.mkdir()
    else:
        owner = tempfile.TemporaryDirectory(prefix='poisson_deeper_pid_archive_', dir='/tmp')
        tmp = Path(owner.name)
    env = os.environ.copy()
    env['DYLD_LIBRARY_PATH'] = str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'
    lib = ROOT/'work/m2-runtime/lib/libflint.24.0.dylib'
    def compile_source(name, source):
        binary = tmp/name
        subprocess.run(['clang++','-O3','-std=c++17',
            '-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),
            '-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),
            '-I/opt/homebrew/include', str(source), str(lib),
            '-L/opt/homebrew/lib','-lgmp','-o',str(binary)], check=True)
        return binary
    sources = {name:HERE/(name+'.cpp') for name in
        ['all_ground_quartics','ground_mod','matrix_input','line_popov','popov_support']}
    for name in ['core_inverse','nuisance_dag','unit_relations']:
        sources[name] = ROOT/'code'/f's61_binary_deeper_quartic_{name}_20261101.cpp'
    sources['cubic_products'] = ROOT/'code/s61_binary_deeper_cubic_products_20261101.cpp'
    if args.exact_inputs:
        sources.update({name:HERE/(name+'.cpp') for name in
            ['ground_exact','ground_denominators','content_units','exact_units']})
    binaries = {name:compile_source(name, source) for name,source in sources.items()}
    def run(name, *paths):
        with (tmp/(name+'.log')).open('w') as log:
            subprocess.run([str(binaries[name]), *map(str,paths)], env=env,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        print(name, round(time.time()-started,3), flush=True)
    data = json.loads(LINEAR.read_text())['linear']
    certificate = json.loads(CERTIFICATE.read_text())
    quad,_ = prepare(data,certificate,tmp/'cubic_input.txt',tmp/'module3.txt')
    prepare_quartic(data,certificate,quad,tmp)
    run('cubic_products',tmp/'cubic_input.txt',tmp/'cubics.txt')
    run('core_inverse',tmp/'core.txt',tmp/'inverse_core.txt')
    run('nuisance_dag',tmp/'outer.txt',tmp/'inverse_core.txt',tmp/'inverse.txt')
    run('unit_relations',tmp/'module3.txt',tmp/'relations.txt')
    run('all_ground_quartics',tmp/'operators.txt',tmp/'inverse.txt',tmp/'relations.txt',tmp/'quartics.txt')
    run('all_ground_quartics',tmp/'operators.txt',tmp/'inverse.txt',HERE/'extra_syzygy.txt',tmp/'extra_quartic.txt')
    run('ground_mod',tmp/'cubics.txt',tmp/'cubic_mod.txt')
    run('ground_mod',tmp/'quartics.txt',tmp/'quartic_mod.txt')
    run('matrix_input',HERE/'seven_quadrics.txt',tmp/'cubic_mod.txt',tmp/'quartic_mod.txt',tmp/'matrix_mod.txt')
    # Matrix_input is an exact rational parser; the following stage reduces
    # its rational quadratic entries and modular higher entries modulo32003.
    run('line_popov',tmp/'matrix_mod.txt',tmp/'line_basis.txt')
    run('popov_support',tmp/'line_basis.txt',HERE/'horizontal_factors.txt',tmp/'line_support.txt')
    support_log = (tmp/'popov_support.log').read_text()
    assert 'DET 430\n' in support_log and 'REMAIN 0\n' in support_log
    assert (tmp/'extra_quartic.txt').read_text().split()[2] == '0'
    if args.exact_inputs:
        run('ground_exact',tmp/'cubics.txt',tmp/'cubic_exact.txt')
        run('ground_exact',tmp/'quartics.txt',tmp/'quartic_exact.txt')
        run('ground_denominators',HERE/'horizontal_factors.txt',tmp/'quartics.txt',tmp/'denominators.txt')
        run('content_units',HERE/'horizontal_factors.txt',tmp/'quartic_exact.txt.contents')
        run('matrix_input',HERE/'seven_quadrics.txt',tmp/'cubic_exact.txt',tmp/'quartic_exact.txt',tmp/'matrix_exact.txt')
        run('exact_units',tmp/'matrix_exact.txt',tmp/'exact_units.txt')
    inputs = sorted(set(sources.values()) | set(HERE.glob('*.hpp')) |
        {HERE/'seven_quadrics.txt',HERE/'horizontal_factors.txt',HERE/'extra_syzygy.txt',Path(__file__),CERTIFICATE,LINEAR})
    report = dict(status='WIP: finite-field one-core slice only',
        prime=32003, fixed_e=2, core_polynomial_variable='d',
        matrix_shape=[266,210], line_rank=210, line_cokernel_length=430,
        known_factor_removed_degrees=[0,42,42,1,2,28,315], remaining_support_degree=0,
        limitations=['This does not prove the QQ(e)[d] statement.',
                     'No new core-stratum or full S61 exhaustion is concluded.',
                     'Generic and exact division-Popov sources are unfinished exploration.'],
        input_sha256={str(x.relative_to(ROOT)):sha(x) for x in inputs},
        large_output_sha256={x.name:sha(x) for x in tmp.glob('*.txt')}, seconds=time.time()-started)
    with output.open('x') as f: json.dump(report,f,indent=2); f.write('\n')
    if owner: owner.cleanup()


if __name__ == '__main__': main()

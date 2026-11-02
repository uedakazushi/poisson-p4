#!/usr/bin/env python3
"""Exact d=3 slice surjectivity and exclusion of unknown vertical core curves.

Regenerates the weighted matrix, proves one minor's degree/known valuations by
constant-QQ kernel peeling, and uses finite Z_(32003)-module Nakayama. The d=4
finite-field full row module excludes the remaining moving quadratic vertical
curve. Isolated core points and S61 saturation are not concluded here.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
import sympy as s
from s61_binary_weighted_seventeen_curves_20261102 import prepare

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / 'research/notes/S61_BINARY_WEIGHTED_VERTICAL_SUPPORT_20261102.md'
NATIVES = [
    's61_binary_vertical_delta_export_20261102',
    's61_binary_constant_kernel_peel_20261102',
    's61_binary_constant_kernel_integer_peel_20261102',
    's61_binary_vertical_local_export_20261102',
    's61_binary_vertical_modular_rowmodule_20261102',
    's61_binary_vertical_determinant_degree_20261102',
    'parametric_implicit_jets_20261101',
    's61_binary_deeper_weighted_extract_20261102',
]
PRIME = 32003


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def peel_result(prefix):
    lines = Path(str(prefix) + '.meta').read_text().splitlines()
    stages = [list(map(int, line.split())) for line in lines if line[0].isdigit()]
    total = int(next(line.split()[1] for line in lines if line.startswith('TOTAL ')))
    assert stages[-1][2] == 0
    assert sum(row[2] for row in stages) == total
    return dict(stages=stages, valuation=total,
                polynomial_kernel_identities_checked=True,
                all_delta_divisions_exact=True,
                final_full_rank=True,
                kernels_sha256=[sha(str(prefix) + '.kernel' + str(i))
                                for i in range(len(stages)-1)],
                final_matrix_sha256=sha(str(prefix) + '.full'))


def modpoly(line, variable):
    data = list(map(int, line.split()))
    assert data[0] == len(data)-1
    return s.Poly.from_list(data[1:][::-1], variable, modulus=PRIME)


def verify(workdir):
    start = time.time()
    wd = Path(workdir).resolve()
    wd.mkdir(parents=True, exist_ok=True)
    library = ROOT / 'work/m2-runtime/lib'
    env = dict(os.environ, DYLD_LIBRARY_PATH=str(library)+':/opt/homebrew/lib')
    flags = ['clang++', '-O3', '-std=c++17',
             '-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),
             '-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),
             '-I/opt/homebrew/include']
    libs = [str(library/'libflint.24.0.dylib'), '-L/opt/homebrew/lib', '-lgmp']
    sources = [ROOT/'code'/f'{name}.cpp' for name in NATIVES]
    binaries = {}
    runs = []
    with tempfile.TemporaryDirectory(prefix='binary-vertical-bin-') as temp:
        for source in sources:
            target = Path(temp)/source.stem
            result = subprocess.run(flags+[str(source)]+libs+['-o', str(target)],
                                    capture_output=True, text=True)
            if result.returncode:
                raise RuntimeError(result.stderr)
            binaries[source.stem] = target

        def run(name, args, timeout=180):
            log = wd/f'run{len(runs):02d}_{name}.log'
            begin = time.time()
            with log.open('x') as stream:
                subprocess.run([str(binaries[name]), *map(str, args)], env=env,
                               stdout=stream, stderr=stream, check=True,
                               timeout=timeout)
            runs.append(dict(name=name, seconds=time.time()-begin,
                             log_sha256=sha(log)))
            return log

        inp, raw, rows, meta, weights, coeff, known = prepare(wd)
        knownfile = wd/'known13'
        with knownfile.open('x') as stream:
            stream.write(str(len(known))+'\n')
            for poly in known:
                stream.write(str(len(poly))+'\n')
                for monomial, coefficient in poly.items():
                    stream.write(str(coefficient)+' '+' '.join(map(str, monomial))+'\n')
        assert len(known) == 13
        run('parametric_implicit_jets_20261101', [inp, raw])
        run('s61_binary_deeper_weighted_extract_20261102', [raw, coeff, rows, meta, weights])
        delta = wd/'delta'
        exportlog = run('s61_binary_vertical_delta_export_20261102',
                        [rows, 3, knownfile, delta])
        assert delta.read_text().splitlines()[0] == 'QQ_DELTA_MINOR 269 HUNGARIAN 1607'
        assert exportlog.read_text().count('REMOVED 0 LEFT 0') == 61
        infinity = wd/'infinity'
        run('s61_binary_constant_kernel_peel_20261102',
            [str(delta)+'.full', 20, infinity])
        infinity_result = peel_result(infinity)
        assert infinity_result['valuation'] == 366
        assert infinity_result['stages'][-1][1] == 269
        local_results = {}
        for name, root, expected in [
                ('r7_24', '7/24', 70), ('r1_10', '1/10', 113),
                ('r1_2', '1/2', 176), ('rm8_15', '-8/15', 197),
                ('quad', 'quad', 184)]:
            local = wd/name
            exportlog = run('s61_binary_vertical_local_export_20261102',
                            [str(delta)+'.full', delta, root, local])
            prefix = wd/(name+'_peel')
            helper = ('s61_binary_constant_kernel_integer_peel_20261102'
                      if name == 'quad' else 's61_binary_constant_kernel_peel_20261102')
            run(helper, [local, 25, prefix])
            result = peel_result(prefix)
            assert result['valuation'] == expected
            assert result['stages'][-1][1] == (538 if name == 'quad' else 269)
            if name == 'quad':
                normline = next(line for line in exportlog.read_text().splitlines()
                                if line.startswith('NORM_POINT_CHECK '))
                values = list(map(int, normline.split()[1:]))
                assert values[0] == values[1]*values[2] % PRIME
                result['independent_modular_norm_evaluation'] = values
                result['original_quadratic_multiplicity'] = expected//2
            local_results[name] = result

        e, delta_symbol = s.symbols('e delta')
        f3 = 270*e**2+321*e-388
        assert s.discriminant(f3, e) == 522081
        assert 722**2 < 522081 < 723**2
        T = s.Matrix([[0, (delta_symbol+388)/270], [1, -s.Rational(321, 270)]])
        assert (270*T*T+321*T-388*s.eye(2)-delta_symbol*s.eye(2)).applyfunc(s.expand) == s.zeros(2)
        det = wd/'mod_minor'
        run('s61_binary_vertical_determinant_degree_20261102', [rows, 3, det])
        detlines = det.read_text().splitlines()
        assert detlines[0] == '3 32003 1607 1241'
        assert detlines[1] == delta.read_text().splitlines()[1]
        P = modpoly(detlines[2], e)
        known5 = [24*e-7, 10*e-1, 2*e-1, 15*e+8, f3]
        powers = [70, 113, 176, 197, 92]
        monics = [s.Poly(poly, e, modulus=PRIME).monic() for poly in known5]
        U = s.Poly(1, e, modulus=PRIME)
        D = s.Poly(1, e, modulus=PRIME)
        for factor, exponent in zip(monics, powers):
            U *= factor**exponent
            D *= factor
        Q, remainder = s.div(P.monic(), U)
        assert remainder.is_zero and Q.degree() == 501
        assert s.gcd(Q, D).degree() == 0
        full_results = {}
        full_polynomials = []
        for d in [3, 4]:
            full = wd/f'full_d{d}'
            run('s61_binary_vertical_modular_rowmodule_20261102', [rows, d, full])
            lines = full.read_text().splitlines()
            header = list(map(int, lines[0].split()))
            assert header[:4] == [d, PRIME, 269, 660]
            H = modpoly(lines[1], e).monic()
            moving = 90*d*e**2+(600-93*d)*e+(4*d-400)
            expected = (monics[0]**65 * monics[1]**113 * monics[2]**138 *
                        monics[3]**160 * s.Poly(moving, e, modulus=PRIME).monic()**92)
            assert H == expected
            full_results[str(d)] = dict(rank=269, determinant_degree=660,
                                        operations=header[4], output_sha256=sha(full))
            full_polynomials.append(H)
        common = s.gcd(*full_polynomials)
        assert common == monics[0]**65*monics[1]**113*monics[2]**138*monics[3]**160
        assert common.degree() == 476
        assert s.gcd(monics[-1], full_polynomials[1]).degree() == 0

    dependencies = [Path(__file__), NOTE, *sources,
                    ROOT/'code/s61_binary_deeper_quartic_field_20261101.hpp',
                    ROOT/'code/s61_binary_weighted_seventeen_curves_20261102.py',
                    ROOT/'code/s61_model.py',
                    ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
                    ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json',
                    ROOT/'reports/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101.json',
                    ROOT/'reports/S61_BINARY_DEEPER_WEIGHTED_FINITE_MODEL_20261102.json']
    return dict(id='S61_BINARY_WEIGHTED_VERTICAL_SUPPORT_20261102', all_passed=True,
                matrix_shape=[390, 269], fixed_d=3, prime=PRIME,
                tropical_degree_upper=1607, infinity=infinity_result,
                exact_minor_degree=1241, local=local_results,
                known_factors=list(map(str, known5)), exact_multiplicities=powers,
                residual_monic_degree=501,
                residual_coprime_to_known5_mod_prime=True,
                full_modular_rowmodules=full_results,
                modular_d3_d4_common_degree=476,
                moving_quadratic_coprime_to_d4_module=True,
                slice_surjectivity='QQ[e,((24e-7)(10e-1)(2e-1)(15e+8)(270e^2+321e-388))^-1]',
                possible_vertical_curve_roots=['7/24', '1/10', '1/2', '-8/15'],
                proof_method='Exact QQ constant-kernel peeling; finite Z_(32003) Nakayama; primitive quadratic divisibility versus d4 finite-field full row module.',
                limitations=['No assertion of surjectivity at every two-parameter core point.',
                             'Isolated core rank-drop points are not excluded.',
                             'Nuisance-zero weighted matrix only; no full S61 saturation or exhaustion.'],
                input_sha256={str(path.relative_to(ROOT)): sha(path) for path in dependencies},
                rebuilt_sha256={'weighted_rows': sha(rows), 'delta_matrix': sha(str(delta)+'.full')},
                runs=runs, seconds=time.time()-start)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workdir', type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    if args.workdir:
        result = verify(args.workdir)
    else:
        with tempfile.TemporaryDirectory(prefix='binary-vertical-') as temp:
            result = verify(temp)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({key: result[key] for key in ['all_passed', 'exact_minor_degree',
                     'exact_multiplicities', 'residual_monic_degree', 'seconds']}))

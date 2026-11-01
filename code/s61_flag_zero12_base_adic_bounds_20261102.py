#!/usr/bin/env python3
"""Finite base-adic cutoff from certified QQ initial relations, not full gr equality."""
import argparse
import hashlib
import json
import os
import subprocess
import tempfile
import time
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def verify():
    started = time.monotonic()
    names = [
        'reports/S61_FLAG_ZERO12_RATIONAL_INITIAL_20261101.json',
        'reports/S61_FLAG_ZERO12_LOW_JETS_20261101.json',
        'reports/S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101.json',
        'research/notes/FINITE_BASE_ADIC_CUTOFF_20261102.md',
        'research/notes/S61_FLAG_ZERO12_BASE_ADIC_BOUNDS_20261102.md',
    ]
    initial, low, quotient = [json.loads((ROOT/n).read_text()) for n in names[:3]]
    assert initial['all_passed'] and low['all_passed'] and quotient['all_passed']
    assert initial['field'] == 'QQ' and initial['geometry']['parameters'] == [0, 1, 3]
    assert low['t_order_at_least'] == 4
    relations = initial['initial_relations']
    assert len(relations) == 64
    expr = []
    for row in relations:
        assert len({len(m) for m, c in row}) == 1
        expr.append('+'.join('('+c+')*'+'*'.join('v_'+str(i) for i in m) for m, c in row))
    expected = [
        [1, 6, 13, 15, 10, 6],
        [1, 9, 31, 52, 55, 31, 14],
        [1, 9, 37, 88, 127, 110, 57, 22],
        [1, 9, 37, 98, 187, 220, 172, 83, 30],
    ]
    program = 'R=QQ[v_0..v_8];\nJ=ideal('+','.join(expr)+');\na=ideal(v_0,v_1,v_3);\n'
    for power, hf in enumerate(expected, 1):
        program += f'K=J+a^{power};\n'
        program += 'assert(apply(toList(0..'+str(len(hf)-1)+'),d->hilbertFunction(d,R/K))=={'+','.join(map(str,hf))+'});\n'
        program += f'assert(dim(R/K)==0); assert(degree(R/K)=={sum(hf)});\n'
        program += f'assert(gens((ideal vars R)^{len(hf)})%K==0);\n'
    program += 'print "BASE ADIC BOUNDS PASS";\n'
    env = os.environ.copy()
    env['DYLD_LIBRARY_PATH'] = str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
    with tempfile.TemporaryDirectory(prefix='zero12-base-adic-bounds-') as temp:
        source = Path(temp)/'verify.m2'
        source.write_text(program)
        run = subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary', '--script', str(source)],
                             cwd=ROOT, env=env, text=True, capture_output=True, timeout=30)
        assert run.returncode == 0 and 'BASE ADIC BOUNDS PASS' in run.stdout, run.stdout+run.stderr
    paths = [Path(__file__)] + [ROOT/n for n in names]
    return dict(
        id='S61_FLAG_ZERO12_BASE_ADIC_BOUNDS_20261102', all_passed=True, field='QQ',
        base_variables=[0, 1, 3], point61=initial['point61'],
        homogeneous_initial_subideal_generator_count=64,
        homogeneous_base_quotient_hilbert_functions=expected,
        homogeneous_base_quotient_lengths=list(map(sum, expected)),
        graded_module_generating_degree_bound=5,
        all_base_powers_cutoff='m^(r+5) S is contained in (v0,v1,v3)^r S for every r>=1',
        exact_finite_computation_degree='r+4 suffices for S/(v0,v1,v3)^r S',
        actual_t_square_base_order_at_least=3,
        kernel_t_square_base_order_at_least=3,
        weighted_first_quotient_column_bounds={str(n):sum(comb(2+b,b)*comb(n+6*(2-b)-1,n) for b in range(2)) for n in (6,7)},
        primary_raw_initial_replayed=False,
        complete_initial_ideal_identified=False, ambient_t_nilpotence_proved=False,
        scope='Specified zero12 fixed-beta point. No actual quotient length equality and no ambient removal.',
        m2_program=program, m2_output=run.stdout+run.stderr,
        proof_note=names[-1],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.monotonic()-started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify()
    with args.output.open('x') as out:
        json.dump(result, out, indent=2)
        out.write('\n')
    print(json.dumps(dict(all_passed=True, seconds=result['seconds'])))

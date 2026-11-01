#!/usr/bin/env python3
"""Independent homogeneous check plus explicit scope of the seventh-step source audit."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]

def verify():
    started = time.time()
    primary = ROOT/'reports/S61_FLAG_ZERO12_SEVENTH_INITIAL_20261101.json'
    initial = ROOT/'reports/S61_FLAG_ZERO12_RATIONAL_INITIAL_20261101.json'
    saved = json.loads(primary.read_text())
    previous = json.loads(initial.read_text())
    assert saved['all_passed'] and previous['all_passed']
    for record in (saved, previous):
        for path, digest in record['input_sha256'].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
    relations = previous['initial_relations']
    assert len(relations) == 64
    assert all(len({len(mon) for mon, coefficient in row if coefficient != '0'}) == 1 for row in relations)
    assert {len(row[0][0]) for row in relations} == {2, 3, 4, 5, 6}
    expressions = ['+'.join('('+str(coefficient)+')'+''.join('*v_'+str(i) for i in mon)
                            for mon, coefficient in row) for row in relations]
    expected = [1, 9, 37, 98, 202, 329, 471, 618]
    program = 'R=QQ[v_0..v_8]; I=ideal('+','.join(expressions)+');\n'
    program += 'G=gb I; HF=apply(toList(0..7),d->hilbertFunction(d,R/I));\n'
    program += 'assert(HF=={'+','.join(map(str, expected))+'});print("INDEPENDENT FULL NINE VARIABLE HF PASS");\n'
    environment = os.environ.copy()
    environment['DYLD_LIBRARY_PATH'] = str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
    with tempfile.TemporaryDirectory(prefix='zero12-seventh-review-') as temp:
        source = Path(temp)/'independent.m2'
        source.write_text(program)
        result = subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary', '--script', str(source)],
                                cwd=ROOT, env=environment, text=True, capture_output=True, check=True, timeout=60)
        assert 'INDEPENDENT FULL NINE VARIABLE HF PASS' in result.stdout
    assert saved['hilbert_function'] == expected
    assert sum(expected) == saved['quotient_by_m8_length'] == 1765
    assert 135+9*saved['lower_trace_rank'] == saved['all_closure_checks'] == 34857
    assert saved['initial_rank_through_seven'] == saved['lower_trace_rank']+6435-618 == 9675
    assert saved['new_seventh_initial_directions'] == 0
    assert saved['full_initial_ideal_identified'] is False and saved['homogenizer_nilpotence_proved'] is False
    paths = [Path(__file__), primary, initial, ROOT/'code/s61_flag_zero12_seventh_initial_20261101.py',
             ROOT/'code/projected_next_initial_replay_20261101.cpp',
             ROOT/'code/s61_flag_zero12_initial_module_20261101.py',
             ROOT/'research/notes/S61_FLAG_ZERO12_SEVENTH_INITIAL_20261101.md']
    return dict(id='S61_FLAG_ZERO12_SEVENTH_INITIAL_REVIEW_20261101', all_passed=True,
                independent_homogeneous_computation='QQ Groebner basis in the full nine-variable polynomial ring; no module projection helper imported.',
                hilbert_function_through_seven=expected, quotient_by_m8_length=1765,
                proof_source_audit=['The projection changes only degree seven, by J7 already in the true ideal modulo m8.',
                                    'Every replay row is constructed from an original relation using earlier rows and variable multiples.',
                                    'Adaptive exact leading-pivot verification catches modularly vanishing coefficients.',
                                    'All original rows and every variable multiple of each basis row reduce to zero, proving complete closure.'],
                not_independently_recomputed=['Original 135 seventh-order implicit jets.', 'All 3858 exact filtered replay rows and 34857 closure reductions.'],
                accepted_scope='True initial degrees through seven at the specified point; no complete initial ideal or homogenizer nilpotence.',
                input_sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
                seconds=time.time()-started)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    report = verify()
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print('ZERO12 SEVENTH REVIEW PASS', report['seconds'])

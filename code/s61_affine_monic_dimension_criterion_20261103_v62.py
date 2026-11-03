#!/usr/bin/env python3
"""Six-parameter exact UNIT criterion for generic affine dimension zero.

The monic coefficient polynomial detects every positive-dimensional component.
This producer checks the complete input and small exact traces; it does not
claim that the full S61 ideal has been solved.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    previous = ROOT / 'reports/S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_v62'
    saved = json.loads((previous / 'report.json').read_text())
    assert saved['all_passed']
    for name, value in saved['input_sha256'].items():
        assert sha(ROOT / name) == value
    for name, value in saved['evidence_sha256'].items():
        assert sha(previous / name) == value
    source = (previous / 'wronski-exact.sing').read_text()
    rows = source.split('ideal I=', 1)[1].split(';\n', 1)[0].split(',')
    assert len(rows) == 136
    rows = rows[:134]
    assert all('zz' not in row and 'cc' not in row for row in rows)
    with gzip.open(previous / 'restricted-original134.json.gz', 'rt') as stream:
        encoded = json.load(stream)
    assert len(encoded) == 134
    counts = Counter()
    for row in encoded:
        for monomial, _ in row:
            counts.update(set(monomial))
    # Only the monic detector's order changes. No coordinate is cut out.
    eliminated = min(range(45), key=lambda i: (counts[i], i))
    assert eliminated == 32
    order = [eliminated] + [i for i in range(45) if i != eliminated]
    assert sorted(order) == list(range(45))
    elimination = 'cc^45-' + '-'.join(f'cc^{j}*x{i}' for j, i in enumerate(order) if j)
    variables = ','.join(f'x{i}' for i in range(45) if i != eliminated)
    parameters = ','.join([*(f'ss{i}' for i in range(1, 6)), 'cc'])
    m2 = f'A=QQ[{parameters}];K=frac A;R=K[{variables}];\nx{eliminated}={elimination};\n'
    m2 += 'I=ideal(' + ','.join(rows) + ');\n'
    m2 += 'print "MONIC6_FULL_EXACT_START";stdio << flush;G=gb(I,ChangeMatrix=>true);print "MONIC6_FULL_EXACT_GB_DONE";stdio << flush;\n'
    m2 += 'if 1_R%G==0 then (C=matrix{{1_R}}//(gens I);assert((gens I)*C==matrix{{1_R}});print "MONIC6_EXACT_UNIT_COFACTOR_PASS";stdio << flush;) else (print("NONUNIT_DIM="|toString dim(R/I));stdio << flush;);exit 0;\n'
    (args.output / 'monic-six-exact.m2').write_text(m2)
    singular = f'ring r=(0,{parameters}),({variables}),dp;\npoly x{eliminated}={elimination};\nideal I=' + ','.join(rows) + ';\n'
    singular += 'print("MONIC6_FULL_EXACT_START");ideal G=slimgb(I);print("MONIC6_FULL_EXACT_GB_DONE");print("DIM="+string(dim(G)));print("VDIM="+string(vdim(G)));\n'
    singular += 'if(size(G)==1 && G[1]==1){matrix C=lift(I,G);matrix M=matrix(I)*C;if(M[1,1]!=1){ERROR("UNIT_IDENTITY_FAILED");quit;}print("MONIC6_EXACT_UNIT_COFACTOR_PASS");}quit;\n'
    (args.output / 'monic-six-exact.sing').write_text(singular)
    toy = f'A=QQ[{parameters}];K=frac A;R=K[x,y];\n'
    toy += 'F=cc^2-x*cc-y;I=ideal(x*x-2,y-3,F);G=gb(I,ChangeMatrix=>true);assert(1_R%G==0);C=matrix{{1_R}}//gens I;assert((gens I)*C==matrix{{1_R}});\n'
    toy += 'J=ideal(x*y-1,F);H=gb(J,ChangeMatrix=>true);assert(1_R%H!=0);D=(gens H)//gens J;assert((gens J)*D==gens H);assert(numColumns basis(R/J)==2);\n'
    toy += 'N=ideal(x*x,y,F);P=gb(N,ChangeMatrix=>true);assert(1_R%P==0);E=matrix{{1_R}}//gens N;assert((gens N)*E==matrix{{1_R}});\n'
    toy += 'print "MONIC_DETECTOR_FINITE_CURVE_AND_NILPOTENT_TRACES_PASS";stdio << flush;exit 0;\n'
    probe = args.output / 'compatibility-probe.m2'
    log = args.output / 'compatibility-probe.log'
    probe.write_text(toy)
    with log.open('w') as stream:
        subprocess.run([shutil.which('M2'), '--script', str(probe.resolve())], stdout=stream,
                       stderr=subprocess.STDOUT, timeout=30, check=True)
    assert 'MONIC_DETECTOR_FINITE_CURVE_AND_NILPOTENT_TRACES_PASS' in log.read_text()
    assert 'error:' not in log.read_text()
    data = dict(original_equations=134, original_coordinates=45,
                monic_degree=45, coefficient_coordinate_order=order,
                eliminated_coordinate=eliminated, eliminated_coefficient=-1,
                substitution=elimination, coordinate_occurrence_terms=dict(sorted(counts.items())),
                coefficient_field='QQ(ss1,...,ss5,cc)', unknowns_after_exact_linear_elimination=44,
                exact_unit_equivalent_to='generic full affine fiber has dimension<=0',
                fixed_pair_or_center_finiteness_assumption=False,
                algebraic_constant_extension=False, generic_to_numeric_specialization=False)
    (args.output / 'monic-model.json').write_text(json.dumps(data, indent=2) + '\n')
    inputs = [Path(__file__), previous / 'report.json', previous / 'wronski-exact.sing',
              previous / 'restricted-original134.json.gz',
              ROOT / 'research/notes/S61_AFFINE_MONIC_DIMENSION_CRITERION_20261103_v62.md']
    report = dict(id='S61_AFFINE_MONIC_DIMENSION_CRITERION_20261103_v62', all_passed=True,
                  transcendental_parameter_count=6, original_equations_preserved=134,
                  exact_unit_criterion_proved=True, small_exact_unit_and_nonunit_traces_passed=True,
                  actual_full_unit_result_obtained=False, full_affine_dimension_zero_claimed=False,
                  five_point_exhaustion_claimed=False,
                  limitations=['No full134-equation GB or UNIT result is computed by this producer.',
                               'Even a future full UNIT result proves dimension0; it does not compute the total affine length or exhaust all points.'],
                  input_sha256={str(p.relative_to(ROOT)): sha(p) for p in inputs},
                  evidence_sha256={p.name: sha(p) for p in args.output.iterdir() if p.is_file()},
                  seconds=time.monotonic() - started)
    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print('PASS six-parameter monic detector; all134 equations; eliminated q32 exactly;', report['seconds'])


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Bounded native FLINT expansion of all frozen 350 weighted residuals."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
PRIME = 32003


def residue(value):
    value = Fraction(value)
    assert value.denominator % PRIME
    return value.numerator * pow(value.denominator, -1, PRIME) % PRIME


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', type=Path, default=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz')
    parser.add_argument('--a', default='2/9')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--timeout', type=int, default=300)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    raw = gzip.decompress(args.model.read_bytes())
    model = json.loads(raw)
    assert len(model['nodes']) == 9647
    a = residue(args.a)
    lines = [f"{len(model['nodes'])} 34 350 {PRIME} {a}",
             ' '.join(map(str, model['input_weights'][1:])),
             ' '.join(map(str, model['residual_weight_bounds'])),
             ' '.join(map(str, model['all350_residual_nodes']))]
    point = [a] + [(37*j*j+19*j+11) % PRIME for j in range(34)]
    expected = []
    for node in model['nodes']:
        op = node[0]
        if op == 'q':
            value = residue(node[1]); lines.append(f'q {value}')
        elif op == 'v':
            value = point[node[1]]; lines.append(f'v {node[1]}')
        elif op == 'inv':
            value = pow(expected[node[1]], -1, PRIME); lines.append(f'i 1 {node[1]}')
        else:
            value = 0 if op == '+' else 1
            for ref in node[1]:
                value = (value + expected[ref] if op == '+' else value * expected[ref]) % PRIME
            lines.append(f"{op} {len(node[1])} " + ' '.join(map(str, node[1])))
        expected.append(value)
    data = args.output/'dag.in'; data.write_text('\n'.join(lines)+'\n')
    native = ROOT/'code/q_weighted_residual_expansion_20261103_v63.cpp'
    executable = args.output/'expand'
    flags = shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],capture_output=True,text=True,check=True).stdout)
    subprocess.run(['c++','-O3','-std=c++17',str(native),*flags,'-o',str(executable)],check=True)
    prefix = args.output/'residuals'
    status = 'completed'
    with (args.output/'native.log').open('w') as log:
        try:
            completed = subprocess.run([str(executable.resolve()),str(data.resolve()),str(prefix.resolve())],stdout=log,stderr=subprocess.STDOUT,timeout=args.timeout)
            if completed.returncode:
                status = f'exit_{completed.returncode}'
        except subprocess.TimeoutExpired:
            status = 'timeout'
    report = dict(id='Q_WEIGHTED_RESIDUAL_EXPANSION_20261103_V63',
                  expansion_completed=status=='completed',status=status,
                  coefficient_field=f'GF({PRIME})',a_rational=args.a,a_mod_prime=a,
                  memory_limit_bytes=6000000000,wall_timeout_seconds=args.timeout,
                  positive_variables=34,all350_equations_retained=True,
                  native_reference_count_freeing=True,
                  finite_field_result_only=True,characteristic_zero_dimension_claimed=False,
                  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__).resolve(),native,args.model.resolve()]})
    if status == 'completed':
        summary=Path(str(prefix)+'.summary').read_text().splitlines()
        lengths=list(map(int,summary[0].split()[1:]))
        actual=list(map(int,summary[1].split()[1:]))
        assert len(lengths)==len(actual)==350
        assert actual==[expected[node] for node in model['all350_residual_nodes']]
        selected=[i for block in model['implicit_blocks'] for i in block['equation_rows']]
        assert len(selected)==74 and all(lengths[i]==0 for i in selected)
        tokens=summary[2].split()
        stats={key:float(value) if key=='SECONDS' else int(value) for key,value in zip(tokens[::2],tokens[1::2])}
        report.update(term_counts=lengths,nonzero_rows=sum(bool(n) for n in lengths),
                      all350_expansions_match_independent_DAG_evaluation=True,
                      all74_selected_polynomials_identically_zero=True,
                      native_statistics=stats,
                      evidence_sha256={p.name:sha(p) for p in [Path(str(prefix)+s) for s in ['.sparse','.m2','.summary']]})
    report['seconds']=time.monotonic()-started
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['term_counts','input_sha256','evidence_sha256']},indent=2))


if __name__=='__main__':
    main()

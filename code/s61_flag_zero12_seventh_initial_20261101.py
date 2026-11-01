#!/usr/bin/env python3
"""Check the actual seventh initial piece at the saved ZERO12 flag point.

Only degree seven is projected by the previously certified initial subideal.
The original six lower degrees and every filtered tail are retained exactly.
"""
import argparse
import hashlib
import json
import os
import subprocess
import tempfile
import time
from itertools import combinations_with_replacement as cwr
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_flag_zero12_initial_module_20261101 import module_program
from truncated_local_ideal_io_20261101 import write_input

ROOT = Path(__file__).resolve().parents[1]


def run(command, **kwargs):
    result = subprocess.run(list(map(str, command)), cwd=ROOT, capture_output=True,
                            text=True, check=True, **kwargs)
    return result


def projection(saved, directory, monomials):
    program, evidence = module_program(saved['initial_relations'])
    program = '\n'.join(program.splitlines()[:2]).replace(
        'S=QQ[x_0..x_3]', 'S=QQ[z,x_0..x_3]')+'\n'
    output = directory/'module-projection.txt'
    program += 'G=gb N;raw=basis(7,F);red=raw%G;assert(numColumns raw==1422);assert(numColumns basis(7,A)==618);\n'
    program += 'out=openOut "'+str(output)+'";out<<numColumns raw<<endl;\n'
    program += ('scan(numColumns raw,j->{src=select(toList(0..6),i->raw_(i,j)!=0);'
                'assert(#src==1);a=src#0;lf=listForm raw_(a,j);assert(#lf==1 and (lf#0)#1==1);'
                'out<<a;scan((lf#0)#0,d->out<<" "<<d);out<<endl;'
                'nr=sum apply(7,i->#listForm red_(i,j));out<<nr<<endl;'
                'scan(7,i->scan(listForm red_(i,j),term->{out<<i<<" "<<toString(term#1);'
                'scan(term#0,d->out<<" "<<d);out<<endl;}));});close out;\n')
    source = directory/'projection.m2'
    source.write_text(program)
    env = os.environ.copy()
    env['DYLD_LIBRARY_PATH'] = str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
    run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary', '--script', source], env=env, timeout=60)
    lines = iter(output.read_text().splitlines())
    count = int(next(lines))
    assert count == 1422
    entries, keys = {}, set()
    for _ in range(count):
        source = list(map(int, next(lines).split()))
        terms = []
        for _ in range(int(next(lines))):
            row = next(lines).split()
            key = int(row[0]), tuple(map(int, row[2:]))
            terms.append((key, s.Rational(row[1])))
            keys.add(key)
        entries[source[0], tuple(source[1:])] = terms
    keys = sorted(keys)
    assert len(keys) == 618
    indices = {key: i for i, key in enumerate(keys)}
    active, radical = [4, 5, 7, 8], [0, 1, 2, 3, 6]
    base = [tuple(x) for x in evidence['algebra_basis']]
    qm = list(cwr(active, 2))
    quadratics = [dict((tuple(m), s.Rational(c)) for m, c in f)
                  for f in saved['initial_relations'] if len(f[0][0]) == 2]
    Q = s.Matrix([[f.get(m, 0) for m in qm] for f in quadratics])
    reduced, pivots = Q.rref()
    free = [j for j in range(10) if j not in pivots]
    assert [qm[j] for j in free] == [(7, 7), (8, 8)]
    rows = []
    for mon in monomials[7]:
        act = tuple(i for i in mon if i in active)
        exponent = tuple(mon.count(i) for i in radical)
        row = {}
        if len(act) < 3:
            if len(act) < 2:
                terms = [(base.index(act), s.S.One)]
            elif qm.index(act) in free:
                terms = [(5+free.index(qm.index(act)), s.S.One)]
            else:
                terms = [(5+j, -reduced[pivots.index(qm.index(act)), col]) for j, col in enumerate(free)]
            for basis, coefficient in terms:
                if coefficient:
                    for key, c in entries[basis, exponent]:
                        column = indices[key]
                        row[column] = row.get(column, 0)+coefficient*c
        rows.append({j: c for j, c in row.items() if c})
    path = directory/'projection.txt'
    seven_index = {mon: j for j, mon in enumerate(monomials[7])}
    with path.open('w') as out:
        out.write(f'618 {len(monomials[7])}\n')
        for row in rows:
            out.write(str(len(row))+''.join(f' {j} {c}' for j, c in row.items())+'\n')
        for variable in range(9):
            out.write(' '.join(str(seven_index[tuple(sorted(mon+(variable,)))])
                               for mon in monomials[6])+'\n')
    return path, hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    started = time.time()
    primary = ROOT/'reports/S61_FLAG_ZERO12_RATIONAL_INITIAL_20261101.json'
    saved = json.loads(primary.read_text())
    assert saved['all_passed'] and saved['full_truncation_degree'] == 6
    assert all(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
               for path, digest in saved['input_sha256'].items())
    model = S61()
    point = s.Matrix(saved['point61']).applyfunc(s.Rational)
    local = calculate(model, True, point)
    J = s.Matrix(local['jacobian']).applyfunc(s.Rational)
    T = s.Matrix(local['tangent_basis']).applyfunc(s.Rational)
    assert J.rank() == 36 and T.shape == (45, 9)
    one = [i for i, w in enumerate(model.reduced['parameter_weights']) if w == 1]
    _, equations = specialize(model, {i: point[i] for i in one}, True)
    columns = local['chart_columns']
    fixed = next(i for i in range(45) if i not in columns)
    terms = []
    for equation in equations:
        row = []
        for mon, coefficient in equation.items():
            if fixed not in mon:
                a, b = [columns.index(i) for i in mon]+[44]*(2-len(mon))
                row.append((a, b, coefficient))
        terms.append(row)
    pivots = list(J.to_DM().rref()[1])
    selected = list(J[:, pivots].T.to_DM().rref()[1])
    U = s.zeros(45, 36)
    inverse = J[selected, pivots].inv(method='DM')
    for i, column in enumerate(pivots):
        U[column, :] = inverse[i, :]
    assert J*T == s.zeros(135, 9) and J[selected, :]*U == s.eye(36)
    factors = [s.ilcm(*[s.denom(c) for c in list(J[i, :])+[c for _, _, c in terms[i]]])
               for i in range(135)]
    Jn = s.diag(*factors)*J
    Un = U*s.diag(*[1/s.Integer(factors[i]) for i in selected])
    denominator = s.ilcm(*[s.denom(c) for c in Un])
    Un = (denominator*Un).applyfunc(s.cancel)
    tangent_denominator = s.ilcm(*[s.denom(c) for c in T])
    Tn = tangent_denominator*T
    monomials = {d: list(cwr(range(9), d)) for d in range(8)}
    indices = {d: {m: i for i, m in enumerate(monomials[d])} for d in range(8)}
    natives = ['integer_implicit_jets_20261101.cpp', 'truncated_local_ideal_filtered_trace_20261101.cpp',
               'projected_next_initial_replay_20261101.cpp']
    with tempfile.TemporaryDirectory(prefix='zero12-seventh-') as tmp:
        directory = Path(tmp)
        binaries = []
        for name in natives:
            binary = directory/Path(name).stem
            run(['c++', '-O3', '-std=c++17', '-I/opt/homebrew/include', ROOT/'code'/name,
                 '-L/opt/homebrew/lib', '-lgmpxx', '-lgmp', '-o', binary], timeout=60)
            binaries.append(binary)
        projected, projection_hash = projection(saved, directory, monomials)
        data = directory/'integer.txt'
        with data.open('w') as out:
            out.write(f'135 9 7 36\n{denominator} {tangent_denominator}\n')
            out.write(' '.join(str(len(monomials[d])) for d in range(8))+'\n')
            out.write(' '.join(map(str, selected))+'\n')
            for matrix in [Jn, Tn, Un]:
                for row in matrix.tolist():
                    out.write(' '.join(map(str, row))+'\n')
            for i, row in enumerate(terms):
                out.write(str(len(row))+' '+' '.join(f'{a} {b} {c*factors[i]}' for a, b, c in row)+'\n')
            for a in range(1, 7):
                for b in range(1, 8-a):
                    out.write(' '.join(str(indices[a+b][tuple(sorted(x+y))])
                                       for x in monomials[a] for y in monomials[b])+'\n')
        coefficients = directory/'qq7.txt'
        native = run([binaries[0], data, coefficients], timeout=120)
        prime = 32003
        dense = [line.split() for line in coefficients.read_text().splitlines()]
        assert len(dense) == 135 and all(len(row) == 11430 for row in dense)
        modular = []
        for row in dense:
            converted = []
            for value in row[:4995]:
                a, b = map(int, value.split('/'))
                converted.append(a*pow(b, -1, prime) % prime)
            modular.append(converted)
        jets, offset = {}, 0
        for d in range(2, 7):
            jets[d] = [row[offset:offset+len(monomials[d])] for row in modular]
            offset += len(monomials[d])
        modular_input = directory/'modular.txt'
        write_input(modular_input, prime, 9, jets, maximum_degree=6)
        trace, order, initial = [directory/name for name in ['trace.txt', 'order.txt', 'initial.txt']]
        guide = run([binaries[1], modular_input, trace, order, initial], timeout=90)
        assert guide.stdout.strip() == '3858'
        new = directory/'new-seven.txt'
        exact = run([binaries[2], modular_input, coefficients, trace, projected, new], timeout=1800)
        assert not new.read_text().strip()
        assert 'PROJECTED SEVENTH CHECKS 34857 NONZERO 0' in exact.stderr
    paths = [str(Path(__file__).relative_to(ROOT)), str(primary.relative_to(ROOT))]
    paths += ['code/'+name for name in natives]
    paths += ['code/s61_model.py', 'code/s61_minimal_infinity_local_data_20261101.py',
              'code/s61_flag_transverse_obstruction_20261101.py',
              'code/s61_flag_zero12_initial_module_20261101.py', 'code/truncated_local_ideal_io_20261101.py']
    return dict(id='S61_FLAG_ZERO12_SEVENTH_INITIAL_20261101', all_passed=True,
                field='QQ', raw135_reconstructed=True, unchanged_lower_degrees=6,
                top_raw_dimension=6435, top_projected_dimension=618,
                lower_trace_rank=3858, all_closure_checks=34857,
                new_seventh_initial_directions=0, initial_rank_through_seven=9675,
                hilbert_function=[1, 9, 37, 98, 202, 329, 471, 618],
                quotient_by_m8_length=1765, full_initial_ideal_identified=False,
                homogenizer_nilpotence_proved=False, projection_sha256=projection_hash,
                native_log=native.stderr, exact_replay_log=exact.stderr,
                limitations=['Only the specified ZERO12 point and initial degrees through seven.',
                             'No complete initial ideal, completed-ring identification, t-nilpotence, or whole-boundary removal.'],
                proof_note='research/notes/S61_FLAG_ZERO12_SEVENTH_INITIAL_20261101.md',
                input_sha256={path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},
                seconds=time.time()-started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify()
    with args.output.open('x') as out:
        json.dump(result, out, indent=2)
        out.write('\n')
    print('ZERO12 EXACT SEVENTH INITIAL PASS', result['seconds'], flush=True)

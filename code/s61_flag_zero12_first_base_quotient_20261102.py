#!/usr/bin/env python3
"""Exact QQ algebra S/(v0,v1,v3)^2 from the original 135 equations."""
import argparse
import hashlib
import json
import subprocess
import tempfile
import time
from fractions import Fraction
from itertools import combinations_with_replacement as cwr
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_flag_transverse_obstruction_20261101 import specialize

ROOT = Path(__file__).resolve().parents[1]
BASE = (0, 1, 3)
NOTE = 'research/notes/S61_FLAG_ZERO12_FIRST_BASE_QUOTIENT_20261102.md'


def verify():
    started = time.monotonic()
    initial_path = ROOT/'reports/S61_FLAG_ZERO12_RATIONAL_INITIAL_20261101.json'
    bounds_path = ROOT/'reports/S61_FLAG_ZERO12_BASE_ADIC_BOUNDS_20261102.json'
    data_path = ROOT/'code/data/s61_flag_zero12_transverse_monomials_20261102.json'
    initial = json.loads(initial_path.read_text())
    bounds = json.loads(bounds_path.read_text())
    assert initial['all_passed'] and bounds['all_passed']
    assert bounds['graded_module_generating_degree_bound'] == 5
    point = s.Matrix(initial['point61']).applyfunc(s.Rational)
    model = S61()
    one = [i for i,w in enumerate(model.reduced['parameter_weights']) if w == 1]
    local = calculate(model, True, point)
    J = s.Matrix(local['jacobian']).applyfunc(s.Rational)
    T = s.Matrix(local['tangent_basis']).applyfunc(s.Rational)
    assert J.rank() == 36 and T.shape == (45, 9)
    _, equations = specialize(model, dict(zip(one, [point[i] for i in one])), True)
    chart = local['chart_columns']
    fixed = next(i for i in range(45) if i not in chart)
    terms = []
    for f in equations:
        row = []
        for mon,c in f.items():
            if fixed not in mon:
                a,b = [chart.index(i) for i in mon]+[44]*(2-len(mon))
                row.append((a,b,c))
        terms.append(row)
    piv = list(J.to_DM().rref()[1])
    selected = list(J[:,piv].T.to_DM().rref()[1])
    U = s.zeros(45,36)
    inverse = J[selected,piv].inv(method='DM')
    for i,col in enumerate(piv): U[col,:] = inverse[i,:]
    assert J*T == s.zeros(135,9) and J[selected,:]*U == s.eye(36)
    scales = [s.ilcm(*[s.denom(c) for c in list(J[i,:])+[c for a,b,c in terms[i]]]) for i in range(135)]
    Jn = s.diag(*scales)*J
    Un = U*s.diag(*[1/s.Integer(scales[i]) for i in selected])
    denominator = s.ilcm(*[s.denom(c) for c in Un])
    Un = (denominator*Un).applyfunc(s.cancel)
    td = s.ilcm(*[s.denom(c) for c in T]); Tn = td*T
    prime, cutoff = 32003, 6
    assert denominator % prime and td % prime
    ms = {d:[m for m in cwr(range(9),d) if sum(i in BASE for i in m) <= 1] for d in range(7)}
    ix = {d:{m:i for i,m in enumerate(mm)} for d,mm in ms.items()}
    mons = sum((ms[d] for d in range(7)),[])
    index = {m:i for i,m in enumerate(mons)}
    assert len(mons) == 2310
    candidates = [tuple(m) for m in json.loads(data_path.read_text())['monomials']]
    assert len(candidates) == len(set(candidates)) == 49
    assert all(all(i not in BASE for i in m) for m in candidates)
    paths = [ROOT/'code/monomial_quotient_implicit_jets_20261102.cpp',
             ROOT/'code/truncated_local_ideal_filtered_trace_20261101.cpp',
             ROOT/'code/s61_flag_first_base_algebra_engine_20261102.cpp']
    logs = []
    with tempfile.TemporaryDirectory(prefix='zero12-first-base-') as temp:
        wd = Path(temp)
        binaries = []
        for path in paths:
            binary = wd/path.stem
            subprocess.run(['c++','-O3','-w','-std=c++17','-I/opt/homebrew/include',str(path),
                            '-L/opt/homebrew/lib','-lgmpxx','-lgmp','-o',str(binary)],check=True)
            binaries.append(binary)
        native_input = wd/'native.txt'
        with native_input.open('w') as out:
            out.write(f'135 9 6 36\n{denominator} {td}\n')
            out.write(' '.join(str(len(ms[d])) for d in range(7))+'\n')
            out.write(' '.join(map(str,selected))+'\n')
            for M in [Jn,Tn,Un]:
                for row in M.tolist(): out.write(' '.join(map(str,row))+'\n')
            for i,row in enumerate(terms):
                out.write(str(len(row))+' '+' '.join(f'{a} {b} {c*scales[i]}' for a,b,c in row)+'\n')
            for a in range(1,6):
                for b in range(1,7-a):
                    out.write(' '.join(str(ix[a+b].get(tuple(sorted(x+y)),-1)) for x in ms[a] for y in ms[b])+'\n')
        qq = wd/'qq.txt'
        run = subprocess.run([str(binaries[0]),str(native_input),str(qq)],capture_output=True,text=True,check=True,timeout=30)
        logs.append(run.stderr)
        rows = [line.split() for line in qq.read_text().splitlines()]
        assert len(rows) == 136 and all(len(row)==2300 for row in rows)
        t = [Fraction(0)]*10+list(map(Fraction,rows[-1]))
        modular = wd/'modular.txt'
        with modular.open('w') as out:
            out.write(f'{prime} 9 2310 135\n')
            out.write(' '.join(str(len(m)) for m in mons)+'\n')
            for j in range(9): out.write(' '.join(str(index.get(tuple(sorted(m+(j,))),-1)) for m in mons)+'\n')
            for row in rows[:135]:
                sparse = []
                for j,c in enumerate(row,10):
                    a,b = map(int,c.split('/')); v = a*pow(b,-1,prime)%prime
                    if v: sparse.append((j,v))
                out.write(str(len(sparse))+' '+' '.join(f'{j} {c}' for j,c in sparse)+'\n')
            out.write(' '.join('0' for _ in mons)+'\n')
        trace, order, first = [wd/n for n in ['trace.txt','order.txt','initial.txt']]
        run = subprocess.run([str(binaries[1]),str(modular),str(trace),str(order),str(first)],capture_output=True,text=True,check=True,timeout=30)
        assert run.stdout.strip() == '2123'
        pivot_indices = {int(line.split()[0]) for line in order.read_text().splitlines()}
        free = [i for i in range(2310) if i not in pivot_indices]
        assert len(free) == 187
        query_mons = candidates+[tuple(sorted(m+(j,))) for j in BASE for m in candidates]
        queries = [{index[m]:Fraction(1)} for m in query_mons]
        def multiply(poly, linear):
            ans = {}
            for m,c in poly.items():
                for j,a in linear.items():
                    k = index.get(tuple(sorted(m+(j,))))
                    if k is not None: ans[k] = ans.get(k,0)+a*c
            return {k:c for k,c in ans.items() if c}
        special = [39,41,42,43,44,45,46,47,48]
        for j in special: queries.append(multiply({candidates[j]:Fraction(1)}, {3:1,1:-6}))
        tpoly = {m:c for m,c in zip(mons,t) if c}
        for linear in [{3:1,1:-6},{6:1,1:3}]+[{j:1} for j in (2,4,5,7,8)]:
            queries.append(multiply(tpoly,linear))
        query_file = wd/'queries.txt'
        query_file.write_text(str(len(queries))+'\n'+'\n'.join(str(len(q))+' '+' '.join(f'{i} {c}' for i,c in sorted(q.items())) for q in queries)+'\n')
        observable = wd/'observable.txt'; observable.write_text(' '.join(map(str,t))+'\n')
        needed = wd/'needed.txt'; needed.write_text('0\n0\n')
        matrix_file, answers = wd/'tmatrix.txt', wd/'answers.txt'
        run = subprocess.run([str(binaries[2]),str(modular),str(qq),str(trace),str(needed),str(wd/'selected.txt'),
                              str(observable),str(matrix_file),str(query_file),str(answers)],
                             capture_output=True,text=True,check=True,timeout=120)
        assert 'COMMUTATORS PASS 6732' in run.stderr
        assert 'QUOTIENT ALGEBRA PASS 187 original 135 all monomials 2310' in run.stderr
        assert 't_squared_zero' in run.stdout
        logs.append(run.stderr)
        def sparse_rows(file):
            lines = file.read_text().splitlines()
            assert int(lines[0]) == len(lines)-1
            result = []
            for line in lines[1:]:
                row = line.split(); assert len(row) == 2*int(row[0])+1
                result.append({int(i):s.Rational(c) for i,c in zip(row[1::2],row[2::2])})
            return result
        ans = sparse_rows(answers); assert len(ans)==212 and all(not row for row in ans[196:])
        Q = s.zeros(187,196)
        for j,row in enumerate(ans[:196]):
            for i,c in row.items(): Q[i,j]=c
        assert Q.to_DM().rank()==187 and Q[:,49:].to_DM().rank()==138
        tm = sparse_rows(matrix_file); assert len(tm)==187
        Mt = s.zeros(187)
        for j,row in enumerate(tm):
            for i,c in row.items(): Mt[i,j]=c
        fi = {k:i for i,k in enumerate(free)}
        tr = [fi[index[()]]]+[fi[index[(j,)]] for j in BASE]
        assert Mt.to_DM().rank()==3 and Mt[:,tr].to_DM().rank()==3
        assert Mt*Mt==s.zeros(187) and Mt[:,tr[3]]==6*Mt[:,tr[2]]
        exact_hf = [sum(len(mons[k])==d for k in free) for d in range(7)]
        assert exact_hf == [1,9,31,52,55,27,12]
        evidence = dict(quotient_basis=[mons[k] for k in free],
                        transverse_generators=candidates,
                        generator_normal_forms=[[[i,str(c)] for i,c in row.items()] for row in ans[:196]],
                        t_multiplication_columns=[[[i,str(c)] for i,c in row.items()] for row in tm],
                        special_linear_relation_generator_indices=special,
                        native_logs=logs, generated_sha256={n:hashlib.sha256(f.read_bytes()).hexdigest() for n,f in
                            [('raw_integer_input',native_input),('raw_qq_coefficients',qq),('guide_trace',trace),('answers',answers),('t_matrix',matrix_file)]})
    paths += [Path(__file__),initial_path,bounds_path,data_path,ROOT/NOTE,
              ROOT/'code/s61_model.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',
              ROOT/'code/s61_flag_transverse_obstruction_20261101.py']
    return dict(id='S61_FLAG_ZERO12_FIRST_BASE_QUOTIENT_20261102',all_passed=True,field='QQ',
                point61=list(map(str,point)),base_variables=list(BASE),base_power=2,
                exact_quotient_length=187,exact_hilbert_function=exact_hf,
                maximal_ideal_power7_zero=True,transverse_generator_count=49,first_base_layer_dimension=138,
                exact_first_linear_relation_count=9,first_linear_relation_factor='v3-6*v1',
                t_nonzero=True,t_square_zero_in_this_quotient=True,t_multiplication_rank=3,
                t_ideal_as_base_module='R/(a^2,v3-6*v1)',
                closed_squarezero_quotient_kernel_mod_base_squared_length=167,
                full_completed_ring_t_nilpotence_proved=False,relative_ambient_removal_proved=False,
                proof_note=NOTE,limitations=['Only the specified fixed-beta quotient S/a^2 is identified.',
                                            'No full initial ideal, full R-adic presentation or full-ring t-nilpotence.'],
                **evidence,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                seconds=time.monotonic()-started)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); result=verify()
    with args.output.open('x') as out:
        json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))

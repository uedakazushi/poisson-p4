#!/usr/bin/env python3
"""Exact rational initial relations at a sparse flag-boundary point.

This proves containment in the actual initial ideal. It does not identify the
complete initial ideal, the completed ring, or nilpotence of the homogenizer.
"""
import argparse
import hashlib
import json
import os
import subprocess
import tempfile
import time
from collections import Counter
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS, sp_basis, infinitesimal_action
from truncated_local_ideal_io_20261101 import write_input
ROOT = Path(__file__).resolve().parents[1]


def geometry(model, case, H, b, T):
    z, bu, beta = primitive_beta(model)
    v = s.Matrix([beta[i, j].diff(z[k]) for (i, j), k in COORDS])
    B = v.jacobian(bu)
    rr = list(B.T.to_DM().rref()[1])
    left = B[rr, :].inv()
    act = s.Matrix.hstack(*[(left*infinitesimal_action(v, A)[rr, :]).applyfunc(s.expand)
                            for A in sp_basis()]).subs(dict(zip(bu, H*b)), simultaneous=True)
    QD = s.Matrix(case['q_orbit_differential']).applyfunc(s.Rational)
    BD = H.row_join(act)
    joint = BD.col_join(s.zeros(44, 15).row_join(QD))
    K = s.Matrix.hstack(*BD.nullspace())
    C = QD*K[15:, :]
    cp = list(C.to_DM().rref()[1])
    C = C[:, cp].col_join(s.zeros(1, len(cp)))
    fr = list(T.T.to_DM().rref()[1])
    CT = T[fr, :].inv()*C[fr, :]
    assert T*CT == C
    pr = list(CT.T.to_DM().rref()[1])
    lv = K[[1], cp]*CT[pr, :].inv()
    L = s.zeros(1, 9)
    for j, col in enumerate(pr):
        L[0, col] = lv[j]
    assert L*CT == K[[1], cp]
    HC = H[:, [j for j in range(15) if j != 1]]
    boundary = HC.row_join(act)
    bj = boundary.col_join(s.zeros(44, 14).row_join(QD))
    assert (BD.rank(), joint.rank(), boundary.rank(), bj.rank(), CT.cols) == (16, 19, 16, 18, 3)
    assert pr == [0, 1, 2] and list(L) == [0, 0, -s.Rational(1, 2), 0, 0, 0, 0, 0, 0]
    enc = lambda M: [[str(v) for v in row] for row in M.tolist()]
    return dict(beta_rank=16, joint_rank=19, boundary_beta_rank=16, boundary_joint_rank=18,
                actual_flag_tangent=enc(CT), parameters=pr, b1_derivative=enc(L))


def verify():
    started = time.time()
    source = ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json'
    saved = json.loads(source.read_text())
    assert saved['all_passed']
    case = saved['flag']
    model = S61()
    one = [i for i, w in enumerate(model.reduced['parameter_weights']) if w == 1]
    H = s.Matrix(case['beta_kernel']).applyfunc(s.Rational)
    b = s.zeros(15, 1)
    b[0] = b[12] = 1
    point = s.Matrix(case['point61']).applyfunc(s.Rational)
    for i, value in zip(one, H*b):
        point[i] = value
    local = calculate(model, True, point)
    assert (local['leading_jacobian_rank'], local['homogenized_jacobian_rank'], local['quadratic_span_rank']) == (35, 36, 9)
    J = s.Matrix(local['jacobian']).applyfunc(s.Rational)
    T = s.Matrix(local['tangent_basis']).applyfunc(s.Rational)
    geo = geometry(model, case, H, b, T)
    two, F = specialize(model, dict(zip(one, H*b)), True)
    cols = local['chart_columns']
    fixed = next(i for i in range(45) if i not in cols)
    terms = []
    for f in F:
        row = []
        for mon, c in f.items():
            if fixed not in mon:
                a, bb = [cols.index(i) for i in mon]+[44]*(2-len(mon))
                row.append((a, bb, c))
        terms.append(row)
    C = list(J.to_DM().rref()[1])
    rows = list(J[:, C].T.to_DM().rref()[1])
    U = s.zeros(45, len(C))
    iv = J[rows, C].inv(method='DM')
    for i, col in enumerate(C):
        U[col, :] = iv[i, :]
    assert J*T == s.zeros(135, 9) and J[rows, :]*U == s.eye(36)
    factors = [s.ilcm(*[s.denom(v) for v in list(J[i, :])+[c for _, _, c in terms[i]]]) for i in range(135)]
    Jn = s.diag(*factors)*J
    Un = U*s.diag(*[1/s.Integer(factors[i]) for i in rows])
    den = s.ilcm(*[s.denom(v) for v in Un])
    Un = (den*Un).applyfunc(s.cancel)
    td = s.ilcm(*[s.denom(v) for v in T])
    Tn = td*T
    prime = 32003
    assert den % prime and td % prime
    ms = {d: list(combinations_with_replacement(range(9), d)) for d in range(7)}
    ix = {d: {m: i for i, m in enumerate(mm)} for d, mm in ms.items()}
    mons = sum((ms[d] for d in range(7)), [])
    stems = ['integer_implicit_jets', 'truncated_local_ideal_filtered_trace',
             'truncated_local_ideal_minimal_generators', 'truncated_local_ideal_adaptive_replay']
    paths = [ROOT/'code'/f'{stem}_20261101.cpp' for stem in stems]
    with tempfile.TemporaryDirectory(prefix='poisson-flag-rational-') as temp:
        temp = Path(temp)
        def run(args, timeout=90):
            r = subprocess.run(list(map(str, args)), cwd=ROOT, capture_output=True, text=True, timeout=timeout)
            assert r.returncode == 0, r.stdout+r.stderr
            return r
        executables = []
        for path in paths:
            binary = temp/path.stem
            run(['c++', '-O3', '-std=c++17', '-I/opt/homebrew/include', path,
                 '-L/opt/homebrew/lib', '-lgmpxx', '-lgmp', '-o', binary], 60)
            executables.append(binary)
        native_input = temp/'integer.txt'
        with native_input.open('w') as out:
            out.write(f'135 9 6 36\n{den} {td}\n')
            out.write(' '.join(str(len(ms[d])) for d in range(7))+'\n')
            out.write(' '.join(map(str, rows))+'\n')
            for M in [Jn, Tn, Un]:
                for row in M.tolist():
                    out.write(' '.join(map(str, row))+'\n')
            for i, ht in enumerate(terms):
                out.write(str(len(ht))+' '+' '.join(f'{a} {bb} {c*factors[i]}' for a, bb, c in ht)+'\n')
            for a in range(1, 6):
                for bb in range(1, 7-a):
                    out.write(' '.join(str(ix[a+bb][tuple(sorted(x+y))]) for x in ms[a] for y in ms[bb])+'\n')
        coefficients = temp/'qq.txt'
        native = run([executables[0], native_input, coefficients])
        assert 'EXACT INTEGER JETS PASS' in native.stderr
        dense = [line.split() for line in coefficients.read_text().splitlines()]
        assert len(dense) == 135 and all(len(row) == 4995 for row in dense)
        modular = []
        for row in dense:
            rr = []
            for value in row:
                a, bb = map(int, value.split('/'))
                rr.append(a*pow(bb, -1, prime) % prime)
            modular.append(rr)
        jets = {}
        offset = 0
        for d in range(2, 7):
            jets[d] = [row[offset:offset+len(ms[d])] for row in modular]
            offset += len(ms[d])
        inp = temp/'modular.txt'
        write_input(inp, prime, 9, jets, maximum_degree=6)
        trace, order, initial = [temp/name for name in ['trace.txt', 'order.txt', 'initial.txt']]
        closure = run([executables[1], inp, trace, order, initial])
        assert closure.stdout.strip() == '3863'
        selected = temp/'selected.txt'
        selection = run([executables[2], inp, initial, selected])
        chosen = list(map(int, selected.read_text().split()))
        assert len(chosen) == 61
        order_rows = [list(map(int, line.split())) for line in order.read_text().splitlines()]
        target = [order_rows[i][1] for i in chosen]
        # The adaptive replay computes dependencies over QQ, including reductions
        # whose coefficients happened to vanish modulo the guide prime.
        needed = temp/'needed.txt'
        needed.write_text('0\n'+str(len(target))+' '+' '.join(map(str, target))+'\n')
        exact = temp/'initial-qq.txt'
        replay = run([executables[3], inp, coefficients, trace, needed, exact, "verify-closure"], 480)
        assert 'EXACT ADAPTIVE REPLAY PASS 61' in replay.stderr
        assert 'EXACT CLOSURE PASS 12105 rank 3863' in replay.stderr
        exact_rows = exact.read_text().splitlines()
        assert int(exact_rows[0]) == len(exact_rows)-1 == 61
        expr, relations, degrees = [], [], []
        for line in exact_rows[1:]:
            row = line.split()
            assert len(row) == 1+2*int(row[0])
            monco = [(mons[int(row[k])], row[k+1]) for k in range(1, len(row), 2)]
            dd = {len(mon) for mon, c in monco}
            assert len(dd) == 1
            degrees.append(dd.pop())
            expr.append('+'.join('('+c+')*'+'*'.join('v_'+str(j) for j in mon) for mon, c in monco))
            relations.append([[list(mon), c] for mon, c in monco])
        assert dict(Counter(degrees)) == {2: 9, 3: 3, 4: 12, 5: 22, 6: 15}
        program = 'R=QQ[v_0..v_8];J=ideal('+','.join(expr)+');\n'
        program += '''C=ideal(v_8,v_7,v_5,v_4,v_3+2*v_6,v_2+4*v_6);
assert(dim(R/J)==3 and degree(R/J)==5 and radical J==C);
B=saturate(J,ideal v_2);E=B/J;A=ann E;
assert(dim(R/B)==3 and degree(R/B)==5 and dim E==2 and degree E==111);
assert(radical A==ideal(v_2,v_3,v_4,v_5,v_6,v_7,v_8));
assert(degree(J+ideal(v_0,v_1,v_2))==50);
hB=hilbertSeries(R/B,Reduce=>true);hJ=hilbertSeries(R/J,Reduce=>true);zz=(ring numerator hB)_0;
assert(numerator hB==1+4*zz and value denominator hB==(1-zz)^3);
assert(numerator hJ==1+6*zz+12*zz^2+15*zz^3+9*zz^4-20*zz^5-8*zz^6-10*zz^7 and value denominator hJ==(1-zz)^3);
print("RATIONAL INITIAL STRUCTURE PASS");stdio<<flush;
'''
        normal_terms = [(j, mon) for j, mon in enumerate(mons[10:]) if all(k > 2 for k in mon)]
        transverse = ['+'.join('('+row[j]+')*'+'*'.join('w_'+str(k) for k in mon)
                              for j, mon in normal_terms if row[j] != '0/1') or '0' for row in dense]
        program += 'S=QQ[w_3..w_8];m=ideal vars S;I=ideal('+','.join(transverse)+');\n'
        program += 'assert(degree(I+m^6)==49 and degree(I+m^7)==49);print("TRANSVERSE LENGTH49 PASS");\n'
        m2file = temp/'verify.m2'
        m2file.write_text(program)
        env = os.environ.copy()
        env['DYLD_LIBRARY_PATH'] = str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
        m2 = subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary', '--script', str(m2file)],
                            cwd=ROOT, env=env, capture_output=True, text=True, timeout=120)
        assert m2.returncode == 0 and 'TRANSVERSE LENGTH49 PASS' in m2.stdout, m2.stdout+m2.stderr
        evidence = dict(native_log=native.stderr, selection_log=selection.stderr, replay_log=replay.stderr,
                        m2_log=m2.stdout+m2.stderr, exact_initial_sha256=hashlib.sha256(exact.read_bytes()).hexdigest(),
                        initial_relations=relations)
    paths += [Path(__file__), source, ROOT/'code/s61_model.py', ROOT/'code/s61_minimal_infinity_local_data_20261101.py',
              ROOT/'code/s61_flag_transverse_obstruction_20261101.py', ROOT/'code/truncated_local_ideal_io_20261101.py',
              ROOT/'code/primitive_jet_elliptic_curve_20261101.py', ROOT/'code/weight1_structure.py',
              ROOT/'research/notes/S61_FLAG_BOUNDARY_RATIONAL_INITIAL_20261101.md']
    return dict(id='S61_FLAG_BOUNDARY_RATIONAL_INITIAL_20261101', all_passed=True, field='QQ',
                point61=list(map(str, point)), H15_coefficients=list(map(str, b)), geometry=geo,
                initial_relation_degrees=dict(Counter(degrees)), initial_subideal_dimension=3,
                initial_subideal_degree=5, full_truncation_degree=6, full_truncation_length=1142, embedded_quotient_dimension=2, embedded_quotient_degree=111,
                transverse_length=49, full_initial_ideal_identified=False, homogenizer_nilpotence_proved=False,
                limitations=['J is contained in the actual initial ideal; equality is not proved.',
                             'Transverse length49 is not the rank of the full ring; the known generic flag rank is5.',
                             'No pointwise or ambient t-nilpotence, and no entire flag-boundary removal.'],
                **evidence, input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                seconds=time.time()-started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = verify()
    with args.output.open('x') as out:
        json.dump(result, out, indent=2)
        out.write('\n')
    print(json.dumps(dict(all_passed=True, seconds=result['seconds'])))

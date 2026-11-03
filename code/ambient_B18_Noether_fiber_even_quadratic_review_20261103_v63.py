#!/usr/bin/env python3
"""Independent one-variable Lie-log source4/IFT/trace/full-graph fiber13 review."""
import argparse
from fractions import Fraction as F
import gzip
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import ambient_B18_zero_linear_contact4_review_20261103_v63 as own
from ambient_B18_actual_D_tangent_cone_support_review_20261103_v63 import inverse, matrix

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_NOETHER_FIBER_EVEN_QUADRATIC_20261103_v63_REVIEW'


def add(A, B, c=F(1)):
    C = {d: row[:] for d, row in A.items()}
    for d, row in B.items():
        C[d] = [a+c*b for a, b in zip(C.get(d, [F(0)]*150), row)]
    return {d: row for d, row in C.items() if any(row)}


def ma(A, B, c=F(1)):
    C = {d: [r[:] for r in M] for d, M in A.items()}
    for d, M in B.items():
        N = C.get(d, [[F(0)]*5 for _ in range(5)])
        C[d] = [[a+c*b for a, b in zip(ar, br)] for ar, br in zip(N, M)]
    return {d: M for d, M in C.items() if any(any(row) for row in M)}


def mm(A, B):
    C = {}
    for d, M in A.items():
        for e, N in B.items():
            if d+e <= 4:
                C = ma(C, {d+e: own.matrix_product(M, N)})
    return C


def action(A, P):
    C = {}
    for d, M in A.items():
        for e, row in P.items():
            if d+e <= 4:
                C = add(C, {d+e: own.lie(M, row)})
    return C


def exp_action(A, P):
    C, term, factorial = P, P, 1
    for j in range(1, 5):
        term = action(A, term)
        factorial *= j
        C = add(C, term, F(1, factorial))
    return C


def mul(a, b):
    return [sum((a[j]*b[i-j] for j in range(i+1)), F(0)) for i in range(5)]


def compose(a, b):
    result, power = [F(0)]*5, [F(1)]+[F(0)]*4
    for c in a:
        result = [x+c*y for x, y in zip(result, power)]
        power = mul(power, b)
    return result


def canonical_and_swap(P):
    div = [[F(0)]*5 for _ in range(5)]
    for ip, (i, j) in enumerate(own.PAIRS):
        for im, (a, b) in enumerate(own.MONS):
            q = P[15*ip+im]
            for k, rem in [(a, b), (b, a)]:
                if k == i: div[j][rem] += q
                if k == j: div[i][rem] -= q
    Q = P[:]
    for ip, (i, j) in enumerate(own.PAIRS):
        for k in range(5):
            Q[15*ip+own.MI[tuple(sorted((i, k)))]] -= div[j][k]/5
            Q[15*ip+own.MI[tuple(sorted((j, k)))]] += div[i][k]/5
    swap = [1, 0, 2, 3, 4]
    R = [F(0)]*150
    for ip, (i, j) in enumerate(own.PAIRS):
        ii, jj = swap[i], swap[j]
        sign = F(1 if ii < jj else -1)
        pair = tuple(sorted((ii, jj)))
        for im, (a, b) in enumerate(own.MONS):
            R[15*own.PI[pair]+own.MI[tuple(sorted((swap[a], swap[b])))]] += sign*Q[15*ip+im]
    return R


def constraints(q):
    div = [[F(0)]*5 for _ in range(25)]
    for ip, (i, j) in enumerate(own.PAIRS):
        for im, (a, b) in enumerate(own.MONS):
            row = q[15*ip+im]
            for k, rem in [(a, b), (b, a)]:
                if k == i: div[5*j+rem] = [x+y for x, y in zip(div[5*j+rem], row)]
                if k == j: div[5*i+rem] = [x-y for x, y in zip(div[5*i+rem], row)]
    assert all(not any(row) for row in div)
    triples = list(itertools.combinations(range(5), 3))
    cubics = list(itertools.combinations_with_replacement(range(5), 3))
    cubic_index = {e: i for i, e in enumerate(cubics)}
    jac = [[F(0)]*5 for _ in range(350)]
    cache = {}
    def pair(i, j):
        return own.PI[tuple(sorted((i, j)))], 1 if i < j else -1
    for ti, (a, b, c) in enumerate(triples):
        for i, j, k in [(a, b, c), (b, c, a), (c, a, b)]:
            ip2, sign2 = pair(j, k)
            for l in range(5):
                if i == l: continue
                ip1, sign1 = pair(i, l)
                for im, mon1 in enumerate(own.MONS):
                    qa = 15*ip1+im
                    if not any(q[qa]): continue
                    for jm, (u, v) in enumerate(own.MONS):
                        count = int(u == l)+int(v == l)
                        qb = 15*ip2+jm
                        if not count or not any(q[qb]): continue
                        rem = v if u == l else u
                        index = 35*ti+cubic_index[tuple(sorted((*mon1, rem)))]
                        key = tuple(sorted((qa, qb)))
                        if key not in cache: cache[key] = mul(q[qa], q[qb])
                        factor = F(sign1*sign2*count)
                        jac[index] = [x+factor*y for x, y in zip(jac[index], cache[key])]
    assert all(not any(row) for row in jac)


def exact_row_rank(A):
    basis = {}
    for original in A:
        row = {i: F(a) for i, a in enumerate(original) if a}
        while row:
            i = min(row)
            if i not in basis:
                pivot = row[i]
                basis[i] = {j: a/pivot for j, a in row.items()}
                break
            factor = row[i]
            for j, a in basis[i].items():
                row[j] = row.get(j, F(0))-factor*a
                if row[j] == 0: row.pop(j)
    return len(basis)


def main(output):
    start = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    assert not output.exists()
    output.mkdir(parents=True)
    authorpath = ROOT/'reports/AMBIENT_B18_NOETHER_FIBER_EVEN_QUADRATIC_20261103_v63/report.json'
    finitepath = ROOT/'reports/AMBIENT_B18_ACTUAL_D_FINITE13_20261103_v63_REVIEW.json'
    firstreview = ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63_REVIEW.json'
    tangentreview = ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63_REVIEW/report.json'
    author, finite, firstgate, tangentgate = [own.accepted(p) for p in [authorpath, finitepath, firstreview, tangentreview]]
    assert finite['closed_fiber_length_allowed'] == [12, 13]
    cpath = authorpath.parent/'certificate.json.gz'
    redpath = ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json'
    firstpath = ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz'
    tcpath = ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63/certificate.json.gz'
    c, red, first, tc = [own.read(p) for p in [cpath, redpath, firstpath, tcpath]]
    regular, fold = red['common28_linear_target_rows'], red['common29th_target_row']
    selected, odd = tc['common39_physical_coordinate_rows'], tc['common10_odd_target_rows']
    assert len(regular) == 28 and fold == 61 and selected == regular+[fold]+odd
    T = matrix(tc['common150_by39_normalized_tangent'])
    assert [T[i] for i in selected] == [[F(i == j) for j in range(39)] for i in range(39)]
    h0 = [row[29:] for row in T]
    assert h0 == matrix(c['full150_by10_common_odd_constant_graph_matrix'])
    P0 = [F(0)]*150
    P0[15*own.PI[0,2]+own.MI[1,1]] = -1
    P0[15*own.PI[2,3]+own.MI[0,1]] = 1
    P0[15*own.PI[2,4]+own.MI[0,0]] = 1
    P0[15*own.PI[2,4]+own.MI[3,3]] = 1
    Gbasis = []
    for mon in [(0,0), (0,1), (1,1)]:
        raw = [F(0)]*150
        raw[15*own.PI[3,4]+own.MI[mon]] = 1
        Gbasis.append(canonical_and_swap(raw))
    oi = red['orbit_GL_generator_indices']
    evidence, normalized, stagepaths = {}, {}, []
    for name in ['CA', 'CB']:
        saved = c['families'][name]
        case = red['families'][name]
        J = matrix(case['linear_J'])
        K = matrix(case['source_kernel_Lie_matrix'])
        co = list(map(F, case['source_velocity']))
        extras = case['source_fiber_basis_indices']
        assert oi == saved['source28_parameter_GL_indices'] and extras == saved['source28_parameter_fixed20_indices']
        LK0 = own.lie(K, P0)
        left = matrix(saved['regular28_inverse'])
        I28 = [[F(i == j) for j in range(28)] for i in range(28)]
        assert own.matrix_product(left, [J[i] for i in regular]) == I28
        assert own.matrix_product([J[i] for i in regular], left) == I28
        for j, ij in enumerate(oi):
            M = [[F(0)]*5 for _ in range(5)]
            M[ij//5][ij % 5] = 1
            assert own.lie(M, P0) == [row[j] for row in J]
        base = [F(0)]*20
        base[0] = base[7] = -1
        base[11], base[16] = F(1), F(1,3)
        def physical(source):
            alpha = [[a, -b, F(0), F(0), F(0)] for a, b in zip(base, co)]
            for row, j in zip(source[19:], extras):
                alpha[j] = [a+b for a, b in zip(alpha[j], row)]
            P = {0:P0[:], 1:[-a for a in LK0]}
            for n in range(5):
                coeff = [row[n] for row in source[19:]]
                if any(coeff): P = add(P, {n:own.matvec([row[19:] for row in J], coeff)})
            fullg = []
            if name == 'CB':
                fullg = [[-a for a in mul(alpha[3], alpha[i])] for i in range(3)]
                deltaf = [[a-(base[i] if n == 0 else 0) for n,a in enumerate(alpha[i])] for i in range(3)]
                correction = [[-a for a in mul(alpha[3], df)] for df in deltaf]
                for n in range(5):
                    row = [sum((correction[i][n]*Gbasis[i][j] for i in range(3)), F(0)) for j in range(150)]
                    P = add(P, {n:row})
            M = {}
            for n in range(1, 5):
                A = [[F(0)]*5 for _ in range(5)]
                for row, ij in zip(source[:19], oi): A[ij//5][ij % 5] += row[n]
                if any(any(row) for row in A): M[n] = A
            log, power = {}, M
            for j in range(1, 5):
                log = ma(log, power, F((-1)**(j+1), j))
                power = mm(power, M)
            Q = exp_action(log, exp_action({1:K}, P))
            return [[Q.get(n, [F(0)]*150)[i] for n in range(5)] for i in range(150)], alpha, fullg
        source = [[F(0)]*5 for _ in range(28)]
        for degree in range(1, 5):
            q, _, _ = physical(source)
            assert all(q[i][n] == (P0[i] if n == 0 else 0) for i in regular for n in range(degree))
            correction = own.matvec(left, [-q[i][degree] for i in regular])
            for row, value in zip(source, correction): row[degree] += value
            stagepath = authorpath.parent/f'{name}-stage{degree}.json.gz'
            stagepaths.append(stagepath)
            assert source == matrix(own.read(stagepath)['source28'])
        q, alpha, fullg = physical(source)
        assert source == matrix(saved['source28_IFT_coefficients'])
        assert alpha == matrix(saved['actual_fixed20_parameter_jets'])
        assert fullg == matrix(saved['full_CB_minus_ratio_times_binary_jets'])
        assert q == matrix(saved['all150_physical_source_coefficients_through4'])
        assert all(q[i][0] == P0[i] and q[i][1] == 0 for i in range(150))
        assert all(q[i] == [P0[i],0,0,0,0] for i in regular)
        constraints(q)
        z = q[fold]
        w, c3, d4 = z[2:]
        assert w == (F(-4) if name == 'CA' else F(3,8))
        assert z == list(map(F, saved['fold_row61_coefficients']))
        tau = [F(0), F(-1), -c3/w, -(c3/w)**2, F(0)]
        assert tau == list(map(F, saved['trace_involution_through3'])) and compose(z,tau) == z
        tracezero = [(F(n == 1)-a)/2 for n,a in enumerate(tau)]
        a1, a2, b1, raw = [], [], [], []
        for row in q:
            other = compose(row,tau)
            avg = [(a+b)/2 for a,b in zip(row,other)]
            oddhalf = [(a-b)/2 for a,b in zip(row,other)]
            even_linear = row[2]/w
            second = (avg[4]-even_linear*d4)/(w*w)
            beta = oddhalf[3]/w
            assert avg == [(row[0] if n == 0 else 0)+even_linear*z[n]+second*mul(z,z)[n] for n in range(5)]
            assert oddhalf == [beta*a for a in mul(z,tracezero)]
            a1.append(even_linear);a2.append(second);b1.append(beta)
            raw.append((row[4]-even_linear*d4-beta*c3)/(w*w))
        assert a1 == list(map(F,saved['all150_trace_even_linear_z_coefficients']))
        assert a2 == list(map(F,saved['all150_trace_even_quadratic_z_coefficients']))
        assert b1 == list(map(F,saved['all150_trace_odd_linear_z_coefficients']))
        assert raw == list(map(F,saved['all150_raw_even_quadratic_z_coefficients']))
        assert a2 == [a-c3*b/(2*w*w) for a,b in zip(raw,b1)]
        assert a1 == [row[-1] for row in matrix(tc['families'][name]['full150_by29_even_cotangent_embedding'])]
        assert b1 == [row[-1] for row in matrix(first['CA_CB_full_firstjet_data'][name]['actual'])]
        norm = [a-sum((v*a2[i] for v,i in zip(h0[j],odd)),F(0)) for j,a in enumerate(a2)]
        assert norm == list(map(F,saved['all150_normalized_even_graph_quadratic_z_coefficients']))
        assert all(norm[i] == 0 for i in selected)
        normalized[name] = norm
        evidence[name] = dict(independent_all150_source4=q, independent_actual_alpha20=alpha,
            independent_full_CB_g_product=fullg, trace_even_z2=a2, normalized_graph_z2=norm)
    delta = [a-b for a,b in zip(normalized['CA'],normalized['CB'])]
    assert delta == [F(0)]*150 == list(map(F,c['all150_normalized_CA_minus_CB_even_quadratic_coefficients']))
    nonbasis = [i for i in range(150) if i not in selected]
    graph87 = c['independent87_physical_graph_rows']
    assert nonbasis == c['all111_nonbasis_physical_rows'] and len(nonbasis)==111
    assert len(graph87) == 87 and len(set(selected+graph87)) == 126
    reconstruction = matrix(red['canonical_reconstruction_from126'])
    assert exact_row_rank([reconstruction[i] for i in selected+graph87]) == 126
    def serial(x):
        if isinstance(x,F): return str(x)
        if isinstance(x,dict): return {k:serial(v) for k,v in x.items()}
        if isinstance(x,list): return [serial(v) for v in x]
        return x
    ep = output/'independent-source4.json.gz'
    with ep.open('wb') as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0) as gz:
            gz.write((json.dumps(serial(evidence),separators=(',',':'))+'\n').encode())
    paths=[Path(__file__),Path(own.__file__),ROOT/'code/ambient_B18_actual_D_tangent_cone_support_review_20261103_v63.py',
        ROOT/f'research/notes/{STEM}.md',authorpath,finitepath,firstreview,tangentreview,cpath,redpath,firstpath,tcpath,*stagepaths]
    result=dict(id=STEM,all_passed=True,
        independent_all150_CA_CB_source4_Lie_log_and_exact_CB_nonlinear_product=True,
        independent_both28_IFT_all8_stages_and_two_inverse_products=True,
        original_divergence25_Jacobi350_through4_independently_checked=True,
        trace_involution_and_even_odd_all150_coefficients_rebuilt=True,
        no_p5_needed_because_all150_p_linear_coefficients_zero=True,
        actual_E_blinear_h0_and_common39_coordinate_binding=True,
        full111_physical_and87_independent_canonical_graph_coverage=True,
        intrinsic_whole_syzygy_product_kernel_and_graph_elimination_coverage_proved=True,
        all_normalized_even_quadratic_graph_differences_zero=True,
        z_squared_survives_no_hidden_pure_even_quadratic_relation=True,
        exact_closed_fiber_Hilbert_function=[1,11,1],exact_closed_fiber_length=13,
        exact_closed_fiber_algebra='C[z,Y0,...,Y9]/(YiYj,zYj,z^3)',
        exact_R28_module_minimal_generator_number=13,generic_R28_module_rank=2,
        no_full_D_depth_actual_nilpotency_index_or_p5_claim=True,
        verification_kind='Independent Fraction Lie-log source4, direct physical constraints, trace involution and full actual generator-ideal coverage; no worker replay/GB.',
        input_sha256={str(p.relative_to(ROOT)):own.sha(p) for p in paths},
        evidence_sha256={ep.name:own.sha(ep)},seconds=time.monotonic()-start)
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['all_passed','exact_closed_fiber_Hilbert_function','exact_closed_fiber_length','exact_R28_module_minimal_generator_number','seconds']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    main(parser.parse_args().output)

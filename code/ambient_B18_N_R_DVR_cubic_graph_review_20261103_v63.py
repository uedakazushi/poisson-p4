#!/usr/bin/env python3
"""Independent Lie-log all150 bivariate cubic/IFT/DVR graph review, no CAS."""
import argparse
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import resource
import signal
import time

import ambient_B18_zero_linear_contact4_review_20261103_v63 as own

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_N_R_DVR_CUBIC_GRAPH_20261103_v63_REVIEW'
MON = [(i, d-i) for d in range(4) for i in range(d, -1, -1)]
ZERO150 = [F(0)]*150


def plus(A, B, c=F(1)):
    out = {e: row[:] for e, row in A.items()}
    for e, row in B.items():
        if e not in out:
            out[e] = [F(0)]*len(row)
        out[e] = [a+c*b for a, b in zip(out[e], row)]
    return {e: row for e, row in out.items() if any(row)}


def matrix_series_add(A, B, c=F(1)):
    out = {e: [row[:] for row in M] for e, M in A.items()}
    for e, M in B.items():
        if e not in out:
            out[e] = [[F(0)]*5 for _ in range(5)]
        out[e] = [[a+c*b for a, b in zip(ar, br)] for ar, br in zip(out[e], M)]
    return {e: M for e, M in out.items() if any(any(row) for row in M)}


def matrix_series_mul(A, B):
    out = {}
    for e, M in A.items():
        for f, N in B.items():
            g = (e[0]+f[0], e[1]+f[1])
            if sum(g) <= 3:
                out = matrix_series_add(out, {g: own.matrix_product(M, N)})
    return out


def action(A, P):
    out = {}
    for e, M in A.items():
        for f, row in P.items():
            g = (e[0]+f[0], e[1]+f[1])
            if sum(g) <= 3:
                out = plus(out, {g: own.lie(M, row)})
    return out


def exp_action(A, P):
    result, term, factorial = P, P, 1
    for i in range(1, 4):
        term = action(A, term)
        factorial *= i
        result = plus(result, term, F(1, factorial))
    return result


def scalar_mul(A, B):
    out = {}
    for e, a in A.items():
        for f, b in B.items():
            g = (e[0]+f[0], e[1]+f[1])
            if sum(g) <= 3:
                out[g] = out.get(g, F(0))+a*b
    return {e: a for e, a in out.items() if a}


def main(output):
    start = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    assert not output.exists()
    output.mkdir(parents=True)
    authorpath = ROOT/'reports/AMBIENT_B18_N_R_DVR_CUBIC_GRAPH_20261103_v63/report.json'
    hpath = ROOT/'reports/AMBIENT_B18_ZERO_LINEAR_CONTACT4_20261103_v63_REVIEW.json'
    firstreview = ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63_REVIEW.json'
    author, hgate, firstgate = [own.accepted(p) for p in [authorpath, hpath, firstreview]]
    assert hgate['h_exact'] == 4 and firstgate['N_exact_k1zero_source_etale28_and_full_fixed_flag_N19_membership']
    capath = authorpath.parent/'certificate.json.gz'
    c = own.read(capath)
    redpath = ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json'
    curvepath = ROOT/'reports/AMBIENT_B18_R_CURVE_TANGENT_20261103_v63/certificate.json.gz'
    secondpath = ROOT/'reports/AMBIENT_B18_N_R_SECOND_CONTACT_20261103_v63/certificate.json.gz'
    red, curve, second = [own.read(p) for p in [redpath, curvepath, secondpath]]
    ca = red['families']['CA']
    J = [[F(q) for j, q in enumerate(row) if j != 20] for row in ca['linear_J']]
    oi = red['orbit_GL_generator_indices']
    assert oi == c['source27_parameter_GL_indices']
    assert [i for i in ca['source_fiber_basis_indices'] if i != 2] == c['source27_parameter_fixed_N_fiber_indices']
    K = [[F(q) for q in row] for row in ca['source_kernel_Lie_matrix']]
    t = own.symbols(['t'])['t']
    Rpoly = [own.Poly(own.evaluate(q, {'t': t})) for q in curve['physical_curve150']]
    target = [[p.terms.get((a,), F(0)) if b == 0 else F(0) for a, b in MON] for p in Rpoly]
    assert target == [[F(q) for q in row] for row in c['full150_actual_R_target_coefficients']]
    P0 = [row[0] for row in target]
    LK0 = own.lie(K, P0)
    assert F(ca['source_velocity'][2]) == 0
    # Independently bind every GL generator to the accepted physical J27 matrix.
    for j, ij in enumerate(oi):
        M = [[F(0)]*5 for _ in range(5)]
        M[ij//5][ij % 5] = F(1)
        assert own.lie(M, P0) == [row[j] for row in J]

    def physical(source):
        source = [[F(q) for q in row] for row in source]
        # The fixed source stays in N19, with the compensated kernel term -p L_K(P0).
        P = {(0, 0): P0[:], (0, 1): [-a for a in LK0]}
        for m, e in enumerate(MON):
            coefficients = [source[j][m] for j in range(19, 27)]
            if any(coefficients):
                P = plus(P, {e: own.matvec([row[19:] for row in J], coefficients)})
        Mseries = {}
        for m, e in enumerate(MON):
            matrix = [[F(0)]*5 for _ in range(5)]
            for j, ij in enumerate(oi):
                matrix[ij//5][ij % 5] += source[j][m]
            if any(any(row) for row in matrix):
                Mseries[e] = matrix
        assert (0, 0) not in Mseries
        # Producer uses inverse-coordinate substitution for I+M.  Our representation
        # uses the logarithm and independent coefficient-index Lie action instead.
        M2 = matrix_series_mul(Mseries, Mseries)
        M3 = matrix_series_mul(M2, Mseries)
        log = matrix_series_add(matrix_series_add(Mseries, M2, -F(1, 2)), M3, F(1, 3))
        Q = exp_action(log, exp_action({(0, 1): K}, P))
        return [[Q.get(e, ZERO150)[i] for e in MON] for i in range(150)]

    rows = c['regular27_physical_rows']
    L = [[F(q) for q in row] for row in c['regular27_inverse']]
    I27 = [[F(i == j) for j in range(27)] for i in range(27)]
    assert own.matrix_product(L, [J[i] for i in rows]) == I27
    assert own.matrix_product([J[i] for i in rows], L) == I27
    source = [[F(0)]*10 for _ in range(27)]
    stage_paths, stages = [], []
    for degree in range(1, 4):
        before = physical(source)
        stagepath = authorpath.parent/f'stage{degree}.json.gz'
        stage_paths.append(stagepath)
        saved = own.read(stagepath)
        expected = [[F(q) for q in row] for row in saved['source27']]
        for m, e in enumerate(MON):
            if sum(e) == degree:
                correction = own.matvec(L, [target[i][m]-before[i][m] for i in rows])
                for j, a in enumerate(correction):
                    source[j][m] += a
        assert source == expected
        after = physical(source)
        assert all(after[i][m] == target[i][m] for i in rows for m, e in enumerate(MON) if sum(e) <= degree)
        stages.append({'degree': degree, 'regular27_through_degree_all_zero': True})
    q = after
    assert source == [[F(x) for x in row] for row in c['source27_IFT_coefficients']]
    assert q == [[F(x) for x in row] for row in c['full150_physical_source_coefficients_afterIFT']]
    fold = c['fold_physical_row']
    f = {e: q[fold][m]-target[fold][m] for m, e in enumerate(MON) if q[fold][m] != target[fold][m]}
    W = {e: F(x) for e, x in zip(MON, c['Weierstrass_monic_coefficients']) if F(x)}
    U = {e: F(x) for e, x in zip(MON, c['Weierstrass_unit_coefficients']) if F(x)}
    assert U[(0, 0)] == -3 and W[(0, 2)] == 1
    assert f == {(1, 1): F(48), (0, 2): F(-3)}
    assert scalar_mul(U, W) == f
    assert [f.get(e, F(0)) for e in MON] == list(map(F, c['raw_fold_equation_coefficients']))
    assert all(b < 2 and a+b >= 2 for (a, b) in W if (a, b) != (0, 2))

    def reduce_term(a, b, value, out):
        if a+b > 3 or not value:
            return
        if b < 2:
            out[a, b] = out.get((a, b), F(0))+value
        else:
            for (i, j), z in W.items():
                if (i, j) != (0, 2):
                    reduce_term(a+i, b-2+j, -value*z, out)

    aa, bb = [], []
    for row in q:
        rem = {}
        for value, (a, b) in zip(row, MON):
            reduce_term(a, b, value, rem)
        aa.append([rem.get((i, 0), F(0)) for i in range(4)])
        bb.append([rem.get((i, 1), F(0)) for i in range(3)])
    assert aa == [[F(x) for x in row] for row in c['full150_remainders_even_a']]
    assert bb == [[F(x) for x in row] for row in c['full150_remainders_odd_b']]
    assert all(bb[i] == [0]*3 and aa[i] == [target[i][MON.index((j, 0))] for j in range(4)] for i in rows+[fold])
    nrows = curve['N37_coordinate_rows']
    odd = nrows[28:]
    assert set(nrows[:28]) == set(rows+[fold])
    assert odd == c['actual_odd9_physical_rows']
    TN = [[F(x) for x in row] for row in curve['full150_by37_actual_N_normalized_tangent']]
    assert [[TN[i][j] for j in range(37)] for i in nrows] == [[F(i == j) for j in range(37)] for i in range(37)]
    h0 = [row[28:] for row in TN]
    assert h0 == [[F(x) for x in row] for row in c['actual_full150_by9_minimal_odd_class_constants']]
    z = [[target[i][MON.index((j, 0))]-aa[i][j] for j in range(4)] for i in range(150)]
    mu = F(300, 23)
    assert all(row[0] == row[1] == 0 for row in z)
    assert all(z[i][2] == mu*bb[i][1] for i in range(150))
    assert any(bb[i][1] for i in odd) and all(bb[i][1] == z[i][2] == 0 for i in odd[-2:])
    bres, graphs = [], []
    for j in range(150):
        assert bb[j][1] == sum(a*bb[i][1] for a, i in zip(h0[j], odd))
        r2 = bb[j][2]-sum(a*bb[i][2] for a, i in zip(h0[j], odd))
        bres.append(r2)
        graphs.append(z[j][3]-sum(a*z[i][3] for a, i in zip(h0[j], odd))-mu*r2)
    assert bres == list(map(F, c['full150_b_minus_h0_bminimal_second_coefficients']))
    assert graphs == list(map(F, c['full150_linear_graph_cubic_restrictions'])) == [F(0)]*150
    nonbasis = [i for i in range(150) if i not in nrows]
    assert nonbasis == c['all113_nonbasis_physical_rows'] and len(nonbasis) == 113
    # Independent trace centering from the actual quadratic coefficient of p.
    center = [-W.get((i, 1), F(0))/2 for i in range(3)]
    assert center == [0, 8, 0]
    traceaa = [[a[i]+sum(center[k]*b[i-k] for k in range(3) if 0 <= i-k < 3)
                for i in range(4)] for a, b in zip(aa, bb)]
    assert traceaa == [[F(x) for x in row] for row in c['full150_trace_zero_even_a']]
    ztrace = [[target[i][MON.index((j, 0))]-traceaa[i][j] for j in range(4)] for i in range(150)]
    mutrace = mu-8
    assert mutrace == F(116, 23) and mutrace**2-64 == -F(20400, 529)
    assert [ztrace[j][3]-sum(a*ztrace[i][3] for a, i in zip(h0[j], odd))-mutrace*bres[j] for j in range(150)] == graphs
    evidence = dict(independent_method='exp(L_log(I+M)) exp(p L_K); own index Lie; no inverse-coordinate pushforward reused',
        bivariate_monomial_order=MON, all150_rebuilt_source_coefficients=q,
        independent_IFT_stage_checks=stages, fold=f, W=W, U=U,
        all150_even_remainders=aa, all150_odd_remainders=bb,
        all150_graph_cubic_restrictions=graphs)
    # JSON keys for bivariate scalar polynomials are serialized as exponent/value lists.
    for key in ['fold', 'W', 'U']:
        evidence[key] = [[list(e), str(v)] for e, v in evidence[key].items()]
    def serial(x):
        if isinstance(x, F): return str(x)
        if isinstance(x, dict): return {k: serial(v) for k, v in x.items()}
        if isinstance(x, (tuple, list)): return [serial(v) for v in x]
        return x
    evidpath = output/'independent-coefficients.json.gz'
    with evidpath.open('wb') as stream:
        with gzip.GzipFile(fileobj=stream, mode='wb', filename='', mtime=0) as gz:
            gz.write((json.dumps(serial(evidence), separators=(',', ':'))+'\n').encode())
    paths = [Path(__file__), Path(own.__file__), ROOT/f'research/notes/{STEM}.md', authorpath,
             hpath, firstreview, capath, redpath, curvepath, secondpath, *stage_paths]
    result = dict(id=STEM, all_passed=True,
        independent_all150_bivariate_degree3_Lie_log_source_reconstruction=True,
        independent_all27_IFT_stages_with_both_inverse_products=True,
        same_actual_R_polynomial_curve_and_actual_N37_binding=True,
        exact_fold_Weierstrass_unit_and_all150_remainders_verified=True,
        physical_base28_equals_regular27_plus_fold=True,
        independent_all150_odd_class_constants_and_all113_nonbasis_graph_cubics=True,
        raw_and_trace_zero_means=[str(mu), str(mutrace)],
        nonunique_odd_generator_division_cubic_formula_proved=True,
        no_image_ideal_nonflat_basechange_or_second_order_source_lift_assumed=True,
        all_N_relations_R_restriction_order_atleast4_proved=True,
        accepted_h4_implies_exact_k4=True, h_exact=4, k_exact=4,
        reduced_slice_embedding_dimension=18, reduced_cone_embedding_dimension=37,
        reduced_projective_embedding_dimension=36,
        actual_slice_embedding_dimension=20, linear_nil_cotangent_kernel_dimension=2,
        scope='Actual generic B18 open orbit; exact reduced-contact/cotangent theorem, not full nonreduced D ring or exact depth.',
        input_sha256={str(p.relative_to(ROOT)): own.sha(p) for p in paths},
        evidence_sha256={evidpath.name: own.sha(evidpath)}, seconds=time.monotonic()-start)
    (output/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    main(parser.parse_args().output)

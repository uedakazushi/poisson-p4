#!/usr/bin/env python3
"""Actual W=1 tangent planes from genuine kernel mixed2 and one kernel3 per source."""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import resource
import signal
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_W1_ACTUAL_TANGENT_INTERSECTION_20261103_v63'


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def enc(A):
    if isinstance(A, (list, tuple)):
        return [enc(a) for a in A]
    return str(A)


def mm(A, B):
    return [[sum((a*b for a, b in zip(row, col) if a and b), F(0)) for col in zip(*B)] for row in A]


def dump(p, value):
    with p.open('wb') as f:
        with gzip.GzipFile(fileobj=f, mode='wb', filename='', mtime=0) as g:
            g.write((json.dumps(value, separators=(',', ':'))+'\n').encode())


def main(out):
    started = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    signal.alarm(60)
    out.mkdir(parents=True, exist_ok=False)
    w_path = ROOT/'reports/AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63/staircase.json.gz'
    allw_path = ROOT/'reports/AMBIENT_B18_ALL_NONZERO_W_FULL_COMPLETIONS_20261103_v63/report.json'
    reviews = [ROOT/'reports/AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63_REVIEW.json',
               ROOT/'reports/AMBIENT_B18_W_UNIFORM_DARBOUX_20261103_v63_REVIEW/report.json',
               ROOT/'reports/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63_REVIEW.json']
    assert read(allw_path)['all_passed'] and all(read(p)['all_passed'] for p in reviews)
    ww = read(w_path)
    h, x, y, z, w, W = s.symbols('h x y z w W')
    X = s.Matrix([h, x, y, z, w])
    pairs = list(itertools.combinations(range(5), 2))
    mons = list(itertools.combinations_with_replacement(range(5), 2))
    pindex, mindex = {p: i for i, p in enumerate(pairs)}, {m: i for i, m in enumerate(mons)}
    def vec(P):
        return [F(s.Poly(P[i, j], *X).coeff_monomial(X[a]*X[b])) for i, j in pairs for a, b in mons]
    def canonical(P):
        d = s.Matrix([sum(s.diff(P[i, j], X[i]) for i in range(5)) for j in range(5)])
        return (P-(X*d.T-d*X.T)/5).applyfunc(s.expand)
    def raw(f, g, H, cb):
        P = s.zeros(5)
        P[1, 2], P[3, 4] = f, g if cb else g*h*h
        P[2, 3] = (-s.diff(f, x)*z/2 if cb else 0)-s.diff(H, w)
        P[2, 4] = (-s.diff(f, x)*w/2 if cb else 0)+s.diff(H, z)
        return P-P.T
    sw = s.eye(5)
    sw[0, 0] = sw[1, 1] = 0
    sw[0, 1] = sw[1, 0] = 1
    change = lambda P: (sw*P.subs({h: x, x: h}, simultaneous=True)*sw.T).applyfunc(s.expand)
    cubics = [s.prod((h, x, z, w)[i] for i in m) for m in itertools.combinations_with_replacement(range(4), 3) if any(i >= 2 for i in m)]
    PW = s.Matrix(ww['CA']['physical_entire_W_curve']).subs(W, 1)
    point = vec(PW)
    assert s.Matrix(ww['CB']['physical_entire_W_curve']).subs(W, 1) == PW
    # A literal rank4 value excludes the whole generic-rank2 component R.
    value = PW.subs({h: 0, x: 1, y: 0, z: 0, w: 0})
    assert value.extract([0, 2, 3, 4], [0, 2, 3, 4]).det() == 1
    def entry(V, i, j, m):
        if i == j:
            return F(0)
        return V[15*pindex[(i, j)]+m] if i < j else -V[15*pindex[(j, i)]+m]
    # Independent coefficient-index Lie action on all150 coordinates.
    def lie(M, V):
        outv = [F(0)]*150
        for k, (i, j) in enumerate(pairs):
            for m in range(15):
                outv[15*k+m] += sum(M[i][a]*entry(V, a, j, m)+M[j][a]*entry(V, i, a, m) for a in range(5))
            for m, (a, b) in enumerate(mons):
                q = V[15*k+m]
                if not q:
                    continue
                for u, other in [(a, b), (b, a)]:
                    for v in range(5):
                        if M[u][v]:
                            outv[15*k+mindex[tuple(sorted((other, v)))]] -= q*M[u][v]
        return outv
    # Independent literal inverse-coordinate and exterior-power source jet.
    def physical(M, delta, Q, order):
        n = order+1
        zero, one = (F(0),)*n, (F(1),)+(F(0),)*(n-1)
        def ja(a, b):
            return tuple(x+y for x, y in zip(a, b))
        def jm(a, b):
            if a == zero or b == zero:
                return zero
            outj = [F(0)]*n
            for i, q in enumerate(a):
                if q:
                    for j in range(n-i):
                        if b[j]:
                            outj[i+j] += q*b[j]
            return tuple(outj)
        def js(values):
            v = zero
            for q in values:
                v = ja(v, q)
            return v
        powers = [[[F(i == j) for j in range(5)] for i in range(5)]]
        for _ in range(order):
            powers.append(mm(powers[-1], M))
        G = [[tuple(powers[k][i][j]/math.factorial(k) for k in range(n)) for j in range(5)] for i in range(5)]
        Gi = [[tuple((-1)**k*powers[k][i][j]/math.factorial(k) for k in range(n)) for j in range(5)] for i in range(5)]
        P = [tuple([a, b, c]+[F(0)]*(n-3)) for a, b, c in zip(point, delta, Q)]
        quad = []
        for a, b in mons:
            rr = []
            for i, j in mons:
                q = jm(Gi[a][i], Gi[b][j])
                if i != j:
                    q = ja(q, jm(Gi[a][j], Gi[b][i]))
                rr.append(q)
            quad.append(rr)
        changed = [[js(jm(P[15*k+l], quad[l][m]) for l in range(15)) for m in range(15)] for k in range(10)]
        wedge = [[tuple(x-y for x, y in zip(jm(G[i][a], G[j][b]), jm(G[i][b], G[j][a]))) for a, b in pairs] for i, j in pairs]
        return [js(jm(wedge[i][k], changed[k][m]) for k in range(10)) for i in range(10) for m in range(15)]
    records, bases, inputs = {}, {}, []
    for name, expected in [('CA', 38), ('CB', 38), ('N', 36)]:
        cb = name == 'CB'
        critical_path = ROOT/f'reports/AMBIENT_B18_{name}_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
        inputs.append(critical_path)
        flags = read(critical_path)['flag_positions']
        generators = []
        for i, j in flags:
            M = [[F(0)]*5 for _ in range(5)]
            M[i][j] = F(1)
            generators.append(M)
        C = [vec(change(canonical(raw(f, 0, 0, cb)))) for f in [h*h, h*x, x*x]]
        CG = [vec(change(canonical(raw(0, f, 0, cb)))) for f in [h*h, h*x, x*x]] if cb else []
        C += [vec(change(canonical(raw(0, h*h if cb else 1, 0, cb))))]
        C += [vec(change(canonical(raw(0, 0, H, cb)))) for H in cubics]
        alpha = [F(0)]*20
        alpha[0] = alpha[7] = -1
        alpha[3] = alpha[11] = 1
        alpha[16] = F(1, 3)
        if cb:
            actual = [sum(alpha[i]*C[i][j] for i in [0, 1, 2, *range(4, 20)])
                      -alpha[3]*sum(alpha[i]*CG[i][j] for i in range(3)) for j in range(150)]
            fixed = [[C[i][j]-alpha[3]*CG[i][j] for j in range(150)] for i in range(3)]
            fixed += [[-sum(alpha[i]*CG[i][j] for i in range(3)) for j in range(150)]]
            fixed += C[4:]
            indices = list(range(20))
        else:
            actual = [sum(a*c[j] for a, c in zip(alpha, C)) for j in range(150)]
            indices = [i for i in range(20) if name != 'N' or i != 2]
            fixed = [C[i] for i in indices]
        assert actual == point
        nf = len(fixed)
        columns = fixed+[lie(M, point) for M in generators]
        J = s.Matrix.hstack(*(s.Matrix(q) for q in columns))
        kernel = J.nullspace()
        assert len(kernel) == 1 and J.rank() == nf+8
        k = kernel[0]
        natural = nf+flags.index([1, 0])
        if k[natural] == 0:
            natural = next(i for i in range(nf, nf+9) if k[i])
        k = k/k[natural]
        kv = list(map(F, k))
        M = [[sum(kv[nf+q]*E[i][j] for q, E in enumerate(generators)) for j in range(5)] for i in range(5)]
        delta = [sum(a*Fv[j] for a, Fv in zip(kv[:nf], fixed)) for j in range(150)]
        assert all(a+b == 0 for a, b in zip(lie(M, point), delta))
        fullvel = [F(0)]*20
        for a, i in zip(kv[:nf], indices):
            fullvel[i] = a
        Q = [-fullvel[3]*sum(fullvel[i]*CG[i][j] for i in range(3)) for j in range(150)] if cb else [F(0)]*150
        lp = lie(M, point)
        l2p = lie(M, lp)
        W2 = [a/2+b+c for a, b, c in zip(l2p, lie(M, delta), Q)]
        mixed = []
        for j, Fv in enumerate(fixed):
            q = lie(M, Fv)
            if cb:
                ii = indices[j]
                if ii < 3:
                    q = [a-fullvel[3]*b for a, b in zip(q, CG[ii])]
                elif ii == 3:
                    q = [a-sum(fullvel[i]*CG[i][l] for i in range(3)) for l, a in enumerate(q)]
            mixed.append(q)
        for E in generators:
            a, b, d = lie(E, lp), lie(M, lie(E, point)), lie(E, delta)
            mixed.append([x/2+y/2+z for x, y, z in zip(a, b, d)])
        # All mixed2 columns are independently checked by direct physical
        # polarization. Only the true kernel curve below has a cubic output.
        for j in range(nf+9):
            vv = list(kv)
            vv[j] += 1
            MM = [[sum(vv[nf+q]*E[i][l] for q, E in enumerate(generators)) for l in range(5)] for i in range(5)]
            dd = [sum(a*Fv[l] for a, Fv in zip(vv[:nf], fixed)) for l in range(150)]
            fv = [F(0)]*20
            for a, i in zip(vv[:nf], indices):
                fv[i] = a
            qq = [-fv[3]*sum(fv[i]*CG[i][l] for i in range(3)) for l in range(150)] if cb else [F(0)]*150
            second = [a[2] for a in physical(MM, dd, qq, 2)]
            pure = ([F(0)]*150 if j < nf else [a/2 for a in lie(generators[j-nf], lie(generators[j-nf], point))])
            assert [a-b-c for a, b, c in zip(second, W2, pure)] == mixed[j]
        H = s.Matrix.hstack(*(s.Matrix(q) for q in mixed))
        S2 = J.row_join(s.Matrix(W2)).row_join(H)
        rank2 = S2.rank()
        assert rank2 in [expected-1, expected]
        cube = None
        if rank2 == expected-1:
            physical_jet = physical(M, delta, Q, 3)
            assert [v[0] for v in physical_jet] == point and all(v[1] == 0 for v in physical_jet)
            assert [v[2] for v in physical_jet] == W2
            cube = [v[3] for v in physical_jet]
            lie_cube = [a/6+b/2+c for a, b, c in zip(lie(M, l2p), lie(M, lie(M, delta)), lie(M, Q))]
            assert cube == lie_cube
            full = S2.row_join(s.Matrix(cube))
        else:
            full = S2
        cols = list(full.rref()[1])
        assert len(cols) == expected
        basis = full[:, cols]
        rows = list(basis.T.rref()[1])
        determinant = s.factor(basis[rows, :].det())
        assert determinant != 0
        bases[name] = basis
        record = dict(source_fixed20_indices=indices, source_flag_positions=flags,
                      source_point20=list(map(str, alpha)), source_Jacobian150_by_n=enc(J.tolist()),
                      source_kernel=list(map(str, k)), source_kernel_GL_matrix=enc(M),
                      source_kernel_fixed_velocity20=enc(fullvel), full_CB_source_quadratic_velocity150=enc(Q),
                      full150_kernel_square=enc(W2), full150_by_n_mixed2=enc(H.tolist()),
                      full150_kernel_cube=None if cube is None else enc(cube),
                      source2_span_rank=rank2, one_kernel_cube_used=cube is not None,
                      actual_tangent_dimension=expected, independent_column_indices=cols,
                      actual_tangent150_by_dimension=enc(basis.tolist()),
                      unit_minor_rows=rows, unit_minor_determinant=str(determinant))
        records[name] = record
        dump(out/f'{name}-tangent.json.gz', record)
    A, B, N = [bases[name] for name in ['CA', 'CB', 'N']]
    assert A.row_join(N).rank() == 38 and B.row_join(N).rank() == 38
    union_rank = A.row_join(B).rank()
    dimD = 76-union_rank
    ar = list(A.T.rref()[1])
    normalizedA = A*A[ar, :].inv()
    rem = B-normalizedA*B[ar, :]
    kernel = rem.nullspace()
    TD = B*s.Matrix.hstack(*kernel)
    assert TD.rank() == dimD
    evidence = dict(physical_point_matrix=enc(PW.tolist()), physical_point150=enc(point),
                    full150_actual_D_tangent_basis=enc(TD.tolist()),
                    full150_CA_by_CB_membership_residual_matrix=enc(rem.tolist()),
                    source_families=records)
    dump(out/'certificate.json.gz', evidence)
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', w_path, allw_path, *reviews, *inputs]
    result = dict(id=STEM, all_passed=True, W=1, source_fields='QQ',
                  actual_CA_CB_N_tangent_dimensions=[38, 38, 36],
                  source2_span_dimensions=[records[name]['source2_span_rank'] for name in ['CA', 'CB', 'N']],
                  independent_physical_polarization_all_kernel_mixed2=True,
                  pure_kernel_cube_count=sum(records[name]['one_kernel_cube_used'] for name in records),
                  all150_true_kernel_cubic_coefficients_independently_physical_and_Lie_match=True,
                  CB_full_nonlinear_ratio_times_f_second_source_term_retained=True,
                  actual_CA_plus_CB_tangent_span_dimension=union_rank,
                  actual_D_tangent_dimension=dimD, actual_N_tangent_dimension=36,
                  actual_D_reduced_germ_equals_N_by_rank4_value_and_global_N_union_R=True,
                  actual_D_nonreduced_at_W1=dimD > 36,
                  nilradical_cotangent_image_dimension=dimD-36,
                  full_D_type_depth_nil_support_dimension_or_new_embedded_component_claimed=False,
                  no_multivariate_source3_or_higher_source_jet_or_GB=True,
                  input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
                  evidence_sha256={p.name: sha(p) for p in [out/'certificate.json.gz', *[out/f'{name}-tangent.json.gz' for name in records]]},
                  seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'source2_span_dimensions',
          'actual_CA_plus_CB_tangent_span_dimension', 'actual_D_tangent_dimension',
          'nilradical_cotangent_image_dimension', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())

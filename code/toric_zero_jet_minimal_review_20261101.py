#!/usr/bin/env python3
"""Independent strict-action and full-matrix degeneration review of M0/M1."""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
import sympy as s
from s61_model import X, quotient_coefficients
from weight1_structure import J

ROOT = Path(__file__).resolve().parents[1]


def verify():
    prior_path = ROOT / 'reports/TORIC_ZERO_JET_DIMENSION_REVIEW_20261101.json'
    source_path = ROOT / 'reports/TORIC_ZERO_JET_MINIMAL_20261101.json'
    prior = json.loads(prior_path.read_text())
    source = json.loads(source_path.read_text())
    h = s.symbols('h0:9')
    x0, x, y, z, w = X
    k, lam, t, alpha = s.symbols('k lam t alpha')
    f0, f1, f2 = s.symbols('f0 f1 f2')
    Y = s.Matrix([x, y, z, w])
    F = f0*x*x+f1*x*w+f2*w*w
    result = {}
    def elem(i, j):
        A = s.zeros(5); A[i, j] = 1
        return A
    for kind, eps in [('M0', 1), ('M1', -1)]:
        # Independent hand-written conformal action from the earlier review:
        # strict action is its subspace p=-m, with no tensor-scaling direction.
        M = s.Matrix(prior['minimal_models'][kind]['residual_potential_action_matrix'])
        N = s.Matrix.hstack(M[:, 0]-M[:, 3], M[:, 1], M[:, 2], M[:, 4], M[:, 5])
        D = N.subs(dict(zip(h[5:], [lam, 0, 0, lam])))
        dm = s.factor(D[[2, 5, 6, 7], [0, 1, 2, 3]].det())
        assert dm != 0 and not (dm.free_symbols - {lam})
        D = N.subs(dict(zip(h[5:], [0, 1, 0, 0])))
        bm = s.factor(D[[2, 3, 6, 7], [0, 1, 3, 4]].det())
        assert bm != 0 and not bm.free_symbols
        D = N.subs(dict(zip(h[5:], [1, 0, 0, 0])))
        rr, cc = [2, 5, 6], [0, 1, 3]
        assert D[rr, cc].det() == -27
        rem = (D[:, [2, 4]]-D[:, cc]*D[rr, cc].inv()*D[rr, [2, 4]]).applyfunc(s.factor)
        gb = s.groebner(list(rem), *h[:5])
        assert list(gb) == [1]
        rank1 = N.subs(dict(zip(h[2:], [1, 0, 0, 0, 0, 0, 0])))
        one_minor = s.factor(rank1[[0, 1, 2, 3], [0, 1, 3, 4]].det())
        assert one_minor != 0 and not one_minor.free_symbols
        semisimple = N.subs(dict(zip(h, [0, 0, 0, k, 0, 0, 0, 0, 0])))
        sm = s.factor(semisimple[[0, 1, 2, 4], [1, 2, 3, 4]].det())
        assert s.cancel(sm/(k*k*(k*k-1))) in [1, -1]
        resonance_minor = None
        if kind == 'M0':
            resonance = N.subs(dict(zip(h, [0, k, 0, 1, 0, 0, 0, 0, 0])))
            found = []
            for rr2 in combinations(range(5), 4):
                for cc2 in combinations(range(5), 4):
                    det = s.factor(resonance[list(rr2), list(cc2)].det())
                    if det and not s.cancel(det/k).free_symbols:
                        found.append((rr2, cc2, det))
            assert found
            resonance_minor = str(found[0][2])
        # Full sl5 kernel, including projective directions: do not infer it from
        # the displayed affine subalgebra alone.
        G = x*z*w if kind == 'M0' else s.S.Zero
        P = s.zeros(5)
        P[1, 2] = x0*x0+x*x
        P[3, 4] = x0*x0+(x*x if kind == 'M1' else 0)
        P[2, 3] = -s.diff(G, w)-(x*z if kind == 'M1' else 0)
        P[2, 4] = s.diff(G, z)-(x*w if kind == 'M1' else 0)
        P -= P.T
        columns = []
        for i in range(5):
            for j in range(5):
                A = elem(i, j)
                columns.append(s.Matrix(quotient_coefficients(A*P+P*A.T-X[j]*P.diff(X[i]))))
        action = s.Matrix.hstack(*columns)
        basis = [s.Matrix(a) for a in source['minimal_families'][kind]['exceptional_strict_basis']]
        expected = 20 if kind == 'M0' else 19
        assert action.rank() == expected
        assert s.Matrix.hstack(*[a.reshape(25, 1) for a in basis+[s.eye(5)]]).rank() == 25-expected
        assert all(action*a.reshape(25, 1) == s.zeros(action.rows, 1) for a in basis+[s.eye(5)])
        if kind == 'M0':
            assert basis[3]*basis[2]-basis[2]*basis[3] == basis[2]/2
        else:
            assert all(basis[i]*basis[j]-basis[j]*basis[i] == s.zeros(5) for i in range(2) for j in range(5))
            assert basis[3]*basis[4]-basis[4]*basis[3] == 2*basis[2]
            assert basis[2]*basis[3]-basis[3]*basis[2] == basis[3]
            assert basis[2]*basis[4]-basis[4]*basis[2] == -basis[4]
        # Reconstruct the entire cotangent-lift matrix in original coordinates
        # and pull back the full AB bivector, retaining all three F coefficients.
        A = s.Matrix([[s.I, 0], [alpha*t, t]]) if kind == 'M0' else s.Matrix([[s.I, 0], [s.I-alpha*t/2, t]])
        pos = A*s.Matrix([x, w]); mom = A.inv().T*s.Matrix([y, -z])
        new = s.Matrix([pos[0], mom[0], pos[1], mom[1]])
        C = new.jacobian(Y)
        assert (C*J*C.T-J).applyfunc(s.simplify) == s.zeros(4)
        AB = J.copy()
        AB[0, 1] -= pos[0]**2; AB[1, 0] += pos[0]**2
        AB[2, 3] -= pos[1]**2; AB[3, 2] += pos[1]**2
        AB[1, 3] += F/A.det(); AB[3, 1] -= F/A.det()
        pulled = (C.inv()*AB*C.inv().T).applyfunc(s.cancel)
        assert all(s.denom(v).subs(t, 0) != 0 for v in pulled)
        target = J.copy()
        target[0, 1] += x*x; target[1, 0] -= x*x
        target[1, 2] -= F; target[2, 1] += F
        target[1, 3] += alpha*x*x; target[3, 1] -= alpha*x*x
        if kind == 'M1':
            target[2, 3] += x*x; target[3, 2] -= x*x
            target[1, 3] -= 2*x*w; target[3, 1] += 2*x*w
        assert (pulled.subs(t, 0)-target).applyfunc(s.expand) == s.zeros(4)
        result[kind] = dict(alternative_strict_action=[[str(v) for v in row] for row in N.tolist()],
            distinct_minor=str(dm), double_minor=str(bm), triple_schur_groebner_basis=['1'],
            rank1_minor=str(one_minor), rank2_minor=str(sm), resonance_minor=resonance_minor,
            full_strict_action_rank=expected, full_symplectic_deformation_checked=True)
    paths = [Path(__file__), prior_path, source_path, ROOT/'code/s61_model.py', ROOT/'code/weight1_structure.py']
    return dict(id='TORIC_ZERO_JET_MINIMAL_REVIEW_20261101', all_passed=True,
        reviewed_cases=result, conclusion='F0 intersect (M0 union M1) is contained in closed AB as reduced sets.',
        proof_note='research/notes/TORIC_ZERO_JET_MINIMAL_REVIEW_20261101.md',
        input_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); result = verify()
    with args.output.open('x') as out:
        json.dump(result, out, indent=2); out.write('\n')
    print(json.dumps(dict(all_passed=True)))


if __name__ == '__main__':
    main()

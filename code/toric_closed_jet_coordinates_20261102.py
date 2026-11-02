#!/usr/bin/env python3
"""Linear Sp4 coordinates on the primitive first jet and the closed projection of the second jet."""
import itertools, json
from pathlib import Path
import sympy as s
ROOT = Path(__file__).resolve().parents[1]
from s61_model import S61, X
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import J
pairs = list(itertools.combinations(range(4), 2))
m2 = list(itertools.combinations_with_replacement(range(4), 2))
V = [(1, 0), (-1, 0), (0, 1), (0, -1)]

def coords():
    """Build the actual S61 jet basis, closed projection, and simple-root actions."""
    model = S61()
    z, u, beta = primitive_beta(model)
    lin = -J * beta * J
    lm = s.Matrix([s.Poly(lin[i, j], *z).coeff_monomial(z[k]) for i, j in pairs for k in range(4)]).jacobian(u)
    li = list(lm.T.to_DM().rref()[1])
    G = lm * lm[li, :].inv()
    q = s.symbols('q:45')
    it = iter(q)
    P = model.reconstruct([0 if w == 1 else next(it) for w in model.reduced['parameter_weights']]).subs(X[0], 1)
    y = X[1:]
    ym = [y[a] * y[b] for a, b in m2]
    E = []
    for i, j in pairs:
        f = s.Poly(P[i + 1, j + 1] - y[i] * P[0, j + 1] + y[j] * P[0, i + 1] - J[i, j], *y)
        E.extend((f.coeff_monomial(m) for m in ym))
    E = s.Matrix(E).jacobian(q)
    sel = json.load(open(ROOT / 'reports/AB_ORIGIN_DEGREE_FOUR_20261101.json'))['coordinate_indices']
    idx = {a: k for k, a in enumerate(pairs)}

    def mapJ(i, j):
        a = next((a for a in range(4) if J[a, i]))
        b = next((b for b in range(4) if J[j, b]))
        return ((min(a, b), max(a, b)), int(J[a, i] * J[j, b]) * (1 if a < b else -1))
    pr = s.zeros(60)
    for h, (i, j) in enumerate(pairs):
        (a, b), sg = mapJ(i, j)
        for m, mon in enumerate(m2):
            for sign, r, l in [(1, a, b), (-1, b, a)]:
                word = sorted((r,) + mon)
                for k in set(word):
                    if k == l:
                        continue
                    rem = word.copy()
                    rem.remove(k)
                    ij, sg2 = mapJ(min(k, l), max(k, l))
                    pr[10 * idx[ij] + m2.index(tuple(rem)), 10 * h + m] += s.Rational(sg * sg2 * sign * (1 if k < l else -1) * word.count(k), 4)
    assert pr * pr == pr and pr * E == E and (pr.rank() == 45)
    Pi = E[sel, :].inv() * pr[sel, :]
    assert Pi * E == s.eye(45) and E * Pi == pr
    roots = []
    for spec in [[(0, 2, 1), (3, 1, -1)], [(2, 0, 1), (1, 3, -1)], [(2, 3, 1)], [(3, 2, 1)]]:
        A = s.zeros(4)
        for a, b, c in spec:
            A[a, b] = c
        assert A * J + J * A.T == s.zeros(4)
        roots.append(A)

    def rep(A, monlist):
        n = len(monlist)
        U = s.zeros(6 * n)
        for h, (i, j) in enumerate(pairs):
            for m, word in enumerate(monlist):
                for u in range(4):
                    for a, b, c in [(u, j, A[u, i]), (i, u, A[u, j])]:
                        if a != b and c:
                            U[n * idx[tuple(sorted((a, b)))] + m, n * h + m] += (1 if a < b else -1) * c
                for p, k in enumerate(word):
                    for l in range(4):
                        if A[k, l]:
                            U[n * h + monlist.index(tuple(sorted(word[:p] + (l,) + word[p + 1:]))), n * h + m] -= A[k, l]
        return U
    acts = []
    wm = lambda indexes: tuple((tuple((V[i][d] + V[j][d] - sum((V[k][d] for k in word)) for d in (0, 1))) for i, j in pairs for word in indexes))
    bw = [wm([(x,) for x in range(4)])[i] for i in li]
    qw = json.load(open(ROOT / 'reports/AB_EXPLICIT_QUADRICS_CUBIC_GENERATION_20261102.json'))['coordinate_weights']
    QW = wm(m2)
    assert all((QW[i] == tuple(qw[j]) for i in range(60) for j in range(45) if E[i, j]))
    for A in roots:
        L = rep(A, [(x,) for x in range(4)])
        Q = rep(A, m2)
        AL = (L * G)[li, :]
        AQ = Pi * Q * E
        assert L * G == G * AL and Q * E == E * AQ and (Q * pr == pr * Q)
        acts.append(s.diag(AL, AQ))
    return dict(li=li, G=G, E=E, Pi=Pi, pr=pr, weights=bw + list(map(tuple, qw)), acts=acts, roots=roots)

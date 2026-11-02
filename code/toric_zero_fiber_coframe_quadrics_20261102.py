#!/usr/bin/env python3
"""Certify all 274 AB lifts on the dense centered coframe using exact Sp4 stability."""
import argparse, hashlib, itertools, json, math, random, time
from pathlib import Path
from collections import defaultdict, deque
import numpy as np
import sympy as s
from sympy.polys.rings import ring
from toric_closed_jet_coordinates_20261102 import coords, pairs, m2
from weight1_structure import J as SJ
ROOT = Path(__file__).resolve().parents[1]
NOTE = 'research/notes/TORIC_ZERO_FIBER_COFRAME_QUADRICS_20261102.md'
J = np.array(SJ.tolist(), dtype=int).astype(object)
comb = itertools.combinations_with_replacement
base = list(comb(range(16), 4)) + [b + (q,) for b in comb(range(16), 2) for q in range(16, 61)]

def group(mons):
    g = defaultdict(list)
    for m in mons:
        g[tuple((sum((W[i][k] for i in m)) for k in (0, 1)))].append(m)
    return g

def action_lists(A):
    return [[(j, A[i, j]) for j in range(61) if A[i, j]] for i in range(61)]

def monact(m, spec):
    b = defaultdict(int)
    for k, v in enumerate(m):
        for j, c in spec[v]:
            b[tuple(sorted(m[:k] + (j,) + m[k + 1:]))] += c
    return {k: v for k, v in b.items() if v}

def highest(M):
    """Compute a basis of the joint kernel of the two raising operators."""
    if not M:
        return []
    cols = []
    for m in M:
        cols.append({(a, j): c for a in (1, 3) for j, c in monact(m, ACT[a]).items()})
    keys = {j: i for i, j in enumerate(sorted(set().union(*(r.keys() for r in cols))))}
    mat = s.MutableSparseMatrix(len(keys), len(M), {(keys[j], k): c for k, r in enumerate(cols) for j, c in r.items()})
    null = mat.to_DM().nullspace().to_Matrix()
    return [{m: c for m, c in zip(M, null.row(i)) if c} for i in range(null.rows)]

def deriv(q, h):
    r = defaultdict(int)
    for m, v in q.items():
        for k, c in monact(m, ACT[h]).items():
            r[k] += c * v
    return {m: v for m, v in r.items() if v}

def addcomb(rows, coef):
    r = defaultdict(int)
    for x, a in zip(rows, coef):
        if a:
            for m, c in x.items():
                r[m] += c * a
    return {m: v for m, v in r.items() if v}

def highest_spaces():
    """Split the AB quadrics and the allowable corrections by highest weight."""
    qr = json.load(open(ROOT / 'reports/AB_EXPLICIT_QUADRICS_CUBIC_GENERATION_20261102.json'))
    rels = [{(i + 16, j + 16): s.Rational(c) for i, j, c in row} for row in qr['quadratic_relations']]
    bg = group(base)
    rr = defaultdict(list)
    for r in rels:
        rr[tuple((sum((W[i][k] for i in next(iter(r)))) for k in (0, 1)))].append(r)
    out = {}
    total = 0
    for w, R in sorted(rr.items()):
        a, b = w
        if not a >= b >= 0:
            continue
        cols = [{(h, k): v for h in (1, 3) for k, v in deriv(r, h).items()} for r in R]
        K = {j: i for i, j in enumerate(set().union(*(x.keys() for x in cols)))}
        A = s.MutableSparseMatrix(len(K), len(R), {(K[j], k): c for k, r in enumerate(cols) for j, c in r.items()})
        kn = A.to_DM().nullspace().to_Matrix()
        if not kn.rows:
            continue
        top = [addcomb(R, kn.row(i)) for i in range(kn.rows)]
        B = highest(bg[w])
        dim = (a - b + 1) * (b + 1) * (a + 2) * (a + b + 3) // 6
        total += dim * len(top)
        out[w] = (B, top)
    assert total == 274, total
    return out

def jets(M, adj, det):
    """Evaluate the first and second jet, scaled by 10 det(R) to their respective weights."""
    A = []
    V = []
    E = np.eye(4, dtype=int).astype(object)
    colsum = M.sum(axis=0)
    for k in range(4):
        ds = M[:, k] * 5 - colsum[k]
        dm = ds[:, None] * M
        V.append(dm.sum(axis=0))
        A.append(adj @ dm - E[:, k, None] * colsum[None, :] * det)
    lt = [T @ J + J @ T.T for T in A]
    ll = [lt[k][i, j] * 2 for i, j in pairs for k in range(4)]
    q = []
    for a, b in m2:
        H = -E[:, a, None] * V[b][None, :]
        T = A[a] @ J @ A[b].T
        if a != b:
            H -= E[:, b, None] * V[a][None, :]
            T += A[b] @ J @ A[a].T
        h = H * det * det
        q.append(T + h @ J + J @ h.T)
    qq = np.array([q[k][i, j] for i, j in pairs for k in range(10)], dtype=object)
    return [ll[i] for i in sel] + list(P4 @ qq)

def ev(f, x):
    return sum((c * math.prod((x[i] for i in m)) for m, c in f.items()))

def build_lifts(HW):
    """Use rational sample evaluations only to propose the fifteen corrected vectors."""
    rng = random.Random(22798)
    pts = []
    while len(pts) < 26:
        mat = s.Matrix(4, 4, [rng.randrange(-2, 3) for _ in range(16)])
        dt = int(mat.det())
        if not dt:
            continue
        arr = lambda m: np.array(m.tolist(), dtype=int).astype(object)
        pts.append(jets(arr(mat), arr(mat.adjugate()), dt))
    SOL = {}
    ranks = {}
    for w, (bs, rs) in HW.items():
        mat = s.Matrix([[ev(f, x) for f in bs + rs] for x in pts])
        rr, piv = mat.to_DM().rref()
        rr = rr.to_Matrix()
        assert all((i < len(bs) for i in piv))
        sol = []
        for j, r in enumerate(rs):
            lift = addcomb([r] + [bs[i] for i in piv], [1] + [-rr[k, len(bs) + j] for k in range(len(piv))])
            assert not deriv(lift, 1) and (not deriv(lift, 3))
            sol.append(lift)
        ranks[str(w)] = len(piv)
        SOL[w] = sol
    return (SOL, ranks)

def orbit(S):
    """Select candidates by independence after restriction to the AB quadratic terms."""
    t = time.monotonic()
    P = 32003
    pe = lambda v: int(s.numer(v)) * pow(int(s.denom(v)), -1, P) % P
    queue = deque((x for rr in S.values() for x in rr))
    basis = {}
    res = []
    while queue:
        f = queue.popleft()
        qr = {m: c for m, c in f.items() if len(m) == 2}
        re = {m: pe(c) for m, c in qr.items() if pe(c)}
        while re:
            m = min(re)
            c = re[m]
            if m not in basis:
                break
            for k, v in basis[m].items():
                a = (re.get(k, 0) - c * v) % P
                if a:
                    re[k] = a
                else:
                    re.pop(k, None)
        if not re:
            continue
        t0 = pow(re[m], -1, P)
        basis[m] = {k: v * t0 % P for k, v in re.items()}
        res.append(f)
        for h in (0, 2):
            queue.append(deriv(f, h))
    assert len(res) == 274
    return res

def chartr():
    R, a, b, c, d, e, f = ring('a,b,c,d,e,f', s.ZZ)
    M = np.array([[1, 0, 0, 0], [0, a, 0, 0], [b, c, 1, 0], [d, e, 0, f]], dtype=object)
    adj = np.array([[a * f, 0, 0, 0], [0, f, 0, 0], [-b * a * f, -c * f, a * f, 0], [-d * a, -e, 0, a]], dtype=object)
    M = np.array([[R(x) for x in row] for row in M], dtype=object)
    adj = np.array([[R(x) for x in row] for row in adj], dtype=object)
    assert all((x == R.zero for x in (M @ adj - np.eye(4, dtype=int).astype(object) * (a * f)).flatten()))
    return (R, jets(M, adj, a * f))

def chart_check(F):
    """Verify polynomial identities coefficient by coefficient on the six-variable chart."""
    R, v = chartr()
    from functools import cache
    vl = [{sum((x << 5 * i for i, x in enumerate(m))): int(c) for m, c in r.items()} for r in v]
    assert all((all((max(m) * 4 < 32 for m in r)) for r in v))

    @cache
    def im(m):
        if len(m) == 1:
            return vl[m[0]]
        w = defaultdict(int)
        for a, c in im(m[:-1]).items():
            for b, d in vl[m[-1]].items():
                w[a + b] += c * d
        return {a: c for a, c in w.items() if c}
    for j, rel in enumerate(F):
        den = math.lcm(*(int(s.denom(x)) for x in rel.values()))
        w = defaultdict(int)
        for m, c in rel.items():
            n = int(c * den)
            for a, b in im(m).items():
                w[a] += b * n
        assert not any(w.values()), (j, [(a, b) for a, b in w.items() if b][:4])
    return im.cache_info().currsize

def primitive(f):
    den = math.lcm(*(int(s.denom(c)) for c in f.values()))
    r = {m: int(c * den) for m, c in f.items() if c}
    g = math.gcd(*r.values())
    return {m: c // g for m, c in r.items()}

def stability(F):
    """Express every root derivative in the selected basis, then recheck every coefficient."""
    B = defaultdict(list)
    T = defaultdict(list)

    def weight(f):
        m = next(iter(f))
        return tuple((sum((W[i][d] for i in m)) for d in (0, 1)))
    for f in F:
        f = primitive(f)
        B[weight(f)].append(f)
        for h in range(4):
            r = deriv(f, h)
            if r:
                T[weight(r)].append(primitive(r))
    n = 0
    for w, Ts in T.items():
        Bs = B[w]
        L = sorted({m for b in Bs for m in b if len(m) == 2})
        C = s.Matrix([[b.get(m, 0) for m in L] for b in Bs])
        piv = C.to_DM().rref()[1]
        assert len(piv) == len(Bs)
        D = s.Matrix([[b.get(L[q], 0) for q in piv] for b in Ts])
        A = (D * C[:, piv].inv()).tolist()
        for g, cs in zip(Ts, A):
            den = math.lcm(*(int(s.denom(c)) for c in cs))
            lhs = defaultdict(int, ((m, c * den) for m, c in g.items()))
            for f, c in zip(Bs, cs):
                if c:
                    c = int(c * den)
                    for m, a in f.items():
                        lhs[m] -= a * c
            assert not any(lhs.values())
        n += len(Ts)
    return n

def verify():
    """Rebuild the relations from the coframe and check specialization, invariance, and vanishing."""
    global W, ACT, P4, sel
    start = time.monotonic()
    d = coords()
    W = d['weights']
    ACT = [action_lists(A) for A in d['acts']]
    P4 = np.array((4 * d['Pi']).tolist(), dtype=int).astype(object)
    assert s.Matrix(P4.tolist()) == 4 * d['Pi']
    sel = d['li']
    H = highest_spaces()
    S, fit = build_lifts(H)
    F = orbit(S)
    stable = stability(F)
    mons = chart_check(F)
    qr = ROOT / 'reports/AB_EXPLICIT_QUADRICS_CUBIC_GENERATION_20261102.json'
    ref = [{(16 + i, 16 + j): s.Rational(c) for i, j, c in r} for r in json.loads(qr.read_text())['quadratic_relations']]
    by = defaultdict(list)
    lhs = defaultdict(list)
    wt = lambda m: tuple((sum((W[i][k] for i in m)) for k in (0, 1)))
    for r in ref:
        by[wt(next(iter(r)))].append(r)
    for f in F:
        r = {m: c for m, c in f.items() if len(m) == 2}
        lhs[wt(next(iter(r)))].append(r)
    for w, rows in by.items():
        cols = sorted(set().union(*(r.keys() for r in rows + lhs[w])))
        rank = lambda rs: s.Matrix([[r.get(m, 0) for m in cols] for r in rs]).to_DM().rank()
        assert rank(rows + lhs[w]) == rank(lhs[w]) == len(rows)
    assert len(F) == 274 and stable == 814 and (sum(map(len, S.values())) == 15)
    files = [Path(__file__), ROOT / 'code/toric_closed_jet_coordinates_20261102.py', ROOT / 'code/primitive_jet_elliptic_curve_20261101.py', ROOT / 'code/s61_model.py', ROOT / 'code/toric_coframe_fiber.py', ROOT / 'reports/TORIC_COFRAME_FIBER_20261031.json', ROOT / 'reports/TORIC_ZERO_JET_EXHAUSTION_20261101.json', qr, ROOT / 'reports/AB_QUADRICS_QUARTIC_MULTIPLICATION_20261102.json', ROOT / 'reports/AB_IDEAL_HIGH_GENERATOR_VANISHING_20261102.json', ROOT / 'reports/AB_ORIGIN_DEGREE_FOUR_20261101.json', ROOT / NOTE]
    return dict(id='TORIC_ZERO_FIBER_COFRAME_QUADRICS_20261102', all_passed=True, field='QQ', highest_space_dimensions={str(w): [len(b), len(r)] for w, (b, r) in H.items()}, sample_ranks_for_candidate_construction_only=fit, certified_independent_quadratic_lifts=len(F), nonzero_simple_root_images_checked_over_QQ=stable, zero_simple_root_images=4 * 274 - stable, exact_symbolic_chart_variables=6, cached_monomial_images=mons, actual_zero_fiber='equals reduced AB as schemes including the origin', global_S61_exhaustion=False, highest_lifts_in_projected_coordinates={str(w): [[[list(m), str(c)] for m, c in sorted(f.items())] for f in fs] for w, fs in S.items()}, input_sha256={str(f.relative_to(ROOT)): hashlib.sha256(f.read_bytes()).hexdigest() for f in files}, seconds=time.monotonic() - start)
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    r = verify()
    with args.output.open('x') as f:
        json.dump(r, f, ensure_ascii=False, separators=(',', ':'))
        f.write('\n')
    print('PASS', r['certified_independent_quadratic_lifts'], r['seconds'])

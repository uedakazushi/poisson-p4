#!/usr/bin/env python3
"""Full A=-1/3 affine centralizer normalization; no tangent/jet computation."""
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json'
DATA = ROOT / 'code/data/classification_nil22_top_sources.json'
OUT = ROOT / 'reports/CLASSIFICATION_NIL22_RESONANCE_NORMALIZATION_20261107.json'


def main():
    started = time.monotonic()
    x, y, z, w = q = S.symbols('x1 x2 x3 x4')
    n = S.symbols('n0:15')
    env = {str(t): t for t in (*q, *n)}
    branch = next(b for b in json.loads(DATA.read_text())['branches'] if b['name'] == 'A_minus_third')
    subst = {env[k]: S.sympify(v, locals=env) for k, v in branch['substitution'].items()}
    source = json.loads(INPUT.read_text())['cases'][4]['full_six_brackets']
    P = S.zeros(4)
    for (i, j), value in zip(itertools.combinations(range(4), 2), source):
        P[i, j] = S.expand(S.sympify(value, locals=env).subs(subst, simultaneous=True))
        P[j, i] = -P[i, j]

    def push(M, shift):
        inv = M.inv() * (S.Matrix(q) - shift)
        return (M * P * M.T).applyfunc(lambda f: S.expand(f.subs(dict(zip(q, inv)), simultaneous=True)))

    u = n[10]
    P = push(S.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, -3*u, 1, 0], [0, 0, 0, 1]]), S.Matrix([0, 0, 0, -3*u]))
    zero = dict.fromkeys(q, 0)
    a, b, c = P[0, 1].subs(zero), S.diff(P[0, 1], w).subs(zero), S.diff(P[0, 1], z).subs(zero)
    d, f = P[1, 2].subs(zero), S.diff(P[1, 2], w).subs(zero)
    e = S.diff(P[0, 2], w, 2)/2
    e = S.expand(e)
    expected = [a+b*w+c*(z-y*w)+(2*x-y*y)/9,
                a*w+d*y+e*w*w+f*y*w+c*z*w+(2*x*w-4*y*z)/9,
                -y*w/3, d+f*w+c*w*w-4*z/9+y*w/9, -w/3, -w*w/3]
    assert all(S.expand(P[i, j]-v) == 0 for (i,j),v in zip(itertools.combinations(range(4),2),expected))
    assert S.expand(e-b-(n[3]-n[6])) == 0
    A = -9*c/2
    B = 9*(A*f-b)
    shift_x = 9*(a-A*d)/2
    delta = 9*d/4
    M = S.Matrix([[1,0,A,B],[0,1,0,A],[0,0,1,-9*f],[0,0,0,1]])
    assert M.det() == 1
    P = push(M, S.Matrix([shift_x,0,-delta,0]))
    lam = 9*(e-b)
    expected = [2*x-y*y, lam*w*w+2*x*w-4*y*z,-3*y*w,-4*z+y*w,-3*w,-3*w*w]
    assert all(S.expand(9*P[i,j]-v)==0 for (i,j),v in zip(itertools.combinations(range(4),2),expected))
    pf = S.expand(P[0,1]*P[2,3]-P[0,2]*P[1,3]+P[0,3]*P[1,2])
    assert S.expand(pf-(n[3]-n[6])*w**3/3) == 0
    report = {
        'status':'PASS', 'scope':'Complete seven-parameter A=-1/3 source, no new componenthood claim',
        'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (INPUT,DATA,Path(__file__))},
        'first_normalized_parameters':dict(zip(['a','b','c','d','e','f'],map(str,[a,b,c,d,e,f]))),
        'second_linear_matrix':[[str(v) for v in row] for row in M.tolist()],
        'second_translation':list(map(str,[shift_x,0,-delta,0])),
        'canonical_six_brackets_scaled_by_9':list(map(str,expected)),
        'lambda':str(lam), 'universal_six_coefficient_identity':True,
        'pfaffian_identity':str(pf), 'all_rank_four_source_is_one_orbit':True,
        'seconds':time.monotonic()-started,
    }
    OUT.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','seconds':report['seconds'],'report':str(OUT.relative_to(ROOT))}))


if __name__ == '__main__':
    main()

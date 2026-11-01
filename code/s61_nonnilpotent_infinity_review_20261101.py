#!/usr/bin/env python3
"""Independent coefficient/zero-scheme audit of the nonnilpotent exclusion."""
from __future__ import annotations
import argparse
import hashlib
import json
from itertools import combinations, product
from pathlib import Path
import sympy as s
from s61_model import S61, X, J as J5

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'reports/S61_NONNILPOTENT_INFINITY_EXCLUSION_20261101.json'


def run():
    data = json.loads(SOURCE.read_text())
    model = S61()
    d, e, a = s.symbols('d e a')
    loc = {str(v): v for v in (*X, d, e, a)}
    matrix = lambda rows: s.Matrix([[s.sympify(v, locals=loc) for v in row] for row in rows])
    weights = model.reduced['parameter_weights']
    first = [i for i, w in enumerate(weights) if w == 1]
    second = [i for i, w in enumerate(weights) if w == 2]
    beta = s.symbols('b0:16')
    scale = s.Symbol('h')

    def rebuild(entry):
        Q = matrix(entry['Q'])
        biv = X[0]**2 * J5
        biv[1:5, 1:5] += Q
        q = model.normalize(biv)['parameters']
        assert all(q[i] == 0 for i in first)
        values = dict(zip(first, beta)) | {i: scale*q[i] for i in second}
        rows = []
        for eq in model.reduced['equations']:
            # Truncate the substitution in h before forming a large symbolic
            # sum. This is independent of the main leading-matrix helper.
            terms = [s.Rational(n, den)*s.prod(values[i] for i in mon)
                     for mon, n, den in eq if sum(i in second for i in mon) == 2]
            f = s.expand(s.Add(*terms))
            lead = f.coeff(scale, 2)
            assert lead.subs(dict.fromkeys(beta, 0)) == 0
            row = [lead.coeff(b) for b in beta]
            assert s.expand(lead - sum(v*b for v, b in zip(row, beta))) == 0
            rows.append(row)
        return Q, s.Matrix(rows)

    def verify_blocks(M, blocks):
        occupied = set()
        for block in blocks:
            rows, cols = block['rows'], block['columns']
            B = matrix(block['matrix']) if rows else s.zeros(0, len(cols))
            assert (M.extract(rows, cols) - B).applyfunc(s.expand) == s.zeros(*B.shape)
            for i in rows:
                for j in cols:
                    occupied.add((i, j))
            assert all(s.expand(B[i,j]*B[k,l]-B[i,l]*B[k,j]) == 0
                       for i,k in combinations(range(B.rows),2)
                       for j,l in combinations(range(B.cols),2))
        assert all(M[i,j] == 0 or (i,j) in occupied
                   for i in range(M.rows) for j in range(M.cols))

    logs = []
    for family in data['log_families']:
        Q, M = rebuild(family)
        verify_blocks(M, family['blocks'])
        points = set()
        audits = []
        for block in family['blocks']:
            B = matrix(block['matrix'])
            gb = s.groebner([v for v in B if v != 0], d, e, domain=s.QQ)
            claimed = {tuple(map(s.Rational, p)) for p in block['zero_points']}
            if gb.contains(s.S.One):
                assert not claimed
                degree = 0
            else:
                assert gb.is_zero_dimensional
                leading = [p.LM(order=gb.order).exponents for p in gb.polys]
                # A pure power in each variable bounds the standard monomials.
                bd = min(i for i,j in leading if j == 0)
                be = min(j for i,j in leading if i == 0)
                degree = sum(not any(i >= u and j >= v for u,v in leading)
                             for i,j in product(range(bd), range(be)))
                assert degree == len(claimed)
                assert all(all(p.as_expr().subs({d:pd,e:pe}) == 0 for p in gb.polys)
                           for pd,pe in claimed)
            # Length equals the number of certified distinct rational points:
            # there are no missing complex points or nonreduced zero blocks.
            points |= claimed
            audits.append({'groebner_basis':list(map(str,gb.polys)), 'length':degree})
        ranks = {(s.Rational(v['d']),s.Rational(v['e'])):v['rank']
                 for v in family['rank_drop_points']}
        assert points == set(ranks)
        assert all(M.subs({d:p[0],e:p[1]}).rank() == rank for p,rank in ranks.items())
        logs.append({'label':family['label'], 'independent_zero_scheme_audit':audits,
                     'exceptional_count':len(points), 'minimum_rank':min(ranks.values())})

    boundary = data['G2_A_zero_boundary']
    Q, M = rebuild(boundary)
    verify_blocks(M, boundary['blocks'])
    roots = set()
    nonzero_blocks = 0
    for block in boundary['blocks']:
        vals = [v for v in matrix(block['matrix']) if v != 0]
        if vals:
            nonzero_blocks += 1
            gb = s.groebner(vals, e, domain=s.QQ)
            g = gb.polys[0].as_expr()
            rr = s.roots(g, e)
            assert sum(rr.values()) == s.degree(g,e)
            roots.update(rr)
    assert nonzero_blocks == 12
    assert roots == {-s.S.One,-s.Rational(1,2),s.S.Zero,s.Rational(1,2)}
    assert all(M.subs(e,r).rank() == 11 for r in roots)

    fixed = []
    for item in data['fixed_representatives']:
        Q, M = rebuild(item)
        assert (M - matrix(item['beta_matrix'])).applyfunc(s.expand) == s.zeros(*M.shape)
        w = item['minor']
        det = M.extract(w['rows_zero_based'],w['columns_zero_based']).det()
        assert det != 0 and det == s.Rational(w['determinant'])
        fixed.append({'label':item['label'],'rank_lower_bound':w['rank'],'determinant':str(det)})
    mixed = []
    for item in data['mixed_lines']:
        Q, M = rebuild(item)
        assert (M - matrix(item['beta_matrix'])).applyfunc(s.expand) == s.zeros(*M.shape)
        det = s.factor(M.extract(item['minor_rows'],item['minor_columns']).det(method='domain-ge'))
        assert det != 0 and not det.has(a) and det == s.sympify(item['constant_minor'])
        mixed.append({'label':item['label'],'uniform_rank_lower_bound':len(item['minor_columns']),
                      'constant_minor':str(det)})
    paths = [Path(__file__), SOURCE, ROOT/'code/s61_nonnilpotent_infinity_exclusion_20261101.py',
             ROOT/'research/notes/S61_NONNILPOTENT_INFINITY_EXCLUSION_20261101.md',
             ROOT/'legacy/v03/proofs/regular_semisimple_theorem.txt',
             ROOT/'legacy/v04/semisimple_boundary.tex',
             ROOT/'legacy/v05/proofs/mixed_modular_theorem.txt']
    return {'id':'S61_NONNILPOTENT_INFINITY_REVIEW_20261101','all_passed':True,
            'field':'QQ','log_families':logs,
            'corrected_A_zero_boundary':{'generic_rank':12,'exceptional_parameters':list(map(str,sorted(roots))),
                                         'exceptional_ranks':[11]*4},
            'fixed_representatives':fixed,'mixed_lines':mixed,
            'scope':'Independent extraction of all h^2 coefficients and exact Groebner zero-scheme audits; historical coverage checked by reading the archived proofs, not recomputed.',
            'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open('x') as f:
        json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'output':str(args.output)}))

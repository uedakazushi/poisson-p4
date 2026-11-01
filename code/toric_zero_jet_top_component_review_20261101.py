#!/usr/bin/env python3
"""Independent action-matrix and divisor audit of the unique top zero-jet component."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
from pathlib import Path
import time
import sympy as s
from s61_model import S61, X, quotient_coefficients

ROOT = Path(__file__).resolve().parents[1]


def verify():
    start = time.monotonic(); h = s.symbols('h0:9')
    prior_path = ROOT/'reports/TORIC_ZERO_JET_DIMENSION_REVIEW_20261101.json'
    top_path = ROOT/'reports/TORIC_ZERO_JET_TOP_COMPONENT_20261101.json'
    prior = json.loads(prior_path.read_text()); top = json.loads(top_path.read_text())
    def schur(M, rows, columns):
        pivot = M.extract(rows, columns)
        assert pivot.det() != 0 and not pivot.det().free_symbols
        rest = [i for i in range(M.cols) if i not in columns]
        return (M[:, rest]-M[:, columns]*pivot.inv()*M[rows, rest]).applyfunc(s.factor)
    R = s.Matrix(prior['r0_action_matrix'])
    double = R.subs(dict(zip(h[5:], [0, 1, 0, 0], strict=True)))
    residual = schur(double, [2, 3, 5, 6, 7], [0, 2, 3, 4, 5])
    GB = s.groebner(list(residual), *h[:5])
    assert list(GB) == [h[0]-h[2]*h[3], h[1]-h[3]**2/4, h[4]]
    assert R.subs(dict(zip(h, [1, 2, 3, 4, 5, 1, 0, 0, 0], strict=True))).rank() == 7
    minimal = {}
    for kind in ['M0', 'M1']:
        M = s.Matrix(prior['minimal_models'][kind]['residual_potential_action_matrix'])
        M = M.subs(dict(zip(h[5:], [1, 0, 0, 0], strict=True)))
        res = schur(M, [2, 5, 6], [0, 1, 4])
        res = res[[i for i in range(res.rows) if any(res[i, :])], :]
        assert res.shape == (4, 3) and s.expand(res.extract([2, 3], [0, 1]).det()-h[4]**2) == 0
        minors = [s.expand(res.extract(rr, cc).det()) for rr in combinations(range(4), 2)
                  for cc in combinations(range(3), 2)]
        gb = s.groebner(minors+[h[4]], *h[:5])
        supports = [{i for i, a in enumerate(p.monoms()[0]) if a} for p in gb.polys]
        dimension = max(len(values) for n in range(6) for values in combinations(range(5), n)
                        if not any(T <= set(values) for T in supports))
        assert dimension == 2
        minimal[kind] = {'alternative_schur_matrix': [[str(v) for v in row] for row in res.tolist()],
            'rank1_minors_plus_h4_groebner_basis': list(map(str, gb)), 'dimension': dimension}
    hist_path = ROOT/'legacy/v07/data/nil22_classification_summary.json'
    history = json.loads(hist_path.read_text()); x, y, z, w = s.symbols('x y z w')
    pfchecks = []
    for recorded in top['nil22_factorization_tests']:
        name = recorded['stratum']; entry = history[name]; q = {}
        for key, value in entry['Q'].items():
            pair = tuple(map(int, key.strip('()').split(','))); q[pair] = s.sympify(value)
        assert s.expand(q[0, 1]+q[2, 3]-x*z) == 0
        pf = s.expand(q[0, 1]*q[2, 3]-q[0, 2]*q[1, 3]+q[0, 3]*q[1, 2])
        assert s.expand(pf-s.sympify(recorded['pfaffian'])) == 0
        passes = False
        if pf and not pf.has(y, w):
            coefficient, factors = s.factor_list(pf, x, z, extension=s.sqrt(3))
            if all(e % 2 == 0 for _, e in factors):
                root = s.sqrt(coefficient)*s.prod(f**(e//2) for f, e in factors)
                determinants = []
                for sign in [1, -1]:
                    p = s.Poly(s.expand(x*z+2*sign*root), x, z)
                    determinants.append(s.simplify(4*p.coeff_monomial(x*x)*p.coeff_monomial(z*z)-p.coeff_monomial(x*z)**2))
                passes = all(v == 0 for v in determinants)
        assert passes == recorded['passes_necessary_complete_factorization'] == (name == 'AB')
        pfchecks.append({'stratum': name, 'passes_necessary_factorization': passes})
    QM = [3*x*x-6*x*z, 0, -2*x*x, 3*x*z-18*z*z, 6*x*x-2*x*y-3*x*w+30*z*w, 6*x*z]
    assert s.expand(QM[0]+QM[5]) == 3*x*x
    assert s.expand(QM[0]*QM[5]-QM[1]*QM[4]+QM[2]*QM[3]) == 12*x**3*z
    x0, x, y, z, w = X
    U = s.Matrix([0, z, 0, 3*x]); V = s.Matrix([-3*x, 5*y, z, w])
    P = s.zeros(5); P[1:, 1:] = V*U.T-U*V.T
    P[1, 2] += x0*x0; P[2, 1] -= x0*x0; P[3, 4] += x0*x0; P[4, 3] -= x0*x0
    model = S61(); model.verify_point(model.normalize(P))
    columns = []
    for i in range(5):
        for j in range(5):
            A = s.zeros(5); A[i, j] = 1
            columns.append(s.Matrix(quotient_coefficients(A*P+P*A.T-X[j]*P.diff(X[i]))))
    columns.append(s.Matrix(quotient_coefficients(P)))
    assert s.Matrix.hstack(*columns).rank() == 22
    paths = [Path(__file__), prior_path, top_path, hist_path, ROOT/'code/s61_model.py',
             ROOT/'legacy/v07/minimal_thickness.tex']
    return {'id': 'TORIC_ZERO_JET_TOP_COMPONENT_REVIEW_20261101', 'all_passed': True,
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'R0_double_alternative_groebner_basis': list(map(str, GB)), 'R0_triple_rank': 7,
        'minimal_triple_alternative_schur_checks': minimal, 'nil22_factorization_checks': pfchecks,
        'third_minimal_representative_trace': '3*x^2', 'third_minimal_representative_Pfaffian': '12*x^3*z',
        'zero_modular_L0_projective_stabilizer_dimension': 3,
        'conclusion': 'AB is the unique dimension10 irreducible component of the reduced actual toric zero-first-jet fiber.',
        'scope': 'Other components of dimension at most9 remain open; historical Q45 exhaustion and AB inclusion remain dependencies.',
        'proof_note': 'research/notes/TORIC_ZERO_JET_TOP_COMPONENT_REVIEW_20261101.md',
        'elapsed_seconds': round(time.monotonic()-start, 3),
        'input_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); report = verify(); args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out: json.dump(report, out, indent=2); out.write('\n')
    print(json.dumps({'all_passed': True, 'seconds': report['elapsed_seconds']}))


if __name__ == '__main__': main()

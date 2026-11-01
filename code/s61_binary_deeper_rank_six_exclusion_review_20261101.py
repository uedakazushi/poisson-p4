#!/usr/bin/env python3
"""Independent divisor/minor and group-weight review of the B-rank<=6 exclusion."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def verify():
    started = time.time()
    primary_path = ROOT/'reports/S61_BINARY_DEEPER_RANK_SIX_EXCLUSION_20261101.json'
    cover_path = ROOT/'reports/S61_BINARY_RANK_DROP_COVER_20261101.json'
    action_path = ROOT/'reports/S61_BINARY_UNIFORM_RANK_20261101.json'
    prior_path = ROOT/'reports/S61_BINARY_DEEPER_REDUCTION_20261101.json'
    primary = json.loads(primary_path.read_text())
    cover = json.loads(cover_path.read_text())
    action = json.loads(action_path.read_text())
    prior = json.loads(prior_path.read_text())
    for report in [primary, cover, action, prior]:
        assert report['all_passed']
    for path, digest in primary['input_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
    u = s.symbols('u1:62')
    a,b,c,d,e,f = [u[i] for i in [3,13,34,35,48,49]]
    core = [a,b,c,d,e,f]
    B = s.Matrix(cover['block_matrices'][1]).applyfunc(s.sympify)
    K,L,Q = 5*c+f, 5*b-e, 600*b*b-255*b*e-63*e*e+100*a*f
    # Q is absolutely irreducible: primitive and linear in a over k[b,d,e,f].
    assert s.gcd(100*f, Q.subs(a,0)) == 1
    chart = B.subs({c:-f/5, a:-Q.subs(a,0)/(100*f)}, simultaneous=True)
    rows, cols = primary['generic_minor_rows'], primary['generic_minor_columns']
    # Clear f first and use a polynomial-domain determinant, independently of
    # the producer's rational-function determinant.
    matrix = (f*chart[rows,cols]).applyfunc(s.cancel)
    assert all(s.denom(x).free_symbols == set() for x in matrix)
    dm = matrix.to_DM()
    determinant = s.cancel(dm.domain.to_sympy(dm.det())/f**7)
    assert s.expand(determinant-s.sympify(primary['generic_seven_minor'])) == 0
    A = 100*b*b-505*b*e+187*e*e
    R = d*A-40*f*f*(2*b-e)
    assert s.gcd(A, 40*f*f*(2*b-e)) == 1
    assert s.degree(R,d) == 1
    assert s.expand(Q.subs(f,0)-3*(5*b-3*e)*(40*b+7*e)) == 0
    factor_product = (5*b-11*e)*(5*b-3*e)*L**2*R
    assert s.cancel(determinant/factor_product) == s.Rational(27,12500000)
    labels = set()
    for case in primary['nonzero_minor_on_each_remaining_divisor']:
        sample = {s.Symbol(k):s.Rational(v) for k,v in case['sample'].items()}
        assert K.subs(sample) == Q.subs(sample) == 0 and L.subs(sample) != 0
        label = case['divisor']; labels.add(label)
        equation = {'5b-11e':5*b-11*e, '5b-3e':5*b-3*e, 'R':R,
                    'f=0,5b-3e=0':5*b-3*e, 'f=0,40b+7e=0':40*b+7*e}[label]
        assert equation.subs(sample) == 0
        if label.startswith('f=0'):
            assert sample[f] == 0
        block = B.subs(sample)[case['rows'],case['columns']]
        value = block.det(method='bareiss')
        assert value == s.Rational(case['minor']) and value != 0
    assert len(labels) == 5
    # Reconstruct the infinitesimal core action from the already reviewed
    # primitive representation, rather than trusting the new report's weights.
    one = [3,4,11,12,13,24,25,26,33,34,35,36,47,48,49,50]
    M = s.Matrix.hstack(*[s.Matrix(v).applyfunc(s.sympify) for v in action['primitive_actions']])
    M = M.subs({u[11]:0,u[12]:0})
    C = M[[one.index(i) for i in [3,13,34,35,48,49]],:]
    assert C[:,[0,1,3]] == s.zeros(6,3)
    assert M[[one.index(11),one.index(12)],:] == s.zeros(2,5)
    for j, multiplier, record in zip([2,4],[1,2],primary['stabilizer_core_torus_weights']):
        weights = record['core_weights']
        assert multiplier*C[:,j] == s.Matrix([w*x for w,x in zip(weights,core)])
        rw = list(map(s.Rational,record['row_weights']))
        cw = list(map(s.Rational,record['column_weights']))
        for i in range(12):
            for k in range(9):
                for mon, coefficient in s.Poly(B[i,k],*core).terms():
                    if coefficient:
                        assert sum(x*w for x,w in zip(mon,weights)) == rw[i]+cw[k]
        for polynomial in [K,L,Q]:
            monomial_weights = {sum(x*w for x,w in zip(mon,weights))
                                for mon,coefficient in s.Poly(polynomial,*core).terms() if coefficient}
            assert len(monomial_weights) == 1
    H = s.Matrix.hstack(*[s.Matrix(h).applyfunc(s.sympify).reshape(16,1)
                         for h in action['stabilizer_generators']])
    assert H.rank() == 5
    assert {r['ideal_index'] for r in prior['excluded_components']} == {0,2,3}
    assert all(r['beta_image_dimension_upper_bound'] <= 15 for r in prior['excluded_components'])
    paths = [Path(__file__), primary_path, cover_path, action_path, prior_path,
             ROOT/'research/notes/S61_BINARY_DEEPER_RANK_SIX_EXCLUSION_REVIEW_20261101.md']
    return dict(id='S61_BINARY_DEEPER_RANK_SIX_EXCLUSION_REVIEW_20261101',
                all_passed=True, accepted=True, field='QQ',
                independent_checks=['Polynomial determinant after clearing f',
                    'Absolute irreducibility by primitive linear polynomials',
                    'All five remaining divisors have a nonzero seven-minor',
                    'Primitive stabilizer actions and all entry weights',
                    'Invariant core dimension at most two and image dimension at most fifteen'],
                beta_image_dimension_upper_bound=15,
                limitations=['Reuses the previously reviewed B block and primitive representation.',
                    'No claim about the remaining B-rank seven locus.'],
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                seconds=time.time()-started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args(); report = verify()
    with args.output.open('x') as out:
        json.dump(report,out,indent=2); out.write('\n')
    print(json.dumps({k:report[k] for k in ['all_passed','beta_image_dimension_upper_bound','seconds']}))

#!/usr/bin/env python3
"""The common missing core curve of the two binary frames has proper beta image."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify():
    started = time.time()
    action_path = ROOT / 'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json'
    old_path = ROOT / 'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json'
    new_path = ROOT / 'reports/S61_BINARY_HYPERPLANE_BFREE_FRAME_20261102.json'
    action = json.loads(action_path.read_text())
    old = json.loads(old_path.read_text())
    new = json.loads(new_path.read_text())
    assert action['all_passed'] and old['all_passed'] and new['all_passed']
    a, e, d, p, q, r, h, k = s.symbols('a e d p q r h k')
    S = 3*e - 2
    F = 500*a + 90*e**2 + 3*e - 24
    B = 50*a*d + 18*d*e**2 - 9*d*e - 2*d + 60*e - 40
    C = 100*a + 18*e**2 - 9*e - 2
    D0 = 50*a + 9*e**2 + 3*e - 6
    av = -3*S*(e+1)/50
    dv = -20/S
    assert s.expand(D0.subs(a, av)) == 0
    assert s.expand(B.subs(a, av) - S*(S*d+20)) == 0
    assert s.expand(C.subs(a, av) + 5*S) == 0
    assert s.expand(F.subs(a, av) - 9*(2-S)) == 0
    assert s.cancel(B.subs({a: av, d: dv}, simultaneous=True)) == 0
    one = [3,4,11,12,13,24,25,26,33,34,35,36,47,48,49,50]
    u = s.symbols('u1:62')
    uu = [u[i] for i in one]
    beta = {v: s.S.Zero for v in uu}
    beta.update({u[3]: av, u[13]: (e+1)/5, u[34]: s.S.One,
                 u[35]: dv, u[48]: e, u[49]: -5,
                 u[24]: p, u[25]: q, u[26]: r, u[36]: h, u[47]: k})
    variables = [e,p,q,r,h,k]
    bv = s.Matrix([beta[v] for v in uu])
    A = s.Matrix(action['primitive_action']).applyfunc(s.sympify)
    assert A.shape == (16,10)
    M = A.subs(beta, simultaneous=True).row_join(bv.jacobian(variables))
    nuisance_rows = [one.index(i) for i in [24,25,26,36,47]]
    assert M[:,11:] == s.eye(16)[:,nuisance_rows]
    keep = [i for i in range(16) if i not in nuisance_rows]
    scale = 50*S**2
    reduced = (M[keep,:11]*scale).applyfunc(lambda f: s.cancel(f).expand())
    domain = s.QQ.poly_ring(*variables)
    determinant = DomainMatrix.from_Matrix(reduced).convert_to(domain).det()
    assert not determinant
    sample = dict(zip(variables, [2,1,2,3,4,5]))
    N = M.subs(sample)
    _, columns = N.to_DM().rref()
    _, rows = N[:,list(columns)].T.to_DM().rref()
    assert len(columns) == len(rows) == 15
    minor = N[list(rows),list(columns)].det()
    assert minor
    point_core = {a: av.subs(sample), e: 2, d: dv.subs(sample)}
    assert all(s.cancel(f.subs(point_core)) for f in [a,F,C])
    note = ROOT / 'research/notes/S61_BINARY_HYPERPLANE_FRAME_COVER_20261102.md'
    paths = [Path(__file__), action_path, old_path, new_path, note]
    hashes = {str(path.relative_to(ROOT)): digest(path) for path in paths}
    for report in [action,old,new]:
        for name, expected in report['input_sha256'].items():
            path = ROOT / name
            assert digest(path) == expected, name
            assert name not in hashes or hashes[name] == expected
            hashes[name] = expected
    return dict(
        id='S61_BINARY_HYPERPLANE_FRAME_COVER_20261102', all_passed=True,
        core_polynomials=dict(F=str(F),B=str(B),C=str(C),D0=str(D0),S=str(S)),
        common_missing_curve=dict(a=str(av),d=str(dv)),
        curve_identities=dict(B='S*(S*d+20)',C='-5*S',F='9*(2-S)'),
        first_jet_parameter_order=[str(v) for v in uu],
        source_parameter_order=list(map(str,variables)),
        source_beta={str(v): str(beta[v]) for v in uu},
        nuisance_identity_rows=nuisance_rows, reduced_rows=keep,
        denominator_clearing_scale=str(scale),
        polynomial_reduced_matrix=[[str(v) for v in row] for row in reduced.tolist()],
        polynomial_determinant_identically_zero=True,
        source_differential_rank_upper=15,
        rank_witness=dict(parameter_values=[2,1,2,3,4,5],
                          rows=list(rows),columns=list(columns),determinant=str(minor)),
        saturated_source_beta_image_dimension=15,
        sufficiently_general_beta_avoids_common_missing_core_curve=True,
        scope='The specified normalized H8 source on aFC!=0. The old B and new D0 core-frame opens cover every general-beta source point. The frames themselves are nuisance-zero models; no relative extension or nilpotence follows here.',
        limitations=['No quartic surjectivity or length65 conclusion.',
                     'No full hyperplane or full S61 exhaustion.'],
        proof_note=str(note.relative_to(ROOT)),input_sha256=hashes,
        seconds=time.time()-started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    result = verify()
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'],
                          rank_witness=result['rank_witness'])))

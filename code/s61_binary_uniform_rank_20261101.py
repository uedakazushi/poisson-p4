#!/usr/bin/env python3
"""Exact generic binary infinity ranks and a dominant nine-dimensional slice."""
import argparse
import hashlib
import json
import time
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import J, COORDS, sp_basis, infinitesimal_action


def verify():
    start = time.time()
    root = Path(__file__).resolve().parents[1]
    model = S61()
    u = s.symbols('u1:62')
    weights = model.reduced['parameter_weights']
    one = [i for i, w in enumerate(weights) if w == 1]
    two = [i for i, w in enumerate(weights) if w == 2]
    index = {i: j for j, i in enumerate(two)}
    beta = {i: (s.S.Zero if i in (11, 12) else u[i]) for i in one}
    extra_path = root / 'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    extra = json.loads(extra_path.read_text())['universal_new_equation']
    equations = model.reduced['equations'] + [extra]
    fixed = index[29]
    q0 = s.zeros(45, 1)
    q0[fixed] = 1  # Q=xw d_y wedge d_z.
    D = s.zeros(135, 45)
    tc = s.zeros(135, 1)
    values = []
    for row, eq in enumerate(equations):
        value = 0
        for mon, n, d in eq:
            qm = [index[i] for i in mon if i in index]
            c = s.Rational(n, d) * s.prod(beta[i] for i in mon if i not in index)
            if len(qm) == 2:
                a, b = qm
                value += c * q0[a] * q0[b]
                D[row, a] += c * q0[b]
                D[row, b] += c * q0[a]
            elif len(qm) == 1:
                tc[row, 0] += c * q0[qm[0]]
        values.append(s.expand(value))
    assert all(v == 0 for v in values)
    K = s.Matrix.hstack(*D[:70, :].nullspace())
    assert K.shape == (45, 29)
    residual = (D[70:, :] * K).applyfunc(s.expand)
    rr, piv = residual.to_DM().rref()
    rr = rr.to_Matrix()
    free = [i for i in range(29) if i not in piv]
    assert len(piv) == 18 and len(free) == 11
    N = s.zeros(29, 11)
    for j, col in enumerate(free):
        N[col, j] = 1
        for row, pivot in enumerate(piv):
            N[pivot, j] = -rr[row, col]
    assert (residual * N).applyfunc(s.cancel) == s.zeros(65, 11)
    assert N[free, :] == s.eye(11)
    sample = {u[i]: s.Integer(j + 2) for j, i in enumerate(one)}
    sample[u[11]] = sample[u[12]] = 0
    Ds = D.subs(sample)
    leading_columns = list(Ds.to_DM().rref()[1])
    leading_rows = list(Ds[:, leading_columns].T.to_DM().rref()[1])
    assert len(leading_columns) == 34
    leading_det = Ds.extract(leading_rows, leading_columns).det(method='domain-ge')
    chart = [i for i in range(45) if i != fixed]
    full = D[:, chart].row_join(tc)
    fulls = full.subs(sample)
    full_columns = list(fulls.to_DM().rref()[1])
    full_rows = list(fulls[:, full_columns].T.to_DM().rref()[1])
    assert len(full_columns) == 35
    full_det = fulls.extract(full_rows, full_columns).det(method='domain-ge')
    assert leading_det != 0 and full_det != 0
    z, primitive_u, primitive = primitive_beta(model)
    Q = s.zeros(4)
    Q[1, 2] = z[0] * z[3]
    Q[2, 1] = -Q[1, 2]
    monomials = [z[i] * z[j] for i, j in combinations_with_replacement(range(4), 2)]
    def vector(M):
        return s.Matrix([s.Poly(s.expand(M[i, j]), *z).coeff_monomial(mon)
            for i, j in combinations(range(4), 2) for mon in monomials])
    sp = sp_basis()
    action_q = [vector(A*Q + Q*A.T - sum(((A*s.Matrix(z))[i]*Q.diff(z[i])
        for i in range(4)), s.zeros(4))) for A in sp]
    stabilizer = s.Matrix.hstack(*action_q, -vector(Q)).nullspace()
    assert len(stabilizer) == 5
    generators = [sum((v[i]*sp[i] for i in range(10)), s.zeros(4)) for v in stabilizer]
    tensor = s.Matrix([primitive[i, j].diff(z[a]) for (i, j), a in COORDS])
    basis = tensor.jacobian(primitive_u)
    basis_rows = list(basis.T.rref()[1])
    left = basis[basis_rows, :].inv()
    actions = [(left*infinitesimal_action(tensor, A)[basis_rows, :]).applyfunc(s.expand)
        for A in generators]
    constraints = {u[11]: 0, u[12]: 0}
    for action in actions:
        assert all(action[primitive_u.index(v)].subs(constraints) == 0 for v in constraints)
    fixed_beta = {u[13]: s.S.One, u[34]: s.S.One, u[4]: s.S.Zero,
                  u[24]: s.S.Zero, u[25]: s.S.Zero}
    slice_sample = sample | fixed_beta
    selected_rows = [primitive_u.index(v) for v in fixed_beta]
    tangent = s.Matrix.hstack(*actions)
    slice_det = tangent[selected_rows, :].subs(slice_sample).det()
    assert slice_det == 553
    def enc(M):
        return [[str(v) for v in row] for row in M.tolist()]
    paths = [Path(__file__), root/'code/s61_model.py',
        root/'code/primitive_jet_elliptic_curve_20261101.py', root/'code/weight1_structure.py', extra_path]
    return dict(id='S61_BINARY_UNIFORM_RANK_20261101', all_passed=True,
        field='QQ(u_i : weight1, u12=u13=0)', normal_form='Q=xw*d_y wedge d_z',
        all135_leading_equations_zero=True, constant_first70_rank=16,
        residual_rank=18, generic_leading_jacobian_rank=34,
        generic_full_homogenized_jacobian_rank=35,
        generic_full_projective_embedding_dimension=10,
        constant_kernel=enc(K), residual_kernel=enc(N), residual_free_columns=free,
        leading_minor_rows=leading_rows, leading_minor_columns=leading_columns,
        leading_minor_at_sample=str(leading_det), full_minor_rows=full_rows,
        full_minor_columns=full_columns, full_minor_at_sample=str(full_det),
        sample={str(k):str(v) for k,v in sample.items()}, stabilizer_generators=[enc(A) for A in generators],
        primitive_actions=[enc(A) for A in actions], slice_fixed_coordinates={str(k):str(v) for k,v in fixed_beta.items()},
        slice_tangent_minor=str(slice_det), slice_dimension=9,
        proof_note='research/notes/S61_BINARY_UNIFORM_RANK_20261101.md',
        limitations=['Generic ranks only; exceptional parameter loci are not exhausted.',
            'Does not prove universal quadratic rank5 or transverse multiplicity65.',
            'No conclusion about t-saturation or affine fiber exhaustion.'],
        source_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.time()-start)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=verify()
    with args.output.open('x') as f:
        json.dump(result,f,indent=2)
        f.write('\n')
    print(json.dumps({k:result[k] for k in ['all_passed','generic_leading_jacobian_rank',
        'generic_full_homogenized_jacobian_rank','slice_tangent_minor','seconds']}))


if __name__=='__main__':
    main()

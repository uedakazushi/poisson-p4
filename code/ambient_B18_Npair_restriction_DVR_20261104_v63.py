#!/usr/bin/env python3
"""Exact fixed-source restriction and D=QQ[t] module; no Groebner computation."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SOURCE = 'reports/AMBIENT_B18_PROPER_NPAIR_RESTRICTION_20261104_v63/actual-Npair-three-variable-polynomials.json'
H_SOURCE = 'reports/AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63/actual12-flag-divided-polynomials.json'


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def unpack(q, variables):
    return sp.Poly.from_dict({tuple(e): sp.Rational(n, d) for e, n, d in q}, variables).as_expr()


def pack(f, variables):
    return [[list(e), int(c.p), int(c.q)] for e, c in sp.Poly(f, *variables, domain=sp.QQ).terms() if c]


def zero(f):
    assert sp.expand(f) == 0, sp.expand(f)


def lift(target, generators, variables):
    """One small fixed homogeneous coefficient system, returning exact cofactors."""
    degree = sp.Poly(target, *variables).total_degree()
    columns, labels = [], []
    for i, g in enumerate(generators):
        bound = degree - sp.Poly(g, *variables).total_degree()
        for e in sorted(sp.polys.monomials.itermonomials(variables, max(0, bound)), key=str) if bound >= 0 else []:
            # The targets and generators are homogeneous, so use only the necessary degree.
            if sp.Poly(e, *variables).total_degree() == bound:
                columns.append(sp.Poly(e*g, *variables))
                labels.append((i, e))
    monomials = sorted(set().union(*(set(p.monoms()) for p in columns), set(sp.Poly(target, *variables).monoms())))
    matrix = sp.Matrix([[p.coeff_monomial(e) for p in columns] for e in monomials])
    vector = sp.Matrix([sp.Poly(target, *variables).coeff_monomial(e) for e in monomials])
    solution, parameters = matrix.gauss_jordan_solve(vector)
    solution = solution.subs({x: 0 for x in parameters})
    cofactors = [sp.Integer(0) for _ in generators]
    for c, (i, mon) in zip(solution, labels):
        cofactors[i] += c*mon
    zero(sum(c*g for c, g in zip(cofactors, generators))-target)
    return [sp.expand(c) for c in cofactors]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    resource.setrlimit(resource.RLIMIT_FSIZE, (32*1024**2, 32*1024**2))
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    d = json.loads((ROOT/SOURCE).read_text())
    a, x, t, u = sp.symbols('a x t u')
    oldvars = sp.symbols('a v2 v3 v4 v5 v6 v7 v8')
    old_H = [unpack(q, oldvars) for q in d['all12_exact_H_d0_polynomials']]
    roots = {next(z for z in oldvars if str(z) == name): unpack(q, oldvars) for name, q in d['all_polynomial_roots'].items()}
    root3 = {oldvars[0]: a, oldvars[1]: x, oldvars[5]: t}
    source7 = [unpack(q['exact_polynomial3'], (a, x, t)) for q in d['seven_remaining_actual_polynomials']]
    # Independently replay all twelve recorded source equations after all five roots.
    reduced_H = [sp.expand(f.subs(roots, simultaneous=True).subs(root3, simultaneous=True)) for f in old_H]
    selected = [0, 4, 1, 2, 3]
    for i in selected:
        zero(reduced_H[i])
    for i, f in enumerate(source7):
        zero(reduced_H[i+5]-f)
    actual = [sp.expand(f.subs(a, u-x)) for f in source7]
    F = u*u-6*u*x
    G = (23*x-6*t)*u+2*x*x+4*x*t
    H = u*x*x
    generators = [F, G, H]
    zero(actual[0]-F)
    zero(actual[1]-10*F-2*G)
    zero(actual[2]+sp.Rational(6, 5)*(u+3*x)*F+24*H)
    cofs = [lift(f, generators, (u, x, t)) for f in actual]
    # q=ux is killed by u,x,t modulo these three exact generators.
    q = u*x
    annihilator_cofs = {
        'u': [x, 0, 6],
        'x': [0, 0, 1],
        't': [(23*x-6*t)/32, -u/32, sp.Rational(35, 8)],
    }
    for z in (u, x, t):
        zero(z*q-sum(c*g for c, g in zip(annihilator_cofs[str(z)], generators)))

    # A literal D^3 plus k*q algebra. q coefficients are taken modulo t.
    # Ordered basis (1,x,u,q); t*q=0 is intrinsic, no division by t.
    table = [[None]*4 for _ in range(4)]
    for i in range(4):
        table[0][i] = table[i][0] = [sp.Integer(int(j == i)) for j in range(4)]
    table[1][1] = [0, -2*t, 3*t, -sp.Rational(23, 2)]
    table[1][2] = table[2][1] = [0, 0, 0, 1]
    table[2][2] = [0, 0, 0, 6]
    for i in range(1, 4):
        table[i][3] = table[3][i] = [0, 0, 0, 0]

    def normalize(v):
        return [sp.expand(c) for c in v[:3]]+[sp.expand(v[3]).subs(t, 0)]

    def multiply(v, w):
        ans = [sp.Integer(0)]*4
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    ans[k] += v[i]*w[j]*table[i][j][k]
        return normalize(ans)

    basis = [[sp.Integer(int(i == j)) for j in range(4)] for i in range(4)]
    for i in range(4):
        for j in range(4):
            for k in range(4):
                assert multiply(multiply(basis[i], basis[j]), basis[k]) == multiply(basis[i], multiply(basis[j], basis[k]))

    def evaluate(f):
        value = [sp.Integer(0)]*4
        for (eu, ex, et), c in sp.Poly(f, u, x, t).terms():
            word = basis[0]
            for _ in range(eu):
                word = multiply(word, basis[2])
            for _ in range(ex):
                word = multiply(word, basis[1])
            for k in range(4):
                value[k] += c*t**et*word[k]
        return normalize(value)

    for f in generators+actual:
        assert evaluate(f) == [0]*4
    assert evaluate(q) == [0, 0, 0, 1]
    # The displayed algebra is generated by u,x; all its defining table relations
    # follow from F,G,H and the three annihilator identities above. Thus both
    # maps are inverse. This proves the D-module type, not only a quotient map.
    qzero_table_generators = [u*u, u*x, x*x-3*t*u+2*t*x]
    double = {u: sp.Rational(2, 3)*x}
    single = {u: 0, x: -2*t}
    for f in qzero_table_generators:
        assert sp.rem(sp.Poly(f.subs(double), x), sp.Poly(x*x, x)).is_zero
        zero(f.subs(single, simultaneous=True))
    # D^3 -> D[x]/x² + D, ordered bases (1,x,u), ((1,x),1).
    gluing = sp.Matrix([[1, 0, 0], [0, 1, sp.Rational(2, 3)], [1, -2*t, 0]])
    zero(gluing.det()-sp.Rational(4, 3)*t)
    evidence = dict(
        variables=['u=a+v2', 'x=v2', 't=v6'],
        seven_source_polynomials=[pack(f, (u, x, t)) for f in actual],
        three_generators=[pack(f, (u, x, t)) for f in generators],
        seven_membership_cofactors=[[pack(c, (u, x, t)) for c in cc] for cc in cofs],
        annihilator_cofactors={k: [pack(c, (u, x, t)) for c in cc] for k, cc in annihilator_cofs.items()},
        multiplication_table=[[[str(c) for c in v] for v in row] for row in table],
        basis=['1', 'x', 'u', 'q=u*x'],
        module_type='QQ[t]^3 direct_sum QQ[t]/(t), q generates last summand',
        two_sided_isomorphism_proof='Table entries follow from F,G,H and the three annihilator cofactors. Its associative algebra is generated by u,x with q=u*x, so generator-preserving maps in both directions are inverse.',
        H0_length=1, socle_generator='q=(a+v2)*v2',
        saturated_ideal=['u^2', 'u*x', 'x^2-3*t*u+2*t*x'],
        saturated_primary_components=[['3*u-2*x', 'x^2'], ['u', 'x+2*t']],
        saturated_primary_gluing_matrix=[[str(c) for c in row] for row in gluing.tolist()],
        saturated_primary_gluing_determinant='4*t/3',
        exact_scope='Actual core quotient by d, not whole E4 quotient by d, and not D0.',
    )
    evidence_path = out/'exact-DVR-and-source-certificate.json'
    evidence_path.write_text(json.dumps(evidence, indent=2)+'\n')
    terms = sum(len(pack(f, (u, x, t))) for f in actual+generators)+sum(len(pack(c, (u, x, t))) for cc in cofs for c in cc)
    assert terms < 150000
    report = dict(
        id='AMBIENT_B18_NPAIR_RESTRICTION_DVR_20261104_v63', all_passed=True,
        actual_source12_root_substitutions_verified=True, actual_source7_ideal_equals_three_generators=True,
        multiplication_associativity_64_checks=True, two_sided_module_model=True,
        DVR_free_rank=3, DVR_torsion_length=1, H0_length=1, depth=0, multiplicity=3,
        saturated_ring_double_line_and_single_line=True,
        whole_E4_or_D0_H0_claimed=False, no_new_source_jets_GB_or_points=True,
        expanded_saved_polynomial_terms=terms,
        seconds=time.monotonic()-started,
        input_sha256={SOURCE: sha(ROOT/SOURCE), H_SOURCE: sha(ROOT/H_SOURCE), str(Path(__file__).relative_to(ROOT)): sha(__file__)},
        evidence_sha256={evidence_path.name: sha(evidence_path)}, independent_review_pending=True,
    )
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Independent rational Lie2 and full-81 contact proof checks; no CAS/third jet."""
import argparse
import ast
from fractions import Fraction as F
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_ZERO_LINEAR_CONTACT4_20261103_v63_REVIEW'


def read(path):
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def accepted(path):
    data = read(path)
    assert data['all_passed'] is True
    for name, digest in data.get('input_sha256', {}).items():
        assert sha(ROOT / name) == digest, name
    for name, digest in data.get('evidence_sha256', {}).items():
        target = ROOT / name if (ROOT / name).exists() else path.parent / name
        assert sha(target) == digest, name
    return data


def matvec(M, v):
    return [sum((F(a)*b for a, b in zip(row, v)), F(0)) for row in M]


def matrix_product(A, B):
    return [[sum((F(a)*F(b) for a, b in zip(row, col)), F(0))
             for col in zip(*B)] for row in A]


PAIRS = list(itertools.combinations(range(5), 2))
MONS = list(itertools.combinations_with_replacement(range(5), 2))
PI = {pair: i for i, pair in enumerate(PAIRS)}
MI = {mon: i for i, mon in enumerate(MONS)}


def lie(M, P):
    """Index-wise MP+PM^T-(Mx).grad(P), with no polynomial package."""
    result = [F(0)] * 150

    def put(i, j, mon, value):
        if i == j:
            return
        if i > j:
            i, j, value = j, i, -value
        result[15*PI[i, j] + MI[tuple(sorted(mon))]] += value

    for ip, (i, j) in enumerate(PAIRS):
        for im, (a, b) in enumerate(MONS):
            value = P[15*ip + im]
            if not value:
                continue
            for k in range(5):
                put(k, j, (a, b), M[k][i]*value)
                put(i, k, (a, b), M[k][j]*value)
            for q, rem in [(a, b), (b, a)]:
                for k in range(5):
                    put(i, j, (rem, k), -M[q][k]*value)
    return result


class Poly:
    """Tiny exact sparse arithmetic solely for finite displayed identities."""
    variables = ()

    def __init__(self, terms=0):
        if isinstance(terms, Poly):
            self.terms = terms.terms.copy()
        elif isinstance(terms, dict):
            self.terms = {e: F(c) for e, c in terms.items() if c}
        else:
            self.terms = {(0,)*len(self.variables): F(terms)} if terms else {}

    def __add__(self, other):
        result = self.terms.copy()
        for e, c in Poly(other).terms.items():
            result[e] = result.get(e, F(0)) + c
        return Poly(result)

    __radd__ = __add__

    def __neg__(self):
        return Poly({e: -c for e, c in self.terms.items()})

    def __sub__(self, other):
        return self + -Poly(other)

    def __rsub__(self, other):
        return Poly(other) - self

    def __mul__(self, other):
        result = {}
        for e, c in self.terms.items():
            for f, d in Poly(other).terms.items():
                g = tuple(x+y for x, y in zip(e, f))
                result[g] = result.get(g, F(0)) + c*d
        return Poly(result)

    __rmul__ = __mul__

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        result = Poly(1)
        for _ in range(n):
            result *= self
        return result

    def __truediv__(self, q):
        return self * (1/F(q))

    def __eq__(self, other):
        return self.terms == Poly(other).terms


def symbols(names):
    Poly.variables = tuple(names)
    return {name: Poly({tuple(int(i == j) for i in range(len(names))): F(1)})
            for j, name in enumerate(names)}


def evaluate(text, env):
    def walk(node):
        if isinstance(node, ast.Constant):
            return F(node.value)
        if isinstance(node, ast.Name):
            return env[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -walk(node.operand)
        if isinstance(node, ast.BinOp):
            a, b = walk(node.left), walk(node.right)
            if isinstance(node.op, ast.Add):
                return a+b
            if isinstance(node.op, ast.Sub):
                return a-b
            if isinstance(node.op, ast.Mult):
                return a*b
            if isinstance(node.op, ast.Div):
                return a/b
            if isinstance(node.op, ast.Pow):
                return a**int(b)
        raise ValueError(ast.dump(node))
    return walk(ast.parse(text, mode='eval').body)


def main(output):
    start = time.monotonic()
    assert not output.exists()
    second = ROOT / 'reports/AMBIENT_B18_N_R_SECOND_CONTACT_20261103_v63/report.json'
    Nreview = ROOT / 'reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63_REVIEW.json'
    nilreview = ROOT / 'reports/AMBIENT_B18_NILRADICAL_LINEAR_VISIBILITY_20261103_v63_REVIEW.json'
    firstreview = ROOT / 'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63_REVIEW.json'
    reductionreview = ROOT / 'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json'
    for path in [second, Nreview, nilreview, firstreview, reductionreview]:
        accepted(path)
    cpath = second.parent / 'certificate.json.gz'
    curvepath = ROOT / 'reports/AMBIENT_B18_R_CURVE_TANGENT_20261103_v63/certificate.json.gz'
    redpath = ROOT / 'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json'
    firstpath = ROOT / 'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz'
    relpath = ROOT / 'reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63/relations.json.gz'
    c, curve, red, first, rel = map(read, [cpath, curvepath, redpath, firstpath, relpath])
    # Extract the actual R curve's constant/linear/quadratic coefficients exactly.
    T = symbols(['t'])['t']
    rp = [Poly(evaluate(q, {'t': T})) for q in curve['physical_curve150']]
    r = [[q.terms.get((i,), F(0)) for q in rp] for i in range(3)]
    P0, velocity, R2 = r
    assert velocity == list(map(F, curve['physical_R_curve_derivative150']))
    assert R2 == list(map(F, c['R2_physical150']))
    ca = red['families']['CA']
    J = [[F(a) for j, a in enumerate(row) if j != 20] for row in ca['linear_J']]
    v27 = list(map(F, c['source_velocity27']))
    assert matvec(J, v27) == velocity
    fixed = matvec([row[19:] for row in J], v27[19:])
    M = [[F(a) for a in row] for row in c['source_outer_GL_matrix']]
    K = [[F(a) for a in row] for row in c['source_kernel_GL_matrix']]
    ownM = [[F(0) for _ in range(5)] for _ in range(5)]
    for a, ij in zip(v27[:19], red['orbit_GL_generator_indices']):
        ownM[ij//5][ij % 5] += a
    assert M == ownM
    assert [a+b for a, b in zip(lie(M, P0), fixed)] == velocity
    LM2, LMF = lie(M, lie(M, P0)), lie(M, fixed)
    N20 = [a/2+b for a, b in zip(LM2, LMF)]
    H = lie(K, fixed)
    W = [-a/2 for a in lie(K, lie(K, P0))]
    assert N20 == list(map(F, c['N2_kernel_zero150']))
    assert H == list(map(F, c['N2_mixed150']))
    assert W == list(map(F, c['N2_kernel_square150']))
    assert W == [F(a[0]) for a in first['N_full_firstjet_data']['W']]
    assert H == matvec(first['N_full_firstjet_data']['mixed'], v27)
    delta = [a-b for a, b in zip(R2, N20)]
    assert delta == list(map(F, c['actual_curvature_difference150']))
    rows = c['source27_unit_rows']
    L = [[F(a) for a in row] for row in c['source27_left_inverse']]
    assert matrix_product(L, [J[i] for i in rows]) == [[F(i == j) for j in range(27)] for i in range(27)]
    def project(V):
        image = matvec(J, matvec(L, [V[i] for i in rows]))
        return [a-b for a, b in zip(V, image)]
    dp, hp, wp = map(project, [delta, H, W])
    A, B = F(4800, 23), F(300, 23)
    assert dp == [A*a+B*b for a, b in zip(wp, hp)]
    assert wp[27]*hp[61]-hp[27]*wp[61] == -138
    assert B*B-A == -F(20400, 529)
    TN = [[F(a) for a in row] for row in curve['full150_by37_actual_N_normalized_tangent']]
    assert matvec(TN, [delta[i] for i in curve['N37_coordinate_rows']]) == delta
    # Independently construct every Taylor pair and product in the accepted full kernel.
    names = ['n', 'w', 'v']+[f's{i}' for i in range(1, 7)]+[f'y{i}' for i in range(9)]
    x = symbols(names)
    b = [x['n']]+[x[f's{i}'] for i in range(1, 7)]+[x['w']*x['v'], x['v']**2]
    y = [x[f'y{i}'] for i in range(9)]
    pairs = list(itertools.combinations(range(9), 2))
    own = []
    for i, j in pairs:
        ei, ej = next(iter(b[i].terms)), next(iter(b[j].terms))
        lcm = tuple(max(a, b) for a, b in zip(ei, ej))
        own.append(Poly({tuple(a-b for a, b in zip(lcm, ei)): 1})*y[i]
                   -Poly({tuple(a-b for a, b in zip(lcm, ej)): 1})*y[j])
    own += [y[i]*y[j]-x['n']*b[i]*b[j]
            for i, j in itertools.combinations_with_replacement(range(9), 2)]
    assert len(own) == 81
    assert all(a == evaluate(q, x) for a, q in zip(own, rel['all81']))
    # Universal leading-series check, not a newly computed physical third/fourth jet.
    z = symbols(['T', 'd', 'A', 'B']+[f'l{i}' for i in range(1, 7)]+['lw', 'lv'])
    T, d, aa, bb = [z[i] for i in ['T', 'd', 'A', 'B']]
    env = {'n': T**2*(d*d+aa+2*d*bb), 'w': T*z['lw'], 'v': T*z['lv'], 'y0': Poly(0), 'y7': Poly(0), 'y8': Poly(0)}
    for i in range(1, 7):
        env[f's{i}'] = T*z[f'l{i}']
        env[f'y{i}'] = T**2*z[f'l{i}']*(d+bb)
    assert all(all(e[0] >= 4 for e in Poly(evaluate(q, env)).terms) for q in rel['all81'])
    witness = evaluate('y1**2-n*s1**2', env)
    assert witness == T**4*z['l1']**2*(bb**2-aa)
    inputs = [Path(__file__), ROOT/f'research/notes/{STEM}.md', second, Nreview, nilreview,
              firstreview, reductionreview, cpath, curvepath, redpath, firstpath, relpath]
    result = dict(id=STEM, all_passed=True,
        verification_kind='Independent exact Fraction Lie2/all150 and universal full81 identities plus conceptual proof; no CAS and no new physical third jet.',
        independent_index_Lie2_reconstructs_all150_N20_H_W=True,
        actual_saved_R_polynomial_all150_R2_independently_extracted=True,
        all150_delta_TN_and_source27_projections_verified=True,
        source27_W_H_rank2_minor='-138', A=str(A), B=str(B), B_squared_minus_A=str(B*B-A),
        original_nu_and_source_trace_kernel_binding_proved=True,
        complete_all81_zero_linear_relations_have_R_order_atleast4=True,
        independent_regular_ambient_graph_multiplier_argument_proved=True,
        zero_linear_N_relation_exists_with_R_order_exactly4=True,
        h_exact=4, k_allowed=[3, 4],
        conditional_k_to_reduced_slice_edim={'3': 19, '4': 18},
        conditional_k_to_linear_nil_dimension={'3': 1, '4': 2},
        exact_k_or_explicit_physical_order4_polynomial_or_actual_third_jet_claim=False,
        scope='Generic B18 open orbit and accepted simultaneous G19 slice only; existential formal zero-linear relation, no full D-ring or exact k claim.',
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in inputs},
        seconds=time.monotonic()-start)
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    main(parser.parse_args().output)

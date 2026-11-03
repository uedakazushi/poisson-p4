#!/usr/bin/env python3
"""Spin4 biform slice and unexpanded mixed-discriminant circuit over QQ."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s

from weight1_structure import COORDS, J, form_matrices, infinitesimal_action, sp_basis

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_FIRSTJET_SPIN4_MIXED_DISCRIMINANT_20261103_v63'


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def rows(M):
    return [[str(x) for x in row] for row in M.tolist()]


class Circuit:
    def __init__(self, inputs):
        self.inputs, self.nodes, self.lookup = list(inputs), [], {}

    def node(self, op, value):
        key = (op, tuple(value) if isinstance(value, list) else value)
        if key not in self.lookup:
            self.lookup[key] = len(self.nodes)
            self.nodes.append([op, value])
        return self.lookup[key]

    def expr(self, expression):
        if expression.is_Rational:
            return self.node('q', str(expression))
        if expression.is_Symbol:
            return self.node('v', self.inputs.index(str(expression)))
        if expression.is_Add:
            return self.node('+', [self.expr(a) for a in expression.args])
        if expression.is_Mul:
            return self.node('*', [self.expr(a) for a in expression.args])
        if expression.is_Pow and expression.exp.is_Integer and expression.exp >= 0:
            return self.node('*', [self.expr(expression.base)] * int(expression.exp))
        raise ValueError(expression)

    def evaluate(self, values):
        computed = []
        for op, value in self.nodes:
            if op == 'q':
                x = Fraction(value)
            elif op == 'v':
                x = Fraction(values[value])
            elif op == '+':
                x = sum((computed[j] for j in value), Fraction(0))
            elif op == '*':
                x = Fraction(1)
                for j in value:
                    x *= computed[j]
            elif op == 'inv':
                x = 1 / computed[value]
            else:
                raise ValueError(op)
            computed.append(x)
        return computed


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    out.mkdir(parents=True, exist_ok=False)
    ip = ROOT / 'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    spinpath = ROOT / 'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    info, previous = read(ip), read(spinpath)
    B = s.Matrix(info['primitive24_by16_basis'])
    gamma = [s.Matrix(a) for a in previous['gamma_matrices']]
    ps = list(itertools.combinations(range(4), 2))[:-1]
    a0, a1, b0, b1, T, e = s.symbols('a0 a1 b0 b1 T e')
    u = s.symbols('u00 u01 u10 u11 u20 u21')
    v = s.symbols('v00 v01 v02 v10 v11 v12')
    variables = (*u, *v)
    U = sum(u[2 * i + j] * a0 ** (2 - i) * a1 ** i * (b0, b1)[j]
            for i in range(3) for j in range(2))
    V = sum(v[3 * i + j] * (a0, a1)[i] * b0 ** (2 - j) * b1 ** j
            for i in range(2) for j in range(3))
    psi = s.zeros(5, 4)
    for i in range(2):
        for j in range(2):
            rr = 1 + 2 * i + j
            psi[rr, 0], psi[rr, 1] = u[2 * i + j] / (1 if i == 0 else 2), u[2 * (i + 1) + j] / (2 if i == 0 else 1)
            psi[rr, 2], psi[rr, 3] = v[3 * i + j] / (1 if j == 0 else 2), v[3 * i + j + 1] / (2 if j == 0 else 1)
    assert sum((gamma[i] * psi[i, :].T for i in range(5)), s.zeros(4, 1)) == s.zeros(4, 1)
    beta = s.Matrix([psi[ps.index(pair), k] if pair != (2, 3) else 0 for pair, k in COORDS])
    E = s.Matrix(24, 12, lambda i, j: beta[i].coeff(variables[j]))
    assert E.rank() == 12
    br = list(B.T.rref()[1])
    E16 = B[br, :].inv() * E[br, :]
    assert B * E16 == E
    gs = sp_basis()
    indices = [i for i, A in enumerate(gs) if A[:2, 2:] == s.zeros(2) and A[2:, :2] == s.zeros(2)]
    assert len(indices) == 6
    action12 = []
    for index in indices:
        A = gs[index]
        aa, bb = A[:2, :2] * s.Matrix([a0, a1]), A[2:, 2:] * s.Matrix([b0, b1])
        dU = s.expand(sum(s.diff(U, x) * y for x, y in zip((a0, a1, b0, b1), (*aa, *bb))))
        dV = s.expand(sum(s.diff(V, x) * y for x, y in zip((a0, a1, b0, b1), (*aa, *bb))))
        target = [dU.coeff(a0, 2 - i).coeff(a1, i).coeff(b0, 1 - j).coeff(b1, j)
                  for i in range(3) for j in range(2)]
        target += [dV.coeff(a0, 1 - i).coeff(a1, i).coeff(b0, 2 - j).coeff(b1, j)
                   for i in range(2) for j in range(3)]
        R = s.Matrix(12, 12, lambda i, j: s.expand(target[i]).coeff(variables[j]))
        assert E * R == s.Matrix.hstack(*(infinitesimal_action(E[:, j], A) for j in range(12)))
        action12.append(R)
    g = s.Matrix.vstack(s.Matrix.hstack(s.zeros(2), s.eye(2)), s.Matrix.hstack(-s.eye(2), s.zeros(2)))
    assert g.T * J * g == J
    finite = []
    for j in range(12):
        forms = form_matrices(E[:, j])
        transformed = [g.T * sum((g[i, k] * forms[i] for i in range(4)), s.zeros(4)) * g for k in range(4)]
        finite.append(s.Matrix([transformed[k][i, j] for (i, j), k in COORDS]))
    swapped = [-v[0], -v[3], -v[1], -v[4], -v[2], -v[5], u[0], u[2], u[4], u[1], u[3], u[5]]
    swapmatrix = s.Matrix(12, 12, lambda i, j: s.expand(swapped[i]).coeff(variables[j]))
    assert s.Matrix.hstack(*finite) == E * swapmatrix
    delta = s.expand(psi[1:, :].det(method='domain-ge'))
    assert s.expand(delta.xreplace(dict(zip(variables, swapped))) + delta) == 0
    vv0, vv1 = sum(v[i] * T ** i for i in range(3)), sum(v[3 + i] * T ** i for i in range(3))
    F = s.Poly(s.expand((u[0] + u[1] * T) * vv1 ** 2 - (u[2] + u[3] * T) * vv0 * vv1
                        + (u[4] + u[5] * T) * vv0 ** 2), T)
    coefficients = [F.nth(i) for i in range(6)]
    R = s.resultant(vv0, vv1, T)
    restriction = s.Matrix(6, 6, lambda i, j: coefficients[i].coeff(u[j]))
    assert s.expand(restriction.det(method='domain-ge') - R ** 3) == 0
    ff = s.symbols('f0:6')
    disc = s.Poly(s.discriminant(sum(ff[i] * T ** i for i in range(6)), T), *ff)
    assert len(disc.terms()) == 59 and set(sum(m) for m, c in disc.terms()) == {8}
    circuit = Circuit(map(str, variables))
    f_nodes = [circuit.expr(f) for f in coefficients]
    R_node = circuit.expr(R)
    terms = []
    for monomial, coefficient in disc.terms():
        factors = [circuit.expr(coefficient)]
        factors += [f_nodes[i] for i, k in enumerate(monomial) for _ in range(k)]
        terms.append(circuit.node('*', factors))
    disc_node = circuit.node('+', terms)
    inv_R = circuit.node('inv', R_node)
    D_node = circuit.node('*', [disc_node, inv_R, inv_R])
    fixed = [2, 3, 5, 7, 11, 13, 1, 2, 3, 2, -1, 1]
    value = circuit.evaluate(fixed)
    swapped_fixed = list(swapmatrix * s.Matrix(fixed))
    value_swap = circuit.evaluate(swapped_fixed)
    assert value[D_node] == value_swap[D_node] == 53268630696
    assert value[R_node] == 50 and value_swap[R_node] == 37
    # A generic-point valuation fixture on the irreducible quadratic-resultant divisor.
    fixture = (1 + T) * (T ** 2 + e) ** 2 - (3 + 2 * T) * T * (T ** 2 + e) + (7 + 5 * T) * T ** 2
    fixture_discriminant = s.discriminant(fixture, T)
    quotient, remainder = s.div(fixture_discriminant, e ** 2, e)
    assert remainder == 0 and quotient.subs(e, 0) == 1466325
    # Bidegrees are propagated through the unexpanded arithmetic circuit.
    degrees = []
    for op, val in circuit.nodes:
        if op == 'q':
            degree = (0, 0)
        elif op == 'v':
            degree = (1, 0) if val < 6 else (0, 1)
        elif op == '*':
            degree = tuple(sum(degrees[j][i] for j in val) for i in range(2))
        elif op == '+':
            nonconstant = {degrees[j] for j in val if circuit.nodes[j] != ['q', '0']}
            assert len(nonconstant) <= 1
            degree = next(iter(nonconstant), (0, 0))
        else:
            degree = tuple(-x for x in degrees[val])
        degrees.append(degree)
    assert degrees[R_node] == (0, 4) and degrees[disc_node] == (8, 16) and degrees[D_node] == (8, 8)
    payload = dict(input_names=list(map(str, variables)), nodes=circuit.nodes,
                   outputs=dict(binary_quintic_coefficients=f_nodes, resultant=R_node,
                                quintic_discriminant=disc_node, mixed_discriminant=D_node),
                   mixed_discriminant_domain='resultant != 0',
                   extension_to_polynomial_is_proved_in_note_not_by_expanding_12_variables=True)
    with gzip.open(out / 'mixed-discriminant-circuit.json.gz', 'wt') as f:
        json.dump(payload, f)
        f.write('\n')
    cert = dict(variables12=list(map(str, variables)), physical_embedding24x12=rows(E),
                accepted_embedding16x12=rows(E16), spin4_sp4_generator_indices=indices,
                six_biform_action_matrices=list(map(rows, action12)), normalizer4=rows(g),
                normalizer_biform12=rows(swapmatrix), normal_component_delta=str(delta),
                normal_component_degree=4, detC_on_slice='delta^2/16',
                U=str(U), V=str(V), quintic_coefficients=list(map(str, coefficients)),
                resultant=str(R), restriction_matrix=rows(restriction),
                restriction_determinant='resultant^3', generic_quintic_discriminant59_terms=[
                    [list(m), str(c)] for m, c in disc.terms()],
                fixed_biform_values=fixed, fixed_swapped_values=list(map(str, swapped_fixed)),
                fixed_discriminant=str(value[D_node]), fixed_resultants=['50', '37'],
                resultant_boundary_fixture_quintic=str(fixture),
                fixture_discriminant_quotient_by_e2=str(quotient),
                fixture_quotient_at_e0='1466325',
                output_bidegrees=dict(quintic_discriminant=[8, 16], resultant=[0, 4], mixed=[8, 8]))
    (out / 'certificate.json').write_text(json.dumps(cert, indent=2) + '\n')
    result = dict(id=STEM, all_passed=True, field='QQ', exact_spin4_slice_dimension=12,
                  full_physical_biform_intertwinings=True, full_normalizer_identity=True,
                  mixed_discriminant_bidegree=[8, 8], mixed_discriminant_total_degree=16,
                  extraneous_resultant_factor_exact_power=2,
                  restriction_linear_map_determinant_is_resultant_cubed=True,
                  circuit_nodes=len(circuit.nodes), no_12_or_16_variable_discriminant_expansion=True,
                  normalizer_invariance_uses_irreducibility_and_one_normalization_fixture=True,
                  descended_rational_invariant_homogeneous_degree=16,
                  possible_polynomial_numerator_degrees='16+8k, k>=0',
                  pole_order_at_detC_or_k_equal_1_proved=False,
                  NEW25_correspondence_or_vanishing_proved=False,
                  no_NEW25_point_evaluation=True, no_degree22_or_higher_invariant_probe=True,
                  input_sha256={str(p.relative_to(ROOT)): sha(p) for p in
                                [Path(__file__), ROOT / 'code/weight1_structure.py', ip, spinpath]},
                  seconds=time.monotonic() - started)
    (out / 'report.json').write_text(json.dumps(result, indent=2) + '\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'circuit_nodes', 'seconds']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    main(parser.parse_args().output.resolve())

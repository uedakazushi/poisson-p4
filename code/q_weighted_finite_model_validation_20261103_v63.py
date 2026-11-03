#!/usr/bin/env python3
"""Independently evaluate the exact finite weighted model over QQ and QQ[eps].

This file imports neither the producer nor its symbolic reviewer.  Physical
divergence and all 350 Jacobi coefficients are computed by differentiation and
multiplication of five-variable exponent dictionaries.  The resonance test is
an identity in QQ[eps], with no truncation or specialization of eps.
"""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
PAIRS = list(itertools.combinations(range(5), 2))
MON2 = list(itertools.combinations_with_replacement(range(5), 2))
MON3 = list(itertools.combinations_with_replacement(range(5), 3))
TRIPLES = list(itertools.combinations(range(5), 3))
EXP2 = [tuple(mon.count(i) for i in range(5)) for mon in MON2]
EXP3 = [tuple(mon.count(i) for i in range(5)) for mon in MON3]
BASIS = list(itertools.product(PAIRS, EXP2))


class RationalDomain:
    zero = Fraction(0)
    one = Fraction(1)

    @staticmethod
    def q(value):
        return Fraction(value)

    @staticmethod
    def add(a, b):
        return a + b

    @staticmethod
    def mul(a, b):
        return a * b

    @staticmethod
    def neg(a):
        return -a

    @staticmethod
    def inv(a):
        assert a
        return 1 / a

    @staticmethod
    def encode(a):
        return str(a)


class EpsilonPolynomialDomain:
    """Univariate polynomials as trimmed coefficient tuples, never truncated."""
    zero = ()
    one = (Fraction(1),)
    eps = (Fraction(0), Fraction(1))

    @staticmethod
    def trim(values):
        values = list(values)
        while values and not values[-1]:
            values.pop()
        return tuple(values)

    @staticmethod
    def q(value):
        value = Fraction(value)
        return (value,) if value else ()

    @classmethod
    def add(cls, a, b):
        result = [Fraction(0)] * max(len(a), len(b))
        for i, value in enumerate(a):
            result[i] += value
        for i, value in enumerate(b):
            result[i] += value
        return cls.trim(result)

    @classmethod
    def mul(cls, a, b):
        if not a or not b:
            return ()
        result = [Fraction(0)] * (len(a) + len(b) - 1)
        for i, left in enumerate(a):
            for j, right in enumerate(b):
                result[i + j] += left * right
        return cls.trim(result)

    @staticmethod
    def neg(a):
        return tuple(-value for value in a)

    @staticmethod
    def inv(a):
        # A-only inverses must be nonzero constants at fixed a=2/9.
        assert len(a) == 1 and a[0], ('nonconstant/nonunit inverse', a)
        return (1 / a[0],)

    @staticmethod
    def encode(a):
        return [str(value) for value in a]


QQ = RationalDomain
QE = EpsilonPolynomialDomain


def digest(value):
    data = json.dumps(value, separators=(',', ':'), sort_keys=True).encode()
    return hashlib.sha256(data).hexdigest()


def evaluate(nodes, inputs, domain):
    values = []
    inverse_count = 0
    for index, node in enumerate(nodes):
        op = node[0]
        if op == 'q':
            value = domain.q(node[1])
        elif op == 'v':
            value = inputs[node[1]]
        elif op in ('+', '*'):
            value = domain.zero if op == '+' else domain.one
            operation = domain.add if op == '+' else domain.mul
            for ref in node[1]:
                assert 0 <= ref < index
                value = operation(value, values[ref])
        elif op == 'inv':
            assert 0 <= node[1] < index
            value = domain.inv(values[node[1]])
            inverse_count += 1
        else:
            raise ValueError(op)
        values.append(value)
    return values, inverse_count


def add_to(poly, exponent, value, domain):
    if value:
        result = domain.add(poly.get(exponent, domain.zero), value)
        if result:
            poly[exponent] = result
        else:
            poly.pop(exponent, None)


def differentiate(poly, direction, domain):
    result = {}
    for exponent, coefficient in poly.items():
        if exponent[direction]:
            reduced = list(exponent)
            reduced[direction] -= 1
            value = domain.mul(domain.q(exponent[direction]), coefficient)
            add_to(result, tuple(reduced), value, domain)
    return result


def add_product(result, left, right, domain):
    for a, c in left.items():
        for b, d in right.items():
            exponent = tuple(a[i] + b[i] for i in range(5))
            add_to(result, exponent, domain.mul(c, d), domain)


def physical_equations(raw150, domain):
    """Use only the physical tensor coefficients, not saved raw equation rows."""
    tensor = {(i, j): {} for i in range(5) for j in range(5)}
    for ((i, j), exponent), coefficient in zip(BASIS, raw150):
        if coefficient:
            tensor[i, j][exponent] = coefficient
            tensor[j, i][exponent] = domain.neg(coefficient)
    derivatives = {(i, j, k): differentiate(tensor[i, j], k, domain)
                   for i in range(5) for j in range(5) for k in range(5)}
    divergence = []
    for j in range(5):
        polynomial = {}
        for i in range(5):
            for exponent, coefficient in derivatives[i, j, i].items():
                add_to(polynomial, exponent, coefficient, domain)
        assert not polynomial, ('physical divergence', j, polynomial)
        divergence.extend([domain.zero] * 5)
    jacobi = []
    for i, j, k in TRIPLES:
        polynomial = {}
        for ell in range(5):
            add_product(polynomial, tensor[i, ell], derivatives[j, k, ell], domain)
            add_product(polynomial, tensor[j, ell], derivatives[k, i, ell], domain)
            add_product(polynomial, tensor[k, ell], derivatives[i, j, ell], domain)
        assert all(sum(exponent) == 3 for exponent in polynomial)
        jacobi.extend(polynomial.get(exponent, domain.zero) for exponent in EXP3)
    assert len(divergence) == 25 and len(jacobi) == 350
    return divergence, jacobi


def raw_point(a, domain):
    result = [domain.zero] * 150
    def put(pair, mon, coefficient):
        exponent = tuple(mon.count(i) for i in range(5))
        result[BASIS.index((pair, exponent))] = coefficient
    put((1, 4), (0, 0), a)
    put((2, 3), (0, 0), domain.q('-2/9'))
    put((2, 4), (0, 1), domain.q('4/3'))
    put((3, 4), (0, 2), domain.q('2/3'))
    put((3, 4), (1, 1), domain.one)
    return result


def recover_and_verify(model, inputs, values, domain):
    reconstructed = [values[index] for index in model['raw150_reconstruction_nodes']]
    residuals = [values[index] for index in model['all350_residual_nodes']]
    assert len(reconstructed) == 150 and len(residuals) == 350
    baseline = raw_point(inputs[0], domain)
    positive = [domain.add(reconstructed[index], domain.neg(baseline[index]))
                for index in model['positive_raw_indices']]
    assert len(positive) == 108
    for n, column in enumerate(model['free_positive_columns']):
        assert positive[column] == inputs[n + 1], ('free input', n, column)
    for index in model['group_cut_raw_indices']:
        assert reconstructed[index] == baseline[index], ('group cut', index)
    # Independently recover every 108-coordinate input from its raw free entry,
    # then test every saved canonical linear reconstruction coefficient.
    coordinates = {-1: domain.one, 0: inputs[0]}
    coordinates.update({i + 1: value for i, value in enumerate(positive)})
    for n, form in enumerate(model['canonical_linear_forms']):
        expected = domain.zero
        for coordinate, coefficient in form:
            expected = domain.add(expected, domain.mul(domain.q(coefficient),
                                                       coordinates[coordinate]))
        assert expected == reconstructed[n], ('canonical reconstruction', n)
    divergence, physical = physical_equations(reconstructed, domain)
    for n, (actual, expected) in enumerate(zip(residuals, physical)):
        assert actual == expected, ('DAG vs physical Jacobi', n)
    selected = [row for block in model['implicit_blocks'] for row in block['equation_rows']]
    assert len(selected) == 74 and len(set(selected)) == 74
    assert all(residuals[row] == domain.zero for row in selected)
    return reconstructed, residuals, dict(
        all150_reconstruction_coefficients_checked=True,
        all108_positive_coordinates_recovered=True,
        all34_free_inputs_recovered=True,
        all17_group_cuts_preserved=True,
        physical_divergence25_zero=True,
        all350_DAG_residuals_equal_independent_physical_Jacobi=True,
        all74_selected_equations_exactly_zero=True,
        nonzero_residual_count=sum(bool(value) for value in residuals),
        reconstruction_sha256=digest([domain.encode(value) for value in reconstructed]),
        residual_sha256=digest([domain.encode(value) for value in residuals]))


def resonance_curve():
    """Explicit Q_(2/9)+eps R-75eps^2 x4^2 d0^d1, from physical entries."""
    curve = raw_point(QE.q('2/9'), QE)
    # This complete physical R is deliberately independent of producer frames.
    terms = {
        (0, 1): [((0, 2), '4'), ((1, 1), '-3')],
        (0, 2): [((0, 3), '4'), ((1, 2), '12')],
        (0, 3): [((0, 4), '4'), ((1, 3), '12'), ((2, 2), '15/2')],
        (0, 4): [((1, 4), '-18')],
        (1, 2): [((0, 4), '-20/3'), ((1, 3), '4'), ((2, 2), '1')],
        (1, 3): [((1, 4), '-6'), ((2, 3), '6')],
        (1, 4): [((2, 4), '-4')],
        (2, 3): [((2, 4), '4'), ((3, 3), '6')],
        (2, 4): [((3, 4), '-4')],
        (3, 4): [((4, 4), '1')],
    }
    for pair, entries in terms.items():
        for mon, coefficient in entries:
            exponent = tuple(mon.count(i) for i in range(5))
            index = BASIS.index((pair, exponent))
            curve[index] = QE.add(curve[index], QE.mul(QE.q(coefficient), QE.eps))
    index = BASIS.index(((0, 1), (0, 0, 0, 0, 2)))
    curve[index] = QE.add(curve[index], QE.mul(QE.q(-75), QE.mul(QE.eps, QE.eps)))
    return curve


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', type=Path, required=True)
    parser.add_argument('--output-directory', type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic()
    args.output_directory.mkdir(parents=True, exist_ok=False)
    args.model = args.model.resolve()
    archive = args.model.read_bytes()
    decoded = gzip.decompress(archive)
    model = json.loads(decoded)
    assert model['input_names'] == ['a'] + [f'z{i}' for i in range(34)]
    assert len(model['nodes']) == 9647
    samples = []
    for sample, a in enumerate([Fraction(2, 9), Fraction(1), Fraction(2)]):
        # Every free variable is nonzero; the three deterministic samples differ
        # in all 34 positive coordinates as well as in the zero-weight a.
        inputs = [a] + [Fraction((-1) ** (i + sample) * (i + 2 + 3 * sample),
                                 2 * i + 3 + 2 * sample) for i in range(34)]
        values, inverses = evaluate(model['nodes'], inputs, QQ)
        reconstructed, residuals, checks = recover_and_verify(model, inputs, values, QQ)
        assert checks['nonzero_residual_count'] > 0  # Dense tests are off the Poisson locus.
        checks.update(a=str(a), input_values=list(map(str, inputs)),
                      inverse_nodes_evaluated=inverses,
                      maximum_output_numerator_bits=max(v.numerator.bit_length()
                                                       for v in reconstructed + residuals),
                      maximum_output_denominator_bits=max(v.denominator.bit_length()
                                                         for v in reconstructed + residuals))
        # Also check that the zero section is exactly the supplied Q_a.
        zero_values, _ = evaluate(model['nodes'], [a] + [QQ.zero] * 34, QQ)
        zero_raw = [zero_values[i] for i in model['raw150_reconstruction_nodes']]
        assert zero_raw == raw_point(a, QQ)
        assert all(zero_values[i] == 0 for i in model['all350_residual_nodes'])
        checks['zero_section_is_exactly_Q_a'] = True
        samples.append(checks)
        print('PASS dense QQ sample', a, '350 physical equations and 74 implicit zeros', flush=True)
    curve = resonance_curve()
    baseline = raw_point(QE.q('2/9'), QE)
    free_raw = [model['positive_raw_indices'][column]
                for column in model['free_positive_columns']]
    inputs = [QE.q('2/9')] + [QE.add(curve[index], QE.neg(baseline[index])) for index in free_raw]
    values, inverses = evaluate(model['nodes'], inputs, QE)
    reconstructed, residuals, curve_checks = recover_and_verify(model, inputs, values, QE)
    assert reconstructed == curve, 'ALL150 resonance coefficients must agree in QQ[eps]'
    assert not any(residuals), 'ALL350 resonance residuals must vanish in QQ[eps]'
    _, direct_jacobi = physical_equations(curve, QE)
    assert not any(direct_jacobi)
    curve_checks.update(
        coefficient_ring='QQ[eps], untruncated',
        all150_coefficients_equal_explicit_resonance_curve_identically=True,
        all350_residuals_zero_identically=True,
        independent_explicit_curve_physical_Jacobi350_zero=True,
        all_DAG_inverses_are_nonzero_constants=True,
        inverse_nodes_evaluated=inverses,
        maximum_reconstruction_eps_degree=max(len(v) - 1 for v in reconstructed),
        input_values=[QE.encode(v) for v in inputs],
        all150_reconstruction_coefficients=[QE.encode(v) for v in reconstructed])
    files = [Path(__file__).resolve(), args.model,
             ROOT / 'code/q_weighted_finite_model_20261103_v63.py',
             ROOT / 'code/q_resonance_integrable_direction_20261103_v62.py',
             ROOT / 'reports/Q_RESONANCE_INTEGRABLE_DIRECTION_20261103_v62.json']
    report = dict(
        id='Q_WEIGHTED_FINITE_MODEL_VALIDATION_20261103_V63',
        all_passed=True,
        no_producer_or_reviewer_imports=True,
        physical_method='5-variable exponent dictionaries with independent diff/mul',
        all_three_dense_QQ_samples_passed=True,
        samples=samples,
        resonance_curve=curve_checks,
        scope='Independent exact evaluations and resonance-curve polynomial identity; universal symbolic raw-row/block proof is the separate reviewer input.',
        universal_symbolic_residual_identity_claimed_from_samples=False,
        archive_sha256=hashlib.sha256(archive).hexdigest(),
        decoded_model_sha256=hashlib.sha256(decoded).hexdigest(),
        input_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in files},
        seconds=time.monotonic() - start)
    path = args.output_directory / 'validation.json'
    path.write_text(json.dumps(report, indent=2) + '\n')
    print('PASS exact QQ[eps] curve, ALL150 coefficients and ALL350 residuals; seconds',
          report['seconds'], flush=True)


if __name__ == '__main__':
    main()

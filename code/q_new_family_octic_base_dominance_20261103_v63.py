#!/usr/bin/env python3
"""One same-point source derivative proves NEW25 dominates the octic spectral base."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import resource
import signal
import sys
import time
from types import SimpleNamespace

import sympy as s
from sympy.polys.matrices import DomainMatrix

from q_new_family_dense_s61_map_20261103_v63 import differentiate
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63'
sys.set_int_max_str_digits(0)


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def rows(M):
    return [[str(x) for x in row] for row in M.tolist()]


def piv(M):
    return list(DomainMatrix.from_Matrix(M).convert_to(s.QQ).rref()[1])


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    out.mkdir(parents=True, exist_ok=False)
    mappath = ROOT / 'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/map.json.gz'
    invpath = ROOT / 'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    pairpath = ROOT / 'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63/certificate.json.gz'
    pointpath = ROOT / 'reports/Q_FIRSTJET_SPIN4_NEW25_POINT_20261103_v63/certificate.json.gz'
    spinpath = ROOT / 'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    payload, info, pair, oldpoint, spin = map(read, [mappath, invpath, pairpath, pointpath, spinpath])
    assert pair['same_point_index'] == oldpoint['sample_index'] == 0
    point = oldpoint['accepted_source_point']['point']
    dag = SimpleNamespace(nodes=payload['nodes'], names=payload['input_names'])
    values, derivatives = differentiate(dag, point)
    output = payload['output_u1_through_u61']
    U61 = s.Matrix([s.Rational(values[i]) for i in output])
    source_Jac = s.Matrix([[s.Rational(derivatives[i].get(j, 0)) for j in range(18)] for i in output])
    assert source_Jac.shape == (61, 18)
    pivot_value = values[payload['normalization']['value_pivot_node']]
    pf_value = values[payload['normalization']['value_pfaffian_node']]
    assert str(pivot_value) == oldpoint['accepted_source_point']['value_pivot']
    assert str(pf_value) == oldpoint['accepted_source_point']['value_pfaffian']
    model = S61()
    z, u, beta = primitive_beta(model)
    one = [i for i, weight in enumerate(model.reduced['parameter_weights']) if weight == 1]
    coefficients = [beta[i, j].coeff(t) for i in range(4) for j in range(i + 1, 4) for t in z]
    betaU = s.Matrix([[s.expand(c).coeff(v) for v in u] for c in coefficients])
    B = s.Matrix(info['primitive24_by16_basis'])
    br = piv(B.T)
    conversion = B[br, :].inv() * betaU[br, :]
    assert B * conversion == betaU
    b = conversion * U61[one, :]
    assert list(map(str, b)) == oldpoint['accepted_source_point']['beta16']
    beta_source_Jac = conversion * source_Jac[one, :]
    C, D = s.Matrix(pair['actual_C4']), s.Matrix(pair['actual_D4'])
    basis = [s.Matrix(A) for A in pair['pair_coordinate_basis4']]
    Jac_pair = s.Matrix(pair['jacobian20x16'])
    words = ['CCCC', 'CCCD', 'CCDD', 'CDCD', 'CDDD']
    gradient = s.zeros(5, 20)
    for k, word in enumerate(words):
        for column, letter in enumerate(['C'] * 10 + ['D'] * 10):
            E = basis[column % 10]
            answer = 0
            for pos, char in enumerate(word):
                if char != letter:
                    continue
                M = s.eye(4)
                for j, char in enumerate(word):
                    M = M * (E if j == pos else {'C': C, 'D': D}[char])
                answer += s.trace(M)
            gradient[k, column] = answer
    gradient_beta = gradient * Jac_pair
    matrix5x18 = gradient_beta * beta_source_Jac
    selected = piv(matrix5x18)
    rank = len(selected)
    assert rank <= 5
    witness = matrix5x18[:, selected]
    determinant = witness.det(method='domain-ge') if rank == 5 else None
    if rank == 5:
        assert determinant
    sp_columns = [i for i, name in enumerate(payload['input_names']) if name.startswith('g')]
    assert len(sp_columns) == 10 and matrix5x18[:, sp_columns] == s.zeros(5, 10)
    # Exact universal linear coordinate change from these five traces to pencil coefficients.
    change = s.Matrix([[0, 0, -s.Rational(1, 2), 0, -1],
                       [-s.Rational(1, 2), -1, 1, 0, 2],
                       [s.Rational(1, 2), 0, -2, -s.Rational(1, 2), -2],
                       [-s.Rational(1, 4), 0, s.Rational(1, 2), 0, 0],
                       [s.Rational(1, 32), 0, s.Rational(1, 16), 0, s.Rational(1, 8)]])
    assert change.det()
    trace_polys = [dict(P) for P in pair['trace_polynomials'][:5]]
    determinant_polys = [dict((m, Fraction(v)) for m, v in P)
                        for P in spin['canonical_adjoint_pencil']['determinant_coefficients']]
    all_keys = set().union(*(P.keys() for P in trace_polys), *(P.keys() for P in determinant_polys))
    for k in range(5):
        assert all(sum(Fraction(str(change[k, j])) * trace_polys[j].get(m, 0) for j in range(5))
                   == 256 * determinant_polys[k].get(m, 0) for m in all_keys)
    cert = dict(same_existing_point_index=0, same_source_point=point, input_names=payload['input_names'],
                source_map_nodes=len(payload['nodes']), all61_output_values=list(map(str, U61)),
                full61x18_source_Jacobian=rows(source_Jac), primitive61_weight1_indices=one,
                primitive_S61_to_beta16_linear_map=rows(conversion), beta16=list(map(str, b)),
                full_beta16x18_source_Jacobian=rows(beta_source_Jac), trace8_basis_words=words,
                trace8_gradient5x20_pair=rows(gradient), trace8_gradient5x16_beta=rows(gradient_beta),
                trace8_source_derivative5x18=rows(matrix5x18), differential_rank=rank,
                selected_source_columns=selected, selected_source_names=[payload['input_names'][j] for j in selected],
                exact5_minor=rows(witness) if rank == 5 else None,
                exact5_minor_determinant=str(determinant) if rank == 5 else None,
                trace8_to_pencil5_linear_change=rows(change), linear_change_determinant=str(change.det()),
                universal_full_octic_coefficient_change_verified=True,
                value_pivot=str(pivot_value), value_pfaffian=str(pf_value))
    with gzip.open(out / 'certificate.json.gz', 'wt') as f:
        json.dump(cert, f, indent=2)
        f.write('\n')
    inputs = [Path(__file__), ROOT / 'code/q_new_family_dense_s61_map_20261103_v63.py',
              ROOT / 'code/primitive_jet_elliptic_curve_20261101.py', ROOT / 'code/s61_model.py',
              mappath, invpath, pairpath, pointpath, spinpath]
    result = dict(id=STEM, all_passed=True, same_existing_point_only=True, new_beta_points=0,
                  full_actual6079_node_source_AD_and_beta_binding=True,
                  all_five_Inv8_source_differential_rank=rank,
                  exact5_minor_nonzero=rank == 5,
                  NEW25_dominates_octic_pencil_base_A5=rank == 5,
                  NEW25_dominates_full_characteristic_base_F4_p0_through_p4=rank == 5,
                  spectral_base_only_proper_divisors_cannot_define_NEW25=rank == 5,
                  no_degree24_or_higher_polynomial_evaluated=True,
                  NEW25_equation_or_multisection_degree_claimed=False,
                  input_sha256={str(p.relative_to(ROOT)): sha(p) for p in inputs},
                  seconds=time.monotonic() - started)
    (out / 'report.json').write_text(json.dumps(result, indent=2) + '\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'all_five_Inv8_source_differential_rank',
          'NEW25_dominates_octic_pencil_base_A5', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output.resolve())

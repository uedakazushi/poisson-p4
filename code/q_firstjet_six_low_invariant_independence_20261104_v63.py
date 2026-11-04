#!/usr/bin/env python3
"""One exact six-invariant Jacobian minor at the same frozen beta point."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_FIRSTJET_SIX_LOW_INVARIANT_INDEPENDENCE_20261104_v63'


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for a in iter(lambda: f.read(1 << 20), b''):
            h.update(a)
    return h.hexdigest()


def rows(M):
    return [[str(q) for q in row] for row in M.tolist()]


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    out.mkdir(parents=True, exist_ok=False)
    fixed = ROOT/'reports/Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63/certificate.json.gz'
    octic = ROOT/'reports/Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63/certificate.json.gz'
    pair = ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63/certificate.json.gz'
    gates = [ROOT/'reports/Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63_REVIEW/report.json',
        ROOT/'reports/Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63_REVIEW/report.json',
        ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63_REVIEW/report.json',
        ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63_REVIEW/report.json']
    assert all(read(p)['all_passed'] for p in gates)
    fc, oc, pc = map(read, [fixed, octic, pair])
    assert pc['same_point_index'] == oc['same_existing_point_index'] == 0
    C, D = s.Matrix(pc['actual_C4']), s.Matrix(pc['actual_D4'])
    G = s.Matrix(fc['universal_six_gradient6x16'])
    assert G.shape == (6, 16) and G[:5, :] == s.Matrix(oc['trace8_gradient5x16_beta'])
    assert G[5, :] == s.Matrix(fc['intrinsic_m_gradient1x16_beta'])
    columns = fc['selected_beta_columns']
    oldminor = G[:, columns]
    olddet = oldminor.det()
    assert olddet == s.Rational(fc['exact6_minor_determinant']) != 0
    change = s.Matrix(oc['trace8_to_pencil5_linear_change'])
    assert change.det() == s.Rational(oc['linear_change_determinant']) != 0
    I8 = s.Matrix([s.trace(C**4), s.trace(C**3*D), s.trace(C**2*D**2), s.trace(C*D*C*D), s.trace(C*D**3)])
    p = change*I8
    p0 = C.det()
    Delta8 = (s.trace(D**2)/2)**2-4*D.det()
    F4 = s.trace((-C+2*D)**2)
    m = s.Rational(fc['actual_same_point_intrinsic_m'])
    N12 = s.trace(C*D**2*C*D**2)-(s.trace(D**2)/2)*s.trace(C**2*D**2)+D.det()*s.trace(C**2)
    assert p[0] == p0 == Delta8 != 0
    assert F4 != 0 and F4**2 == 400*p0+1600*p[4]
    assert N12 == 2*p0*m
    # Globally: p0=Delta8 by the already-verified f3/Newton polynomial
    # identity, and F4²=400p0+1600p4 by the verified full pencil identity.
    # Differentiate these IDENTITIES, not numerical equalities at one point.
    outer = s.zeros(6)
    outer[0, 0], outer[0, 4] = 200/F4, 800/F4
    for i in range(1, 5):
        outer[i, i] = 1
    outer[5, 0], outer[5, 5] = 2*m, 2*p0
    first = s.diag(change, s.ones(1, 1))
    chain = outer*first
    assert chain.det() == change.det()*400*p0/F4 != 0
    new = chain*G
    assert new[:1, :] == (200/F4)*(change[0, :]+4*change[4, :])*G[:5, :]
    assert new[5, :] == 2*m*change[0, :]*G[:5, :]+2*p0*G[5, :]
    minor = new[:, columns]
    determinant = minor.det()
    assert determinant == chain.det()*olddet != 0
    evidence = dict(same_existing_point_index=0, selected_beta_columns=columns,
        invariant_order=['F4', 'p1', 'p2', 'p3', 'p4', 'N12'],
        old_order=['I8_0', 'I8_1', 'I8_2', 'I8_3', 'I8_4', 'm'],
        F4=str(F4), p0=str(p0), Delta8=str(Delta8), m=str(m), N12=str(N12),
        old_exact6minor=rows(oldminor), old_exact6minor_determinant=str(olddet),
        I8_to_pencil5=rows(change), derivative_chain6x6=rows(chain), derivative_chain_determinant=str(chain.det()),
        actual6x16_gradient=rows(new), actual6minor=rows(minor), actual6minor_determinant=str(determinant),
        global_derivative_identities=['dF4=(200/F4)(dp0+4dp4)', 'dN12=2m dp0+2p0 dm'],
        global_identity_dependencies=['detC=(tr(D^2)/2)^2-4detD', 'F4^2=400p0+1600p4', 'N12=2Delta8*m'])
    (out/'certificate.json').write_text(json.dumps(evidence, indent=2)+'\n')
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', fixed, octic, pair, *gates]
    result = dict(id=STEM, all_passed=True,
        exact_selected_F4_p1_p2_p3_p4_N12_are_algebraically_independent=True,
        one_exact6_by16_QQ_minor_at_same_frozen_beta_point=True,
        all_derivatives_obtained_from_frozen_full_gradients_by_verified_global_chain=True,
        F4_p0_nonzero_and_all6_minor_entries_saved=True,
        no_new_point_native_probe_or_high_degree_test=True,
        NEW_equation_or_degree24_relation_claimed=False,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'certificate.json': sha(out/'certificate.json')},
        seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'exact_selected_F4_p1_p2_p3_p4_N12_are_algebraically_independent', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())

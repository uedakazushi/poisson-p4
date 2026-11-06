#!/usr/bin/env python3
"""G: exact sparse direction incidence and bounded normal-elimination cost.

No Groebner basis or exhaustion claim. The numeric cost comparison uses the
archived T01 fibre over GF(32003). The incidence identities are checked over QQ
against all original 134 equations before any modular calculation.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction as F
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

import numpy as np
import sympy as sp
from sympy.polys.rings import ring

ROOT = Path(__file__).resolve().parents[1]
PRIME = 32003
FRAME = 'reports/S61_AFFINE_RESIDUAL42_20261103_v63/linear129x42-QQ.json.gz'
CROSS = 'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63/cross-input.txt.gz'
MODEL = 'legacy/v02/data/symplectic_slice_reduced_equations.json'
POINT = 'reports/S61_KNOWN_FIBER_AUDIT_20261031.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)


def mod(value):
    value = F(str(value))
    assert value.denominator % PRIME
    return value.numerator * pow(value.denominator, -1, PRIME) % PRIME


def serialize(poly):
    return [[list(mon), str(c)] for mon, c in sorted(poly.items())]


def add(poly, mon, c):
    c = (poly.get(mon, 0) + c) % PRIME
    if c:
        poly[mon] = c
    else:
        poly.pop(mon, None)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--source', type=Path, required=True)
    args = ap.parse_args()
    resource.setrlimit(resource.RLIMIT_AS, (1536 * 1024**2,) * 2)
    resource.setrlimit(resource.RLIMIT_CPU, (90, 91))
    started = time.monotonic()
    data, model, records = (load(ROOT / p) for p in (FRAME, MODEL, POINT))
    point = list(map(sp.QQ, next(p for p in records['points'] if p['component_id'] == 'T01')['S61_coordinates_u1_through_u61']))
    qids = [i for i, w in enumerate(model['parameter_weights']) if w == 2]
    qi = {old: j for j, old in enumerate(qids)}
    W = [[sp.QQ(c) for c in row] for row in data['embedding']]
    Q, *x = ring(','.join('x'+str(i) for i in range(45)), sp.QQ)
    forms = [Q(point[old]) + sum((W[j][h] * x[h] for h in range(3)), Q.zero)
             + (x[j] if j >= 3 else Q.zero) for j, old in enumerate(qids)]
    physical = []
    for eq in model['equations']:
        poly = Q.zero
        for mon, numerator, denominator in eq:
            term = Q(sp.QQ(numerator, denominator))
            for old in mon:
                term *= forms[qi[old]] if old in qi else point[old]
            poly += term
        physical.append(poly)
    assert len(physical) == 134
    cross, linear, quadratic, incidence = [], [], [], []
    normal_mons = [(j,) for j in range(42)] + list(itertools.combinations_with_replacement(range(42), 2))
    pure_changes = [[sp.QQ(c) for c in row] for row in data['pure_plane_subtraction_coefficients']]
    for i in range(129):
        residual = physical[i+5] - sum((c * physical[j] for j, c in enumerate(pure_changes[i])), Q.zero)
        given = Q.zero
        for j, terms in enumerate(data['all129_linear_factor_rows'][i]):
            for mon, c in terms:
                given += x[j+3] * sp.QQ(c) * (x[mon[0]] if mon else Q.one)
        assert residual == given
        row_cross = [sp.QQ.zero] * 126
        row_linear, row_quadratic = {}, {}
        for exp, coeff in residual.items():
            zword = tuple(h for h in range(3) for _ in range(exp[h]))
            nword = tuple(j for j in range(42) for _ in range(exp[j+3]))
            assert nword and len(zword)+len(nword) <= 2
            if zword:
                assert len(zword) == len(nword) == 1
                row_cross[42*zword[0]+nword[0]] = coeff
            elif len(nword) == 1:
                row_linear[nword] = coeff
            else:
                row_quadratic[nword] = coeff
        cross.append(row_cross)
        linear.append(row_linear)
        quadratic.append(row_quadratic)
        incidence.append({
            'z_columns3_linear_in_v': [serialize({(j,): row_cross[42*h+j] for j in range(42) if row_cross[42*h+j]}) for h in range(3)],
            'r_column_quadratic_in_v': serialize(row_quadratic),
            'constant_column_linear_in_v': serialize(row_linear),
        })
    rest = [[row1.get(mon, row2.get(mon, sp.QQ.zero)) for mon in normal_mons] for row1, row2 in zip(linear, quadratic)]
    raw = ('129 126 945\n' + '\n'.join(' '.join(str(c) for c in row) for row in cross+rest) + '\n').encode()
    assert raw == gzip.decompress((ROOT/CROSS).read_bytes())
    # Evaluate the SAME cross solve only to quantify the density it creates.
    C = np.array([[mod(c) for c in row] for row in cross], dtype=np.int64)
    D = np.array([[mod(c) for c in row] for row in rest], dtype=np.int64)
    aug = np.concatenate((C.copy(), np.eye(129, dtype=np.int64)), axis=1)
    pivots = []
    for j in range(126):
        pivot = next(k for k in range(j, 129) if aug[k, j])
        aug[[j, pivot]] = aug[[pivot, j]]
        aug[j] = aug[j] * pow(int(aug[j, j]), -1, PRIME) % PRIME
        for k in range(129):
            if k != j and aug[k, j]:
                aug[k] = (aug[k] - aug[k, j] * aug[j]) % PRIME
        pivots.append(pivot)
    assert np.array_equal(aug[:126,:126], np.eye(126, dtype=np.int64))
    assert not aug[126:,:126].any()
    transform = aug[:,126:]
    P = (-transform[:126] @ D) % PRIME
    H = transform[126:] @ D % PRIME
    assert np.array_equal(transform[:126] @ C % PRIME, np.eye(126, dtype=np.int64))
    assert not (transform[126:] @ C % PRIME).any()
    residual_matrix = (C @ P + D) % PRIME
    assert not (transform[:126] @ residual_matrix % PRIME).any()
    assert np.array_equal(transform[126:] @ residual_matrix % PRIME, H)
    polynomials = [{m: int(c) for m, c in zip(normal_mons, row) if c} for row in P]
    moments = []
    for h in range(3):
        m = dict(polynomials[42*h])
        for mon, c in polynomials[42*h+34].items():
            add(m, mon, c)
        moments.append(m)
    compat = []
    for h in range(3):
        for j in range(42):
            if j == 34:
                continue  # C_h34 = -C_h0 for ell=n0+n34.
            poly = {}
            for mon, c in polynomials[42*h+j].items():
                for k in (0, 34):
                    add(poly, tuple(sorted(mon+(k,))), c)
            for mon, c in moments[h].items():
                add(poly, tuple(sorted(mon+(j,))), -c)
            compat.append(poly)
    assert len(compat) == 123
    source = {
        'task_id': 'G', 'coefficient_field': 'QQ at archived T01 beta',
        'normalization': 'n=r*v; ell(v)=v0+v34=1; r inverted through 1-u*r',
        'matrix129x5_row_order': ['z0','z1','z2','r','constant'],
        'matrix129x5': incidence,
        'selected_F5_in_z_n': [serialize(p) for p in physical[:5]],
        'selected_F5_encoding': '45-entry exponent tuple: first3 z then42 n; substitute n_j=r*v_j',
        'all134_identity': 'F_(i+5)=sum_j pure_change[i,j]*F_j + r*matrix_row_i(v)*(z0,z1,z2,r,1)^t',
        'pure_changes129x5': [[str(c) for c in row] for row in pure_changes],
        'maps': {'forward': 'n=r*v', 'inverse': 'r=n0+n34; v=n/r; u=1/r; z unchanged'},
        'all134_QQ_reconstructed_and129_residual_identities_checked': True,
        'rank_statement': 'First3 columns rank3 for every nonzero v, by injective X129x126. Rank of full matrix<=4 is only necessary without the affine kernel condition and F5.',
    }
    args.source.parent.mkdir(parents=True, exist_ok=True)
    with args.source.open('x') as f:
        json.dump(source, f, indent=2); f.write('\n')
    report = {
        'task_id': 'G', 'all_passed': True,
        'all134_QQ_reconstructed_and129_residual_identities_checked': True,
        'full_cross_input_reconstruction_byte_equal': True,
        'fixed_beta_ONLY': True, 'modular_cost_field': 'GF(32003)',
        'cross_matrix_rank126': True,
        'source_incidence_terms': sum(sum(bool(c) for c in row) for row in cross+rest),
        'original134_in_z_n_terms': sum(len(p) for p in physical),
        'inverse_P126_terms': sum(map(len, polynomials)),
        'normal_H3_terms': int(np.count_nonzero(H)),
        'principal_eta_M3_terms': list(map(len, moments)),
        'normal_compat123_terms': sum(map(len, compat)),
        'normal_compat123_degrees': dict(Counter(len(m) for p in compat for m in p)),
        'normal_only_before_five_quartics_terms': int(np.count_nonzero(H))+sum(map(len, compat))+3,
        'quartics_expanded': False,
        'no_unit_no_GB_no_generic_exhaustion': True,
        'decision': 'Keep sparse direction incidence; do not launch dense normal-only GB solely because it removes three variables.',
        'source_path': str(args.source), 'source_sha256': sha(args.source),
        'input_sha256': {p: sha(ROOT/p) for p in (FRAME,CROSS,MODEL,POINT,'code/classification_normal_incidence.py')},
        'seconds': time.monotonic()-started,
        'resource_limit': {'cpu_seconds':90,'address_space_mib':1536},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as f:
        json.dump(report, f, indent=2); f.write('\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()

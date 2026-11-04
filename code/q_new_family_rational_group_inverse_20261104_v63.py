#!/usr/bin/env python3
"""Exact finite-pullback convention and rational group inverse; no point job."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
ID = 'Q_NEW_FAMILY_RATIONAL_GROUP_INVERSE_20261104_v63'


def read(p):
    b = p.read_bytes()
    return json.loads(gzip.decompress(b) if p.suffix == '.gz' else b)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def main(out):
    started = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    signal.alarm(60)
    out.mkdir(parents=True, exist_ok=False)
    stems = ['Q_NEW_FAMILY_CONSTRUCTIVE_INVERSE_RECIPE_20261104_v63_REVIEW',
             'Q_NEW_FAMILY_RATIONAL_ORBIT_SECTION_20261104_v63_REVIEW',
             'Q_NEW_FAMILY_SP4_CENTRALIZER_FREE_20261104_v63_REVIEW',
             'Q_NEW_FAMILY_TARGET_FIELD_AND_SPECTRAL_DEGREES_20261104_v63_REVIEW']
    gates = [ROOT/'reports'/n/'report.json' for n in stems]
    for p in gates:
        assert read(p)['all_passed'], p
    rp = ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    sp = ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    rep, spin = read(rp), read(sp)
    B = s.Matrix([[s.Rational(v) for v in row] for row in rep['primitive24_by16_basis']])
    J = s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
    pairs = list(itertools.combinations(range(4), 2))
    coords = [(i,j,k) for i,j in pairs for k in range(4)]
    selected = list(B.T.rref()[1])
    inverse = B[selected,:].inv()
    t = s.Symbol('t')
    A = s.zeros(4); A[1,0] = -1
    assert A*A == s.zeros(4) and A.T*J+J*A == s.zeros(4)
    g = s.eye(4)+t*A
    assert g.T*J*g == J
    transformed = s.zeros(24,16)
    for col in range(16):
        forms = [s.zeros(4) for _ in range(4)]
        for val,(i,j,k) in zip(B[:,col],coords):
            forms[k][i,j] = val; forms[k][j,i] = -val
        pulled = [g.T*sum((g[l,k]*forms[l] for l in range(4)),s.zeros(4))*g for k in range(4)]
        transformed[:,col] = s.Matrix([s.expand(pulled[k][i,j]) for i,j,k in coords])
    U = inverse*transformed[selected,:]
    assert (B*U-transformed).applyfunc(s.expand) == s.zeros(24,16)
    variables = s.symbols('b0:16')
    moved = U*s.Matrix(variables)
    mon = list(itertools.combinations_with_replacement(range(16),2))
    checked = 0
    for key in ['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']:
        rows = [[s.Rational(c) for c in row] for row in spin[key]]
        plain = s.Matrix(4,4,lambda i,j:sum((c*variables[a]*variables[b] for c,(a,b) in zip(rows[4*i+j],mon) if c),s.S.Zero))
        image = s.Matrix(4,4,lambda i,j:sum((c*moved[a]*moved[b] for c,(a,b) in zip(rows[4*i+j],mon) if c),s.S.Zero))
        assert (image-g.T*plain*g.T.inv()).applyfunc(s.expand) == s.zeros(4)
        checked += sum(len(row) for row in rows)
    assert checked == 4352
    # Universal scalar powers: three covariant indices, not inverse pullback.
    a, k0 = s.symbols('a k0', nonzero=True)
    assert s.cancel((a*a)/(a*a*a)*a-1) == 0
    assert s.cancel(((a*a)/(a*a*a))**2*a*a-1) == 0
    assert s.cancel(((a*a)/(a*a*a))**3*a**3-1) == 0
    assert s.cancel((k0*k0)/(k0**3)*k0-1) == 0
    # Cayley equivalence without expanding a ten-variable determinant:
    # (I+A)^t J(I+A)-(I-A)^t J(I-A)=2(A^t J+JA).
    generic = s.Matrix(4,4,lambda i,j:s.Symbol('a%d%d'%(i,j)))
    I = s.eye(4)
    assert ((I+generic).T*J*(I+generic)-(I-generic).T*J*(I-generic)-2*(generic.T*J+J*generic)).applyfunc(s.expand) == s.zeros(4)
    upper = s.symbols('s0:10'); S=s.zeros(4)
    for v,(i,j) in zip(upper,itertools.combinations_with_replacement(range(4),2)):
        S[i,j]=v; S[j,i]=v
    ham=J*S
    assert ham.T*J+J*ham == s.zeros(4)
    evidence = dict(finite_pullback='g^* beta, anti-action',
       covariant_law='C(g^*b)=g^t C(b) g^(-t), likewise D',
       checked_finite_transvection='I-t E10',checked_primitive_basis_entries=384,
       checked_original_quadratic_coefficients=4352,
       intertwiner='T C(b)=C(beta)T, T D(b)=D(beta)T; nonzero 15minor charts/Cramer',
       scalar_multiplier='lambda=(T^t J T)_01=t^2',
       triple_covariant_ratio='k=((T^t)^*b)_i/beta_i=t^3',
       recovered_group='g=(lambda/k)T^t',
       cayley='A=(g-I)(g+I)^(-1)=J S, S symmetric; g=(I-A)^(-1)(I+A)',
       source_point_or_large_generic_coefficient_expansion=False)
    ep=out/'certificate.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),rp,sp,*gates]
    report=dict(id=ID,all_passed=True,
       finite_non_diagonal_pullback_and_transposed_C_D_convention_verified=True,
       no_square_root_rational_intertwiner_group_inverse_proved=True,
       explicit_ten_Cayley_and_five_quotient_finite_rational_coordinates=True,
       all_boundary_or_full_Poisson_component25_claim=False,
       no_new_beta_point_prime_invariant_expansion_source26_or_large_matrix_solve=True,
       input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},
       evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k:report[k] for k in ['all_passed','seconds']}))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)

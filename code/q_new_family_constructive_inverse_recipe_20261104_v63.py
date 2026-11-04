#!/usr/bin/env python3
"""Finite inverse operator recipe; only re-use the saved 100-entry old matrix.

No polynomial coefficient expansion, new point, norm, or source26 evaluation.
The emitted JSON embeds all polynomial inputs for the finite operator recipe.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_NEW_FAMILY_CONSTRUCTIVE_INVERSE_RECIPE_20261104_v63'


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def solve_two(M, p):
    """Solve exactly M*a=e_r and M*b=e_1 in the SAME saved F_p matrix."""
    n = len(M)
    A = [[x % p for x in row] + [int(i == 1), int(i == 0)]
         for i, row in enumerate(M)]
    pivots = []
    for k in range(n):
        j = next(i for i in range(k, n) if A[i][k])
        if j != k:
            A[j], A[k] = A[k], A[j]
        u = A[k][k]
        pivots.append(u)
        A[k] = [x * pow(u, -1, p) % p for x in A[k]]
        for i in range(n):
            if i != k:
                z = A[i][k]
                A[i] = [(x - z*y) % p for x, y in zip(A[i], A[k])]
    a, b = [row[n] for row in A], [row[n+1] for row in A]
    images = [[sum(M[i][j]*v[j] for j in range(n)) % p
               for i in range(n)] for v in (a, b)]
    assert images == [[int(i == 1) for i in range(n)],
                      [int(i == 0) for i in range(n)]]
    return dict(a=a, b=b, M_a=images[0], M_b=images[1],
                all20_coefficient_identities_verified=True,
                elimination_pivots=pivots)


def j14_operator(beta, monomials, C_coeff, D_coeff, omega_inverse):
    """Exact finite field operations for {tr(C^4),tr(C^3 D)} at beta.

    Valid over ANY characteristic-zero field containing beta (including the
    rational shape field).  This function is deliberately NOT called on the
    symbolic source: the inverse packet stores its finite algorithm instead.
    Coefficients must first be coerced into the chosen field.
    """
    z = beta[0] * 0

    def mm(A, B):
        return [[sum((A[i][k]*B[k][j] for k in range(4)), z)
                 for j in range(4)] for i in range(4)]

    def tr(A):
        return sum((A[i][i] for i in range(4)), z)

    matrices, derivatives = [], []
    for coefficients in (C_coeff, D_coeff):
        v = [sum((c*beta[i]*beta[j] for c, (i, j)
                  in zip(row, monomials)), z) for row in coefficients]
        matrices.append([v[4*i:4*i+4] for i in range(4)])
        deriv = []
        for k in range(16):
            w = [sum((c*((beta[j] if i == k else z) +
                         (beta[i] if j == k else z))
                      for c, (i, j) in zip(row, monomials)), z)
                 for row in coefficients]
            deriv.append([w[4*i:4*i+4] for i in range(4)])
        derivatives.append(deriv)
    C, D = matrices
    dC, dD = derivatives
    C2 = mm(C, C)
    C3 = mm(C2, C)
    pieces = [mm(C2, D), mm(mm(C, D), C), mm(D, C2)]
    K = [[sum((A[i][j] for A in pieces), z) for j in range(4)]
         for i in range(4)]
    grad0 = [4*tr(mm(A, C3)) for A in dC]
    grad1 = [tr(mm(A, K)) + tr(mm(C3, B)) for A, B in zip(dC, dD)]
    return sum((grad0[i]*omega_inverse[i][j]*grad1[j]
                for i in range(16) for j in range(16)), z)


def main(out):
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    signal.alarm(60)
    out.mkdir(parents=True, exist_ok=False)
    stems = {
        'krylov': 'Q_NEW_FAMILY_ZN_PRIMITIVE_KRYLOV_20261104_v63',
        'krylov_review': 'Q_NEW_FAMILY_ZN_PRIMITIVE_KRYLOV_20261104_v63_REVIEW',
        'content_review': 'Q_NEW_FAMILY_TWO_SHAPE_CONTENT_PRIMITIVITY_20261104_v63_REVIEW',
        'scale_review': 'Q_NEW_FAMILY_GENERAL_SOURCE_SCALE_20261104_v63_REVIEW',
        'lower26_review': 'Q_NEW_FAMILY_FULL_DEGREE24_SOURCE_RELATION_20261104_v63_REVIEW',
        'inv14_review': 'Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63_REVIEW',
        'shape_review': 'Q_NEW_FAMILY_TWO_SHAPE_RECOVERY_20261104_v63_REVIEW',
        'source_review': 'Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW',
        'compression_review': 'Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63_REVIEW',
    }
    reports = {k: ROOT/'reports'/v/'report.json' for k, v in stems.items()}
    for p in reports.values():
        assert read(p)['all_passed'], p
    krylov = ROOT/'reports'/stems['krylov']/'certificate.json'
    spin = ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    inv = ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    G = ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_CONTENT_GCD_20261104_v63/primitive.txt'
    C1 = ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_PRS_CONTENT_PREFLIGHT_20261104_v63/C1-primitive.txt.gz'
    shape = ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_HOMOGENIZATION_20261104_v63/certificate.json'
    orbit = ROOT/'reports/Q_NEW_FAMILY_RATIONAL_ORBIT_SECTION_20261104_v63/certificate.json'
    circuit = ROOT/'reports/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63/source-five-circuit.json.gz'
    invdefs = ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63/definitions.json'
    minors = ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63/minor-witnesses.json'
    k, sd, iv, sh, orb, ci, defs, mi = map(read, [krylov, spin, inv, shape, orbit, circuit, invdefs, minors])
    p = k['prime']
    assert p == 32003 and k['oldpoint_index'] == 0 and k['Krylov_det_mod32003'] == 21349
    solves = solve_two(k['Krylov_matrix_rows'], p)
    # Literal saved G and BOTH full C1 coefficients, never re-expanded.
    lines = G.read_text().splitlines()
    assert list(map(int, lines[0].split())) == [6, 891]
    Gterms = [list(map(int, line.split())) for line in lines[1:]]
    assert len(Gterms) == 891 and all(v[1] == 0 for v in Gterms)
    with gzip.open(C1, 'rt') as f:
        assert f.readline().split() == ['5', '2']
        Cterms = []
        for n in (10314, 10698):
            assert int(f.readline()) == n
            Cterms.append([list(map(int, f.readline().split())) for _ in range(n)])
        assert not f.read().strip()
    O = S.Matrix(iv['symplectic_form16'])
    Oi = O.inv()
    assert O*Oi == S.eye(16)
    assert defs['Inv8_trace_words'][:2] == ['0000', '0001']
    assert mi['14']['columns'][0] == 0 and mi['14']['minor_mod32003'] == 7584
    assert defs['degree14'][0] == dict(kind='Poisson_Inv8_Inv8', left=0, right=1)
    assert len(orb['beta_bar16']) == 16 and all(z['quadratic_remainder_odd_part_zero'] for z in orb['all16_checks'])
    assert ci['input_names'] == ['u', 'c', 'X', 'W', 'B']
    operators = [
        dict(name='base', outputs=['R','Z0','Z1','Z2','eta'],
             expressions=['(p3+4*p4)/p0','p0/F4^2','p1/F4^2','p2/F4^2','N12/F4^3']),
        dict(name='quotient_field', coefficient_field='K=QQ(R,Z0,Z1,Z2)',
             algebra='L=K[r]/(monic(G))', basis=['r^'+str(i) for i in range(10)],
             inputs='G_terms_embedded, C1_terms_embedded, N79_N41_Fn_J_embedded'),
        dict(name='balanced_pair', C1='A(r)*s+B(r); constant FIRST, linear SECOND',
             hom='H_d(P)=sum p_ijk*r^i*R^k*A^(d-j)*(-B)^j, for s-degree j',
             U='Z0*r*(r-1)^2*H1(J)*H3(N79)+25*Z0^2*H2(Fn)*H2(N41)',
             V='20*r*(r-1)^2*H2(Fn)*H1(J)^2',
             reduction='all products reduced modulo monic(G), using Euclidean division'),
        dict(name='matrix', columns='coefficient vector of U^j*V^(9-j), j=0,...,9',
             rows='1,r,...,r^9', determinant='Delta=det(M), nonzero by accepted Krylov certificate'),
        dict(name='two_linear_solves', a='adj(M)*e_r/Delta', b='adj(M)*e_1/Delta',
             e_r=[0,1]+[0]*8, e_1=[1]+[0]*9,
             adjugate='adj(M)_{ij}=(-1)^(i+j)*det(M deleting row j,column i); det via signed permutation sum'),
        dict(name='recover_r', expression='(sum a_j*eta^j)/(sum b_j*eta^j)',
             denominator_nonzero='sum b_j*eta^j=V^(-9) in generic field'),
        dict(name='recover_s', expression='-B(r)/A(r)',
             unit='A(r) nonzero by accepted birational degree3/5 linear C1 recovery'),
        dict(name='recover_q', T='(r-1)*R-2*r', Lq='8*R*r*(r-1)-6*r^2-11*r+1',
             Nq='(s-1)*Lq-(r-1)', q='Nq/(40*r*T)'),
        dict(name='recover_B', h='s^2-25*r*q^2',
             F='100*r*(r-1)*q^2-4*r*s+s^2+4*r',
             D='2359296000*r^4*(r-1)*q^2*(-s+5*q+1)', B='Z0*h*F^2/D'),
        dict(name='h14', J14='gradient_beta(tr(C^4))^T*Omega^(-1)*gradient_beta(tr(C^3*D))',
             derivative_order='differentiate all16 ambient beta coordinates BEFORE source substitution',
             beta='embedded rational orbit beta_bar16 at xi=1',
             output='h14(r,s,q,B)=J14(beta_bar16(1,r,s,q,B))',
             gradient_algorithm='j14_operator in frozen producer code; quadratic C/D arrays and Omega inverse embedded',
             fallback='original embedded source-five circuit at u^2=1600r,X=s/(48u),c=B*u/40,W=(q+1)/3072; trace/2 in quadratic etale extension yields same h14 by known tau Sp4 action',
             nonzero='J14 universal nonzero from Inv14 minor7584; NEW defining degree>=26, so J14 restriction cannot vanish'),
        dict(name='recover_xi', f4='(2/405)*F', expression='(J14_target/F4_target^3)*f4^3/h14',
             proof='J14_source=xi^7*h14 and F4_source=xi^2*f4 by accepted all16 cocharacter scaling'),
    ]
    payload = dict(
        coefficient_field='QQ', operator_program=operators,
        target_required_values=['F4','p0','p1','p2','p3','p4','N12','J14'],
        generic_open='source chart plus F4*p0*Delta*A(r)*V*h14 nonzero; rational recipe, no boundary extension',
        G_terms_format=['integer','s','r','R','Z0','Z1','Z2'], G_terms_embedded=Gterms,
        C1_terms_format=['integer','Z2','Z1','Z0','R','r'],
        C1_terms_constant_then_linear_embedded=Cterms,
        N79_N41_Fn_J_embedded={
            'N79': sh['rational_coefficient_pairs']['E1']['numerator'],
            'N41': sh['rational_coefficient_pairs']['E2']['numerator'],
            'Fn': sh['factors']['Fn'], 'J': sh['factors']['J']},
        global_quadratic_monomials=sd['quadratic_monomials'],
        global_C_coefficients16x136=sd['spin_adjoint_C_coefficients16x136'],
        global_D_coefficients16x136=sd['vector_adjoint_D_coefficients16x136'],
        global_Omega=iv['symplectic_form16'],
        global_Omega_inverse=[[str(Oi[i,j]) for j in range(16)] for i in range(16)],
        rational_orbit_section=orb,
        original_source_five_circuit=ci,
        word_trace_derivative_formulas=dict(
            T0='tr(C^4)', T1='tr(C^3 D)',
            dT0='4 tr(dC C^3)',
            dT1='tr(dC (C^2 D+C D C+D C^2))+tr(C^3 dD)',
            bracket='sum_i,j dT0_i Omega_inverse_ij dT1_j'),
        no_generic_expanded_adjugate_h14_or_defining_equation_computation=True)
    encoded=json.dumps(payload,separators=(',',':')).encode()
    assert len(encoded)<32*1024*1024
    with gzip.open(out/'recipe.json.gz','wb') as f:
        f.write(encoded+b'\n')
    certificate=dict(prime=p, oldpoint_index=0,
                     provenance='UNCHANGED matrix from accepted Krylov certificate; no beta or target ratio re-evaluation',
                     saved_matrix=k['Krylov_matrix_rows'], saved_determinant=21349,
                     **solves)
    (out/'old_matrix_solves.json').write_text(json.dumps(certificate,indent=2)+'\n')
    inputs=[Path(__file__),krylov,spin,inv,G,C1,shape,orbit,circuit,invdefs,minors,*reports.values()]
    report=dict(id=STEM,all_passed=True,complete_finite_rational_inverse_operator_recipe=True,
                generic_coefficient_expansions_performed=False,
                same_saved_old_matrix_two_solves_verified=True,
                new_beta_points=0, new_target_points=0, new_primes=0,
                no_source26_Inv14_source_expansion_norm_or_defining_equation_job=True,
                h14_direct_rational_orbit_recipe_with_prior_quadratic_trace_fallback=True,
                target_field_equality_proof='separate manuscript companion; this packet supplies its finite inverse operators',
                uncompressed_recipe_bytes=len(encoded),
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs},
                evidence_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()},
                seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k:report[k] for k in ('all_passed','complete_finite_rational_inverse_operator_recipe','uncompressed_recipe_bytes','seconds')},indent=2))


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output.resolve())

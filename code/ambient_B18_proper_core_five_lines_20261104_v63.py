#!/usr/bin/env python3
"""Fixed polynomial support and generic Artin models of the actual proper core.

No Groebner basis, decomposition search, source solve, point or prime is used.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
ID = 'AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    resource.setrlimit(resource.RLIMIT_FSIZE, (32*1024*1024, 32*1024*1024))
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    hashes = {}

    def read(p):
        hashes[str(p.relative_to(ROOT))] = sha(p)
        return json.loads(p.read_text())

    def gate(p):
        record = read(p)
        assert record['all_passed']
        for name, digest in record.get('input_sha256', {}).items():
            assert sha(ROOT/name) == digest, name
            hashes[name] = digest
        for name, digest in record.get('evidence_sha256', {}).items():
            path = p.parent/name
            assert sha(path) == digest, name
            hashes[str(path.relative_to(ROOT))] = digest
        return record

    b = ROOT/'reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63'
    gate(b.with_name(b.name+'_REVIEW')/'report.json')
    primary = read(b/'report.json')
    assert primary['all_passed']
    prep = read(b/'preparation.json')
    old = read(ROOT/prep['source'])
    assert sha(ROOT/prep['source']) == prep['source_sha256']
    a, d, x, y = variables = sp.symbols('a d v2 v6')
    t, u, v, w = sp.symbols('t u v w')
    P, D = d*d+5, 1-d*d
    r, s = x-a-d, a+x

    def decode(rows):
        return sp.Poly.from_dict({tuple(m):sp.Rational(n,q) for m,n,q in rows}, variables, domain=sp.QQ).as_expr()

    def enc(expr, vs=variables):
        return [[list(m), str(c)] for m,c in sp.Poly(sp.expand(expr), *vs, domain=sp.QQ).terms() if c]

    def exact(left, right):
        assert sp.Poly(sp.expand(left-right), *variables, domain=sp.QQ).is_zero

    Q0, Q1, f = [decode(q) for q in prep['generators3']]
    K = [decode(rec['actual_core_polynomial']) for rec in old['remaining7_actual_core_polynomials']]
    assert len(enc(Q0)) == 11 and len(enc(Q1)) == 27
    exact(K[0], D*Q0)
    exact(K[1], D*Q1)
    exact(f, s*r*r)
    q2, rem = sp.div(-K[2], D*r, *variables)
    assert rem == 0
    exact(Q0-2*q2, -P*s*r)
    assert prep['initial_regular_sequence_2_2_3_proved']

    preflight = dict(source_terms=[11,27,8], generic_models=5,
        maximum_IFT_unknowns=2, IFT_order=2, largest_local_initial_space=2,
        only_2x2_or_3x3_fixed_determinants=True,
        no_GB_decomposition_search_new_source_point_prime_or_higher_source_jet=True,
        seconds_cap=60, address_bytes_cap=2_000_000_000,
        output_bytes_cap=32*1024*1024, expanded_terms_cap=150000)
    (out/'preflight.json').write_text(json.dumps(preflight, indent=2)+'\n')

    # Exhaustion by elementary prime cases, keeping the only formal units D,P.
    qr0 = sp.expand(Q0.subs(x,a+d))
    qr1 = sp.expand(Q1.subs(x,a+d))
    h = -a*a-3*a*d+d*y
    j = a*a*d+2*a*d*d-5*a+5*y
    exact(qr0, 8*h)
    exact(qr1, -8*a*j)
    exact(d*j-5*h, a*P*(a+2*d))
    exact(d*qr1+5*a*qr0, -8*a*a*P*(a+2*d))
    qs0 = sp.expand(Q0.subs(x,-a))
    qs1 = sp.expand(Q1.subs(x,-a))
    hs = -a+2*d+2*y
    exact(qs0, 4*d*hs)
    qsquot, qsrem = sp.div(qs1, hs, a, d, y)
    assert sp.expand(qsrem-2*d*(d+2*y)*P) == 0
    exact(qs1, qsquot*hs+2*d*(d+2*y)*P)
    exact(qs1.subs(d,0), 10*a*(a-2*y))

    lines = {
      'N': {a:0,d:0,x:0,y:t},
      'C': {a:0,d:t,x:t,y:0},
      'R': {a:-2*t,d:t,x:-t,y:-2*t},
      'P': {a:2*t,d:0,x:-2*t,y:t},
      'L': {a:t,d:t,x:-t,y:-t/2},
    }
    line_primes = {
      'N':[a,d,x], 'C':[a,y,x-d],
      'R':[a+2*d,x+d,y+2*d],
      'P':[d,a-2*y,x+2*y], 'L':[a-d,x+d,2*y+d],
    }
    for line in lines.values():
        assert all(sp.expand(q.subs(line)) == 0 for q in [Q0,Q1,f])
    support = dict(
        r_restriction_Q0=enc(qr0), r_restriction_Q1=enc(qr1),
        r_elimination_h=enc(h), r_elimination_j=enc(j),
        r_literal_elimination=enc(-8*a*a*P*(a+2*d)),
        s_restriction_Q0=enc(qs0), s_restriction_Q1=enc(qs1),
        s_linear_factor=enc(hs), s_division_quotient=enc(qsquot),
        s_division_remainder=enc(qsrem), s_d_zero=enc(qs1.subs(d,0)),
        primes={k:[enc(z) for z in values] for k,values in line_primes.items()},
        all5_literal_line_substitutions_zero=True,
        proof='Any prime over (Q0,Q1,sr^2) contains r or s. At r=0: d=0 forces a=0; d notin prime lets the elimination force a=0 or a+2d=0, giving C or R. At s=0: d=0 gives a=0 or a-2v6=0; d notin prime forces -a+2d+2v6=0 and d+2v6=0, giving L. P=d^2+5 is a unit. The converse is the literal five-line substitution.')

    field = sp.QQ.frac_field(t)

    def polynomial(expr, vs):
        return sp.Poly(sp.expand(expr), *vs, domain=field)

    def rat(expr):
        return str(sp.factor(sp.cancel(expr)))

    def det_fixed(matrix):
        n = len(matrix)
        ans = 0
        for perm in itertools.permutations(range(n)):
            inv = sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))
            ans += (-1)**inv*sp.prod(matrix[i][perm[i]] for i in range(n))
        return sp.factor(ans)

    def solve2(A, rhs):
        det = det_fixed(A)
        assert det != 0
        sol = [sp.cancel((rhs[0]*A[1][1]-A[0][1]*rhs[1])/det),
               sp.cancel((A[0][0]*rhs[1]-rhs[0]*A[1][0])/det)]
        assert all(sp.cancel(sum(A[i][j]*sol[j] for j in range(2))-rhs[i]) == 0 for i in range(2))
        return sol

    def twojet(label, free, eliminated):
        center = lines[label]
        A = [[sp.diff(q,z).subs(center) for z in eliminated] for q in [Q0,Q1]]
        base = dict(center)
        base[free] = center[free]+u
        at_free = [polynomial(q.subs(base), (u,)) for q in [Q0,Q1]]
        slopes = solve2(A, [-p.nth(1) for p in at_free])
        first = dict(base)
        for z,c in zip(eliminated,slopes): first[z] = center[z]+c*u
        p2 = [polynomial(q.subs(first), (u,)) for q in [Q0,Q1]]
        corrections = solve2(A, [-p.nth(2) for p in p2])
        jets = dict(first)
        for z,c,e in zip(eliminated,slopes,corrections): jets[z] = center[z]+c*u+e*u*u
        values = [polynomial(q.subs(jets), (u,)) for q in [Q0,Q1]]
        assert all(p.nth(k) == 0 for p in values for k in range(3))
        ff = polynomial(f.subs(jets), (u,))
        return dict(center={str(z):rat(center[z]) for z in variables},
            free_variable=str(free), eliminated=list(map(str,eliminated)),
            IFT_Jacobian=[[rat(c) for c in row] for row in A],
            IFT_determinant=rat(det_fixed(A)),
            first_coefficients=list(map(rat,slopes)), second_coefficients=list(map(rat,corrections)),
            all4_coordinate_jets={str(z):rat(jets[z]) for z in variables},
            Q0_Q1_all_coefficients_through2_zero=True,
            f_first5_coefficients=[rat(ff.nth(k)) for k in range(5)]), jets, ff

    N, njets, nf = twojet('N', a, [d,x])
    assert sp.cancel(njets[d]-4*u*u/t) == 0
    assert sp.cancel(njets[x]+3*u-2*u*u/t) == 0
    assert all(nf.nth(k)==0 for k in range(3)) and nf.nth(3)==-32
    N.update(length=3, embedding_dimension=1, HF=[1,1,1],
        Artin_model='QQ(t)[[u]]/(u^3)', nilpotence_index=3,
        d_order=2, d_leading=rat(4/t), d_squared_zero=True)

    # C has a genuine quadratic, two-normal local complete intersection.
    csub={a:u,d:t,x:t+u+v}
    cQ0=sp.expand(Q0.subs(csub))
    assert sp.diff(cQ0,y)==8*t
    croot=sp.cancel(-cQ0.subs(y,0)/(8*t))
    assert sp.expand(cQ0-8*t*(y-croot)) == 0
    cQ1=polynomial(Q1.subs(csub).subs(y,croot), (u,v))
    assert all(cQ1.coeff_monomial(m)==0 for m in [1,u,v])
    ch2=sum(c*u**i*v**j for (i,j),c in cQ1.terms() if i+j==2)
    assert polynomial(ch2,(u,v)).coeff_monomial(u*u)==-16*(t*t+5)
    assert sp.cancel(cQ1.as_expr().subs(v,0)+8*(t*t+5)*u*u*(u+2*t)/t)==0
    C=dict(center={str(z):rat(lines['C'][z]) for z in variables},
        normal_coordinates=['u=a','v=r=v2-a-d'], Q0_y_unit=rat(8*t),
        exact_y_root=rat(croot), exact_Q1_after_y=[(list(m),rat(c)) for m,c in cQ1.terms()],
        quadratic_initial=rat(ch2), coefficient_u_squared=rat(-16*(t*t+5)),
        other_relation='v^2 (since s=t+2u+v is a unit)',
        initial_regular_sequence_reason='The first quadratic has nonzero u^2 coefficient and is coprime to v^2.',
        length=4, embedding_dimension=2, HF=[1,2,1], nilpotence_index=3,
        r_divisor_length=2, s_is_unit=True, d_is_unit=True)

    def normal_jac(label, last):
        normal = [a,x,y]
        jac = [[sp.factor(sp.diff(q,z).subs(lines[label])) for z in normal] for q in [Q0,Q1,last]]
        det = det_fixed(jac)
        assert det != 0
        return dict(center={str(z):rat(lines[label][z]) for z in variables},
            parameter='d=t', normal_variables=list(map(str,normal)),
            exact_normal_Jacobian=[[rat(c) for c in row] for row in jac],
            determinant=rat(det), source_normal_coordinates=['Q0','Q1',str(last)])

    R=normal_jac('R',r)
    assert sp.expand(s.subs(lines['R'])) == -3*t
    R.update(length=2,embedding_dimension=1,HF=[1,1],nilpotence_index=2,
        Artin_model='QQ(t)[[r]]/(r^2)',s_is_unit=True,d_is_unit=True)
    PP,pjets,pf=twojet('P',d,[a,x])
    assert sp.cancel((pjets[a]+pjets[x]).coeff(u,1)) == 0
    assert polynomial(pjets[a]+pjets[x],(u,)).nth(2) == -1/(3*t)
    assert pf.nth(0)==pf.nth(1)==0 and sp.cancel(pf.nth(2)+16*t/3)==0
    PP.update(length=2,embedding_dimension=1,HF=[1,1],nilpotence_index=2,
        Artin_model='QQ(t)[[d]]/(d^2)',d_order=1,d_nonzero=True,d_squared_zero=True,
        s_order=2,s_leading=rat(-1/(3*t)),r_is_unit=True)
    L=normal_jac('L',s)
    assert sp.expand(r.subs(lines['L'])) == -3*t
    L.update(length=1,embedding_dimension=0,HF=[1],nilpotence_index=1,
        Artin_model='QQ(t)',r_is_unit=True,d_is_unit=True)
    models=dict(N=N,C=C,R=R,P=PP,L=L)
    assert [models[k]['length'] for k in ['N','C','R','P','L']] == [3,4,2,2,1]
    assert sum(m['length'] for m in models.values()) == 12

    # Explicit finite checks suffice for the filtration/associated-prime proof.
    core=dict(reduced_support_exact_five_lines=True, geometric_generic_lengths=[3,4,2,2,1],
        generic_HFs={k:v['HF'] for k,v in models.items()},
        divisor_lengths_r=dict(N=1,C=2,R=1),divisor_lengths_s=dict(N=1,P=2,L=1),
        associativity_multiplicity_sum=12, CM_and_no_embedded_primes_from_accepted_actual_CI=True,
        nilradical_index_exact=3,
        nilradical_proof='The one-dimensional CM actual core has only its five minimal associated primes. It injects into the product of those generic localizations. Their maximal ideals have nilpotence indices 3,3,2,2,1; thus Nil^3=0 globally, while Nil^2 is nonzero at N and C.')
    whole=dict(ring='QQ[[u1..u5,a,d,v2,v6]] / d*(u1..u5,Q0,Q1,s*r^2)',
        Ass=['(d)']+[k for k in ['N','C','R','P','L']],
        minimal_Ass=['(d)','C','R','L'], embedded_Ass=['N','P'],
        Ass_proof='The exact injection core --times d--> E4 and quotient E4/(d)=S/(d) bound Ass(E4) by these six primes and put all five core primes in Ass(E4); (d) is minimal and hence associated.',
        at_N_d_in_core='4*u^2/t mod u^3',at_P_d_in_core='d mod d^2',
        nilradical_index_exact=3,
        nilradical_proof='At N and P the only minimal prime after localization is (d), so Nil(E4)=(d)/(dJ). Via multiplication by d its square is d*core and its cube is d^2*core=0; d*core is nonzero by the exact generic models. At C,R,L the element d is a unit and E4 is the corresponding core algebra, with nilpotence at most3. At (d) the ring is generically reduced. A submodule vanishing at every associated prime is zero, proving Nil(E4)^3=0 globally, with square nonzero.',
        full_dimension=8,depth=1,
        no_original_D0_D_depth_nilradical_or_actual_Poisson_component_identification=True)

    evidence=dict(source_Q0=enc(Q0),source_Q1=enc(Q1),source_f=enc(f),
        source_K0_D_Q0_K1_D_Q1_K2_relation_checked=True,
        support=support,generic_models=models,actual_core=core,normalized_E4=whole,
        no_CAS_GB_decomposition_search_numeric_sourcepoint_prime_or_new_source_jet=True)
    terms=sum(len(sp.Poly(q,*variables).terms()) for q in [Q0,Q1,f,qr0,qr1,qs0,qs1,h,j,hs,qsquot,qsrem])
    assert terms < 150000
    (out/'five-line-and-generic-models.json').write_text(json.dumps(evidence,indent=2)+'\n')
    proofnote=ROOT/f'research/notes/{ID}.md'
    for p in [Path(__file__),proofnote]: hashes[str(p.relative_to(ROOT))]=sha(p)
    report=dict(id=ID,all_passed=True,source_bound_actual_full_core=True,
        exact_reduced_core_support_five_lines=True,core_generic_lengths=[3,4,2,2,1],
        actual_core_nilradical_index_exact3=True,normalized_E4_nilradical_index_exact3=True,
        normalized_E4_associated_primes_six_embedded_N_P=True,
        no_original_D0_or_full_D_nil_depth_claim=True,
        no_GB_new_point_prime_sourcejet_or_decomposition_search=True,
        preflight_expanded_identity_terms=terms,input_sha256=hashes,
        evidence_sha256={p.name:sha(p) for p in out.iterdir()},
        seconds=time.monotonic()-started,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        frozen=True,independent_review_pending=True)
    path=out/'report.json'
    path.write_text(json.dumps(report,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir()) < 32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],report_sha256=sha(path))))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)

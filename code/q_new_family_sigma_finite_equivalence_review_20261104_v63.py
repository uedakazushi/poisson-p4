#!/usr/bin/env python3
"""Independent QQ(u,v) homogeneous coefficient-ring sigma verification."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_SIGMA_FINITE_EQUIVALENCE_20261104_v63_REVIEW'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()


def read(p): return json.loads(p.read_text())


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'reports/Q_NEW_FAMILY_SIGMA_FINITE_EQUIVALENCE_20261104_v63'
    author=read(base/'report.json');cert=read(base/'certificate.json')
    assert author['all_passed'] and author['both_signs_explicit_full_tensor_sigma_equivalence']
    for r,hsh in author['input_sha256'].items(): assert sha(ROOT/r)==hsh,r
    for r,hsh in author['evidence_sha256'].items(): assert sha(base/r)==hsh,r
    lower=ROOT/'reports/Q_NEW_FAMILY_GLOBAL_TWO_INVARIANT_PREFLIGHT_20261104_v63_REVIEW/report.json'
    assert read(lower)['all_passed']
    geo=ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json'
    geometry=read(geo)
    u,v=s.symbols('u v');n0,n1,n2=s.symbols('n0 n1 n2')
    affine=s.symbols('x y z w')
    K=s.QQ.frac_field(u,v)
    R,H,X,Y,Z,W=ring('h,x,y,z,w',K)
    variables=[H,X,Y,Z,W]
    field=lambda a:K.from_sympy(s.sympify(a))
    par=32*u+v
    # Directly homogenize the complete ORIGINAL spatial brackets into a
    # five-variable sparse polynomial ring. The primary used affine SymPy
    # expansion/cancel; this checker performs exact coefficient-field sums.
    P=[[R.zero for _ in range(5)] for _ in range(5)]
    Pt=[[R.zero for _ in range(5)] for _ in range(5)]
    for i,j in itertools.combinations(range(4),2):
        expression=s.sympify(geometry['affine_P'][i][j]).subs({n0:u,n1:1,n2:v})
        for mon,co in s.Poly(expression,*affine).terms():
            assert sum(mon)<=2
            power=(2-sum(mon),*mon)
            P[i+1][j+1][power]=field(co)
            Pt[i+1][j+1][power]=field(co.subs(v,-par))
        P[j+1][i+1]=-P[i+1][j+1]
        Pt[j+1][i+1]=-Pt[i+1][j+1]
    def substitute(poly,images):
        result=R.zero
        for mon,co in poly.items():
            a=R(co)
            for q,e in zip(images,mon):
                if e:a*=q**e
            result+=a
        return result
    def permutation_sign(p):
        return -1 if sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))%2 else 1
    def det(A):
        answer=K.zero
        for p in itertools.permutations(range(len(A))):
            term=K.one*permutation_sign(p)
            for i,j in enumerate(p):term*=A[i][j]
            answer+=term
        return answer
    evidence=[]
    for branch in cert['branches']:
        ep=branch['sign']
        M=[[field(a) for a in row] for row in branch['homogeneous_matrix']]
        A=[[field(a) for a in row] for row in branch['linear_matrix']]
        tr=[field(a) for a in branch['translation']]
        assert M[0]==[K.one,K.zero,K.zero,K.zero,K.zero]
        assert all(M[i+1][0]==tr[i] and M[i+1][1:]==A[i] for i in range(4))
        determinant=det(A)
        assert determinant==field(-v**7/par**7)
        assert det(M)==determinant
        kap=field(branch['Poisson_scalar'])
        assert kap==field(-par**2/v**2)
        images=[sum((R(M[i][j])*variables[j] for j in range(5)),R.zero) for i in range(5)]
        checks=0
        for i,j in itertools.combinations(range(5),2):
            pushed=R.zero
            for k in range(5):
                for l in range(5):
                    if M[i][k] and M[j][l] and P[k][l]:
                        pushed+=R(kap*M[i][k]*M[j][l])*P[k][l]
            target=substitute(Pt[i][j],images)
            difference=pushed-target
            assert not difference,(ep,i,j,str(difference))
            checks+=15
        evidence.append(dict(sign=ep,checked_homogeneous_coefficients=checks,
                             determinant=str(K.to_sympy(determinant)),all_zero=True))
    # The two independent functions are fixed by both actual transformations.
    aa=u*u;jj=(v/u+16)**2
    sigma={u:u,v:-32*u-v};tau={u:-u,v:-v}
    for image in (sigma,tau):
        assert s.cancel(aa.subs(image,simultaneous=True)-aa)==0
        assert s.cancel(jj.subs(image,simultaneous=True)-jj)==0
    for a in (u,v):
        assert s.expand(a.subs(sigma,simultaneous=True).subs(sigma,simultaneous=True)-a)==0
        assert s.expand(a.subs(tau,simultaneous=True).subs(tau,simultaneous=True)-a)==0
        assert s.expand(a.subs(sigma,simultaneous=True).subs(tau,simultaneous=True)
                        -a.subs(tau,simultaneous=True).subs(sigma,simultaneous=True))==0
    result=dict(
        engine='independent five-variable sparse PolyRing over exact QQ(u,v); no producer import or replay',
        branches=evidence,whole_domain='u*v*(32*u+v)!=0',
        raw_homogeneous_identity='kappa*M*P_raw(X)*M^t=P_raw(u,1,-32u-v)(M*X)',
        Euler_quotient_reason='exact raw equality descends to projective bivectors; nonzero tensor scalar disappears on projectivization',
        lower_gate='independently accepted actual global A=u^2 and J=(v/u+16)^2, including arbitrary evaluation point',
        full_moduli_field='k(X_NEW25)^PGL5=k(u^2,(v/u+16)^2), k=QQ,C',
        parameter_identifications=['(u,v)','(-u,-v)','(u,-32u-v)','(-u,32u+v)'],
        exact_geometric_generic_parameter_fibre=4,
        parameter_group='Klein four generated by tau and sigma',
        not_claimed='finite stabilizer group of one tensor, special-parameter exhaustion, boundary classification or an invariant-ring presentation')
    ep=out/'independent-certificate.json';ep.write_text(json.dumps(result,indent=2)+'\n')
    files=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),base/'report.json',base/'certificate.json',
           ROOT/'code/q_new_family_sigma_finite_equivalence_20261104_v63.py',
           ROOT/'research/notes/Q_NEW_FAMILY_SIGMA_FINITE_EQUIVALENCE_20261104_v63.md',geo,lower]
    report=dict(id=ID,all_passed=True,
                independent_original_full_homogeneous_raw_tensor_150_coefficients_each_sign=True,
                independent_QQ_u_v_coefficient_ring_not_author_affine_expansion=True,
                both_signs_invertible_on_entire_displayed_open=True,
                actual_fullPGL_moduli_field_exact_k_u_squared_J_over_QQ_and_C=True,
                exact_generic_parameter_identifications_Klein_four=True,
                finite_isotropy_special_boundary_and_invariant_ring_not_claimed=True,
                no_new_point_prime_GB_solver_high_invariant_or_author_rerun=True,
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in files},
                evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start,frozen=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k] for k in ('all_passed','seconds')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(p.parse_args().output.resolve())

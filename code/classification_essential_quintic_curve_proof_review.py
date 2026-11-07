#!/usr/bin/env python3
"""Finite polynomial identities for the separate eight-curve necessity proof.
No CAS/producer imports. Geometry and all open complements are in the review.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,resource,signal,time
from classification_minimal_quartic_review import LaurentRing
ROOT=Path(__file__).resolve().parents[1]
JOB='research/classification/jobs/essential_three_quintic_curve_proof_review.json'
OUT='reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_CURVE_PROOF_REVIEW_20261107.json'
def main():
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,)*2)
    j=json.loads((ROOT/JOB).read_text());assert not(ROOT/OUT).exists()
    for b in j['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    # All A12 five-parameter odd-vanishing conditions after B=2C,D=2A-C^3,E=-AC^2.
    R=LaurentRing(('A','C'),('r',));p=R.parse
    cut=lambda P:{e:c for e,c in P.items()if e[-1]<14}
    mul=lambda *Ps:cut(R.mul(*Ps))
    den=p('1+2*C*r+(2*A-C**3)*r**3-A*C**2*r**5');v=R.add(p('1'),R.scale(den,-1));inv=p('1');term=p('1')
    for k in range(1,14):term=mul(term,v);inv=R.add(inv,term)
    assert cut(R.mul(den,inv))==p('1')
    x=mul(p('r**4'),inv);y=mul(p('r**2*(1+C*r+A*r**3)'),inv)
    residual=R.add(x,R.scale(mul(y,y),-1),mul(p('C**2'),y,y,y),R.scale(mul(p('C**4-2*A*C'),y,y,y,y),-1),R.scale(mul(p('-A**2+6*A*C**3-C**6'),y,y,y,y,y),-1))
    coeff=lambda P,k:{e[:-1]+(0,):c for e,c in P.items()if e[-1]==k}
    assert all(not coeff(residual,k)for k in (5,7,9))
    assert coeff(residual,11)==p('-A**2*C')
    # Bicuspidal exact necessary conditions, all A!=0 and remaining charts.
    S=LaurentRing(('A',),(),invert=('A',));p=S.parse
    C=p('-2-1/A');K=S.add(C,p('-A'));D=p('(A+1)**2*(A-1)/A**3')
    assert S.add(S.mul(p('A'),C),p('2*A+1'))=={}
    assert S.add(S.mul(K,S.add(C,p('3'))),S.mul(p('A'),D))=={}
    L2=S.add(S.mul(D,S.add(S.mul(p('A'),C),p('1'))),S.scale(K,-3))
    L3=S.add(S.mul(D,S.add(S.scale(C,2),p('-A'))),K)
    assert S.add(L3,S.mul(S.add(C,p('4')),L2))==p('2*(A+1)**2/A**4')
    C=p('-1/A');D=p('-1/A**3');assert S.add(S.mul(C,C),S.mul(p('A'),D))=={}
    assert S.mul(p('2'),C,D)==p('2/A**4')
    # Both Cremona maps have explicit birational inverses on their stated open.
    T=LaurentRing((),('x','y','z'));p=T.parse;x,y,z=[p(n)for n in ('x','y','z')]
    X,Y,Z=p('x**2'),p('x*y'),p('y**2-x*z')
    assert [T.mul(X,X),T.mul(X,Y),T.add(T.mul(Y,Y),T.scale(T.mul(X,Z),-1))]==[p('x**4'),p('x**3*y'),p('x**3*z')]
    X,Y,Z=p('x*y'),p('y**2'),p('x*z')
    assert [T.mul(X,X),T.mul(X,Y),T.mul(Y,Z)]==[p('x**2*y**2'),p('x*y**3'),p('x*y**2*z')]
    # Unique A6 quartic and its actual dual sections; no quoted list required.
    B=LaurentRing((),('s','t'));p=B.parse
    X,Y,Z=p('s**4-s*t**3'),p('s**2*t**2'),p('t**4');H=B.add(B.mul(Y,Y),B.scale(B.mul(X,Z),-1))
    assert B.add(B.mul(H,H),B.scale(B.mul(Y,Z,Z,Z),-1))=={}
    gradient=[B.scale(B.mul(Z,H),-2),B.add(B.scale(B.mul(Y,H),4),B.scale(B.mul(Z,Z,Z),-1)),B.add(B.scale(B.mul(X,H),-2),B.scale(B.mul(Y,Z,Z),-3))]
    dual=[p('-2*s*t**4'),p('t**2*(4*s**3-t**3)'),p('-s**2*(2*s**3+t**3)')]
    assert gradient==[B.mul(p('t**7'),a)for a in dual]
    # All eight explicit cusp/source certificates were independently established
    # before the stronger necessity proof; do not use necessity circularly.
    old=json.loads((ROOT/'reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_EIGHT_SOURCE_REVIEW_20261107.json').read_text());assert old['all_passed']
    records=[dict(id=r['id'],signature=r['intrinsic_cusp_recognition']['signature'])for r in old['all8_explicit_sources']]
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,
        all8_necessary_curve_classes_pureproof_independently_reviewed=True,
        A12_general_odd5_7_9_11_conditions=True,bicuspidal_all_generic_and_tangent_charts=True,
        both_Cremona_actual_inverses=True,A6_quartic_and_actual_dual_identity=True,
        independent_explicit_source_receipt_used_without_necessity_circularity=records,
        geometry_scope='Full proof including all source-section charts/basepoints/RH/Plucker/dual genus and projective normalizations is in ESSENTIAL_THREE_QUINTIC_CURVE_PROOF_SEPARATE_REVIEW.md. This finite arithmetic alone is not a curve-classification proof.',
        no_Namba_or_unread_primary_classification_dependency=True,no_CAS_or_producer_imports=True,
        inputs=j['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (ROOT/OUT).write_text(json.dumps(result,indent=2)+'\n');print('PASS',OUT,round(result['seconds'],3),flush=True)
if __name__=='__main__':main()

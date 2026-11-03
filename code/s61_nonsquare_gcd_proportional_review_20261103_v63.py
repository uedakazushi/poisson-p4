#!/usr/bin/env python3
"""Independent exact-QQ proportionality exclusion on a*P*J!=0; no GB."""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import tarfile
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_GCD_GATE_PROPORTIONAL_EXCLUSION_20261103_v63_REVIEW'


def read(path):return json.loads(path.read_text())
def sha(path):
    h=hashlib.sha256()
    with path.open('rb')as f:
        for c in iter(lambda:f.read(1048576),b''):h.update(c)
    return h.hexdigest()
def decode(rows):return{tuple(ex):Q(c)for ex,c in rows if Q(c)}
def mul(a,b):
    out={}
    for ex,c in a.items():
        for ey,d in b.items():
            mon=tuple(x+y for x,y in zip(ex,ey));out[mon]=out.get(mon,Q(0))+c*d
    return{mon:c for mon,c in out.items()if c}
def add(a,b):
    out=a.copy()
    for mon,c in b.items():out[mon]=out.get(mon,Q(0))+c
    return{mon:c for mon,c in out.items()if c}


def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(120)
    parent=ROOT/'reports/S61_NONSQUARE_GCD_GATE_PROPORTIONAL_EXCLUSION_20261103_v63';report=read(parent/'report.json');assert report['all_exact_checks_passed']and report['frozen']
    archive=ROOT/report['archive'];assert sha(archive)==report['archive_sha256'];cert=None;memberhashes={}
    with tarfile.open(archive)as tar:
        for m in tar:
            assert m.isfile();h=hashlib.sha256()
            with tar.extractfile(m)as file:
                content=file.read();h.update(content)
            memberhashes[m.name]=h.hexdigest()
            if m.name.endswith('/certificate.json'):cert=json.loads(content)
    assert cert is not None
    for rel,h in report['input_sha256'].items():assert (memberhashes[rel]if rel in memberhashes else sha(ROOT/rel))==h
    source=ROOT/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json';data=read(source);assert data['all_passed']and sha(source)=='80f7c46b1bd4ffef4a78f492ad4f7a2bb40b2ffa27cddcd23d5b5d660db4a7e0'
    a,e,P=s.symbols('a e P');xx=s.symbols('x0:6');r,u,v,t,R=s.symbols('r u v t R');ys=(r,u,v,t,R)
    quads=[];terms=0
    for i,row in enumerate(data['cleared_qfour_sparse']):
        q=0
        for mon,coeffs in row:
            assert len(mon)==7 and mon[6]==0 and sum(mon)==2
            param=sum(s.Rational(c)*a**(ex[0]-ex[2]+[0,2,2,1][i])*e**ex[1]*P**ex[2]for ex,c in coeffs);terms+=len(coeffs)
            q+=param*s.prod(x**power for x,power in zip(xx,mon[:6]))
        quads.append(q)
    assert terms==804
    J=5*P-6*e-4;Lp=20*P+3*e+20;Rc=15*P**2-6*P*e+6*P+80*a;Tc=-5*P+6*e-2;Sc=21*P*e+16*P+140*a;Vc=9*P*e+4*P+60*a;Kc=85*P**2-12*P*e+58*P+600*a
    B=30*a*Rc*xx[0]-9*P*Sc*xx[1]+90*P*a*Tc*xx[2]-18*P*Vc*xx[3]-2*a*Kc*xx[4];Af=(B+10*a*P*J*xx[4])/3
    assert s.expand(quads[1].subs(xx[5],0)-Af*B/6075)==0
    assert s.cancel(((quads[2]-e*quads[1])/(a*a*J)).subs(xx[5],0)-40*xx[0]*B/6075)==0
    assert s.expand(3*Sc-7*Vc-20*P)==0
    inverse={xx[0]:r,xx[1]:v,xx[2]:R/a-u-6*(e+1)*v/(5*a)-(9*e+4)*t/(10*a),xx[3]:(t-v)/2,xx[4]:5*u+5*r-4*R/a+24*(e+1)*v/(5*a)+2*(9*e+4)*t/(5*a),xx[5]:0}
    change=s.Matrix([inverse[x]for x in xx[:5]]).jacobian(ys);change_det=s.factor(change.det());assert s.Poly(s.cancel(change_det*a),a,e,P).total_degree()==0 and change_det!=0
    assert s.cancel(4*quads[0].subs(inverse,simultaneous=True)-u*R+5*v*t)==0
    # Coefficient extraction, not the author's Hessian call.
    q=s.Poly(s.cancel(a*a*quads[3].subs(inverse,simultaneous=True)),*ys)
    H=s.Matrix(5,5,lambda i,j:s.Poly(q.coeff_monomial(ys[i]*ys[j])/(1 if i==j else 2),a,e,P,domain=s.QQ).as_expr())
    literal=lambda terms,variables:sum(s.Rational(c)*s.prod(x**n for x,n in zip(variables,ex))for ex,c in terms)
    savedH=s.Matrix([[literal(q,(a,e,P))for q in row]for row in cert['global_H_a2Q3']]);assert H==savedH
    bn=s.Poly(s.cancel(B.subs(inverse,simultaneous=True)),*ys);bc=[s.Poly(bn.coeff_monomial(x),a,e,P,domain=s.QQ).as_expr()for x in ys]
    assert bc==[literal(q,(a,e,P))for q in cert['global_B_normal_coefficients']]
    DB=20*P**2+3*P*e+20*P+180*a;br=bc[0];assert s.expand(br+20*a*DB)==0
    cleared=H[0,0]*bc[1]**2-2*br*H[0,1]*bc[1]+br**2*H[1,1]
    assert s.Poly(cleared+a**5*br**2*J/12,a,e,P,domain=s.QQ).is_zero
    # DB=0: exact polynomial cross-minor identities, dividing only recorded P,Lp,J units.
    ad=-P*Lp/180;assert s.cancel(H[0,0].subs(a,ad))==0
    hs=[s.Poly(s.expand(H[0,i].subs(a,ad)),P,e,domain=s.QQ)for i in range(1,5)];bs=[s.Poly(s.expand(c.subs(a,ad)),P,e,domain=s.QQ)for c in bc[1:]]
    Gi=[decode(q)for q in cert['normalized_minors']];cofs=[decode(q)for q in cert['Bezout_coefficients']]
    assert len(Gi)==len(cofs)==5
    actualcross=[]
    for i,((left,right),powers,scalar,g)in enumerate(zip(cert['corresponding_cross_minor_indices'],cert['known_unit_powers_P_L_J'],cert['normalization_QQ_scalars'],Gi)):
        assert 0<=left<right<4 and len(powers)==3 and min(powers)>=0
        raw=hs[left]*bs[right]-hs[right]*bs[left];gexpr=sum(s.Rational(c.numerator,c.denominator)*P**ex[0]*e**ex[1]for ex,c in g.items());known=P**powers[0]*Lp**powers[1]*J**powers[2]
        assert raw==s.Poly(known*gexpr/s.Rational(scalar),P,e,domain=s.QQ);actualcross.append(str(raw.as_expr()))
    total={}
    for c,g in zip(cofs,Gi):total=add(total,mul(c,g))
    assert total=={(0,0):Q(1)}
    finite=ROOT/'reports/S61_NONSQUARE_FINITE_GCD_EXCEPTION_NAKAYAMA_20261103_v63_REVIEW/report.json';finiteproof=read(finite);assert finiteproof['all_passed']and finiteproof['finite_V_cokernel_Nakayama_QQ_surjectivity_on_entire_finite_Sigma_reviewed']
    evidence=dict(split_coordinate_determinant=str(change_det),global_split_Q0='4Q0=u*R-5*v*t',radical_coefficient=str(br),
        DB_nonzero_cleared_U2_identity_verified=True,unscaled_DB_nonzero_Q3_U2_coefficient='-a^3*J/12',
        exact_global_H_a2Q3=[[str(c)for c in row]for row in H.tolist()],exact_B_normal_coefficients=list(map(str,bc)),
        DB_zero_cross_minors=actualcross,only_stripped_units=['P','20P+3e+20','5P-6e-4'],
        DB_zero_Lp_equals_minus180a_overP=True,known_unit_powers=cert['known_unit_powers_P_L_J'],five_polynomial_literal_QQ_Bezout_equals1=True)
    (out/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    result=dict(id=STEM,all_passed=True,own_actual804_quadric_coefficients_rebuilt=True,actual_A_B_Q2_factorization_verified=True,
        split_rank4_coordinate_change_determinant=str(change_det),restricted_Q0_rank_at_least2_globally=True,
        DB_nonzero_actual_unit_coefficient_minus_a3_J_over12_verified=True,DB_zero_all_actual_cross_minor_unit_factorizations_and_QQ_Bezout_verified=True,
        no_proportionality_on_a_P_J_open=True,FULL_U_highrank_gcd_failure_excluded=True,
        all_gcd_gate_failure_points_on_FULL_U_raw_good_by_finite_Sigma_review=True,
        pair_or_five_gate_remaining_and_no_whole_product41_claim=True,seconds=time.monotonic()-start)
    paths=[Path(__file__),source,parent/'report.json',archive,ROOT/'code/s61_nonsquare_gcd_gate_proportional_exclusion_20261103_v63.py',finite];note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['evidence_sha256']={'independent-certificate.json':sha(out/'independent-certificate.json')};(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

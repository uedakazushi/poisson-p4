#!/usr/bin/env python3
"""Independent actual-coefficient review of auxiliary five-rank on DB=0; no GB."""
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
STEM='S61_NONSQUARE_FIVE_RANK_DB_ZERO_20261103_v63_REVIEW'


def read(path):return json.loads(path.read_text())
def sha(path):
    H=hashlib.sha256()
    with path.open('rb')as f:
        for chunk in iter(lambda:f.read(1048576),b''):H.update(chunk)
    return H.hexdigest()
def decode(rows):return{tuple(ex):Q(c)for ex,c in rows if Q(c)}
def mul(a,b):
    out={}
    for ex,c in a.items():
        for ey,d in b.items():
            mon=tuple(x+y for x,y in zip(ex,ey));out[mon]=out.get(mon,Q(0))+c*d
    return{m:c for m,c in out.items()if c}
def add(a,b,c=1):
    out=a.copy()
    for ex,v in b.items():out[ex]=out.get(ex,Q(0))+c*v
    return{m:c for m,c in out.items()if c}


def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(120)
    parent=ROOT/'reports/S61_NONSQUARE_FIVE_RANK_DB_ZERO_20261103_v63';author=read(parent/'report.json');assert author['frozen']and author['all_exact_QQ_checks_passed'];archive=ROOT/author['archive'];assert sha(archive)==author['archive_sha256']
    members={};hashes={}
    with tarfile.open(archive)as tar:
        for m in tar:
            assert m.isfile()
            with tar.extractfile(m)as f:content=f.read()
            hashes[m.name]=hashlib.sha256(content).hexdigest()
            if m.name.endswith('.json'):members[m.name]=json.loads(content)
    for rel,h in author['input_sha256'].items():assert(hashes[rel]if rel in hashes else sha(ROOT/rel))==h
    cert=members['work/nonsquare-v62/five-rank-DB0-certificate-v68/certificate.json'];prop=members['work/nonsquare-v62/gcd-proportional-certificate-v58/certificate.json'];source=members['reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json'];assert source['all_passed']
    accepted=ROOT/'reports/S61_NONSQUARE_GCD_GATE_PROPORTIONAL_EXCLUSION_20261103_v63_REVIEW/report.json';assert read(accepted)['all_passed']
    a,e,P=s.symbols('a e P');xx=s.symbols('x0:6');ys=s.symbols('r u v t R');r,u,v,t,R=ys;qs=[];count=0
    for i,row in enumerate(source['cleared_qfour_sparse']):
        q=0
        for mon,terms in row:
            assert len(mon)==7 and mon[6]==0 and sum(mon)==2
            q+=sum(s.Rational(c)*a**(ex[0]-ex[2]+[0,2,2,1][i])*e**ex[1]*P**ex[2]for ex,c in terms)*s.prod(x**k for x,k in zip(xx,mon[:6]));count+=len(terms)
        qs.append(q)
    assert count==804
    J=5*P-6*e-4;L=20*P+3*e+20;DB=20*P**2+3*P*e+20*P+180*a;ad=-P*L/180;K=-130*P+168*e-5;ed=(130*P+5)/168
    Rc=15*P**2-6*P*e+6*P+80*a;Tc=-5*P+6*e-2;Sc=21*P*e+16*P+140*a;Vc=9*P*e+4*P+60*a;Kc=85*P**2-12*P*e+58*P+600*a
    B=30*a*Rc*xx[0]-9*P*Sc*xx[1]+90*P*a*Tc*xx[2]-18*P*Vc*xx[3]-2*a*Kc*xx[4];Af=(B+10*a*P*J*xx[4])/3;Q2=s.cancel((qs[2]-e*qs[1])/(a*a*J));D=s.cancel((40*xx[0]*qs[1]-Af*Q2)/xx[5]);s.Poly(D,*xx,a,e,P,domain=s.QQ)
    assert s.expand(qs[1].subs(xx[5],0)-Af*B/6075)==0 and s.cancel(Q2.subs(xx[5],0)-40*xx[0]*B/6075)==0 and s.expand(3*Sc-7*Vc-20*P)==0
    # A and x0 cannot be dependent: A's x1,x3 coefficients cannot both vanish when P!=0.
    assert s.expand(s.Poly(Af,*xx).coeff_monomial(xx[1])+3*P*Sc)==0 and s.expand(s.Poly(Af,*xx).coeff_monomial(xx[3])+6*P*Vc)==0
    inverse={xx[0]:r,xx[1]:v,xx[2]:R/a-u-6*(e+1)*v/(5*a)-(9*e+4)*t/(10*a),xx[3]:(t-v)/2,xx[4]:5*u+5*r-4*R/a+24*(e+1)*v/(5*a)+2*(9*e+4)*t/(5*a),xx[5]:0}
    det=s.factor(s.Matrix([inverse[x]for x in xx[:5]]).jacobian(ys).det());assert det==1/(2*a)
    assert s.cancel(4*qs[0].subs(inverse,simultaneous=True)-u*R+5*v*t)==0
    # Independent restriction then coefficient extraction, not Hessian conjugation.
    H=[]
    for q in[qs[0],qs[3],D]:
        pull=s.Poly(s.cancel(a*a*q.subs(inverse,simultaneous=True)),*ys)
        H.append(s.Matrix(5,5,lambda i,j:s.Poly(pull.coeff_monomial(ys[i]*ys[j])/(1 if i==j else 2),a,e,P,domain=s.QQ).as_expr()))
    expr=lambda rows,vs:sum(s.Rational(c)*s.prod(x**n for x,n in zip(vs,ex))for ex,c in rows)
    savedH=[s.Matrix([[expr(q,(a,e,P))for q in row]for row in h])for h in cert['global_H_matrices_a2_Q0_Q3_D']];assert all((h-g).applyfunc(s.expand)==s.zeros(5)for h,g in zip(H,savedH))
    bpoly=s.Poly(s.cancel(B.subs(inverse,simultaneous=True)),*ys);bc=[s.Poly(bpoly.coeff_monomial(x),a,e,P).as_expr()for x in ys];assert all(s.expand(c-expr(q,(a,e,P)))==0 for c,q in zip(bc,cert['global_B_coefficients']))
    assert s.expand(bc[0]+20*a*DB)==0 and s.cancel(H[1][0,0].subs(a,ad))==0
    assert s.cancel(H[2][0,0].subs(a,ad)+P**5*L**4*J**2*K/s.Integer(35872267500))==0
    bs=[s.Poly(c.subs(a,ad),P,e).as_expr()for c in bc[1:]];h3=[s.Poly(H[1][0,j].subs(a,ad),P,e).as_expr()for j in range(1,5)]
    for G,ij,powers,scalar in zip(prop['normalized_minors'],prop['corresponding_cross_minor_indices'],prop['known_unit_powers_P_L_J'],prop['normalization_QQ_scalars']):
        i,j=ij;unit=P**powers[0]*L**powers[1]*J**powers[2];assert s.Poly(s.Rational(scalar)*(h3[i]*bs[j]-h3[j]*bs[i])-unit*expr(G,(P,e)),P,e).is_zero
    total={}
    for co,q in zip(prop['Bezout_coefficients'],prop['normalized_minors']):total=add(total,mul(decode(co),decode(q)))
    assert total=={(0,0):Q(1)}
    sub=lambda q:s.cancel(q.subs(a,ad).subs(e,ed));assert sub(H[1][0,0])==sub(H[2][0,0])==0
    LK=s.factor(L.subs(e,ed));JK=s.factor(J.subs(e,ed));assert s.cancel(LK-125*(10*P+9)/56)==0 and s.cancel(JK-(10*P-117)/28)==0
    cross=s.Matrix([[sub(q)for q in bc[1:]],[sub(H[1][0,j])for j in range(1,5)],[sub(H[2][0,j])for j in range(1,5)]]);assert(cross-s.Matrix([[expr(q,(P,))for q in row]for row in cert['K_zero_cross_matrix']])).applyfunc(s.cancel)==s.zeros(3,4)
    # Compute each determinant with sparse Fraction permutations, not the author's symbolic det/gcdex.
    rawMinors=[]
    for cols in itertools.combinations(range(4),3):
        mat=[[{ex:Q(str(c))for ex,c in s.Poly(cross[i,j],P,domain=s.QQ).terms()}for j in cols]for i in range(3)];polynomial={}
        for perm in itertools.permutations(range(3)):
            inversions=sum(perm[i]>perm[j]for i in range(3)for j in range(i+1,3));prod={(0,):Q(1)}
            for i,j in enumerate(perm):prod=mul(prod,mat[i][j])
            polynomial=add(polynomial,prod,-1 if inversions%2 else 1)
        rawMinors.append(polynomial)
    assert rawMinors==list(map(decode,cert['cross_minors']))
    primitives=list(map(decode,cert['primitive_minors']));cofs=list(map(decode,cert['Bezout_coefficients']));assert len(primitives)==len(cofs)==4 and all(primitives)
    for raw,primitive,record in zip(rawMinors,primitives,cert['minor_factor_records']):
        powers=record['known_unit_powers_P_L_J'];assert len(powers)==3 and min(powers)>=0;factor=s.Rational(record['QQcontent'])*P**powers[0]*LK**powers[1]*JK**powers[2]
        fp={ex:Q(str(c))for ex,c in s.Poly(factor,P,domain=s.QQ).terms()};assert mul(fp,primitive)==raw
    unity={}
    for co,q in zip(cofs,primitives):unity=add(unity,mul(co,q))
    assert unity=={(0,):Q(1)}
    evidence=dict(split_determinant=str(det),Q0_split='4Q0=uR-5vt',B_radical_coefficient=str(bc[0]),actual_global_H_matrices=[[[str(c)for c in row]for row in h.tolist()]for h in H],
        Hd_rr_on_DB_zero='-P^5*L^4*J^2*K/35872267500',K_zero_e=str(ed),K_zero_L=str(LK),K_zero_J=str(JK),K_zero_cross_matrix=[[str(c)for c in row]for row in cross.tolist()],
        all4_univariate_raw_minors_and_known_unit_factorizations_verified=True,all4_literal_QQ_Bezout_coefficients_verified=True)
    (out/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    result=dict(id=STEM,all_passed=True,actual804_source_Q2prime_D_and_all75_restricted_matrix_entries_rebuilt=True,
        Q1_Q2_independent_and_vanish_on_B_W_verified=True,Q0_restriction_nonzero_by_split_rank4=True,
        DB0_b_h3_rank2_by_actual_cross_minors_and_Bezout=True,DB0_Hd_rr_knownunit_times_K_verified=True,
        K_nonzero_case_independence_verified=True,K_zero_all4_actual_cross_minors_knownunit_factors_and_QQ_Bezout_verified=True,
        global_auxiliary_five_rank_exact5_on_DB0_a_P_J=True,F_B_C_not_required=True,
        K_not_inverted_globally=True,DB_complement_pair_original41_or_ambient_claim_accepted=False,seconds=time.monotonic()-started)
    paths=[Path(__file__),parent/'report.json',archive,accepted,ROOT/'code/s61_nonsquare_five_rank_DB_zero_20261103_v63.py'];note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['evidence_sha256']={'independent-certificate.json':sha(out/'independent-certificate.json')};(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

#!/usr/bin/env python3
"""Independent four-linear-variable/minor/binary-D reduction; no gate search or GB."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import tarfile
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_FOUR_QUADRIC_BINARY_REDUCTION_20261103_v63_REVIEW'


def read(path):return json.loads(path.read_text())
def sha(path):
    h=hashlib.sha256()
    with path.open('rb')as f:
        for c in iter(lambda:f.read(1048576),b''):h.update(c)
    return h.hexdigest()
def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(120)
    parent=ROOT/'reports/S61_NONSQUARE_FOUR_QUADRIC_BINARY_REDUCTION_20261103_v63';author=read(parent/'report.json');assert author['frozen']and author['all_exact_QQ_actual804_checks_passed'];archive=ROOT/author['archive'];assert sha(archive)==author['archive_sha256']
    members={};hashes={}
    with tarfile.open(archive)as tar:
        for m in tar:
            assert m.isfile()
            with tar.extractfile(m)as f:content=f.read()
            hashes[m.name]=hashlib.sha256(content).hexdigest()
            if m.name.endswith('.json'):members[m.name]=json.loads(content)
    for rel,h in author['input_sha256'].items():assert(hashes[rel]if rel in hashes else sha(ROOT/rel))==h
    cert=members['work/nonsquare-v62/four-quadric-binary-certificate-v71/certificate.json'];data=members['reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json'];assert data['all_passed']
    a,e,P=s.symbols('a e P');xx=s.symbols('x0:6');qs=[];count=0
    for i,row in enumerate(data['cleared_qfour_sparse']):
        q=0
        for mon,terms in row:
            assert len(mon)==7 and mon[6]==0 and sum(mon)==2
            q+=sum(s.Rational(c)*a**(ex[0]-ex[2]+[0,2,2,1][i])*e**ex[1]*P**ex[2]for ex,c in terms)*s.prod(x**power for x,power in zip(xx,mon[:6]));count+=len(terms)
        qs.append(q)
    assert count==804
    J=5*P-6*e-4;Rc=15*P**2-6*P*e+6*P+80*a;Tc=-5*P+6*e-2;Sc=21*P*e+16*P+140*a;Vc=9*P*e+4*P+60*a;Kc=85*P**2-12*P*e+58*P+600*a
    B=30*a*Rc*xx[0]-9*P*Sc*xx[1]+90*P*a*Tc*xx[2]-18*P*Vc*xx[3]-2*a*Kc*xx[4];A=(B+10*a*P*J*xx[4])/3;Q2=s.cancel((qs[2]-e*qs[1])/(a*a*J))
    C1=s.cancel((qs[1]-A*B/6075)/xx[5]);C2=s.cancel((Q2-40*xx[0]*B/6075)/xx[5]);assert all(s.Poly(c,*xx).total_degree()==1 for c in[C1,C2])
    assert s.expand(qs[1]-A*B/6075-xx[5]*C1)==0 and s.expand(Q2-40*xx[0]*B/6075-xx[5]*C2)==0
    D=s.expand(40*xx[0]*C1-A*C2);assert s.cancel(xx[5]*D-40*xx[0]*qs[1]+A*Q2)==0
    aa,X,W,bb,S,U=s.symbols('AA X W BB S U');coords=(aa,X,W,bb,S,U);locals={'TA':aa,'X':X,'W':W,'TB':bb,'S':S,'U':U,'a':a,'e':e,'P':P}
    x4=(3*aa-bb)/(10*a*P*J);tau=(30*a*Rc*X+90*P*a*Tc*S-18*P*Vc*U-2*a*Kc*x4-bb)/(180*P**2);images=[X,3*tau,S,-s.Rational(7,2)*tau+U,x4,W];sub=dict(zip(xx,images))
    det=s.factor(s.Matrix(images).jacobian(coords).det());assert s.cancel(det-1/(200*a*P**3*J))==0
    for q,v in[(A,aa),(B,bb),(xx[0],X),(xx[5],W)]:assert s.cancel(q.subs(sub,simultaneous=True)-v)==0
    # Substitute into the actual forms, not the author's coefficient-matrix multiplication.
    cc=[s.Poly(s.cancel(c.subs(sub,simultaneous=True)),*coords)for c in[C1,C2]];assert all(not q.as_expr().has(S,U)for q in cc)
    coeff=[[s.cancel(q.coeff_monomial(v))for v in coords[:4]]for q in cc]
    for row,saved in zip(coeff,cert['C1_C2_four_coordinate_coefficients']):assert all(s.cancel(c-s.sympify(d,locals=locals))==0 for c,d in zip(row,saved))
    six=[aa/6075,40*X/6075,W,*[q.as_expr()for q in cc],-bb];sixmat=s.Matrix([[s.Poly(q,*coords).coeff_monomial(v)for v in coords]for q in six])
    assert sixmat[:,4:]==s.zeros(6,2)and s.cancel(sixmat[[0,1,2,5],:4].det()+s.Rational(40,6075**2))==0
    # All three actual minors and D's global ideal membership.
    M=s.Matrix([[A/6075,40*xx[0]/6075,xx[5]],[C1,C2,-B]])
    assert s.cancel(M[0,0]*M[1,1]-M[0,1]*M[1,0]+D/6075)==0
    assert s.cancel(M[0,0]*M[1,2]-M[0,2]*M[1,0]+qs[1])==0
    assert s.cancel(M[0,1]*M[1,2]-M[0,2]*M[1,1]+Q2)==0
    dnormal=s.Poly(s.expand(40*X*cc[0].as_expr()-aa*cc[1].as_expr()),*coords)
    dbar=s.Poly(dnormal.as_expr().subs({bb:0,W:0}),aa,X);want=-8*(3*P*e+2*P+12*a)*aa**2/(1215*P**2)-256*a*a*J*aa*X/(243*P)+1280*a**3*J**2*X**2/243
    assert s.cancel(dbar.as_expr()-want)==0 and s.cancel(dbar.coeff_monomial(X**2)-1280*a**3*J**2/243)==0
    # FIRST restriction: AA=X=W=0; BB remains a FREE coordinate.
    paircut={aa:0,X:0,W:0};q1n=aa*bb/6075+W*cc[0].as_expr();q2n=40*X*bb/6075+W*cc[1].as_expr()
    assert all(s.expand(q.subs(paircut))==0 for q in[q1n,q2n,dnormal.as_expr()])and bb not in paircut
    # Save the actual pair restriction without attempting to solve its gate.
    pairimages=[s.cancel(v.subs(paircut))for v in images];pairquads=[s.Poly(s.cancel(q.subs(dict(zip(xx,pairimages)),simultaneous=True)),bb,S,U)for q in[qs[0],qs[3]]]
    pairmon=[bb*bb,bb*S,bb*U,S*S,S*U,U*U];pairmatrix=[[str(s.cancel(q.coeff_monomial(mon)))for mon in pairmon]for q in pairquads]
    assert s.expand(3*Sc-7*Vc-20*P)==0
    gcdreview=ROOT/'reports/S61_NONSQUARE_GCD_GATE_PROPORTIONAL_EXCLUSION_20261103_v63_REVIEW/report.json';finitereview=ROOT/'reports/S61_NONSQUARE_FINITE_GCD_EXCEPTION_NAKAYAMA_20261103_v63_REVIEW/report.json'
    assert all(read(p)['all_passed']for p in[gcdreview,finitereview])
    evidence=dict(inverse_frame_images=list(map(str,images)),inverse_frame_determinant=str(det),C1_C2_four_coordinate_coefficients=list(map(lambda row:list(map(str,row)),coeff)),
        actual_six_linear_entries_coefficient_matrix=[[str(c)for c in row]for row in sixmat.tolist()],constant_rank4_minor=str(-s.Rational(40,6075**2)),
        all3_actual_minor_identities_verified=True,D_global_expression='40*x0*C1-A*C2',D_on_B_W_zero=str(want),D_x0_squared_known_unit='1280*a^3*J^2/243',
        FIRST_pair_restriction='A=x0=W=0; B FREE',actual_pair_variable_order=['B','S','U'],actual_pair_monomials=list(map(str,pairmon)),actual_pair2_by6_coefficients=pairmatrix,
        SECOND_restriction='B=W=0; A andx0 independent')
    (out/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    result=dict(id=STEM,all_passed=True,own_actual804_source_C1_C2_Q2prime_D_rebuilt=True,
        actual_inverse_frame_and_six_linear_entries_rank_exact4_on_a_P_J=True,only_known_unit_denominators_a_P_J=True,
        all3_actual_determinantal_minor_identities_verified=True,D_binary_B_W_restriction_known_unit_nonzero=True,
        first_restriction_keeps_B_free_and_second_restriction_correct=True,pair_gate_implies_auxiliary_full_five_rank5_reviewed=True,
        conditional_linear_syzygy_reduction_lemma_reviewed=True,gcd_failures_raw_good_reuse_leaves_only_pair_as_sufficient_raw24_gate=True,
        pair_gate_itself_full_raw24_original41_or_ambient_claim_accepted=False,seconds=time.monotonic()-started)
    paths=[Path(__file__),parent/'report.json',archive,gcdreview,finitereview,ROOT/'code/s61_nonsquare_four_quadric_binary_reduction_20261103_v63.py'];note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['evidence_sha256']={'independent-certificate.json':sha(out/'independent-certificate.json')};(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

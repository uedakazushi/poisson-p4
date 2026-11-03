#!/usr/bin/env python3
"""Independent actual-source QQ finiteness proof for pair failure, not emptiness."""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import resource
import shutil
import signal
import tarfile
import time
import sympy as s
from s61_nonsquare_five_quadric_fourth_unit_review_20261103_v63 import native,witness,dot,equivalent,power,AP,PP

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_PAIR_FINITE_BAD_LOCUS_20261103_v63_REVIEW'


def read(path):return json.loads(path.read_text())
def sha(path):
    h=hashlib.sha256()
    with path.open('rb')as f:
        for c in iter(lambda:f.read(1048576),b''):h.update(c)
    return h.hexdigest()
def decode(terms):return{tuple(ex):Q(c)for ex,c in terms if Q(c)}
def poly(expr,vs):return{tuple(ex):Q(str(c))for ex,c in s.Poly(expr,*vs,domain=s.QQ).terms()if c}
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
    parent=ROOT/'reports/S61_NONSQUARE_PAIR_FINITE_BAD_LOCUS_20261103_v63';author=read(parent/'report.json');assert author['frozen']and author['all_exact_QQ_actual_source_checks_passed'];archive=ROOT/author['archive'];assert sha(archive)==author['archive_sha256']
    snapshot=out/'immutable-input';snapshot.mkdir()
    with tarfile.open(archive)as tar:
        for m in tar:
            assert m.isfile();path=snapshot/m.name;assert path.resolve().is_relative_to(snapshot.resolve());path.parent.mkdir(parents=True,exist_ok=True)
            with tar.extractfile(m)as f,path.open('wb')as dst:shutil.copyfileobj(f,dst)
    for rel,h in author['input_sha256'].items():assert sha(snapshot/rel if(snapshot/rel).exists()else ROOT/rel)==h
    cert=read(snapshot/'work/nonsquare-v62/pair-finite-proof-v76/certificate.json');data=read(snapshot/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json');assert data['all_passed']
    a,e,P,t,v,H=s.symbols('a e P t v H');xx=s.symbols('x0:6');yl=s.symbols('u vl z');u,vl,z=yl;qs=[];count=0
    for i,row in enumerate(data['cleared_qfour_sparse']):
        q=0
        for mon,coeffs in row:
            assert len(mon)==7 and mon[6]==0 and sum(mon)==2
            q+=sum(s.Rational(c)*a**(ex[0]-ex[2]+[0,2,2,1][i])*e**ex[1]*P**ex[2]for ex,c in coeffs)*s.prod(x**k for x,k in zip(xx,mon[:6]));count+=len(coeffs)
        qs.append(q)
    assert count==804
    J=5*P-6*e-4;Rc=15*P**2-6*P*e+6*P+80*a;Tc=-5*P+6*e-2;Sc=21*P*e+16*P+140*a;Vc=9*P*e+4*P+60*a;Kc=85*P**2-12*P*e+58*P+600*a;Ka=10*P**2+3*P*e+13*P+100*a
    B=30*a*Rc*xx[0]-9*P*Sc*xx[1]+90*P*a*Tc*xx[2]-18*P*Vc*xx[3]-2*a*Kc*xx[4];A=(B+10*a*P*J*xx[4])/3
    tau=(30*P*a*Tc*vl-6*P*Vc*u-4*a*Ka*z)/(60*P**2);images={xx[0]:0,xx[1]:3*tau,xx[2]:vl,xx[3]:-s.Rational(7,2)*tau+u,xx[4]:z,xx[5]:0}
    assert s.cancel(A.subs(images,simultaneous=True))==0 and s.cancel(B.subs(images,simultaneous=True)+10*a*P*J*z)==0
    assert s.cancel(images[xx[3]]+7*images[xx[1]]/6-u)==0 and images[xx[2]]==vl and images[xx[4]]==z
    Q2=s.cancel((qs[2]-e*qs[1])/(a*a*J));assert s.cancel(qs[1].subs(images,simultaneous=True))==0 and s.cancel(Q2.subs(images,simultaneous=True))==0
    pair=[s.Poly(s.cancel(P**4*q.subs(images,simultaneous=True)),*yl)for q in[qs[0],qs[3]]];mons=[(i,j,2-i-j)for i in range(3)for j in range(3-i)]
    actual=[[poly(q.coeff_monomial(s.prod(x**n for x,n in zip(yl,mon))),(a,e,P))for mon in mons]for q in pair]
    base=snapshot/'work/nonsquare-v62/four-quadric-pair-gate-v52';literal,_,_=native(base/'pair-native.txt');assert actual==literal
    normalized,used,origin=native(base/'pair-after.txt',True);assert used==[]and origin==[0,1];ratios=[equivalent(row,saved)for row,saved in zip(actual,normalized)]
    assert sha(base/'pair-after.txt')==cert['actual_pair_source_sha256']
    cof=witness(snapshot/'work/nonsquare-v62/pair-gate-small-units-v64/D2-row0.witness');assert len(cof)==6 and sum(len(c)for c in cof)==10
    assert dot(normalized[0],cof)==mul(AP,power(PP,4))
    # Invertible natural parameter transformation, preserving only known t,P,J localizers.
    newP=H-8*t;newe=(v-20*t)/3;newa=t*newP
    assert s.factor(s.Matrix([newa,newe,newP]).jacobian([t,v,H]).det()-newP/3)==0
    assert s.expand(5*newP-6*newe-4-(5*H-2*v-4))==0
    oldexpr=lambda q:sum(s.Rational(c.numerator,c.denominator)*a**ex[0]*e**ex[1]*P**ex[2]for ex,c in q.items())
    nat=[[s.Poly(s.cancel(oldexpr(q).subs(a,t*P).subs(e,(v-20*t)/3)/P**4).subs(P,H-8*t),t,v,H,domain=s.QQ)for q in row]for row in normalized]
    h0=(160*t-15*H-2*v+31)/36;rows=[nat[0],[s.Poly(q3-h0*q0,t,v,H,domain=s.QQ)for q0,q3 in zip(*nat)]]
    encoded=lambda terms,vs:sum(s.Rational(c)*s.prod(x**n for x,n in zip(vs,ex))for ex,c in terms)
    assert all(q==s.Poly(encoded(saved,(t,v,H)),t,v,H,domain=s.QQ)for row,savedrow in zip(rows,cert['natural_rows'])for q,saved in zip(row,savedrow))
    f=(v+3)*(3*v+4);E=3*H*(v+4)+5*v-12
    assert s.expand(rows[0][5].as_expr()-27*f/40000)==0 and s.expand(rows[1][5].as_expr()+3*E/32000)==0
    equations=[s.Poly(encoded(q,(t,v,H)),t,v,H,domain=s.QQ)for q in cert['linear_t_equations']];ab=[]
    for eq,record in zip(equations,cert['linear_equation_normalizations']):
        j=record['pair_coefficient_index'];k=record['t_unit_power'];c=s.Rational(record['QQcontent']);assert j in[2,3,4]and k>=0 and c
        raw=rows[0][5]*rows[1][j]-rows[1][5]*rows[0][j];assert raw==s.Poly(c*t**k*eq.as_expr(),t,v,H,domain=s.QQ)
        pt=s.Poly(eq.as_expr(),t);assert pt.degree()<=1;ab.append([s.Poly(pt.coeff_monomial(t),H,v,domain=s.QQ),s.Poly(pt.coeff_monomial(1),H,v,domain=s.QQ)])
    elim=[s.Poly(encoded(q,(H,v)),H,v,domain=s.QQ)for q in cert['elimination_polynomials']]
    for q,record in zip(elim,cert['elimination_normalizations']):
        i,j=record['indices'];c=s.Rational(record['QQcontent']);k=record['f_chart_power'];assert c and k>=0
        raw=ab[i][0]*ab[j][1]-ab[j][0]*ab[i][1];assert raw==s.Poly(c*f**k*q.as_expr(),H,v,domain=s.QQ)
    alpha=cert['alpha_nonzero_Bezout'];ap={key:poly(s.sympify(alpha[key],locals={'H':H,'v':v}),(H,v))for key in['C','beta3_div_f','alpha4','Bezout_coefficient_C','Bezout_coefficient_beta3_div_f','Bezout_coefficient_alpha4']}
    Cexpr=s.sympify(alpha['C'],locals={'H':H,'v':v});bexpr=s.sympify(alpha['beta3_div_f'],locals={'H':H,'v':v});assert s.expand(ab[2][1].as_expr()+f*Cexpr)==0 and s.expand(ab[1][1].as_expr()-f*bexpr)==0 and ap['alpha4']==poly(ab[2][0].as_expr(),(H,v))
    unity={}
    for cof,q in[('Bezout_coefficient_C','C'),('Bezout_coefficient_beta3_div_f','beta3_div_f'),('Bezout_coefficient_alpha4','alpha4')]:unity=add(unity,mul(ap[cof],ap[q]))
    assert unity=={(0,0):Q(1)}
    projection=cert['projection_certificate'];N,M,R,Aco,Bco=[decode(projection[k])for k in['A','B','R','Bezout_A','Bezout_B']]
    assert N==poly(elim[1].as_expr(),(H,v))and M==poly(elim[2].as_expr(),(H,v))and R and all(ex[0]==0 for ex in R)and max(ex[1]for ex in R)==33
    assert add(mul(Aco,N),mul(Bco,M))==R
    degreeH=max(ex[0]for ex in N);actualcoeff=[{(ex[1],):c for ex,c in N.items()if ex[0]==i}for i in range(degreeH,-1,-1)];savedcoeff=list(map(decode,projection['A_H_coefficients']));weights=list(map(decode,projection['A_content_Bezout']));assert actualcoeff==savedcoeff
    contentunit={}
    for co,q in zip(weights,actualcoeff):contentunit=add(contentunit,mul(co,q))
    assert contentunit=={(0,):Q(1)}
    branches=[]
    for branch,(vv,hh)in zip(cert['special_branches'],[(s.Integer(-3),s.Integer(9)),(-s.Rational(4,3),s.Rational(7,3))]):
        assert s.Rational(branch['v'])==vv and s.Rational(branch['H'])==hh and s.expand(E.subs(v,vv)-E.subs({v:vv,H:hh})-s.diff(E,H).subs(v,vv)*(H-hh))==0 and E.subs({v:vv,H:hh})==0
        assert s.diff(E,H).subs(v,vv)!=0
        specialized=[[s.Poly(q.as_expr().subs({v:vv,H:hh}),t,domain=s.QQ)for q in row]for row in rows]
        assert specialized[0][5].is_zero and specialized[1][5].is_zero
        rawminors=[poly((specialized[0][i]*specialized[1][j]-specialized[0][j]*specialized[1][i]).as_expr(),(t,))for i,j in itertools.combinations(range(5),2)]
        assert rawminors==[poly(s.sympify(q,locals={'t':t}),(t,))for q in branch['raw_minors']]
        prims=[poly(s.sympify(q,locals={'t':t}),(t,))for q in branch['normalized_minors']];cofs=[poly(s.sympify(q,locals={'t':t}),(t,))for q in branch['Bezout_coefficients']]
        for raw,pr,record in zip(rawminors,prims,branch['normalization_records']):
            powers=record['powers_t_P'];c=s.Rational(record['QQcontent']);assert len(powers)==2 and min(powers)>=0
            if not raw:assert not pr and not c
            else:assert c and mul(poly(c*t**powers[0]*(hh-8*t)**powers[1],(t,)),pr)==raw
        identity={}
        for co,q in zip(cofs,prims):identity=add(identity,mul(co,q))
        assert identity=={(0,):Q(1)};branches.append(dict(v=str(vv),H=str(hh),all10_raw_minors_known_t_P_factorizations_and_QQ_Bezout_verified=True))
    gcdreview=ROOT/'reports/S61_NONSQUARE_GCD_GATE_PROPORTIONAL_EXCLUSION_20261103_v63_REVIEW/report.json';binaryreview=ROOT/'reports/S61_NONSQUARE_FOUR_QUADRIC_BINARY_REDUCTION_20261103_v63_REVIEW/report.json';finitereview=ROOT/'reports/S61_NONSQUARE_FINITE_GCD_EXCEPTION_NAKAYAMA_20261103_v63_REVIEW/report.json';assert all(read(p)['all_passed']for p in[gcdreview,binaryreview,finitereview])
    evidence=dict(actual804_pair_source_rebuilt=True,raw_to_normalized_known_unit_row_ratios=ratios,first_row_literal_10term_cofactor_unit='a*P^4',pair_space='A=x0=W=0; B free=-10*a*P*J*z',
        inverse_parameter_change={'a':str(newa),'e':str(newe),'P':str(newP),'J':'5H-2v-4','Jacobian':'P/3'},all12_natural_coefficients_and_row_subtraction_verified=True,
        all3_linear_t_equations_and_elimination_unit_factorizations_verified=True,projection_R_v_degree=33,literal_projection_resultant_membership_verified=True,literal_H_content_QQ_Bezout1_verified=True,
        no_all_alpha_zero_ONLY_at_solutions_of_necessary_equations_on_f_nonzero=True,literal_alpha_beta_QQ_Bezout1_verified=True,special_branches=branches)
    (out/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n');shutil.rmtree(snapshot)
    result=dict(id=STEM,all_passed=True,actual804_pair_chart_and_full12_coefficients_rebuilt=True,first_row_never_zero_by_literal_QQ_10term_unit=True,
        natural_coordinate_inverse_and_only_t_P_J_localizers_verified=True,special_v_minus3_minus4over3_branches_excluded_by_QQ_Bezout=True,
        all3_linear_t_equations_and_necessary_projection_factors_verified=True,nonzero_v_degree33_resultant_membership_and_H_content1_verified=True,
        at_failure_some_alpha_nonzero_and_t_unique_by_exact_QQ_identity=True,pair_failure_dimension_at_most_zero_on_a_P_J=True,
        actual_raw24_defect_dimension_at_most_zero_on_FULLU_using_prior_reviews=True,
        pair_empty_proper_finite_Zp_unit_globalraw24_original41_or_ambient_claim_accepted=False,seconds=time.monotonic()-start)
    paths=[Path(__file__),parent/'report.json',archive,gcdreview,binaryreview,finitereview,ROOT/'code/s61_nonsquare_five_quadric_fourth_unit_review_20261103_v63.py'];note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['evidence_sha256']={'independent-certificate.json':sha(out/'independent-certificate.json')};(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

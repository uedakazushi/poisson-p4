#!/usr/bin/env python3
"""Independent finite-Z_(32003) Sigma raw-module proof review; no GB construction."""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import tarfile
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1];PCHAR=32003
STEM='S61_NONSQUARE_FINITE_GCD_EXCEPTION_NAKAYAMA_20261103_v63_REVIEW'


def read(path):return json.loads(path.read_text())
def sha(path):
    h=hashlib.sha256()
    with path.open('rb')as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return h.hexdigest()
def qqpoly(expr,variables):return{tuple(ex):Q(str(c))for ex,c in s.Poly(expr,*variables,domain=s.QQ).terms()}
def decode(terms):return{tuple(ex):Q(c)for ex,c in terms if Q(c)}
def mod(q):
    out={}
    for ex,c in q.items():
        c=Q(c);assert c.denominator%PCHAR;v=c.numerator*pow(c.denominator,-1,PCHAR)%PCHAR
        if v:out[ex]=v
    return out
def add(a,b,c=1):
    out=a.copy()
    for ex,v in b.items():
        w=(out.get(ex,0)+c*v)%PCHAR
        if w:out[ex]=w
        elif ex in out:del out[ex]
    return out
def mul(a,b):
    out={}
    for x,c in a.items():
        for y,d in b.items():
            ex=tuple(i+j for i,j in zip(x,y));out[ex]=(out.get(ex,0)+c*d)%PCHAR
    return{ex:v for ex,v in out.items()if v}
def shift(q,ex,c=1):return{tuple(i+j for i,j in zip(x,ex)):v*c%PCHAR for x,v in q.items()if v*c%PCHAR}
def qmul(a,b):
    out={}
    for x,c in a.items():
        for y,d in b.items():
            ex=tuple(i+j for i,j in zip(x,y));out[ex]=out.get(ex,Q(0))+c*d
    return{ex:v for ex,v in out.items()if v}
def qadd(a,b):
    out=a.copy()
    for ex,v in b.items():out[ex]=out.get(ex,Q(0))+v
    return{ex:v for ex,v in out.items()if v}
def native(path,nvars):
    tokens=iter(path.read_text().split());n,m,flag=[int(next(tokens))for _ in range(3)];assert flag==0;rows=[]
    for _ in range(n):
        row=[]
        for _ in range(m):
            q={}
            for _ in range(int(next(tokens))):
                c=Q(next(tokens));ex=tuple(int(next(tokens))for _ in range(nvars));assert ex not in q;q[ex]=c
            row.append(q)
        rows.append(row)
    assert next(tokens,None)is None;return rows


def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(120)
    parent=ROOT/'reports/S61_NONSQUARE_FINITE_GCD_EXCEPTION_NAKAYAMA_20261103_v63';author=read(parent/'report.json');assert author['author_exact_proof_complete']and author['frozen']
    archive=ROOT/author['archive'];assert sha(archive)==author['archive_sha256'];snapshot=out/'immutable-input';snapshot.mkdir()
    manifest={v['path']:v for v in author['members']}
    with tarfile.open(archive)as tar:
        assert set(tar.getnames())==set(manifest)
        for member in tar.getmembers():
            assert member.isfile();path=snapshot/member.name;assert path.resolve().is_relative_to(snapshot.resolve());path.parent.mkdir(parents=True,exist_ok=True)
            with tar.extractfile(member)as src,path.open('wb')as dst:shutil.copyfileobj(src,dst)
            assert path.stat().st_size==manifest[member.name]['bytes']and sha(path)==manifest[member.name]['sha256']
    for rel,h in author['input_sha256'].items():assert sha(snapshot/rel if rel in manifest else ROOT/rel)==h
    data=snapshot/'work/nonsquare-v62/rank2-finite-model-v59';old=read(snapshot/'work/nonsquare-v62/gcd-gate-exception-v54/certificate.json')
    src=read(snapshot/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json');assert src['all_passed']
    assert sha(snapshot/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json')=='80f7c46b1bd4ffef4a78f492ad4f7a2bb40b2ffa27cddcd23d5b5d660db4a7e0'
    # Reconstruct every literal row directly from804 actual quadratic parameter terms, not a Schur source.
    quads=[];qcount=0;shifts=[]
    for row in src['cleared_qfour_sparse']:
        sh=max(0,max(ex[2]-ex[0]for mon,terms in row for ex,c in terms));shifts.append(sh);quad={}
        for mon,terms in row:
            assert len(mon)==7 and mon[-1]==0 and sum(mon)==2;co={}
            for ex,c in terms:
                c=Q(c);assert c.denominator%PCHAR;new=(ex[0]-ex[2]+sh,ex[1],ex[2]);assert min(new)>=0;assert new not in co;co[new]=c;qcount+=1
            quad[tuple(mon[:6])]=co
        quads.append(quad)
    assert shifts==[0,2,2,1]and qcount==804
    monomials=list(itertools.combinations_with_replacement(range(6),3));monomials=sorted(tuple(word.count(i)for i in range(6))for word in monomials);assert len(monomials)==56
    products=[]
    for quad in quads:
        for k in range(6):products.append({tuple(e+int(i==k)for i,e in enumerate(mon)):co for mon,co in quad.items()})
    rebuilt=[[row.get(mon,{})for mon in monomials]for row in products]
    literalpath=snapshot/'work/nonsquare-v62/raw-fourquad-GF-v36/M-QQ-native.txt';literal=native(literalpath,3);assert literal==rebuilt and sum(len(q)for row in literal for q in row)==4824
    assert sha(literalpath)=='f9590655c56c4c7790d43bc780950ff9e7bf4a3339a4a0879e2901369a3b78aa'
    # Bind N and the restrictedQ3 minor to the actual four quadrics, including DB elimination.
    a,e,P=s.symbols('a e P');xx=s.symbols('x0:6');zz=s.symbols('z0:4');variables=(P,e)
    expression=lambda q:sum(s.Rational(c.numerator,c.denominator)*a**ex[0]*e**ex[1]*P**ex[2]for ex,c in q.items())
    QQ=[sum(expression(co)*s.prod(v**i for v,i in zip(xx,mon))for mon,co in q.items())for q in quads]
    Rc=15*P**2-6*P*e+6*P+80*a;Tc=-5*P+6*e-2;Sc=21*P*e+16*P+140*a;Vc=9*P*e+4*P+60*a;Kc=85*P**2-12*P*e+58*P+600*a;J=5*P-6*e-4;Lp=20*P+3*e+20
    B=30*a*Rc*xx[0]-9*P*Sc*xx[1]+90*P*a*Tc*xx[2]-18*P*Vc*xx[3]-2*a*Kc*xx[4]
    Afactor=(B+10*a*P*J*xx[4])/3
    assert s.expand(QQ[1].subs(xx[5],0)-Afactor*B/6075)==0
    assert s.cancel(((QQ[2]-e*QQ[1])/(a*a*J)).subs(xx[5],0)-40*xx[0]*B/6075)==0
    assert s.expand(3*Sc-7*Vc-20*P)==0
    r,u,v,t,R=s.symbols('r u v t R');sub={xx[0]:r,xx[1]:v,xx[2]:R/a-u-6*(e+1)*v/(5*a)-(9*e+4)*t/(10*a),xx[3]:(t-v)/2,xx[4]:5*u+5*r-4*R/a+24*(e+1)*v/(5*a)+2*(9*e+4)*t/(5*a),xx[5]:0}
    assert s.cancel(4*QQ[0].subs(sub,simultaneous=True)-u*R+5*v*t)==0
    bpoly=s.Poly(s.cancel(B.subs(sub,simultaneous=True)),r,u,v,t,R);bc=[bpoly.coeff_monomial(x)for x in[r,u,v,t,R]];DB=20*P**2+3*P*e+20*P+180*a
    assert s.expand(bc[0]+20*a*DB)==0;aelim=-P*Lp/180;norm=bc[1]*bc[4]-bc[2]*bc[3]/5;N=s.Poly(s.cancel(-s.Rational(10125,4)*norm.subs(a,aelim)/P**2),P,e,domain=s.QQ)
    cut=(30*a*Rc*zz[0]+90*P*a*Tc*zz[1]-18*P*Vc*zz[3]-2*a*Kc*zz[2])/(180*P**2)
    q3=s.Poly(s.cancel(P**4*QQ[3].subs({xx[0]:zz[0],xx[1]:3*cut,xx[2]:zz[1],xx[3]:-s.Rational(7,2)*cut+zz[3],xx[4]:zz[2],xx[5]:0},simultaneous=True).subs(a,aelim)),*zz)
    H=s.Matrix(4,4,lambda i,j:s.cancel(q3.coeff_monomial(zz[i]*zz[j])/(1 if i==j else 2)))
    load=lambda terms:sum(s.Rational(c)*P**ex[0]*e**ex[1]for ex,c in terms)
    assert H==s.Matrix([[load(q)for q in row]for row in old['global_H3_entries']])and N==s.Poly(load(old['N']),P,e,domain=s.QQ)
    model=read(data/'primitive-input.json');assert N==s.Poly(s.sympify(model['N'],locals={'P':P,'e':e}),P,e,domain=s.QQ);M=s.Poly(s.sympify(model['primitive_minor'],locals={'P':P,'e':e}),P,e,domain=s.QQ)
    Hsmall=H.applyfunc(lambda q:s.cancel(q/(P**6*Lp**2*J)));assert all(s.Poly(q,P,e,domain=s.QQ)is not None for q in Hsmall)
    content=Q(model['minor_QQ_content']);assert content.numerator%PCHAR and content.denominator%PCHAR
    assert s.Poly(Hsmall[:3,:3].det(method='domain-ge'),P,e,domain=s.QQ)==s.Poly(P*Lp*s.Rational(content.numerator,content.denominator)*M.as_expr(),P,e,domain=s.QQ)
    # Check the literal QQ Bezout identity using independent sparse Fraction arithmetic.
    membership=read(data/'resultant-membership-certificate.json');nq=qqpoly(N.as_expr(),variables);mq=qqpoly(M.as_expr(),variables);rq,Aq,Bq=[decode(membership[k])for k in['R','A','B']]
    assert qadd(qmul(Aq,nq),qmul(Bq,mq))==rq
    for q in[nq,mq,rq,Aq,Bq]:assert all(c.denominator%PCHAR for c in q.values())
    assert max(ex[1]for ex in nq)==4 and nq[0,4]==72442188 and nq[0,4].numerator%PCHAR
    assert all(ex[1]==0 for ex in rq)and max(ex[0]for ex in rq)==28 and rq[28,0].numerator%PCHAR
    U=aelim*P*(500*aelim+90*e*e+3*e-24)*((50*aelim+18*e*e-9*e-2)*P+(60*e-40)*aelim)*(100*aelim+18*e*e-9*e-2)*J;uq=qqpoly(U,variables)
    models=[mod(q)for q in[nq,mq,uq]];assert models==[mod(q)for q in native(data/'finite-model-polys-native2.txt',2)[0]]
    G=[mod(q)for q in native(data/'finite-G-native2.txt',2)[0]];lifts=[[mod(q)for q in row]for row in native(data/'finite-G-lift-native2.txt',2)];wit=[[mod(q)for q in row]for row in native(data/'finite-U-bezout-native2.txt',2)]
    assert len(G)==5 and len(lifts)==2 and len(wit)==3;order=lambda ex:(sum(ex),ex[0]);heads=[max(g,key=order)for g in G];assert all(g[ex]==1 for g,ex in zip(G,heads))
    def nf(q):
        work=q.copy();rem={}
        while work:
            mon=max(work,key=order);c=work[mon]
            for g,lm in zip(G,heads):
                if mon[0]>=lm[0]and mon[1]>=lm[1]:work=add(work,shift(g,(mon[0]-lm[0],mon[1]-lm[1]),c),-1);break
            else:rem[mon]=c;del work[mon]
        return rem
    for j,g in enumerate(G):assert add(mul(models[0],lifts[0][j]),mul(models[1],lifts[1][j]))==g
    assert all(not nf(q)for q in models[:2])
    for i,j in itertools.combinations(range(5),2):
        common=tuple(max(a,b)for a,b in zip(heads[i],heads[j]));S=add(shift(G[i],tuple(a-b for a,b in zip(common,heads[i]))),shift(G[j],tuple(a-b for a,b in zip(common,heads[j]))),-1);assert not nf(S)
    for i,g in enumerate(G):assert all(not all(x>=y for x,y in zip(mon,lm))for mon in g for j,lm in enumerate(heads)if j!=i)
    unit={}
    for model,cof in zip(models,wit):unit=add(unit,mul(model,cof[0]))
    assert unit=={(0,0):1};uinv=nf(wit[2][0]);assert nf(add(mul(models[2],uinv),{(0,0):1},-1))=={}
    limP=min(x for x,y in heads if y==0);lime=min(y for x,y in heads if x==0);basis=sorted([(i,j)for i in range(limP)for j in range(lime)if not any(i>=a and j>=b for a,b in heads)],key=lambda ex:(sum(ex),*ex));assert len(basis)==28
    assert basis==list(map(tuple,read(data/'algebra-verification.json')['basis']));index={ex:i for i,ex in enumerate(basis)};operators=[]
    for inc in[(1,0),(0,1)]:
        op=[[0]*28 for _ in range(28)]
        for j,ex in enumerate(basis):
            red=nf({(ex[0]+inc[0],ex[1]+inc[1]):1});assert set(red)<=set(index)
            for mon,c in red.items():op[index[mon]][j]=c
        operators.append(op)
    encoded=list(map(int,(data/'multiplication-matrices.txt').read_text().split()));assert encoded[:2]==[28,PCHAR]and encoded[2:]==[c for op in operators for row in op for c in row]
    opsfile=out/'own-operators.txt';opsfile.write_text('28 32003\n'+'\n'.join(' '.join(map(str,row))for op in operators for row in op)+'\n')
    rawfile=out/'own-raw-modp.txt'
    with rawfile.open('w')as f:
        f.write('24 56 3\n')
        for row in literal:
            for q in row:
                f.write(str(len(q))+'\n')
                for ex,c in sorted(q.items()):f.write(str(c.numerator*pow(c.denominator,-1,PCHAR)%PCHAR)+' '+' '.join(map(str,ex))+'\n')
    polyfile=out/'own-model-uinverse.txt'
    with polyfile.open('w')as f:
        f.write('4 2\n')
        for q in[*models,uinv]:
            f.write(str(len(q))+'\n')
            for ex,c in sorted(q.items()):f.write(str(c)+' '+' '.join(map(str,ex))+'\n')
    selected=read(data/'raw-block-rank.json')['selected_columns'];assert len(selected)==672 and len(set(selected))==672
    choose=out/'selected-columns.txt';choose.write_text('672\n'+' '.join(map(str,selected))+'\n')
    cpp=ROOT/'code/s61_nonsquare_finite_sigma_module_review_20261103_v63.cpp';binary=out/'own-native-bin';subprocess.run(['g++','-O3','-std=c++17',str(cpp),'-o',str(binary)],check=True,timeout=30)
    nativeout=out/'own-native.json';subprocess.run([str(binary),str(opsfile),str(rawfile),str(polyfile),str(choose),str(data/'selected-square-inverse.txt'),str(nativeout)],check=True,timeout=30);nr=read(nativeout);assert nr['all_passed']
    # All immutable inputs remain in the author's hash-bound archive; do not publish ELF or duplicate snapshots.
    binary.unlink();shutil.rmtree(snapshot)
    result=dict(id=STEM,all_passed=True,own_actual804_quadratic_terms_and_all4824_raw_product_coefficients_match=True,
        own_actual_DB_dual_norm_N_and_restricted_Q3_minor_M_verified=True,exact_QQ_R_equals_A_N_plus_B_M_verified=True,
        p_integral_model_and_e4_P28_leading_coefficients_units=True,finite_V_algebra_without_flatness_reviewed=True,
        own_GF_GB_both_inclusions_all10_S_pairs_reduced_basis28_verified=True,FULLU_unit_on_ENTIRE_special_fiber_and_finite_V_model_reviewed=True,
        all1568_multiplication_operator_entries_reconstructed=True,independent_plain_integer_native_full672x1568_module_both_inverse_products_verified=True,
        finite_V_cokernel_Nakayama_QQ_surjectivity_on_entire_finite_Sigma_reviewed=True,
        all_U_global_GF_promotion_CRT_or_original41_claim_accepted=False,seconds=time.monotonic()-started)
    note=ROOT/f'research/notes/{STEM}.md';paths=[Path(__file__),cpp,parent/'report.json',archive]
    if note.exists():paths.append(note)
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['evidence_sha256']={p.name:sha(p)for p in[opsfile,rawfile,polyfile,choose,nativeout]};(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);main(parser.parse_args().output)

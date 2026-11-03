#!/usr/bin/env python3
"""Independent QQ/integer polynomial audit of four AUXILIARY five-quadric pivots."""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
import math
from pathlib import Path
import resource
import shutil
import signal
import tarfile
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_FIVE_QUADRIC_FOURTH_UNIT_20261103_v63_REVIEW'
Z=(0,0,0);JA={(0,0,1):5,(0,1,0):-6,Z:-4};AP={(1,0,0):1};PP={(0,0,1):1}


def read(path):return json.loads(path.read_text())
def sha(path):
    h=hashlib.sha256()
    with path.open('rb')as f:
        for c in iter(lambda:f.read(1048576),b''):h.update(c)
    return h.hexdigest()
def add(a,b,c=1):
    out=a.copy()
    for ex,v in b.items():out[ex]=out.get(ex,0)+c*v
    return{ex:c for ex,c in out.items()if c}
def scale(q,c):return{ex:v*c for ex,v in q.items()if v*c}
def mul(a,b):
    out={}
    for ex,c in a.items():
        for ey,d in b.items():
            m=tuple(x+y for x,y in zip(ex,ey));out[m]=out.get(m,0)+c*d
    return{ex:c for ex,c in out.items()if c}
def power(q,n):
    out={Z:1}
    for _ in range(n):out=mul(out,q)
    return out
def parse_poly(tokens):
    out={}
    for _ in range(int(next(tokens))):
        c=Q(next(tokens));ex=tuple(int(next(tokens))for _ in range(3));assert ex not in out;out[ex]=c
    return out
def native(path,metadata=False):
    tokens=iter(path.read_text().split());n,m,k=[int(next(tokens))for _ in range(3)];used=[];orig=[]
    if metadata:used=[int(next(tokens))for _ in range(k)];orig=[int(next(tokens))for _ in range(n)]
    else:assert k==0
    rows=[[parse_poly(tokens)for _ in range(m)]for _ in range(n)];assert next(tokens,None)is None
    return rows,used,orig
def witness(path):
    tokens=iter(path.read_text().split());n=int(next(tokens));out=[parse_poly(tokens)for _ in range(n)];assert next(tokens,None)is None;return out
def clear(row):
    D=1
    for q in row:
        for c in q.values():D=math.lcm(D,Q(c).denominator)
    return[{ex:Q(c).numerator*(D//Q(c).denominator)for ex,c in q.items()}for q in row],D
def divideJ(q):
    if not q:return{}
    q=q.copy();out={}
    while q:
        ex=max(q);c=q[ex]
        if ex[1]==0 or c%(-6):return None
        shift=(ex[0],ex[1]-1,ex[2]);v=c//(-6);out[shift]=out.get(shift,0)+v
        for m,t in JA.items():
            key=tuple(a+b for a,b in zip(shift,m));q[key]=q.get(key,0)-v*t
            if not q[key]:del q[key]
    return{ex:c for ex,c in out.items()if c}
def strip(row):
    row=[q.copy()for q in row];powers=[]
    for pos in[0,2]:
        n=min(ex[pos]for q in row for ex in q);powers.append(n)
        if n:row=[{tuple(v-n if i==pos else v for i,v in enumerate(ex)):c for ex,c in q.items()}for q in row]
    n=0
    while True:
        divided=[divideJ(q)for q in row]
        if any(q is None for q in divided):break
        row=divided;n+=1
    powers.append(n);content=0
    for q in row:
        for c in q.values():content=math.gcd(content,abs(c))
    first=next(q for q in row if q);content*=1 if first[max(first)]>0 else-1
    assert content
    return[{ex:c//content for ex,c in q.items()}for q in row],content,powers
def primitive(row):
    z,D=clear(row);z,c,p=strip(z);return[{ex:Q(v)for ex,v in q.items()}for q in z],dict(content=str(Q(c,D)),powers_a_P_J=p)
def equivalent(row,expected):
    ra,da=clear(row);rb,db=clear(expected);a,ca,pa=strip(ra);b,cb,pb=strip(rb);assert a==b
    return dict(QQ_ratio=str(Q(ca*db,da*cb)),powers_a_P_J=[x-y for x,y in zip(pa,pb)])
def dot(row,cofs):
    out={}
    for q,c in zip(row,cofs):out=add(out,mul(q,c))
    return out


def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(120)
    stems=['S61_NONSQUARE_FOUR_QUADRIC_STRUCTURE_20261103_v63','S61_NONSQUARE_FIVE_QUADRIC_RANK_GATE_FRONTIER_20261103_v63','S61_NONSQUARE_FIVE_QUADRIC_FOURTH_UNIT_20261103_v63'];reports=[];snapshot=out/'immutable-input';snapshot.mkdir()
    archives=[]
    for stem in stems:
        r=read(ROOT/'reports'/stem/'report.json');assert r['frozen'];archive=ROOT/r['archive'];assert sha(archive)==r['archive_sha256'];archives.append(archive);reports.append(r)
        with tarfile.open(archive)as tar:
            for m in tar:
                assert m.isfile();path=snapshot/m.name;assert path.resolve().is_relative_to(snapshot.resolve());path.parent.mkdir(parents=True,exist_ok=True)
                with tar.extractfile(m)as src:
                    if path.exists():
                        h=hashlib.sha256()
                        for chunk in iter(lambda:src.read(1048576),b''):h.update(chunk)
                        assert sha(path)==h.hexdigest()
                    else:
                        with path.open('wb')as dst:shutil.copyfileobj(src,dst)
    for r in reports:
        for rel,h in r['input_sha256'].items():assert sha(snapshot/rel if(snapshot/rel).exists()else ROOT/rel)==h
    base=snapshot/'work/nonsquare-v62/raw-fourquad-structure-v49';third=snapshot/'work/nonsquare-v62/five-quadric-minors-v50';fourth=snapshot/'work/nonsquare-v62/five-quadric-degree23-v60';last=snapshot/'work/nonsquare-v62/five-rank-last-row-v61'
    assert(base/'aPJ-units.txt').read_text().splitlines()==['a','d','5*d-6*e-4']
    source=read(snapshot/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json');a,e,P=s.symbols('a e P');xx=s.symbols('x0:6');qs=[]
    for i,row in enumerate(source['cleared_qfour_sparse']):qs.append(sum(sum(s.Rational(c)*a**(ex[0]-ex[2]+[0,2,2,1][i])*e**ex[1]*P**ex[2]for ex,c in terms)*s.prod(x**n for x,n in zip(xx,mon[:6]))for mon,terms in row))
    J=5*P-6*e-4;Rc=15*P**2-6*P*e+6*P+80*a;Tc=-5*P+6*e-2;Sc=21*P*e+16*P+140*a;Vc=9*P*e+4*P+60*a;Kc=85*P**2-12*P*e+58*P+600*a
    B=30*a*Rc*xx[0]-9*P*Sc*xx[1]+90*P*a*Tc*xx[2]-18*P*Vc*xx[3]-2*a*Kc*xx[4];A=(B+10*a*P*J*xx[4])/3;Q2=s.cancel((qs[2]-e*qs[1])/(a*a*J));D=s.cancel((40*xx[0]*qs[1]-A*Q2)/xx[5]);assert not s.denom(D).has(a,e,P,*xx)
    polys=[s.Poly(q,*xx)for q in[qs[0],qs[1],Q2,qs[3],D]];mons=sorted(set().union(*(q.monoms()for q in polys)));assert len(mons)==21
    rebuilt=[[{ex:Q(str(c))for ex,c in s.Poly(q.coeff_monomial(s.prod(x**n for x,n in zip(xx,mon))),a,e,P,domain=s.QQ).terms()}if q.coeff_monomial(s.prod(x**n for x,n in zip(xx,mon)))else{}for mon in mons]for q in polys]
    original,_,_=native(base/'five-native.txt');assert original==rebuilt
    # First pivot literal(a/4) at row0,col2. Independent exact remaining-row normalization.
    assert original[0][2]=={(1,0,0):Q(1,4)};columns=[j for j in range(21)if j!=2]
    stage=[]
    for i in range(1,5):
        row=[add(mul(original[0][2],original[i][j]),mul(original[i][2],original[0][j]),-1)for j in columns];row,_=primitive(row);stage.append(row)
    # Exact QQ constant-span second pivot, on the original Q2prime row.
    target=mul(PP,JA);monomials=sorted(set(target).union(*(set(q)for q in stage[1])));M=s.Matrix([[s.Rational(q.get(mon,0))for q in stage[1]]for mon in monomials]);rhs=s.Matrix([target.get(mon,0)for mon in monomials])
    pc=list(DomainMatrix.from_Matrix(M).convert_to(s.QQ).rref()[1]);pr=list(DomainMatrix.from_Matrix(M[:,pc].T).convert_to(s.QQ).rref()[1]);solution=M[pr,pc].inv()*rhs[pr,:];cofs=[Q(0)]*20
    for j,c in zip(pc,solution):cofs[j]=Q(str(c))
    assert M*s.Matrix(list(map(s.Rational,cofs)))==rhs;replace=columns.index(3);assert cofs[replace]
    newcol=[dot(row,[{Z:c}if c else{}for c in cofs])for row in stage];assert newcol[1]==target
    remaining=[]
    for i in[0,2,3]:remaining.append([add(mul(target,stage[i][j]),mul(newcol[i],stage[1][j]),-1)for j in range(20)if j!=replace])
    saved3,u3,o3=native(base/'five-after.txt',True);assert u3==[2,3]and o3==[1,3,4];first2ratios=[equivalent(r,q)for r,q in zip(remaining,saved3)]
    # Third pivot: explicit invertible constant-column replacement, no search.
    cw3=witness(third/'primitive-row1.witness');assert len(cw3)==19;g3=mul(power(AP,3),mul(PP,power(JA,2)));assert dot(saved3[1],cw3)==g3
    assert cw3[5]=={Z:Q(55000,3)}and cw3[7]=={Z:Q(-11000,3)}and sum(bool(q)for q in cw3)==2
    v3=[dot(row,cw3)for row in saved3];after3=[[add(mul(g3,saved3[i][j]),mul(v3[i],saved3[1][j]),-1)for j in range(19)if j!=5]for i in[0,2]]
    saved2,u2,o2=native(third/'after-primitive.txt',True);assert u2==[5]and o2==[0,2];thirdratios=[equivalent(r,q)for r,q in zip(after3,saved2)]
    # Fourth: integer-cleared175-term dot and fraction-free Schur, independently of FLINT.
    rowfile,_,_=native(fourth/'source-row0.txt');assert rowfile==[saved2[0]];cw4=witness(fourth/'D2-row0.witness');assert len(cw4)==18 and all(cw4)and sum(len(q)for q in cw4)==175
    assert max(sum(ex)for q in cw4 for ex in q)==4
    ai0,d0=clear(saved2[0]);ai1,d1=clear(saved2[1]);ci,dc=clear(cw4);g4=mul(power(AP,3),mul(power(PP,3),power(JA,2)));v0=dot(ai0,ci);v1=dot(ai1,ci);assert v0==scale(g4,d0*dc)
    rawSchur=[add(scale(mul(g4,q1),d0*dc),mul(v1,q0),-1)for q0,q1 in zip(ai0,ai1)];denSchur=d0*d1*dc
    # Strip only known units and one QQ scalar, retaining every unknown polynomial factor.
    rawPrimitive,rawContent,rawPowers=strip(rawSchur);savedLast,usedLast,originLast=native(fourth/'last-row.txt',True);assert usedLast==[18]and originLast==[1]and len(savedLast)==1
    li,dl=clear(savedLast[0]);lastPrimitive,lastContent,lastPowers=strip(li);assert rawPrimitive==lastPrimitive
    fourthratio=dict(QQ_ratio=str(Q(rawContent*dl,denSchur*lastContent)),powers_a_P_J=[x-y for x,y in zip(rawPowers,lastPowers)])
    assert sum(len(q)for q in savedLast[0])==8475 and max(sum(ex)for q in savedLast[0]for ex in q)==19
    prim,_,_=native(last/'primitive-row.txt');assert len(prim)==1;scalings=read(last/'column-scalings.json');assert len(scalings)==18
    for j,datum in enumerate(scalings):
        assert datum['column']==j;coeff=Q(datum['QQ_content']);pa,pp,pj=datum['powers_a_P_J'];assert coeff and min(pa,pp,pj)>=0 and all(c.denominator==1 for c in prim[0][j].values())
        factor=mul(power(AP,pa),mul(power(PP,pp),power(JA,pj)));leftInt,ld=clear([savedLast[0][j]]);right=mul(factor,{ex:int(c)for ex,c in prim[0][j].items()})
        assert scale(right,ld*coeff.numerator)==scale(leftInt[0],coeff.denominator)
    assert sum(len(q)for q in prim[0])==7865 and max(sum(ex)for q in prim[0]for ex in q)==18
    evidence=dict(second_pivot_constant_cofactors=list(map(str,cofs)),second_pivot_known_unit='P*J',first2_output_row_known_unit_ratios=first2ratios,
        third_original_row_identity='(55000/3)M[1,5]-(11000/3)M[1,7]=a^3*P*J^2',third_output_row_known_unit_ratios=thirdratios,
        fourth_literal_integer_dot_identity='a^3*P^3*J^2',fourth_output_row_known_unit_ratio=fourthratio,all18_primitive_column_scalings_verified=True,
        saved_expanded_output_sha256=sha(fourth/'last-row.txt'),saved_primitive_output_sha256=sha(last/'primitive-row.txt'))
    (out/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n');shutil.rmtree(snapshot)
    result=dict(id=STEM,all_passed=True,actual5x21_from804_qfour_coefficients_rebuilt=True,first_two_known_unit_pivots_independently_verified=True,
        third_exact_constant_witness_and_2x18_output_verified=True,fourth_all175_QQ_cofactor_terms_integer_dot_verified=True,
        safe_append18_unitpivot_and_all1x18_Schur_row_identities_verified=True,all18_primitive_known_unit_column_identities_and_output_hashes_verified=True,
        auxiliary_rank_at_least4_on_a_P_J_open=True,only_inverted_polynomials=['a','P','5P-6e-4'],main_constructive34_unchanged=True,
        full_aux_rank5_pair_all_QQ_product41_or_ambient_claim_accepted=False,seconds=time.monotonic()-started)
    paths=[Path(__file__),*[ROOT/'reports'/stem/'report.json'for stem in stems],*archives];note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['evidence_sha256']={'independent-certificate.json':sha(out/'independent-certificate.json')};(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

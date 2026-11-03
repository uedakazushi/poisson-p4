#!/usr/bin/env python3
"""Independent companion-module reconstruction and selected670-minor audit."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_AFFINE_NORMAL_QUINTIC_CONSTANT_BARRIER_20261103_v63_REVIEW'
P=32003


def sha(path):
    h=hashlib.sha256()
    with path.open('rb')as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def read(path):
    with (gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)


def decode(poly):
    answer={}
    for degree,terms in poly:
        assert degree not in answer
        coefficients={}
        for mon,text in terms:
            mon=tuple(mon);assert mon==tuple(sorted(mon))and mon not in coefficients
            value=Fraction(text);assert value;coefficients[mon]=value
        if coefficients:answer[degree]=coefficients
    return answer


def add_scaled(target,source,coefficient):
    for mon,value in source.items():
        result=target.get(mon,Fraction(0))+coefficient*value
        if result:target[mon]=result
        else:target.pop(mon,None)


def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports/S61_AFFINE_NORMAL_QUINTIC_CONSTANT_BARRIER_20261103_v63'
    parent=ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63'
    parentreview=ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63_REVIEW.json'
    docs=[(folder/'report.json',read(folder/'report.json')),(parentreview,read(parentreview))]
    checked={}
    for path,d in docs:
        assert d['all_passed']
        for name,want in d.get('input_sha256',{}).items():assert sha(ROOT/name)==want;checked[name]=want
        for name,want in d.get('evidence_sha256',{}).items():assert sha(path.parent/name)==want;checked[str((path.parent/name).relative_to(ROOT))]=want
    accepted=docs[1][1]
    assert accepted['independent_all134_fiber_centering_and_first5_actual_plane_binding']
    assert accepted['independent_flat46symbol_sparse_all134_substitution_division_identities']
    assert accepted['generic_beta_only_Macaulay_spreading_and_two_inverse_algebra_maps_accepted']
    data=read(parent/'presentation.json.gz');certificate=read(folder/'certificate.json')
    assert len(data['coefficient_symbols'])==45 and data['coefficient_symbols'][-3:]==['S1(n)','S2(n)','H(n)']
    G=decode(data['monic_quintic_G']);assert G[5]=={():Fraction(1)}
    lower=[]
    for degree in range(5):
        coeff=dict(G.get(degree,{}))
        if degree==0:assert coeff.pop((44,))==Fraction(-1)
        assert set(coeff)<=set([()]);lower.append(coeff.get((),Fraction(0)))
    residuals=[decode(row)for row in data['all134_remainders']];assert len(residuals)==134
    rows=[]
    for original in residuals:
        assert max(original,default=0)<=4
        vector=[dict(original.get(degree,{}))for degree in range(5)]
        for j in range(5):
            rows.append(vector)
            if j==4:break
            # Independent h companion action; no producer long-division functions imported.
            top=vector[4];following=[{}for _ in range(5)]
            for degree in range(5):
                if degree:following[degree]=dict(vector[degree-1])
                add_scaled(following[degree],top,-lower[degree])
            for mon,value in top.items():
                shifted=tuple(sorted(mon+(44,)));following[0][shifted]=following[0].get(shifted,Fraction(0))+value
                if not following[0][shifted]:following[0].pop(shifted)
            vector=following
    assert len(rows)==670
    monomials=sorted({mon for row in rows for poly in row for mon in poly})
    assert ()not in monomials and all(1<=len(mon)<=2 for mon in monomials)
    quad=[(degree,mon)for degree in range(5)for mon in monomials if len(mon)==2]
    linear=[(degree,mon)for degree in range(5)for mon in monomials if len(mon)==1]
    assert len(quad)==5145 and len(linear)==225
    labels=quad+linear;index={label:j for j,label in enumerate(labels)};entries={}
    for i,row in enumerate(rows):
        for degree,poly in enumerate(row):
            for mon,coefficient in poly.items():
                assert coefficient.denominator%P
                value=coefficient.numerator*pow(coefficient.denominator,-1,P)%P
                if value:entries[i,index[degree,mon]]=value
    assert len(entries)==237811
    # Every original input entry is independently matched, including all linear columns.
    with gzip.open(folder/'input-mod32003.txt.gz','rt')as f:
        assert tuple(map(int,f.readline().split()))==(670,5145,225,237811)
        seen=set()
        for line in f:
            if not line.strip():continue
            i,j,c=map(int,line.split());assert(i,j)not in seen;seen.add((i,j));assert entries.get((i,j))==c
        assert set(entries)==seen
    selected=certificate['minor_columns_zero_based'];assert len(selected)==len(set(selected))==670
    assert certificate['minor_column_labels']==[[degree,list(mon)]for degree,mon in(quad[j]for j in selected)]
    assert certificate['nonzero_QQ_punit_minor_rows']==list(range(670))
    minor=out/'minor.txt'
    with minor.open('w')as f:
        f.write('670 32003\n')
        for i in range(670):f.write(' '.join(str(entries.get((i,j),0))for j in selected)+'\n')
    cpp=ROOT/'code/s61_affine_normal_quintic_barrier_review_20261103_v63.cpp';binary=out/'check'
    subprocess.run(['c++','-std=c++17','-O2',str(cpp),'-o',str(binary)],check=True,capture_output=True,text=True)
    run=subprocess.run([str(binary),str(minor)],check=True,capture_output=True,text=True,timeout=30)
    assert run.stdout.strip()=='INDEPENDENT_INTEGER_MINOR670_DETERMINANT 4745'
    assert certificate['minor_determinant_mod_prime']==4745
    (out/'native.log').write_text(run.stdout+run.stderr);binary.unlink()
    with minor.open('rb')as src,gzip.open(out/'minor.txt.gz','wb')as dst:
        import shutil
        shutil.copyfileobj(src,dst,1<<20)
    minor.unlink()
    # Universal plane-origin change: a beta-only h translation changes both module
    # bases by unipotent binomial matrices, hence preserves parameter-degree ranks.
    c,h,H=s.symbols('c h H');g=s.symbols('g0:5');poly=h**5+sum(g[i]*h**i for i in range(5))-H
    T=s.zeros(5)
    for j in range(4):T[j+1,j]=1
    for i in range(5):T[i,4]=-s.expand(poly).coeff(h,i)
    C=s.Matrix(5,5,lambda i,j:s.binomial(j,i)*(-c)**(j-i)if i<=j else 0)
    Ci=s.Matrix(5,5,lambda i,j:s.binomial(j,i)*c**(j-i)if i<=j else 0)
    assert(C*Ci).applyfunc(s.expand)==s.eye(5)and C.det()==1
    shifted=s.Poly(s.expand(poly.subs(h,h+c)),h);Ts=s.zeros(5)
    for j in range(4):Ts[j+1,j]=1
    for i in range(5):Ts[i,4]=-shifted.coeff_monomial(h**i)
    assert(Ci*(T-c*s.eye(5))*C-Ts).applyfunc(s.expand)==s.zeros(5)
    paths=[Path(__file__),cpp,ROOT/f'research/notes/{STEM}.md',folder/'report.json',folder/'certificate.json',folder/'input-mod32003.txt.gz',parentreview,parent/'presentation.json.gz']
    result=dict(id=STEM,all_passed=True,current_dependency_hashes_checked=len(checked),
        full_original134_binding_inherited_from_pinned_parent_independent_review=True,
        independent_companion_reconstruction_all134_times5_module_rows=True,
        all670_formal45_parameter_rows_degree_atmost2_and_no_constant=True,
        all237811_modular_entries_including225_linear_columns_independently_equal=True,
        all_QQ_row_coefficient_denominators_32003_units=True,
        selected670_column_labels_bound=True,independent_integer_Gaussian670minor_determinant=4745,
        exact_QQ_quadratic_rank=670,generic_beta_quadratic_rank=670,
        generic_spreading_is_nonzero_punit_minor_lower_bound_not_modular_columnspan_lift=True,
        universal_binomial_plane_origin_translation_and_companion_conjugation_verified=True,
        parameter_constant_purelinear_consequence_span_zero_scope='free formal45 ring K[n42,S1,S2,H], coefficients K=C(beta)',
        nonlinear_normal42_basechange_preserves_this_negative_claim=False,
        parameter_dependent_cofactor_impossibility_or_unit_or_exhaustion_claimed=False,
        no_new_Macaulay_or_GB_or_higher_cofactor_search=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={name:sha(out/name)for name in ['native.log','minor.txt.gz']},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('INDEPENDENT QUINTIC FORMAL45 CONSTANT BARRIER /670MINOR4745 PASS',result['seconds'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

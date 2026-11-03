#!/usr/bin/env python3
"""Independent full physicalQQ point-Artin5 audit, without generic-flatness inference."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix

from q_positive_weight6_small_model_review_20261103_v63 import physical_coefficients,own_action,read,qm,mm,rank,limits

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def multiply(A,words,degree):
    factors=list(itertools.combinations_with_replacement(range(2),degree))
    target=list(itertools.combinations_with_replacement(range(2),len(words[0])+degree));idx={word:i for i,word in enumerate(target)}
    out=s.zeros(A.rows*len(factors),len(target))
    for row in range(A.rows):
        for k,factor in enumerate(factors):
            for mon,c in zip(words,A.row(row)):
                if c:out[row*len(factors)+k,idx[tuple(sorted(mon+factor))]]+=c
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    packet=ROOT/'reports/Q_NEW_FAMILY_TRANSVERSE_ARTIN_20261103_v63'
    dependencies=[packet/'report.json',packet/'point111/report.json',packet/'point111-through4/report.json']
    for path in dependencies:
        r=read(path);assert r['all_passed']
        for p,digest in r.get('input_sha256',{}).items():assert sha(ROOT/p)==digest
        for p,digest in r.get('evidence_sha256',{}).items():assert sha(path.parent/p)==digest
    point=read(packet/'point111/certificate.json.gz');jet=read(packet/'point111-through4/certificate.json.gz')
    family=read(ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz')
    canonical=read(ROOT/'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63/certificate.json.gz')
    C=qm(canonical['canonical_basis']);free=canonical['canonical_free_raw_indices']
    ns=s.symbols('n0:3');subs=dict(zip(ns,[1,1,1]));expressions=list(map(s.sympify,family['canonical150']))
    q=s.Matrix([v.subs(subs)for v in expressions]);assert q==s.Matrix(list(map(s.Rational,point['physical_center150'])))
    N=own_action(q);M=N.row_join(s.Matrix([[s.diff(v,n).subs(subs)for n in ns]for v in expressions]))
    assert rank(N)==24 and rank(M)==26
    selected=M[:,point['family_image_columns']];cuts=point['family_image_cut_canonical_rows']
    minor=selected[[free[i]for i in cuts],:]
    assert DomainMatrix.from_Matrix(minor).convert_to(s.QQ).det()==s.Rational(point['family_image26_minor'])!=0
    S=C[:,point['remaining100_canonical_columns']]
    pairs,D=physical_coefficients();assert D*C==s.zeros(25,126)
    L=s.zeros(350,150)
    for row,poly in enumerate(pairs):
        value=0
        for(a,b),co in poly.items():
            value+=co*q[a]*q[b];L[row,a]+=co*q[b];L[row,b]+=co*q[a]
        assert value==0
    assert L*M==s.zeros(350,29)
    J=mm(L,S);assert J==qm(point['full350_by100_linearization'])
    rows=point['implicit98_rows'];columns=point['implicit98_columns'];inverse=qm(point['implicit98_inverse'])
    assert J[rows,columns]*inverse==s.eye(98)
    R1=qm(point['tangent150_by2']);T=R1[free,:][point['remaining100_canonical_columns'],:]
    assert T[point['normal2_columns'],:]==s.eye(2)and mm(J,T)==s.zeros(350,2)and mm(S,T)==R1
    dependent=S[:,columns];H2=qm(point['implicit98_by3_quadratic_correction'])
    Qraw=qm(point['full350_by3_uncorrected_quadratic']);F2=qm(point['full350_by3_true_quadratic_obstructions'])
    assert F2==Qraw+mm(J[:,columns],H2)and F2[rows,:]==s.zeros(98,3)and rank(F2)==1
    R2=mm(dependent,H2);R3=mm(dependent,qm(jet['implicit98_cubic']));R4=mm(dependent,qm(jet['implicit98_quartic']))
    F3=qm(jet['full350_cubic']);F4=qm(jet['full350_quartic'])
    assert F3[rows,:]==s.zeros(98,4)and F4[rows,:]==s.zeros(98,5)
    words=[list(itertools.combinations_with_replacement(range(2),d))for d in range(5)]
    native=ROOT/'code/q_positive_weight6_small_model_review_20261103_v63.cpp';binary=ROOT/'work/q-new-artin-point-review-native'
    flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
    subprocess.run(['g++','-O2','-std=c++17',str(native),*flags,'-o',str(binary)],check=True,timeout=30)
    inputpath=args.output/'physical-input.txt'
    def writepoly(stream,row,monomials):
        terms=[(word,c)for word,c in zip(monomials,row)if c];stream.write(str(len(terms))+'\n')
        for word,c in terms:stream.write(str(c)+' '+' '.join(str(word.count(i))for i in range(8))+'\n')
    with inputpath.open('w')as stream:
        for matrix,mons in zip([q,R1,R2,R3,R4],words):
            for row in matrix.tolist():writepoly(stream,row,mons)
        for row in pairs:
            stream.write(str(len(row))+'\n')
            for(a,b),co in sorted(row.items()):stream.write(f'{a} {b} {co}\n')
        for row in D.tolist():
            terms=[(i,int(c))for i,c in enumerate(row)if c];stream.write(str(len(terms))+'\n')
            for i,c in terms:stream.write(f'{i} {c}\n')
        for matrix,mons in zip([s.zeros(350,1),s.zeros(350,2),F2,F3,F4],words):
            for row in matrix.tolist():writepoly(stream,row,mons)
    with(args.output/'physical-replay.log').open('w')as log:
        subprocess.run([str(binary),str(inputpath)],stdout=log,stderr=subprocess.STDOUT,timeout=60,preexec_fn=limits,check=True)
    assert 'ALL25_DIVERGENCE_ALL5GRADES_PASS'in(args.output/'physical-replay.log').read_text()
    chosen=jet['quadratic_selected_original_row'];assert len(chosen)==1
    A2,A3,A4=F2[chosen,:],F3[chosen,:],F4[chosen,:]
    coefficients=qm(jet['constant_quadratic_combination350']);assert mm(coefficients,A2)==F2
    remaining=[i for i in range(350)if i not in chosen]
    G3=(F3-mm(coefficients,A3))[remaining,:];G4=(F4-mm(coefficients,A4))[remaining,:]
    assert G3==qm(jet['full349_order3_cubics'])and G4==qm(jet['full349_order3_quartics'])
    source=multiply(A2,words[2],1).col_join(G3)
    assert source==qm(jet['complete_degree3_multiplier_source'])and source.shape==(351,4)and rank(source)==4
    kernel=qm(jet['complete_degree3_kernel'])
    assert kernel.shape==(347,351)and rank(kernel)==347 and mm(kernel,source)==s.zeros(347,4)
    carry=mm(kernel,multiply(A3,words[3],1).col_join(G4));assert carry==qm(jet['complete_kernel_carry4'])
    u,v=s.symbols('u v');quadratic=s.expand(sum(c*s.prod([u,v][i]for i in mon)for c,mon in zip(A2.row(0),words[2])))
    discriminant=s.factor(s.Poly(quadratic,u,v).coeff_monomial(u*v)**2-4*s.Poly(quadratic,u,v).coeff_monomial(u*u)*s.Poly(quadratic,u,v).coeff_monomial(v*v))
    assert discriminant==s.Rational(read(packet/'report.json')['exact_binary_quadratic_discriminant'])!=0
    cp=list(map(s.sympify,jet['new_genuine_cubic_polynomials']))
    C3=s.Matrix([[s.Poly(p,u,v).coeff_monomial(s.prod([u,v][i]for i in word))for word in words[3]]for p in cp])
    assert C3.rows==2 and rank(multiply(A2,words[2],1).col_join(C3))==4
    assert rank(multiply(A2,words[2],2).col_join(multiply(C3,words[3],1)).col_join(carry))==5
    script='R=QQ[u,v];I=ideal('+','.join(str(p).replace('**','^')for p in[quadratic,*cp])+');G=gb I;\n'
    script+='assert(dim(R/I)==0);assert(degree(R/I)==5);assert(gens((ideal(u,v))^3)%G==0);print("EXACT_POINT_ARTIN5_M3_ZERO_PASS");\n'
    (args.output/'own-closed-fiber.m2').write_text(script);exe=shutil.which('M2');assert exe
    with(args.output/'own-closed-fiber.log').open('w')as log:
        subprocess.run([exe,'--script',str(args.output/'own-closed-fiber.m2')],stdout=log,stderr=subprocess.STDOUT,timeout=30,check=True)
    assert 'EXACT_POINT_ARTIN5_M3_ZERO_PASS'in(args.output/'own-closed-fiber.log').read_text()
    paths=[Path(__file__),native,*dependencies,packet/'point111/certificate.json.gz',packet/'point111-through4/certificate.json.gz',
        ROOT/'code/q_positive_weight6_small_model_review_20261103_v63.py',
        ROOT/'research/notes/Q_NEW_FAMILY_TRANSVERSE_ARTIN_20261103_v63_REVIEW.md']
    result=dict(id='Q_NEW_FAMILY_TRANSVERSE_ARTIN_20261103_v63_REVIEW',all_passed=True,
        actual_full126canonical_and350physical_directions=True,
        own_actual26family_cut_and_not_group_only=True,
        given98inverse_verified_not_reconstructed=True,
        native_full_QQ_two_variable_coefficients_through4_all350=True,
        complete_cubic351source_rank4_and347kernel=True,
        closed_point_fiber_HF=[1,2,2],closed_point_fiber_length=5,closed_point_maximal_ideal_cubed_zero=True,
        geometric_closed_point_fiber='C[u,v]/(u*v,u^3,v^3)',
        full_completion_finite_over26smooth_branch_parameters_generated_atmost5=True,
        no26group_product_or_generic_length_or_flatness_inference=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('INDEPENDENT POINT ARTIN5 REVIEW PASS',result['seconds'],flush=True)


if __name__=='__main__':main()

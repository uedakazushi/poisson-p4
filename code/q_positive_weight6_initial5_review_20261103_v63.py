#!/usr/bin/env python3
"""Independent all350 exact physical fifth coefficients and complete initial5 audit."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import time

import sympy as s

from q_positive_weight6_small_model_review_20261103_v63 import physical_coefficients,read,mm,rank,qm,multiply,limits,own_action

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    quadratic=ROOT/'reports/Q_POSITIVE_WEIGHT6_QUADRATIC_MODEL_20261103_v63'
    fourth=ROOT/'reports/Q_POSITIVE_WEIGHT6_INITIAL4_20261103_v63'
    fifth=ROOT/'reports/Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63'
    dependency_paths=[ROOT/'reports/Q_POSITIVE_WEIGHT6_SMALL_MODEL_20261103_v63_REVIEW/report.json',fifth/'report.json']
    for path in dependency_paths:
        r=read(path);assert r['all_passed']
        for p,digest in r.get('input_sha256',{}).items():assert sha(ROOT/p)==digest
        for p,digest in r.get('evidence_sha256',{}).items():assert sha(path.parent/p)==digest
    Q=read(quadratic/'certificate.json.gz');F4=read(fourth/'certificate.json.gz')
    J5=read(fifth/'fifth-jets.json.gz');F5=read(fifth/'fifth-complete-spanning.json.gz')
    selected=read(fifth/'selected-fifth.json')
    canonical=read(ROOT/Q['canonical_basis_dependency']);C=qm(canonical['canonical_basis'])
    dependent=C[:,Q['remaining103_canonical_indices']][:,Q['implicit_coordinate_columns']]
    R1=qm(Q['tangent150_by8']);R2=qm(F4['raw150_homogeneous_implicit2']);R3=qm(F4['raw150_homogeneous_implicit3'])
    R4=mm(dependent,qm(F4['dependent95_homogeneous_implicit4']));R5=mm(dependent,qm(J5['dependent95_fifth_coefficients']))
    q=s.Matrix([s.Rational(v)for v in Q['physical_center150']])
    # Explicit projective stabilizer: two commuting square-zero directions,
    # plus the GL center.  No connected torus can grade this fixed point.
    U=s.zeros(5);V=s.zeros(5);U[3,0]=1;V[2,0]=s.Rational(4,3);V[3,1]=1
    stabilizer=s.Matrix.hstack(*[s.Matrix(list(A)+[0])for A in[U,V,s.eye(5)]])
    action=own_action(q);assert rank(action)==23 and rank(stabilizer)==3
    assert action*stabilizer==s.zeros(150,3)
    assert U*U==V*V==U*V==V*U==s.zeros(5)
    Q2=qm(Q['full350_by36_corrected_obstructions']);Q3=qm(F4['full350_corrected_cubic'])
    Q4=qm(F4['full350_corrected_quartic']);Q5=qm(J5['full350_fifth_residual'])
    words=[list(itertools.combinations_with_replacement(range(8),degree))for degree in range(6)]
    assert words[5]==[tuple(m)for m in J5['normal_monomials5']]
    assert Q5[Q['implicit_equation_rows'],:]==s.zeros(95,792)
    pairs,D=physical_coefficients()
    native=ROOT/'code/q_positive_weight6_initial5_review_20261103_v63.cpp';binary=ROOT/'work/q-positive-initial5-review-native'
    flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
    subprocess.run(['g++','-O2','-std=c++17',str(native),*flags,'-o',str(binary)],check=True,timeout=30)
    inputpath=args.output/'physical-input.txt'
    def writepoly(stream,row,monomials):
        nonzero=[(word,c)for word,c in zip(monomials,row)if c];stream.write(str(len(nonzero))+'\n')
        for word,c in nonzero:stream.write(str(c)+' '+' '.join(str(word.count(i))for i in range(8))+'\n')
    with inputpath.open('w')as stream:
        for matrix,mons in zip([q,R1,R2,R3,R4,R5],words):
            for row in matrix.tolist():writepoly(stream,row,mons)
        for row in pairs:
            stream.write(str(len(row))+'\n')
            for(a,b),c in sorted(row.items()):stream.write(f'{a} {b} {c}\n')
        for row in D.tolist():
            terms=[(i,int(c))for i,c in enumerate(row)if c];stream.write(str(len(terms))+'\n')
            for i,c in terms:stream.write(f'{i} {c}\n')
        for matrix,mons in zip([s.zeros(350,1),s.zeros(350,8),Q2,Q3,Q4,Q5],words):
            for row in matrix.tolist():writepoly(stream,row,mons)
    with(args.output/'physical-replay.log').open('w')as log:
        subprocess.run([str(binary),str(inputpath)],stdout=log,stderr=subprocess.STDOUT,timeout=120,preexec_fn=limits,check=True)
    assert 'ALL25_DIVERGENCE_ALL6GRADES_PASS'in(args.output/'physical-replay.log').read_text()
    print('OWN ALL350 PHYSICAL QQ ALLDEGREE5 COEFFICIENTS PASS',time.monotonic()-started,flush=True)
    # Exact polynomial row operations isolate the minimal2,3,4 generators.
    chosen=F4['independent_quadratic_original_rows'];remaining=F4['remaining345_original_rows']
    coefficients=qm(F4['quadratic_combination350_by5']);A2=Q2[chosen,:];A3=Q3[chosen,:]
    A4=Q4[chosen,:];A5=Q5[chosen,:]
    G3=qm(F4['order3_generators_cubic345']);G4=qm(F4['order3_generators_quartic345'])
    G5=(Q5-mm(coefficients,A5))[remaining,:]
    ci=F5['cubic_G_indices'];C3,C4,C5=G3[ci,:],G4[ci,:],G5[ci,:]
    M3=multiply(A2,words[2],1)[0].col_join(C3);assert M3.shape==(48,120)and rank(M3)==48
    HC=qm(F5['G_cubic_representation345_by48']);assert mm(HC,M3)==G3
    actual4=multiply(A3,words[3],1)[0].col_join(C4)
    actual5=multiply(A4,words[4],1)[0].col_join(C5)
    hi=F5['remaining_H_indices'];assert hi==[i for i in range(345)if i not in ci]
    H4=(G4-mm(HC,actual4))[hi,:];H5=(G5-mm(HC,actual5))[hi,:]
    assert H4==qm(F5['quartic_H337'])and H5==qm(F5['fifth_H337'])
    di=F5['quartic_H_indices'];D4,D5=H4[di,:],H5[di,:]
    M4=multiply(A2,words[2],2)[0].col_join(multiply(C3,words[3],1)[0]).col_join(D4)
    assert M4.shape==(252,330)and M4==qm(F5['minimal_quartic_source252'])and rank(M4)==204
    K4=qm(F5['full_quartic_kernel48']);assert K4.shape==(48,252)and rank(K4)==48
    assert mm(K4,M4)==s.zeros(48,330)
    actual5=multiply(A3,words[3],2)[0].col_join(multiply(C4,words[4],1)[0]).col_join(D5)
    assert actual5==qm(F5['quartic_source_fifth252'])
    carry=mm(K4,actual5);assert carry==qm(F5['complete_quartic_kernel_fifth_carries48'])
    independent=F5['quartic_independent_source_rows'];assert len(independent)==204 and rank(M4[independent,:])==204
    EC=qm(F5['H_quartic_representation337_by204']);assert mm(EC,M4[independent,:])==H4
    E5=H5-mm(EC,actual5[independent,:]);assert E5[di,:]==s.zeros(8,792)
    ei=F5['remaining_E_indices'];assert ei==[i for i in range(337)if i not in di]
    E5=E5[ei,:];assert E5==qm(F5['complete_fifth_order_new329'])
    lower=multiply(A2,words[2],3)[0].col_join(multiply(C3,words[3],2)[0]).col_join(multiply(D4,words[4],1)[0])
    full=lower.col_join(E5).col_join(carry);assert full.shape==(1329,792)
    assert full==qm(F5['complete_degree5_spanning1329'])
    z=s.symbols('z0:8');index={word:i for i,word in enumerate(words[5])}
    def vector(expression):
        poly=s.Poly(s.sympify(expression),*z,domain=s.QQ);out=[s.S.Zero]*792
        for exponent,c in poly.terms():
            word=tuple(i for i,e in enumerate(exponent)for _ in range(e));assert len(word)==5
            out[index[word]]=c
        return out
    selected_rows=selected['selected_full_spanning_rows'];assert len(selected_rows)==15
    assert s.Matrix([vector(p)for p in selected['selected_QQ_fifth_polynomials']])==full[selected_rows,:]
    print('OWN COMPLETE48 CARRY +329 ORDER5 ROWS PASS',time.monotonic()-started,flush=True)
    # All377 exact reductions supply the QQ upper bound; no modular rank upper.
    def polynomial(row,monomials):
        return '+'.join('('+str(c)+')*'+'*'.join('z'+str(i)for i in word)for c,word in zip(row,monomials)if c)or'0'
    lowerpolys=[polynomial(row,mons)for matrix,mons in[(A2,words[2]),(C3,words[3]),(D4,words[4])]for row in matrix.tolist()]
    selectedpolys=[polynomial(row,words[5])for row in full[selected_rows,:].tolist()]
    allnew=[polynomial(row,words[5])for row in E5.col_join(carry).tolist()]
    m2='R=QQ['+','.join(map(str,z))+'];\nI=ideal('+','.join(p.replace('**','^')for p in lowerpolys+selectedpolys)+');G=gb I;\n'
    m2+='ALL=matrix{{'+','.join(allnew)+'}};assert(ALL%G==0);\n'
    m2+='assert(hilbertFunction(5,R/I)==171);assert(dim(R/I)==3);assert(degree(R/I)==5);\n'
    m2+='print("ALL377_EXACT_QQ_REDUCTIONS_AND_HF5_171_PASS");\n'
    m2+='A=QQ[u0,u2,u3,u4,u6];fa=map(A,R,{u0,0,u2,u3,u4,0,u6,0});IA=fa I;RA=radical IA;KA=ideal(u6,u4,u2*u3);\n'
    m2+='assert(gens RA%gb KA==0);assert(gens KA%gb RA==0);\n'
    m2+='B=QQ[v0,v2,v3,v5,v6];fb=map(B,R,{v0,0,v2,v3,-2*v6,v5,v6,0});IB=fb I;RB=radical IB;KB=ideal(v6,v2*v3,v2*v5);\n'
    m2+='assert(gens RB%gb KB==0);assert(gens KB%gb RB==0);\n'
    m2+='print("TWO_QUADRATIC_PLANES_RADICAL_SPLIT_PASS");\n'
    (args.output/'own-fifth-cone.m2').write_text(m2);exe=shutil.which('M2');assert exe
    with(args.output/'own-fifth-cone.log').open('w')as log:
        subprocess.run([exe,'--script',str(args.output/'own-fifth-cone.m2')],stdout=log,stderr=subprocess.STDOUT,timeout=60,preexec_fn=limits,check=True)
    assert 'TWO_QUADRATIC_PLANES_RADICAL_SPLIT_PASS'in(args.output/'own-fifth-cone.log').read_text()
    paths=[Path(__file__),native,*dependency_paths,fifth/'fifth-jets.json.gz',fifth/'fifth-complete-spanning.json.gz',
        ROOT/'code/q_positive_weight6_small_model_review_20261103_v63.py',
        ROOT/'research/notes/Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63_REVIEW.md']
    report=dict(id='Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63_REVIEW',all_passed=True,
        independent_native_full350_QQ_8variable_coefficients_through5=True,
        independent_physical_divergence_all25_through5=True,
        all_polynomial_row_eliminations_reconstructed=True,
        full48_kernel_of252rank204_verified=True,
        complete_fifth_initial1329_by792_sources_verified=True,
        all377_reductions_exact_QQ_and_reverse_subset_checked=True,
        exact_QQ_degree5_initial_rank=621,actual_HF_through5=[1,8,31,72,126,171],
        exact_fifth_cone_dimension3_degree5=True,
        radical_two_linear_components_dimensions=[3,2],
        radical='(z1,z4,z6,z7,z2*z3,z2*z5)',
        projective_connected_stabilizer_Ga2_verified=True,
        stabilizer_generators=['E30','(4/3)*E20+E31'],
        actual_projective_dimension_upper25_only=True,
        actual_dimension_equality_or_branches_or_fullinitial_claimed=False,
        no115second_producer_rerun=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('INDEPENDENT FIFTH REVIEW PASS',report['seconds'],flush=True)


if __name__=='__main__':main()

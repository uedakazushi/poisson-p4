#!/usr/bin/env python3
"""Independent physical, formal8, and complete initial4 review for the positive ray."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import shutil
import subprocess
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(5),2));MON2=list(itertools.combinations_with_replacement(range(5),2))
MON3=list(itertools.combinations_with_replacement(range(5),3));TRIPLES=list(itertools.combinations(range(5),3))
BASIS=list(itertools.product(PAIRS,MON2));BI={v:i for i,v in enumerate(BASIS)}


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(gzip.decompress(path.read_bytes())if path.suffix=='.gz'else path.read_bytes())
def mm(A,B):return (DomainMatrix.from_Matrix(A).convert_to(s.QQ)*DomainMatrix.from_Matrix(B).convert_to(s.QQ)).to_Matrix()
def rank(A):return len(DomainMatrix.from_Matrix(A).rref()[1])
def qm(rows):return s.Matrix([[s.Rational(c)for c in row]for row in rows])


def physical_coefficients():
    """Build physical x5 differentiation into constant raw150 coefficient pairs."""
    result=[defaultdict(int)for _ in range(350)];position={m:i for i,m in enumerate(MON3)}
    for rr,triple in enumerate(TRIPLES):
        for a,b,c in[triple,(triple[1],triple[2],triple[0]),(triple[2],triple[0],triple[1])]:
            for ell in range(5):
                if a==ell or b==c:continue
                left=(min(a,ell),max(a,ell));right=(min(b,c),max(b,c))
                sign=(1 if a<ell else-1)*(1 if b<c else-1)
                for mon in MON2:
                    i=BI[left,mon]
                    for other in MON2:
                        if ell not in other:continue
                        rest=list(other);rest.remove(ell);j=BI[right,other]
                        row=rr*35+position[tuple(sorted(mon+tuple(rest)))]
                        result[row][tuple(sorted((i,j)))]+=sign*other.count(ell)
    result=[{p:c for p,c in row.items()if c}for row in result]
    D=s.zeros(25,150)
    for col,((i,j),mon)in enumerate(BASIS):
        for direction,component,sign in[(i,j,1),(j,i,-1)]:
            if direction in mon:
                remaining=list(mon);remaining.remove(direction)
                D[component*5+remaining[0],col]+=sign*mon.count(direction)
    return result,D


def own_action(q):
    result=s.zeros(150,26)
    def coefficient(i,j,mon):return q[BI[(i,j),mon]]if i<j else(-q[BI[(j,i),mon]]if j<i else 0)
    for a in range(5):
        for b in range(5):
            col=5*a+b
            for row,((i,j),mon)in enumerate(BASIS):
                if i==a:result[row,col]+=coefficient(b,j,mon)
                if j==a:result[row,col]+=coefficient(i,b,mon)
                for occurrence,k in enumerate(mon):
                    if k==a:result[BI[(i,j),tuple(sorted((mon[1-occurrence],b)))],col]-=q[row]
    result[:,25]=q
    return result


def multiply(A,words,degree):
    factors=list(itertools.combinations_with_replacement(range(8),degree))
    target=list(itertools.combinations_with_replacement(range(8),len(words[0])+degree));idx={m:i for i,m in enumerate(target)}
    out=s.zeros(A.rows*len(factors),len(target))
    for row in range(A.rows):
        coefficients=[(mon,A[row,j])for j,mon in enumerate(words)if A[row,j]]
        for k,factor in enumerate(factors):
            for mon,c in coefficients:out[row*len(factors)+k,idx[tuple(sorted(mon+factor))]]+=c
    return out,target


def limits():resource.setrlimit(resource.RLIMIT_AS,(2200*1024*1024,2200*1024*1024))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    quadratic=ROOT/'reports/Q_POSITIVE_WEIGHT6_QUADRATIC_MODEL_20261103_v63'
    initial=ROOT/'reports/Q_POSITIVE_WEIGHT6_INITIAL4_20261103_v63'
    dependencies=[]
    for folder in[quadratic,initial]:
        r=read(folder/'report.json');assert r['all_passed']
        for path,digest in r['input_sha256'].items():assert sha(ROOT/path)==digest
        for path,digest in r.get('evidence_sha256',{}).items():assert sha(folder/path)==digest
        assert sha(folder/'certificate.json.gz')==r['certificate_sha256']
        dependencies.append(folder/'report.json')
    Q=read(quadratic/'certificate.json.gz');T4=read(initial/'certificate.json.gz')
    canonical=read(ROOT/Q['canonical_basis_dependency']);C=qm(canonical['canonical_basis']);free=Q['canonical_free_raw_indices']
    pairs,D=physical_coefficients();assert D.shape==(25,150)and rank(D)==24
    assert C.shape==(150,126)and C[free,:]==s.eye(126)and D*C==s.zeros(25,126)
    q=s.Matrix([s.Rational(c)for c in Q['physical_center150']]);N=own_action(q)
    for index,value in enumerate(q):
        if value:
            pair,mon=BASIS[index];power=sum(pair)-sum(mon)-5
            assert power==(-6 if index in[64,83,116]else 0)
    oldN=read(ROOT/'reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63/certificate.json')
    assert N==qm(oldN['GL25_plus_scalar_matrix'])and D*N==s.zeros(25,26)
    assert rank(N)==23 and rank(N[:,:25])==22
    cuts=Q['orbit_cut_canonical_indices'];ocols=Q['orbit_columns_zero_based'];O=N[:,ocols]
    assert O[[free[i]for i in cuts],:].det()==s.Rational(Q['orbit_cut_minor'])!=0
    S=C[:,Q['remaining103_canonical_indices']];L=s.zeros(350,150)
    for row,poly in enumerate(pairs):
        value=0
        for(a,b),co in poly.items():
            value+=co*q[a]*q[b];L[row,a]+=co*q[b];L[row,b]+=co*q[a]
        assert value==0
    assert L*N==s.zeros(350,26)
    J=mm(L,S);assert J==qm(Q['full350_by103_linear_Jacobian'])
    rows=Q['implicit_equation_rows'];cols=Q['implicit_coordinate_columns'];normal=Q['normal8_coordinate_columns']
    B=J[rows,cols];inverse=qm(Q['implicit_inverse95']);assert B*inverse==s.eye(95)
    raw1=qm(Q['tangent150_by8']);kernel=raw1[free,:][Q['remaining103_canonical_indices'],:]
    assert kernel[normal,:]==s.eye(8)and J*kernel==s.zeros(350,8)and S*kernel==raw1
    assert raw1[[free[i]for i in cuts],:]==s.zeros(23,8)
    H2=qm(Q['implicit95_by36_second_order_coefficients']);Qraw=qm(Q['full350_by36_raw_quadratic'])
    Q2=qm(Q['full350_by36_corrected_obstructions']);JD=J[:,cols]
    assert B*H2==-Qraw[rows,:]and Q2==Qraw+mm(JD,H2)
    assert rank(Q2)==5
    words=[[( )],[(i,)for i in range(8)],list(itertools.combinations_with_replacement(range(8),2)),
        list(itertools.combinations_with_replacement(range(8),3)),list(itertools.combinations_with_replacement(range(8),4))]
    raw2=mm(S[:,cols],H2);raw3=mm(S[:,cols],qm(T4['dependent95_homogeneous_implicit3']))
    raw4=mm(S[:,cols],qm(T4['dependent95_homogeneous_implicit4']))
    assert raw2==qm(T4['raw150_homogeneous_implicit2'])and raw3==qm(T4['raw150_homogeneous_implicit3'])
    Q3=qm(T4['full350_corrected_cubic']);Q4=qm(T4['full350_corrected_quartic'])
    assert Q3[rows,:]==s.zeros(95,120)and Q4[rows,:]==s.zeros(95,330)
    # Native review uses only physical x5 exponent coefficients generated above.
    native=ROOT/'code/q_positive_weight6_small_model_review_20261103_v63.cpp';binary=ROOT/'work/q-positive-small-review-native-final'
    flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
    subprocess.run(['g++','-O2','-std=c++17',str(native),*flags,'-o',str(binary)],check=True,timeout=30)
    inputpath=args.output/'physical-input.txt'
    def writepoly(stream,row,monomials):
        nonzero=[(word,c)for word,c in zip(monomials,row)if c];stream.write(str(len(nonzero))+'\n')
        for word,c in nonzero:stream.write(str(c)+' '+' '.join(str(word.count(i))for i in range(8))+'\n')
    with inputpath.open('w')as stream:
        for matrix,mons in zip([q,raw1,raw2,raw3,raw4],words):
            for row in matrix.tolist():writepoly(stream,row,mons)
        for row in pairs:
            stream.write(str(len(row))+'\n')
            for(a,b),co in sorted(row.items()):stream.write(f'{a} {b} {co}\n')
        for row in D.tolist():
            entries=[(i,int(c))for i,c in enumerate(row)if c];stream.write(str(len(entries))+'\n')
            for i,c in entries:stream.write(f'{i} {c}\n')
        for matrix,mons in zip([s.zeros(350,1),s.zeros(350,8),Q2,Q3,Q4],words):
            for row in matrix.tolist():writepoly(stream,row,mons)
    with(args.output/'physical-replay.log').open('w')as log:
        subprocess.run([str(binary),str(inputpath)],stdout=log,stderr=subprocess.STDOUT,timeout=120,preexec_fn=limits,check=True)
    assert 'ALL25_DIVERGENCE_ALL5GRADES_PASS'in(args.output/'physical-replay.log').read_text()
    print('OWN FULL PHYSICAL350 ALLQQ COEFFICIENTS THROUGH4 PASS',time.monotonic()-start,flush=True)
    # Complete degree3/4 cancellation audit, without importing producer routines.
    chosen=T4['independent_quadratic_original_rows'];coeff=qm(T4['quadratic_combination350_by5'])
    A2,A3,A4=Q2[chosen,:],Q3[chosen,:],Q4[chosen,:]
    assert len(chosen)==rank(A2)==5 and mm(coeff,A2)==Q2
    rem=T4['remaining345_original_rows'];assert rem==[i for i in range(350)if i not in chosen]
    G3=(Q3-mm(coeff,A3))[rem,:];G4=(Q4-mm(coeff,A4))[rem,:]
    assert G3==qm(T4['order3_generators_cubic345'])and G4==qm(T4['order3_generators_quartic345'])
    W3,_=multiply(A2,words[2],1);all3=W3.col_join(G3)
    assert all3==qm(T4['complete_degree3_source_matrix'])and rank(all3)==48
    K=qm(T4['complete_degree3_kernel']);assert K.shape==(337,385)and rank(K)==337
    assert mm(K,all3)==s.zeros(337,120)
    source4=multiply(A3,words[3],1)[0].col_join(G4);carry=mm(K,source4)
    assert carry==qm(T4['complete_degree3_kernel_carry4'])
    cubic=G3[T4['cubic_minimal_row_indices'],:];assert cubic.rows==8 and rank(W3)==40
    lower4=multiply(A2,words[2],2)[0].col_join(multiply(cubic,words[3],1)[0])
    assert rank(lower4)==196 and rank(lower4.col_join(carry))==204
    quartic=carry[T4['quartic_minimal_carry_indices'],:];assert quartic.rows==8 and rank(lower4.col_join(quartic))==204
    executable=shutil.which('M2');assert executable
    for folder,name in[(quadratic,'quadrics'),(initial,'initial4')]:
        with(args.output/(name+'.log')).open('w')as log:
            subprocess.run([executable,'--script',str(folder/(name+'.m2'))],stdout=log,stderr=subprocess.STDOUT,timeout=30,check=True)
    assert 'EXACT_RADICAL_TWO_FIVE_PLANES_PASS'in(args.output/'quadrics.log').read_text()
    m2=(args.output/'initial4.log').read_text();assert 'DIM4=3'in m2 and'DEG4=6'in m2
    inputs=[Path(__file__),native,*dependencies,quadratic/'certificate.json.gz',initial/'certificate.json.gz',
        ROOT/'research/notes/Q_POSITIVE_WEIGHT6_SMALL_MODEL_20261103_v63_REVIEW.md']
    report=dict(id='Q_POSITIVE_WEIGHT6_SMALL_MODEL_20261103_v63_REVIEW',all_passed=True,
        own_canonical_divergence_rank24_kernel126=True,own_full_GL25_plus_scalar_action=True,
        exact_orbit23_and_projective22_factor=True,own_physical350_by150_linearization=True,
        exact_implicit95_and_formal8_presentation=True,
        independent_native_all350_QQ_coefficients_degrees0_through4=True,
        independent_all25_divergence_coefficients_degrees0_through4=True,
        complete_QQ_degree3_kernel_rank337=True,all_degree4_carries_recomputed=True,
        actual_HF_through4=[1,8,31,72,126],new_cubic_count=8,new_quartic_count=8,
        exact_small_M2_radical_and_dim3_degree6_replayed=True,
        actual_transverse_dimension_upper3_projective_upper25_only=True,
        actual_dimension_equality_or_allbranches_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('INDEPENDENT SMALL MODEL REVIEW PASS',report['seconds'],flush=True)


if __name__=='__main__':main()

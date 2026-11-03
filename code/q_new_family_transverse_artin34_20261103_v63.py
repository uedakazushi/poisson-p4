#!/usr/bin/env python3
"""Exact closed two-normal fiber through4 at the saved26-cut new-family point."""
import argparse
import itertools
from pathlib import Path
import json
import shutil
import subprocess
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from q_positive_weight6_quadratic_model_20261103_v63 import ROOT,read,matrix,sha,pivot_columns
from q_positive_weight6_initial4_20261103_v63 import mm,bilinear,multiply_words,write_gzip


def multiply2(A,words,degree):
    full,allwords=multiply_words(A,words,degree)
    wanted=[i for i,word in enumerate(allwords)if set(word)<={0,1}]
    # multiply_words makes8 factors; keep only actual2-variable factors too.
    factorwords=list(itertools.combinations_with_replacement(range(8),degree))
    keep=[row*len(factorwords)+i for row in range(A.rows)for i,word in enumerate(factorwords)if set(word)<={0,1}]
    return full[keep,wanted],[allwords[i]for i in wanted]


def pair(U,mu,V,mv):
    value,words=bilinear(U,mu,V,mv)
    keep=[i for i,word in enumerate(words)if set(word)<={0,1}]
    return value[:,keep],[words[i]for i in keep]


def main(source,output):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    d=read(source/'certificate.json.gz');qreport=read(source/'report.json');assert qreport['all_passed']
    C=s.Matrix(read(ROOT/'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63/certificate.json.gz')['canonical_basis'])
    D=C[:,d['remaining100_canonical_columns']][:,d['implicit98_columns']]
    B=s.Matrix(d['implicit98_inverse']);J=s.Matrix(d['full350_by100_linearization'])[:,d['implicit98_columns']]
    rows=d['implicit98_rows'];R1=s.Matrix(d['tangent150_by2']);mu1=[(0,),(1,)]
    H2=s.Matrix(d['implicit98_by3_quadratic_correction']);mu2=[(0,0),(0,1),(1,1)];R2=mm(D,H2)
    a,mu3=pair(R1,mu1,R2,mu2);b,check=pair(R2,mu2,R1,mu1);assert check==mu3
    raw3=a+b;H3=-mm(B,raw3[rows,:]);F3=raw3+mm(J,H3);assert F3[rows,:]==s.zeros(98,4)
    R3=mm(D,H3)
    a,mu4=pair(R1,mu1,R3,mu3);b,check=pair(R2,mu2,R2,mu2);assert check==mu4
    c,check=pair(R3,mu3,R1,mu1);assert check==mu4
    raw4=a+b+c;H4=-mm(B,raw4[rows,:]);F4=raw4+mm(J,H4);assert F4[rows,:]==s.zeros(98,5)
    F2=s.Matrix(d['full350_by3_true_quadratic_obstructions'])
    chosen=pivot_columns(F2.T);assert len(chosen)==1
    A2=F2[chosen,:];A3=F3[chosen,:];A4=F4[chosen,:];qc=pivot_columns(A2)
    coefficients=mm(F2[:,qc],A2[:,qc].inv());assert mm(coefficients,A2)==F2
    remaining=[i for i in range(350)if i not in chosen]
    G3=(F3-mm(coefficients,A3))[remaining,:];G4=(F4-mm(coefficients,A4))[remaining,:]
    W3,check=multiply2(A2,mu2,1);assert check==mu3
    source3=W3.col_join(G3);piv3=pivot_columns(source3.T);r3=len(piv3)
    ci=[i-W3.rows for i in piv3 if i>=W3.rows];C3=G3[ci,:]
    kernel=DomainMatrix.from_Matrix(source3.T).nullspace().to_Matrix()
    assert mm(kernel,source3)==s.zeros(kernel.rows,4)
    high3,check=multiply2(A3,mu3,1);assert check==mu4
    carry=mm(kernel,high3.col_join(G4))
    low4,check=multiply2(A2,mu2,2);assert check==mu4
    cx,check=multiply2(C3,mu3,1);assert check==mu4
    all4=low4.col_join(cx).col_join(carry);r4=len(pivot_columns(all4.T))
    u,v=s.symbols('u v')
    def polynomials(M,words):return [str(s.factor(sum(coefficient*s.prod([u,v][i]for i in word)for coefficient,word in zip(row,words)if coefficient)))for row in M.tolist()]
    qp=polynomials(A2,mu2);cp=polynomials(C3,mu3)
    text='R=QQ[u,v];\nI=ideal('+','.join(x.replace('**','^')for x in qp+cp)+');\nG=gb I;\n'
    text+='print("DIM="|toString dim(R/I));print("DEG="|toString degree(R/I));print("HS="|toString hilbertSeries(R/I));\nprint("GB="|toString gens G);\n'
    (output/'closed-fiber.m2').write_text(text)
    executable=shutil.which('M2');assert executable
    with (output/'closed-fiber.log').open('w')as log:subprocess.run([executable,'--script',str(output/'closed-fiber.m2')],stdout=log,stderr=subprocess.STDOUT,timeout=25,check=True)
    hf=[1,2,2,4-r3,5-r4];length=sum(hf)if r4==5 else None
    payload=dict(point=d['point'],quadratic_selected_original_row=chosen,normal_monomials3=mu3,normal_monomials4=mu4,
        full350_cubic=matrix(F3),full350_quartic=matrix(F4),implicit98_cubic=matrix(H3),implicit98_quartic=matrix(H4),
        constant_quadratic_combination350=matrix(coefficients),full349_order3_cubics=matrix(G3),full349_order3_quartics=matrix(G4),
        complete_degree3_multiplier_source=matrix(source3),complete_degree3_kernel=matrix(kernel),complete_kernel_carry4=matrix(carry),
        quadratic_polynomials=qp,new_genuine_cubic_polynomials=cp)
    write_gzip(output/'certificate.json.gz',payload)
    paths=[Path(__file__),source/'report.json',source/'certificate.json.gz',
        ROOT/'code/q_positive_weight6_initial4_20261103_v63.py',ROOT/'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63/certificate.json.gz']
    report=dict(id='Q_NEW_FAMILY_TRANSVERSE_ARTIN34_20261103_v63',all_passed=True,point=d['point'],
        full350_physical_closed_fiber_orders3_and4=True,all98_implicit_rows_zero=True,
        complete_cubic_initial_rank=r3,complete_quartic_initial_rank=r4,actual_closed_fiber_HF_through4=hf,
        actual_closed_fiber_length_exact=length,actual_closed_fiber_maximal_ideal_fourth_power_zero=r4==5,
        quadratic_polynomials=qp,new_genuine_cubic_polynomials=cp,
        actual_completion_finite_over26_smooth_branch_parameters=r4==5,
        generic_flatness_or_generic_length_or_generic_nonreducedness_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={p.name:sha(p)for p in output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('CLOSED NEW-C TRANSVERSE2 TRUEHF',hf,'LENGTH',length,'SEC',report['seconds'],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();main(args.source,args.output)

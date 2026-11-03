#!/usr/bin/env python3
"""Complete exact fifth initial spanning data using48 minimal fourth syzygies."""
import argparse
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from q_positive_weight6_quadratic_model_20261103_v63 import ROOT,read,matrix,sha,pivot_columns
from q_positive_weight6_initial4_20261103_v63 import mm,bilinear,multiply_words,write_gzip


def main(output,rank_program,run_m2=False):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    qfolder=ROOT/'reports/Q_POSITIVE_WEIGHT6_QUADRATIC_MODEL_20261103_v63'
    ffolder=ROOT/'reports/Q_POSITIVE_WEIGHT6_INITIAL4_20261103_v63'
    q=read(qfolder/'certificate.json.gz');f=read(ffolder/'certificate.json.gz')
    for folder in [qfolder,ffolder]:
        report=read(folder/'report.json');assert report['all_passed']
        for name,digest in report['input_sha256'].items():assert sha(ROOT/name)==digest
    base=read(ROOT/q['canonical_basis_dependency']);C=s.Matrix(base['canonical_basis'])
    D=C[:,q['remaining103_canonical_indices']][:,q['implicit_coordinate_columns']]
    B=s.Matrix(q['implicit_inverse95']);J=s.Matrix(q['full350_by103_linear_Jacobian'])[:,q['implicit_coordinate_columns']]
    R1=s.Matrix(q['tangent150_by8']);m1=[(i,)for i in range(8)]
    R2=s.Matrix(f['raw150_homogeneous_implicit2']);m2=[tuple(w)for w in q['normal_quadratic_monomials']]
    R3=s.Matrix(f['raw150_homogeneous_implicit3']);m3=[tuple(w)for w in f['normal_monomials3']]
    R4=mm(D,s.Matrix(f['dependent95_homogeneous_implicit4']));m4=[tuple(w)for w in f['normal_monomials4']]
    terms=[]
    for U,mu,V,mv in [(R1,m1,R4,m4),(R2,m2,R3,m3),(R3,m3,R2,m2),(R4,m4,R1,m1)]:
        value,m5=bilinear(U,mu,V,mv);terms.append(value)
    raw5=sum(terms,s.zeros(350,792));H5=-mm(B,raw5[q['implicit_equation_rows'],:]);F5=raw5+mm(J,H5)
    assert F5[q['implicit_equation_rows'],:]==s.zeros(95,792)
    write_gzip(output/'fifth-jets.json.gz',dict(normal_monomials5=m5,full350_fifth_residual=matrix(F5),dependent95_fifth_coefficients=matrix(H5)))
    print('SAVED all350 exact fifth jets',time.monotonic()-started,flush=True)
    chosen=f['independent_quadratic_original_rows'];remaining=f['remaining345_original_rows']
    coefficients=s.Matrix(f['quadratic_combination350_by5'])
    A2=s.Matrix(q['full350_by36_corrected_obstructions'])[chosen,:]
    A3=s.Matrix(f['full350_corrected_cubic'])[chosen,:]
    A4=s.Matrix(f['full350_corrected_quartic'])[chosen,:];A5=F5[chosen,:]
    G3=s.Matrix(f['order3_generators_cubic345']);G4=s.Matrix(f['order3_generators_quartic345'])
    G5=(F5-mm(coefficients,A5))[remaining,:]
    ci=f['cubic_minimal_row_indices'];C3=G3[ci,:];C4=G4[ci,:];C5=G5[ci,:]
    W3,_=multiply_words(A2,m2,1);M3=W3.col_join(C3)
    pc3=pivot_columns(M3);assert len(pc3)==48
    inverse3=DomainMatrix.from_Matrix(M3[:,pc3]).inv().to_Matrix()
    HC=mm(G3[:,pc3],inverse3);assert mm(HC,M3)==G3
    L4,_=multiply_words(A3,m3,1);L4=L4.col_join(C4)
    L5,_=multiply_words(A4,m4,1);L5=L5.col_join(C5)
    H4=G4-mm(HC,L4);HH5=G5-mm(HC,L5)
    assert H4[ci,:]==s.zeros(8,330) and HH5[ci,:]==s.zeros(8,792)
    hi=[i for i in range(345)if i not in ci];H4=H4[hi,:];HH5=HH5[hi,:]
    lower4,_=multiply_words(A2,m2,2);cx4,_=multiply_words(C3,m3,1);lower4=lower4.col_join(cx4)
    p4=pivot_columns(lower4.col_join(H4).T)
    di=[i-lower4.rows for i in p4 if i>=lower4.rows];assert len(di)==8
    D4=H4[di,:];D5=HH5[di,:];M4=lower4.col_join(D4)
    independent=pivot_columns(M4.T);assert len(independent)==204
    K4=DomainMatrix.from_Matrix(M4.T).nullspace().to_Matrix()
    assert K4.rows==48 and mm(K4,M4)==s.zeros(48,330)
    A_actual5,_=multiply_words(A3,m3,2);C_actual5,_=multiply_words(C4,m4,1)
    actual5=A_actual5.col_join(C_actual5).col_join(D5)
    carry=mm(K4,actual5)
    basis=M4[independent,:];pc4=pivot_columns(basis);assert len(pc4)==204
    inverse4=DomainMatrix.from_Matrix(basis[:,pc4]).inv().to_Matrix()
    EC=mm(H4[:,pc4],inverse4);assert mm(EC,basis)==H4
    E5=HH5-mm(EC,actual5[independent,:]);assert E5[di,:]==s.zeros(8,792)
    ei=[i for i in range(337)if i not in di];E5=E5[ei,:]
    low5,_=multiply_words(A2,m2,3);cx5,_=multiply_words(C3,m3,2);dx5,_=multiply_words(D4,m4,1)
    low5=low5.col_join(cx5).col_join(dx5);assert low5.rows==952
    full5=low5.col_join(E5).col_join(carry);assert full5.shape==(1329,792)
    packet=dict(normal_monomials5=m5,quadratic_original_rows=chosen,cubic_G_indices=ci,remaining_H_indices=hi,
        G_cubic_representation345_by48=matrix(HC),quartic_H337=matrix(H4),fifth_H337=matrix(HH5),
        quartic_H_indices=di,minimal_quartic_source252=matrix(M4),full_quartic_kernel48=matrix(K4),
        quartic_source_fifth252=matrix(actual5),complete_quartic_kernel_fifth_carries48=matrix(carry),
        quartic_independent_source_rows=independent,H_quartic_representation337_by204=matrix(EC),
        remaining_E_indices=ei,complete_fifth_order_new329=matrix(E5),complete_degree5_spanning1329=matrix(full5))
    write_gzip(output/'fifth-complete-spanning.json.gz',packet)
    print('SAVED COMPLETE exact fifth spanning;48 syzygies+329 order5 generators',time.monotonic()-started,flush=True)
    # The modular step selects candidate new generators, never a QQ rank upper.
    prime=32003
    if any(v.q%prime==0 for v in full5):prime=1000003
    assert all(v.q%prime for v in full5)
    inverses={}
    def mod(v):
        denominator=int(v.q)
        if denominator not in inverses:inverses[denominator]=pow(denominator,-1,prime)
        return int(v.p)*inverses[denominator]%prime
    path=output/'modular-rank-input.txt'
    with path.open('w')as out:
        out.write(f'{full5.cols} {full5.rows} {prime}\n')
        for i in range(full5.cols):out.write(' '.join(str(mod(full5[j,i]))for j in range(full5.rows))+'\n')
    run=subprocess.run([str(rank_program),str(path)],capture_output=True,text=True,timeout=25,check=True)
    (output/'modular-rank.log').write_text(run.stdout)
    lines=run.stdout.splitlines();rank=int(lines[0]);pivots=list(map(int,lines[1].split()));assert len(pivots)==rank
    selected=[i for i in pivots if i>=952]
    z=s.symbols('z0:8')
    def polynomials(A,words):
        return [str(s.factor(sum(v*s.prod(z[j]for j in word)for v,word in zip(row,words)if v)))for row in A.tolist()]
    lower_polys=polynomials(A2,m2)+polynomials(C3,m3)+polynomials(D4,m4)
    selected_polys=polynomials(full5[selected,:],m5)
    (output/'selected-fifth.json').write_text(json.dumps(dict(prime=prime,modular_rank_lower=rank,selected_full_spanning_rows=selected,
        selected_QQ_fifth_polynomials=selected_polys,lower21_polynomials=lower_polys),indent=2)+'\n')
    text='R=QQ['+','.join(map(str,z))+'];\nI=ideal('+','.join(p.replace('**','^')for p in lower_polys+selected_polys)+');\n'
    text+='print("DIM5="|toString dim(R/I));print("DEG5="|toString degree(R/I));print("HS5="|toString hilbertSeries(R/I));\n'
    (output/'selected-fifth.m2').write_text(text)
    gb_status='not_run'
    if run_m2:
        executable=shutil.which('M2');assert executable
        try:
            with (output/'selected-fifth.log').open('w')as log:
                completed=subprocess.run([executable,'--script',str(output/'selected-fifth.m2')],stdout=log,stderr=subprocess.STDOUT,timeout=40)
            gb_status='completed'if completed.returncode==0 else f'exit{completed.returncode}'
        except subprocess.TimeoutExpired:gb_status='timeout40s'
    paths=[Path(__file__),ROOT/'code/q_positive_weight6_initial5_rank_20261103_v63.cpp',
        qfolder/'report.json',qfolder/'certificate.json.gz',ffolder/'report.json',ffolder/'certificate.json.gz',
        ROOT/'code/q_positive_weight6_initial4_20261103_v63.py']
    note=ROOT/'research/notes/Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63.md'
    if note.exists():paths.append(note)
    report=dict(id='Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63',all_passed=True,
        exact_QQ_all350_fifth_formal_residuals_and_selected95_zero=True,
        complete_degree5_initial_spanning_QQ_rows=1329,lower_degree_multiple_rows=952,
        complete_order5_generators=329,complete_degree4_syzygy_carries=48,
        all_lower_degree_cancellations_accounted_for=True,
        modular_prime=prime,modular_degree5_initial_rank_lower=rank,selected_new_modular_classes=len(selected),
        QQ_rank_upper_or_Hilbert_function_at5_claimed=False,
        selected_QQ_GB_status=gb_status,full_initial_or_actual_branches_exhaustion_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={p.name:sha(p)for p in output.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('COMPLETE INITIAL5 GF LOWER',rank,'NEW',len(selected),'GB',gb_status,'SEC',report['seconds'],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--rank-program',type=Path,required=True);ap.add_argument('--m2',action='store_true')
    args=ap.parse_args();main(args.output,args.rank_program,args.m2)

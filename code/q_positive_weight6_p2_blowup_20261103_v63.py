#!/usr/bin/env python3
"""Bounded actual P2-direction blowup initial through degree3, all six normals allowed."""
import argparse
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import time
import sympy as s
from q_positive_weight6_quadratic_model_20261103_v63 import ROOT, read, matrix, sha, pivot_columns
from q_positive_weight6_initial4_20261103_v63 import mm, write_gzip


def main(out, run_m2):
    start = time.monotonic()
    out.mkdir(parents=True, exist_ok=False)
    f2 = ROOT/'reports/Q_POSITIVE_WEIGHT6_QUADRATIC_MODEL_20261103_v63'
    f4 = ROOT/'reports/Q_POSITIVE_WEIGHT6_INITIAL4_20261103_v63'
    f5 = ROOT/'reports/Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63'
    q, c, e = [read(p) for p in [f2/'certificate.json.gz', f4/'certificate.json.gz', f5/'fifth-jets.json.gz']]
    data = [(q['full350_by36_corrected_obstructions'], q['normal_quadratic_monomials']),
            (c['full350_corrected_cubic'], c['normal_monomials3']),
            (c['full350_corrected_quartic'], c['normal_monomials4']),
            (e['full350_fifth_residual'], e['normal_monomials5'])]
    normal = {i:j for j,i in enumerate([1,3,4,5,6,7])}
    words = [list(itertools.combinations_with_replacement(range(7), d)) for d in range(4)]
    pos = [{w:i for i,w in enumerate(level)} for level in words]
    jets = [s.zeros(350, len(level)) for level in words]
    for degree, (rows, level) in enumerate(data, 2):
        for col, word in enumerate(level):
            target = tuple(sorted([normal[i] for i in word if i in normal] + [6]*(degree-2)))
            if len(target)>3:
                continue
            j = pos[len(target)][target]
            for row in range(350):
                val = rows[row][col]
                if val != '0':
                    jets[len(target)][row,j] += s.Rational(val)
    assert jets[0] == s.zeros(350,1)
    assert jets[1] == s.zeros(350,7)
    Q, C = jets[2:]
    qr = pivot_columns(Q.T)
    A = Q[qr,:]
    qc = pivot_columns(A)
    coeff = mm(Q[:,qc], A[:,qc].inv())
    assert mm(coeff,A) == Q
    residual = C - mm(coeff,C[qr,:])
    cr = pivot_columns(residual.T)
    selected = residual[cr,:]
    variables = s.symbols('n1 n3 n4 n5 n6 n7 r')
    def polys(M, level):
        return [str(s.factor(sum(v*s.prod(variables[i] for i in w) for v,w in zip(row,level) if v)))
                for row in M.tolist()]
    qp = polys(A,words[2]); cp = polys(selected,words[3])
    # The degree3 initial source consists of S1*Q plus the constant-cancellation cubics.
    rows3=[]
    index3=pos[3]
    for row in A.tolist():
        for variable in range(7):
            result=[s.S.Zero]*len(words[3])
            for value,word in zip(row,words[2]):
                if value:result[index3[tuple(sorted(word+(variable,)))]] += value
            rows3.append(result)
    source3=s.Matrix(rows3).col_join(selected)
    rank3=len(pivot_columns(source3.T))
    payload=dict(direction=[1,0,1,0,0,0,0,0],direction_chart='z2=1',
                 direction_base_parameter='z0-1, set to0 in this closed fiber',
                 normal_coordinates=[1,3,4,5,6,7],radial_coordinate=6,
                 full350_quadratic28=matrix(Q),full350_cubic84=matrix(C),
                 independent_quadratic_rows=qr,quadratic_combination350=matrix(coeff),
                 cubic_residual350=matrix(residual),independent_cubic_rows=cr,
                 quadratic_polynomials=qp,cubic_polynomials=cp,
                 complete_degree3_initial_rank=rank3)
    write_gzip(out/'certificate.json.gz',payload)
    m2='R=QQ[n1,n3,n4,n5,n6,n7,r];\nI=ideal('+','.join(p.replace('**','^')for p in qp+cp)+');\n'
    m2+='print("DIM="|toString dim(R/I));print("DEG="|toString degree(R/I));print("GB="|toString gens gb I);\n'
    m2+='K=ideal(n1,n7,n5*r,r*(15*n3-5*n6-10368*r),n5*(n4+2*n6));\n'
    m2+='assert(dim(R/I)==3);assert(degree(R/I)==6);J=radical I;assert(gens J%gb K==0);assert(gens K%gb J==0);\n'
    m2+='L=ideal(n1,n5,n7,15*n3-5*n6-10368*r);KR=saturate(K,ideal(r));assert(gens KR%gb L==0);assert(gens L%gb KR==0);\n'
    m2+='assert(r%gb J!=0);print("RAD="|toString J);print("RADIAL_SATURATED_INITIAL_RADICAL="|toString L);print("ALL_EXACT_ASSERTIONS_PASS");\n'
    (out/'initial3.m2').write_text(m2)
    status='not_run'
    if run_m2:
        with (out/'initial3.log').open('w') as log:
            job=subprocess.run([shutil.which('M2'),'--script',str(out/'initial3.m2')],stdout=log,stderr=subprocess.STDOUT,timeout=35)
        assert job.returncode==0
        status='passed'
    paths=[Path(__file__),f2/'certificate.json.gz',f4/'certificate.json.gz',f5/'fifth-jets.json.gz',
           ROOT/'code/q_positive_weight6_initial4_20261103_v63.py']
    note=ROOT/'research/notes/Q_POSITIVE_WEIGHT6_P2_BLOWUP_20261103_v63.md'
    if note.exists():paths.append(note)
    report=dict(id='Q_POSITIVE_WEIGHT6_P2_BLOWUP_20261103_v63',all_passed=True,
                exact_actual_closed_blowup_direction=[1,0,1,0,0,0,0,0],
                all6_off_P2_normals_allowed=True,full350_linear_rank=0,
                complete_quadratic_rank=len(qr),complete_cubic_initial_rank=rank3,
                new_cubic_constant_cancellation_rank=len(cr),
                degree3_initial_ideal_dimension=3,degree3_initial_ideal_degree=6,
                degree3_initial_radical='(n1,n7,n5*r,r*(15*n3-5*n6-10368*r),n5*(n4+2*n6))',
                radial_not_in_initial_radical=True,
                radial_saturated_initial_radical='(n1,n5,n7,15*n3-5*n6-10368*r)',
                full_formal_r_saturation_or_actual_branch_claimed=False,
                quadratic_polynomials=qp,cubic_polynomials=cp,tiny_GB_status=status,
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
                evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},
                seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('P2 actual blowup linear0/Qrank',len(qr),'newCrank',len(cr),'initial3rank',rank3,'seconds',report['seconds'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--m2',action='store_true')
    args=p.parse_args();main(args.output,args.m2)

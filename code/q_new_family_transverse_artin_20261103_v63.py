#!/usr/bin/env python3
"""Closed transverse fiber at a specified new-family point; no generic inference."""
import argparse
import json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from q_positive_weight6_quadratic_model_20261103_v63 import ROOT,read,matrix,sha,pivot_columns,physical_linearization,physical_quadratic
from q_positive_weight6_initial4_20261103_v63 import mm,write_gzip


def main(output,point):
    import time
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    family_folder=ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63'
    data=read(family_folder/'certificate.json.gz');ns=s.symbols('n0:3');sub=dict(zip(ns,point))
    q=s.Matrix([s.sympify(v).subs(sub)for v in data['canonical150']])
    # Retain the FULL canonical126 physical directions. The26 cuts are
    # transverse to the ACTUAL family+group image, not all group directions.
    basis_path=ROOT/'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63/certificate.json.gz'
    base=read(basis_path);C=s.Matrix(base['canonical_basis']);free=base['canonical_free_raw_indices']
    assert C*q[free,:]==q
    key=str(tuple(point));M=s.Matrix(data['rank_certificates'][key]['full150_by29_matrix'])
    assert C*M[free,:]==M
    image_columns=pivot_columns(M);assert len(image_columns)==26
    N=M[:,image_columns];cuts=pivot_columns(N[free,:].T);assert len(cuts)==26
    image_minor=DomainMatrix.from_Matrix(N[[free[i]for i in cuts],:]).convert_to(s.QQ).det();assert image_minor
    remaining=[i for i in range(126)if i not in cuts];S=C[:,remaining]
    L=physical_linearization(q);assert L*M==s.zeros(350,29)
    J=L*S;dep=pivot_columns(J);assert len(dep)==98
    rows=pivot_columns(J[:,dep].T);assert len(rows)==98
    B=J[rows,dep];inverse=DomainMatrix.from_Matrix(B).inv().to_Matrix();assert B*inverse==s.eye(98)
    normal=[i for i in range(100)if i not in dep];assert len(normal)==2
    T=s.zeros(100,2)
    for j,k in enumerate(normal):T[k,j]=1
    values=-mm(inverse,J[rows,normal])
    for j,k in enumerate(dep):T[k,:]=values[j,:]
    raw=mm(S,T);assert L*raw==s.zeros(350,2)
    Q,monomials=physical_quadratic(raw);assert monomials==[(0,0),(0,1),(1,1)]
    H2=-mm(inverse,Q[rows,:]);O=Q+mm(J[:,dep],H2)
    assert O[rows,:]==s.zeros(98,3)
    qr=pivot_columns(O.T);rank=len(qr)
    minor=O[qr,:].det()if rank==3 else None
    payload=dict(id='Q_NEW_FAMILY_TRANSVERSE_ARTIN_20261103_v63',point=list(map(str,point)),physical_center150=list(map(str,q)),
        family_image_columns=image_columns,family_image_cut_canonical_rows=cuts,family_image_cut_raw_rows=[free[i]for i in cuts],
        family_image26_minor=str(image_minor),remaining100_canonical_columns=remaining,
        implicit98_rows=rows,implicit98_columns=dep,implicit98_inverse=matrix(inverse),normal2_columns=normal,
        normal2_raw_indices=[free[remaining[i]]for i in normal],tangent150_by2=matrix(raw),
        full350_by100_linearization=matrix(J),full350_by3_uncorrected_quadratic=matrix(Q),
        implicit98_by3_quadratic_correction=matrix(H2),full350_by3_true_quadratic_obstructions=matrix(O),
        selected_quadratic_rows=qr,full_m2_quadratic_minor=str(minor)if minor else None)
    write_gzip(output/'certificate.json.gz',payload)
    paths=[Path(__file__),family_folder/'report.json',family_folder/'certificate.json.gz',
        ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63_REVIEW/report.json',basis_path,
        ROOT/'code/q_positive_weight6_quadratic_model_20261103_v63.py',ROOT/'code/q_positive_weight6_initial4_20261103_v63.py']
    result=dict(id=payload['id'],all_passed=True,point=list(map(str,point)),full_canonical126_physical350_retained=True,
        actual_family_plus_group_rank=26,closed_slice_variables=100,exact_linear_implicit_rank=98,
        closed_transverse_embedding_dimension=2,exact_corrected_quadratic_rank=rank,
        closed_transverse_Artin_length_exact=3 if rank==3 else None,
        closed_transverse_square_zero_maximal_ideal=rank==3,
        actual_completion_finiteness_proved_by_quadrics_alone=rank==3,
        group_action_product_with26_claimed=False,generic_rank98_flatness_or_nonreducedness_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={p.name:sha(p)for p in output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('NEW FAMILY CLOSED TRANSVERSE2 Qrank',rank,'Artinlength',result['closed_transverse_Artin_length_exact'],'SEC',result['seconds'],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--point',nargs=3,default=['1','1','1'])
    args=ap.parse_args();main(args.output,list(map(s.Rational,args.point)))

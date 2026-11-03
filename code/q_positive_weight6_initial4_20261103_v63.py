#!/usr/bin/env python3
"""Exact eight-normal implicit3/4jets and complete initial cancellations through4."""
import argparse
from collections import defaultdict
import gzip
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from q_positive_weight6_quadratic_model_20261103_v63 import ROOT,BASIS,TRIPLES,MON3,read,matrix,sha,pivot_columns


def mm(A,B):
    return (DomainMatrix.from_Matrix(A).convert_to(s.QQ)*DomainMatrix.from_Matrix(B).convert_to(s.QQ)).to_Matrix()


def bilinear(U,mu,V,mv):
    def encode(M,words):
        result={}
        for row,(pair,mon)in enumerate(BASIS):
            poly={word:s.QQ.convert(M[row,j])for j,word in enumerate(words)if M[row,j]}
            if poly:result.setdefault(pair,{})[mon]=poly
        for i,j in list(result):result[j,i]={m:{w:-c for w,c in p.items()}for m,p in result[i,j].items()}
        return result
    left,right=encode(U,mu),encode(V,mv)
    words=list(itertools.combinations_with_replacement(range(8),len(mu[0])+len(mv[0])))
    index={m:i for i,m in enumerate(words)};xm={m:i for i,m in enumerate(MON3)}
    result=[defaultdict(lambda:s.QQ.zero)for _ in range(350)]
    for r,(i,j,k)in enumerate(TRIPLES):
        for b,c,d in [(i,j,k),(j,k,i),(k,i,j)]:
            for ell in range(5):
                for mon,ul in left.get((b,ell),{}).items():
                    for other,vr in right.get((c,d),{}).items():
                        if ell not in other:continue
                        rem=list(other);rem.remove(ell)
                        row=35*r+xm[tuple(sorted(mon+tuple(rem)))];multiplicity=other.count(ell)
                        for u,a in ul.items():
                            for v,z in vr.items():result[row][tuple(sorted(u+v))]+=multiplicity*a*z
    out=s.zeros(350,len(words))
    for i,poly in enumerate(result):
        for mon,coefficient in poly.items():
            if coefficient:out[i,index[mon]]=s.QQ.to_sympy(coefficient)
    return out,words


def multiply_words(A,monomials,degree):
    factors=list(itertools.combinations_with_replacement(range(8),degree))
    words=list(itertools.combinations_with_replacement(range(8),len(monomials[0])+degree))
    index={word:i for i,word in enumerate(words)}
    result=s.zeros(A.rows*len(factors),len(words))
    for i in range(A.rows):
        coefficients=[(m,A[i,j])for j,m in enumerate(monomials)if A[i,j]]
        for k,factor in enumerate(factors):
            for mon,c in coefficients:result[i*len(factors)+k,index[tuple(sorted(mon+factor))]]+=c
    return result,words


def write_gzip(path,payload):
    with path.open('wb')as f:
        with gzip.GzipFile(fileobj=f,mode='wb',filename='',mtime=0)as out:out.write((json.dumps(payload,separators=(',',':'))+'\n').encode())


def validate_physical(q,stages,monomials,residuals):
    """A separate scalar-series convolution of the physical tensor at3 QQ jets.

    The all-coefficient proof is bilinear(); this independent dense evaluation
    checks signs, degree placement and every one of350 recorded residuals.
    """
    tested=[]
    for point in [[1]*8,[1,-1,2,-2,3,-3,4,-4],[2,3,5,7,11,13,17,19]]:
        def evaluate(A,words):
            vector=s.Matrix([s.prod(point[j]for j in word)for word in words])
            return [s.QQ.convert(v)for v in mm(A,vector)]
        values=[list(map(s.QQ.convert,q)),*[evaluate(A,words)for A,words in zip(stages,monomials)]]
        tensors={}
        for row,(pair,mon)in enumerate(BASIS):
            coeff=[value[row]for value in values]
            if any(coeff):tensors.setdefault(pair,{})[mon]=coeff
        for i,j in list(tensors):tensors[j,i]={m:[-v for v in c]for m,c in tensors[i,j].items()}
        positions={m:i for i,m in enumerate(MON3)}
        actual=[[s.QQ.zero]*5 for _ in range(350)]
        for row,(i,j,k)in enumerate(TRIPLES):
            for b,c,d in [(i,j,k),(j,k,i),(k,i,j)]:
                for ell in range(5):
                    for mon,left in tensors.get((b,ell),{}).items():
                        for other,right in tensors.get((c,d),{}).items():
                            if ell not in other:continue
                            rem=list(other);rem.remove(ell)
                            target=35*row+positions[tuple(sorted(mon+tuple(rem)))]
                            for degree in range(5):
                                actual[target][degree]+=other.count(ell)*sum((left[r]*right[degree-r]for r in range(degree+1)),s.QQ.zero)
        expected=[evaluate(A,words)for A,words in zip(residuals,monomials[1:])]
        assert all(actual[i][:2]==[0,0]and actual[i][2:]==[values[i]for values in expected]for i in range(350))
        tested.append(point)
    return tested


def main(output,run_m2=False):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports/Q_POSITIVE_WEIGHT6_QUADRATIC_MODEL_20261103_v63'
    old=read(folder/'report.json');assert old['all_passed']
    for filename,digest in old['input_sha256'].items():assert sha(ROOT/filename)==digest
    data=read(folder/'certificate.json.gz')
    base=read(ROOT/data['canonical_basis_dependency'])
    C=s.Matrix(base['canonical_basis'])
    S=C[:,data['remaining103_canonical_indices']]
    columns=data['implicit_coordinate_columns'];rows=data['implicit_equation_rows']
    D=S[:,columns];B=s.Matrix(data['implicit_inverse95'])
    J=s.Matrix(data['full350_by103_linear_Jacobian']);JD=J[:,columns]
    R1=s.Matrix(data['tangent150_by8']);mu1=[(i,)for i in range(8)]
    H2=s.Matrix(data['implicit95_by36_second_order_coefficients']);mu2=[tuple(m)for m in data['normal_quadratic_monomials']]
    Q2=s.Matrix(data['full350_by36_corrected_obstructions'])
    R2=mm(D,H2)
    C12,mu3=bilinear(R1,mu1,R2,mu2);C21,check=bilinear(R2,mu2,R1,mu1);assert mu3==check
    C3=C12+C21;H3=-mm(B,C3[rows,:]);Q3=C3+mm(JD,H3)
    assert Q3[rows,:]==s.zeros(95,len(mu3))
    R3=mm(D,H3)
    C13,mu4=bilinear(R1,mu1,R3,mu3);C22,check=bilinear(R2,mu2,R2,mu2);assert check==mu4
    C31,check=bilinear(R3,mu3,R1,mu1);assert check==mu4
    C4=C13+C22+C31;H4=-mm(B,C4[rows,:]);Q4=C4+mm(JD,H4)
    assert Q4[rows,:]==s.zeros(95,len(mu4))
    print('FULL350 implicit3/4 QQ jets',time.monotonic()-started,flush=True)
    # An invertible constant350-row change isolates5 generators A with
    # independent quadratic parts, and345 generators G of order>=3.
    chosen=pivot_columns(Q2.T);assert len(chosen)==5
    A2=Q2[chosen,:];A3=Q3[chosen,:];A4=Q4[chosen,:]
    qc=pivot_columns(A2);assert len(qc)==5
    coefficients=mm(Q2[:,qc],DomainMatrix.from_Matrix(A2[:,qc]).inv().to_Matrix())
    assert mm(coefficients,A2)==Q2
    G3=Q3-mm(coefficients,A3);G4=Q4-mm(coefficients,A4)
    remaining=[i for i in range(350)if i not in chosen]
    G3=G3[remaining,:];G4=G4[remaining,:]
    W3,check=multiply_words(A2,mu2,1);assert check==mu3
    all3=W3.col_join(G3)
    piv3=pivot_columns(all3.T);r3=len(piv3)
    cubic_rows=[i-W3.rows for i in piv3 if i>=W3.rows]
    cubic=G3[cubic_rows,:]
    # ALL degree3 cancellations: linear coefficients of the five A plus
    # constant coefficients of all345 G. No uncomputed syzygy is omitted.
    kernel=DomainMatrix.from_Matrix(all3.T).nullspace().to_Matrix()
    assert kernel.rows==all3.rows-r3 and mm(kernel,all3)==s.zeros(kernel.rows,120)
    source4,check=multiply_words(A3,mu3,1);assert check==mu4
    source4=source4.col_join(G4)
    carry=mm(kernel,source4)
    lower4,check=multiply_words(A2,mu2,2);assert check==mu4
    cubic_multiples,check=multiply_words(cubic,mu3,1);assert check==mu4
    lower4=lower4.col_join(cubic_multiples)
    all4=lower4.col_join(carry)
    piv4=pivot_columns(all4.T);r4=len(piv4)
    quartic_rows=[i-lower4.rows for i in piv4 if i>=lower4.rows]
    quartic=carry[quartic_rows,:]
    R4=mm(D,H4)
    checks=validate_physical(s.Matrix(data['physical_center150']),[R1,R2,R3,R4],[mu1,mu2,mu3,mu4],[Q2,Q3,Q4])
    z=s.symbols('z0:8')
    def polys(A,words):return [str(s.factor(sum(c*s.prod(z[j]for j in mon)for c,mon in zip(row,words)if c)))for row in A.tolist()]
    qp=polys(A2,mu2);cp=polys(cubic,mu3);dp=polys(quartic,mu4)
    payload=dict(id='Q_POSITIVE_WEIGHT6_INITIAL4_20261103_v63',quadratic_dependency=str((folder/'certificate.json.gz').relative_to(ROOT)),
        normal_monomials3=mu3,normal_monomials4=mu4,
        raw150_homogeneous_implicit2=matrix(R2),raw150_homogeneous_implicit3=matrix(R3),
        dependent95_homogeneous_implicit3=matrix(H3),dependent95_homogeneous_implicit4=matrix(H4),
        full350_corrected_cubic=matrix(Q3),full350_corrected_quartic=matrix(Q4),
        independent_quadratic_original_rows=chosen,quadratic_combination350_by5=matrix(coefficients),
        remaining345_original_rows=remaining,order3_generators_cubic345=matrix(G3),order3_generators_quartic345=matrix(G4),
        complete_degree3_source_matrix=matrix(all3),complete_degree3_kernel=matrix(kernel),
        complete_degree3_kernel_carry4=matrix(carry),cubic_minimal_row_indices=cubic_rows,quartic_minimal_carry_indices=quartic_rows,
        quadratic_generators=qp,cubic_new_generators=cp,quartic_new_generators=dp,
        homogeneous_degree3_rank=r3,homogeneous_degree4_rank=r4,
        independent_scalar_physical_QQ_jet_checks=checks)
    archive=output/'certificate.json.gz';write_gzip(archive,payload)
    text='R=QQ['+','.join(map(str,z))+'];\nI3=ideal('+','.join(p.replace('**','^')for p in qp+cp)+');\n'
    text+='I4=ideal('+','.join(p.replace('**','^')for p in qp+cp+dp)+');\n'
    text+='print("DIM3="|toString dim(R/I3));print("DEG3="|toString degree(R/I3));print("HS3="|toString hilbertSeries(R/I3));\n'
    text+='print("DIM4="|toString dim(R/I4));print("DEG4="|toString degree(R/I4));print("HS4="|toString hilbertSeries(R/I4));\n'
    (output/'initial4.m2').write_text(text)
    if run_m2:
        executable=shutil.which('M2');assert executable
        with (output/'initial4.log').open('w')as log:
            subprocess.run([executable,'--script',str(output/'initial4.m2')],stdout=log,stderr=subprocess.STDOUT,timeout=40,check=True)
        log=(output/'initial4.log').read_text()
        assert all(marker in log for marker in ['DIM3=5','DEG3=1','DIM4=3','DEG4=6'])
    paths=[Path(__file__),folder/'report.json',folder/'certificate.json.gz',ROOT/data['canonical_basis_dependency'],
           ROOT/'code/q_positive_weight6_quadratic_model_20261103_v63.py']
    note=ROOT/'research/notes/Q_POSITIVE_WEIGHT6_INITIAL4_20261103_v63.md'
    if note.exists():paths.append(note)
    result=dict(id=payload['id'],all_passed=True,all350_original_formal_residuals_retained=True,
        exact_implicit_orders3_and4_over_QQ=True,selected95_residuals_zero_through4=True,
        independent_quadratic_generators=5,new_cubic_generator_count=len(cp),new_quartic_generator_count=len(dp),
        degree3_source_rows=all3.rows,degree3_source_rank=r3,complete_degree3_kernel_dimension=kernel.rows,
        exact_QQ_degree4_initial_rank=r4,full_actual_transverse_Hilbert_function_through4=[1,8,31,120-r3,330-r4],
        complete_true_initial_ideal_through_degree4_proved=True,
        independent_scalar_physical_QQ_jet_checks=checks,
        tiny8variable_QQ_initial_ideal_GB_run=run_m2,
        through3_ideal_dimension=5 if run_m2 else None,through4_ideal_dimension=3 if run_m2 else None,
        through4_ideal_degree=6 if run_m2 else None,
        actual_transverse_local_dimension_upper_bound=3 if run_m2 else None,
        actual_projective_Poisson_local_dimension_upper_bound=25 if run_m2 else None,
        full_initial_ideal_actual_dimension_or_all_branches_exhaustion_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},certificate_sha256=sha(archive),seconds=time.monotonic()-started)
    result['evidence_sha256']={p.name:sha(p)for p in output.iterdir()if p.is_file()and p.name!='report.json'}
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('TRUEINITIAL4',len(cp),'CUBICS',len(dp),'QUARTICS; HF',result['full_actual_transverse_Hilbert_function_through4'],'SEC',result['seconds'],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--m2',action='store_true',help='Run only the tiny8-variable exact initial3/4 ideals')
    args=ap.parse_args();main(args.output,args.m2)

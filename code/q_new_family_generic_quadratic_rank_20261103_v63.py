#!/usr/bin/env python3
"""Generic intrinsic2-normal quadratic obstruction rank1 by full augmented rank99."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import re
import subprocess
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring
from q_positive_weight6_quadratic_model_20261103_v63 import physical_linearization,physical_quadratic

ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_bytes())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
def build(args):
    start=time.monotonic();n0,n1,n2=s.symbols('n0 n1 n2');params=[n0,n2]
    family_path=ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz';family=load(family_path)
    smith_path=ROOT/'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63/certificate.json.gz';smith=load(smith_path)
    normal_path=ROOT/'reports/Q_NEW_FAMILY_POLYNOMIAL_NORMAL_KERNEL_20261103_v63/certificate.json.gz';normal=load(normal_path)
    C=s.Matrix([[s.Rational(v)for v in row]for row in smith['canonical_basis']]);q=s.Matrix([s.sympify(v).subs(n1,1)for v in family['canonical150']])
    R0=s.QQ.poly_ring(*params)
    J=(DomainMatrix.from_Matrix(physical_linearization(q)).convert_to(R0)*DomainMatrix.from_Matrix(C).convert_to(R0)).to_Matrix()
    T=s.Matrix([[s.sympify(v)for v in row]for row in normal['extra_normal126_by2']]);physicalT=C*T
    assert (DomainMatrix.from_Matrix(J).convert_to(R0)*DomainMatrix.from_Matrix(T).convert_to(R0)).to_Matrix()==s.zeros(350,2)
    H,monomials=physical_quadratic(physicalT);assert monomials==[(0,0),(0,1),(1,1)]
    assert H==s.Matrix([[s.sympify(v)for v in row]for row in normal['full350_by3_Hessian_columns']])
    M=J.row_join(H);assert M.shape==(350,129)
    R,*nv=ring('n0,n2',s.QQ);A=[[R.from_expr(s.expand(v))for v in row]for row in M.tolist()]
    witness=load(args.witness);steps=witness.get('steps',witness.get('rational_unit_pivots'));assert len(steps)==79
    rows=list(range(350));cols=list(range(129));verified=[]
    for oi,oj,c in steps:
        i=rows.index(oi);j=cols.index(oj);constant=R.domain.convert(s.Rational(c));assert A[i][j]==constant
        prow=A[i];A=[[v-(row[j]/constant)*prow[k]for k,v in enumerate(row)if k!=j]for l,row in enumerate(A)if l!=i]
        verified.append([oi,oj,str(constant)]);del rows[i];del cols[j]
    assert len(A)==271 and len(A[0])==50
    rr=[i for i,row in enumerate(A)if any(row)];cc=[j for j in range(50)if any(row[j]for row in A)]
    zerorows=[rows[i]for i,row in enumerate(A)if not any(row)];zerocols=[cols[j]for j in range(50)if not any(row[j]for row in A)]
    assert len(zerorows)==98 and len(zerocols)==8
    residual=[[A[i][j]for j in cc]for i in rr];assert (len(residual),len(residual[0]))==(173,42)
    args.work.mkdir(parents=True,exist_ok=True);input_path=args.work/'matrix.txt'
    with input_path.open('w')as stream:
        stream.write('2 173 42\n')
        for row in residual:
            for value in row:
                den=s.ilcm(1,1,*(s.denom(s.Rational(c))for c in value.values()));stream.write(str(len(value))+'\n')
                for exponent,c in value.items():stream.write(str(int(s.Rational(c)*den))+' '+' '.join(map(str,exponent))+'\n')
                stream.write('1\n'+str(den)+' 0 0\n')
    if args.native:
        finished=subprocess.run([str(args.native.resolve()),str(input_path.resolve())],capture_output=True,text=True,check=True,timeout=180);native_log=finished.stdout;reused=False
    else:
        assert args.completed_input and args.completed_log
        assert input_path.read_bytes()==args.completed_input.read_bytes()
        native_log=args.completed_log.read_text();reused=True
    assert re.search(r'EXACT_FUNCTION_FIELD_RANK 20 ROWS 173 COLUMNS 42',native_log)
    sample={n0:1,n2:1};Mp=M.subs(sample);pivcols=list(DomainMatrix.from_Matrix(Mp).rref()[1]);pivrows=list(DomainMatrix.from_Matrix(Mp[:,pivcols].T).rref()[1]);assert len(pivcols)==99
    det99=s.Rational(DomainMatrix.from_Matrix(Mp[pivrows,pivcols]).convert_to(s.QQ).det());assert det99
    # A point nondegenerate quadratic row certifies the generic open condition.
    Jp=J.subs(sample);Jcols=list(DomainMatrix.from_Matrix(Jp).rref()[1]);Jrows=list(DomainMatrix.from_Matrix(Jp[:,Jcols].T).rref()[1]);assert len(Jcols)==98
    Hp=H.subs(sample);correction=Jp[:,Jcols]*Jp[Jrows,Jcols].inv()*Hp[Jrows,:];coker=Hp-correction
    assert DomainMatrix.from_Matrix(coker).rank()==1
    qi=next(i for i in range(350)if any(coker[i,:]));binary=list(coker[qi,:]);disc=s.factor(binary[1]**2-4*binary[0]*binary[2]);assert disc!=0
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'native-input.txt').write_bytes(input_path.read_bytes());(args.output/'native-log.txt').write_text(native_log)
    certificate=dict(normalized_parameters=['n0','n2'],full350_by3_quadratic_Hessian=mat(H),rational_unit_pivots=verified,
       after79_shape=[271,50],all98_identically_zero_original_rows=zerorows,all8_identically_zero_original_columns=zerocols,
       remaining_original_rows=[rows[i]for i in rr],remaining_original_columns=[cols[j]for j in cc],full173_by42_polynomial_residual=[[str(v)for v in row]for row in residual],
       lower99_sample=[1,1],lower99_rows=pivrows,lower99_columns=pivcols,lower99_minor=str(det99),
       nondegenerate_binary_QQ_row_index=qi,nondegenerate_binary_QQ_coefficients=list(map(str,binary)),nondegenerate_binary_discriminant=str(disc))
    archive=args.output/'certificate.json.gz'
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write((json.dumps(certificate,separators=(',',':'))+'\n').encode())
    inputs=[Path(__file__),ROOT/'research/notes/Q_NEW_FAMILY_GENERIC_QUADRATIC_RANK_20261103_v63.md',family_path,smith_path,normal_path,
       normal_path.parent/'report.json',ROOT/'code/q_positive_weight6_quadratic_model_20261103_v63.py',ROOT/'code/rational_function_matrix_rank_20261103_v62.cpp',
       ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/transform_audit_homogeneous.json']
    result=dict(id='Q_NEW_FAMILY_GENERIC_QUADRATIC_RANK_20261103_v63',all_passed=True,field='QQ(n0,n2), dense corrected-shear normalization',
       actual_Jacobian350_by126_and_all_three_physical_Hessian_columns_rebuilt=True,all350_Jacobian_times_normal_vectors_zero=True,
       augmented_matrix_shape=[350,129],rational_unit_pivots=79,identically_zero_rows_removed=98,identically_zero_columns_removed=8,
       exact_residual_field_rank=20,residual_shape=[173,42],all_final_native_zero_rows_checked=True,
       exact_generic_augmented_rank=99,exact_generic_Jacobian_rank_from_universal28_kernel=98,
       intrinsic_generic_two_normal_quadratic_obstruction_rank=1,quadratic_form_generically_nondegenerate=True,
       lower99_QQ_minor=str(det99),sample_binary_discriminant=str(disc),
       descent_mod_actual26='Universal JT=0 differentiated along actual family, group and scalar shows mixed Hessian terms lie in imJ; actual family Hessians also lie in imJ. Thus the two complementary universal kernel sections compute the intrinsic normal quadratic map.',
       generic_Artin_length_or_full_flat_model_claimed_in_this_packet=False,
       completed_native_run_reused=reused,native_run_seconds_if_reused=77.9066,
       no_generic_upper_inferred_from_QQ_samples_or_failed_degree2_syzygies=True,
       input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},evidence_sha256={p.name:sha(p)for p in[archive,args.output/'native-input.txt',args.output/'native-log.txt']},seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('EXACT GENERIC AUGMENTED99 / INTRINSIC QUADRATIC RANK1 PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--work',type=Path,required=True);parser.add_argument('--witness',type=Path,required=True);parser.add_argument('--native',type=Path);parser.add_argument('--completed-input',type=Path);parser.add_argument('--completed-log',type=Path);build(parser.parse_args())

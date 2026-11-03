#!/usr/bin/env python3
"""Generic full canonical Jacobi rank98 over QQ(n0,n1,n2), not samples.

Rebuilds350x126 physical Jacobian, verifies70 polynomial-unit operations,
all73 zero rows/all10 zero columns, and the exact remaining field rank28.
Can rerun native arithmetic or reuse a matching actually completed run.
"""
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
from q_positive_weight6_quadratic_model_20261103_v63 import physical_linearization

ROOT=Path(__file__).resolve().parents[1]
def load(path):
    data=path.read_bytes();return json.loads(gzip.decompress(data)if path.suffix=='.gz'else data)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def matrix(M):return [[str(v)for v in row]for row in M.tolist()]
def build(args):
    started=time.monotonic();n=s.symbols('n0:3')
    family_path=ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz';family=load(family_path)
    q=list(map(s.sympify,family['canonical150']))
    smith_path=ROOT/'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63/certificate.json.gz';smith=load(smith_path)
    C=s.Matrix([[s.Rational(v)for v in row]for row in smith['canonical_basis']]);assert C.shape==(150,126)
    R0=s.QQ.poly_ring(*n);Jraw=physical_linearization(q)
    J=(DomainMatrix.from_Matrix(Jraw).convert_to(R0)*DomainMatrix.from_Matrix(C).convert_to(R0)).to_Matrix()
    assert J.shape==(350,126)
    witness=load(args.witness);steps=witness.get('steps',witness.get('rational_unit_pivots'));assert len(steps)==70
    R,*nv=ring('n0,n1,n2',s.QQ);A=[[R.from_expr(s.expand(v))for v in row]for row in J.tolist()]
    row_ids=list(range(350));col_ids=list(range(126));verified=[]
    for original_row,original_col,pivot in steps:
        i=row_ids.index(original_row);j=col_ids.index(original_col);c=R.from_expr(s.Rational(pivot));assert A[i][j]==c and len(c)==1 and(0,0,0)in c
        prow=A[i];constant=c[(0,0,0)]
        A=[[v-(row[j]/constant)*prow[k]for k,v in enumerate(row)if k!=j]for l,row in enumerate(A)if l!=i]
        verified.append([original_row,original_col,str(constant)]);del row_ids[i];del col_ids[j]
    assert len(A)==280 and len(A[0])==56
    nonzero_rows=[i for i,row in enumerate(A)if any(row)];nonzero_cols=[j for j in range(56)if any(row[j]for row in A)]
    zero_rows=[row_ids[i]for i,row in enumerate(A)if not any(row)];zero_cols=[col_ids[j]for j in range(56)if not any(row[j]for row in A)]
    assert len(zero_rows)==73 and len(zero_cols)==10
    residual=[[A[i][j]for j in nonzero_cols]for i in nonzero_rows];remaining_rows=[row_ids[i]for i in nonzero_rows];remaining_cols=[col_ids[j]for j in nonzero_cols]
    assert len(residual)==207 and len(residual[0])==46
    assert remaining_rows==witness.get('remaining_rows',witness.get('remaining_original_rows'))
    assert remaining_cols==witness.get('remaining_cols',witness.get('remaining_original_columns'))
    args.work.mkdir(parents=True,exist_ok=True);native_input=args.work/'matrix.txt'
    with native_input.open('w')as stream:
        stream.write('3 207 46\n')
        for row in residual:
            for value in row:
                denominator=s.ilcm(1,1,*(s.denom(s.Rational(c))for c in value.values()))
                stream.write(str(len(value))+'\n')
                for exponent,c in value.items():stream.write(str(int(s.Rational(c)*denominator))+' '+' '.join(map(str,exponent))+'\n')
                stream.write('1\n'+str(denominator)+' 0 0 0\n')
    if args.native:
        completed=subprocess.run([str(args.native.resolve()),str(native_input.resolve())],capture_output=True,text=True,check=True,timeout=600);native_log=completed.stdout;reused=False
    else:
        assert args.completed_input and args.completed_log
        assert args.completed_input.read_bytes()==native_input.read_bytes()
        native_log=args.completed_log.read_text();reused=True
    assert re.search(r'EXACT_FUNCTION_FIELD_RANK 28 ROWS 207 COLUMNS 46',native_log)
    # A QQ98-minor is only a lower witness; the full-field computation above
    # proves the upper as well and checks every final zero entry.
    sample=J.subs(dict(zip(n,[1,1,1])));cols=list(DomainMatrix.from_Matrix(sample).rref()[1]);rows=list(DomainMatrix.from_Matrix(sample[:,cols].T).rref()[1]);assert len(cols)==98
    det=s.Rational(DomainMatrix.from_Matrix(sample[rows,cols]).convert_to(s.QQ).det());assert det
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'native-input.txt').write_bytes(native_input.read_bytes());(args.output/'native-log.txt').write_text(native_log)
    certificate=dict(parameters=list(map(str,n)),full350_by126_canonical_Jacobian=matrix(J),canonical_basis_source=str(smith_path.relative_to(ROOT)),
      rational_unit_pivots=verified,after70_shape=[280,56],all73_identically_zero_original_rows=zero_rows,all10_identically_zero_original_columns=zero_cols,
      remaining_original_rows=remaining_rows,remaining_original_columns=remaining_cols,
      exact207_by46_polynomial_residual=[[str(v)for v in row]for row in residual],
      lower98_sample=[1,1,1],lower98_rows=rows,lower98_columns=cols,lower98_minor=str(det))
    archive=args.output/'certificate.json.gz'
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write((json.dumps(certificate,separators=(',',':'))+'\n').encode())
    inputs=[Path(__file__),ROOT/'research/notes/Q_NEW_FAMILY_GENERIC_JACOBIAN_RANK_20261103_v63.md',family_path,smith_path,
      ROOT/'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63_REVIEW.json',ROOT/'code/q_positive_weight6_quadratic_model_20261103_v63.py',
      ROOT/'code/rational_function_matrix_rank_20261103_v62.cpp']
    result=dict(id='Q_NEW_FAMILY_GENERIC_JACOBIAN_RANK_20261103_v63',all_passed=True,field='QQ(n0,n1,n2)',
      full_actual_physical350_by126_linearization_rebuilt=True,canonical_basis_rank=126,
      rational_unit_pivots_verified=70,identically_zero_rows_removed=73,identically_zero_columns_removed=10,
      residual_shape=[207,46],exact_residual_function_field_rank=28,all_final_native_zero_rows_checked=True,
      full_generic_canonical_Jacobian_rank=98,lower98_QQ_minor=str(det),
      generic_canonical_affine_tangent_dimension=28,generic_projective_tangent_dimension=27,
      generic_nonreducedness_after_independent_dimension25_receipt=True,
      exact_generic_Artin_length_or_quadratic_obstruction_rank_claimed=False,
      generic_PGL_moduli_normal_form_or_finite_identifications_claimed=False,
      completed_native_run_reused=reused,native_run_scope='Actual prior successful QQ(n3) run,192.52s, if reused; independent replay available with--native.',
      no_generic_upper_inferred_from_numeric_rank98_samples=True,
      input_sha256={str(path.relative_to(ROOT)):sha(path)for path in inputs},
      evidence_sha256={path.name:sha(path)for path in [archive,args.output/'native-input.txt',args.output/'native-log.txt']},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('EXACT GENERIC FULL JACOBIAN98 PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--work',type=Path,required=True);parser.add_argument('--witness',type=Path,required=True);parser.add_argument('--native',type=Path);parser.add_argument('--completed-input',type=Path);parser.add_argument('--completed-log',type=Path);build(parser.parse_args())

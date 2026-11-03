#!/usr/bin/env python3
"""Exact generic firstjet rank of the actual new three-parameter family.

Raw value+linear-jet rank25 is certified over QQ(n0,n1,n2,p0,...,p3),
with20 polynomial unit pivots and a frozen exact native field-rank engine.
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

ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def strings(M):return [[str(v)for v in row]for row in M.tolist()]

def build(args):
    start=time.monotonic();source=ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz'
    saved=json.loads(gzip.decompress(source.read_bytes()));n=s.symbols('n0:3');p=s.symbols('p0:4');x=s.Matrix(s.symbols('x y z w'));t=s.Matrix(s.symbols('t0:4'))
    variables=[*n,*p];pairs=list(itertools.combinations(range(4),2));P=s.Matrix([[s.sympify(v)for v in row]for row in saved['full_affine_P']])
    C=P.subs(dict(zip(x,p)),simultaneous=True)
    L=sum((P.diff(x[k]).subs(dict(zip(x,p)),simultaneous=True)*t[k]for k in range(4)),s.zeros(4))
    def jet(A,B):return s.Matrix([A[i,j]for i,j in pairs]+[s.diff(B[i,j],v)for i,j in pairs for v in t])
    columns=[];labels=[]
    for v in variables:columns.append(jet(C.diff(v),L.diff(v)));labels.append(str(v))
    for i,j in itertools.product(range(4),repeat=2):
        A=s.zeros(4);A[i,j]=1
        columns.append(jet(A*C+C*A.T,A*L+L*A.T-sum(((A*t)[k]*L.diff(t[k])for k in range(4)),s.zeros(4))));labels.append(f'GL4_{i}_{j}')
    for i in range(4):
        al=s.eye(4)[:,i];v=C*al
        columns.append(jet(s.zeros(4),-2*(al.dot(t))*C+t*v.T-v*t.T));labels.append(f'perspective_{i}')
    columns.append(jet(C,L));labels.append('P_scalar');M=s.Matrix.hstack(*columns);assert M.shape==(30,28)
    R,*rv=ring('n0,n1,n2,p0,p1,p2,p3',s.QQ)
    A=[[R.from_expr(s.expand(v))for v in row]for row in M.tolist()];units=[];row_ids=list(range(30));column_ids=list(range(28))
    while A and A[0]:
        candidates=[]
        for i,row in enumerate(A):
            for j,f in enumerate(row):
                if len(f)==1 and (0,)*7 in f:
                    cost=sum(len(v)for v in row)*sum(len(rr[j])for rr in A)
                    candidates.append((cost,i,j,f[(0,)*7]))
        if not candidates:break
        cost,i,j,c=min(candidates);pivot=A[i]
        units.append(dict(current_row=i,current_column=j,original_row=row_ids[i],original_column=column_ids[j],pivot=str(c)))
        A=[[v-(row[j]/c)*pivot[k]for k,v in enumerate(row)if k!=j]for l,row in enumerate(A)if l!=i]
        del row_ids[i];del column_ids[j]
    assert len(units)==20 and (len(A),len(A[0]))==(10,8)
    # Every step used a rational UNIT. All identities hold in the polynomial
    # ring itself, not merely at a generic parameter specialization.
    args.work.mkdir(parents=True,exist_ok=True);input_path=args.work/'matrix.txt'
    with input_path.open('w')as stream:
        stream.write('7 10 8\n')
        for row in A:
            for value in row:
                denominator=s.ilcm(1,1,*(s.denom(s.Rational(c))for c in value.values()))
                stream.write(str(len(value))+'\n')
                for exponent,c in value.items():stream.write(str(int(s.Rational(c)*denominator))+' '+' '.join(map(str,exponent))+'\n')
                stream.write('1\n'+str(denominator)+' 0 0 0 0 0 0 0\n')
    completed=subprocess.run([str(args.native.resolve()),str(input_path.resolve())],check=True,capture_output=True,text=True,timeout=60)
    assert re.search(r'EXACT_FUNCTION_FIELD_RANK 5 ROWS 10 COLUMNS 8',completed.stdout)
    total_rank=20+5;assert total_rank==25
    sample=dict(zip(variables,[1,1,1,0,0,0,0]));sampleM=M.subs(sample)
    cols=list(sampleM.rref()[1]);rows=list(sampleM[:,cols].T.rref()[1]);assert len(cols)==25
    determinant=s.Rational(DomainMatrix.from_Matrix(sampleM[rows,cols]).convert_to(s.QQ).det());assert determinant
    Cs=C.subs(sample);pf=s.factor(Cs[0,1]*Cs[2,3]-Cs[0,2]*Cs[1,3]+Cs[0,3]*Cs[1,2]);assert pf==-s.Rational(16,81)
    perspective=sampleM[:,23:27];assert perspective[:6,:]==s.zeros(6,4)and perspective.rank()==4
    assert sampleM[:6,7:23].rank()==6
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'native-input.txt').write_bytes(input_path.read_bytes())
    (args.output/'native-log.txt').write_text(completed.stdout)
    certificate=dict(field_variables=list(map(str,variables)),raw_jet_row_order='sixC entries, then24linearL entries (pair first,variable second)',
      source_column_labels=labels,full30_by28_polynomial_matrix=strings(M),
      twenty_rational_unit_pivots=units,remaining_original_rows=row_ids,remaining_original_columns=column_ids,
      residual10_by8_polynomial_matrix=[[str(v)for v in row]for row in A],
      sample={str(v):str(value)for v,value in sample.items()},sample_value_Pfaffian=str(pf),lower25_minor_rows=rows,lower25_minor_columns=cols,lower25_minor_determinant=str(determinant))
    archive=args.output/'certificate.json.gz'
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write((json.dumps(certificate,separators=(',',':'))+'\n').encode())
    native_source=ROOT/'code/rational_function_matrix_rank_20261103_v62.cpp'
    inputs=[Path(__file__),ROOT/'research/notes/Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63.md',source,source.parent/'report.json',native_source,
      ROOT/'research/notes/KNOWN51_FIRSTJET_IMAGE_DIMENSIONS_20261103_v62.md',ROOT/'reports/KNOWN51_FIRSTJET_IMAGE_DIMENSIONS_20261103_v62_REVIEW.json']
    result=dict(id='Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63',all_passed=True,field='QQ(n0,n1,n2,p0,p1,p2,p3)',
      full_actual_three_parameter_family_and_all_point_directions_retained=True,full_GL4_16_and_perspective4_and_P_scalar_retained=True,
      raw_value_and_linear_jet_generic_image_dimension=25,polynomial_unit_pivots=20,residual_native_function_field_rank=5,
      all_final_native_zero_rows_checked_exactly=True,lower25_QQ_minor=str(determinant),sample_rank4_value_Pfaffian=str(pf),
      independent_value_nuisance_dimension=6,independent_perspective_nuisance_dimension=4,
      canonical_primitive_firstjet_image_dimension=15,does_not_dominate_W16=True,
      image_codimension_in_W16=1,dominant_general_firstjet_counterexample_claimed=False,
      generic_S61_fiber_dimension_if_saturated_image_is_component25=1,
      no_upper_inferred_from_finite_samples=True,
      input_sha256={str(path.relative_to(ROOT)):sha(path)for path in inputs},
      evidence_sha256={path.name:sha(path)for path in [archive,args.output/'native-input.txt',args.output/'native-log.txt']},seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('EXACT GENERIC RAWJET25 / PRIMITIVE BETA15 PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--work',type=Path,required=True);parser.add_argument('--native',type=Path,required=True);build(parser.parse_args())

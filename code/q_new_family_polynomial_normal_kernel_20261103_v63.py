#!/usr/bin/env python3
"""Two universal linear polynomial normals, giving a short generic rank98 proof."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from q_positive_weight6_new_family_20261103_v63 import vector,lie
from q_positive_weight6_quadratic_model_20261103_v63 import physical_linearization,physical_quadratic

ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_bytes())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
def build(args):
    start=time.monotonic();n0,n1,n2=s.symbols('n0 n1 n2');variables=[n0,n2]
    family_path=ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz';family=load(family_path)
    smith_path=ROOT/'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63/certificate.json.gz';smith=load(smith_path)
    C=s.Matrix([[s.Rational(v)for v in row]for row in smith['canonical_basis']]);free=smith['canonical_free_raw_indices'];assert C[free,:]==s.eye(126)
    q=s.Matrix([s.sympify(v).subs(n1,1)for v in family['canonical150']]);P=s.Matrix([[s.sympify(v).subs(n1,1)for v in row]for row in family['full_canonical_homogeneous_P']])
    R=s.QQ.poly_ring(*variables)
    J=(DomainMatrix.from_Matrix(physical_linearization(q)).convert_to(R)*DomainMatrix.from_Matrix(C).convert_to(R)).to_Matrix()
    witness=load(args.witness);T=s.Matrix([[s.sympify(v)for v in row]for row in witness.get('T126x2',witness.get('extra_normal126_by2'))]);assert T.shape==(126,2)
    assert all(s.Poly(v,*variables).total_degree()<=1 for v in T)
    actual=s.Matrix.hstack(*[vector(lie(P,i,j))for i,j in itertools.product(range(5),repeat=2)],q,q.diff(n0),q.diff(n2))[free,:]
    sample={n0:1,n2:1};Ap=actual.subs(sample);chosen=list(DomainMatrix.from_Matrix(Ap).rref()[1]);assert len(chosen)==26
    K=actual[:,chosen].row_join(T)
    product=(DomainMatrix.from_Matrix(J).convert_to(R)*DomainMatrix.from_Matrix(K).convert_to(R)).to_Matrix();assert product==s.zeros(350,28)
    Kp=K.subs(sample);rows=list(DomainMatrix.from_Matrix(Kp.T).rref()[1]);assert len(rows)==28
    det28=s.Rational(DomainMatrix.from_Matrix(Kp[rows,:]).convert_to(s.QQ).det());assert det28
    Jp=J.subs(sample);cols98=list(DomainMatrix.from_Matrix(Jp).rref()[1]);rows98=list(DomainMatrix.from_Matrix(Jp[:,cols98].T).rref()[1]);assert len(cols98)==98
    det98=s.Rational(DomainMatrix.from_Matrix(Jp[rows98,cols98]).convert_to(s.QQ).det());assert det98
    physicalT=C*T;Hess,monomials=physical_quadratic(physicalT);assert monomials==[(0,0),(0,1),(1,1)]
    args.output.mkdir(parents=True,exist_ok=False)
    certificate=dict(normalized_parameters=['n0','n2'],normalization='n1=1 on dense corrected-shear chart D(n0*n1)',
       extra_normal126_by2=mat(T),extra_normal_physical150_by2=mat(physicalT),actual_family_selected26_columns=chosen,
       full126_by28_polynomial_kernel=mat(K),all350_by28_Jacobian_products_zero=True,
       sample=[1,1],kernel28_rows=rows,kernel28_minor=str(det28),Jacobian98_rows=rows98,Jacobian98_columns=cols98,Jacobian98_minor=str(det98),
       full350_by3_Hessian_columns=mat(Hess),quadratic_monomials=monomials)
    archive=args.output/'certificate.json.gz'
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write((json.dumps(certificate,separators=(',',':'))+'\n').encode())
    inputs=[Path(__file__),ROOT/'research/notes/Q_NEW_FAMILY_POLYNOMIAL_NORMAL_KERNEL_20261103_v63.md',family_path,smith_path,
       ROOT/'code/q_positive_weight6_quadratic_model_20261103_v63.py',ROOT/'code/q_positive_weight6_new_family_20261103_v63.py',
       ROOT/'code/r4_positive_weight6_transform_audit_20261103_v63.py',ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/transform_audit_homogeneous.json']
    result=dict(id='Q_NEW_FAMILY_POLYNOMIAL_NORMAL_KERNEL_20261103_v63',all_passed=True,field='QQ(n0,n2), equivalent dense full component by corrected shear normalization',
       two_extra_universal_degree1_polynomial_normals=True,all350_by28_full_physical_Jacobian_kernel_identities_exact=True,
       genuine_family_and_GL_and_P_scalar_kernel_rank=26,full_polynomial_kernel_rank=28,kernel28_QQ_minor=str(det28),
       generic_Jacobian_upper_from_kernel=98,point98_lower_minor=str(det98),exact_generic_canonical_Jacobian_rank=98,
       native_rational_function_rank_not_needed_for_this_proof=True,
       full350_by3_quadratic_Hessian_saved=True,generic_quadratic_cokernel_rank_or_exact_Artin_length_claimed=False,
       input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},certificate_sha256=sha(archive),seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('UNIVERSAL POLYNOMIAL28 KERNEL / GENERIC JAC98 SHORT PROOF PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--witness',type=Path,required=True);build(parser.parse_args())

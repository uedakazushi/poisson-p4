#!/usr/bin/env python3
"""Exact11-dimensional intrinsic beta orbit excludes five small image closures."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import J

ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def strings(M):return [[str(v)for v in row]for row in M.tolist()]
def build(out):
    start=time.monotonic();eta=s.Symbol('eta');certpath=ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63/certificate.json.gz'
    point=json.loads(gzip.decompress(certpath.read_bytes()));u=list(map(s.sympify,point['S61']))
    z,us,beta=primitive_beta(S61());B=beta.subs({v:u[int(str(v)[1:])-1]for v in us})
    pairs=list(itertools.combinations(range(4),2));vec=lambda T:s.Matrix([s.diff(T[i,j],v)for i,j in pairs for v in z])
    variables=s.symbols('m0:16');M=s.Matrix(4,4,variables);linear=s.linear_eq_to_matrix(list(M*J+J*M.T),variables)[0]
    lie=[s.Matrix(4,4,v)for v in linear.nullspace()];assert len(lie)==10
    columns=[]
    for A in lie:
        assert A*J+J*A.T==s.zeros(4)
        action=-A.T*B-B*A-sum(((A*s.Matrix(z))[i]*B.diff(z[i])for i in range(4)),s.zeros(4))
        columns.append(vec(action))
    O=s.Matrix.hstack(*columns,vec(B));rows=list(O.subs(eta,1).T.rref()[1]);assert len(rows)==11
    R=s.QQ.poly_ring(eta);det=s.factor(R.to_sympy(DomainMatrix.from_Matrix(O[rows,:]).convert_to(R).det()))
    polynomial=s.Poly(det,eta);assert len(polynomial.terms())==1 and polynomial.terms()[0][1]!=0
    assert polynomial.terms()[0][0][0]>0
    tablepath=ROOT/'reports/KNOWN51_FIRSTJET_DIMENSION_TABLE_20261103_v62.json';table=json.loads(tablepath.read_text());assert table['all_passed']
    excluded=[(r['component_id'],r['firstjet_image_dimension'])for r in table['records']if r['firstjet_image_dimension']<11]
    assert excluded==[('T21',10),('R4_CA',9),('R4_CB',9),('R4_ZERO_LOG_5H',10),('R4_CMINUS_MAG',10)]
    reviewpath=ROOT/'reports/KNOWN51_FIRSTJET_IMAGE_DIMENSIONS_20261103_v62_REVIEW.json';review=json.loads(reviewpath.read_text());assert review['all_passed']
    out.mkdir(parents=True,exist_ok=False)
    certificate=dict(primitive_beta=str(B),coordinate_variables=list(map(str,z)),sp4_basis=[strings(A)for A in lie],
       full24_by11_beta_orbit_and_scalar_matrix=strings(O),minor_rows=rows,all11_columns=list(range(11)),exact_minor=str(det),excluded_components_and_closed_image_dimensions=excluded)
    (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'research/notes/Q_POSITIVE_WEIGHT6_BETA_ORBIT_20261103_v63.md',certpath,tablepath,reviewpath,
       ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'code/s61_model.py']
    result=dict(id='Q_POSITIVE_WEIGHT6_BETA_ORBIT_20261103_v63',all_passed=True,field='QQ[eta]',
       complete_sp4_lie_algebra_and_full24_coefficient_action_rebuilt=True,
       exact_minor=str(det),orbit_dimension_Sp4_times_Cstar=11,parameter_scope='Every eta!=0',
       all_nonzero_axis_points_outside_entire_known_component_closures=[name for name,dim in excluded],
       proof='Every known closed firstjet image is Sp4 and weighted-Cstar invariant. Containment of beta would contain its11D orbit and contradict the certified image dimension<11.',
       other_component_memberships_or_global_new_component_claimed=False,
       input_sha256={str(path.relative_to(ROOT)):sha(path)for path in inputs},certificate_sha256=sha(out/'certificate.json'),seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('POSITIVE BETA ORBIT11 AND T21/FOUR IMAGE EXCLUSIONS PASS',det,result['seconds'],flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);build(parser.parse_args().output)

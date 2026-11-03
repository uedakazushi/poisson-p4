#!/usr/bin/env python3
"""Independent intrinsic two-form Lie action for the positive beta orbit."""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic();folder=ROOT/'reports/Q_POSITIVE_WEIGHT6_BETA_ORBIT_20261103_v63';r=json.loads((folder/'report.json').read_text());cert=json.loads((folder/'certificate.json').read_text());assert r['all_passed']
 for name,want in r['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==want
 universal=ROOT/'reports/PRIMITIVE_JET_ELLIPTIC_CURVE_20261101.json';j=json.loads(universal.read_text());assert j['all_passed']
 source=ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63/certificate.json.gz';ray=json.loads(gzip.decompress(source.read_bytes()));eta=s.Symbol('eta')
 sub={s.Symbol('u'+str(i+1)):s.sympify(v)for i,v in enumerate(ray['S61'])}
 B=s.Matrix(j['primitive_beta']).subs(sub).applyfunc(s.expand);assert B==s.sympify(cert['primitive_beta'])
 z=list(map(s.Symbol,cert['coordinate_variables']));pairs=list(itertools.combinations(range(4),2));J=s.zeros(4);J[0,1]=J[2,3]=1;J=J-J.T
 lie=list(map(s.Matrix,cert['sp4_basis']));assert len(lie)==10 and s.Matrix.hstack(*[A.reshape(16,1)for A in lie]).rank()==10
 columns=[]
 for A in lie:
  assert A*J+J*A.T==s.zeros(4)
  transformed=(-A.T*B-B*A-sum(((A*s.Matrix(z))[k]*B.diff(z[k])for k in range(4)),s.zeros(4))).applyfunc(s.expand)
  columns.append(s.Matrix([s.diff(transformed[i,j],v)for i,j in pairs for v in z]))
 columns.append(s.Matrix([s.diff(B[i,j],v)for i,j in pairs for v in z]));M=s.Matrix.hstack(*columns)
 assert M==s.Matrix(cert['full24_by11_beta_orbit_and_scalar_matrix'])
 det=s.factor(M[cert['minor_rows'],:].det());assert det==s.sympify(cert['exact_minor'])and det!=0 and s.Poly(det,eta).terms()[0][0][0]>0
 tablepath=ROOT/'reports/KNOWN51_FIRSTJET_DIMENSION_TABLE_20261103_v62.json';table=json.loads(tablepath.read_text());small=[(c['component_id'],c['firstjet_image_dimension'])for c in table['records']if c['firstjet_image_dimension']<11];assert [list(v)for v in small]==cert['excluded_components_and_closed_image_dimensions']
 paths=[Path(__file__),folder/'report.json',folder/'certificate.json',universal,source,tablepath,ROOT/'research/notes/Q_POSITIVE_WEIGHT6_BETA_ORBIT_20261103_v63_REVIEW.md']
 result=dict(id='Q_POSITIVE_WEIGHT6_BETA_ORBIT_20261103_v63_REVIEW',all_passed=True,universal_primitive_two_form_specialization_reused=True,own_full_covariant_Sp4_Lie_action_and_scalar=True,exact_minor=str(det),nonzero_eta_orbit_dimension=11,whole_component_T21_exclusion_by_closed_invariant_image_dimension=True,other_small_image_exclusions=[c for c,d in small],new_component_or_exhaustion_inferred=False,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('INDEPENDENT POSITIVE BETA ORBIT11/T21 REVIEW PASS',result['seconds'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

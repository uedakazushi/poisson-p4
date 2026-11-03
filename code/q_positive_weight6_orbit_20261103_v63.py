#!/usr/bin/env python3
"""Exact projective orbit of the positive weight-six ray and four exclusions."""
import argparse,gzip,hashlib,json,time
from fractions import Fraction as Q
from itertools import combinations,combinations_with_replacement
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 folder=ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63'
 cert=json.loads(gzip.decompress((folder/'certificate.json.gz').read_bytes()))
 eta=s.Symbol('eta');raw=list(map(s.sympify,cert['raw150']))
 pairs=list(combinations(range(5),2));mons=list(combinations_with_replacement(range(5),2))
 basis=[(i,j,m)for i,j in pairs for m in mons]
 expected={64:-s.Rational(256,3)*eta,83:-64*eta,90:s.Rational(8,9),105:-s.Rational(2,9),116:64*eta,121:s.Rational(4,3),137:s.Rational(2,3),140:s.Integer(1)}
 assert {i:v for i,v in enumerate(raw)if v}==expected
 v=[Q(str(x.subs(eta,1)))for x in raw]
 P={(i,j,m):c for (i,j,m),c in zip(basis,v)if c}
 def get(i,j,m):return P.get((i,j,m),0)if i<j else -P.get((j,i,m),0)
 cols=[]
 for a in range(5):
  for b in range(5):
   d={}
   def add(k,c):
    if c:d[k]=d.get(k,Q(0))+c
   for i,j in pairs:
    for m in mons:
     if i==a:add((i,j,m),get(b,j,m))
     if j==a:add((i,j,m),get(i,b,m))
     c=get(i,j,m)
     for pos,k in enumerate(m):
      if k==a:add((i,j,tuple(sorted((m[1-pos],b)))),-c)
   cols.append([d.get(k,0)for k in basis])
 M=s.Matrix(cols).T;N=M.row_join(s.Matrix(v))
 rows=list(N.T.rref()[1]);columns=list(N[rows,:].rref()[1])
 assert len(rows)==len(columns)==23
 minor=N.extract(rows,columns).det();assert minor
 kernel=N.nullspace();assert len(kernel)==3
 K=s.Matrix.hstack(*kernel);assert N*K==s.zeros(150,3)and K.rank()==3
 assert M.rank()==22
 # diag(1,t,t²,t³,t⁴), followed by scalar t^-5, sends eta to t^-6 eta.
 weights=list(range(5));exponents={}
 for k,c in enumerate(raw):
  if not c:continue
  i,j,m=basis[k];power=weights[i]+weights[j]-sum(weights[a]for a in m)-5
  assert power==(-6 if eta in c.free_symbols else 0)
  exponents[k]=power
 catalogpath=ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json'
 catalog=json.loads(catalogpath.read_text());ids=['R4_DOUBLE_SMOOTH_QUADRIC','R4_TRIPLE_HYPERPLANE_QUADRIC']
 for name in ids:
  record=next(r for r in catalog['records']if r['component_id']==name)
  assert record['expected_s61_dimension']+9==22
  assert record['full_s61_component_dominance']
 membership=ROOT/'reports/Q_EXACT_KNOWN_MPS_MEMBERSHIP_20261103_v63.json'
 md=json.loads(membership.read_text());assert md['all_passed']
 for name,digest in md['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
 evidence=dict(nonzero_minor_rows_zero_based=rows,nonzero_minor_columns_zero_based=columns,nonzero_minor=str(minor),GL25_plus_scalar_matrix=[[str(x)for x in row]for row in N.tolist()],kernel_basis=[[str(x)for x in row]for row in K.tolist()],weight_action_exponents=exponents,dimension_exclusions=ids)
 (out/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
 paths=[Path(__file__),folder/'report.json',folder/'certificate.json.gz',catalogpath,membership,ROOT/'legacy/v03/p4_poisson_completeness_v3.tex',ROOT/'legacy/v05/triple_hyperplane.tex',ROOT/'research/notes/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63.md']
 result=dict(id='Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63',all_passed=True,
  exact_GL25_rank=22,exact_GL25_plus_scalar_rank=23,projective_PGL5_orbit_dimension=22,
  every_nonzero_eta_is_in_one_projective_orbit=True,
  explicit_action='diag(1,t,t^2,t^3,t^4), scalar t^-5: eta maps to t^-6 eta',
  orbit_closure_contains_physical_Qa_eight_over_nine=True,
  all_nonzero_ray_points_outside_T39_T40=True,
  all_nonzero_ray_points_outside_two_dimension22_quadric_components=True,
  quadric_exclusion_uses_equal_dimension_orbit_closure_and_closed_pure_fifth_power_locus=True,
  new_component_or_exhaustion_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
 print('POSITIVE RAY PROJECTIVE ORBIT22 AND FOUR EXCLUSIONS PASS',result['seconds'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

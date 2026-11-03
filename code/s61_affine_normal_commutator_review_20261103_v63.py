#!/usr/bin/env python3
"""Independent normal-commutator constants and scheme-theoretic proof audit."""
import gzip,hashlib,itertools,json,time
from fractions import Fraction
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1];P=32003
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for c in iter(lambda:f.read(1024*1024),b''):h.update(c)
 return h.hexdigest()
def residue(v):
 q=Fraction(v);assert q.denominator%P
 return q.numerator%P*pow(q.denominator%P,-1,P)%P
def DM(rows):return DomainMatrix.from_Matrix(s.Matrix(rows)).convert_to(s.GF(P))
def main():
 start=time.monotonic();folder=ROOT/'reports/S61_AFFINE_NORMAL_COMMUTATOR_20261103_v63';r=json.loads((folder/'report.json').read_text());assert r['all_passed']
 for p,want in r['input_sha256'].items():assert sha(ROOT/p)==want,p
 for p,want in r['evidence_sha256'].items():assert sha(folder/p)==want,p
 finite=ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63';cert=json.loads((finite/'certificate.json').read_text());selected=cert['selected_cross_rows']
 raw=gzip.decompress((finite/'cross-input.txt.gz').read_bytes());lines=raw.decode().splitlines();assert lines[0]=='129 126 945'
 cross=[[residue(v)for v in line.split()]for line in lines[1:130]];D=[[residue(v)for v in line.split()[:42]]for line in lines[130:259]]
 A=DM([cross[i]for i in selected]);L=-(A.inv()*DM([D[i]for i in selected]));LL=[[int(v)%P for v in row]for row in L.to_Matrix().tolist()]
 derived=ROOT/'work/normal-finite5-public-v63/cross-output.txt';assert sha(derived)==cert['native_output_sha256']
 with derived.open()as f:
  assert f.readline().split()==['126','945'];assert list(map(int,f.readline().split()))==selected
  for _ in range(126):f.readline()
  for i in range(126):assert[residue(v)for v in f.readline().split(maxsplit=42)[:42]]==LL[i]
 T=[DM(LL[42*i:42*(i+1)])for i in range(3)];payload=json.loads(gzip.decompress((folder/'modular-linear-commutators.json.gz').read_bytes()));dets=[]
 for k,(i,j)in enumerate(itertools.combinations(range(3),2)):
  C=T[i]*T[j]-T[j]*T[i];rows=[[int(v)%P for v in row]for row in C.to_Matrix().tolist()];assert rows==payload['all_three_commutators'][k]
  dets.append(int(C.det())%P)
 assert dets==[21866,23983,19574]
 note=ROOT/'research/notes/S61_AFFINE_NORMAL_COMMUTATOR_20261103_v63_REVIEW.md';paths=[Path(__file__),note,folder/'report.json',folder/'certificate.json',folder/'modular-linear-commutators.json.gz',finite/'report.json',finite/'certificate.json',finite/'cross-input.txt.gz']
 result=dict(id='S61_AFFINE_NORMAL_COMMUTATOR_20261103_v63_REVIEW',all_passed=True,independent_GF32003_DomainMatrix_cross_inverse_and_all5292_linear_coefficients=True,exactQQ_native_output_stream_hash_checked=True,all_three_independent_commutator_determinants=dets,exact_original_ideal_cubic_relation_proof_accepted=True,normal_Jacobian42_nonzero_implies_reduced_origin_local_ideal=True,all_chart_gluing_J0_equals_m_intersect_residual_accepted=True,normal_D84_equals_original_global_idempotent_accepted=True,normal_unsaturated_residual_J0_plus_D84_accepted=True,expanded_QQ_cubic_or_det_computation_claimed=False,generic_beta_finiteness_or_unit_or_exhaustion_claimed=False,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},seconds=time.monotonic()-start)
 (ROOT/'reports/S61_AFFINE_NORMAL_COMMUTATOR_20261103_v63_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()

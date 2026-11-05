#!/usr/bin/env python3
"""Independent fixed-certificate verification; no Gaussian/rank search."""
import hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_T2_FULL_G3_LINEAR_KERNEL_20261105_v63_REVIEW'
PRIMARY=STEM.removesuffix('_REVIEW')
GATE=PRIMARY.replace('_KERNEL_','_KERNEL_INPUT_GATE_')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def sp(a):
 d={int(i):Q(c) for i,c in a}
 assert len(d)==len(a) and all(d.values())
 return d
def image(A,w):
 d={}
 for j,c in w.items():
  for i,a in A[j].items():
   t=d.get(i,Q(0))+c*a
   if t:d[i]=t
   else:d.pop(i,None)
 return d
def main():
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 paths={
  'primary':ROOT/'reports'/PRIMARY/'report.json',
  'fixed_certificate':ROOT/'reports'/PRIMARY/'actual-full-G3linear-active-fixed-rank-and-kernel.json',
  'source_gate':ROOT/'reports'/GATE/'report.json',
  'gate_data':ROOT/'reports'/GATE/'actual-full-G3-linear-zero3-split-and-active-kernel-input.json',
  'source_peer':ROOT/'reports'/(GATE+'_REVIEW')/'report.json'}
 P=read(paths['primary']);G=read(paths['source_gate']);peer=read(paths['source_peer'])
 assert P['all_passed'] and G['all_passed'] and peer['all_passed']
 for k,v in P['input_sha256'].items():assert sha(ROOT/k)==v
 assert sha(paths['fixed_certificate'])==P['evidence_sha256'][paths['fixed_certificate'].name]
 assert sha(ROOT/'code/ambient_B18_ell2_T2_full_G3_linear_kernel_20261105_v63.py')==P['code_sha256']
 assert sha(ROOT/'research/notes'/f'{PRIMARY}.md')==P['note_sha256']
 assert sha(paths['gate_data'])==G['evidence_sha256'][paths['gate_data'].name]
 D=read(paths['gate_data']);C=read(paths['fixed_certificate'])
 A=[sp(c) for c in D['actual_active_G3linear252_by238_columns']]
 assert len(A)==238 and all(0<=i<252 for c in A for i in c)
 assert D['all418_fixed_diagonal_projection_identities']
 assert D['zero_variable_coordinates3']==[2,3,6]
 selected=C['selected_original_source_columns'];S=set(selected)
 assert len(S)==62 and selected==sorted(selected)
 pivots=[];product=Q(1);basis=[]
 for j,record in zip(selected,C['normalized_basis_and_original_source_representations'],strict=True):
  receipt=record['receipt'];p=receipt['pivot_row'];c=Q(receipt['divide_by'])
  v=sp(record['normalized']);w=sp(record['original238_source_representation'])
  assert receipt['source_column']==j and c and min(v)==p and v[p]==1
  assert p not in pivots and set(w)<=S and max(w)==j and w[j]==1/c
  assert image(A,w)==v
  pivots.append(p);product*=c;basis.append(v)
 assert pivots==C['fixed_original_minor_rows']
 assert [str(Q(q['receipt']['divide_by'])) for q in C['normalized_basis_and_original_source_representations']]==C['fixed_minor_determinant_product']
 # In increasing leading-row order the normalized minor is triangular with
 # diagonal 1. Simultaneously reorder its columns, so its determinant stays 1.
 order=sorted(range(62),key=lambda j:pivots[j])
 assert all(basis[order[j]].get(pivots[order[i]],0)==0 for j in range(62) for i in range(j))
 free=[]
 for record in C['full_constant_kernel_columns']:
  j=record['free_source_index'];w=sp(record['original238_source_combination'])
  assert j not in S and w.get(j)==1 and set(w)<=(S|{j})
  assert max(w)==j and not image(A,w)
  free.append(j)
 assert len(free)==176 and len(set(free))==176 and S|set(free)==set(range(238))
 assert C['active_rank']==P['active_rank']==62
 assert C['full_G3linear_rank']==P['full_G3linear_rank']==180+62==242
 assert C['full_constant_kernel_dimension']==P['full_constant_kernel_dimension']==176
 assert not C['Taylor2_coefficients_queried'] and not C['actual_T2_HF2_or_total_length_computed']
 evidence={
  'all62_exact_basis_memberships':True,'all176_exact_kernel_identities':True,
  'source_partition_triangular_unit_change':True,'leading_minor_and_original_minor_product_verified':True,
  'original62_minor_determinant':str(product),'active_rank':62,'full_rank':242,'full_kernel_dimension':176,
  'full_kernel_embedding':'active first238 coordinates, then180 zeros; source-peer fixed diagonal forces this',
  'no_rank_search_or_producer_import':True,'no_Taylor2_or_HF2_or_total_length_claim':True}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 ep=out/'independent-fixed-source-rank-and-kernel.json';ep.write_text(json.dumps(evidence,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'active_rank':62,'full_G3linear_rank':242,'full_constant_kernel_dimension':176,
  'no_new_rank_search_or_Taylor2':True,'input_sha256':{k:{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for k,p in paths.items()},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
 print(json.dumps({'all_passed':True,'active_rank':62,'full_rank':242,'kernel':176,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()

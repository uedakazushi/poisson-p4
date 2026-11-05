#!/usr/bin/env python3
"""Independent fixed rank/kernel postcheck for the 105-row condition."""
import hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_full_G3_linear_kernel_review_20261105_v63 import sp,image
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_T2_G3_PLANE_KERNEL_20261105_v63';STEM=A+'_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 pp=ROOT/'reports'/A/'report.json';P=json.loads(pp.read_text());assert P['all_passed'];inputs={'primary':{'path':str(pp.relative_to(ROOT)),'sha256':sha(pp)}}
 for p,h in P['input_sha256'].items():assert sha(ROOT/p)==h;inputs[p]={'path':p,'sha256':h}
 assert sha(ROOT/'code/ambient_B18_ell2_T2_G3_plane_kernel_20261105_v63.py')==P['code_sha256']
 assert sha(ROOT/'research/notes'/f'{A}.md')==P['note_sha256']
 ep=pp.parent/'actual-G3-plane-condition-fixed-rank-and-kernel.json';assert sha(ep)==P['evidence_sha256'][ep.name]
 inputs['certificate']={'path':str(ep.relative_to(ROOT)),'sha256':sha(ep)};C=json.loads(ep.read_text())
 gate=ROOT/'reports/AMBIENT_B18_ELL2_T2_G3_PLANE_KERNEL_INPUT_GATE_20261105_v63/actual378-annihilation-checks-and105by238-condition.json'
 D=json.loads(gate.read_text());Acols=list(map(sp,D['condition105_by238_columns']))
 sourcepeer=ROOT/'reports/AMBIENT_B18_ELL2_T2_G3_PLANE_KERNEL_INPUT_GATE_20261105_v63_REVIEW/report.json';assert json.loads(sourcepeer.read_text())['all_passed']
 assert len(Acols)==238 and C['condition_shape']==[105,238]
 selected=C['selected_original_source_columns'];S=set(selected);assert len(S)==20 and selected==sorted(selected)
 pivots=[];norms=[];det=Q(1)
 for j,record in zip(selected,C['normalized_basis_and_original_source_representations'],strict=True):
  rec=record['receipt'];p=rec['pivot_row'];c=Q(rec['divide_by']);v=sp(record['normalized']);w=sp(record['original238_source_representation'])
  assert rec['source_column']==j and c and min(v)==p and v[p]==1 and p not in pivots
  assert set(w)<=S and max(w)==j and w[j]==1/c and image(Acols,w)==v
  pivots.append(p);norms.append(v);det*=c
 assert pivots==C['fixed_original_minor_rows']
 assert [r['receipt']['divide_by'] for r in C['normalized_basis_and_original_source_representations']]==C['fixed_minor_determinant_product']
 order=sorted(range(20),key=lambda i:pivots[i]);assert all(not norms[order[j]].get(pivots[order[i]],0) for j in range(20) for i in range(j))
 free=[]
 for record in C['full_constant_kernel_columns']:
  j=record['free_source_index'];w=sp(record['original238_source_combination'])
  assert j not in S and w.get(j)==1 and set(w)<=(S|{j}) and max(w)==j
  assert not image(Acols,w);free.append(j)
 assert len(free)==len(set(free))==218 and S|set(free)==set(range(238))
 assert C['condition_rank']==P['condition_rank']==20 and C['enlarged_constant_kernel_dimension']==P['enlarged_constant_kernel_dimension']==218
 assert not C['Taylor2_coefficients_queried'] and not C['actual_T2_HF2_or_total_length_computed']
 data={'all20_source_membership_identities':True,'all218_exact_kernel_identities':True,'source_partition_and_triangular_unit_change':True,
       'normalized_leading_minor_determinant':1,'original20_minor_determinant':str(det),'rank':20,'kernel_dimension':218,
       'scope':'Only the active105-condition family. Actual quadratic classes need the separate corrected182-span source formula; no HF2/total length or unrestricted393 conclusion.',
       'no_rank_search_or_producer_import':True}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);path=out/'independent-plane-condition-rank20-and-kernel218.json';path.write_text(json.dumps(data,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'condition_rank':20,'enlarged_constant_kernel_dimension':218,'rank_search_or_quadratic_queries':False,
         'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
         'evidence_sha256':{path.name:sha(path)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
 print(json.dumps({'all_passed':True,'rank':20,'kernel':218,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()

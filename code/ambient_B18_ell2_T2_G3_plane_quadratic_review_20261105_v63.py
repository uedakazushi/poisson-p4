#!/usr/bin/env python3
"""Independent all218 signed-vector postcheck; no new rank search."""
import json,resource,signal,time
from pathlib import Path
from ambient_B18_ell2_T2_G3_quadratic_source_review_20261105_v63 import sha,sp,add,scale
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_T2_G3_PLANE_QUADRATIC_20261105_v63';STEM=A+'_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 pp=ROOT/'reports'/A/'report.json';P=json.loads(pp.read_text());assert P['all_passed'];inputs={'primary':{'path':str(pp.relative_to(ROOT)),'sha256':sha(pp)}}
 for p,h in P['input_sha256'].items():assert sha(ROOT/p)==h;inputs[p]={'path':p,'sha256':h}
 assert sha(ROOT/'code/ambient_B18_ell2_T2_G3_plane_quadratic_20261105_v63.py')==P['code_sha256']
 assert sha(ROOT/'research/notes'/f'{A}.md')==P['note_sha256']
 ep=pp.parent/'actual218-G3-quadratic-initial-classes-and-fixed-rank.json';assert sha(ep)==P['evidence_sha256'][ep.name];E=json.loads(ep.read_text());inputs['certificate']={'path':str(ep.relative_to(ROOT)),'sha256':sha(ep)}
 ip=ROOT/'reports/AMBIENT_B18_ELL2_T2_G3_PLANE_QUADRATIC_INPUT_GATE_20261105_v63/actual-enlarged-G3quadratic-factored-source-and-known-NF-input.json';D=json.loads(ip.read_text())
 Acols=list(map(sp,D['actual_enlarged_selected126_correction_linear_columns']));linearP=list(map(sp,D['actual_selected126_rawD31_Pi42_linear_columns']));NF=list(map(sp,D['PiY_L48_degree2_NF42_times20_columns']));pairindex={tuple(p):j for j,p in enumerate(D['degree2_active_pair_order20'])}
 sums=[];signedterms=0
 for correction in Acols:
  v={}
  for k,x in correction.items():
   source,u=divmod(k,6)
   for i,y in linearP[source].items():
    row,w=divmod(i,6);pair=tuple(sorted((u,w)))
    if pair==(5,5):continue
    term=NF[20*row+pairindex[pair]];signedterms+=len(term)
    v=add(v,term,-x*y)
  assert not v;sums.append(v)
 assert len(sums)==218 and signedterms==D['signed_product_summand_count']==6
 assert list(map(sp,E['actual218_projected_quadratic_initial_relation_columns']))==sums
 assert E['rank']==P['known_quadratic_initial_relation_rank']==0
 assert E['known_quadratic_cover_standard_dimension']==P['known_quadratic_cover_standard_dimension']==6
 assert all(not E[k] for k in ['selected_original_kernel_relation_indices','fixed_original_minor_rows','fixed_minor_determinant_product','fixed_rank_receipts','normalized_basis','basis_source_representations'])
 assert E['kernel_free_source_indices']==list(range(218)) and all(not q for q in E['all218_factorizations_in_selected_original_relation_basis'])
 assert E['actual_HF2_or_total_length_not_computed'] and not E['new_physical_or_implicit_coefficients']
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);path=out/'independent-all218-quadratic-zero-and-scope.json'
 path.write_text(json.dumps({'all218_signed_vectors_independently_zero':True,'signed_nonzero_factor_terms_combined':6,'rank':0,'all218_factorizations_and_identity_kernel_recipe':True,
                             'scope':'Stop this sufficient active family. No HF2, full quadratic initial generation, total length or H0 inference.'},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'all218_actual_quadratic_classes_zero':True,'actual_HF2_or_total_length_computed':False,
         'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
         'evidence_sha256':{path.name:sha(path)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
 print(json.dumps({'all_passed':True,'rank':0,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()

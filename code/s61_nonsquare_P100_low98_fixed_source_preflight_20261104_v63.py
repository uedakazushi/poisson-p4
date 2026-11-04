#!/usr/bin/env python3
"""Read-only fixed source inventory; no scalar coefficient or minor evaluation."""
import hashlib,json,collections
from pathlib import Path
R=Path(__file__).resolve().parents[1];S='S61_NONSQUARE_P100_LOW98_FIXED_SOURCE_PREFLIGHT_20261104_v63';O=R/'reports'/S
P=R/'reports/S61_NONSQUARE_P100_H15_FULL239_FINITE_TEST_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main():
 O.mkdir(exist_ok=True);inputs={};fs=[]
 for i in range(7):
  p=P/f'factor{i}.json';inputs[str(p.relative_to(R))]=sha(p);fs.append(json.loads(p.read_text()))
 rows=fs[0]['selected_source319_rows'];assert len(rows)==98 and len(set(rows))==98 and all(f['selected_source319_rows']==rows and f['rank']==98 for f in fs)
 p=R/'reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/certificate.json';inputs[str(p.relative_to(R))]=sha(p);c=json.loads(p.read_text());labels=[c['source_rows'][c['retained_row_indices'][i]] for i in rows];counts=dict(collections.Counter(x[0] for x in labels));assert counts=={'cubic_normalized':87,'quadric':11}
 assert all(x[0]!='quartic' for x in labels)
 light={0,1,2,3,4,6,7}
 for label in labels:
  if label[0]=='quadric':assert sum(label[2])==2 or label[2]==[0,0,0,0,0,1,0,0]
  else:assert label[-1] in light
 evidence={'selected_source319_rows_common_to_all7_factors':rows,'source_labels':labels,'source_counts':counts,'selected_quadric_rows':[i for i,x in zip(rows,labels) if x[0]=='quadric'],'global_lazy_matrix':'B(z)=transpose(Gamma319x98)[all98targets,selected98sources]','global_lazy_D':'det(B(z)); actual source functions, not arbitrary lifts of modp values','D_mod_all7_H_factors_nonzero_from_accepted_fixed_minors':True,'new_global_coefficient_or_determinant_values_computed':False,'no_whole239_from_finite_data':True,'method_frontier':'Need actual global D factor/unit coverage or other residual89 generation certificate; generically zero curve coker may be finite torsion away fromH15. PID kernel freeness alone does not prevent rank drop.'}
 e=O/'source-inventory.json';e.write_text(json.dumps(evidence,indent=2)+'\n')
 report={'id':S,'source_only':True,'all_static_source_checks_passed':True,'input_sha256':inputs,'evidence_sha256':{e.name:sha(e)},'source_sha256':sha(Path(__file__)),'note_sha256':sha(R/'research/notes'/f'{S}.md'),'new_coefficient_NF_rank_minor_or_source_query':False,'wholeP100_239_not_proved':True};(O/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()

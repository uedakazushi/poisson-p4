#!/usr/bin/env python3
"""Read-only TIME frontier and retained H-slice sizes, no profile replay."""
from pathlib import Path
import argparse,hashlib,json,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_ROW6_PROFILE_FRONTIER_20261105_v64'
OLD='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_ROW6_COEFFICIENT_PROFILE_20261105_v64'
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def main(a):
 assert a.freeze
 start=time.monotonic();out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);old=ROOT/'reports'/OLD;inputs={}
 def read(p):
  inputs[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 failure=read(old/'failure.json');assert failure['forecast']=='hard60s'and len(failure['completed_phase_files'])==8
 phasefacts=[];common={}
 for rec in failure['completed_phase_files']:
  p=old/rec['file'];assert sha(p)==rec['sha256'];z=read(p);n=z['nodes'][z['accumulator']];assert sum(c for h,c in n['H_slice_sizes'])==n['terms_upper']
  phasefacts.append({'file':p.name,'literal_ops':len(z['ops']),'accumulator':z['accumulator'],'support_upper':n['terms_upper'],'max_H_slice':max((c for h,c in n['H_slice_sizes']),default=0)})
  for k in ['UF10','C10_Ni56','C9_C10_alpha_i']:
   if k in z['nodes']:common[k]=z['nodes'][k]
 partial=read(old/'partial-profile.json')
 assert 'new_416_sum'in partial['nodes'];n=partial['nodes']['new_416_sum'];assert n['terms_upper']==36620
 Hmax=max(c for h,c in n['H_slice_sizes']);assert Hmax==2970
 assert{key:val['terms_upper']for key,val in common.items()}=={'UF10':951,'C10_Ni56':1762,'C9_C10_alpha_i':5794}
 controller=ROOT/'code/s61_sparse_generic5_constrained83_col53_row6_coefficient_profile_20261105_v64.py'
 helper=ROOT/'code/s61_sparse_generic5_source_packed5_profile_20261105_v64.py'
 assert sha(controller)=='374bb437964ce7706b9e262ae0da675e41ceef2dc7a20c6cddd5279ea6bec92e'
 assert sha(helper)=='66384df036330d38e8b9f50c1ac12393742dcc3906ea9b44752e2708de7c881e'
 for p in[controller,helper,Path(__file__),ROOT/f'research/notes/{ID}.md']:inputs[str(p.relative_to(ROOT))]=sha(p)
 for p,h in partial['input_sha256'].items():assert sha(ROOT/p)==h;inputs[p]=h
 ev={'all_passed':True,'read_only_no_source_or_VALUE_replay':True,'SOURCE_worker_completed':False,'fixed_row6_col53_only':True,'completed_raw_components':8,'total_raw_components':11,'completed_phase_profiles':phasefacts,'current_largest_complete_accumulator':{'node':'new_416_sum','terms_upper':36620,'H_slice_sizes':n['H_slice_sizes'],'max_H_slice':2970},'common_multiplier_profiles':{k:{'terms_upper':v['terms_upper'],'max_H_slice':max(c for h,c in v['H_slice_sizes']),'exponent_box5':v['exponent_box5']}for k,v in common.items()},'whole_VALUE_unscaled_output_profile_forecast':7*36620,'only_current_output_H_slice_7numeric_forecast':7*2970,'not_a_whole_worker_budget':True,'future_SOURCE_method_only':{'intermediate_hash':'Do not tuple-sort every intermediate; record exact cardinality,BOX,Hcounts and DAG. Hash source callback2 boundaries and each11 raw component boundary.','live_snapshot':'At every phase stream ALL live packed support sets to immutable uncompressed5byte integers, with q/M/node IDs/index; no source arithmetic lost.','prior_BYTE_history_retained':True,'old_completed_support_recovery_requires_distinct_code_and_source_budget':True,'no_new_SOURCE_or_VALUE_GO_from_this_report':True},'all167x14_87plus66_RHS_oneP290_F9_F10_complements_preserved':True}
 ep=out/'saved-phase-H-fibres-and-changed-SOURCE-method-proposal.json';ep.write_bytes(enc(ev))
 raw=sum((ROOT/p).stat().st_size for p in inputs)+ep.stat().st_size+65536;assert raw<32*1024*1024
 report={'id':ID,'all_passed':True,'source_only':True,'old_hard60_SOURCE_stop_preserved':True,'no_rerun_or_new_matrix_values':True,'not_completed_col53_profile':True,'eight_of_eleven_components_saved':True,'current_positive_support':36620,'current_max_H_slice':2970,'common_positive_supports':[951,1762,5794],'whole_VALUE_or_recovery_worker_NOT_authorized':True,'READ_and_output_RAW_upper':raw,'input_sha256':inputs,'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start}
 (out/'report.json').write_bytes(enc(report));print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',action='store_true');main(p.parse_args())

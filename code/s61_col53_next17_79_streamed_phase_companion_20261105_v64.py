"""Metadata-only complete normal/failure lifetime companion; no VALUE maps."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
ID='S61_COL53_NEXT17_79_STREAMED_PHASE_COMPANION_20261105_v64'
FINAL='reports/S61_COL53_NEXT17_79_FINALIZED_SOURCE_20261105_v64'
SOURCE='reports/S61_COL53_NEXT17_79_SOURCE_PROJECTION_20261105_v64'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def enc(d):return(json.dumps(d,separators=(',',':'))+'\n').encode()
def main():
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 pins={'code/s61_col53_next17_79_streamed_value_20261105_v64.py':'a946d904a20a78bc7f9a211ee010620ac0280035bce6565eb7b729ecbaffe808',
 'reports/S61_COL53_NEXT17_79_STREAMED_VALUE_GATE_20261105_v64/report.json':'c5823b8df2afda1ff7c4285889da3f45e06ace7d00d89fc27c72133464fa98e3',
 'reports/S61_COL53_NEXT17_79_STREAMED_VALUE_GATE_20261105_v64/streamed-whole-lifetime-and-raw-budget.json':'4ad988d00463fe00796be7103496c8ebe4d6dfe1a5d135ef597146e222145f59',
 'reports/S61_COL53_NEXT17_79_STREAM_VALUE_METHOD_SOURCE_V2_20261105_v64/report.json':'e15323af46f399a232e5ad482439b7c1e910e012ca6f85bcb6fa572cb2338c99',
 FINAL+'/report.json':'548ebf5f4078707d6ee2bb6dd1867c7161c2d32bd553ebef4ce454600c7db877',
 SOURCE+'/report.json':'8cba44d77933886756700a63fb250058256db3310272b6ca626936116c785d60',
 'reports/S61_COL53_NEXT17_79_STREAMED_VALUE_CUSTODY_20261105_v64/report.json':'e75cf2b3576812d20022139e41f7096b701830e04703265f0255094dc9d6975f'}
 for f,h in pins.items():assert sha(ROOT/f)==h
 paths=[FINAL+'/compiled-full-SLP-and-fixed-projection-plan.json',SOURCE+'/all-projected-SOURCE-node-counts.json']
 for f in paths:
  receipt=json.loads((ROOT/(FINAL if f.startswith(FINAL+'/')else SOURCE)/'report.json').read_text())
  pins[f]=sha(ROOT/f);assert pins[f]==receipt['evidence_sha256'][Path(f).name]
 D=json.loads((ROOT/paths[0]).read_text());C=dict(json.loads((ROOT/paths[1]).read_text()));assert len(C)==1884
 rows=[];maxima={};META=118000
 for k in sorted(C):
  n=D['nodes'][k];children=list(dict.fromkeys(n['args']if n['op']in['copy','multiply']else D['active_tied_HL_children'][str(k)]))
  ct=sum(C[c]for c in children);ni=C[k]
  phases={'evaluation':META+7*(ct+ni)+4096,
   'child_decode_including_2048':META+7*ct+8*max([C[c]for c in children]+[0])+2048+4096,
   'normal_codec_children_still_held':META+7*(ct+ni)+3*ni+4096}
  if n['op']=='input'and n['q']!='0':
   phases['normal_input']=META+16152+7*ni+4096
   phases['leaf_STOP_codec_retaining_input_decode']=META+16152+7*ni+3*ni+4096
  if k==D['FULL_F11_original_root']:
   assert k==1202 and children==[1196,1201]and ct==120 and ni==62
   phases['anchor_STOP_codec_retaining_children_and_anchor_decode']=META+7*(ct+ni)+8*411+2048+3*ni+4096
  assert max(phases.values())<150000
  for name,v in phases.items():
   if v>maxima.get(name,{'slots':0})['slots']:maxima[name]={'slots':v,'node':k}
  rows.append({'node':k,'children_source_terms':ct,'output_source_terms':ni,'phases':phases})
 expected={'evaluation':124805,'child_decode_including_2048':126982,'normal_codec_children_still_held':125192,'normal_input':138535,'leaf_STOP_codec_retaining_input_decode':138658,'anchor_STOP_codec_retaining_children_and_anchor_decode':128892}
 assert{name:r['slots']for name,r in maxima.items()}==expected
 normal=max(v for name,v in expected.items()if '_STOP_'not in name);prefix=max(expected.values());assert(normal,prefix)==(138535,138658)
 proof={'metadata_only_no_signed_value_queries':True,'metadata_base_unchanged':META,'all1884_node_phase_bounds':rows,'phase_maxima':maxima,'normal_success_whole_upper':normal,'complete_guard_prefix_whole_upper':prefix,'child_decoder_extra_2048_charged':True,'anchor_retained_120_child_terms_extra_slots':840,'all_strictly_below_same150000':True,'frozen_worker_gate_method_unchanged':True}
 (out/'all-normal-and-guard-prefix-phase-bounds.json').write_bytes(enc(proof))
 own='code/s61_col53_next17_79_streamed_phase_companion_20261105_v64.py';note='research/notes/'+ID+'.md';pins[own]=sha(ROOT/own);pins[note]=sha(ROOT/note)
 (out/'custody.json').write_bytes(enc({'exact_input_sha256':pins,'no_chosen_VALUE_execution':True,'old_frozen_nominal138535_not_rewritten':True,'fresh_RAW_within_old_control_reserve':True,'companion_is_metadata_only_not_new_codec_generation':True}))
 report={'id':ID,'all_passed':True,'METADATA_SOURCE_only_no_chosen_signed_VALUE':True,'normal_success_whole_upper':normal,'guard_prefix_whole_upper':prefix,'corrected_child_decode_upper':126982,'normal_codec_upper':125192,'anchor_prefix_upper_with_retained_children':128892,'global_same_guard':150000,'old_worker_gate_source_method_unchanged':True,'input_sha256':pins,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()}}
 (out/'report.json').write_bytes(enc(report));print('report_sha256',sha(out/'report.json'));print('phase_maxima',maxima)
if __name__=='__main__':main()

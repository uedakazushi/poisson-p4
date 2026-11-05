#!/usr/bin/env python3
"""SOURCE-only exact thin controls and one-leaf custody for staged col56."""
from pathlib import Path
import argparse,json,hashlib,time,signal,resource,re,gc
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_THIN_CONTROL_SOURCE_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def digest(v):return hashlib.sha256(enc(v)).hexdigest()
def scalars(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,str):return 3 if re.fullmatch(r'[+-]?[0-9]+(?:/[0-9]+)?',v)else 0
 if isinstance(v,list):return sum(map(scalars,v))
 if isinstance(v,dict):return sum(map(scalars,v.values()))
 return 0
def main(a):
 assert a.build
 begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);(out/'leaves').mkdir();(out/'phases').mkdir()
 rp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64/report.json';assert sha(rp)=='e332bea912ec182d8c6f7ed614a1cb7a25cc0dda4ab13f5d6033c85c9200b5ab';r=read(rp);assert r['all_passed']and r['actual_new_matrix_coefficients']==0
 inputs={ROOT/p:h for p,h in r['input_sha256'].items()};inputs[rp]=sha(rp)
 for name,h in r['evidence_sha256'].items():inputs[rp.parent/name]=h
 inputs.update({Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md')})
 for p,h in inputs.items():assert sha(p)==h
 pp=rp.parent/'staged-literal-col56-program-and-positive-profiles.json';program=read(pp);ops=program['literal_ops'];profiles=program['positive_support_and_height_profiles'];outs=program['row15_output_nodes'];assert len(ops)==700 and len(outs)==15;leafrefs={};maxleaf=0;custody=[];sourcepeak=0
 base=25000+scalars(program)+scalars(r)+scalars(outs);assert base<150000
 for i,op in enumerate(ops):
  if op['op']!='input':continue
  bind=op['binding'];parent=rp.parent/bind['scope'];assert parent in inputs;data=read(parent)
  if bind.get('entire_operand'):z=data
  else:
   pool=data['common_operands']if'common_operands'in data else data['row_operands'];z=pool[bind['label']]
  leaf={'all_passed':True,'parent_path':str(parent.relative_to(ROOT)),'parent_sha256':sha(parent),'original_binding':bind,'original_input_op_sha256':digest(op),'operand':z};file=out/'leaves'/f'input-{i:03d}.json';file.write_bytes(enc(leaf));back=read(file);assert back==leaf and back['operand']==z;maxleaf=max(maxleaf,scalars(leaf));leafrefs[i]={'file':str(file.relative_to(ROOT)),'sha256':sha(file),'bytes':file.stat().st_size,'leaf_parse_scalars':scalars(leaf)};custody.append({'op':i,'node':op['output'],'parent':str(parent.relative_to(ROOT)),'parent_sha256':sha(parent),'leaf':str(file.relative_to(ROOT)),'leaf_sha256':sha(file),'exact_operand_object_equal':True});sourcepeak=max(sourcepeak,base+scalars(data)+scalars(leaf)+scalars(back)+scalars(custody)+4096);assert sourcepeak<150000;del data,z,leaf,back;gc.collect()
 # Derive exact release instructions, checking every old frozen last-use flag.
 uses={k:0 for k in profiles}
 for op in ops:
  for key in('A','B'):
   if key in op:uses[op[key]]+=1
 for o in outs:uses[o['node']]+=1
 schedules=[];old=program['lastuse_schedule']
 for i,op in enumerate(ops):
  dest=op['output'];reuse=op['op']=='add'and uses[op['A']]==1;assert reuse==old[i]['reuse_left_lastuse_buffer'];release=[]
  if uses[dest]==0:release.append(dest)
  for key in('A','B'):
   if key in op:
    uses[op[key]]-=1
    if uses[op[key]]==0:release.append(op[key])
  output=next((o for o in outs if o['node']==dest),None)
  if output:uses[dest]-=1;release.append(dest)
  schedules.append({'op':i,'reuse_left_last_use_buffer':reuse,'release_after_arithmetic_before_encoding':[k for k in release if k!=dest or output is None],'release_output_after_encoding':dest if output else None,'original_release_forecast_sha256':digest(old[i])})
 assert all(v==0 for v in uses.values())
 first=outs[0]['phase_begin_op'];ranges=[{'kind':'shared','start':0,'end':first-1}]+[{'kind':'row','row':o['row'],'start':o['phase_begin_op'],'end':o['phase_end_op'],'output_node':o['node']}for o in outs];assert ranges[0]['start']==0 and ranges[-1]['end']==699
 for a,b in zip(ranges,ranges[1:]):assert a['end']+1==b['start']
 phaseindex=[];fragmentcounts=[];flatten=[]
 for phase in ranges:
  records=[]
  for i in range(phase['start'],phase['end']+1):
   op=ops[i];rec={'index':i,'original_op':op,'original_op_sha256':digest(op),'output_profile':profiles[op['output']],'output_profile_sha256':digest(profiles[op['output']]),'last_use':schedules[i]}
   if i in leafrefs:rec['input_leaf']=leafrefs[i]
   records.append(rec);flatten.append(rec['original_op_sha256'])
  frag={'all_passed':True,'source_only':True,'phase':phase,'parent_program_sha256':sha(pp),'literal_records':records};file=out/'phases'/('shared.json'if phase['kind']=='shared'else f"row-{phase['row']}.json");file.write_bytes(enc(frag));assert read(file)==frag;count=scalars(frag);fragmentcounts.append(count);phaseindex.append({'phase':phase,'file':str(file.relative_to(ROOT)),'sha256':sha(file),'bytes':file.stat().st_size,'parse_scalars':count});sourcepeak=max(sourcepeak,base+scalars(frag)+scalars(phaseindex)+scalars(custody)+4096);assert sourcepeak<150000;del frag,records;gc.collect()
 assert flatten==[digest(op)for op in ops]and len(flatten)==700
 # VALUE keeps shared metadata plus ONE phase and ONE leaf, not700records.
 controller={'all_passed':True,'source_only':True,'parent_program':str(pp.relative_to(ROOT)),'parent_program_sha256':sha(pp),'phase_index16':phaseindex,'shared_N8_35_56_node':program['literal_shared_N8_35_56_node'],'output_nodes15':outs,'future_VALUE_no_full700_program_JSON_parse':True,'future_VALUE_one_input_leaf_parse':True,'full168x15_other14_87plus66_RHS_oneP290_all_complements_kept':True};controlfile=out/'thin-controller.json';controlfile.write_bytes(enc(controller));controls=scalars(controller);metarecords=[];META=0
 for phase,count in zip(ranges,fragmentcounts):
  allowance=fragmentcounts[0]+(0 if phase['kind']=='shared'else count)+maxleaf+controls+10000;META=max(META,allowance);metarecords.append({'phase':phase,'shared_fragment_scalars':fragmentcounts[0],'current_phase_scalars':0 if phase['kind']=='shared'else count,'one_leaf_scalars_max':maxleaf,'controller_scalars':controls,'helper_runtime_and_parser_extra':10000,'whole_fixed_metadata_upper':allowance})
 META=((META+999)//1000)*1000;assert META<50000;sourcepeak=max(sourcepeak,base+scalars(controller)+scalars(metarecords)+scalars(custody)+4096);assert sourcepeak<150000
 valueopmax=max(z['whole_coefficient_owned_upper']for z in old)-50000+META;valueencmax=max(z['planned_streaming_binary_owned_upper']for z in program['streaming_binary_phase_forecasts'])-50000+META;valuepeak=max(valueopmax,valueencmax);assert valuepeak<150000
 rawcaps=sum(z['binary_RAW_upper']for z in program['streaming_binary_phase_forecasts']);inputRAW=sum(p.stat().st_size for p in inputs);thinRAW=sum(p.stat().st_size for p in out.rglob('*')if p.is_file());futureRAW=inputRAW+thinRAW+rawcaps+2*1024*1024;assert futureRAW<32*1024*1024
 proof={'all_passed':True,'source_only':True,'700_original_ops_and_profiles_exact':True,'input_leaf_custody':custody,'all16_metadata_phase_bounds':metarecords,'future_fixed_metadata_upper':META,'future_whole_coefficient_owned_upper':valuepeak,'future_operation_owned_upper':valueopmax,'future_stream_encoder_owned_upper':valueencmax,'max_single_leaf_parse_scalars':maxleaf,'raw_output_upper':rawcaps,'future_RAW_upper_including_all_READ_and2MiB_next_gate_receipts':futureRAW,'no_old_guard_or_program_changed':True,'SOURCE_no_coefficients_evaluated':True,'source_metadata_integer_strings_charged_three_scalars':True,'independent_source_gate_fresh_worker_custody_REQUIRED_before_coefficients':True};(out/'thin-custody-lifetime-and-RAW-proof.json').write_bytes(enc(proof));report={'id':ID,'all_passed':True,'source_only':True,'exact_thin_700_ops':700,'exact_input_leaves':len(leafrefs),'metadata_phases':16,'future_fixed_metadata_upper':META,'future_whole_coefficient_owned_upper':valuepeak,'future_RAW_upper':futureRAW,'SOURCE_owned_peak':sourcepeak,'zero_new_matrix_coefficients':True,'future_worker_and_independent_gate_REQUIRED':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',action='store_true');main(p.parse_args())

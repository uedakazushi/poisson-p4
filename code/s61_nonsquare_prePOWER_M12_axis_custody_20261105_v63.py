#!/usr/bin/env python3
"""Immutable literal-only M12 source custody. NO signed new axis evaluation."""
import argparse,gc,hashlib,json,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_M12_AXIS_WORKER_GATE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def key_start(s,key):
 tag=json.dumps(key)+':';assert s.count(tag)==1;i=s.index(tag)+len(tag)
 while s[i].isspace():i+=1
 return i
def skip(s,i):
 while s[i].isspace():i+=1
 if s[i]in'[{':
  depth=0;quoted=False;escape=False
  for j in range(i,len(s)):
   c=s[j]
   if quoted:
    if escape:escape=False
    elif c=='\\':escape=True
    elif c=='"':quoted=False
   elif c=='"':quoted=True
   elif c in'[{':depth+=1
   elif c in']}':
    depth-=1
    if not depth:return j+1
 return json.JSONDecoder().raw_decode(s,i)[1]
def select(s,start,indices):
 i=start
 for index in indices:
  assert s[i]=='[';i+=1
  for j in range(index+1):
   while s[i].isspace()or s[i]==',':i+=1
   if j==index:break
   i=skip(s,i)
 return json.JSONDecoder().raw_decode(s,i)[0]
def main(out):
 started=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={};texts={}
 def load(rel):p=ROOT/rel;pins[rel]=sha(p);return json.loads(p.read_text())
 def text(rel):
  if rel not in texts:p=ROOT/rel;pins[rel]=sha(p);texts[rel]=p.read_text()
  return texts[rel]
 st='reports/S61_NONSQUARE_PREPOWER_M12_AXIS_SOURCE_20261105_v63/'
 r=load(st+'report.json');s=load(st+'added-source-profiles.json');f=load(st+'whole-axis-source-worker-forecast.json');assert r['all_passed']
 ordered=load('reports/S61_NONSQUARE_PREPOWER_M12_AXIS_ORDERED_L1_GATE_20261105_v63/report.json');assert ordered['completed_source_bound_bits']==1984
 frac=load('reports/S61_NONSQUARE_PREPOWER_M12_AXIS_FRACTION_SOURCE_20261105_v63/report.json');assert frac['all_unreduced_Fraction_temporary_bits']==238 and frac['whole_prep_numeric_upper']==19212
 peer=load('reports/S61_NONSQUARE_PREPOWER_TARGET17_THREE_CRAMER_ENTRY_20261105_v63_REVIEW/report.json');assert peer['all_passed']
 assert all(f[k]==0 for k in['new_row_a_shift','new_column_a_shift','new_row_B_clear','new_column_B_clear'])
 records=[]
 with(out/'literal-atoms67.jsonl').open('w')as w:
  for i,a in enumerate(s['all_atom_profiles']):
   ref=a['reference']
   if 'constant'in ref:v=[[[[0,0,0],str(ref['constant'])]],[[[0,0,0],'1']]]
   elif 'text_span'in ref:
    ss=text(ref['path']);start,end=ref['text_span'];v,e=json.JSONDecoder().raw_decode(ss,start);assert e==end
   elif 'json_pointer'in ref:
    ss=text(ref['path']);index=int(ref['json_pointer'].split('/')[-1]);v=select(ss,key_start(ss,'new_unit_atoms_num_den'),[index])
   elif 'target'in ref:
    ss=text(ref['path']);labels=json.JSONDecoder().raw_decode(ss,key_start(ss,'original_target_labels'))[0];num=select(ss,key_start(ss,'source_h9'),[labels.index(ref['target'])]);v=[num,[[[0,0,0],'1']]]
   else:
    ss=text(ref['path']);num=select(ss,key_start(ss,'raw_compact14x5'),[ref['compact'],ref['heavy_column']]);v=[num,[[[0,0,0],a['denominator_QQ_unit']]]]
   assert all(m[1]>=0 and m[2]>=0 for m,c in v[0]);w.write(json.dumps(dict(atom=i,literal_num_den=v,source_reference=ref,expected=a))+'\n');records.append(ref)
   del v
 assert len(records)==67
 qr='reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json';qt=text(qr);roots=json.JSONDecoder().raw_decode(qt,key_start(qt,'remaining15x22_scalar_roots'))[0];cols=json.JSONDecoder().raw_decode(qt,key_start(qt,'new22_target_original32_columns'))[0]
 def atomid(definition):
  op=definition[0]
  if op=='checkpoint':matches=[i for i,v in enumerate(records)if v.get('row')==definition[1]and v.get('column')==definition[2]]
  elif op=='new_unit_atom':matches=[i for i,v in enumerate(records)if v.get('json_pointer')==f'/scalar_DAG/new_unit_atoms_num_den/{definition[1]}']
  elif op=='one':matches=[i for i,v in enumerate(records)if v.get('constant')==1]
  else:return None
  assert matches;return matches[0]
 nodes=[dict(node=v['node'],definition=v['definition'],axis_bound=v['bound']['axis'],atom=atomid(v['definition']))for v in s['finite_DAG_profiles']]
 h9idx=next(i for i,v in enumerate(records)if v.get('target')==17);h3idx=[next(i for i,v in enumerate(records)if v.get('compact')==13 and v.get('heavy_column')==j)for j in[1,3,2]]
 hr='reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json';rho=json.JSONDecoder().raw_decode(text(hr),key_start(text(hr),'known_row_units_rho14'))[0]
 for i in[6,8,9,10,11,13]:assert rho[i]['a_d_F_B_C_J_powers']==[1,0,1,0,2,0]
 assert rho[12]['a_d_F_B_C_J_powers']==[0,0,1,0,2,0]
 assert rho[13]['constant']=='8000000/16119'
 axrel='reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_SOURCE_V2_20261105_v63/whole-axis-literal-source.json';old=load(axrel);assert sum(len(v['axis_integer_polynomial'])for v in old['axis_coefficients'])==579
 dump(out/'old-axis579.json',dict(axis_coefficients=[dict(row=v['row'],column=v['column'],polynomial=v['axis_integer_polynomial'])for v in old['axis_coefficients']],old_row_LCD=old['axis_row_denominators'],old_reference=axrel,old_sha256=pins[axrel]))
 tr='reports/S61_NONSQUARE_PREPOWER_K92_ORIGINAL_TERMINAL_FACE_SOURCE_20261105_v63/terminal-face-and-normalization.json';tt=text(tr);rowdual=json.JSONDecoder().raw_decode(tt,key_start(tt,'row_dual'))[0];coldual=json.JSONDecoder().raw_decode(tt,key_start(tt,'column_dual'))[0]
 recipe=dict(nodes=nodes,current15_roots=roots,current22_original_labels=cols,h9_atom=h9idx,H3_atoms_in_order=h3idx,constants={str(i):rho[i]['constant']for i in[6,8,9,10,11,12]},
  rowdual=rowdual,coldual=coldual,extra_row_QQ_integer_scales=f['old11_row_extra_QQ_integer_scales'],new_row_LCD=f['new_row_LCD'],generic_e_d_regular_after_aB_units=True,
  old_row_operation_also_on_newcolumn='row1+2row2',exact_H3carry='rho13/L=8000000/16119; cols28,30,29 ->compact13[1,3,2]',all_minor_H=f['all_minor_L1_upper'],native_DB=108,native_CB=1096,native_RB=2048,native_owned=33275,
  source_and_extra_axis_bounds=dict(addedcol=[v['axis']for v in s['added_column11']],bottom=[v['axis']for v in s['added_G134_row11']],corner=s['corner']['axis']),
  no_kappa3_denominator=True,no_signed_new_axis_values_yet=True)
 dump(out/'recipe.json',recipe);del texts,s,old,roots,rho;gc.collect()
 cpp=ROOT/'code/s61_nonsquare_M12_axis_binary_bareiss_20261105_v63.cpp';parser=ROOT/'code/s61_nonsquare_M12_binary_inspect_20261105_v63.py';driver=ROOT/'code/s61_nonsquare_prePOWER_M12_axis_entry_20261105_v63.py'
 for p in[cpp,parser,driver,Path(__file__)]:pins[str(p.relative_to(ROOT))]=sha(p)
 exe=ROOT/'work/d-zero-v63/M12axis/worker';exe.parent.mkdir(parents=True,exist_ok=True);c=subprocess.run(['g++','-std=c++17','-O2',str(cpp),'-lflint','-lgmp','-o',str(exe)],capture_output=True,text=True,timeout=20);(out/'compile.stdout').write_text(c.stdout);(out/'compile.stderr').write_text(c.stderr);assert c.returncode==0
 h=subprocess.run([str(exe),'--help'],capture_output=True,text=True,timeout=3);assert h.returncode==0;(out/'native-help.txt').write_text(h.stdout)
 fixture=out/'format-fixture.bin';z=subprocess.run([str(exe),'--binary-fixture',str(fixture)],capture_output=True,text=True,timeout=3);assert z.returncode==0
 z=subprocess.run([str(ROOT/'.venv/bin/python'),str(parser),str(fixture)],capture_output=True,text=True,timeout=3);assert z.returncode==0;(out/'format-fixture-parse.json').write_text(z.stdout)
 note=ROOT/'research/notes/S61_NONSQUARE_PREPOWER_M12_AXIS_ENTRY_20261105_v63.md';pins[str(note.relative_to(ROOT))]=sha(note)
 prepraw=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);guardraw=3000000;mainraw=7899784;total=external+prepraw+guardraw+mainraw;assert total<=22239243
 gate=dict(id=ID,all_passed=True,frozen=True,mathematical_axis_job_NOT_run=True,source_literal_atoms=67,source_DAG_nodes=214,known579_literal_reuse=True,
  Fraction_source_bound238=True,Bareiss_source_bound1984=True,hard_bit_guard2048=True,whole_source_prep_numeric_upper=19212,native_whole_owned_upper=33275,
  ALL_numeric_key_metadata_cap=33275,source_runtime_metadata_reserve=11000,native_ALL_metadata_reserve=15000,main_RAW_upper=mainraw,native_guard_BINARY_upper=guardraw,
  actual_external_input_RAW=external,prepared_source_RAW=prepraw,TOTAL_RAW_upper=total,unchanged_RAW_cap=22239243,
  whole_source67_e_d_nonnegative=True,H3_carry_rho13_over_L='8000000/16119',no_kappa3_inverse=True,
  fixed12_axis_only=True,soft45_hard60seconds=True,memory_bytes=2000000000,
  binary_parser_fixture_NEG_ZERO_2048_PASS=True,source_input_sha256={p.name:sha(p)for p in out.iterdir()},input_sha256=pins,executable_sha256=sha(exe),seconds=time.monotonic()-started)
 dump(out/'report.json',gate);print(json.dumps({k:v for k,v in gate.items()if k not in['source_input_sha256','input_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

#!/usr/bin/env python3
"""Freeze the unrun HB_pivot03 worker, import-only checks and whole budget."""
import ast,hashlib,importlib.util,json,re,subprocess,time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_FULL_HB_PIVOT03_VALUE_CUSTODY_20261105_v63'
WORKER='code/ambient_B18_C_full_HB_pivot03_cache_values_20261105_v63.py'
SOURCE='reports/AMBIENT_B18_C_FULL_HB_PIVOT03_CACHE_SOURCE_20261105_v63/report.json'
COMPILED='reports/AMBIENT_B18_C_FULL_HB_PIVOT03_CACHE_SOURCE_20261105_v63/compiled-six-pivot-even-roots-and-pure-ell4-inputs.json'
PROFILE='reports/AMBIENT_B18_C_FULL_HB_PIVOT03_CACHE_SOURCE_20261105_v63/full-HB_pivot03-source-frame-and-whole-operational-profile.json'
FORMAL='reports/AMBIENT_B18_C_FULL_HB_PIVOT03_SOURCE_20261105_v63/full-HB6-by2-affine-source-operator-and-old27-leaves.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,(int,float)):return 1
 if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.values()))
 return 0

def main():
 t=time.monotonic();worker=ROOT/WORKER;tree=ast.parse(worker.read_text())
 source=json.loads((ROOT/SOURCE).read_text());c=json.loads((ROOT/COMPILED).read_text());p=json.loads((ROOT/PROFILE).read_text());f=json.loads((ROOT/FORMAL).read_text())
 spec=importlib.util.spec_from_file_location('unrun_HB_pivot03_worker',worker);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 assert m.SOURCE_SHA==sha(ROOT/SOURCE)
 assert not(ROOT/'reports'/m.STEM).exists()
 # Unrelated guard fixtures: no source root is composed here.
 x=list(map(m.parse,['1','2/3','-3/7','4/11']));v=m.inverse(x)
 assert m.conv(x,v)==[1,0,0,0]==m.conv(v,x)
 assert m.parse('1/4')!=m.parse('1')
 try:m.times(m.Q(1<<2040),m.Q(1<<20))
 except AssertionError:pass
 else:raise AssertionError('before-operation hard guard failed')
 help_result=subprocess.run([str(ROOT/'.venv/bin/python'),str(worker),'--help'],capture_output=True,text=True,timeout=10)
 assert help_result.returncode==0 and'--source-peer-report'in help_result.stdout
 # The eventual peer numeric body is bounded independently at5000 slots;
 # dynamic worker repeats the complete exact count before any root values.
 owned=count(source)+count(c)+count(p)+count(f)+5000+2*7095+525+4*162+16+14000
 assert owned<75000<100000<150000
 raw=sum((ROOT/k).stat().st_size for k in[SOURCE,COMPILED,PROFILE,FORMAL])
 assert raw+80000+150000<1000000
 data={
  'worker_sha256':sha(worker),'source_sha256':sha(ROOT/SOURCE),'compiled_sha256':sha(ROOT/COMPILED),'profile_sha256':sha(ROOT/PROFILE),'formal_operator_sha256':sha(ROOT/FORMAL),
  'unrun_value_folder_absent':True,'import_only_help_passed':True,'unrelated_inverse_numerator_denominator_and_overflow_fixtures':True,
  'requested_HB_pivot03_values_computed':False,
  'source_proven_before_operation_bits':p['max_operational_before_add_mul_bound'],'root_preop_bits':p['root_phase_before_operation_bound'],'post_preop_bits':p['post_affine_rank_inverse_operational_bound'],'hard_bits':2048,
  'owned_with5000_future_peer_numeric_slots':owned,'static_numeric':75000,'hard_numeric':100000,'global_numeric':150000,
  'existing_RAW_bytes':raw,'future_peer_RAW_reserve':80000,'new_output_RAW_reserve':150000,'total_fresh_RAW_cap':1000000,
  'seconds_cap':60,'memory_bytes_cap':2147483648,
  'AST_functions':[n.name for n in tree.body if isinstance(n,ast.FunctionDef)],
  'future_command':'.venv/bin/python '+WORKER+' --source-peer-report reports/AMBIENT_B18_C_FULL_HB_PIVOT03_CACHE_SOURCE_20261105_v63_REVIEW/report.json',
  'one_same_cache_job_only':'After fresh SOURCE/code peer and parent conditional GO. No IFT, gauge, new physical source or jet.',
  'scope':'FullHB6x2 knownsubmodule: rank2 givesHF4=4/TorM2(2,2)/e3or4; rank1 onlypartialnegative. No literalell3Hcut/newsource3/delta4.'}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'unrun-worker-import-guard-and-whole-budget-custody.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'UNRUN_NEW_HB_PIVOT03_WORKER_FROZEN_NO_VALUES',
  'input_sha256':{k:sha(ROOT/k)for k in[WORKER,SOURCE,COMPILED,PROFILE,FORMAL]},'code_sha256':sha(Path(__file__)),
  'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()

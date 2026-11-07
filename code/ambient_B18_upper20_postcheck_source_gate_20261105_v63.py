#!/usr/bin/env python3
"""Counts/custody only for independent fixed20 postcheck; never executes it."""
import ast,hashlib,json,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_UPPER20_POSTCHECK_INPUT_GATE_20261105_v63'
WORKER='code/ambient_B18_upper20_independent_postcheck_20261105_v63.py'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple,set)):return sum(map(count,v))
 if isinstance(v,str):return int(all(x.lstrip('-').isdigit() for x in v.split('/')) and len(v.split('/'))<=2)
 return int(isinstance(v,(int,float)))
def main():
 start=time.monotonic();out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);pins={};raw=0
 def read(p):
  nonlocal raw
  q=ROOT/p;pins[p]=sha(q);raw+=2*q.stat().st_size;return json.loads(q.read_text())
 data='reports/AMBIENT_B18_ELL2_T2_UPPER20_WHOLE_INPUT_GATE_V4_20261105_v63/whole20-coupled-pruned-cache-input-and-literal-owner-height.json'
 pp='reports/AMBIENT_B18_ELL2_T2_UPPER20_VALUES_V4_20261105_v63/'
 d=read(data);pr=read(pp+'report.json');sv=read(pp+'actual20-full-components-upper2-six-and-shared-tables.json');rank=read(pp+'actual20-fixed-integer-rank-all-factors-and-kernel.json')
 assert pins[data]=='dd8c9d761d6954cee49150166e11fa76c50071702f6450857c1523decf7e9342' and pins[pp+'report.json']=='253aff768ad00dfe9bdbfb5771bc2151fa2006ea472567697dbbdf8912a943a0'
 assert pr['all_passed'] and len(sv['all20'])==20 and count(d)==11845 and count(sv)==3144
 assert d['committed_prior']==354 and d['temporary_prior']==1063
 pins['reports/AMBIENT_B18_ELL2_T2_UPPER20_WHOLE_V4_20261105_v63_REVIEW/report.json']=sha(ROOT/'reports/AMBIENT_B18_ELL2_T2_UPPER20_WHOLE_V4_20261105_v63_REVIEW/report.json')
 base=2*count(d)+2*count(sv);assert base==29978
 native={'Nparsed8_by7needed_by9_by4_plus56keys':2072,'G576_nu36':612,'Gamma_sparse36_Ddense1120':1156,'L90dense32_plus90tuplekeys4':3240,'oneSc18_by32':576,'nativebasis_indices_and_tables':666,'aux16_points16_kernel25':177,'current_work':500,'metadata':1000}
 old={'leaf48_dense21':1008,'T118_dense21_andkeys':2596,'OLD21_two_tables':126,'B21indices_table':840,'GammaD':1156,'native_retained':177,'work':500,'metadata':1000}
 first={'GammaD':1156,'parsed_L0Le':375,'FIRST62_dense56_plus_tuplekeys':3596,'raw21_components':189,'all1742_FIRST_and21keys':1763,'OLDtables':126,'native_retained':177,'B21indices':840,'results20':420,'Q18':18,'metadata':1000}
 phase=max(sum(native.values()),sum(old.values()),sum(first.values()));assert phase==9999
 # R/N/G/nu are released afternative phase;C/T/leaves afterOLD phase.
 # Maximum NEW guard state occursbeforefirst.clear:3596+1763+189+
 # native177+OLD126+overcountresults420+work500=6771. Extra2codecs
 # plusmetadata1000 fit15000. SUCCESS/GUARD madeexclusive byTBrelease.
 guard_new=3596+1763+189+177+126+420+500;assert guard_new==6771 and 2*guard_new+1000<15000
 forecast=base+phase+15000+3000;assert forecast==57977<60000
 assert raw+2_500_000<4_000_000<32*1024*1024
 tree=ast.parse((ROOT/WORKER).read_text());allowed={'argparse','gc','hashlib','json','resource','signal','time','traceback','fractions','pathlib'}
 for node in ast.walk(tree):
  if isinstance(node,ast.Import):assert all(x.name in allowed for x in node.names)
  elif isinstance(node,ast.ImportFrom):assert node.module in allowed
 run=subprocess.run([sys.executable,str(ROOT/WORKER),'--help'],capture_output=True,timeout=10,cwd=ROOT);assert run.returncode==0 and not run.stderr
 evidence={'SOURCE_counts_only_no_signed_postcheck':True,'input_numeric_whole':count(d),'provided_primary_numeric':count(sv),'base_two_copies':base,'native_phase_ledger':native,'OLD_phase_ledger':old,'FIRST_phase_ledger':first,'explicit_release_order':['R/N/G/nu/Sc afterphase1native','C/T/leaves afterphase2OLD','FIRST62 afterraw21phase'],
 'maximum_phase_derived':phase,'NEW_guard_max6771':guard_new,'two_extra_codec_andmetadata_upper':2*guard_new+1000,'codec_reserve':15000,'output_metadata_reserve':3000,'forecast_numeric':forecast,'private_numeric_cap':60000,'RAW_cap':4000000,
 'coefficient_prior354_temp1063_inherited':True,'L90_extra22_point_terms':'68pointfreeL plusGammaL17 andD10L2point terms, bothchannels. Required to matchFORMALraw16chosenpartialN beforetruekernelrestriction; nounknownN3values.',
 'no_actual_raw_aux_claim':True,'no_new_cases_or_rank_search':'Onlyprovided20andfixedrank0cert arepostverified; unexpectednonzerostops. Remaining12opaque24seconds notqueried.',
 'scientific_scope':'all64FORMALaux scalars,15actual+10SOURCEzero kernels,21FORMALOLDpartials,1742FIRST,20fivepieces+full6 andprovidedzero2x20/I20kernel',
 'runtime':'55s math;clearframes/context/cause/gc then5s prefix;ROOTexternalhardtimeout60 andAS2GB required.',
 'fresh_ROOT_GO_before_worker':True}
 ep='source-lifetime-and-fixed-input-counts-no-responses.json';(out/ep).write_text(json.dumps(evidence,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'SOURCE_only_no_postcheck_execution':True,'private_numeric_cap':60000,'derived_forecast57977':forecast,'RAW_cap':4000000,'prior':[354,1063],'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'worker_sha256':sha(ROOT/WORKER),'code_sha256':sha(Path(__file__)),'seconds':time.monotonic()-start};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
 cp=ROOT/'reports'/'AMBIENT_B18_UPPER20_POSTCHECK_PREPARED_CUSTODY_20261105_v63';cp.mkdir(exist_ok=False)
 custody={'id':cp.name,'all_passed':True,'worker_sha256':sha(ROOT/WORKER),'source_gate_sha256':sha(out/'report.json'),'input_sha256':{**pins,str((out/'report.json').relative_to(ROOT)):sha(out/'report.json'),WORKER:sha(ROOT/WORKER)},'safe_help_returncode':run.returncode,'safe_help_stdout_sha256':hashlib.sha256(run.stdout).hexdigest(),'stdlib_AST_PASS':True,'private_numeric60000_forecast57977':True,'RAW_cap4000000':True,'source354_temp1063':True,'all_signed_postcheck_operations_UNRUN':True,'ROOTfreshGO_required':True};(cp/'report.json').write_text(json.dumps(custody,separators=(',',':'))+'\n')
 print(json.dumps({'gate_sha256':sha(out/'report.json'),'worker_sha256':sha(ROOT/WORKER),'custody_sha256':sha(cp/'report.json'),'forecast':forecast,'RAW':4000000}))
if __name__=='__main__':main()

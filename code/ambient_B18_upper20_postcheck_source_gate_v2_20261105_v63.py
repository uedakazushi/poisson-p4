#!/usr/bin/env python3
"""Counts/custody only for independent fixed20 postcheck; never executes it."""
import ast,hashlib,json,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_UPPER20_POSTCHECK_INPUT_GATE_V2_20261105_v63'
WORKER='code/ambient_B18_upper20_independent_postcheck_v2_20261105_v63.py'
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
 old_pins={'code/ambient_B18_upper20_independent_postcheck_20261105_v63.py':'2c5b6d0fbd4dd2bdbd1479d6dd48cefa476175a0d918dd0937287d50adcd79e3','reports/AMBIENT_B18_UPPER20_POSTCHECK_INPUT_GATE_20261105_v63/report.json':'80d37a75fcd6851f034c16fb90087ca9859661757f7b93c31c898364768285d0','reports/AMBIENT_B18_UPPER20_POSTCHECK_PREPARED_CUSTODY_20261105_v63/report.json':'b8fef8d2171f41148495b231a1f2394210fc4997cd05a99975f6722e9d924078'}
 for p,h in old_pins.items():assert sha(ROOT/p)==h;pins[p]=h;raw+=2*(ROOT/p).stat().st_size
 # This is a structural count, not any FIRST response pairing.
 pair_order=d['OLD']['pair_order10'];factor_counts=[]
 for case,entries in d['OLD']['original_full_a1_raw21']:
  size=sum(int(slot%9==x)+int(slot%9==y and x!=y) for slot,q in entries for x,y in pair_order);factor_counts.append([case,size])
 assert max(n for case,n in factor_counts)==83
 control_toy={'last_binary':'1','dot':'1','dot_cursor':{'index':55,'z':'1'},'poly':['1']*21,'current_vector':['1']*56,'env':['1',['1','1']],'science':{'OLD_Gamma':'1','OLD_D':'1','FIRSTcarry':'1','outer':'1','FIRSTfv':['1']*83}}
 assert count(control_toy)<200<500
 assert json.loads(json.dumps(control_toy))==control_toy
 assert raw+2_500_000<4_000_000<32*1024*1024
 tree=ast.parse((ROOT/WORKER).read_text());assert not any(isinstance(v,ast.Name) and isinstance(v.ctx,ast.Store) and v.id=='gc' for v in ast.walk(tree));allowed={'argparse','gc','hashlib','json','resource','signal','time','traceback','fractions','pathlib'}
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
 'fresh_ROOT_GO_before_worker':True,'V2_control_only':True,'all_V1_bytes_preserved_UNRUN':old_pins,'maximum_current_FIRSTfv83':factor_counts,'current_work_count_toy_upper':count(control_toy),'work500_still_covers_all_new_science_state':True,'control_changes':['current_R_z_and_inner_N_nu_cursor','native_Gamma_D_accumulators_and_dot_index','intrinsic_globalSc_cursor_no_duplicate576','actual_kernel_point_and_Gamma_D_partials','FIRST_recipe_partial56','rawcase_four_accumulators_and_FIRSTfv83','current_C0_TYPED_vector_and_product_coordinates','corrected20_partialvectors_and_nativeaccumulator','rename_OLD_Gamma_local_gc_to_gcols_to_preserve_garbage_collector_in_failure']}
 ep='source-lifetime-and-fixed-input-counts-no-responses.json';(out/ep).write_text(json.dumps(evidence,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'SOURCE_only_no_postcheck_execution':True,'private_numeric_cap':60000,'derived_forecast57977':forecast,'RAW_cap':4000000,'prior':[354,1063],'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'worker_sha256':sha(ROOT/WORKER),'code_sha256':sha(Path(__file__)),'seconds':time.monotonic()-start};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
 cp=ROOT/'reports'/'AMBIENT_B18_UPPER20_POSTCHECK_V2_PREPARED_CUSTODY_20261105_v63';cp.mkdir(exist_ok=False)
 custody={'id':cp.name,'all_passed':True,'worker_sha256':sha(ROOT/WORKER),'source_gate_sha256':sha(out/'report.json'),'input_sha256':{**pins,str((out/'report.json').relative_to(ROOT)):sha(out/'report.json'),WORKER:sha(ROOT/WORKER)},'safe_help_returncode':run.returncode,'safe_help_stdout_sha256':hashlib.sha256(run.stdout).hexdigest(),'stdlib_AST_PASS':True,'private_numeric60000_forecast57977':True,'RAW_cap4000000':True,'source354_temp1063':True,'all_signed_postcheck_operations_UNRUN':True,'ROOTfreshGO_required':True};(cp/'report.json').write_text(json.dumps(custody,separators=(',',':'))+'\n')
 print(json.dumps({'gate_sha256':sha(out/'report.json'),'worker_sha256':sha(ROOT/WORKER),'custody_sha256':sha(cp/'report.json'),'forecast':forecast,'RAW':4000000}))
if __name__=='__main__':main()

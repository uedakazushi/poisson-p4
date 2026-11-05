#!/usr/bin/env python3
"""Freeze prepared HEAD23 worker custody, never run its actual main."""
import ast,hashlib,json,subprocess,sys,time
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_CONNECTING_HEAD23_VALUE_CUSTODY_20261105_v63'
W='code/ambient_B18_C_connecting_head23_cached_values_20261105_v63.py'
C='reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHE_SOURCE_20261105_v63'
H='reports/AMBIENT_B18_C_CONNECTING_HEAD23_SOURCE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();worker=ROOT/W;ast.parse(worker.read_text());actual=ROOT/'reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63';assert not actual.exists()
 help=subprocess.run([sys.executable,str(worker),'--help'],capture_output=True,text=True,check=True)
 fixtures=subprocess.run([sys.executable,str(worker),'--self-test'],capture_output=True,text=True,check=True)
 assert 'PURE_GUARD_AND_B6_INVERSE_FIXTURES_PASS_NO_ACTUAL_SOURCE_VALUES'in fixtures.stdout and not actual.exists()
 paths=[f'{C}/report.json',f'{C}/compiled-single-bB63-root-and-provided-joint105-cache.json',f'{C}/single-root-head23-coherent-positive-height-and-whole-budget.json',
  f'{H}/report.json',f'{H}/full-coherent-head23-two-socle-cache-polynomials-and-low-two-row-source.json']
 d=[json.loads((ROOT/p).read_text())for p in paths];profile=d[2]
 assert profile['whole_before_operation_bound']==207 and profile['future_static_owned_upper']==31087
 base=S.count(d)+2000+6000+1000 # two peer reports + own custody fields
 root_live=2*2194+105*7+105*7+7*122+500
 root_upper=base+root_live+200 # phase1 snapshot/serialization
 post_terms=sum(len(p)for p in d[4]['head1_genuine_high_A_B_UNPAIRED']+d[4]['head1_canonical_socle_coordinates_rg_k2_UNPAIRED']+[p for pair in d[4]['five_head_low_rows_rH_ellH_UNPAIRED']for p in pair])
 post_upper=base+S.count(d[4]['all_input_values_already_cached'])*2+12*post_terms+1000
 bound=max(root_upper,post_upper);assert bound<=40000
 raw_base=sum((ROOT/p).stat().st_size for p in paths)+worker.stat().st_size
 raw_upper=raw_base+200000+100000+50000;assert raw_upper<2000000
 record={'id':STEM,'all_passed':True,'scope':'PREPARED/unrun worker code andwhole-lifetime custody only, not source mathematical orvalue acceptance',
  'prepared_worker_path':W,'prepared_worker_sha256':sha(worker),'source_code_helpers':'stdlibonly in VALUE worker',
  'input_sha256':{p:sha(ROOT/p)for p in paths},'source_root_records':2194,'lastuse_nodes':120,
  'strict_preoperation_bounds':{'root':83,'post':207,'hard':2048},'source_stored_bounds':{'root':41,'post_including_minor':103},
  'base_numeric_including_two_peer_and_custody_reserves':base,'root_numeric_upper':root_upper,'post_numeric_upper':post_upper,
  'whole_positive_owned_upper':bound,'static_owned_cap':40000,'hard_owned_cap':60000,'global_cap':150000,
  'RAW_base_files':raw_base,'whole_fresh_RAW_positive_upper':raw_upper,'fresh_RAW_cap':2000000,'global_RAW_cap':33554432,
  'seconds':60,'memory_bytes':2147483648,'all14_outputs_termwise_retained':True,
  'help_result':help.stdout,'guard_inverse_fixture_result':fixtures.stdout,'actual_value_folder_absent':True,
  'required_peers':['freshfullHEAD23coherentsourcepeer','freshsinglecache/code/custodywholepeer'],'ROOT_ONE_GO_required':True,
  'failure_semantics':'Strictbefore-operationQQSTOP; ifpostfails,theactualrootphase1 remainsfrozen.Noautomaticretryorrecompute.',
  'positive_scope':'Ifsimpletargeth0/h1minor!=0,ker(M7→G4)J2+k. Otherheadhigheredgecoords/originalH0/nilcubestayformal.',
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'seconds_for_custody':time.monotonic()-t}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);(out/'prepared-worker-code-source-custody.json').write_text(json.dumps(record,indent=2)+'\n')
 (out/'report.json').write_text(json.dumps({'id':STEM,'all_passed':True,'prepared_worker_sha256':sha(worker),
  'custody_sha256':sha(out/'prepared-worker-code-source-custody.json'),'code_sha256':sha(Path(__file__)),'actual_values_not_run':True},indent=2)+'\n')
 print(json.dumps({'custody_sha256':sha(out/'prepared-worker-code-source-custody.json'),'whole_positive_owned_upper':bound,'RAW_upper':raw_upper}))
if __name__=='__main__':main()

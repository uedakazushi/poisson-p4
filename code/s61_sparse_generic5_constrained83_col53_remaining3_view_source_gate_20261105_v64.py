#!/usr/bin/env python3
"""Fresh SOURCE gate for snapshot07 resume, support views and fused SOURCE union."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,ast,hashlib,json,os,resource,signal,subprocess,sys,time
from s61_sparse_generic5_source_view_profile_20261105_v64 import Engine
from s61_sparse_generic5_source_packed5_profile_v2_20261105_v64 import pack,unpack,rawindex
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_REMAINING3_VIEW_SOURCE_GATE_20261105_v64'
WORKER='code/s61_sparse_generic5_constrained83_col53_remaining3_view_profile_20261105_v64.py'
HELPER='code/s61_sparse_generic5_source_view_profile_20261105_v64.py'
PREV='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_ROW6_COEFFICIENT_PROFILE_V2_20261105_v64'
SCALAR='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_REMAINING3_SCALAR_SOURCE_V2_20261105_v64'
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def main(a):
 assert a.freeze
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.signal(signal.SIGALRM,lambda s,f:(_ for _ in()).throw(TimeoutError('gate60s')));signal.alarm(60)
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={}
 def add(p):
  p=Path(p);inputs[str(p.relative_to(ROOT))]=sha(p);return p
 def read(p):return json.loads(add(p).read_text())
 try:
  sp=add(ROOT/'reports'/SCALAR/'report.json');assert sha(sp)=='d3da153ee4aa9276783c947e7634740c65a73139ecd6d29f9dc3deec66834664';sr=read(sp);assert sr['all_passed']
  scalar=read(sp.parent/'old-live07-and-remaining93-scalar-H-majorants.json')
  assert scalar['max_positive_M_bits_upper']==404 and scalar['max_QQ_field_bits_upper']==80 and scalar['max_H_BOX_width_upper']==26
  previous=ROOT/'reports'/PREV;fail=read(previous/'failure.json');assert fail['forecast']==155866 and len(fail['completed_phase_files'])==8
  phase=read(previous/'phase-07.json');shot=previous/phase['complete_live_snapshot']['file'];assert sha(add(shot))==phase['complete_live_snapshot']['sha256']
  assert phase['complete_live_snapshot']['total_keys']==44555 and phase['complete_live_snapshot']['maps']==22
  for p,h in read(previous/'partial-profile.json')['input_sha256'].items():assert sha(add(ROOT/p))==h
  # Exact history closure: no preceding source code, snapshot, failed prefix or
  # actual oldREAD file is silently dropped from this fresh RAW budget.
  for stem in[PREV,'S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_ROW6_COEFFICIENT_PROFILE_20261105_v64','S61_SPARSE_GENERIC5_CONSTRAINED83_TENTH_COL53_STAGED_SOURCE_20261105_v64']:
   for p in(ROOT/'reports'/stem).rglob('*'):
    if p.is_file():add(p)
  for p,h in sr['input_sha256'].items():assert sha(add(ROOT/p))==h
  for f in[WORKER,HELPER,'code/s61_sparse_generic5_source_packed5_profile_v2_20261105_v64.py','code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py',str(Path(__file__).relative_to(ROOT)),f'research/notes/{ID}.md',f'research/notes/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_REMAINING3_VIEW_PROFILE_20261105_v64.md']:add(ROOT/f)
  for f in[WORKER,HELPER]:ast.parse((ROOT/f).read_text())
  helpcheck=subprocess.run([sys.executable,str(ROOT/WORKER),'--help'],capture_output=True,text=True,timeout=10,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'})
  assert helpcheck.returncode==0 and '--root-review-sha'in helpcheck.stdout
  assert not(ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_REMAINING3_VIEW_PROFILE_20261105_v64').exists()
  # Pure exponent/QQ-majorant implementation fixtures, never beta/source points.
  e=Engine(metadata_base=25000);SA={pack(x)for x in[(0,0,0,0,0),(1,0,0,0,0),(0,1,0,0,0)]};SB={pack(x)for x in[(0,0,0,0,0),(0,0,1,0,0)]};SC={pack((0,0,0,0,0))}
  e.save('A',SA,Q(7,11),17,{'op':'formal_fixture'});e.save('B',SB,Q(-13,17),19,{'op':'formal_fixture'});e.save('C',SC,Q(19,23),23,{'op':'formal_fixture'})
  e.copy('view','A',Q(3),(1,0,0,0,-1));assert e.live['view']['S']is e.live['A']['S']and e.live['view']['off']==rawindex((1,0,0,0,-1))
  expected=set(SC)|{x+e.live['view']['off']+y-64*sum(256**j for j in range(5))for x in SA for y in SB}
  e.fused_union_product('F','C','view','B','literal_dot')
  assert e.live['F']['S']==expected and e.live['F']['off']==0
  assert e.nodes['literal_dot']['support_profile'].startswith('UNMATERIALIZED')
  op=e.ops[-1];qnew=e.live['F']['q'];assert qnew*int(op['integer_weight_old'])==Q(19,23)and qnew*int(op['integer_weight_dot'])==Q(7,11)*3*Q(-13,17)
  h=e.boundary_hash('F');manual=hashlib.sha256(enc([list(x)for x in sorted(map(unpack,expected))])).hexdigest();assert h==manual
  vshot=e.failure_view_snapshot(out/'formal-view-storage-roundtrip.q5view')
  assert vshot['maps']==4 and vshot['backings']==3 and vshot['all_backing_keys_header_and_view_offsets_stream_roundtrip']
  e.drop('A','view','B')
  flat=e.snapshot(out/'formal-boundary-roundtrip.q5supports');assert flat['maps']==1 and flat['all_keys_exact_stream_roundtrip']
  e.save('empty',set(),Q(0),0,{'op':'empty_fixture'});e.copy('empty_view','empty',Q(3),(3,0,0,0,0));assert e.nodes['empty_view']['exponent_box5']==[[0,0]]*5
  # Tight scalar facts are already independently readable from the exact
  #93-term SOURCE preflight, not guessed from the old guard.
  boundary_cap=12+1024*39+256+5*125000;failure_cap=12+1024*64+4096+5*125000
  assert boundary_cap==665204 and failure_cap==694644
  term_cap=49152;phase_cap=65536
  generated=93*term_cap+3*phase_cap+3*boundary_cap+failure_cap+2*262144+65536
  assert generated==8047824 and generated<8*1024*1024
  inputbytes=sum((ROOT/p).stat().st_size for p in inputs);gate_reserve=2*65536+vshot['bytes']+flat['bytes'];root_reserve=65536
  fresh=inputbytes+gate_reserve+root_reserve+generated
  assert fresh<32*1024*1024
  budget={'id':ID,'source_only':True,'worker':WORKER,'worker_sha256':inputs[WORKER],'helper':HELPER,'helper_sha256':inputs[HELPER],'resume_snapshot':str(shot.relative_to(ROOT)),'resume_snapshot_sha256':sha(shot),'resume_phase07':True,'no_old8_arithmetic_replay':True,'remaining_components':[8,9,10],'remaining_terms':93,
   'SOURCE_metadata_base_upper':50000,'same_SOURCE_owned_limit':150000,'failure_numeric_reserve_ALWAYS_charged':4096,'failure_header_index_numeric_upper':2048,'stored_M_bit_guard':512,'stored_QQ_num_den_bit_guard':128,'positive_M_proven_upper_bits':404,'QQ_proven_upper_bits':80,'H_BOX_width_proven_upper':26,'SOURCE_temporary_integer_bit_upper':16384,'no_VALUE_bit_policy_change':True,
   'snapshot_boundary_maps_upper':39,'snapshot_unique_backing_key_upper':125000,'one_boundary_snapshot_RAW_upper':boundary_cap,'failure_view_snapshot_RAW_upper':failure_cap,'term_receipt_RAW_upper':term_cap,'phase_receipt_RAW_upper':phase_cap,'generated_RAW_upper':generated,'actual_READ_files_RAW':inputbytes,'gate_receipt_RAW_reserve':gate_reserve,'ROOT_receipt_RAW_reserve':root_reserve,'fresh_RAW_upper':fresh,'seconds_hard':60,'math_alarm':55,'checkpoint_deadline':52,'AS_bytes':2000000000,
   'ownership_proof':['A copy creates only immutable view metadata and an exact packed offset; it never copies or mutates backing support.','Every multiply reads views with offsets and creates one new true packed set after a pair-BOX proof.','Fused SOURCE product-union mutates only a uniquely-owned offset0 accumulator; it never allocates the dot support.','Unique backing sets plus4numeric view fields,4index fields per backing,all metadata andscratch are charged.4096 numeric slots are ALWAYS reserved for failure serialization.','The literal dot is UNMATERIALIZED; its card/Hcounts are explicitly safe pair-count uppers, not exact. Accumulator card/Hcounts remain exact positive-support statistics.','All metadata/ops are streamed after EVERY term. Only currently live node profiles and small receipt index remain resident.','A boundary has exactly38function mapstorages at most plusone accumulator, all offsets0 and all backingsdistinct; snapshot output is byte/value lossless.'],
   'RAW_proof':['Each93term has<=7function occurrences: copy+6mul+shiftview+2fusedprovenance+<=7newleaf profiles=<=17nodes.','Each node/ops record<=2048bytes using QQ128,M512,Hwidth26,paircount<=125000*6,80charids.17records+refs fit49152bytes. Runtime capchecked BEFOREwrite.','Each boundary snapshot keys use5bytes; header<=1024bytes/map. Failure view snapshot stores shared backings ONCE plusoffset metadata, not logical duplicatedsupports.','All oldREAD and all failure histories are in the exactinput ledger; allnew receipts,snapshots,partialfailure andreport are simultaneously charged. No gzip or retiredhistory deletion.'],
   'no_coefficient_VALUE_GO':True,'all167x14_87plus66_RHS_oneP290_F9_F10_complements':True,'formal_fixture_snapshot_records':[vshot,flat]}
  bp=out/'resume-views-fused-budget.json';bp.write_bytes(enc(budget));assert bp.stat().st_size<=65536
  report={'id':ID,'all_passed':True,'source_only':True,'profile_not_executed':True,'resume_last_atomic_snapshot_only':True,'view_and_fused_union_fixtures_passed':True,'all_snapshot_fixtures_roundtrip':True,'remaining_terms':93,'generated_RAW_upper':generated,'fresh_RAW_upper':fresh,'same_owned_limit':150000,'stored_bits_M_QQ':[512,128],'SOURCE_temp16384_unchanged':True,'prepared_worker_help_passed':True,'input_sha256':inputs,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start}
  rp=out/'report.json';rp.write_bytes(enc(report));assert rp.stat().st_size<=65536
  signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(rp))
 except Exception as exc:
  (out/'failure.json').write_bytes(enc({'all_passed':False,'source_only':True,'error':repr(exc),'no_profile_or_VALUE_executed':True,'no_automatic_retry':True}));raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',action='store_true');main(p.parse_args())

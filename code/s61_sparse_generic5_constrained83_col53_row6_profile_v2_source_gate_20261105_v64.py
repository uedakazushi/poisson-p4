#!/usr/bin/env python3
"""SOURCE/code/RAW gate for ONE boundary-hash and full-live-snapshot profile."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,ast,hashlib,json,os,resource,signal,subprocess,sys,time
from s61_sparse_generic5_source_packed5_profile_v2_20261105_v64 import Engine,pack,unpack,count_numeric
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_ROW6_PROFILE_V2_SOURCE_GATE_20261105_v64'
WORKER='code/s61_sparse_generic5_constrained83_col53_row6_coefficient_profile_v2_20261105_v64.py'
HELPER='code/s61_sparse_generic5_source_packed5_profile_v2_20261105_v64.py'
TIME='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_ROW6_COEFFICIENT_PROFILE_20261105_v64'
PACKSTOP='S61_SPARSE_GENERIC5_CONSTRAINED83_TENTH_COL53_STAGED_SOURCE_20261105_v64'
OLD='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'
PLAN='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_COEFFICIENT_SIDE_SOURCE_20261105_v64'
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def main(a):
 assert a.freeze
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 signal.signal(signal.SIGALRM,lambda s,f:(_ for _ in()).throw(TimeoutError('gate60s')));signal.alarm(60)
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={}
 def add(p):
  p=Path(p);inputs[str(p.relative_to(ROOT))]=sha(p);return p
 def read(p):return json.loads(add(p).read_text())
 try:
  previous=ROOT/'reports'/TIME
  failure=read(previous/'failure.json');assert failure['forecast']=='hard60s'and len(failure['completed_phase_files'])==8
  partial=read(previous/'partial-profile.json')
  for p,h in partial['input_sha256'].items():assert sha(add(ROOT/p))==h
  # ALL source artifacts from BOTH failed methods remain byte-bound, not
  # deleted or relabeled successful. No old source job is executed.
  for stem in[TIME,PACKSTOP]:
   for p in(ROOT/'reports'/stem).rglob('*'):
    if p.is_file():add(p)
  for f in['code/s61_sparse_generic5_constrained83_col53_row6_coefficient_profile_20261105_v64.py','code/s61_sparse_generic5_source_packed5_profile_20261105_v64.py','code/s61_sparse_generic5_constrained83_tenth_col53_staged_source_20261105_v64.py']:add(ROOT/f)
  old=read(ROOT/'reports'/OLD/'staged-literal-col56-program-and-positive-profiles.json')
  for op in old['literal_ops']:
   if op['op']=='input':add(ROOT/'reports'/OLD/op['binding']['scope'])
  plan=read(ROOT/'reports'/PLAN/'exact-11-row-coefficient-operator-and-callback-owners.json')
  token={x for coefficient in plan['raw_coefficient_side11']for factors,q in coefficient for x,n in factors if x not in['H','s','q','t']}
  maps=len(token)+1;assert maps==39
  for f in[WORKER,HELPER,'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py',str(Path(__file__).relative_to(ROOT)),f'research/notes/{ID}.md',f'research/notes/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_ROW6_COEFFICIENT_PROFILE_V2_20261105_v64.md']:add(ROOT/f)
  for f in[WORKER,HELPER]:ast.parse((ROOT/f).read_text())
  assert not(ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_ROW6_COEFFICIENT_PROFILE_V2_20261105_v64').exists()
  helpcheck=subprocess.run([sys.executable,str(ROOT/WORKER),'--help'],capture_output=True,text=True,timeout=10,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'})
  assert helpcheck.returncode==0 and '--root-review-sha'in helpcheck.stdout
  # Tiny exact formal exponent/metadata fixtures, not beta/source samples.
  S={pack(e)for e in[(-64,-64,-64,-64,-64),(191,191,191,191,191),(0,0,0,0,0),(-1,2,-3,4,-5),(3,0,2,0,1),(3,1,2,0,1)]}
  assert all(pack(unpack(k))==k for k in S)
  e=Engine(metadata_base=25000)
  e.save('fixture',S,Q(7,11),17,{'op':'formal_exponent_fixture'})
  actual=e.boundary_hash('fixture')
  expected=hashlib.sha256(enc([list(x)for x in sorted(map(unpack,S))])).hexdigest();assert actual==expected
  e.save('max4096_header',set(S),Q((1<<4095)-1,(1<<4094)+1),(1<<4095)-1,{'op':'SOURCE_metadata_height_fixture'})
  fixture=e.snapshot(out/'formal-exponent-and-header-roundtrip.q5supports')
  assert fixture['all_keys_exact_stream_roundtrip']and fixture['bytes']<10000
  # Source fields up to4096bits can temporarily multiply to<=16384bits:
  # M-product<=8209; LCM/QQ intermediate<=12288; sum<=12289.
  maxkeys=150000-25000
  header_per_map=4096;headerfixed=256
  one_snapshot=12+header_per_map*maps+headerfixed+5*maxkeys
  assert one_snapshot==785012
  phasecap=262144
  generated=11*one_snapshot+12*phasecap+262144
  assert generated==12043004
  # Include every new gate artifact in the futureREAD set. Their combined
  # record sizes are separately reserved (report/proof max64KiB each).
  inputbytes=sum((ROOT/p).stat().st_size for p in inputs)
  gate_receipts=2*65536+fixture['bytes'];root_receipt=65536
  fresh=inputbytes+gate_receipts+root_receipt+generated
  assert inputbytes<18*1024*1024 and fresh<32*1024*1024
  budget={'id':ID,'source_only':True,'worker':WORKER,'worker_sha256':inputs[WORKER],'helper':HELPER,'helper_sha256':inputs[HELPER],'actual_input_RAW_bytes':inputbytes+gate_receipts,'future_ROOT_review_RAW_reserve':root_receipt,'max_boundary_live_maps':maps,'function_tokens38':sorted(token),'max_total_snapshot_keys':maxkeys,'stored_q_M_bit_limit':4096,'SOURCE_integer_temporary_bit_upper':16384,'one_snapshot_header_RAW_upper':4096*maps+256,'one_snapshot_RAW_upper':one_snapshot,'phase_JSON_RAW_upper':phasecap,'all_generated_RAW_upper':generated,'fresh_RAW_upper':fresh,'SOURCE_numeric_limit':150000,'ROOT_review_metadata_limit':2000,'gate_metadata_limit':4000,'seconds_hard':60,'source_math_alarm_seconds':55,'phase_deadline_seconds':52,'AS_bytes':2000000000,
   'source_change':['Recover the previous8 boundaries ONCE under the changed hash policy; compare their canonical hashes EXACTLY.','Check frozen old callback cardinality/BOX/QQ/M fields at every source node and canonical hashes at BOTH target callback boundaries.','New intermediate support stores exact card/BOX/Hcounts with literal DAG; canonical tuplelex hash only at each11raw component boundary.','One pass obtains BOX/Hcounts; no sorting or full alias at intermediate nodes.','Stream ALL live supports at each boundary as exact5byte packed integers, then readback entry-by-entry in unchanged set order.','Snapshot READ for future resume must use the exact recorded key count/fresh budget; no claim every theoretical125k-key snapshot resumes under150k.'],
   'snapshot_losslessness_proof':['At boundary live map names equal exactly common function-registry values plus ONE accumulator;38tokens+1=39maps.','Every support key satisfies0<=k<2^40 by checked signed5-coordinate BOX, so5bytes is exact.','Writer and reader traverse the SAME unchanged sets in the SAME order, compare EVERY5byte key, header and EOF. No copied set, private coefficient array or full alias.','Header stores id,normalized QQ scaling,positive M and exact term count for every live map. No support or source function is dropped.','Source q/M fields have<=4096bits; each map header<=4096bytes including80charid; all map headers<=39*4096+256.','At every source allocation the SAME150k guard includes base>=25k, metadata and scratch, hence live support keys<=125000.'],
   'no_coefficient_VALUE_GO':True,'all167x14_87plus66_RHS_oneP290_F9_F10_complements':True,'formal_fixture':fixture,'gate_owns_peak':e.peak}
  bp=out/'whole-SOURCE-and-snapshot-budget.json';bp.write_bytes(enc(budget));assert bp.stat().st_size<=65536
  report={'id':ID,'all_passed':True,'source_only':True,'no_actual_support_profile_or_VALUE_executed':True,'prepared_changed_SOURCE_worker_help_passed':True,'old_TWO_stops_and_all_HISTORY_READ_bytes_preserved':True,'formal_snapshot_fixture_roundtrip':True,'max_boundary_maps':39,'one_snapshot_RAW_upper':one_snapshot,'all_generated_RAW_upper':generated,'fresh_RAW_upper':fresh,'same_owned_limit':150000,'SOURCE_temporary_integer_bits_upper':16384,'input_sha256':inputs,'evidence_sha256':{bp.name:sha(bp),'formal-exponent-and-header-roundtrip.q5supports':sha(out/'formal-exponent-and-header-roundtrip.q5supports')},'seconds':time.monotonic()-start}
  rp=out/'report.json';rp.write_bytes(enc(report));assert rp.stat().st_size<=65536
  assert count_numeric(report)+count_numeric(budget)<=4000
  signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(rp))
 except Exception as exc:
  (out/'failure.json').write_bytes(enc({'all_passed':False,'source_only':True,'error':repr(exc),'no_actual_SOURCE_profile_or_VALUE_run':True,'no_automatic_retry':True}));raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',action='store_true');main(p.parse_args())

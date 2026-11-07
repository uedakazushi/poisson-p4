"""SOURCE-only budget proposal: immutable low bands plus new high-band codecs."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time
from s61_sparse_generic5_source_view_profile_20261105_v64 import count_numeric
ROOT=Path(__file__).resolve().parents[1]
ID='S61_COL53_NEXT17_79_FULL_HL_EXTENSION_BUDGET_SOURCE_20261105_v64'
FINAL='reports/S61_COL53_NEXT17_79_FINALIZED_SOURCE_20261105_v64'
FULL='reports/S61_COL53_NEXT17_79_FULL_HL_SUPPORT_20261105_v64'
LOW='reports/S61_COL53_NEXT17_79_SOURCE_PROJECTION_20261105_v64'
OLD='reports/S61_COL53_NEXT17_79_STREAMED_VALUE_20261105_v64'
PEER='reports/S61_COL53_NEXT17_79_STREAMED_VALUE_20261105_v64_LOCAL_IDENTITY_REVIEW/report.json'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def enc(d):return(json.dumps(d,separators=(',',':'))+'\n').encode()
def ui(n):return max(1,(n.bit_length()+6)//7)
def si(bits):n=(bits+7)//8;return ui(2*n+1)+n
def main():
 start=time.monotonic();out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={}
 reports={FINAL+'/report.json':'548ebf5f4078707d6ee2bb6dd1867c7161c2d32bd553ebef4ce454600c7db877',FULL+'/report.json':'4488f4a9425da0c0c81b93b5546411f9f03edceaaf7767be19ed4a6277350ace',LOW+'/report.json':'8cba44d77933886756700a63fb250058256db3310272b6ca626936116c785d60',OLD+'/report.json':'4b1ed0275771c6f13a8f250b14a71a125277d3b8e6f17fd5995142b8d4dc0c54',PEER:'2b0b2212d91ad5f74481ffc5428a9c33bd1d6cfb018a78c7a3cea7130747dd88'}
 for f,h in reports.items():assert sha(ROOT/f)==h;pins[f]=h
 assert json.loads((ROOT/PEER).read_text())['all_passed']
 def evidence(folder,name):
  r=json.loads((ROOT/folder/'report.json').read_text());f=folder+'/'+name;pins[f]=sha(ROOT/f);assert pins[f]==r['evidence_sha256'][name];return ROOT/f
 D=json.loads(evidence(FINAL,'compiled-full-SLP-and-fixed-projection-plan.json').read_text());Dnumeric=count_numeric(D);assert Dnumeric==65701
 F=dict(json.loads(evidence(FULL,'all1884-exact-complete-face-source-counts.json').read_text()));L=dict(json.loads(evidence(LOW,'all-projected-SOURCE-node-counts.json').read_text()));A={}
 index=evidence(OLD,'node-values-index.jsonl');spool=evidence(OLD,'node-values.q5stream')
 with index.open('rb')as f:
  for line in f:r=json.loads(line);A[r['node']]=r['actual_nonzero_slots']
 assert len(F)==len(L)==len(A)==1884 and F.keys()==L.keys()==A.keys()
 H={k:F[k]-L[k]for k in F};V={k:H[k]+A[k]for k in F};assert min(H.values())>=0 and all(A[k]<=L[k]for k in F)
 normal={};prefix={};records=[];onegen=0;weights=0;Sourcepeak=0;record_numeric=0
 for k in sorted(F):
  n=D['nodes'][k];q=Q(n['q']);M=int(n['M']);args=n['args'];child=list(dict.fromkeys(args if n['op']in['copy','multiply']else[a for a in args if D['nodes'][a]['q']!='0'and D['nodes'][a]['lead5'][:2]==n['lead5'][:2]]))
  ct=sum(V[c]for c in child);hh=H[k];ww=[];bound=0
  if not q:assert M>=0
  elif n['op']=='input':bound=M
  elif n['op']=='copy':
   c=args[0];w=Q(D['nodes'][c]['q'])*Q(n['scalar'])/q;assert w.denominator==1;ww=[w.numerator];bound=abs(w.numerator)*int(D['nodes'][c]['M'])
  elif n['op']=='multiply':
   a,b=args;w=Q(D['nodes'][a]['q'])*Q(D['nodes'][b]['q'])*Q(n['scalar'])/q;assert w.denominator==1;ww=[w.numerator];bound=min(F[a],F[b])*int(D['nodes'][a]['M'])*int(D['nodes'][b]['M'])*abs(w.numerator)
  else:
   assert n['op']=='add'
   for j,c in enumerate(args):
    if c not in child:continue
    w=Q(D['nodes'][c]['q'])*(Q(n['scalar'])if j else Q(1))/q;assert w.denominator==1;ww.append(w.numerator);bound+=abs(w.numerator)*int(D['nodes'][c]['M'])
  assert bound<=M,(k,'full-face positive partial majorant');weights+=len(ww)
  phases={'higher_band_arithmetic':118000+7*(ct+hh)+4096,'two_band_child_decode':118000+7*ct+8*max([max(A[c],H[c])for c in child]+[0])+2048+4096,'codec_after_child_release':118000+10*hh+4096}
  exceptional={}
  if n['op']=='input'and q:
   phases['input_decode']=118000+16152+7*hh+4096;exceptional['leaf_STOP_codec']=118000+16152+10*hh+4096
  if k==D['FULL_F11_original_root']:
   phases['anchor_after_child_release']=118000+7*hh+8*411+2048+4096;exceptional['anchor_STOP_codec']=118000+10*hh+8*411+2048+4096
  for name,v in phases.items():
   if v>normal.get(name,{'slots':0})['slots']:normal[name]={'slots':v,'node':k}
  for name,v in exceptional.items():
   if v>prefix.get(name,{'slots':0})['slots']:prefix[name]={'slots':v,'node':k}
  assert max([*phases.values(),*exceptional.values()])<150000
  # Every SOURCE support coordinate was checked in[-64,191]. All5 signed
  # shift varints and all5 relative exponents therefore need at most2 bytes.
  cell=9 if not hh else 9+si(abs(q.numerator).bit_length()+M.bit_length())+si(q.denominator.bit_length())+10+ui(hh)+hh*(si(M.bit_length())+10)
  onegen+=cell
  records.append({'node':k,'old_positive_slots':L[k],'old_actual_slots':A[k],'full_positive_slots':F[k],'new_positive_slots':hh,'full_actual_slots_upper':V[k],'new_codec_RAW_upper':cell,'integer_weights':list(map(str,ww)),'all_partial_source_bounds_passed':True})
  # This SOURCE budget preparation owns scalar metadata only, never maps.
  record_numeric+=count_numeric(records[-1]);Sourcepeak=max(Sourcepeak,25000+Dnumeric+10*1884+record_numeric+4096);assert Sourcepeak<150000
 generated=2*onegen+2048*1884+524288+1048576
 # The future program may reuse the exact prior physical read closure. New
 # code/gate/plan documents are an explicit reserve, not an implemented worker.
 oldreport=json.loads((ROOT/OLD/'report.json').read_text())
 future_READ={**oldreport['input_sha256'],**pins}
 future_READ['code/s61_col53_next17_79_full_HL_extension_budget_source_20261105_v64.py']=sha(__file__)
 for f,h in future_READ.items():assert sha(ROOT/f)==h
 reserved_new_code_plan_gate=786432;reserved_new_peer_ROOT=131072;control=262144
 raw=sum((ROOT/f).stat().st_size for f in future_READ)+generated+reserved_new_code_plan_gate+reserved_new_peer_ROOT+control
 assert raw<32*1024*1024
 proof={'SOURCE_BUDGET_PROPOSAL_only_no_extension_worker_implemented':True,'all1884_source_weight_and_full_Cauchy_checks_passed':True,'unchanged_global_M_max_bits':max(int(n['M']).bit_length()for n in D['nodes']),'new_positive_slots_sum':sum(H.values()),'old_actual_slots_sum':sum(A.values()),'root_new_positive_slots':H[D['minor_B_root']],'normal_phase_maxima':normal,'guard_prefix_phase_maxima':prefix,'META118000_representation_obligation':'replace prior2N source-count dictionary by2N extra paired-index numeric fields; stream old manifest/current per-node plan, no all-plan parse; paired index6N numeric, same explicit baseline111753 plus original6247 reserve','new_codec_one_generation_RAW_upper':onegen,'TWO_new_highband_generations_RAW_upper':2*onegen,'ALL_generated_RAW_upper':generated,'future_fresh_RAW_upper':raw,'future_fixed_READ_sha256':future_READ,'new_code_plan_gate_reserved_BYTES':reserved_new_code_plan_gate,'new_peer_ROOT_reserved_BYTES':reserved_new_peer_ROOT,'control_reserved_BYTES':control,'all_node_records':records,'old_prefix_immutable_and_not_recomputed':True,'mixed_oldlow_times_oldlow_new_output_terms_retained':True,'full5coords_VALUE7_unchanged':True,'runtime_worker_and_whole_lifetime_audit_still_required':True}
 (out/'source-cauchy-two-band-lifetime-and-raw-proposal.json').write_bytes(enc(proof))
 pins['code/s61_col53_next17_79_full_HL_extension_budget_source_20261105_v64.py']=sha(__file__)
 report={'id':ID,'all_passed':True,'SOURCE_PROPOSAL_ONLY_no_new_signed_values_no_VALUE_GO':True,'complete_leading_HL_face_S_range':[7,38],'same1884_nodes':True,'new_high_positive_slots':sum(H.values()),'root_new_positive_slots':H[D['minor_B_root']],'all1884_full_Cauchy_passed':True,'normal_upper':max(r['slots']for r in normal.values()),'guard_prefix_upper':max(r['slots']for r in prefix.values()),'fresh_RAW_upper':raw,'Source_preparation_owned_upper':Sourcepeak,'future_worker_CODE_not_yet_implemented':True,'input_sha256':pins,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start}
 (out/'report.json').write_bytes(enc(report));print('report_sha256',sha(out/'report.json'));print('normal/prefix/RAW',report['normal_upper'],report['guard_prefix_upper'],raw)
if __name__=='__main__':main()

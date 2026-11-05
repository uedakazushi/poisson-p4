"""Prepared streamed signed needed-projection evaluation. Requires fresh peer+ROOT pins."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
from s61_col53_projection_leaf_io_20261105_v64 import read_leaf,leaf_paths
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
from s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64 import write_cell,validate_cell
from s61_sparse_generic5_source_view_profile_20261105_v64 import count_numeric
ID='S61_COL53_NEXT17_79_STREAMED_VALUE_20261105_v64'
GATE='S61_COL53_NEXT17_79_STREAMED_VALUE_GATE_20261105_v64'
FINAL='reports/S61_COL53_NEXT17_79_FINALIZED_SOURCE_20261105_v64'
SOURCE='reports/S61_COL53_NEXT17_79_SOURCE_PROJECTION_20261105_v64'
METHOD='reports/S61_COL53_NEXT17_79_STREAM_VALUE_METHOD_SOURCE_V2_20261105_v64'
def enc(d):return(json.dumps(d,separators=(',',':'))+'\n').encode()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def load_data():
 p=ROOT/FINAL/'compiled-full-SLP-and-fixed-projection-plan.json';D=json.loads(p.read_text())
 p2=ROOT/SOURCE/'all-projected-SOURCE-node-counts.json';C=dict(json.loads(p2.read_text()))
 return D,C
def recover_committed(out):
 # Only complete manifest records whose exact stream bytes hash correctly count.
 result={'nodes_committed':0,'last_node':None,'stream_prefix_bytes':0,'manifest_prefix_bytes':0}
 path=out/'node-values-index.jsonl'
 if not path.exists():return result
 with path.open('rb')as f,(out/'node-values.q5stream').open('rb')as s:
  while True:
   line=f.readline()
   if not line:break
   try:r=json.loads(line)
   except json.JSONDecodeError:break
   s.seek(r['stream_offset']);body=s.read(r['bytes'])
   if len(body)!=r['bytes']or hashlib.sha256(body).hexdigest()!=r['payload_sha256']:break
   assert r['stream_offset']==result['stream_prefix_bytes']
   result={'nodes_committed':result['nodes_committed']+1,'last_node':r['node'],'stream_prefix_bytes':r['stream_offset']+r['bytes'],'manifest_prefix_bytes':f.tell()}
 return result
def prepare():
 out=ROOT/'reports'/GATE;out.mkdir(exist_ok=False);pins={}
 expected={FINAL+'/report.json':'548ebf5f4078707d6ee2bb6dd1867c7161c2d32bd553ebef4ce454600c7db877',SOURCE+'/report.json':'8cba44d77933886756700a63fb250058256db3310272b6ca626936116c785d60',METHOD+'/report.json':'e15323af46f399a232e5ad482439b7c1e910e012ca6f85bcb6fa572cb2338c99'}
 for f,h in expected.items():assert sha(ROOT/f)==h;pins[f]=h
 D,C=load_data();assert len(C)==1884 and sum(C.values())==71584 and max(C.values())==129
 for f in[FINAL+'/compiled-full-SLP-and-fixed-projection-plan.json',SOURCE+'/all-projected-SOURCE-node-counts.json',SOURCE+'/boundary-2997.json']:
  pins[f]=sha(ROOT/f)
  owner=FINAL if f.startswith(FINAL+'/')else SOURCE
  receipt=json.loads((ROOT/owner/'report.json').read_text())
  assert pins[f]==receipt['evidence_sha256'][Path(f).name]
 # Actual byte-pinned existing leaf physical data only; no new signed query.
 maxdecode=0
 for k in C:
  n=D['nodes'][k]
  if n['op']=='input'and n['q']!='0'and C[k]:
   for f in leaf_paths(n):pins[f]=sha(ROOT/f)
   p,shift,d=read_leaf(n);maxdecode=max(maxdecode,d+2048)
 assert maxdecode<=16152
 for f in['code/s61_col53_next17_79_streamed_value_20261105_v64.py','code/s61_col53_projection_leaf_io_20261105_v64.py','code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py','code/s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64.py','code/s61_sparse_generic5_source_view_profile_20261105_v64.py','code/s61_sparse_generic5_source_packed5_profile_v2_20261105_v64.py']:
  pins[f]=sha(ROOT/f)
 anchor='reports/S61_COL53_H_FACES_VALUES_20261105_v64/H-minus3-face.bin';pins[anchor]=sha(ROOT/anchor)
 assert pins[anchor]=='cc9165dc8ec941004451f7924944e839c36db3e0b4e939be17b4a8f9afd88bce'
 method=json.loads((ROOT/METHOD/'report.json').read_text());assert method['whole_VALUE_forecast']==138535 and method['generated_RAW_forecast']==23968000
 # Method proof is hash-read, never retained as a second table in the VALUE job.
 proof=METHOD+'/streamed-residency-height-codec-storage-source-proof.json';pins[proof]=sha(ROOT/proof)
 assert pins[proof]==method['evidence_sha256'][Path(proof).name]
 fixture='reports/S61_COL53_NEXT17_79_STREAM_CODEC_FIXTURES_20261105_v64/report.json'
 fr=json.loads((ROOT/fixture).read_text());assert fr['all_passed']and fr['synthetic_only_no_source_value_queries']
 pins[fixture]=sha(ROOT/fixture)
 for f,h in fr['input_sha256'].items():assert sha(ROOT/f)==h;pins[f]=h
 raw=sum((ROOT/f).stat().st_size for f in pins)+23968000+512*1024+2*65536+262144
 assert raw<32*1024*1024
 B={'worker_sha256':sha(__file__),'metadata_base':118000,'whole_numeric_forecast':138535,'global_numeric_guard':150000,'maxdecode':16152,'commit_max_bits':1024,'actual_source_M_max_bits':824,'temporary_guard_bits':1921,'hard_VALUE_bits':2048,'codec_naive_q_di_G_temp_bits':1676,'one_generation_RAW_upper':9268352,'ALL_generated_RAW_upper':23968000,'manifest_line_RAW_max':2048,'fresh_RAW_upper':raw,'math_s':50,'signal_s':55,'hard_s':60,'AS_bytes':2000000000,'Root_or_peer_document_bytes':65536,'Root_or_peer_document_numeric':2000,'node_file_two_generations_preserved':True,'fixed_q_NOT_primitive':True,'zero_math_G':0,'existing_codec_zero_G_sentinel':1}
 (out/'streamed-whole-lifetime-and-raw-budget.json').write_bytes(enc(B))
 d={'id':GATE,'all_passed':True,'SOURCE_PREPARATION_only_no_new_signed_values':True,'whole_numeric_forecast':138535,'fresh_RAW_upper':raw,'two_codec_generations_all_counted':True,'zero_content_and_G_conventions_explicit':True,'nonempty_qStored_div_qNode_positive_integer_even_if_qNode_negative':True,'FULL_F11_and_all5coefficientkeys_kept':True,'no_matrix_RHS_update':True,'WHOLE_CODE_PEER_and_ROOT_GO_REQUIRED':True,'input_sha256':pins,'evidence_sha256':{'streamed-whole-lifetime-and-raw-budget.json':sha(out/'streamed-whole-lifetime-and-raw-budget.json')}}
 (out/'report.json').write_bytes(enc(d));print('report_sha256',sha(out/'report.json'));print('whole/RAW',138535,raw)
def value(a):
 assert all([a.gate_file,a.gate_sha,a.source_peer_file,a.source_peer_sha,a.root_review_file,a.root_review_sha])
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 tmp=ROOT/'work'/(ID+'-canonical-cells');tmp.mkdir(exist_ok=False)
 B={'ALL_generated_RAW_upper':23968000};META=118000;pins={};index={};held=[];work={};transient=0;peak=0;context={};cursor={'nodes_committed':0,'last_node':None,'stream_prefix_bytes':0,'manifest_prefix_bytes':0};tempbytes=0;stream=None;manifest=None;root_records={};updates=0;D=None
 def put(n,d):
  b=enc(d);assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())+sum(p.stat().st_size for p in tmp.iterdir()if p.is_file())+len(b)<=B['ALL_generated_RAW_upper'];(out/n).write_bytes(b)
 def alarm(s,f):raise TimeoutError('signedVALUE55s; preserve all completed canonical segments')
 signal.signal(signal.SIGALRM,alarm);signal.alarm(55)
 def own(extra=4096,new=False):
  nonlocal peak
  maps=[work]+held;seen=set();terms=0
  for m in maps:
   if id(m)not in seen:seen.add(id(m));terms+=len(m)
  n=META+7*(terms+int(new))+transient+extra;peak=max(peak,n);assert n<=150000,('prospective_stream_VALUE_owned',n,context)
 def raw(v):assert abs(v).bit_length()<=1921;return v
 def update(e,v,M):
  nonlocal updates
  assert len(e)==5 and all(isinstance(x,int)for x in e)and abs(v)<=M and abs(v).bit_length()<=1024
  if v:
   if e not in work:own(new=True)
   work[e]=v
  else:work.pop(e,None)
  updates+=1
 def child(k):
  nonlocal transient
  rec=index[k]
  path=tmp/('node-'+str(k)+'.bin');body=path.read_bytes();assert len(body)==rec[1]and hashlib.sha256(body).hexdigest()==rec[3]
  if rec[2]==0:
   assert body==b'Q5CELL1\n\0';return{}
  qi=Q(D['nodes'][k]['q']);assert qi
  transient=8*C[k]+2048;own();assert body[:9]==b'Q5CELL1\n\1';p=decode(body[9:])
  ratio=Q(p['content'])/qi;assert ratio.denominator==1 and ratio.numerator>0;G=ratio.numerator
  assert G.bit_length()<=824;vals={};held.append(vals);M=int(D['nodes'][k]['M'])
  for t in p['terms']:
   e=tuple(p['shift'][j]+int(t[j+1])for j in range(5));c=raw(int(t[0])*G);assert c and abs(c)<=M and e not in vals
   # held contains the map under construction, count its prospective new key.
   own(new=True);vals[e]=c
  assert len(vals)==rec[2]<=C[k];del p,body;transient=0;own();return vals
 def anchor_check(vals,qi,lead,degree):
  nonlocal transient
  transient=8*411+2048;own();body=(ROOT/'reports/S61_COL53_H_FACES_VALUES_20261105_v64/H-minus3-face.bin').read_bytes();p=decode(body[9:]);G=Q(p['content'])/qi;assert G.denominator==1 and G>0;expected={}
  for t in p['terms']:
   e=tuple(p['shift'][j]+int(t[j+1])for j in range(5))
   if e[0]==lead[0]and e[1]==lead[1]and 0<=e[2]-lead[2]<=degree:expected[e]=int(t[0])*G.numerator
  assert vals==expected;del expected,p,body;transient=0;own()
 try:
  gp=ROOT/a.gate_file;assert sha(gp)==a.gate_sha;g=json.loads(gp.read_text());assert g['all_passed'];pins[str(gp.relative_to(ROOT))]=sha(gp)
  bp=gp.parent/'streamed-whole-lifetime-and-raw-budget.json';assert sha(bp)==g['evidence_sha256'][bp.name];B=json.loads(bp.read_text());pins[str(bp.relative_to(ROOT))]=sha(bp);assert sha(__file__)==B['worker_sha256']
  for filename,want in[(a.source_peer_file,a.source_peer_sha),(a.root_review_file,a.root_review_sha)]:
   p=ROOT/filename;assert sha(p)==want and p.stat().st_size<=65536;r=json.loads(p.read_text());assert r['all_passed']and count_numeric(r)<=2000
   assert r['input_sha256'][a.gate_file]==a.gate_sha and r['input_sha256']['code/s61_col53_next17_79_streamed_value_20261105_v64.py']==B['worker_sha256']
   pins[filename]=want
  for f,h in g['input_sha256'].items():assert sha(ROOT/f)==h;pins[f]=h
  D,C=load_data();nodes=D['nodes'];order=sorted(C);degree=D['needed_node_degree'];active=D['active_tied_HL_children'];root=D['minor_B_root'];assert len(order)==1884 and C[root]==129
  stream=(out/'node-values.q5stream').open('w+b');manifest=(out/'node-values-index.jsonl').open('w+b')
  root_support={tuple(e)for e in json.loads((ROOT/SOURCE/'boundary-2997.json').read_text())};assert len(root_support)==129
  for pos,k in enumerate(order):
   n=nodes[k];qi=Q(n['q']);M=int(n['M']);lead=n['lead5'];cap=degree[str(k)];context={'node':k,'order_position':pos,'phase':'start_node','fixed_SOURCE_lead':lead};memo={};held=[];work={};own()
   if not qi or not C[k]:pass
   elif n['op']=='input':
    transient=16152;own();p,shift,decoded=read_leaf(n);assert decoded+2048<=16152
    for ti,t in enumerate(p['terms']):
     context['term_cursor']=ti;e=tuple(p['shift'][j]+int(t[j+1])+shift[j]for j in range(5))
     if e[0]==lead[0]and e[1]==lead[1]and 0<=e[2]-lead[2]<=cap:update(e,raw(int(t[0])),M)
    del p;transient=0
   else:
    args=n['args'];wanted=args if n['op']in['copy','multiply']else active[str(k)]
    for c in dict.fromkeys(wanted):memo[c]=child(c)
    if n['op']=='copy':
     c=args[0];w=Q(nodes[c]['q'])*Q(n['scalar'])/qi;assert w.denominator==1
     for ti,(e,v)in enumerate(memo[c].items()):
      context['term_cursor']=ti;ee=tuple(e[j]+n['shift5'][j]for j in range(5))
      if ee[0]==lead[0]and ee[1]==lead[1]and 0<=ee[2]-lead[2]<=cap:update(ee,raw(v*w.numerator),M)
    elif n['op']=='multiply':
     a0,b0=args;w=Q(nodes[a0]['q'])*Q(nodes[b0]['q'])*Q(n['scalar'])/qi;assert w.denominator==1;la=nodes[a0]['lead5'];lb=nodes[b0]['lead5']
     for ai,(e,v)in enumerate(memo[a0].items()):
      for bi,(f,z)in enumerate(memo[b0].items()):
       context['pair_cursor']=[ai,bi]
       if(e[2]-la[2])+(f[2]-lb[2])>cap:continue
       ee=tuple(e[j]+f[j]+n['shift5'][j]for j in range(5));term=raw(raw(v*z)*w.numerator);update(ee,raw(work.get(ee,0)+term),M)
    else:
     assert n['op']=='add'
     for j,c in enumerate(args):
      if c not in memo:continue
      w=Q(nodes[c]['q'])*(Q(n['scalar'])if j else Q(1))/qi;assert w.denominator==1
      for ti,(e,v)in enumerate(memo[c].items()):
       context['sum_cursor']=[j,ti]
       if e[0]==lead[0]and e[1]==lead[1]and 0<=e[2]-lead[2]<=cap:update(e,raw(work.get(e,0)+raw(v*w.numerator)),M)
   assert len(work)<=C[k]
   if k==D['FULL_F11_original_root']:anchor_check(work,qi,lead,cap)
   if k==root:assert set(work)<=root_support
   context['phase']='canonical_codec_and_stream_commit'
   own(extra=3*len(work)+4096)
   upper=512+116*C[k];assert tempbytes+upper+stream.seek(0,2)+upper+manifest.tell()+2048+1572864<=B['ALL_generated_RAW_upper']
   path=tmp/('node-'+str(k)+'.bin');rec=write_cell(path,qi,work);assert rec['bytes']<=upper and validate_cell(path,qi,work,rec)
   actualG=0 if not work else rec['G']
   if not work:assert rec['content']=='0'and rec['terms']==0 and rec['G']==1
   else:
    stored=Q(rec['content']);ratio=stored/qi;assert ratio.denominator==1 and ratio.numerator==actualG>0
    assert(stored>0)==(qi>0)
   tempbytes+=rec['bytes'];offset=stream.seek(0,2)
   with path.open('rb')as t:
    for body in iter(lambda:t.read(65536),b''):stream.write(body)
   stream.flush();end=stream.tell();assert end-offset==rec['bytes'];stream.seek(offset);body=stream.read(rec['bytes']);assert hashlib.sha256(body).hexdigest()==rec['sha256'];stream.seek(end)
   record={'node':k,'SOURCE_positive_slots':C[k],'actual_nonzero_slots':len(work),'stream_offset':offset,'bytes':rec['bytes'],'payload_sha256':rec['sha256'],'codec':rec,'math_G':actualG,'codec_zero_G_sentinel':1 if not work else None,'fixed_SOURCE_q':str(qi),'fixed_SOURCE_lead_UNCHANGED':lead,'projection_zero_only_NOT_whole_function_zero':not work and bool(qi),'two_byteidentical_generations_preserved':True}
   line=enc(record);assert len(line)<=2048;manifest.write(line);manifest.flush();index[k]=(offset,rec['bytes'],len(work),rec['sha256']);cursor={'nodes_committed':pos+1,'last_node':k,'stream_prefix_bytes':end,'manifest_prefix_bytes':manifest.tell()}
   if k in D['four_callback_roots'].values()or k==root:root_records[str(k)]=record
   if k==root:
    target=tuple(D['target_full5']);targetq=qi*work.get(target,0);witness=None
    if work:
     e=min(work);witness={'exponent5':list(e),'coefficient_QQ':str(qi*work[e])}
    put('actual-needed-minor-projection-summary.json',{'target_exponent5':list(target),'target_coefficient_QQ':str(targetq),'actual_nonzero_slots':len(work),'nonzero_witness':witness,'root_node':root,'fixed_q':str(qi),'full_F11_not_replaced':True,'only_needed_HL_s2_projection':True,'all_matrix_RHS_and_complements_UNCHANGED':True})
   memo.clear();held.clear();work={};assert time.monotonic()-start<50
  stream.close();manifest.close();stream=None;manifest=None
  committed=recover_committed(out);assert committed==cursor and committed['nodes_committed']==1884
  generated=sum(p.stat().st_size for p in out.iterdir()if p.is_file())+sum(p.stat().st_size for p in tmp.iterdir()if p.is_file());assert generated<=B['ALL_generated_RAW_upper']
  rawtotal=sum((ROOT/f).stat().st_size for f in pins)+generated+262144;assert rawtotal<=B['fresh_RAW_upper']
  put('byteidentity-storage-and-cursor-receipt.json',{'cursor':cursor,'all1884_segments_HASH_and_canonical_temp_readback':True,'all_temporary_bytes_preserved_under_work':str(tmp.relative_to(ROOT)),'all_temporary_bytes_reconstructible_from_manifest_stream_segments':True,'ALL_two_generations_RAW_counted':True,'generated_RAW_bytes':generated,'root_cell_records':root_records})
  report={'id':ID,'all_passed':True,'only_needed_fixed_minor_projection_evaluated':True,'nodes_committed':1884,'target_coefficient_QQ':str(targetq),'nonzero_projected_minor_witness':witness,'actual_whole_owned_peak':peak,'source_whole_forecast':138535,'global_guard':150000,'all_terms_full5coords_pluscoeff7':True,'fixed_q_never_assumed_primitive':True,'zero_mathG0_codec_sentinelG1_sign_rules_checked':True,'new_matrix_RHS_update_NOT_performed':True,'fresh_RAW_upper_including_control_reserve':rawtotal,'input_sha256':pins,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start}
  put('report.json',report);signal.alarm(0);print('report_sha256',sha(out/'report.json'));print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}))
 except Exception as exc:
  if stream is not None:stream.flush()
  if manifest is not None:manifest.flush()
  recovered=recover_committed(out)if(out/'node-values.q5stream').exists()else cursor
  if D is not None and work:
   k=context['node'];own(extra=3*len(work)+4096)
   rec=write_cell(out/'partial-current-accumulator.bin',Q(D['nodes'][k]['q']),work);validate_cell(out/'partial-current-accumulator.bin',Q(D['nodes'][k]['q']),work,rec)
   put('partial-current-accumulator-status.json',{'node':k,'context':context,'complete_polynomial_NOT_claimed':True,'resume_only_from_completed_node_cursor_without_new_verified_plan':True,'codec':rec})
  put('failure.json',{'all_passed':False,'error':repr(exc),'context':context,'committed_cursor':recovered,'no_automatic_retry':True,'all_completed_segments_and_tempfiles_preserved':True})
  raise
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--prepare-gate',action='store_true');p.add_argument('--value',action='store_true');p.add_argument('--gate-file');p.add_argument('--gate-sha');p.add_argument('--source-peer-file');p.add_argument('--source-peer-sha');p.add_argument('--root-review-file');p.add_argument('--root-review-sha');a=p.parse_args()
 if a.prepare_gate:prepare()
 else:assert a.value;value(a)

"""Prepared ONE SOURCE-only leadingHL/s2 support profile; requires fresh ROOT pin."""
from pathlib import Path
import argparse,hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
from s61_col53_projection_leaf_io_20261105_v64 import leaf_paths,read_leaf
from s61_sparse_generic5_source_packed5_profile_v2_20261105_v64 import pack,unpack,Engine
from s61_sparse_generic5_source_view_profile_20261105_v64 import count_numeric
ID='S61_COL53_NEXT17_79_SOURCE_PROJECTION_20261105_v64'
GATE='S61_COL53_NEXT17_79_SOURCE_PROJECTION_GATE_20261105_v64'
FINAL='reports/S61_COL53_NEXT17_79_FINALIZED_SOURCE_20261105_v64'
def enc(d):return(json.dumps(d,separators=(',',':'))+'\n').encode()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def prepare():
 out=ROOT/'reports'/GATE;out.mkdir(exist_ok=False);pins={}
 p=ROOT/FINAL/'report.json';assert sha(p)=='548ebf5f4078707d6ee2bb6dd1867c7161c2d32bd553ebef4ce454600c7db877';r=json.loads(p.read_text());assert r['all_passed'];pins[str(p.relative_to(ROOT))]=sha(p)
 p=ROOT/FINAL/'compiled-full-SLP-and-fixed-projection-plan.json';assert sha(p)==r['evidence_sha256'][p.name];D=json.loads(p.read_text());pins[str(p.relative_to(ROOT))]=sha(p)
 nodes=D['nodes'];N=len(D['needed_node_degree']);assert N==1884 and len(nodes)==2998
 meta=25000+count_numeric(D)+8*N+4096;assert meta==109869
 peakdecode=0;maxterms=0;inputrecords=[]
 for kstr in D['needed_node_degree']:
  k=int(kstr);n=nodes[k]
  if n['op']!='input'or n['q']=='0':continue
  for f in leaf_paths(n):pins[f]=sha(ROOT/f)
  p,shift,numeric=read_leaf(n);peakdecode=max(peakdecode,numeric+2048);maxterms=max(maxterms,len(p['terms']))
  assert meta+numeric+2048+4096<=150000
  inputrecords.append({'node':k,'decoded_numeric':numeric,'existing_terms':len(p['terms']),'all_source_QQ_transforms_exact':True})
 assert peakdecode<=16152
 helpers=['code/s61_col53_projection_leaf_io_20261105_v64.py','code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py','code/s61_sparse_generic5_source_view_profile_20261105_v64.py','code/s61_sparse_generic5_source_packed5_profile_v2_20261105_v64.py','code/s61_col53_next17_79_source_projection_20261105_v64.py']
 for f in helpers:pins[f]=sha(ROOT/f)
 # SOURCE sets use one exact40bit encoded key. Every product decodes both
 # operands and checks all5 coords before pack; there is no packed addition.
 # VALUE representations remain5coordinates+coefficient7, not this packing.
 raw=sum((ROOT/f).stat().st_size for f in pins)+4*1024*1024+512*1024+65536+262144;assert raw<32*1024*1024
 b={'worker_sha256':sha(__file__),'metadata_base':meta,'whole_SOURCE_guard':150000,'decoded_leaf_upper':16152,'input_max_terms':maxterms,'fixed_source_packed_key_bits':40,'SOURCE_only_one_owned_int_per_key_full5_exact':True,'future_VALUE_7_count_NOT_changed':True,'generated_RAW_upper':4*1024*1024,'fresh_RAW_upper':raw,'math_s':50,'signal_s':55,'hard_s':60,'AS_bytes':2000000000,'ROOT_doc_bytes':65536,'ROOT_doc_numeric':2000,'no_completion_promise_if_prospective_support_guard_hits':True}
 (out/'SOURCE-support-lifetime-budget.json').write_bytes(enc(b));(out/'exact-input-transform-and-decode-profiles.json').write_bytes(enc(inputrecords))
 d={'id':GATE,'all_passed':True,'SOURCE_support_JOB_NOT_RUN':True,'no_signed_VALUE_GO':True,'fixed_candidate':[17,79],'functional_target':[-5,0,9,-3,19],'metadata_base':meta,'max_leaf_numeric_upper':peakdecode,'fresh_RAW_upper':raw,'full5coords_checked_before_packing':True,'ROOT_review_required_before_source_profile':True,'input_sha256':pins,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()}}
 (out/'report.json').write_bytes(enc(d));print('report_sha256',sha(out/'report.json'));print('meta/decode/RAW',meta,peakdecode,raw)
def execute(a):
 assert a.gate_file and a.gate_sha and a.root_review_file and a.root_review_sha
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);context={};cache={};work=set();counts=[];pins={};transient=0;peak=0;pairs=0;B={'generated_RAW_upper':65536};META=0
 def put(n,d):
  b=enc(d);assert sum(p.stat().st_size for p in out.iterdir())+len(b)<=B['generated_RAW_upper'];(out/n).write_bytes(b)
 def alarm(s,f):raise TimeoutError('SOURCEprojection55s; no retry')
 signal.signal(signal.SIGALRM,alarm);signal.alarm(55)
 def own(extra=2048,new=False):
  nonlocal peak
  value=META+sum(map(len,cache.values()))+len(work)+int(new)+transient+extra;peak=max(peak,value);assert value<=150000,('prospective_SOURCE_packed_owned',value,context)
 def emit(e):
  code=pack(e);assert unpack(code)==tuple(e)
  if code not in work:own(new=True);work.add(code)
 def boundary(k,S):
  own(extra=len(S)+2048);keys=Engine._lex_alias(S);p=out/('boundary-'+str(k)+'.json');h=hashlib.sha256()
  with p.open('wb')as f:
   f.write(b'[');h.update(b'[')
   for i,v in enumerate(keys):
    body=(b','if i else b'')+json.dumps(list(unpack(v)),separators=(',',':')).encode();f.write(body);h.update(body)
   f.write(b']\n');h.update(b']\n')
  del keys;return{'node':k,'positive_support_terms':len(S),'file':p.name,'sha256':h.hexdigest()}
 try:
  gp=ROOT/a.gate_file;assert sha(gp)==a.gate_sha;g=json.loads(gp.read_text());assert g['all_passed'];pins[str(gp.relative_to(ROOT))]=sha(gp)
  bp=gp.parent/'SOURCE-support-lifetime-budget.json';assert sha(bp)==g['evidence_sha256'][bp.name];B=json.loads(bp.read_text());pins[str(bp.relative_to(ROOT))]=sha(bp);assert sha(__file__)==B['worker_sha256']
  rp=ROOT/a.root_review_file;assert sha(rp)==a.root_review_sha and rp.stat().st_size<=65536;r=json.loads(rp.read_text());assert r['all_passed']and count_numeric(r)<=2000;pins[str(rp.relative_to(ROOT))]=sha(rp)
  for f,h in g['input_sha256'].items():assert sha(ROOT/f)==h;pins[f]=h
  D=json.loads((ROOT/FINAL/'compiled-full-SLP-and-fixed-projection-plan.json').read_text());nodes=D['nodes'];degree=D['needed_node_degree'];active=D['active_tied_HL_children'];order=sorted(map(int,degree));META=B['metadata_base'];assert META==109869
  uses={k:0 for k in order};root=D['minor_B_root']
  for k in order:
   for c in active[str(k)]:uses[c]+=1
  uses[root]+=1;boundaries=[]
  for k in order:
   n=nodes[k];l=n['lead5'];d=degree[str(k)];op=n['op'];work=set();context={'node':k,'op':op,'relative_s_degree':d};own()
   if n['q']=='0':pass
   elif op=='input':
    transient=16152;own();p,shift,decoded=read_leaf(n);assert decoded+2048<=16152
    for t in p['terms']:
     e=tuple(p['shift'][j]+int(t[j+1])+shift[j]for j in range(5));assert tuple(e[:3])>=tuple(l[:3])
     if e[0]==l[0]and e[1]==l[1]and 0<=e[2]-l[2]<=d:emit(e)
    del p;transient=0
   elif op=='copy':
    c=n['args'][0];shift=n['shift5'];assert l[:3]==[nodes[c]['lead5'][j]+shift[j]for j in range(3)]
    for v in cache[c]:
     e=unpack(v)
     if e[2]-nodes[c]['lead5'][2]<=d:emit(tuple(e[j]+shift[j]for j in range(5)))
   elif op=='multiply':
    a0,b0=n['args'];la=nodes[a0]['lead5'];lb=nodes[b0]['lead5'];shift=n['shift5'];assert l[:3]==[la[j]+lb[j]+shift[j]for j in range(3)]
    for v in cache[a0]:
     e=unpack(v)
     for w in cache[b0]:
      f=unpack(w)
      if(e[2]-la[2])+(f[2]-lb[2])>d:continue
      pairs+=1;emit(tuple(e[j]+f[j]+shift[j]for j in range(5)))
   else:
    assert op=='add'
    for c in active[str(k)]:
     for v in cache[c]:
      e=unpack(v)
      if e[0]==l[0]and e[1]==l[1]and 0<=e[2]-l[2]<=d:emit(e)
   cache[k]=work;work=set();counts.append([k,len(cache[k])])
   for c in active[str(k)]:
    uses[c]-=1
    if not uses[c]:cache.pop(c,None)
   if not uses[k]:cache.pop(k,None)
   if k in[D['FULL_F11_original_root'],root]:boundaries.append(boundary(k,cache[k]))
   own();assert time.monotonic()-start<50
  S=cache[root];target=tuple(D['target_full5']);target_in=pack(target)in S
  put('all-projected-SOURCE-node-counts.json',counts)
  put('root-support-and-source-custody.json',{'boundaries':boundaries,'target_exponent5':list(target),'target_in_positive_support':target_in,'full_F11_preserved':True,'no_signed_coefficient_values':True,'all5coords_retained_SOURCE_packing_only':True})
  raw=sum((ROOT/f).stat().st_size for f in pins)+sum(p.stat().st_size for p in out.iterdir())+262144;assert raw<=B['fresh_RAW_upper']
  report={'id':ID,'all_passed':True,'SOURCE_support_only_no_signed_values':True,'positive_root_slots':len(S),'target_in_positive_support':target_in,'source_owned_peak':peak,'SOURCE_global_guard':150000,'packed_SOURCE_only_future_VALUE7_unchanged':True,'product_pairs_after_relative_degree_filter':pairs,'fresh_RAW_upper':raw,'no_nonzero_rank_pivot_claim':True,'input_sha256':pins,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};put('report.json',report);signal.alarm(0);print('report_sha256',sha(out/'report.json'));print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}))
 except Exception as exc:
  # Save sets without sort/alias allocation; arbitrary-order exact packed40keys.
  p=out/'live-SOURCE-sets.q5stream'
  with p.open('wb')as f:
   for k,S in cache.items():
    f.write(enc({'node':k,'terms':len(S),'complete_node':True,'next_bytes':5*len(S)}))
    for code in S:f.write(code.to_bytes(5,'big'))
   if work:
    f.write(enc({'node':context.get('node'),'terms':len(work),'complete_node':False,'next_bytes':5*len(work)}))
    for code in work:f.write(code.to_bytes(5,'big'))
  put('failure.json',{'all_passed':False,'error':repr(exc),'context':context,'counts':counts,'no_automatic_retry':True,'no_signed_values':True});raise
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--prepare-gate',action='store_true');p.add_argument('--source-profile',action='store_true');p.add_argument('--gate-file');p.add_argument('--gate-sha');p.add_argument('--root-review-file');p.add_argument('--root-review-sha');a=p.parse_args()
 if a.prepare_gate:prepare()
 else:assert a.source_profile;execute(a)

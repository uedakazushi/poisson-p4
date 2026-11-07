"""Prepared complete leading-HL SOURCE support stream; no new signed VALUE."""
from pathlib import Path
import argparse,hashlib,json,resource,signal,time
from s61_col53_projection_leaf_io_20261105_v64 import read_leaf,leaf_paths
from s61_sparse_generic5_source_packed5_profile_v2_20261105_v64 import pack,unpack,Engine
from s61_sparse_generic5_source_view_profile_20261105_v64 import count_numeric
ROOT=Path(__file__).resolve().parents[1]
ID='S61_COL53_NEXT17_79_FULL_HL_SUPPORT_20261105_v64'
GATE='S61_COL53_NEXT17_79_FULL_HL_SUPPORT_GATE_20261105_v64'
FACE='reports/S61_COL53_NEXT17_79_FULL_HL_FACE_SOURCE_20261105_v64'
FINAL='reports/S61_COL53_NEXT17_79_FINALIZED_SOURCE_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def enc(d):return(json.dumps(d,separators=(',',':'))+'\n').encode()
def prepare():
 out=ROOT/'reports'/GATE;out.mkdir(exist_ok=False);pins={}
 for f,h in[(FACE+'/report.json','c7e13694fe97d78c590b894a27c85231bf4397f00a470d26f83af4567412f332'),(FINAL+'/report.json','548ebf5f4078707d6ee2bb6dd1867c7161c2d32bd553ebef4ce454600c7db877')]:assert sha(ROOT/f)==h;pins[f]=h
 for folder,name in[(FACE,'complete-leading-HL-extrema-and-band-reuse.json'),(FINAL,'compiled-full-SLP-and-fixed-projection-plan.json')]:
  r=json.loads((ROOT/folder/'report.json').read_text());f=folder+'/'+name;pins[f]=sha(ROOT/f);assert pins[f]==r['evidence_sha256'][name]
 D=json.loads((ROOT/FINAL/'compiled-full-SLP-and-fixed-projection-plan.json').read_text());P=json.loads((ROOT/FACE/'complete-leading-HL-extrema-and-band-reuse.json').read_text());nodes=D['nodes'];U=[0]*len(nodes)
 for row in P['all_node_metadata']:U[row['node']]=row['positive_support_cardinality_UPPER_not_exact']
 order=[row['node']for row in P['all_node_metadata']];assert len(order)==1884 and order==sorted(map(int,D['needed_node_degree']))
 assert P['exact_positive_face_S_min_max']==[7,38]and P['positive_root_support_UPPER_not_exact']==1472
 del P
 phases=[];peak=0;pairupper=0;decodepeak=0
 for k in order:
  n=nodes[k];ch=n['args']if n['op']in['copy','multiply']else[c for c in n['args']if nodes[c]['q']!='0'and nodes[c]['lead5'][:2]==n['lead5'][:2]]
  ct=sum(U[c]for c in set(ch));u=U[k]
  states={'operation':118000+ct+u+4096,'child_decode_stream':118000+ct+4096,'canonical_output_after_child_release':118000+2*u+4096}
  if n['op']=='input'and n['q']!='0':
   for f in leaf_paths(n):pins[f]=sha(ROOT/f)
   p,shift,dec=read_leaf(n);assert dec+2048<=16152;decodepeak=max(decodepeak,dec+2048);del p
   states['input_decode']=118000+16152+u+4096
  if n['op']=='multiply':pairupper+=U[n['args'][0]]*U[n['args'][1]]
  peak=max(peak,*states.values());assert max(states.values())<150000
  phases.append([k,states])
 # Mathematical scalar inventory: source graph, order, bound vector, thin
 # source-cell index and counts. Both review docs/control fit the remaining7187.
 exactmeta=25000+65701+1884+2998+4*1884+2*1884+4096;assert exactmeta==110983<118000
 generated=5*sum(U)+512*1884+64*1472+524288+65536+5*1472
 for f in['code/s61_col53_next17_79_full_HL_support_20261105_v64.py','code/s61_col53_projection_leaf_io_20261105_v64.py','code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py','code/s61_sparse_generic5_source_view_profile_20261105_v64.py','code/s61_sparse_generic5_source_packed5_profile_v2_20261105_v64.py']:pins[f]=sha(ROOT/f)
 raw=sum((ROOT/f).stat().st_size for f in pins)+generated+524288+65536+262144;assert raw<32*1024*1024
 (out/'support-cardinality-upper-vector.json').write_bytes(enc(U))
 B={'worker_sha256':sha(__file__),'metadata_base':118000,'explicit_scalar_inventory':exactmeta,'whole_SOURCE_upper':peak,'same_global_guard':150000,'decode_upper':16152,'generated_RAW_upper':generated,'fresh_RAW_upper':raw,'pair_count_UPPER':pairupper,'SOURCE_only_packed40_full5_each_product_checked':True,'future_VALUE7_unchanged':True,'math_seconds':50,'signal_seconds':55,'hard_seconds':60,'AS_bytes':2000000000,'ROOT_doc_bytes':65536,'ROOT_doc_numeric':2000}
 (out/'whole-SOURCE-support-stream-budget.json').write_bytes(enc(B));(out/'all1884-phase-forecasts.json').write_bytes(enc(phases))
 d={'id':GATE,'all_passed':True,'SOURCE_PREPARATION_ONLY_support_profile_UNRUN':True,'no_signed_VALUE_authorization':True,'full_leading_HL_face_S_range':[7,38],'same1884_original_nodes':True,'SOURCE_numeric_upper':peak,'fresh_RAW_upper':raw,'input_sha256':pins,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()}}
 (out/'report.json').write_bytes(enc(d));print('report_sha256',sha(out/'report.json'));print('source/RAW/pairs',peak,raw,pairupper)
def execute(a):
 assert a.gate_file and a.gate_sha and a.root_review_file and a.root_review_sha
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 pins={};index={};counts=[];held=[];work=set();transient=0;peak=0;pairs=0;context={};spool=None;manifest=None;B={'generated_RAW_upper':65536}
 def put(name,d):(out/name).write_bytes(enc(d))
 def alarm(s,f):raise TimeoutError('complete-HL SOURCE55s; preserve prefix, no retry')
 signal.signal(signal.SIGALRM,alarm);signal.alarm(55)
 def own(extra=4096,new=False):
  nonlocal peak
  seen={id(s):len(s)for s in[work,*held]};v=118000+sum(seen.values())+transient+extra+int(new);peak=max(peak,v);assert v<=150000,('prospective_SOURCE_stream_owner',v,context)
 def emit(e):
  z=pack(e);assert unpack(z)==tuple(e)
  if z not in work:own(new=True);work.add(z)
 def child(k):
  off,nb,card,h=index[k];s=set();held.append(s);spool.seek(off);digest=hashlib.sha256()
  for j in range(card):
   body=spool.read(5);assert len(body)==5;digest.update(body);z=int.from_bytes(body,'big');e=unpack(z);assert pack(e)==z and z not in s;own(new=True);s.add(z)
  assert nb==5*card and digest.hexdigest()==h;return s
 try:
  gp=ROOT/a.gate_file;assert sha(gp)==a.gate_sha;g=json.loads(gp.read_text());assert g['all_passed'];pins[a.gate_file]=a.gate_sha
  bp=gp.parent/'whole-SOURCE-support-stream-budget.json';assert sha(bp)==g['evidence_sha256'][bp.name];B=json.loads(bp.read_text());pins[str(bp.relative_to(ROOT))]=sha(bp);assert B['worker_sha256']==sha(__file__)
  rp=ROOT/a.root_review_file;assert sha(rp)==a.root_review_sha and rp.stat().st_size<=65536;r=json.loads(rp.read_text());assert r['all_passed']and count_numeric(r)<=2000;pins[a.root_review_file]=a.root_review_sha
  assert r['input_sha256'][a.gate_file]==a.gate_sha and r['input_sha256']['code/s61_col53_next17_79_full_HL_support_20261105_v64.py']==B['worker_sha256']
  for f,h in g['input_sha256'].items():assert sha(ROOT/f)==h;pins[f]=h
  uv=gp.parent/'support-cardinality-upper-vector.json';assert sha(uv)==g['evidence_sha256'][uv.name];U=json.loads(uv.read_text());pins[str(uv.relative_to(ROOT))]=sha(uv)
  D=json.loads((ROOT/FINAL/'compiled-full-SLP-and-fixed-projection-plan.json').read_text());ns=D['nodes'];order=sorted(map(int,D['needed_node_degree']));root=D['minor_B_root']
  spool=(out/'all-node-positive-supports.q5keys').open('w+b');manifest=(out/'source-node-index.jsonl').open('wb')
  for pos,k in enumerate(order):
   n=ns[k];op=n['op'];lead=n['lead5'];held=[];work=set();memo={};context={'node':k,'order_position':pos,'phase':'SOURCE_operation'};own()
   if n['q']=='0':pass
   elif op=='input':
    transient=16152;own();p,shift,dec=read_leaf(n);assert dec+2048<=16152
    for t in p['terms']:
     e=tuple(p['shift'][j]+int(t[j+1])+shift[j]for j in range(5))
     if e[:2]==tuple(lead[:2]):emit(e)
    del p;transient=0
   else:
    ch=n['args']if op in['copy','multiply']else[c for c in n['args']if ns[c]['q']!='0'and ns[c]['lead5'][:2]==lead[:2]]
    for c in dict.fromkeys(ch):memo[c]=child(c)
    if op=='copy':
     for z in memo[n['args'][0]]:emit(tuple(v+s for v,s in zip(unpack(z),n['shift5'])))
    elif op=='multiply':
     aa,bb=n['args']
     for z in memo[aa]:
      e=unpack(z)
      for zz in memo[bb]:
       f=unpack(zz);pairs+=1;emit(tuple(e[j]+f[j]+n['shift5'][j]for j in range(5)))
    else:
     assert op=='add'
     for c in ch:
      for z in memo[c]:emit(unpack(z))
   assert len(work)<=U[k];memo.clear();held.clear();context['phase']='SOURCE_atomic_segment';own(extra=len(work)+4096)
   keys=sorted(work);spool.seek(0,2);off=spool.tell();hh=hashlib.sha256()
   assert off+5*len(keys)+manifest.tell()+512+64*1472+589824<=B['generated_RAW_upper']
   for z in keys:b=z.to_bytes(5,'big');spool.write(b);hh.update(b)
   spool.flush();rec={'node':k,'offset':off,'bytes':5*len(keys),'positive_terms':len(keys),'sha256':hh.hexdigest()};line=enc(rec);assert len(line)<=512;manifest.write(line);manifest.flush();index[k]=(off,rec['bytes'],len(keys),rec['sha256']);counts.append([k,len(keys)]);del keys;work=set()
   assert time.monotonic()-start<50
  S=child(root);work=S;held.clear();own(extra=len(S)+4096);keys=Engine._lex_alias(S);path=out/'complete-leading-HL-positive-support.json'
  with path.open('wb')as f:
   f.write(b'[')
   for i,z in enumerate(keys):f.write((b','if i else b'')+json.dumps(list(unpack(z)),separators=(',',':')).encode())
   f.write(b']\n')
  del keys;assert min(unpack(z)[2]for z in S)==7 and max(unpack(z)[2]for z in S)==38
  assert all(unpack(z)[:2]==(-5,0)and unpack(z)[4]==19 for z in S)
  spool.close();manifest.close();spool=None;manifest=None;put('all1884-exact-complete-face-source-counts.json',counts)
  generated=sum(p.stat().st_size for p in out.iterdir());assert generated<=B['generated_RAW_upper'];raw=sum((ROOT/f).stat().st_size for f in pins)+generated+262144;assert raw<=B['fresh_RAW_upper']
  report={'id':ID,'all_passed':True,'complete_leading_HL_SOURCE_only_no_signed_values':True,'complete_nodes':1884,'positive_root_terms':len(S),'positive_S_range':[7,38],'source_owned_peak':peak,'pair_operations':pairs,'fresh_RAW_upper':raw,'old_signed_prefix_cells_unchanged_reusable':True,'no_face_zero_nonzero_or_pivot_claim':True,'input_sha256':pins,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};put('report.json',report);signal.alarm(0);print('report_sha256',sha(out/'report.json'));print('root_SOURCE_terms',len(S))
 except Exception as exc:
  if spool is not None:spool.flush()
  if manifest is not None:manifest.flush()
  with(out/'incomplete-current-positive-support.q5keys').open('wb')as f:
   for z in work:f.write(z.to_bytes(5,'big'))
  put('failure.json',{'all_passed':False,'SOURCE_only_no_signed_values':True,'error':repr(exc),'context':context,'completed_source_counts':counts,'complete_segments_and_current_partial_preserved':True,'no_automatic_retry':True});raise
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--prepare-gate',action='store_true');p.add_argument('--source-profile',action='store_true');p.add_argument('--gate-file');p.add_argument('--gate-sha');p.add_argument('--root-review-file');p.add_argument('--root-review-sha');a=p.parse_args()
 if a.prepare_gate:prepare()
 else:assert a.source_profile;execute(a)

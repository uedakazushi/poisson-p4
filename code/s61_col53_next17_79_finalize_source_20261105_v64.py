"""Changed metadata-only finalization from immutable completed compiler prefix."""
from pathlib import Path
import argparse,hashlib,json,re,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_COL53_NEXT17_79_FINALIZED_SOURCE_20261105_v64'
GATE='S61_COL53_NEXT17_79_FINALIZE_SOURCE_GATE_20261105_v64'
PREFIX='reports/S61_COL53_NEXT17_79_COMPILED_SOURCE_20261105_v64/partial-compiled-source.json'
PREFIX_SHA='9a20b33d3af0ec833be03abf2f88dd2ef002ad70f5865c9db55871c9569b1a62'
PLAN='reports/S61_COL53_NEXT17_79_FUNCTIONAL_SOURCE_20261105_v64/chosen-functional-and-thin-source-inputs.json'
from s61_sparse_generic5_source_view_profile_20261105_v64 import count_numeric
def enc(x):return(json.dumps(x,separators=(',',':'))+'\n').encode()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def census(text):
 token=re.compile(r'"(?:\\.|[^"\\])*"|-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?|true|false|null');digits=re.compile(r'[+-]?[0-9]+(?:/[0-9]+)?');n=0
 for m in token.finditer(text):
  t=m.group()
  if t[0]=='"':n+=3 if digits.fullmatch(t[1:-1])else 0
  elif t not in['true','false','null']:n+=1
 return n
def gate():
 out=ROOT/'reports'/GATE;out.mkdir(exist_ok=False);p=ROOT/PREFIX;assert sha(p)==PREFIX_SHA
 text=p.read_text();n=census(text);assert n==86903
 nb=0;nonnull=0
 for m in re.finditer(r'"box5":(null|\[\[.*?\]\])',text):nb+=1;nonnull+=m.group(1)!='null'
 nn=sum(1 for _ in re.finditer(r'"n_upper":-?\d+',text));assert(nb,nonnull,nn)==(2998,2320,2998)
 compact=n-10*nonnull-nn;assert compact==60705
 cold=25000+n+4096;mode=25000+compact+8*nb+12000+4000
 assert cold==115999 and mode==125689 and max(cold,mode)<126000<150000
 # Preserve the ancestor read-only hash closure without parsing the graph array.
 marker='"input_sha256":';i=text.index(marker)+len(marker);closure=json.JSONDecoder().raw_decode(text,i)[0]
 closure[PREFIX]=PREFIX_SHA
 closure['code/s61_col53_next17_79_compile_source_20261105_v64.py']=sha(ROOT/'code/s61_col53_next17_79_compile_source_20261105_v64.py')
 closure['reports/S61_COL53_NEXT17_79_COMPILED_SOURCE_20261105_v64/failure.json']=sha(ROOT/'reports/S61_COL53_NEXT17_79_COMPILED_SOURCE_20261105_v64/failure.json')
 closure[PLAN]=sha(ROOT/PLAN)
 for f in['code/s61_col53_next17_79_finalize_source_20261105_v64.py','code/s61_sparse_generic5_source_view_profile_20261105_v64.py','code/s61_sparse_generic5_source_packed5_profile_v2_20261105_v64.py']:
  closure[f]=sha(ROOT/f)
 for f,h in closure.items():assert sha(ROOT/f)==h
 raw=sum((ROOT/f).stat().st_size for f in closure)+4*1024*1024+512*1024+65536+262144
 assert raw<32*1024*1024
 budget={'worker_sha256':sha(Path(__file__)),'prefix_sha256':PREFIX_SHA,'nodes':2998,'raw_parse_numeric_upper':86903,'cold_with25k_and4096':cold,'compact_graph_numeric_upper':compact,'max_mode_index_and_stack':8*nb+12000,'literal_reroute_and_extra_metadata':4000,'whole_numeric_upper':126000,'precise_mode_upper':mode,'ROOT_doc_bytes':65536,'ROOT_doc_numeric':2000,'generated_RAW_upper':4*1024*1024,'fresh_RAW_upper':raw,'no_Source_arithmetic_replay':True,'no_support_or_signed_values':True}
 (out/'read-only-finalization-budget.json').write_bytes(enc(budget))
 d={'id':GATE,'all_passed':True,'SOURCE_gate_only_no_finalization_executed':True,'preserves_old157050_stop':True,'changed_method':'new worker parses saved2998-node prefix, in-place drops auxiliaryBOX/n; no oldcompiledscalaroperations replay/CSE/parse-doc retention','worker_sha256':budget['worker_sha256'],'input_sha256':closure,'evidence_sha256':{'read-only-finalization-budget.json':sha(out/'read-only-finalization-budget.json')},'whole_numeric_upper':126000,'fresh_RAW_upper':raw,'ROOT_review_GO_required_before_finalize':True}
 (out/'report.json').write_bytes(enc(d));print('report_sha256',sha(out/'report.json'));print('whole/RAW',126000,raw)
def finalize(a):
 assert a.gate_file and a.gate_sha and a.root_review_file and a.root_review_sha
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);nodes=None;pins={};context={};peak=0
 def put(n,d):
  b=enc(d);assert sum(p.stat().st_size for p in out.iterdir())+len(b)<=4*1024*1024
  (out/n).write_bytes(b)
 def alarm(s,f):raise TimeoutError('SOURCEmetadata55s; preserve prefix; noautomaticretry')
 signal.signal(signal.SIGALRM,alarm);signal.alarm(55)
 try:
  gp=ROOT/a.gate_file;assert sha(gp)==a.gate_sha;g=json.loads(gp.read_text());assert g['all_passed'];pins[str(gp.relative_to(ROOT))]=a.gate_sha
  bp=gp.parent/'read-only-finalization-budget.json';assert sha(bp)==g['evidence_sha256'][bp.name];B=json.loads(bp.read_text());assert sha(__file__)==B['worker_sha256'];pins[str(bp.relative_to(ROOT))]=sha(bp)
  rp=ROOT/a.root_review_file;assert sha(rp)==a.root_review_sha and rp.stat().st_size<=65536;review=json.loads(rp.read_text());assert review['all_passed']and count_numeric(review)<=2000;pins[str(rp.relative_to(ROOT))]=a.root_review_sha
  for f,h in g['input_sha256'].items():assert sha(ROOT/f)==h;pins[f]=h
  p=ROOT/PREFIX;D=json.loads(p.read_text());assert count_numeric(D)<=B['raw_parse_numeric_upper'];peak=B['cold_with25k_and4096'];nodes=D['nodes_completed'];assert len(nodes)==2998
  for k,n in enumerate(nodes):
   assert all(isinstance(j,int)and 0<=j<k for j in n['args'])
   n.pop('box5');n.pop('n_upper')
  assert count_numeric(D)<=B['compact_graph_numeric_upper']
  # Pure byte-equal input alias to the existing copied pack; no old14MB reread.
  P=json.loads((ROOT/PLAN).read_text());bysha={d['payload_sha256']:d for d in P['selected48_literal_payloads']if'payload_sha256'in d};rerouted=0
  thinfile='reports/S61_COL53_NEXT17_79_FUNCTIONAL_SOURCE_20261105_v64/selected48-existing-cell-payloads.bin'
  for n in nodes:
   if n['op']!='input':continue
   src=n.get('source',{});old=src.get('legacy_source',{})
   if old.get('file','').endswith('current-matrix.bin'):
    d=bysha[old['payload_sha256']];assert old['bytes']==d['bytes']
    src['same_original_accepted_source']=old
    src['legacy_source']={'thin_payload_file':thinfile,'offset':d['new_offset'],'bytes':d['bytes'],'payload_sha256':d['payload_sha256']}
    rerouted+=1
  root=len(nodes)-1;assert nodes[root]['name']=='B17_79'and nodes[root]['op']=='add'and nodes[root]['scalar']=='-1'
  diag,cross=nodes[root]['args'];assert nodes[diag]['op']==nodes[cross]['op']=='multiply'
  F11,droot=nodes[diag]['args'];broot,croot=nodes[cross]['args'];assert F11==1202 and nodes[F11]['name']=='new_1013_sum'
  assert nodes[root]['lead5']==[-5,0,7,-3,19]
  target=P['new_functional_request']['exponent5'];assert target==[-5,0,9,-3,19]
  degree={};active={}
  def own():
   nonlocal peak
   numeric=25000+count_numeric(D)+count_numeric(P)+8*len(nodes)+12000
   peak=max(peak,numeric);assert numeric<=B['whole_numeric_upper']<=150000
  def visit(k,d):
   if degree.get(k,-1)>=d:return
   degree[k]=d;n=nodes[k];children=[]
   if n['q']=='0':active[k]=[];return
   if n['op']in['copy','multiply']:children=[(c,d)for c in n['args']]
   elif n['op']=='add':
    for c in n['args']:
     lead=nodes[c]['lead5'];parent=n['lead5']
     if lead is not None and lead[:2]==parent[:2]and lead[2]-parent[2]<=d:
      assert lead[2]>=parent[2];children.append((c,d-(lead[2]-parent[2])))
   active[k]=[c for c,dd in children]
   for c,dd in children:visit(c,dd)
  own();visit(root,2);own()
  maxM=max(int(n['M']).bit_length()for n in nodes);assert maxM<=2048
  graph={'nodes':nodes,'FULL_F11_original_root':F11,'four_callback_roots':{'6,53':F11,'17,79':droot,'17,53':broot,'6,79':croot},'minor_B_root':root,'actual_C_factor':'(P7*F9)^2','target_full5':target,'max_relative_s_degree':2,'needed_node_degree':{str(k):v for k,v in degree.items()},'active_tied_HL_children':{str(k):v for k,v in active.items()},'full_original_prefix_preserved':{'file':PREFIX,'sha256':PREFIX_SHA},'auxiliaryBOX_n_only_omitted_from_resident_view':True,'old157050_stop_unchanged':True,'whole_F11_not_face411':True,'no_support_or_signed_values_evaluated':True}
  put('compiled-full-SLP-and-fixed-projection-plan.json',graph)
  raw=sum((ROOT/f).stat().st_size for f in pins)+sum(p.stat().st_size for p in out.iterdir())+262144;assert raw<=B['fresh_RAW_upper']
  r={'id':ID,'all_passed':True,'metadata_finalization_only_no_Source_arithmetic_replay':True,'full_nodes':len(nodes),'needed_projection_nodes':len(degree),'full_F11_preserved':True,'new_raw_payload_aliases_exact':rerouted,'source_M_bits_upper':maxM,'whole_numeric_upper':B['whole_numeric_upper'],'actual_prospective_numeric_peak':peak,'fresh_RAW_upper':raw,'projection_SOURCE_support_and_signed_VALUE_STILL_UNRUN':True,'input_sha256':pins,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};put('report.json',r);signal.alarm(0);print('report_sha256',sha(out/'report.json'));print(json.dumps({k:v for k,v in r.items()if k not in('input_sha256','evidence_sha256')}))
 except Exception as exc:
  if nodes is not None:put('metadata-prefix-on-stop.json',{'nodes':nodes,'complete_finalization_NOT_claimed':True,'input_sha256':pins})
  put('failure.json',{'all_passed':False,'metadata_only':True,'error':repr(exc),'context':context,'no_automatic_retry':True,'no_support_signed_values':True});raise
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--prepare-gate',action='store_true');p.add_argument('--finalize',action='store_true');p.add_argument('--gate-file');p.add_argument('--gate-sha');p.add_argument('--root-review-file');p.add_argument('--root-review-sha');a=p.parse_args()
 if a.prepare_gate:gate()
 else:assert a.finalize;finalize(a)

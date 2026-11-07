#!/usr/bin/env python3
"""Read-only local identity postpeer. Signed verification REQUIRES fresh ROOT GO."""
import argparse,hashlib,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRIMARY='S61_COL53_NEXT17_79_STREAMED_VALUE_20261105_v64'
PREP=PRIMARY+'_LOCAL_IDENTITY_POSTPEER_GATE'
ID=PRIMARY+'_LOCAL_IDENTITY_REVIEW'
CODE='code/s61_col53_next17_79_local_identity_postpeer_20261105_v64.py'
FINAL='reports/S61_COL53_NEXT17_79_FINALIZED_SOURCE_20261105_v64/compiled-full-SLP-and-fixed-projection-plan.json'
COUNTS='reports/S61_COL53_NEXT17_79_SOURCE_PROJECTION_20261105_v64/all-projected-SOURCE-node-counts.json'
TMP='work/'+PRIMARY+'-canonical-cells'
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def qparse(s,cap=1024):
 a=str(s).split('/');assert len(a)<=2
 n=int(a[0]);d=int(a[1])if len(a)==2 else 1;assert d>0 and max(abs(n).bit_length(),d.bit_length())<=cap
 return Q(n,d)
def decode_cell(data,qi,M,C):
 """Own bounded canonical parser, returning fixed-qi INTEGER coefficients."""
 assert data[:8]==b'Q5CELL1\n'and len(data)>=9
 if data[8]==0:assert len(data)==9;return{}
 assert data[8]==1 and qi;at=9
 def u():
  nonlocal at
  start=at;n=0;factor=1
  while True:
   assert at<len(data)and at-start<5;c=data[at];at+=1;n+=(c%128)*factor
   if c<128:break
   factor*=128
  assert at-start==max(1,(n.bit_length()+6)//7);return n
 def integer(bitcap):
  nonlocal at
  tag=u();size=tag//2;assert size<=(bitcap+7)//8 and at+size<=len(data)
  b=data[at:at+size];at+=size;assert not b or b[0];n=int.from_bytes(b,'big');assert n.bit_length()<=bitcap and(n or not tag%2)
  return-n if tag%2 else n
 n=integer(852);d=integer(164);assert n and d>0 and math.gcd(n,d)==1
 assert abs(n).bit_length()+qi.denominator.bit_length()<=1921 and d.bit_length()+abs(qi.numerator).bit_length()<=1921
 ratio=Q(n*qi.denominator,d*qi.numerator);assert ratio.denominator==1 and 0<ratio.numerator<=M;G=ratio.numerator
 shift=[]
 for k in range(5):z=u();shift.append(z//2 if z%2==0 else-(z//2)-1)
 N=u();assert 0<N<=C;v={};last=None;gc=0
 for k in range(N):
  c=integer(824);assert c;gc=math.gcd(gc,abs(c));e=tuple(shift[j]+u()for j in range(5));assert all(-64<=x<=191 for x in e)
  assert last is None or(sum(last),last)>(sum(e),e);last=e;assert e not in v
  assert abs(c).bit_length()+G.bit_length()<=1921;x=c*G;assert abs(x)<=M;v[e]=x
 assert gc==1 and at==len(data);return v
def weights(n,nodes):
 qi=qparse(n['q']);a=n['args']
 if not qi:return[]
 if n['op']=='copy':xs=[qparse(nodes[a[0]]['q'])*qparse(n['scalar'])/qi]
 elif n['op']=='multiply':xs=[qparse(nodes[a[0]]['q'])*qparse(nodes[a[1]]['q'])*qparse(n['scalar'])/qi]
 elif n['op']=='add':xs=[qparse(nodes[c]['q'])*(qparse(n['scalar'])if j else Q(1))/qi for j,c in enumerate(a)]
 else:return[]
 assert all(x.denominator==1 for x in xs);return[x.numerator for x in xs]
def tree_digest(order):
 h=hashlib.sha256();total=0
 for k in order:
  p=ROOT/TMP/f'node-{k}.bin';z=sha(p);n=p.stat().st_size;total+=n;h.update(f'{k} {n} {z}\n'.encode())
 return h.hexdigest(),total
def prepare():
 start=time.monotonic();out=ROOT/'reports'/PREP;out.mkdir(exist_ok=False);pins={}
 def read(p,h=None):
  z=sha(ROOT/p)
  if h:assert z==h
  pins[p]=z;return json.loads((ROOT/p).read_text())
 pp=f'reports/{PRIMARY}/report.json';r=read(pp,'4b1ed0275771c6f13a8f250b14a71a125277d3b8e6f17fd5995142b8d4dc0c54');assert r['all_passed']and r['nodes_committed']==1884
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 for p,h in r['evidence_sha256'].items():p=f'reports/{PRIMARY}/{p}';assert sha(ROOT/p)==h;pins[p]=h
 peer=read('reports/S61_COL53_NEXT17_79_STREAMED_VALUE_GATE_20261105_v64_REVIEW/report.json','012f7d578a5240351e40451479f4caeb22997b85247a9e2956de00f572cab27a');assert peer['all_passed']
 D=read(FINAL);C=dict(read(COUNTS));order=sorted(C);assert len(C)==1884 and max(C.values())==129 and count(D)==51115
 maxraw=0
 for k in order:
  n=D['nodes'][k];M=int(n['M']);assert M.bit_length()<=824
  if not qparse(n['q'])or not C[k]:continue
  w=weights(n,D['nodes']);args=n['args']
  if n['op']=='copy':assert abs(w[0])*int(D['nodes'][args[0]]['M'])<=M
  elif n['op']=='multiply':assert min(C[a]for a in args)*abs(w[0])*math.prod(int(D['nodes'][a]['M'])for a in args)<=M
  elif n['op']=='add':assert sum(abs(x)*int(D['nodes'][a]['M'])for x,a in zip(w,args))<=M
  maxraw=max(maxraw,2*824,852+164)
 tree,treebytes=tree_digest(order);stream=ROOT/'reports'/PRIMARY/'node-values.q5stream';index=ROOT/'reports'/PRIMARY/'node-values-index.jsonl'
 assert treebytes==stream.stat().st_size and treebytes<=9268352 and index.stat().st_size<=2048*1884
 assert sha(ROOT/CODE);pins[CODE]=sha(ROOT/CODE)
 oldread=sum((ROOT/p).stat().st_size for p in r['input_sha256'])
 RAW=oldread+23968000+2_000_000+2*65536;assert RAW<32*1024*1024
 # The new layout has one extra saved-current map besides the two children
 # and expected accumulator; prospective C bounds charge all full5 keys.
 phase={};at={}
 for k in order:
  n=D['nodes'][k];ni=C[k];cs=list(dict.fromkeys(n['args']if n['op']in['copy','multiply']else D['active_tied_HL_children'][str(k)]if n['op']=='add'else[]));ct=sum(C[c]for c in cs);mc=max([C[c]for c in cs]+[ni])
  for name,z in [('maps',118000+7*(ct+2*ni)+4096),('oracle_decode',118000+7*(ct+ni)+8*mc+2048+4096),('root_copy',118000+7*(ct+3*ni)+4096)]:
   if z>phase.get(name,0):phase[name]=z;at[name]=k
  if n['op']=='input':
   for name,z in[('leaf',118000+21000+14*ni+4096),('leaf_STOP',118000+21000+17*ni+4096)]:
    if z>phase.get(name,0):phase[name]=z;at[name]=k
 assert max(phase.values())<150000
 fixtures=out/'unrelated-bounded-parser';fixtures.mkdir()
 # Only pre-existing unrelated codec fixture bytes, no selected node cells.
 fp='reports/S61_COL53_NEXT17_79_STREAMED_VALUE_GATE_20261105_v64_REVIEW/unrelated-codec-fixtures'
 for k,qi,M in[(0,Q(-7,30),18),(1,Q(11,37),25),(2,Q(-13,21),1),(3,Q(0),1)]:
  p=ROOT/fp/f'case-{k}.bin';pins[str(p.relative_to(ROOT))]=sha(p);v=decode_cell(p.read_bytes(),qi,M,2);assert len(v)==(2 if k==0 else 1 if k==1 else 0)
 meta_bound=count(D)+2*1884+1884+4*1884+3*1884+25000+2*2000+3000
 assert meta_bound<=118000
 e={'all_passed':True,'SOURCE_PREP_ONLY_no_signed_selected_node_arithmetic':True,'whole_numeric_prior':max(phase.values()),'META':118000,'independent_metadata_bound_including1884_report_rows':meta_bound,'hard_numeric_guard':150000,'phase_bounds':phase,'peak_nodes':at,
  'source_M_bits824_fixedqi28_164':True,'source_commit_bits':824,'max_raw_integer_prior':1921,'hard_VALUE_bits':2048,
  'canonical_tree_digest':tree,'canonical_tree_bytes':treebytes,'canonical_files':1884,'thin_index_4N_not_full_manifest':True,'both_generations_compared_at_runtime':True,
  'source_C_cards_full5_VALUE7':True,'input_single_literal_representation_transient':21000,'source_input_decoded14104_implies_atmost2349_terms':True,'all1884_report_rows_charged_META':True,'whole_RAW_upper':RAW,'generated_output_reserve':2000000,'AS_bytes':2000000000,'seconds_guard':55,
  'local_identity_method':'Read savedcurrent and savedchildren with own bounded canonical parser, restore fixedqi integer coefficients; independently evaluate COPY/ADD/MUL in reverse term order; check exact coefficient dictionary equality at each1884node. Leaves rebound to original immutable physical data. No producer evaluator/replay.',
  'scope':'Read-only actual postpeer, ROOTGO still required. Target/full projection only, no wholeminorzero/rank/matrixupdate.'}
 ep=out/'positive-local-identity-prior-and-custody.json';save(ep,e)
 rep={'id':PREP,'all_passed':True,'SOURCE_PREPARATION_ONLY':True,'input_sha256':pins,'code_sha256':sha(ROOT/CODE),'evidence_sha256':{ep.name:sha(ep)},'signed_postpeer_UNRUN':True,'seconds':time.monotonic()-start};save(out/'report.json',rep);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'numeric':max(phase.values()),'RAW':RAW}))
def verify(args):
 start=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};context={};work={};rows=[]
 def read(p,h=None):
  z=sha(ROOT/p)
  if h:assert z==h
  pins[p]=z;return json.loads((ROOT/p).read_text())
 gate=read(f'reports/{PREP}/report.json',args.gate_sha);assert gate['all_passed']and gate['code_sha256']==sha(ROOT/CODE)
 bname=next(iter(gate['evidence_sha256']));B=read(f'reports/{PREP}/{bname}',gate['evidence_sha256'][bname])
 assert(ROOT/args.root_report).stat().st_size<=65536
 go=read(args.root_report,args.root_sha);assert count(go)<=2000 and go['all_passed']and go['input_sha256'][CODE]==gate['code_sha256']and go['input_sha256'][f'reports/{PREP}/report.json']==args.gate_sha
 for p,h in gate['input_sha256'].items():assert sha(ROOT/p)==h
 D=read(FINAL);nodes=D['nodes'];C=dict(read(COUNTS));order=sorted(C);index={};stream=ROOT/'reports'/PRIMARY/'node-values.q5stream';offset=0
 with(ROOT/'reports'/PRIMARY/'node-values-index.jsonl').open()as f:
  for line in f:
   assert len(line)<=2048;r=json.loads(line);k=r['node'];assert k not in index and r['stream_offset']==offset and r['actual_nonzero_slots']<=C[k]
   index[k]=(r['stream_offset'],r['bytes'],r['actual_nonzero_slots'],r['payload_sha256']);offset+=r['bytes']
 assert sorted(index)==order and offset==stream.stat().st_size
 htree=hashlib.sha256();peak=0;held=[];oracle={};work={};transient=0
 def own(extra=4096,new=False):
  nonlocal peak
  n=118000+7*(sum(len(v)for v in held)+len(work)+len(oracle)+int(new))+transient+extra;peak=max(peak,n);assert n<=150000,(context,n)
 def signed(x):assert abs(x).bit_length()<=1921;return x
 def load(k,f):
  nonlocal transient
  off,size,n,hs=index[k];p=ROOT/TMP/f'node-{k}.bin';body=p.read_bytes();assert len(body)==size and hashlib.sha256(body).hexdigest()==hs
  f.seek(off);assert f.read(size)==body
  transient=8*C[k]+2048;own();v=decode_cell(body,qparse(nodes[k]['q']),int(nodes[k]['M']),C[k]);assert len(v)==n;transient=0;return v
 def source_leaf(n):
  nonlocal transient
  src=n['source'];s=src.get('legacy_source',src);shift=[0]*5;scale=Q(1);transient=21000;own()
  if'binding'in s:
   a=s['binding'];doc=json.loads((ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'/a['scope']).read_text());v=(doc['common_operands']if'common_operands'in doc else doc['row_operands'])[a['label']];p=v['original_literal'];scale=qparse(v['QQ_transform_scalar']);shift=v['transform_monomial_exponent5']
  elif'thin_payload_file'in s or('file'in s and'offset'in s):
   filename=s.get('thin_payload_file',s.get('file'))
   with(ROOT/filename).open('rb')as f:f.seek(s['offset']);body=f.read(s['bytes'])
   assert hashlib.sha256(body).hexdigest()==s['payload_sha256'];p=payload_to_literal(body)
  elif'cell_file'in s or'file'in s:
   filename=s.get('cell_file',s.get('file'));body=(ROOT/filename).read_bytes();assert body[:9]==b'Q5CELL1\n\1';p=payload_to_literal(body[9:])
   if'legacy_source'in src:
    prof=src['literal_profile'];low=[p['shift'][j]+min(int(t[j+1])for t in p['terms'])for j in range(5)];shift=[prof['exponent_box5'][j][0]-low[j]for j in range(5)];scale=qparse(prof['known_QQ_content'])/qparse(p['content'])
    assert(scale,shift)==((Q(353894400000),[1,0,0,0,1])if'col56-ninth-cleared'in filename else(Q(1),[0]*5))
  elif s.get('QQ_monomial'):
   prof=src.get('literal_profile');m=s.get('exponent5',[b[0]for b in prof['exponent_box5']]if prof else[0]*5);p={'content':n['q'],'shift':m,'terms':[['1',0,0,0,0,0]]}
  else:assert s.get('known_zero')or s.get('known_source_zero');p={'content':'0','shift':[0]*5,'terms':[]}
  assert qparse(p['content'])*scale==qparse(n['q']);return p,shift
 def payload_to_literal(data):
  # Directly construct only ONE six-number literal term representation.
  # Do not first build a full-key map/list and then duplicate its terms.
  at=0
  def u():
   nonlocal at
   start=at;n=0;factor=1
   while True:
    assert at<len(data)and at-start<5;c=data[at];at+=1;n+=(c%128)*factor
    if c<128:break
    factor*=128
   assert at-start==max(1,(n.bit_length()+6)//7);return n
  def integer(cap):
   nonlocal at
   tag=u();size=tag//2;assert size<=(cap+7)//8 and at+size<=len(data);b=data[at:at+size];at+=size;assert not b or b[0]
   n=int.from_bytes(b,'big');assert n.bit_length()<=cap and(n or not tag%2);return-n if tag%2 else n
  n=integer(852);d=integer(512);assert n and d>0 and math.gcd(n,d)==1;content=Q(n,d);shift=[]
  for j in range(5):z=u();shift.append(z//2 if z%2==0 else-(z//2)-1)
  N=u();assert 0<N<=2349;terms=[];previous=None;gc=0
  for k in range(N):
   c=integer(824);assert c;gc=math.gcd(gc,abs(c));offsets=[u()for j in range(5)];e=tuple(shift[j]+offsets[j]for j in range(5))
   assert previous is None or(sum(previous),previous)>(sum(e),e);previous=e;terms.append([str(c),*offsets])
  assert at==len(data)and gc==1;return{'content':str(content),'shift':shift,'terms':terms}
 try:
  with stream.open('rb')as f:
   for k in order:
    context={'node':k,'verified_nodes':len(rows)};held=[];oracle={};work={};transient=0;n=nodes[k];qi=qparse(n['q']);M=int(n['M']);lead=n['lead5'];cap=D['needed_node_degree'][str(k)]
    oracle=load(k,f);held=[];args0=n['args']
    def put(e,x):
     assert len(e)==5 and abs(x)<=M and all(-64<=z<=191 for z in e)
     if x:
      if e not in work:own(new=True)
      work[e]=x
     else:work.pop(e,None)
     assert len(work)<=C[k]
    if qi and C[k]:
     if n['op']=='input':
      p,shift=source_leaf(n)
      for term in reversed(p['terms']):
       e=tuple(p['shift'][j]+int(term[j+1])+shift[j]for j in range(5))
       if e[:2]==tuple(lead[:2])and 0<=e[2]-lead[2]<=cap:put(e,int(term[0]))
      del p;transient=0
     else:
      wanted=args0 if n['op']in['copy','multiply']else D['active_tied_HL_children'][str(k)];memo={}
      for c in reversed(list(dict.fromkeys(wanted))):v=load(c,f);memo[c]=v;held.append(v)
      ws=weights(n,nodes)
      if n['op']=='copy':
       for e,x in reversed(list(memo[args0[0]].items())):
        ee=tuple(e[j]+n['shift5'][j]for j in range(5))
        if ee[:2]==tuple(lead[:2])and 0<=ee[2]-lead[2]<=cap:put(ee,signed(x*ws[0]))
      elif n['op']=='multiply':
       a,b=args0;la=nodes[a]['lead5'];lb=nodes[b]['lead5']
       for e,x in reversed(list(memo[b].items())):
        for g,y in reversed(list(memo[a].items())):
         if(e[2]-lb[2])+(g[2]-la[2])>cap:continue
         ee=tuple(e[j]+g[j]+n['shift5'][j]for j in range(5));put(ee,signed(work.get(ee,0)+signed(signed(x*y)*ws[0])))
      else:
       for j in reversed(range(len(args0))):
        c=args0[j]
        if c not in memo:continue
        for e,x in reversed(list(memo[c].items())):
         if e[:2]==tuple(lead[:2])and 0<=e[2]-lead[2]<=cap:put(e,signed(work.get(e,0)+signed(x*ws[j])))
      memo.clear()
    assert work==oracle,('local_identity',k)
    off,size,nn,hs=index[k];htree.update(f'{k} {size} {hs}\n'.encode());rows.append({'node':k,'actual_terms':len(oracle),'full_coefficient_identity_zero':True,'payload_sha256':hs})
    if k==D['minor_B_root']:rootvals=dict(oracle)
    held.clear();oracle={};work={};own();assert time.monotonic()-start<50
  assert htree.hexdigest()==B['canonical_tree_digest']and len(rows)==1884
  assert not rootvals;target=tuple(D['target_full5']);assert rootvals.get(target,0)==0
  ep=out/'all1884-independent-local-coefficient-identities.json';save(ep,rows)
  report={'id':ID,'all_passed':True,'independent_read_only_signed_local_postpeer':True,'all1884_local_identities_from_trueleaves_and_savedchildren':True,'both_canonical_generations_byteidentical':True,'needed_projection_actual_empty':True,'target_coefficient_QQ':'0','whole_minor_zero_or_new_rank_claimed':False,'actual_numeric_peak':peak,'SOURCE_prior':B['whole_numeric_prior'],'input_sha256':pins,'code_sha256':sha(ROOT/CODE),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start};save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'peak':peak,'seconds':report['seconds']}))
 except BaseException as error:
  save(out/'STOP-only-local-proof-prefix.json',{'all_passed':False,'error':repr(error),'context':context,'all_completed_local_identities':rows,'partial_expected_coefficients':[[list(e),str(x)]for e,x in sorted(work.items())],'no_automatic_retry_or_new_query':True});raise
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--verify',action='store_true');p.add_argument('--gate-sha');p.add_argument('--root-report');p.add_argument('--root-sha');a=p.parse_args()
 if a.prepare:prepare()
 else:assert a.verify and a.gate_sha and a.root_report and a.root_sha;verify(a)

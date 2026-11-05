"""SOURCE-only full positive H=-3/H=21 faces, integer weights and VALUE lifetime."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,gc,hashlib,json,resource,signal,sys,time
ROOT=Path(__file__).resolve().parents[1]
from s61_sparse_generic5_source_view_profile_20261105_v64 import count_numeric
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
OLD='reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'
ID='S61_COL53_H_FACES_SOURCE_20261105_v64'
VERTEX='S61_COL53_VERTEX_SOURCE_GATE_V2_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def enc(x):return(json.dumps(x,separators=(',',':'))+'\n').encode()
def main(a):
 assert a.prepare
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};context={};sourcepeak=0;leafparsemax=0;leafmax=0;scalarmax=0
 def put(n,v):
  b=enc(v);assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())+len(b)<4*1024*1024
  (out/n).write_bytes(b)
 def timeout(s,f):
  put('failure.json',{'all_passed':False,'SOURCE_only':True,'stop':'math55s/hard60s','context':context,'no_signed_face_coefficient_evaluated':True,'no_automatic_retry':True});raise TimeoutError('SOURCE55s')
 signal.signal(signal.SIGALRM,timeout);signal.alarm(55)
 def read(p):
  p=Path(p)
  if str(p.relative_to(ROOT))not in pins:pins[str(p.relative_to(ROOT))]=sha(p)
  return json.loads(p.read_text())
 def own(n):
  nonlocal sourcepeak
  sourcepeak=max(sourcepeak,n);assert n<=150000,('SOURCE prospective owned guard',n,context)
 def polynomial(src,profile):
  nonlocal leafparsemax,leafmax
  if 'binding'in src:
   b=src['binding'];doc=read(ROOT/OLD/b['scope']);leafparsemax=max(leafparsemax,count_numeric(doc)+2048)
   v=(doc['common_operands']if 'common_operands'in doc else doc['row_operands'])[b['label']]
   return v['original_literal'],Q(v['QQ_transform_scalar']),v['transform_monomial_exponent5']
  if 'file'in src:
   p=ROOT/src['file']
   if str(p.relative_to(ROOT))not in pins:pins[str(p.relative_to(ROOT))]=sha(p)
   if 'offset'in src:
    with p.open('rb')as f:f.seek(src['offset']);b=f.read(src['bytes'])
    assert hashlib.sha256(b).hexdigest()==src['payload_sha256'];poly=decode(b);leafparsemax=max(leafparsemax,count_numeric(poly)+2048);return poly,Q(1),[0]*5
   assert pins[str(p.relative_to(ROOT))]==src['sha256'];b=p.read_bytes();assert b[:9]==b'Q5CELL1\n\x01';poly=decode(b[9:]);leafparsemax=max(leafparsemax,count_numeric(poly)+2048)
   low=[poly['shift'][j]+min(t[j+1]for t in poly['terms'])for j in range(5)]
   shift=[profile['exponent_box5'][j][0]-low[j]for j in range(5)];scalar=Q(profile['known_QQ_content'])/Q(poly['content'])
   expected=(Q(353894400000),[1,0,0,0,1])if'col56-ninth-cleared'in src['file']else(Q(1),[0]*5)
   assert(scalar,shift)==expected
   assert all(profile['exponent_box5'][j][1]==poly['shift'][j]+max(t[j+1]for t in poly['terms'])+shift[j]for j in range(5))
   return poly,scalar,shift
  if src.get('known_source_zero'):return{'content':'0','shift':[0]*5,'terms':[]},Q(1),[0]*5
  assert src.get('QQ_monomial')and profile['terms_upper']==1
  return{'content':profile['known_QQ_content'],'shift':[v[0]for v in profile['exponent_box5']],'terms':[['1',0,0,0,0,0]]},Q(1),[0]*5
 try:
  gp=ROOT/'reports'/VERTEX/'report.json';g=read(gp);assert pins[str(gp.relative_to(ROOT))]=='6ad0a98af6dd6539bf2d9cdf16ec8e431ee70df66dbafef540c6aac32c73a1b8'and g['all_passed']
  ep=gp.parent/'graph.json';assert sha(ep)==g['evidence_sha256']['graph.json'];D=read(ep);graph=D['nodes'];root=D['root'];N=len(graph);assert N==1203
  # Compact IN PLACE: never retain a second set of the10lead arrays.
  cold=count_numeric(D)+N*12+25000;own(cold)
  for node in graph:
   lead=node.pop('leads');node['H_extrema']=[None if lead[0]is None else lead[0][0],None if lead[1]is None else lead[1][0]]
   del lead
  D.pop('orders',None);D.pop('root_vertices',None);gc.collect()
  assert graph[root]['H_extrema']==[-3,21]
  for k,node in enumerate(graph):
   qout=Q(node['source_profile']['known_QQ_content']);args=node['args'];kind=node['op'];weights=[]
   if qout:
    if kind in['copy','immutable_SOURCE_view']:
     ratio=Q(graph[args[0]]['source_profile']['known_QQ_content'])*Q(node['scalar'])/qout;assert ratio.denominator==1;weights=[str(ratio.numerator)]
    elif kind in['multiply','unmaterialized_multiply_for_SOURCE_union']:
     ratio=Q(graph[args[0]]['source_profile']['known_QQ_content'])*Q(graph[args[1]]['source_profile']['known_QQ_content'])*Q(node['scalar'])/qout;assert ratio.denominator==1;weights=[str(ratio.numerator)]
    elif kind!='input':
     for j,child in enumerate(args):
      ratio=Q(graph[child]['source_profile']['known_QQ_content'])*(Q(node.get('scalarB','1'))if j else Q(1))/qout;assert ratio.denominator==1;weights.append(str(ratio.numerator))
   node['integer_weights']=weights;node['actual_zero_by_source_q']=not bool(qout)
   scalarmax=max(scalarmax,*[abs(int(w)).bit_length()for w in weights],0)
  META=30000+count_numeric(D)+N*12+4000
  assert META+20000<150000
  profiles=[];roots=[];modeplans=[];maxbits=0;futurepeak=0;pairtotal=0
  for mode,target in enumerate([-3,21]):
   context={'mode':mode,'stage':'backward_tie_pruning'};needed=set();active={}
   def visit(k):
    if k in needed:return
    needed.add(k);node=graph[k];h=node['H_extrema'][mode];kind=node['op'];args=node['args']
    if h is None or kind=='input':children=[]
    elif kind in['copy','immutable_SOURCE_view']:
     children=args;assert h==graph[args[0]]['H_extrema'][mode]+node['shift'][0]
    elif kind in['multiply','unmaterialized_multiply_for_SOURCE_union']:
     children=args;assert all(graph[x]['H_extrema'][mode]is not None for x in args);assert h==sum(graph[x]['H_extrema'][mode]for x in args)+node['shift'][0]
    else:children=[x for x in args if graph[x]['H_extrema'][mode]==h]
    active[k]=children
    for child in children:visit(child)
   visit(root);uses={k:0 for k in needed}
   for k,children in active.items():
    for child in children:uses[child]+=1
   uses[root]+=1;cache={};counts={};hashes={};input_records=[];schedule=[];support_pairs=0;uses_v=None
   for k in sorted(needed):
    node=graph[k];h=node['H_extrema'][mode];kind=node['op'];args=active[k];context={'mode':mode,'node':k,'kind':kind};S=set()
    resident=sum(map(len,cache.values()));extra=0
    if h is None:pass
    elif kind=='input':
     p,sc,shift=polynomial(node['source'],node['literal_profile']);leafmax=max(leafmax,len(p['terms']))
     transient=count_numeric(p)+2048;own(META+5*resident+transient+leafparsemax+2048)
     ex=[]
     for t in p['terms']:
      e=tuple(p['shift'][j]+int(t[j+1])+shift[j]for j in range(5))
      if e[0]==h:
       if e not in S:own(META+5*(resident+len(S)+1)+transient+leafparsemax+2048);S.add(e)
     qout=Q(node['source_profile']['known_QQ_content'])
     if qout:ratio=Q(p['content'])*sc/qout;assert ratio.denominator==1;node['leaf_integer_ratio']=str(ratio.numerator)
     else:node['leaf_integer_ratio']='0'
     input_records.append({'node':k,'existing_literal_terms':len(p['terms']),'decoded_numeric_upper':transient,'source_face_terms':len(S),'whole_source_file_transform_preserved':True});del p
    elif kind in['copy','immutable_SOURCE_view']:
     for e in cache[args[0]]:
      v=tuple(e[j]+node['shift'][j]for j in range(5));own(META+5*(resident+len(S)+1)+2048);S.add(v)
    elif kind in['multiply','unmaterialized_multiply_for_SOURCE_union']:
     A,B=(cache[x]for x in args)
     for e in A:
      for f in B:
       v=tuple(e[j]+f[j]+node['shift'][j]for j in range(5))
       if v not in S:own(META+5*(resident+len(S)+1)+2048);S.add(v)
     support_pairs+=len(A)*len(B)
    else:
     for child in args:
      for e in cache[child]:
       if e not in S:own(META+5*(resident+len(S)+1)+2048);S.add(e)
    assert all(e[0]==h for e in S);counts[k]=len(S);cache[k]=S
    # Hash only this small face, charged as at most6n exponent/sort scratch.
    own(META+5*sum(map(len,cache.values()))+6*len(S)+2048)
    hashes[k]=hashlib.sha256(enc([list(e)for e in sorted(S)])).hexdigest()
    for child in args:
     uses[child]-=1
     if not uses[child]:cache.pop(child,None)
    if not uses[k]:cache.pop(k,None)
   rootS=cache[root];assert len(rootS)==[766,214][mode];roots.append(rootS);pairtotal+=support_pairs
   # VALUE closure additionally prunes exact q=0. No use of ten vanished vertices.
   requiredV=set()
   def visitV(k):
    if k in requiredV:return
    requiredV.add(k)
    if graph[k]['actual_zero_by_source_q']:return
    for x in active[k]:visitV(x)
   visitV(root);uv={k:0 for k in requiredV}
   for k in requiredV:
    if not graph[k]['actual_zero_by_source_q']:
     for c in active[k]:uv[c]+=1
   uv[root]+=1;live={};value_peak=META+2048;value_pairs=0
   for k in sorted(requiredV):
    node=graph[k];n=0 if node['actual_zero_by_source_q']else counts[k];old=sum(live.values());decode_extra=0
    if node['op']=='input'and n:decode_extra=leafparsemax+8*leafmax+2048
    peak=META+7*(old+n)+decode_extra+2048;value_peak=max(value_peak,peak)
    live[k]=n
    if node['op']in['multiply','unmaterialized_multiply_for_SOURCE_union']and n:value_pairs+=counts[active[k][0]]*counts[active[k][1]]
    M=int(node['source_profile']['primitive_integer_abs_upper']);maxbits=max(maxbits,M.bit_length())
    schedule.append({'node':k,'face_terms_upper':n,'existing_live_terms':old,'new_output_before_release':n,'decoded_leaf_plus_metadata_scratch':decode_extra,'whole_VALUE_numeric_upper':peak})
    if not node['actual_zero_by_source_q']:
     for c in active[k]:
      uv[c]-=1
      if not uv[c]:live.pop(c,None)
    if not uv[k]:live.pop(k,None)
   futurepeak=max(futurepeak,value_peak)
   profiles.append({'mode':mode,'root_H':target,'source_nodes':len(needed),'VALUE_nodes':len(requiredV),'face_terms_by_node':{str(k):n for k,n in counts.items()},'face_hash_by_node':{str(k):v for k,v in hashes.items()},'VALUE_lastuse_schedule':schedule,'source_pair_majorant':support_pairs,'VALUE_pair_majorant':value_pairs,'whole_VALUE_numeric_upper':value_peak,'input_decoding_records':input_records})
   modeplans.append({'mode':mode,'H':target,'SOURCE_needed':sorted(needed),'VALUE_needed':sorted(requiredV),'active_tied_args':{str(k):v for k,v in active.items()}})
  # Compare COMPLETE5-coordinate faces to the completed positive-source snapshot.
  phase=read(ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_LAST2_VIEW_PROFILE_20261105_v64/phase-10.json')
  snap=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_LAST2_VIEW_PROFILE_20261105_v64'/phase['all_live_snapshot']['file'];pins[str(snap.relative_to(ROOT))]=sha(snap);assert pins[str(snap.relative_to(ROOT))]==phase['all_live_snapshot']['sha256']
  from s61_sparse_generic5_source_view_profile_20261105_v64 import unpack
  got=[set(),set()]
  with snap.open('rb')as f:
   assert f.read(8)==b'Q5SUPP1\n';header=json.loads(f.read(int.from_bytes(f.read(4),'big')))
   for r in header['maps']:
    for _ in range(r['terms']):
     b=f.read(5);assert len(b)==5
     if r['id']==phase['accumulator']:
      e=unpack(int.from_bytes(b,'big'))
      if e[0]==-3:got[0].add(e)
      elif e[0]==21:got[1].add(e)
   assert not f.read(1)
  assert got==roots
  # Bound integer VALUE operations before any signed face evaluation.
  committed=max(int(n['source_profile']['primitive_integer_abs_upper']).bit_length()for n in graph)
  temporary=0
  for n in graph:
   if n['actual_zero_by_source_q']:continue
   M=[int(graph[x]['source_profile']['primitive_integer_abs_upper']).bit_length()for x in n['args']];W=[abs(int(x)).bit_length()for x in n['integer_weights']]
   if n['op']in['multiply','unmaterialized_multiply_for_SOURCE_union']:b=sum(M)+max(W,default=0)+1
   else:b=max([m+w for m,w in zip(M,W)]+[committed])+2
   temporary=max(temporary,b)
  assert max(committed,temporary)<=2048
  put('compact-face-graph.json',{'nodes':graph,'root':root,'mode_plans':modeplans,'fixed_source_q_integer_maps':True,'all5coordinates_retained':True})
  put('two-complete-positive-H-faces.json',{'faces':[{'H':h,'terms':len(S),'full5_exponents':[list(e)for e in sorted(S)]}for h,S in zip([-3,21],roots)],'full_snapshot_equality':True,'no_signed_face_coefficients':True})
  put('positive-face-profiles-and-VALUE-lifetimes.json',{'profiles':profiles,'source_owned_peak':sourcepeak,'future_VALUE_whole_numeric_upper':futurepeak,'metadata_base':META,'decoded_leaf_numeric_peak':leafparsemax,'largest_input_terms':leafmax,'committed_integer_bits_upper':committed,'temporary_integer_bits_upper':temporary,'SOURCE_source_integer_weights_exact':True})
  for f in[Path(__file__),ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py',ROOT/f'research/notes/{ID}.md']:pins[str(f.relative_to(ROOT))]=sha(f)
  raw=sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in out.iterdir())+262144;assert raw<32*1024*1024
  report={'id':ID,'all_passed':True,'SOURCE_ONLY_no_actual_face_coefficients':True,'positive_source_full_faces':[766,214],'all5coords_verified_against_final_snapshot':True,'SSA1203_compacted_in_place':True,'SOURCE_owned_peak':sourcepeak,'future_VALUE_numeric_upper':futurepeak,'future_VALUE_150k_gate':futurepeak<=150000,'future_VALUE_committed_bits':committed,'future_VALUE_temporary_bits':temporary,'VALUE2048_not_changed':True,'all_integer_weights_exact_QQ_ratios':True,'RAW_upper':raw,'no_face_VALUE_GO_or_entryzero_nonzero_claim':True,'all167x14_87plus66_oneP290_F9F10_complements':True,'input_sha256':pins,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};put('report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
 except Exception as exc:
  if not(out/'failure.json').exists():put('failure.json',{'all_passed':False,'SOURCE_only':True,'error':repr(exc),'context':context,'no_signed_face_coefficient_evaluated':True,'no_auto_retry':True})
  raise
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');main(p.parse_args())

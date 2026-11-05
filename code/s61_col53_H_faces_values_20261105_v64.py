"""Prepared exact integer full H-face coefficients; requires SOURCE and ROOT pins."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,gc,hashlib,json,resource,signal,time
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
from s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64 import write_cell,validate_cell
from s61_sparse_generic5_source_view_profile_20261105_v64 import count_numeric
ROOT=Path(__file__).resolve().parents[1]
ID='S61_COL53_H_FACES_VALUES_20261105_v64'
SOURCE='S61_COL53_H_FACES_SOURCE_20261105_v64'
OLD='reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def main(a):
 assert a.execute and a.gate_file and a.gate_sha and a.root_review_file and a.root_review_sha
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 inputs={};budget={'generated_RAW_upper':65536};context={};completed=[];cache={};work={};phasepeak=0;updates=0;last_uses=0;META=0;transient=0
 def put(name,v):
  b=enc(v);assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())+len(b)<=budget['generated_RAW_upper']
  p=out/name;assert not p.exists();p.write_bytes(b);return p
 def stop(s,f):
  checkpoint_maps()
  put('failure.json',{'all_passed':False,'context':context,'completed_faces':completed,'no_automatic_retry':True,'no_partial_face_claim':True});raise TimeoutError('VALUE55s; preserve atomic faces')
 signal.signal(signal.SIGALRM,stop);signal.alarm(55)
 def read(p):
  p=Path(p);inputs[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def own(extra=2048,newterm=False):
  nonlocal phasepeak
  n=META+7*(sum(map(len,cache.values()))+len(work)+(1 if newterm else 0))+transient+extra
  phasepeak=max(phasepeak,n);assert n<=budget['whole_VALUE_numeric_upper']<=150000,('VALUE prospective whole-owned guard',n,context)
 def commit(e,v,M):
  nonlocal updates
  assert len(e)==5 and all(isinstance(x,int)for x in e)and abs(v)<=M and abs(v).bit_length()<=640
  updates+=1
  if v:
   if e not in work:own(newterm=True)
   work[e]=v
  else:work.pop(e,None)
 def rawtemp(v):assert abs(v).bit_length()<=1921;return v
 def checkpoint_maps():
  nonlocal META
  if getattr(checkpoint_maps,'busy',False):return
  checkpoint_maps.busy=True
  try:
   if not META:return
   graph_local=graph
  except NameError:return
  reserve=budget['failure_numeric_reserve'];META-=reserve
  try:
   path=out/'failure-cache-index.jsonl';assert not path.exists()
   def maps_to_save():
    for k,A in cache.items():yield k,A
    if work and 'node'in context and context['node']not in cache:yield context['node'],work
   with path.open('wb')as index:
    for k,A in maps_to_save():
     own(extra=3*len(A)+4096)
     assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())+budget['one_cell_RAW_upper']+65536<=budget['generated_RAW_upper']
     name='failure-cache-'+str(k)+('-partial'if A is work else'')+'.bin'
     rec=write_cell(out/name,Q(graph_local[k]['source_profile']['known_QQ_content']),A)
     assert validate_cell(out/name,Q(graph_local[k]['source_profile']['known_QQ_content']),A,rec)
     index.write(enc({'node':k,'partial_node_not_complete':A is work,'file':name,**rec}))
  finally:META+=reserve

 def polynomial(src,profile):
  nonlocal transient
  transient=budget['leaf_decode_numeric_upper'];own()
  if 'binding'in src:
   b=src['binding'];z=read(ROOT/OLD/b['scope']);assert count_numeric(z)<=budget['leaf_decode_numeric_upper']
   p=(z['common_operands']if'common_operands'in z else z['row_operands'])[b['label']]
   return p['original_literal'],Q(p['QQ_transform_scalar']),p['transform_monomial_exponent5']
  if'file'in src:
   p=ROOT/src['file']
   if str(p.relative_to(ROOT))not in inputs:inputs[str(p.relative_to(ROOT))]=sha(p)
   if'offset'in src:
    with p.open('rb')as f:f.seek(src['offset']);b=f.read(src['bytes'])
    assert hashlib.sha256(b).hexdigest()==src['payload_sha256'];return decode(b),Q(1),[0]*5
   assert inputs[str(p.relative_to(ROOT))]==src['sha256'];b=p.read_bytes();assert b[:9]==b'Q5CELL1\n\x01';poly=decode(b[9:])
   lo=[poly['shift'][j]+min(t[j+1]for t in poly['terms'])for j in range(5)]
   shift=[profile['exponent_box5'][j][0]-lo[j]for j in range(5)];scalar=Q(profile['known_QQ_content'])/Q(poly['content'])
   expected=(Q(353894400000),[1,0,0,0,1])if'col56-ninth-cleared'in src['file']else(Q(1),[0]*5)
   assert(scalar,shift)==expected and all(profile['exponent_box5'][j][1]==poly['shift'][j]+max(t[j+1]for t in poly['terms'])+shift[j]for j in range(5))
   return poly,scalar,shift
  if src.get('known_source_zero'):return{'content':'0','shift':[0]*5,'terms':[]},Q(1),[0]*5
  assert src.get('QQ_monomial')
  return{'content':profile['known_QQ_content'],'shift':[b[0]for b in profile['exponent_box5']],'terms':[['1',0,0,0,0,0]]},Q(1),[0]*5
 try:
  gp=ROOT/a.gate_file;assert sha(gp)==a.gate_sha;g=read(gp);assert g['all_passed']
  bp=gp.parent/'whole-integer-face-VALUE-budget.json';assert sha(bp)==g['evidence_sha256'][bp.name];budget=read(bp)
  assert sha(Path(__file__))==budget['worker_sha256']
  rp=ROOT/a.root_review_file;assert sha(rp)==a.root_review_sha and rp.stat().st_size<=65536
  rootdoc=read(rp);assert rootdoc['all_passed']and count_numeric(rootdoc)<=2000
  for p,h in g['input_sha256'].items():assert sha(ROOT/p)==h;inputs[p]=h
  src=ROOT/'reports'/SOURCE;D=read(src/'compact-face-graph.json');P=read(src/'positive-face-profiles-and-VALUE-lifetimes.json');F=read(src/'two-complete-positive-H-faces.json')
  graph=D['nodes'];root=D['root'];META=budget['VALUE_metadata_base'];assert len(graph)==1203
  Mlimits=[int(n['source_profile']['primitive_integer_abs_upper'])for n in graph]
  assert max(m.bit_length()for m in Mlimits)<=640
  put('before-actual-face-values.json',{'ONE_bounded_face_VALUE_job':True,'source':SOURCE,'gate_sha':a.gate_sha,'root_review_sha':a.root_review_sha,'metadata_base':META,'full5coordinates_7numeric_per_term':True,'whole_numeric_upper':budget['whole_VALUE_numeric_upper'],'committed640_temp1921_within2048':True})
  for mode,target in enumerate([-3,21]):
   context={'mode':mode,'H':target,'stage':'start'};plan=D['mode_plans'][mode];required=set(plan['VALUE_needed']);active={int(k):v for k,v in plan['active_tied_args'].items()};counts={int(k):v for k,v in P['profiles'][mode]['face_terms_by_node'].items()};allowed={tuple(e)for e in F['faces'][mode]['full5_exponents']}
   uses={k:0 for k in required}
   for k in required:
    if not graph[k]['actual_zero_by_source_q']:
     for c in active[k]:uses[c]+=1
   uses[root]+=1;cache={};work={};own()
   for k in sorted(required):
    node=graph[k];h=node['H_extrema'][mode];kind=node['op'];context={'mode':mode,'H':target,'node':k,'kind':kind};work={};M=Mlimits[k];own()
    if node['actual_zero_by_source_q']or h is None:pass
    elif kind=='input':
     p,sc,shift=polynomial(node['source'],node['literal_profile']);assert count_numeric(p)<=budget['leaf_decode_numeric_upper']
     qout=Q(node['source_profile']['known_QQ_content']);ratio=Q(p['content'])*sc/qout;assert ratio.denominator==1 and ratio.numerator==int(node['leaf_integer_ratio'])
     w=ratio.numerator
     for term in p['terms']:
      e=tuple(p['shift'][j]+int(term[j+1])+shift[j]for j in range(5))
      if e[0]==h:
       assert e not in work;v=rawtemp(int(term[0])*w);commit(e,v,M)
     del p;transient=0
    elif kind in['copy','immutable_SOURCE_view']:
     w=int(node['integer_weights'][0]);child=active[k][0]
     for e,c in cache[child].items():
      v=tuple(e[j]+node['shift'][j]for j in range(5));assert v[0]==h;commit(v,rawtemp(c*w),M)
    elif kind in['multiply','unmaterialized_multiply_for_SOURCE_union']:
     left,right=active[k];w=int(node['integer_weights'][0])
     for e,c in cache[left].items():
      for f,d in cache[right].items():
       v=tuple(e[j]+f[j]+node['shift'][j]for j in range(5));assert v[0]==h
       term=rawtemp(rawtemp(c*d)*w);value=rawtemp(work.get(v,0)+term);commit(v,value,M)
    else:
     assert kind in['add','SOURCE_fused_union_product']
     for j,child in enumerate(node['args']):
      if child not in active[k]:continue
      w=int(node['integer_weights'][j])
      for e,c in cache[child].items():
       assert e[0]==h;commit(e,rawtemp(work.get(e,0)+rawtemp(c*w)),M)
    assert len(work)<=counts[k]and all(e[0]==h for e in work)
    # Move the one actual map into cache, never duplicate or repack4coordinates.
    cache[k]=work;work={};own()
    if not node['actual_zero_by_source_q']:
     for c in active[k]:
      uses[c]-=1;last_uses+=1
      if not uses[c]:cache.pop(c,None)
    if not uses[k]:cache.pop(k,None)
    assert time.monotonic()-start<50
   result=cache[root];assert set(result)<=allowed and all(e[0]==target for e in result);qroot=Q(graph[root]['source_profile']['known_QQ_content'])
   # Integer codec primitives/content and streaming readback have separate
   # prospective headroom; keep full5-coordinate keys unchanged.
   own(extra=3*len(result)+4096)
   name='H-minus3-face.bin'if mode==0 else'H-plus21-face.bin';path=out/name
   assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())+budget['face_binary_RAW_upper'][mode]+65536<=budget['generated_RAW_upper']
   rec=write_cell(path,qroot,result);assert validate_cell(path,qroot,result,rec)and rec['bytes']<=budget['face_binary_RAW_upper'][mode]
   completed.append({'H':target,'positive_source_terms':len(allowed),'actual_nonzero_coefficients':len(result),'file':name,**rec,'all_values_exact_integer_FIXED_source_q':True,'source_q':str(qroot),'full5_coordinate_stream_roundtrip':True})
   put('completed-face-'+str(mode)+'.json',completed[-1]);cache={};work={};del result;gc.collect()
  raw=sum((ROOT/p).stat().st_size for p in inputs)+sum(p.stat().st_size for p in out.iterdir()if p.is_file())+262144
  assert raw<=budget['fresh_RAW_upper']
  report={'id':ID,'all_passed':True,'two_actual_full_H_faces_evaluated':True,'faces':completed,'nonzero_face_count':sum(bool(v['actual_nonzero_coefficients'])for v in completed),'entry_nonzero_proved':any(v['actual_nonzero_coefficients']for v in completed),'both_zero_does_NOT_prove_entry_zero':True,'full5_coords_VALUE7term_charge':True,'actual_whole_owned_peak':phasepeak,'whole_numeric_upper':budget['whole_VALUE_numeric_upper'],'integer_updates':updates,'last_use_releases':last_uses,'committed640_temporary1921_within_VALUE2048':True,'actual_READ_output_RAW_upper':raw,'scope':'fixedrow6/col53 numeratoronly; full167x14/87+66/oneP290/F9F10/allcomplements retained; no all14/fullinverse/exhaustion','input_sha256':inputs,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};put('report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256','faces')}));print('report_sha256',sha(out/'report.json'))
 except Exception as exc:
  if not(out/'failure.json').exists():
   checkpoint_maps()
   put('failure.json',{'all_passed':False,'error':repr(exc),'context':context,'completed_faces':completed,'no_automatic_retry':True,'no_partial_face_claim':True})
  raise
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--execute',action='store_true');p.add_argument('--gate-file');p.add_argument('--gate-sha');p.add_argument('--root-review-file');p.add_argument('--root-review-sha');main(p.parse_args())

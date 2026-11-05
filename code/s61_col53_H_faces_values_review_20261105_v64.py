#!/usr/bin/env python3
"""Independent reversed signed integer FULL2 H-face source replay/readback."""
import hashlib,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
from s61_col53_vertex_source_values_review_20261105_v64 import decoder,literal_record
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='S61_COL53_H_FACES_VALUES_20261105_v64';SRC='S61_COL53_H_FACES_SOURCE_20261105_v64';GATE='S61_COL53_H_FACES_VALUE_GATE_20261105_v64'
ID='S61_COL53_H_FACES_VALUES_20261105_v64_REVIEW'
OLD='reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main():
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 def timeout(s,f):raise TimeoutError('Independentfacepeer55s; preserveprefix')
 signal.signal(signal.SIGALRM,timeout);signal.alarm(55)
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};RAW=0
 def pin(p,h=None,parse=True):
  nonlocal RAW
  p=pathlib.Path(p);k=str(p.relative_to(ROOT))
  if k not in pins:pins[k]=sha(p);RAW+=p.stat().st_size
  if h:assert pins[k]==h,(k,pins[k],h)
  return json.loads(p.read_text()) if parse else None
 pr=pin(ROOT/'reports'/PRI/'report.json','b8d4309eaddecdd9b362934f7d2f8bde1973487d071ca507b730f426db8fd983')
 for p,h in pr['input_sha256'].items():pin(ROOT/p,h,False)
 sr=pin(ROOT/'reports'/SRC/'report.json','5e66838fb845def468a9a1c883fd08a5836f21fdfd3a88de8f4d62109013d718')
 gr=pin(ROOT/'reports'/GATE/'report.json','fe261a87001601168e3f52eae66fe4bd5e1e80ddda810836de7b656736ac09fe')
 peer=pin(ROOT/'reports/S61_COL53_H_FACES_SOURCE_CODE_20261105_v64_REVIEW/report.json','7c6e258796e3a5cd5f4b6a5d612c56c2def02b4ee03fd587fbefa11a7ff13783');assert peer['all_passed']
 def doc(st,n,r):return pin(ROOT/'reports'/st/n,r['evidence_sha256'][n])
 D=doc(SRC,'compact-face-graph.json',sr);P=doc(SRC,'positive-face-profiles-and-VALUE-lifetimes.json',sr);F=doc(SRC,'two-complete-positive-H-faces.json',sr)
 B=doc(GATE,'whole-integer-face-VALUE-budget.json',gr);G=D['nodes'];root=D['root'];assert len(G)==1203 and root==1202
 metadata=B['VALUE_metadata_base']+4096;normal_limit=B['whole_VALUE_numeric_upper']+4096
 assert metadata==102267 and normal_limit==130358<150000
 inputs=[n for n in G if n['op']=='input']
 bmax=max(n['literal_profile']['terms_upper'] for n in inputs if 'binding' in n['source'])
 fmax=max(n['literal_profile']['terms_upper'] for n in inputs if 'file' in n['source'])
 assert max(14*bmax+1024,6*fmax+128)<=16152
 assert B['committed_SOURCE_integer_bits']==401 and B['normal_commit_bits']==640 and B['normal_temp_bits']==1921
 save(out/'positive-before-independent-signed-twoface-replay.json',{'SOURCE_peer':pins['reports/S61_COL53_H_FACES_SOURCE_CODE_20261105_v64_REVIEW/report.json'],
  'same_fixed_q_integer_maps_and_source_majorants':True,'reverse_order_integer_arithmetic':True,'metadata_with4096peer_extra':metadata,'normal_owned_forecast':normal_limit,
  'leaf_decode_upper':16152,'strict5coordinate_term_charge':7,'commit_limit':640,'temporary_limit':1921,'VALUE_hard':2048,'hard_owned':150000,
  'no_new_H_or_parameterpoint_or_prime':True,'all411_expected_values_NOT_used_for_pruning':True})
 cache={};work={};transient=0;peak=metadata;tempbits=0;stored=0;updates=0;completed=[];context={}
 def own(extra=2048,new=False):
  nonlocal peak
  used=metadata+7*(sum(map(len,cache.values()))+len(work)+int(new))+transient+extra
  peak=max(peak,used);assert used<=normal_limit,(used,context)
 def mul(a,b):
  nonlocal tempbits
  h=abs(a).bit_length()+abs(b).bit_length();tempbits=max(tempbits,h);assert h<=1921
  v=a*b;assert abs(v).bit_length()<=1921;return v
 def add(a,b):
  nonlocal tempbits
  h=max(abs(a).bit_length(),abs(b).bit_length())+1;tempbits=max(tempbits,h);assert h<=1921
  return a+b
 def put(e,v,M):
  nonlocal stored,updates
  assert len(e)==5 and all(isinstance(x,int) for x in e) and abs(v)<=M and abs(v).bit_length()<=640
  stored=max(stored,abs(v).bit_length());updates+=1
  if v:
   if e not in work:own(new=True)
   work[e]=v
  else:work.pop(e,None)
 def leaf(n):
  nonlocal transient
  transient=16152;own();s=n['source'];p=n['literal_profile']
  if 'binding' in s:
   b=s['binding'];path=ROOT/OLD/b['scope'];pin(path,pr['input_sha256'][str(path.relative_to(ROOT))],False)
   rec=literal_record(path,b['label']);poly=rec['original_literal'];q=Q(poly['content']);off=poly['shift'];terms=[(int(t[0]),tuple(t[1:])) for t in poly['terms']]
   scale=Q(rec['QQ_transform_scalar']);shift=rec['transform_monomial_exponent5']
  else:
   path=ROOT/s['file'];pin(path,pr['input_sha256'][s['file']],False)
   with path.open('rb') as f:
    if 'offset' in s:f.seek(s['offset']);raw=f.read(s['bytes']);assert hashlib.sha256(raw).hexdigest()==s['payload_sha256']
    else:raw=f.read();assert raw[:9]==b'Q5CELL1\n\x01';raw=raw[9:]
   q,off,terms=decoder(raw);scale=Q(1);shift=[0]*5
   if 'offset' not in s:
    low=[off[j]+min(e[j] for c,e in terms) for j in range(5)];shift=[p['exponent_box5'][j][0]-low[j] for j in range(5)]
    scale=Q(p['known_QQ_content'])/q
    assert (scale,shift)==((Q(353894400000),[1,0,0,0,1]) if 'col56-ninth-cleared' in s['file'] else (Q(1),[0]*5))
  ratio=q*scale/Q(n['source_profile']['known_QQ_content'])
  assert ratio.denominator==1 and str(ratio.numerator)==n['leaf_integer_ratio']
  return off,terms,shift,ratio.numerator
 try:
  for mode,target in enumerate([-3,21]):
   plan=D['mode_plans'][mode];req=set(plan['VALUE_needed']);active={int(k):v for k,v in plan['active_tied_args'].items()};counts={int(k):v for k,v in P['profiles'][mode]['face_terms_by_node'].items()}
   uses={k:0 for k in req}
   for k in req:
    if not G[k]['actual_zero_by_source_q']:
     for child in active[k]:uses[child]+=1
   uses[root]+=1;cache={};work={};own()
   for k in sorted(req):
    n=G[k];kind=n['op'];h=n['H_extrema'][mode];M=int(n['source_profile']['primitive_integer_abs_upper']);context={'mode':mode,'node':k,'kind':kind};work={};own()
    if n['actual_zero_by_source_q'] or h is None:pass
    elif kind=='input':
     off,terms,shift,w=leaf(n)
     for c,e in reversed(terms):
      key=tuple(off[j]+e[j]+shift[j] for j in range(5))
      if key[0]==h:assert key not in work;put(key,mul(c,w),M)
     del terms;transient=0
    elif kind in ['copy','immutable_SOURCE_view']:
     child=cache[active[k][0]];w=int(n['integer_weights'][0])
     for e in reversed(child):put(tuple(e[j]+n['shift'][j] for j in range(5)),mul(child[e],w),M)
    elif kind in ['multiply','unmaterialized_multiply_for_SOURCE_union']:
     left,right=[cache[j] for j in active[k]];w=int(n['integer_weights'][0])
     for e in reversed(left):
      for f in reversed(right):
       key=tuple(e[j]+f[j]+n['shift'][j] for j in range(5));term=mul(mul(left[e],right[f]),w);put(key,add(work.get(key,0),term),M)
    else:
     assert kind in ['add','SOURCE_fused_union_product']
     for j in reversed(range(len(n['args']))):
      child=n['args'][j]
      if child not in active[k]:continue
      w=int(n['integer_weights'][j]);a=cache[child]
      for e in reversed(a):put(e,add(work.get(e,0),mul(a[e],w)),M)
    assert len(work)<=counts[k] and all(e[0]==h for e in work)
    cache[k]=work;work={};own()
    if not n['actual_zero_by_source_q']:
     for child in active[k]:
      uses[child]-=1
      if not uses[child]:cache.pop(child,None)
    if not uses[k]:cache.pop(k,None)
    assert time.monotonic()-start<50
   result=cache[root];qroot=Q(G[root]['source_profile']['known_QQ_content'])
   assert set(result)<={tuple(e) for e in F['faces'][mode]['full5_exponents']} and all(e[0]==target for e in result)
   # OnlyNOW read the saved actual cell; canonical signed integer support
   # is checked term-by-term with this independently replayed full result.
   rec=pr['faces'][mode];path=ROOT/'reports'/PRI/rec['file'];pin(path,pr['evidence_sha256'][rec['file']],False);raw=path.read_bytes()
   assert raw[:8]==b'Q5CELL1\n'
   if not result:
    assert raw==b'Q5CELL1\n\x00' and rec['terms']==0
    content=Q(0);gcd=1;encoded_terms=0;digest=hashlib.sha256(b'0\n').hexdigest();first=None
   else:
    assert raw[8]==1;content,shift,terms=decoder(raw[9:]);own(extra=6*len(terms)+2048)
    assert len(terms)==len(result)==rec['terms'];ratio=content/qroot;assert ratio.denominator==1 and ratio.numerator>0
    gcd=0
    for c in result.values():gcd=math.gcd(gcd,abs(c))
    assert ratio.numerator==gcd;previous=None;primitive_gcd=0;digest_obj=hashlib.sha256((str(len(terms))+'\n').encode())
    for c,off in terms:
     e=tuple(shift[j]+off[j] for j in range(5));key=(sum(e),e)
     assert previous is None or previous>key;previous=key
     assert result[e]==c*gcd;primitive_gcd=math.gcd(primitive_gcd,abs(c))
     qq=qroot*result[e];assert max(abs(qq.numerator).bit_length(),qq.denominator.bit_length())<=2048
     digest_obj.update((str(qq)+' '+' '.join(map(str,e))+'\n').encode())
    assert primitive_gcd==1 and all(min(off[j] for c,off in terms)==0 for j in range(5))
    digest=digest_obj.hexdigest();first=[list(shift[j]+terms[0][1][j] for j in range(5)),str(content*terms[0][0])];encoded_terms=len(terms);del terms
   assert digest==rec['value_sha256'] and str(content)==rec['content']
   witness=(-3,0,6,-1,10)
   if mode==0:assert qroot*result[witness]==Q(241815,134217728)
   completed.append({'H':target,'source_positive_terms':[766,214][mode],'independent_actual_nonzero_terms':len(result),'canonical_content':str(content),'canonical_gcd':str(gcd),'full_cell_value_hash':digest,'first_canonical_term':first,'all5coordinates_readback':True})
   save(out/f'completed-independent-face-{mode}.json',completed[-1]);cache={};work={};del result
  assert [x['independent_actual_nonzero_terms'] for x in completed]==[411,0]
  assert all(sha(ROOT/p)==h for p,h in pins.items()) and RAW+1000000<32<<20
  proof={'all_original_signed_source_operations_replayed_independently':True,'reverse_integer_input_product_and_sum_order':True,
   'faces':completed,'all411_actual_nonzero_coefficients_and_zero_face_match':True,
   'explicit_exact_entry_nonzero_witness':{'exponent5':[-3,0,6,-1,10],'coefficient':'241815/134217728'},
   'same9byte_input_cell_and_scalar_monomial_transforms':True,'no_vertex_zero_or_saved_face_zero_pruning':True,
   'entry_row6_col53_cleared_numerator_nonzero':True,
   'scope':'ONLYH=-3andH=21fullfaces. H=-3nonzeroprovesfixedentrynonzero. NootherH/fullentrymap/pivot/rank/fullinverse/167x14exhaustion/parameterpointorprime computed.',
   'committed640_temporary1921_VALUE2048_unchanged':True,'independent_metadata_base':metadata,'positive_normal_owned':normal_limit,'actual_owned_peak':peak,
   'actual_integer_preoperation_bits':tempbits,'actual_stored_integer_bits':stored,'integer_updates':updates,'unique_RAW_inputs':RAW}
  ep=out/'independent-reversed-signed-twoface-source-replay-and-all-cell-values.json';save(ep,proof)
  report={'id':ID,'all_passed':True,'frozen':True,'independent_actual_signed_SOURCE_VALUE':True,'input_sha256':pins,'code_sha256':sha(__file__),
   'evidence_sha256':{ep.name:sha(ep)},'all411_minH_coefficients_and_zero_maxH_face_PASS':True,'fixed_entry_nonzero':True,'new_rank_claim':False,'seconds':time.monotonic()-start}
  save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned_peak':peak,'prebits':tempbits,'seconds':report['seconds']}))
 except BaseException as exc:
  with (out/'failure-live-prefix.jsonl').open('w') as f:
   for k,A in cache.items():
    for e,c in A.items():f.write(json.dumps([k,list(e),str(c)],separators=(',',':'))+'\n')
   if context.get('node') not in cache:
    for e,c in work.items():f.write(json.dumps(['partial',context.get('node'),list(e),str(c)],separators=(',',':'))+'\n')
  save(out/'failure.json',{'error':repr(exc),'context':context,'completed':completed,'no_automatic_retry':True,'partial_not_face_claim':True});raise
if __name__=='__main__':main()

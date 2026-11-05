#!/usr/bin/env python3
"""Independent original-SSA/signed-Laurent-vertex SOURCE and coefficient peer."""
import hashlib,itertools,json,math,pathlib,re,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1]
SRC='S61_COL53_VERTEX_SOURCE_GATE_V2_20261105_v64';VAL='S61_COL53_VERTEX_VALUES_V2_20261105_v64'
ID='S61_COL53_VERTEX_SOURCE_VALUES_V2_20261105_v64_REVIEW'
OLD='reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'
ORDERS=[(tuple((i+j)%5 for j in range(5)),sgn) for i in range(5) for sgn in [1,-1]]
PRE=0
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def ratio(s):
 q=Q(s);assert max(abs(q.numerator).bit_length(),q.denominator.bit_length())<=640;return q
def safe(a,b,multiply=False):
 global PRE
 if a and b:
  na,da=abs(a.numerator).bit_length(),a.denominator.bit_length();nb,db=abs(b.numerator).bit_length(),b.denominator.bit_length()
  bound=max(na+nb,da+db) if multiply else max(max(na+db,nb+da)+1,da+db)
  PRE=max(PRE,bound);assert bound<=1921<2048
 return a*b if multiply else a+b
def decoder(data):
 """Own canonical unsigned-varint/signed-length integer decoder."""
 pos=0
 def uint():
  nonlocal pos
  start=pos;v=0;power=0
  while True:
   assert pos<len(data);x=data[pos];pos+=1;v+=(x&127)<<power
   if x<128:break
   power+=7
  assert pos-start==max(1,(v.bit_length()+6)//7)
  return v
 def sint():
  nonlocal pos
  t=uint();n=t//2;assert pos+n<=len(data);v=int.from_bytes(data[pos:pos+n],'big')
  assert not n or data[pos]!=0;pos+=n;assert v or not t%2
  return -v if t%2 else v
 num,den=sint(),sint();assert num and den>0 and math.gcd(num,den)==1
 shifts=[]
 for _ in range(5):
  x=uint();shifts.append(-(x//2)-1 if x%2 else x//2)
 n=uint();assert n>0;terms=[]
 for _ in range(n):
  c=sint();assert c;terms.append((c,tuple(uint() for _ in range(5))))
 assert pos==len(data)
 return Q(num,den),tuple(shifts),terms
def literal_record(path,label):
 # Exact named binding; only that literal record is decoded from JSON.
 text=path.read_text();hits=list(re.finditer(r'"'+re.escape(label)+r'"\s*:',text));assert len(hits)==1
 at=hits[0].end()
 while text[at].isspace():at+=1
 return json.JSONDecoder().raw_decode(text,at)[0]
def orderkey(m,k):
 p,s=ORDERS[k];return tuple(s*m[i] for i in p)
def main():
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};RAW=0
 def pin(path,h=None,parse=True):
  nonlocal RAW
  p=pathlib.Path(path);got=sha(p)
  if h:assert got==h,(str(p),got,h)
  name=str(p.relative_to(ROOT))
  if name not in pins:RAW+=p.stat().st_size
  pins[name]=got
  return json.loads(p.read_text()) if parse else None
 sr=pin(ROOT/'reports'/SRC/'report.json','6ad0a98af6dd6539bf2d9cdf16ec8e431ee70df66dbafef540c6aac32c73a1b8')
 vr=pin(ROOT/'reports'/VAL/'report.json','ef9cb1c8175ce59df1cd95274bcc92df719036a61b4d7d41d080237e4b25ea20')
 for p,h in sr['input_sha256'].items():pin(ROOT/p,h,False)
 graphdoc=pin(ROOT/'reports'/SRC/'graph.json',sr['evidence_sha256']['graph.json']);G=graphdoc['nodes']
 pin(ROOT/'code/s61_col53_vertex_coefficients_v2_20261105_v64.py',sr['code_sha256'],False)
 assert vr['code_sha256']==sr['code_sha256'] and len(G)==1203 and graphdoc['root']==1202
 assert graphdoc['orders']==[[list(p),s] for p,s in ORDERS]
 # Original eight phases +45+48 literal receipts, independently SSA-renamed.
 phasebase=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_ROW6_COEFFICIENT_PROFILE_V2_20261105_v64'
 records=[phasebase/f'phase-{i:02d}.json' for i in range(8)]
 for stem,cs in [('S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_REMAINING3_VIEW_PROFILE_20261105_v64',[8]),('S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_LAST2_VIEW_PROFILE_20261105_v64',[9,10])]:
  for c in cs:records+=sorted((ROOT/'reports'/stem).glob(f'term-{c}-*.json'))
 assert len(records)==101
 full=[];names={};receipt_peak=0;raw_versions=[]
 for path in records:
  rec=pin(path,sr['input_sha256'][str(path.relative_to(ROOT))]);receipt_peak=max(receipt_peak,count(rec))
  ops=rec.get('ops',rec.get('literal_ops'));assert ops is not None
  assert len({x['output'] for x in ops})==len(ops)
  for op in ops:
   name=op['output'];kind=op['op'];prof=rec['nodes'][name]
   n={'name':name,'op':kind,'source_profile':{k:prof[k] for k in ['known_QQ_content','primitive_integer_abs_upper']},'args':[],
      'scalar':op.get('scalar','1'),'shift':op.get('shift5',[0]*5)}
   if kind=='input':
    n['source']=op['source'];n['literal_profile']={k:prof[k] for k in ['terms_upper','known_QQ_content','primitive_integer_abs_upper','exponent_box5']}
   elif kind in ['copy','immutable_SOURCE_view']:n['args']=[names[op['A']]]
   elif kind in ['multiply','unmaterialized_multiply_for_SOURCE_union']:n['args']=[names[op['A']],names[op['B']]]
   elif kind=='add':n['args']=[names[op['A']],names[op['B']]];n['scalarB']=op.get('scalarB','1')
   else:
    assert kind=='SOURCE_fused_union_product'
    dot=names[op['literal_dot']];assert full[dot]['op']=='unmaterialized_multiply_for_SOURCE_union'
    assert full[dot]['args']==[names[op['A']],names[op['B']]]
    n['args']=([names[op['C']]] if op['C'] is not None else [])+[dot];n['scalarB']='1'
   if name.startswith('raw53'):raw_versions.append([len(full),name])
   names[name]=len(full);full.append(n)
 endpoint=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_LAST2_VIEW_PROFILE_20261105_v64/phase-10.json'
 end=pin(endpoint,sr['input_sha256'][str(endpoint.relative_to(ROOT))]);root=names[end['accumulator']]
 required=set();todo=[root]
 while todo:
  i=todo.pop()
  if i not in required:required.add(i);todo.extend(full[i]['args'])
 ren={old:new for new,old in enumerate(sorted(required))}
 own=[{**full[i],'args':[ren[j] for j in full[i]['args']]} for i in sorted(required)]
 assert len(own)==1203 and ren[root]==1202
 for a,b in zip(own,G):assert a=={k:v for k,v in b.items() if k!='leads'}
 del full,required,ren,rec,ops
 # Own original source parsing, known scalar/monomial transforms, and input
 # vertices. Store literal source coefficients only; target values wait.
 seeds={};input_proof=[];poly_peak=0;leaf_height=0
 for i,n in enumerate(own):
  if n['op']!='input':continue
  src=n['source'];profile=n['literal_profile'];sc=Q(1);extra=(0,)*5
  if 'binding' in src:
   b=src['binding'];path=ROOT/OLD/b['scope'];pin(path,sr['input_sha256'][str(path.relative_to(ROOT))],False)
   record=literal_record(path,b['label']);poly=record['original_literal']
   content=ratio(poly['content']);shift=tuple(poly['shift']);terms=[(int(t[0]),tuple(map(int,t[1:]))) for t in poly['terms']]
   sc=ratio(record['QQ_transform_scalar']);extra=tuple(record['transform_monomial_exponent5']);poly_peak=max(poly_peak,count(record))
  else:
   path=ROOT/src['file'];pin(path,sr['input_sha256'][src['file']],False)
   with path.open('rb') as f:
    if 'offset' in src:f.seek(src['offset']);raw=f.read(src['bytes']);assert hashlib.sha256(raw).hexdigest()==src['payload_sha256']
    else:raw=f.read();assert sha(path)==src['sha256'];assert raw[:8]==b'Q5CELL1\n' and raw[8]==1;raw=raw[9:]
   content,shift,terms=decoder(raw);poly_peak=max(poly_peak,count(terms)+count(shift)+2)
   if 'offset' not in src:
    low=[shift[j]+min(e[j] for c,e in terms) for j in range(5)]
    extra=tuple(profile['exponent_box5'][j][0]-low[j] for j in range(5));sc=ratio(profile['known_QQ_content'])/content
    expected=(Q(353894400000),(1,0,0,0,1)) if 'col56-ninth-cleared' in src['file'] else (Q(1),(0,)*5)
    assert (sc,extra)==expected
  actual=[(tuple(shift[j]+e[j]+extra[j] for j in range(5)),c) for c,e in terms]
  assert len({e for e,c in actual})==len(actual)
  assert all(c for e,c in actual)
  if actual:
   lo=[min(e[j] for e,c in actual) for j in range(5)];hi=[max(e[j] for e,c in actual) for j in range(5)]
   assert [[a,b] for a,b in zip(lo,hi)]==profile['exponent_box5']
   vertices=[min((e for e,c in actual),key=lambda e:orderkey(e,k)) for k in range(10)]
   coeffs=dict(actual);rawcoeff=[coeffs[e] for e in vertices]
  else:vertices=[None]*10;rawcoeff=[0]*10
  assert [list(e) if e is not None else None for e in vertices]==G[i]['leads']
  seeds[i]=(vertices,rawcoeff,content,sc)
  leaf_height=max(leaf_height,max(abs(content.numerator).bit_length(),content.denominator.bit_length())+max(abs(sc.numerator).bit_length(),sc.denominator.bit_length())+max([abs(c).bit_length() for c in rawcoeff] or [0]))
  input_proof.append([i,n['name'],len(terms),str(sc),list(extra)])
  del terms,actual
 # All internal positive vertices, before signed target coefficient work.
 for direction in range(10):
  leads=[]
  for i,n in enumerate(own):
   kind=n['op'];args=n['args'];shift=n['shift']
   if kind=='input':e=seeds[i][0][direction]
   elif kind in ['copy','immutable_SOURCE_view']:
    p=leads[args[0]];e=None if p is None else tuple(x+y for x,y in zip(p,shift))
   elif kind in ['multiply','unmaterialized_multiply_for_SOURCE_union']:
    a,b=[leads[j] for j in args];e=None if a is None or b is None else tuple(a[k]+b[k]+shift[k] for k in range(5))
   else:
    assert shift==[0]*5 and ratio(n['scalar'])==1
    poss=[leads[j] for j in args if leads[j] is not None];e=min(poss,key=lambda p:orderkey(p,direction)) if poss else None
   assert e==(None if G[i]['leads'][direction] is None else tuple(G[i]['leads'][direction]));leads.append(e)
 # Own streamed positive snapshot parser; no support sets or producer import.
 snap=endpoint.parent/end['all_live_snapshot']['file'];pin(snap,sr['input_sha256'][str(snap.relative_to(ROOT))],False)
 best=[None]*10;snapshot_count=0
 with snap.open('rb') as f:
  assert f.read(8)==b'Q5SUPP1\n';head=json.loads(f.read(int.from_bytes(f.read(4),'big')))
  for rec in head['maps']:
   for _ in range(rec['terms']):
    raw=f.read(5);assert len(raw)==5
    if rec['id']==end['accumulator']:
     key=int.from_bytes(raw,'big');e=tuple(((key>>(8*j))&255)-64 for j in range(5));snapshot_count+=1
     for k in range(10):
      if best[k] is None or orderkey(e,k)<orderkey(best[k],k):best[k]=e
  assert not f.read(1)
 assert [list(e) for e in best]==graphdoc['root_vertices']
 # Positive profiles bound EACH source coefficient; no unknown entry values.
 heights=[max(abs(ratio(n['source_profile']['known_QQ_content']).numerator).bit_length()+int(n['source_profile']['primitive_integer_abs_upper']).bit_length(),ratio(n['source_profile']['known_QQ_content']).denominator.bit_length()) for n in own]
 assert max(heights)==sr['coefficient_prior_bits']==402 and leaf_height<=640
 owned=count([sr,vr,graphdoc,own,seeds,input_proof,pins,head,end])+receipt_peak+2*poly_peak+20000
 assert owned<150000 and RAW+1000000<32<<20
 prior={'SOURCE_PASS_before_independent_target_values':True,'all1203_SSA_nodes_and_original101_receipts':True,
  'all91_literal_bindings_and_transforms':True,'all12030_positive_vertex_coordinates_checked':True,'all_snapshot_root_keys':snapshot_count,
  'coefficient_positive_prior':402,'committed_coefficient_limit':640,'committed_all_partial_temporary_limit':1921,'hard_bits':2048,
  'largest_input_three_factor_bound':leaf_height,'whole_owned_including_decoded_list_labels':owned,'RAW_inputs_plus_output_reserve':RAW+1000000,
  'tuple_list_accounting':'Independent count includes ALLdecoded exponentlabels, rather than relying on pre-JSONtuplecount. The resulting positive wholebound remainsbelow150k.'}
 save(out/'SOURCE-and-positive-budget-before-values.json',prior)
 # Independent one-direction-at-a-time projection, never full polynomial
 # expansion and never all1203x10 QQ rows retained simultaneously.
 result=[];valuehash=hashlib.sha256()
 for direction in range(10):
  vals=[]
  for i,n in enumerate(own):
   kind=n['op'];args=n['args']
   if kind=='input':
    vv=seeds[i];q=safe(safe(vv[2],vv[3],True),Q(vv[1][direction]),True)
   elif kind in ['copy','immutable_SOURCE_view']:q=safe(vals[args[0]],ratio(n['scalar']),True)
   elif kind in ['multiply','unmaterialized_multiply_for_SOURCE_union']:q=safe(safe(vals[args[0]],vals[args[1]],True),ratio(n['scalar']),True)
   else:
    q=Q(0);target=G[i]['leads'][direction]
    for j,parent in enumerate(args):
     if G[parent]['leads'][direction]==target:
      c=ratio(n['scalarB']) if j else Q(1);q=safe(q,safe(vals[parent],c,True))
   assert max(abs(q.numerator).bit_length(),q.denominator.bit_length())<=640
   vals.append(q);valuehash.update(str(q).encode()+b'\n')
  result.append({'monomial':graphdoc['root_vertices'][direction],'coefficient':str(vals[-1]),'nonzero':bool(vals[-1])})
 observed=pin(ROOT/'reports'/VAL/'ten-actual-vertex-coefficients.json',vr['evidence_sha256']['ten-actual-vertex-coefficients.json'])
 assert result==observed and all(x['coefficient']=='0' for x in result)
 # Preserve the first reader failure; it cannot authorize a broad replay.
 hist=ROOT/'reports/S61_COL53_VERTEX_SOURCE_GATE_20261105_v64_FIRST_CELL_HEADER_READER'
 failure=pin(hist/'failure.json');pin(hist/'failed-source.py',None,False);pin(hist/'failure.log',None,False)
 assert failure['source_reader_failure_only'] and failure['no_new_col53_coefficient_values']
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'all_source_and_value_checks_PASS':True,'ten_actual_coefficients':result,'all_signed_values_zero':True,
  'projection_theorem':[
   'Each fixed signedcycliclex is a translation-invariant totalorder onZ5. A finite positive support has a minimum. All monomials ofa child are at leastthat childminimum, evenwhen its actualcoefficient cancels tozero.',
   'For a product, the coefficient at the sumofsource minima receives onlythe minimum/minimum pair: twononnegative orderedgroup differences cannot sumtozero unlessbothzero. It equals the productofthose two coefficients, times the literalQQscalar. At a sumunionminimum, onlychildren withthat samepositive minimumcontribute.',
   'Induction remainsvalid if an intermediate extremal coefficient cancels. We retain its ORIGINALpositive-support exponent, do not replace it byan actualnext vertex, and do not inferthat the wholechildpolynomial vanishes. Shifts translateexponents; QQscalars altercoefficients only.',
   'All ten coordinates are therefore exactordinary Laurent coefficient functionals ofthe SAMErow6/col53 clearednumerator in the originalsourceframe. Repeatedraw53names were distinctSSAassignments, not one mutableinput.'
  ],
  'all_original_SSA_nodes':1203,'source_receipts':101,'input_transforms':input_proof,
  'all_snapshot_keys_streamed':snapshot_count,'all_positive_vertex_coordinates':12030,
  'independent_value_stream_hash':valuehash.hexdigest(),'maximum_actual_before_operation_bits':PRE,
  'source_positive_budget':prior,'NO_full_entry_zero_nonzero_rank_pivot_or_exhaustion_claim':True,
  'retained_scope':'Full167x14,87additionalRHS,66constraints,oneP290,F9/F10andALLcomplements retained. Onlythe ten fixedcoefficient values arezero. No fullentry, face, point, prime, inverse orranksearch was computed.'}
 ep=out/'independent-original-SSA-and-ten-signed-vertex-coefficients.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'independent_SOURCE_and_VALUE':True,'input_sha256':pins,'code_sha256':sha(__file__),
  'evidence_sha256':{ep.name:sha(ep)},'all10_exact_vertex_coefficients_zero':True,'whole_entry_zero_claim':False,'seconds':time.monotonic()-started}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned,'prebits':PRE,'seconds':report['seconds']}))
if __name__=='__main__':main()

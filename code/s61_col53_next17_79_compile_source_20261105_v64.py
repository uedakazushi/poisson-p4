"""Compile full source bridge for ONE minor projection. No supports or signed values."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd,prod
import argparse,hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
from s61_sparse_generic5_source_view_profile_20261105_v64 import count_numeric
ID='S61_COL53_NEXT17_79_COMPILED_SOURCE_20261105_v64'
PREF='reports/S61_COL53_NEXT17_79_FUNCTIONAL_SOURCE_20261105_v64'
OLD='reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'
SIDE='reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_COEFFICIENT_SIDE_SOURCE_20261105_v64'
def enc(x):return(json.dumps(x,separators=(',',':'))+'\n').encode()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def main(a):
 assert a.prepare
 resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));start=time.monotonic();out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 pins={};docs={};context={};peak=0
 def put(n,x):
  b=enc(x);assert sum(p.stat().st_size for p in out.iterdir())+len(b)<=4*1024*1024
  (out/n).write_bytes(b)
 def load(f):
  if f not in docs:
   p=ROOT/f;pins[f]=sha(p);docs[f]=json.loads(p.read_text())
  return docs[f]
 def own(extra=12000):
  nonlocal peak
  index_numeric=len(names)+sum(9+len(k[1])for k in dedup)+len(leafcache)+8*len(raw)
  n=25000+sum(count_numeric(d)for d in docs.values())+node_numeric+index_numeric+extra
  peak=max(peak,n);assert n<=150000,('prospective_SOURCE_compile_owned',n,context)
  assert time.monotonic()-start<50,('SOURCE_compile_time',context)
 def qm(a,b):
  assert max(abs(a.numerator).bit_length()+abs(b.numerator).bit_length(),a.denominator.bit_length()+b.denominator.bit_length())<=2048
  return a*b
 def qg(a,b):
  den=a.denominator*b.denominator//gcd(a.denominator,b.denominator);assert den.bit_length()<=2048
  aa=a.numerator*(den//a.denominator);bb=b.numerator*(den//b.denominator);return Q(gcd(abs(aa),abs(bb)),den)
 def primitive(p,scale=Q(1),shift=(0,0,0,0,0)):
  q=qm(Q(p['content']),scale)
  if not p['terms']or not q:return q,0,None,None,0
  lo=[10**6]*5;hi=[-10**6]*5;lead=None;M=0
  for term in p['terms']:
   e=tuple(p['shift'][j]+int(term[j+1])+shift[j]for j in range(5));lead=e if lead is None else min(lead,e)
   M=max(M,abs(int(term[0])))
   for j in range(5):lo[j]=min(lo[j],e[j]);hi[j]=max(hi[j],e[j])
  return q,M,lead,[[lo[j],hi[j]]for j in range(5)],len(p['terms'])
 nodes=[];node_numeric=0;names={};dedup={};oldmap={};leafcache={};raw={}
 def add_record(name,op,args,q,M,lead,box,n,scalar='1',shift=(0,0,0,0,0),source=None,old_origin=None):
  nonlocal node_numeric
  assert all(0<=j<len(nodes)for j in args)
  assert max(abs(q.numerator).bit_length(),q.denominator.bit_length(),int(M).bit_length(),int(n).bit_length())<=2048
  if box:assert all(-64<=lo<=hi<=191 for lo,hi in box)
  rec={'name':name,'op':op,'args':list(args),'q':str(q),'M':str(M),'lead5':list(lead)if lead is not None else None,'box5':box,'n_upper':n,'scalar':str(scalar),'shift5':list(shift)}
  if source is not None:rec['source']=source
  if old_origin is not None:rec['old_origin']=old_origin
  own(count_numeric(rec)+12000);k=len(nodes);nodes.append(rec);node_numeric+=count_numeric(rec);names[name]=k;return k
 def newop(op,args,scalar=Q(1),shift=(0,0,0,0,0),name=None):
  key=(op,tuple(args),str(scalar),tuple(shift))
  if key in dedup:return dedup[key]
  A=nodes[args[0]];qa=Q(A['q']);Ma=int(A['M']);ba=A['box5'];la=A['lead5'];na=A['n_upper']
  if not qa:Ma=0;ba=None;la=None;na=0
  if op=='copy':
   q=qm(qa,scalar);M=Ma;box=None if ba is None else[[v+shift[j]for v in b]for j,b in enumerate(ba)];lead=None if la is None else tuple(la[j]+shift[j]for j in range(5));n=na
  else:
   B=nodes[args[1]];qb=Q(B['q']);Mb=int(B['M']);bb=B['box5'];lb=B['lead5'];nb=B['n_upper']
   if not qb:Mb=0;bb=None;lb=None;nb=0
   if op=='multiply':
    q=qm(qm(qa,qb),scalar)
    if ba is None or bb is None:box=None;lead=None;n=0;M=0
    else:
     box=[[ba[j][i]+bb[j][i]+shift[j]for i in range(2)]for j in range(5)];lead=tuple(la[j]+lb[j]+shift[j]for j in range(5))
     assert Ma.bit_length()+Mb.bit_length()+min(na,nb).bit_length()<=2048
     M=Ma*Mb*min(na,nb);n=min(na*nb,prod(hi-lo+1 for lo,hi in box))
   else:
    assert op=='add';qb=qm(qb,scalar);q=qg(qa,qb)
    if not qb:Mb=0;bb=None;lb=None;nb=0
    if q:wa=qa/q;wb=qb/q;assert wa.denominator==wb.denominator==1;M=abs(wa.numerator)*Ma+abs(wb.numerator)*Mb
    else:M=0
    boxes=[b for b in[ba,bb]if b is not None];leads=[tuple(v)for v in[la,lb]if v is not None]
    box=None if not boxes else[[min(b[j][0]for b in boxes),max(b[j][1]for b in boxes)]for j in range(5)]
    lead=min(leads)if leads else None;n=0 if box is None else min(na+nb,prod(hi-lo+1 for lo,hi in box))
  if not q:M=0;lead=None;box=None;n=0
  k=add_record(name or('N'+str(len(nodes))),op,args,q,M,lead,box,n,scalar,shift);dedup[key]=k;return k
 def binding(scope,label):
  key=(scope,label)
  if key in leafcache:return leafcache[key]
  d=load(OLD+'/'+scope);v=(d['common_operands']if'common_operands'in d else d['row_operands'])[label]
  q,M,l,b,n=primitive(v['original_literal'],Q(v['QQ_transform_scalar']),v['transform_monomial_exponent5'])
  k=add_record('B_'+scope+'_'+label,'input',(),q,M,l,b,n,source={'binding':{'scope':scope,'label':label}})
  leafcache[key]=k;return k
 def cached(stem,rec):
  key=(stem,rec['file'])
  if key in leafcache:return leafcache[key]
  p=ROOT/stem/rec['file'];pins[str(p.relative_to(ROOT))]=sha(p);assert pins[str(p.relative_to(ROOT))]==rec['sha256'];body=p.read_bytes();assert body[:9]==b'Q5CELL1\n\1'
  q,M,l,b,n=primitive(decode(body[9:]));k=add_record('CACHE_'+str(len(nodes)),'input',(),q,M,l,b,n,source={'cell_file':str(p.relative_to(ROOT)),'sha256':pins[str(p.relative_to(ROOT))]})
  leafcache[key]=k;return k
 def alarm(s,f):raise TimeoutError('SOURCE55s; no automatic retry')
 signal.signal(signal.SIGALRM,alarm);signal.alarm(55)
 try:
  pre=load(PREF+'/report.json');assert pins[PREF+'/report.json']=='359acee48739dc5e99b14bd24d9a65791fd9704510b6cc8ad308b7309b9b870b'and pre['all_passed']
  plan=load(PREF+'/chosen-functional-and-thin-source-inputs.json')
  # In-place compact the accepted full original graph before loading old700.
  vp='reports/S61_COL53_VERTEX_SOURCE_GATE_V2_20261105_v64/graph.json';D=load(vp);assert pins[vp]=='26a7c0109615e546767a97091f9ff290868fdffccc12af297fd77415562be42d'
  original=D.pop('nodes');assert len(original)==1203
  # Charge cold original graph separately before compacting it.
  cold=25000+sum(count_numeric(d)for d in docs.values())+count_numeric(original)+4096;assert cold<=150000;peak=max(peak,cold)
  for k,old in enumerate(original):
   context={'stage':'old1203_compaction','node':k}
   op=old['op'];op={'immutable_SOURCE_view':'copy','unmaterialized_multiply_for_SOURCE_union':'multiply','SOURCE_fused_union_product':'add'}.get(op,op)
   q=Q(old['source_profile']['known_QQ_content']);M=int(old['source_profile']['primitive_integer_abs_upper']);lead=old['leads'][0];args=old['args'];shift=old.get('shift',[0]*5)
   if op=='input':box=old['literal_profile']['exponent_box5']if old['literal_profile']['terms_upper']else None;n=old['literal_profile']['terms_upper']
   elif op=='copy':
    b=nodes[args[0]]['box5'];box=None if b is None else[[x+shift[j]for x in v]for j,v in enumerate(b)];n=nodes[args[0]]['n_upper']
   elif op=='multiply':
    a,b=(nodes[j]['box5']for j in args);box=None if a is None or b is None else[[a[j][i]+b[j][i]+shift[j]for i in range(2)]for j in range(5)];n=0 if box is None else min(nodes[args[0]]['n_upper']*nodes[args[1]]['n_upper'],prod(hi-lo+1 for lo,hi in box))
   else:
    boxes=[nodes[j]['box5']for j in args if nodes[j]['box5']is not None];box=None if not boxes else[[min(b[j][0]for b in boxes),max(b[j][1]for b in boxes)]for j in range(5)];n=0 if box is None else min(sum(nodes[j]['n_upper']for j in args),prod(hi-lo+1 for lo,hi in box))
   # No second copy of the legacy ten-lead arrays is retained.
   rec={'name':old['name'],'op':op,'args':args,'q':str(q),'M':str(M),'lead5':lead,'box5':box,'n_upper':n,'scalar':old.get('scalarB','1')if op=='add'else old.get('scalar','1'),'shift5':shift}
   if op=='input':rec['source']={'legacy_source':old['source'],'literal_profile':old['literal_profile']}
   old.clear();old.update(rec);nodes.append(old);names[old['name']]=k
   assert max(q.denominator.bit_length(),abs(q.numerator).bit_length(),M.bit_length())<=2048
  del original;docs.pop(vp);D.clear();del D
  node_numeric=sum(count_numeric(n)for n in nodes)
  own()
  fullF11=1202;assert nodes[fullF11]['name']=='new_1013_sum'
  prog=load(OLD+'/staged-literal-col56-program-and-positive-profiles.json');opmap={v['output']:v for v in prog['literal_ops']};oldprofiles=prog['positive_support_and_height_profiles']
  def oldnode(name):
   if 'O_'+name in names:return names['O_'+name]
   op=opmap[name];kind=op['op']
   if kind=='input':k=binding(op['binding']['scope'],op['binding']['label'])
   elif kind=='copy':k=newop('copy',[oldnode(op['A'])],Q(op.get('scalar','1')),op.get('shift5',[0]*5),'O_'+name)
   elif kind=='multiply':k=newop('multiply',[oldnode(op['A']),oldnode(op['B'])],Q(op.get('scalar','1')),op.get('shift5',[0]*5),'O_'+name)
   else:k=newop('add',[oldnode(op['A']),oldnode(op['B'])],Q(op.get('scalarB','1')),name='O_'+name)
   names['O_'+name]=k;profile=oldprofiles[name]
   # The old profile is an already certified tighter coefficient majorant.
   assert Q(nodes[k]['q'])==Q(profile['known_QQ_content'])
   nodes[k]['M']=profile['primitive_integer_abs_upper'];nodes[k]['n_upper']=profile['terms_upper'];return k
  neededold=[oldnode('row'+str(i)+'_N8_sum_m')for i in[6,17,58]]
  for rr in[6,17,35,58]:
   for label in['U','V','C','F','N','M']:oldnode(label+str(rr))
  for label in['Q','P7','E']:oldnode(label)
  # All later oldnode lookups now return at the name guard. Release the whole
  # parse-only700-op/profile documents before any new coefficient-side graph.
  for f in list(docs):
   if f.startswith(OLD+'/'):docs.pop(f)
  del prog,opmap,oldprofiles
  # Cache physical input payloads in the copied48 thin pack only.
  pack=ROOT/PREF/'selected48-existing-cell-payloads.bin';pins[str(pack.relative_to(ROOT))]=sha(pack)
  desc={(x['row'],x['column']):x for x in plan['selected48_literal_payloads']};raw={}
  def rawnode(row,col):
   key=(row,col)
   if key in raw:return raw[key]
   d=desc[key]
   if d.get('known_existing_zero_by_complete_index'):q,M,l,b,n=Q(0),0,None,None,0;src={'known_zero':True}
   else:
    with pack.open('rb')as f:f.seek(d['new_offset']);body=f.read(d['bytes'])
    assert hashlib.sha256(body).hexdigest()==d['payload_sha256'];q,M,l,b,n=primitive(decode(body));src={'thin_payload_file':str(pack.relative_to(ROOT)),'offset':d['new_offset'],'bytes':d['bytes'],'payload_sha256':d['payload_sha256'],'original_cell':d['old']}
   k=add_record('RAW_'+str(row)+'_'+str(col),'input',(),q,M,l,b,n,source=src);raw[key]=k;return k
  # FullF11 remains all1203 unpruned nodes; new expressions use uniform117 terms.
  side=load(SIDE+'/exact-11-row-coefficient-operator-and-callback-owners.json')
  a15=load('reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64/actual15-col56-and-once-shared-row35.json');I={r['row']:r for r in a15['actual15_ninth_col56_numerators']}
  a16=load('reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64/actual-cleared-col82-and-commonQ3RE-representation.json');B={r['row']:r for r in a16['actual16_col82_numerator_values']}
  one=add_record('ONE','input',(),Q(1),1,(0,0,0,0,0),[[0,0]]*5,1,source={'QQ_monomial':True,'coefficient':'1','exponent5':[0]*5})
  mono_cache={}
  def monoshift(name,power):
   key=(name,power)
   if key not in mono_cache:
    j={'H':0,'s':2,'q':3,'t':4}[name];v=[0]*5;v[j]=power;mono_cache[key]=newop('copy',[one],shift=v)
   return mono_cache[key]
  bb58=cached('reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64',B[58])
  ni58=cached('reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64',I[58])
  UF10=newop('copy',[ni58],Q(353894400000),[1,0,0,0,1])
  roots={'6,53':fullF11}
  for row,columns in[(6,[79]),(17,[53,79])]:
   V={name:oldnode(name)for name in['Q','P7','E']}
   for tag in['i','35','58']:
    rr=row if tag=='i'else int(tag)
    for f in['U','V','C','F','N','M']:V[f+tag]=oldnode(f+str(rr))
   for tag in['i','58','31','83','14','35','104']:
    rr=row if tag=='i'else int(tag)
    for c in[49,51]:V['A'+tag+'_'+str(c)]=rawnode(rr,c)
   ni=cached('reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64',I[row]);bi=cached('reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64',B[row])
   V['UF10']=UF10;V['C10_Ni56']=newop('copy',[ni],Q(353894400000),[1,0,0,0,1])
   aa=newop('multiply',[bb58,oldnode('row'+str(row)+'_N8_sum_m')]);ab=newop('multiply',[bi,oldnode('row58_N8_sum_m')])
   alpha=newop('add',[aa,ab],Q(-1))
   V['C9_C10_alpha_i']=newop('copy',[alpha],Q(62208000000)*Q(353894400000),[1,0,1,1,3])
   weights=[]
   for rawrow,polynomial in zip(side['raw_row_ids11'],side['raw_coefficient_side11']):
    terms=[]
    for mon,c in polynomial:
     factors=[]
     for name,power in mon:
      if name in['H','s','q','t']:factors.append(monoshift(name,power))
      else:
       assert power>=0;factors.extend([V[name]]*power)
     term=one
     for f in factors:term=newop('multiply',[term,f])
     term=newop('copy',[term],Q(c));terms.append(term)
    weights.append((rawrow,terms))
   for col in columns:
    terms=[]
    for rawrow,wt in weights:
     rr=row if rawrow=='i'else int(rawrow);v=rawnode(rr,col)
     for w in wt:terms.append(newop('multiply',[w,v]))
    ans=terms[0]
    for t in terms[1:]:ans=newop('add',[ans,t])
    roots[str(row)+','+str(col)]=ans
  diag=newop('multiply',[fullF11,roots['17,79']])
  cross=newop('multiply',[roots['17,53'],roots['6,79']])
  minor=newop('add',[diag,cross],Q(-1),name='B17_79')
  assert nodes[minor]['lead5']==[-5,0,7,-3,19],nodes[minor]['lead5']
  target=plan['new_functional_request']['exponent5'];assert target==[-5,0,9,-3,19]
  # SOURCE propagation budget only. The actual support profile is a later ROOT GO.
  maxdegree={};active={}
  def need(k,d):
   if Q(nodes[k]['q'])==0:maxdegree[k]=max(maxdegree.get(k,-1),d);active[k]=[];return
   if maxdegree.get(k,-1)>=d:return
   maxdegree[k]=d;n=nodes[k];op=n['op'];children=[]
   if op in['copy','multiply']:children=[(c,d)for c in n['args']]
   elif op=='add':
    for c in n['args']:
     cl=nodes[c]['lead5'];pl=n['lead5']
     if cl is not None and cl[:2]==pl[:2]and cl[2]-pl[2]<=d:children.append((c,d-(cl[2]-pl[2])))
   active[k]=[c for c,dd in children]
   for c,dd in children:need(c,dd)
  need(minor,2);own(20000)
  for f in[Path(__file__),ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py',ROOT/'code/s61_sparse_generic5_source_view_profile_20261105_v64.py',ROOT/'code/s61_sparse_generic5_source_packed5_profile_v2_20261105_v64.py']:pins[str(f.relative_to(ROOT))]=sha(f)
  # Drop all bulky parse-only old700/operand documents before writing the result.
  maxbits=max(int(n['M']).bit_length()for n in nodes);assert maxbits<=2048
  graph={'nodes':nodes,'FULL_F11_original_root':fullF11,'four_callback_roots':roots,'minor_B_root':minor,'actual_C_factor':'(P7*F9)^2','target_full5':target,'max_relative_s_degree':2,'needed_node_degree':{str(k):d for k,d in maxdegree.items()},'active_tied_HL_children':{str(k):v for k,v in active.items()},'all_old_and_new_complements_and_RHS_kept':True,'no_actual_support_or_signed_values_evaluated':True}
  put('compiled-full-SLP-and-fixed-projection-plan.json',graph)
  raw=sum((ROOT/f).stat().st_size for f in pins)+4*1024*1024+262144;assert raw<32*1024*1024
  report={'id':ID,'all_passed':True,'SOURCE_COMPILATION_ONLY_no_support_or_signed_VALUE':True,'full_nodes':len(nodes),'needed_projected_nodes':len(maxdegree),'M_positive_bits_upper':maxbits,'SOURCE_compile_peak_numeric':peak,'fresh_RAW_upper':raw,'future_projection_SOURCE_worker_gate_pending':True,'input_sha256':pins,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};put('report.json',report)
  print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'));signal.alarm(0)
 except Exception as exc:
  if not(out/'partial-compiled-source.json').exists():put('partial-compiled-source.json',{'nodes_completed':nodes,'whole_prefix_only':True,'completed_norm_profiles_not_coefficient_values':True,'input_sha256':pins})
  if not(out/'failure.json').exists():put('failure.json',{'all_passed':False,'SOURCE_compile_only':True,'error':repr(exc),'context':context,'nodes_completed':len(nodes),'no_support_or_signed_values':True,'no_automatic_retry':True})
  raise
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');main(p.parse_args())

#!/usr/bin/env python3
"""Independent sparse packed-exponent polynomial replay; two heads only."""
import hashlib,json,pathlib,re,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1];PRI='AMBIENT_B18_C_CONNECTING_HEADS2930_CACHED_VALUES_20261105_v63';SRC='AMBIENT_B18_C_CONNECTING_HEADS2930_CACHE_SOURCE_V2_20261105_v63';SP='AMBIENT_B18_C_CONNECTING_HEADS2930_CACHE_SOURCE_CODE_20261105_v63_REVIEW';ID=PRI+'_REVIEW'
PRE=0;STORED=0;LIMIT=103
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def hb(q):return max(abs(q.numerator).bit_length(),q.denominator.bit_length())
def keep(q):
 global STORED
 STORED=max(STORED,hb(q));assert hb(q)<=384;return q
def parse(s):
 global PRE
 assert isinstance(s,(str,int));a=str(s).split('/');n=int(a[0]);d=int(a[1])if len(a)>1 else 1;assert d
 h=max(abs(n).bit_length(),abs(d).bit_length());PRE=max(PRE,h);assert h<=LIMIT and h<=1153
 return keep(Q(n,d))
def calc(a,b,m=False):
 global PRE
 if not a or not b:return keep(a*b if m else a+b)
 na,da=abs(a.numerator).bit_length(),a.denominator.bit_length();nb,db=abs(b.numerator).bit_length(),b.denominator.bit_length()
 h=max(na+nb,da+db)if m else max(na+db,nb+da)+1;PRE=max(PRE,h);assert h<=LIMIT and h<=1153
 return keep(a*b if m else a+b)
def main():
 global LIMIT
 t=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};RAW=0
 def pin(p,h=None,decode=True):
  nonlocal RAW
  p=pathlib.Path(p);k=str(p.relative_to(ROOT));s=sha(p)
  if h:assert h==s,(k,h,s)
  if k not in pins:RAW+=p.stat().st_size
  pins[k]=s;assert RAW<2_000_000
  return json.loads(p.read_text())if decode else None
 pr=pin(ROOT/'reports'/PRI/'report.json','29133a93c9c3f24d1bf228cacd4f49f8fd5233c255a3b85a159598da21cc297d')
 sr=pin(ROOT/'reports'/SRC/'report.json','fabde8ddc9d5f8c5b321c73c634c3094ba4bf9e520b7ddb47a75a10a20467b37')
 peer=pin(ROOT/'reports'/SP/'report.json','985e63285949c4503136288971c13494f898ae15e6e97325b0d17d5b87e6ec5c');assert peer['all_passed']
 cp=ROOT/'reports'/SRC/'compiled-ten-roots-mixed63-plane10-and-actual-provided-lower-inputs.json';C=pin(cp,sr['evidence_sha256'][cp.name])
 pp=ROOT/'reports'/SRC/'positive-before-actual-twohead-cache-and-whole-budget.json';profile=pin(pp,sr['evidence_sha256'][pp.name])
 assert profile['root_stored_bound']==51 and profile['root_before_operation_bound']==103 and profile['post_partial_stored_bound']==118 and profile['post_before_operation_bound']==237
 mons={63:list(map(tuple,C['basis63_monomials'])),10:list(map(tuple,C['plane10_monomials']))}
 def code(m):return sum(5**i for i in m)
 codes={p:list(map(code,mm))for p,mm in mons.items()};allowed={p:set(v)for p,v in codes.items()};degree={code(m):len(m)for m in mons[63]}
 assert len(set(codes[63]))==63 and len(set(codes[10]))==10
 inputs=C['provided105_lower_extension63'];terms=sum(q!='0'for i,a in inputs for q in a)
 metadata=count([pr,sr,peer,C,profile,codes,degree,pins])+16000
 forecast=metadata+2*terms+2*105+2*profile['lastuse_peak_QQ_slots']+2*206+6*2*63+6000
 assert forecast<150000 and RAW+450000<2_000_000
 save(out/'positive-before-independent-signed-packed-polynomial-replay.json',{'fresh_prior_source_peer':pins[str((ROOT/'reports'/SP/'report.json').relative_to(ROOT))],'packed_exponent_radix5_with_totaldegree_cut_before_add':True,'all_cache_term_keys_and_values_charged':2,'whole_numeric_forecast':forecast,'root_stored_preop':[51,103],'post_stored_preop':[118,237],'no_new_implicit_physical_gauge_jet_or_heads':True})
 cache={};outputs={};work={};powers=[];given={};peak=metadata;cursor={};root_hash=hashlib.sha256();units=[]
 def own(extra=0):
  nonlocal peak
  n=metadata+2*sum(map(len,given.values()))+2*len(given)+2*sum(map(len,cache.values()))+len(cache)+2*sum(map(len,outputs.values()))+2*len(work)+sum(2*len(p)for p in powers)+extra
  peak=max(peak,n);assert n<=forecast and n<150000
 def project(a,p):return a if p==63 else{k:q for k,q in a.items()if k in allowed[10]}
 def polyadd(a,b):
  z=dict(a)
  for k in sorted(b,reverse=True):
   q=calc(z.get(k,Q(0)),b[k])
   if q:z[k]=q
   else:z.pop(k,None)
  return z
 def scale(a,q):return{k:calc(v,q,True)for k,v in a.items()if v and q}
 def multiply(a,b,p):
  z={};maxdegree=4 if p==63 else 3
  for e in sorted(a,reverse=True):
   for f in sorted(b,reverse=True):
    if degree[e]+degree[f]>maxdegree:continue
    k=e+f # degree<=4 guarantees no radix5 carry in any coordinate
    if k not in allowed[p]:continue
    q=calc(z.get(k,Q(0)),calc(a[e],b[f],True))
    if q:z[k]=q
    else:z.pop(k,None)
  return z
 try:
  for node,arr in inputs:given[node]={k:parse(q)for k,q in zip(codes[63],arr)if q!='0'}
  entries=C['original_records_precision_and_retirement'];assert len(entries)==8081
  roots=C['delta_safe_root_pair']+[r['original_node']for r in C['odd_trusted_root_requests']];groots=set(C['delta_safe_root_pair'])
  for rec,p,gone in entries:
   i,op=rec[:2];cursor={'phase':'original_scalar_sparse_replay','node':i,'precision':p,'op':op};powers=[]
   if op=='input':work=project(given[i],p)
   elif op=='c':
    q=parse(str(Q(rec[2],rec[3])));work={0:q}if q else{}
   elif op=='a':work=polyadd(project(cache[rec[2]],p),project(cache[rec[3]],p))
   elif op=='n':work=scale(project(cache[rec[2]],p),Q(-1))
   elif op=='m':work=multiply(project(cache[rec[2]],p),project(cache[rec[3]],p),p)
   elif op=='i':
    a=project(cache[rec[2]],p);assert a.get(0,Q(0))==1
    nil=dict(a);nil.pop(0);power={0:Q(1)}
    for k in range(1,5):power=multiply(power,nil,p);powers.append(scale(power,Q((-1)**k)))
    work={0:Q(1)}
    for power in reversed(powers):work=polyadd(work,power)
    assert multiply(a,work,p)=={0:Q(1)}==multiply(work,a,p);units.append([i,rec[2],'1',p])
   else:raise AssertionError(op)
   cache[i]=work;work={};own(2000)
   vector=[str(cache[i].get(k,Q(0)))for k in codes[p]]
   root_hash.update((json.dumps([i,p,vector],separators=(',',':'))+'\n').encode())
   if i in roots:outputs[i]=project(cache[i],63 if i in groots else 10).copy()
   for j in gone:del cache[j]
   assert time.monotonic()-t<48
  assert not cache
  delta=polyadd(outputs[20144],scale(outputs[40922],Q(-1)))
  assert STORED<=51
  # Read the primary root certificate ONLY after independent completion.
  p1name='phase1-safe-delta63-and-eight-trusted-bB-cache-readbacks.json';P1=pin(ROOT/'reports'/PRI/p1name,pr['evidence_sha256'][p1name])
  assert P1['all_node_chosen_lower_extension_value_hash']==root_hash.hexdigest()
  assert P1['all_unit_both_products']==units
  assert [str(delta.get(k,Q(0)))for k in codes[63]]==P1['safe_delta63_coefficients']
  expectedodd={node:(row,v)for row,node,v in P1['actual_bB_plane10_roots']}
  for req in C['odd_trusted_root_requests']:
   node=req['original_node'];assert set(outputs[node])<=allowed[10]
   assert expectedodd[node][0]==req['physical_row']and[str(outputs[node].get(k,Q(0)))for k in codes[10]]==expectedodd[node][1]
   for m,q in req['root_lower_checks']:assert outputs[node].get(code(m),Q(0))==parse(q)
  assert delta.get(code([7,7,7,7]),Q(0))==parse(C['known_delta63_plane4_r4_r3ell'][0])
  assert delta.get(code([7,7,7,8]),Q(0))==parse(C['known_delta63_plane4_r4_r3ell'][1])
  assert outputs[38648].get(code([7,8,8]),Q(0))==parse(C['known_bB63_r_ell2'])
  assert delta.get(0,Q(0))==0 and all(delta.get(code([i]),Q(0))==0 for i in range(9))
  # All14 literal nonzz quadratic coefficients, not just the one needed for
  # the inactive beta3 clarification, are checked from independent output.
  c2path=ROOT/'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json';pin(c2path,sr['input_sha256'][str(c2path.relative_to(ROOT))],False)
  phys=next(r for r in field(c2path,'physical12_coefficients')if r['physical_row']==63)
  rawpath=ROOT/'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json';pin(rawpath,decode=False)
  rb=field(rawpath,'actual33_basis');checks=[]
  for i in range(7):
   for v in[7,8]:
    slot=rb.index([i+1,int(v==7),int(v==8)]);q=delta.get(code([i,v]),Q(0));assert q==parse(dict(phys['delta']).get(slot,'0'))==0;checks.append([i,v,str(q)])
  assert delta.get(code([1,8]),Q(0))==0
  root_summary={'all8081_original_node_values_hash':root_hash.hexdigest(),'full_safe_delta63_matches':True,'all_eight_bB_plane10_arrays_match':True,'all14_nonzz_quadratics':checks,'original_root_unit_both_products':units,'root_max_stored':STORED,'root_max_pre':PRE}
  save(out/'completed-independent-root-cache-checkpoint.json',root_summary)
  headpaths=C['head29_head30_paths'];requests=C['odd_trusted_root_requests'];rowroot={r['physical_row']:r['original_node']for r in requests}
  del entries,inputs,given,cache,powers,work,P1
  C.pop('original_records_precision_and_retirement');C.pop('provided105_lower_extension63')
  LIMIT=237;leaf={};heads=[]
  for path in headpaths:
   pin(ROOT/path,sr['input_sha256'][path],False)
   h={k:field(ROOT/path,k)for k in['known_cached_inputs','four_genuine_and_canonical_high_polynomials_UNPAIRED']};heads.append(h)
   for label,q,locator in h['known_cached_inputs']:
    k=tuple(label);v=parse(q)
    if k in leaf:assert leaf[k]==v
    leaf[k]=v
  for key,m in[(('NEW_delta63_plane4','r4'),[7,7,7,7]),(('NEW_delta63_plane4','r3ell'),[7,7,7,8]),(('NEW_delta63_plane4','r2ell2'),[7,7,8,8])]:leaf[key]=delta.get(code(m),Q(0))
  leaf[('NEW_CACHED_bB63_ell3',)]=outputs[38648].get(code([8,8,8]),Q(0))
  for m,k in zip(mons[63],codes[63]):leaf[('SAFE_delta63',)+m]=delta.get(k,Q(0))
  for row,node in rowroot.items():leaf[('TRUSTED_plane3','bB',row,1,2)]=outputs[node].get(code([7,8,8]),Q(0))
  values=[];termwise=[]
  for hi,h in enumerate(heads):
   vv=[];tt=[]
   for j,poly in enumerate(h['four_genuine_and_canonical_high_polynomials_UNPAIRED']):
    cursor={'phase':'reversed_postoperator','head':hi,'coordinate':j};terms=[None]*len(poly);result=Q(0)
    for i in reversed(range(len(poly))):
     labels,weight=poly[i];term=parse(weight)
     for label in reversed(labels):term=calc(term,leaf[tuple(label)],True)
     result=calc(result,term);terms[i]=str(term)
    vv.append(result);tt.append(terms)
   assert vv[2]==calc(parse('5/3'),vv[0],True)
   assert vv[3]==calc(vv[1],calc(parse('-121/48'),vv[0],True))
   values.append(vv);termwise.append(tt)
  assert STORED<=118
  p2name='phase2-actual-full-head29-head30-and-four-Kbar-parameters.json';P2=pin(ROOT/'reports'/PRI/p2name,pr['evidence_sha256'][p2name])
  assert list(map(str,values[0]))==P2['head29_A_B_rg_k2']and list(map(str,values[1]))==P2['head30_A_B_rg_k2']
  assert termwise==P2['all_eight_full_post_termwise_products']
  assert len(leaf)<400 and sum(len(poly)for h in heads for poly in h['four_genuine_and_canonical_high_polynomials_UNPAIRED'])==2673
  assert RAW+60000<2_000_000 and all(sha(ROOT/p)==h for p,h in pins.items())
  proof={'independent_signed_sparse_packed_exponent_polynomial_algebra':True,'all8081_node_hash_matches_after_independent_replay':root_hash.hexdigest(),
   'inverse_method':'nilpotent signed geometric series, reversed power summation; reversed sparse product/add order; source original n/d parsed exactly',
   'all_safe_delta63_and_eight_projected_bB_plane10_arrays_match':True,'shared38648_odd_capture_projected63_to10':True,
   'all2673_postproducts_and_eight_operator_outputs_match':True,'head29_A_B_rg_k2':list(map(str,values[0])),'head30_A_B_rg_k2':list(map(str,values[1])),
   'explicit_actual_delta_z1ell_zero_and_all14_nonzz_quadratics':checks,
   'generic_q7_completeness_source':'Own295b/a5ff independently derived the coherent operators. Actual omitted inactive-beta3 term is zero by this fresh safe-delta output; seven bA plane3 cancellations are literal source aliases. All postleaves now come from this independent cache or pinned prior inputs.',
   'lower_extension_scope':'Only DELTA difference true on full63 by all105 gradient0. Individual even lower-skeleton values and full63 bB outputs are not asserted actual; bB truth onlyplane10. No implicitPhi4=0 assertion.',
   'only_four_new_Kbar_source_basis_parameters':{k:str(v)for k,v in zip(['p2','b2','p3','b3'],values[0][2:]+values[1][2:])},
   'no_other_head_high_Kbar_actions_extension_originalH0_or_nil_inference':True,
   'fresh_forecast_numeric':forecast,'actual_tracked_numeric_peak':peak,'actual_max_stored_bits':STORED,'actual_max_preoperation_bits':PRE,'RAW_unique':RAW}
  ep=out/'independent-sparse-twohead-all-root-coefficients-and-all-termwise-values.json';save(ep,proof)
  report={'id':ID,'all_passed':True,'frozen':True,'independent_actual_cached_SOURCE_VALUE':True,'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},'full_HEAD29_HEAD30_values_PASS':True,'new_source_solve_or_jet':False,'seconds':time.monotonic()-t}
  save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'forecast':forecast,'peak':peak,'stored':STORED,'pre':PRE,'seconds':report['seconds']}))
 except BaseException as exc:
  with(out/'failure-live-prefix.jsonl').open('w')as f:
   for node,v in cache.items():
    for k,q in v.items():f.write(json.dumps([node,k,str(q)],separators=(',',':'))+'\n')
  save(out/'failure.json',{'error':repr(exc),'cursor':cursor,'no_automatic_retry_or_capraise':True,'root_checkpoint_exists':(out/'completed-independent-root-cache-checkpoint.json').exists()});raise
if __name__=='__main__':main()

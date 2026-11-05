#!/usr/bin/env python3
"""Independent mixed63/10 source/custody audit; no signed root response."""
import ast,collections,hashlib,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1];SRC='AMBIENT_B18_C_CONNECTING_HEADS2930_CACHE_SOURCE_V2_20261105_v63';CUS='AMBIENT_B18_C_CONNECTING_HEADS2930_VALUE_CUSTODY_20261105_v63';ID='AMBIENT_B18_C_CONNECTING_HEADS2930_CACHE_SOURCE_CODE_20261105_v63_REVIEW'
SOURCE_HASH='fabde8ddc9d5f8c5b321c73c634c3094ba4bf9e520b7ddb47a75a10a20467b37';WORKER_HASH='faf2c7f4a25c844ca3ab6e45be05f2b6d2c30be43a69a2e591815318d5c1fbd1'
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
class Majorants:
 def __init__(self):self.H=0
 def track(self,v):
  h=max(abs(v[0].numerator).bit_length(),v[0].denominator.bit_length())if isinstance(v[0],Q)else(v[0]*v[1]).bit_length()+1
  self.H=max(self.H,h);assert self.H<=384 and 2*self.H+1<=1153;return v
 def bound(self,q):return self.track(((abs(q.numerator)+q.denominator-1)//q.denominator,q.denominator))
 def plus(self,a,b):
  if isinstance(a[0],Q)and isinstance(b[0],Q):return self.track((a[0]+b[0],1))
  if isinstance(a[0],Q):a=self.bound(a[0])
  if isinstance(b[0],Q):b=self.bound(b[0])
  return self.track((a[0]+b[0],math.lcm(a[1],b[1])))
 def times(self,a,b):
  if isinstance(a[0],Q)and isinstance(b[0],Q):return self.track((a[0]*b[0],1))
  if isinstance(a[0],Q):a=self.bound(a[0])
  if isinstance(b[0],Q):b=self.bound(b[0])
  return self.track((a[0]*b[0],a[1]*b[1]))if a[0]and b[0]else(0,1)
def main():
 t=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};RAW=0
 def pin(p,h=None,parse=True):
  nonlocal RAW
  p=pathlib.Path(p);k=str(p.relative_to(ROOT));s=sha(p)
  if h:assert h==s,(k,h,s)
  if k not in pins:RAW+=p.stat().st_size
  pins[k]=s
  return json.loads(p.read_text())if parse else None
 sr=pin(ROOT/'reports'/SRC/'report.json',SOURCE_HASH)
 for p,h in sr['input_sha256'].items():
  pin(ROOT/p,h,False)
  if p.endswith('_REVIEW/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
 producer=ROOT/'code/ambient_B18_C_connecting_heads2930_cache_source_v2_20261105_v63.py';pin(producer,sr['code_sha256'],False)
 profpath=ROOT/'reports'/SRC/'positive-before-actual-twohead-cache-and-whole-budget.json';profile=pin(profpath,sr['evidence_sha256'][profpath.name])
 cust=pin(ROOT/'reports'/CUS/'report.json');cn=next(iter(cust['evidence_sha256']));CUST=pin(ROOT/'reports'/CUS/cn,cust['evidence_sha256'][cn])
 assert cust['all_passed']and CUST['all_passed']and CUST['code_sha256']==WORKER_HASH
 for p,h in CUST['input_sha256'].items():pin(ROOT/p,h,False)
 worker=ROOT/CUST['code_path'];pin(worker,WORKER_HASH,False)
 def src(s):return ROOT/next(p for p in sr['input_sha256']if p.endswith(s))
 original=field(src('compiled-original164-DAG-and-provided-pure-ell-inputs.json'),'operator_records_originalID_op_args')
 om={r[0]:r for r in original};phys={r['physical_row']:r for r in field(src('actual-physical12-and-C2-right-exact-quotient.json'),'physical12_coefficients')}
 oddrows=[18,19,28,29,60,63,67,72];groots=[20144,40922];oroots=[phys[i]['roots']['bB_node']for i in oddrows];roots=groots+oroots;assert len(set(roots))==10
 def closure(roots):
  wanted=set();todo=list(roots)
  while todo:
   i=todo.pop()
   if i in wanted:continue
   wanted.add(i);r=om[i]
   if r[1]not in['c','input']:todo.extend(r[2:])
  return wanted
 gneed=closure(groots);need=closure(roots);expected=[om[i]for i in sorted(need)]
 phase0=count([original,sr,profile,cust,CUST,phys,pins])+len(om)+len(gneed)+len(need)+len(expected)+15000
 assert phase0<150000 and RAW+2000000<32<<20
 save(out/'positive-before-original-source-compaction.json',{'SOURCE_only':True,'no_signed_root_or_head_value':True,'all_original_record_numeric_labels_and_maps_charged':True,'phase0_upper':phase0})
 del original,om,need
 cp=ROOT/'reports'/SRC/'compiled-ten-roots-mixed63-plane10-and-actual-provided-lower-inputs.json';C=pin(cp,sr['evidence_sha256'][cp.name])
 entries=C['original_records_precision_and_retirement'];assert len(entries)==len(expected)==8081
 precision={r[0]:p for r,p,gone in entries}
 assert sum(p==63 for p in precision.values())==4387 and sum(p==10 for p in precision.values())==3694
 assert groots==C['delta_safe_root_pair']and [r['original_node']for r in C['odd_trusted_root_requests']]==oroots
 uses=collections.Counter(roots);active=set();maxlive=0;maxQQ=0
 for r in expected:
  if r[1]not in['c','input']:uses.update(r[2:])
 for (r,p,gone),r0 in zip(entries,expected):
  assert r==r0 and p==(63 if r[0]in gneed else 10)
  i,op=r[:2];active.add(i);maxlive=max(maxlive,len(active));maxQQ=max(maxQQ,sum(precision[j]for j in active))
  if i in roots:uses[i]-=roots.count(i)
  retire=[]
  if op not in['c','input']:
   for j in r[2:]:
    assert p==10 or precision[j]==63
    uses[j]-=1
    if not uses[j]:active.remove(j);retire.append(j)
  if not uses[i]:active.remove(i);retire.append(i)
  assert retire==gone
 assert not active and maxlive==241==profile['lastuse_peak_nodes']and maxQQ==8662==profile['lastuse_peak_QQ_slots']
 assert sum(len(gone)for r,p,gone in entries)==8081
 assert precision[38648]==63 and [38648,10]in C['captured_output_precisions']
 assert C['captured_output_precisions']==[[i,63 if i in groots else 10]for i in roots]
 del expected,gneed,uses,active
 mon63=list(map(tuple,C['basis63_monomials']));mon10=list(map(tuple,C['plane10_monomials']));mons={63:mon63,10:mon10};mi={p:{m:i for i,m in enumerate(mm)}for p,mm in mons.items()};i10=[mi[63][m]for m in mon10]
 assert mon10==[m for m in mon63 if all(i>=7 for i in m)and len(m)<=3]
 B63=src('finite63-socle-safe-output-theorem-and-actual-lower-input-overlaps.json')
 assert C['provided105_lower_extension63']==field(B63,'provided105_lower_extension63')
 assert C['basis63_monomials']==field(B63,'basis63_monomials')
 pp={p:[(i,j,ii[tuple(sorted(a+b))])for i,a in enumerate(mm)for j,b in enumerate(mm)if tuple(sorted(a+b))in ii]for p,mm in mons.items()for ii in[mi[p]]}
 arithmetic_prior=count([sr,profile,cust,CUST,C,pp,mon63,mon10,pins,phys])+2*105*63+2*maxQQ+18000
 assert arithmetic_prior<150000
 save(out/'positive-before-independent-majorants.json',{'whole_SOURCE_arithmetic_prior':arithmetic_prior,'different_inverse_bound_algorithm':'nilpotent geometric series, no signed coefficients','odd_capture_project63_to10_required':True,'no_requested_values':True})
 E=Majorants()
 def zero(p):return[(Q(0),1)]+[(0,1)]*(p-1)
 def one(p):return[(Q(1),1)]+[(0,1)]*(p-1)
 def restrict(v,p):
  assert len(v)==p or(p==10 and len(v)==63)
  return v if len(v)==p else[v[j]for j in i10]
 def conv(a,b,p):
  v=zero(p)
  for i,j,k in pp[p]:v[k]=E.plus(v[k],E.times(a[i],b[j]))
  return v
 given={}
 for node,arr in C['provided105_lower_extension63']:given[node]=[(Q(q),1)if i==0 else E.bound(Q(q))for i,q in enumerate(arr)]
 assert len(given)==105
 cache={};outputs={};majorhash=hashlib.sha256();units=[];constants=[]
 for r,p,gone in entries:
  i,op=r[:2]
  if op=='input':v=restrict(given[i],p)
  elif op=='c':v=zero(p);v[0]=(Q(r[2],r[3]),1);constants.append(r)
  elif op=='a':v=[E.plus(a,b)for a,b in zip(restrict(cache[r[2]],p),restrict(cache[r[3]],p))]
  elif op=='n':v=restrict(cache[r[2]],p).copy();v[0]=(-v[0][0],1)
  elif op=='m':v=conv(restrict(cache[r[2]],p),restrict(cache[r[3]],p),p)
  elif op=='i':
   a=restrict(cache[r[2]],p);assert a[0][0]==1
   nil=a.copy();nil[0]=(Q(0),1);v=one(p);power=one(p)
   # Positive geometric series counts all compositions by degree. No actual
   # inverse response or coefficient sign is used.
   for degree in range(1,5):power=conv(power,nil,p);v=[E.plus(x,y)for x,y in zip(v,power)]
   conv(a,v,p);conv(v,a,p);units.append([i,r[2],'1'])
  else:raise AssertionError(op)
  cache[i]=v;majorhash.update((json.dumps([i,p,[[str(a),str(d)]for a,d in v]],separators=(',',':'))+'\n').encode())
  if i in roots:outputs[i]=restrict(v,63 if i in groots else 10).copy()
  for j in gone:del cache[j]
 assert not cache and majorhash.hexdigest()==profile['all_node_positive_majorant_sha256']
 assert units==C['unit_inverse_centers']and constants==profile['literal_n_d_constants']
 delta=[E.plus(a,b)for a,b in zip(outputs[groots[0]],outputs[groots[1]])]
 assert E.H<=profile['root_stored_bound']==51 and 2*E.H+1<=profile['root_before_operation_bound']==103
 # Non-DAG parse/check scalars must obey the SAME root-phase contract.
 root_literals=[Q(q)for req in C['odd_trusted_root_requests']for m,q in req['root_lower_checks']]+list(map(Q,C['known_delta63_plane4_r4_r3ell']))+[Q(C['known_bB63_r_ell2'])]
 assert max(max(abs(q.numerator).bit_length(),q.denominator.bit_length())for q in root_literals)<=51
 assert max(max(abs(r[2]).bit_length(),abs(r[3]).bit_length())for r in constants)<=51
 del cache,given,entries
 heads=[];leaf={}
 for path in C['head29_head30_paths']:
  assert pins[path]==sr['input_sha256'][path]
  h={k:field(ROOT/path,k)for k in['known_cached_inputs','four_genuine_and_canonical_high_polynomials_UNPAIRED']};heads.append(h)
  for label,q,loc in h['known_cached_inputs']:
   b=E.bound(Q(q));key=tuple(label)
   if key in leaf:assert leaf[key]==b
   leaf[key]=b
 for key,m in[(('NEW_delta63_plane4','r4'),(7,7,7,7)),(('NEW_delta63_plane4','r3ell'),(7,7,7,8)),(('NEW_delta63_plane4','r2ell2'),(7,7,8,8))]:leaf[key]=delta[mi[63][m]]
 root63=phys[63]['roots']['bB_node'];assert root63==38648 and len(outputs[root63])==10
 leaf[('NEW_CACHED_bB63_ell3',)]=outputs[root63][mi[10][(8,8,8)]]
 for m,j in mi[63].items():leaf[('SAFE_delta63',)+m]=delta[j]
 for row,node in zip(oddrows,oroots):leaf[('TRUSTED_plane3','bB',row,1,2)]=outputs[node][mi[10][(7,8,8)]]
 post=[];diagnostics=[]
 for h in heads:
  result=[]
  for poly in h['four_genuine_and_canonical_high_polynomials_UNPAIRED']:
   value=(0,1)
   for labels,q in poly:
    term=(Q(q),1)
    for label in labels:term=E.times(term,leaf[tuple(label)])
    value=E.plus(value,term)
   result.append(value)
  post.append([[str(u),str(d)]for u,d in result])
  diagnostics.append([E.times((Q(5,3),1),result[0]),E.times((Q(-121,48),1),result[0]),E.plus(result[1],E.times((Q(-121,48),1),result[0]))])
 assert post==[v for tag,v in profile['post_output_U_D_bounds']]
 assert E.H<=profile['post_partial_stored_bound']==118 and 2*E.H+1<=profile['post_before_operation_bound']==237
 # Fresh literal worker audit and unrelated guard fixtures, no main/import.
 text=worker.read_text();tree=ast.parse(text);funcs={x.name:x for x in tree.body if isinstance(x,ast.FunctionDef)}
 assert {x.module for x in tree.body if isinstance(x,ast.ImportFrom)}=={'fractions','pathlib'}
 assert {a.name for x in tree.body if isinstance(x,ast.Import)for a in x.names}<={'argparse','hashlib','json','re','signal','time'}
 code=ast.Module(body=[funcs[k]for k in['hb','pre','keep','parse','add','times','divide']],type_ignores=[])
 env={'Q':Q,'re':__import__('re'),'MAXPRE':0,'MAXSTORED':0,'LIMIT':2048};exec(compile(code,str(worker),'exec'),env)
 assert env['parse']('-4/5')==Q(-4,5)and env['add'](Q(1,3),Q(2,5))==Q(11,15)
 assert env['times'](Q(-3,7),Q(14,5))==Q(-6,5)and env['divide'](Q(4,9),Q(-2,3))==Q(-2,3)
 for name,a,b in[('add',Q(1<<20),Q(1)),('times',Q(1<<20),Q(1<<20))]:
  env['LIMIT']=10
  try:env[name](a,b)
  except AssertionError:pass
  else:raise AssertionError('guard did not stop')
 env['LIMIT']=2048
 try:env['parse']('3/0')
 except AssertionError:pass
 else:raise AssertionError('zero denominator')
 required=["if i in roots:outputs[i]=restrict(v,63 if i in groots else 10).copy()",'for j in gone:del cache[j]',
  "for i,s in enumerate(arr):arr[i]=parse(s)","del compiled,records,given,cache,products,unitchecks,checks,first,requests",
  "p1hash=write(p1,first)","LIMIT=profile['post_before_operation_bound']",'assert conv(a,v,p)==one==conv(v,a,p)',
  "assert MAXSTORED<=profile['root_stored_bound']","assert MAXSTORED<=profile['post_partial_stored_bound']",
  "'preserve_phase1_no_automatic_retry':True","assert n<=100000 and n<=150000",'assert h<=LIMIT and h<=1153 and h<=2048']
 for item in required:assert item in text,item
 assert 'Counter'not in text and not(ROOT/'reports/AMBIENT_B18_C_CONNECTING_HEADS2930_CACHED_VALUES_20261105_v63').exists()
 # Rebuild phase1 components from the exact runtime layout, including all
 # numeric labels in records, retirement lists, inputIDs and product tables.
 entry=C['original_records_precision_and_retirement']
 runtime1={'records_precision_and_precompiled_retirement':count(entry),'runtime_dependency_Counter':0,'given105_vectors_and_IDs':105*63+105,
  'mixed_precision_lastuse_QQ_and_keys':maxQQ+2*maxlive,'root_outputs_and_product_tables':sum(map(len,outputs.values()))+count(pp),
  'all_lowercheck_values_and_metadata':count(C['odd_trusted_root_requests'])+1,'controllers_copies_guard_scratch':10000}
 assert runtime1==profile['future_runtime1_numeric_components']and sum(runtime1.values())==67465
 assert max(sum(profile['future_runtime1_numeric_components'].values()),sum(profile['future_runtime2_numeric_components'].values()))==67465<=75000<100000
 assert CUST['future_static_VALUE_numeric_upper']==67465 and CUST['root_stored_preop']==[51,103]and CUST['post_stored_preop']==[118,237]
 assert CUST['whole_fresh_RAW_upper']==CUST['actual_input_RAW']+CUST['future_output_RAW_bound']+CUST['future_peer_controller_RAW_reserve']==1064308<2000000
 assert sum(len(h['four_genuine_and_canonical_high_polynomials_UNPAIRED'][j])for h in heads for j in range(4))==CUST['known_term_products']==2673
 assert (2*36+2)*2673+50000<CUST['future_output_RAW_bound']
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'all8081_literal_original_DAG_records_and_root_names_verified':True,'precisions63_10':[4387,3694],
  'all_original_retirement_and_source_closure_lastuses_rebuilt':[maxlive,maxQQ],
  'shared_bB38648_computed63_but_captured_projected10_and_all_odd_outputs10':True,
  'given105_lower_extension_equals_accepted63_SOC_source':True,
  'safe_delta_only_other_odd_outputs_trusted_plane10':True,
  'all8081_majorant_fingerprint_independently_rebuilt_with_geometric_inverse':majorhash.hexdigest(),
  'both_inverse_diagnostics_and_all_root_lower_parse_checks_in_root51_103':True,
  'all_four_plus_four_postoperator_majorants_and_extra_canonical_diagnostics_in118_237':True,
  'source_gates_and_inactivebeta3_zero_companion_pinned':True,
  'standard_library_worker_guard_fixtures_and_literal_phase1_lifetime_PASS':True,
  'runtime1_owned_components':runtime1,'future_whole_static_VALUE_upper':67465,'hard_owned':100000,
  'source_review_phase0_upper':phase0,'source_review_arithmetic_prior':arithmetic_prior,
  'fresh_VALUE_RAW_upper':1064308,'review_unique_RAW':RAW,
  'code_caps':'root51/pre103,post118/pre237,retained384/1153,global2048/150k/32MiB,time60s,memory2GB; no capraise',
  'no_signed_root_or_head_response_composed':True,'worker_currently_UNRUN':True,
  'ROOT_actual_ONE_GO_still_required':'This report establishes readiness only; no actual value permission inferred. Execute exact frozen worker and this source peer argument under the ROOT condition.',
  'no_original_H0_Kbar_actions_extension_nil_inference':True}
 ep=out/'independent-mixed-precision-source-all-majorants-retirement-and-unrun-code-custody.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'SOURCE_CODE_ONLY':True,'input_sha256':pins,'source_report_sha256':SOURCE_HASH,'worker_sha256':WORKER_HASH,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},'worker_values_UNRUN':True,'seconds':time.monotonic()-t}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'source_phase0':phase0,'source_arithmetic':arithmetic_prior,'value_static':67465,'seconds':report['seconds']}))
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Independent prior/source, absolute-majorant and unrun ETA4 worker audit."""
import ast,hashlib,json,math,pathlib,re,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_FULL_ETA4_CACHE_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
SOURCE='5052f2af263e5168bdb09c73ec4ec9affdac04548713b5f49f76b3acf05fe0e7'
WORKER='code/ambient_B18_C_full_eta4_cache_values_20261105_v63.py'
WORKER_SHA='c6e00c676b82e268a066269937b6e9f3c66911dbc58b5900ff8539c208c7d0ff'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def nc(x):
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(list,tuple)):return sum(map(nc,x))
 if isinstance(x,dict):return sum(map(nc,x.values()))
 return 0
def enc(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def ceiling(q):return (abs(q.numerator)+q.denominator-1)//q.denominator
height=0
def mark(t):
 global height
 h=max(abs(t[0].numerator).bit_length(),t[0].denominator.bit_length()) if isinstance(t[0],Q) else (t[0].bit_length()+t[1].bit_length() if t[0] else 0)
 height=max(height,h);assert height<=2048
 return t
def sum_bound(a,b):
 if isinstance(a[0],Q):
  assert isinstance(b[0],Q);return mark((a[0]+b[0],1))
 assert not isinstance(b[0],Q)
 return mark((a[0]+b[0],math.lcm(a[1],b[1])))
def product_bound(a,b):
 if isinstance(a[0],Q) and isinstance(b[0],Q):return mark((a[0]*b[0],1))
 if isinstance(a[0],Q):
  return mark((ceiling(a[0])*b[0],a[0].denominator*b[1])) if a[0] and b[0] else mark((0,1))
 if isinstance(b[0],Q):return product_bound(b,a)
 return mark((a[0]*b[0],a[1]*b[1])) if a[0] and b[0] else mark((0,1))
def zseries():return [(Q(0),1)]+[(0,1)]*3
def product_series(a,b):
 v=zseries();v[0]=product_bound(a[0],b[0])
 # Reverse-order summation gives the same final U,D, and checks every
 # absolute summand and partial sum without evaluating any positive root.
 for n in range(1,4):
  for k in reversed(range(n+1)):v[n]=sum_bound(v[n],product_bound(a[k],b[n-k]))
 return v
def inverse_bounds(x):
 assert x[0][0]==1
 # Closed positive envelopes for1/(1+x1 ell+x2 ell²+x3 ell³).
 x1,x2,x3=x[1:];square=product_bound(x1,x1)
 v=[(Q(1),1),x1,sum_bound(square,x2),
    sum_bound(sum_bound(product_bound(square,x1),product_bound((Q(2),1),product_bound(x1,x2))),x3)]
 for q in v:mark(q)
 return v
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};rawbytes=0
 def pin(p,h=None):
  nonlocal rawbytes
  p=pathlib.Path(p);actual=sha(p)
  if h:assert actual==h,(str(p),actual,h)
  key=str(p.relative_to(ROOT))
  if key not in pins:rawbytes+=p.stat().st_size
  pins[key]=actual
  return json.loads(p.read_text()) if p.suffix=='.json' else p.read_text()
 s=pin(ROOT/'reports'/PRI/'report.json',SOURCE);assert s['all_passed'] and not s['requested_values_evaluated']
 for p,h in s['input_sha256'].items():
  assert sha(ROOT/p)==h;pins[p]=h
 c=pin(ROOT/'reports'/PRI/'compiled-eight-free-even-roots-and-pure-ell4-inputs.json',s['evidence_sha256']['compiled-eight-free-even-roots-and-pure-ell4-inputs.json'])
 prof=pin(ROOT/'reports'/PRI/'full-eta4-source-frame-and-whole-operational-profile.json',s['evidence_sha256']['full-eta4-source-frame-and-whole-operational-profile.json'])
 pin(ROOT/'code/ambient_B18_C_full_eta4_cache_source_20261105_v63.py',s['code_sha256'])
 pin(ROOT/'research/notes'/f'{PRI}.md',s['note_sha256'])
 worker=pin(ROOT/WORKER,WORKER_SHA)
 custody=pin(ROOT/'reports/AMBIENT_B18_C_FULL_ETA4_VALUE_CUSTODY_20261105_v63/report.json',
  'eb37e44ba8ed684172b5ed93c211cfd526362f18cd88a929133cecd1895080ae')
 for p,h in custody['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in custody['evidence_sha256'].items():cw=pin(ROOT/'reports/AMBIENT_B18_C_FULL_ETA4_VALUE_CUSTODY_20261105_v63'/p,h)
 theorem=pin(ROOT/'reports/AMBIENT_B18_C_FULL_CUBIC_LINE_SOURCE_20261105_v63_REVIEW/report.json',
  '57db669d19310ff27229dd025ed04d02ed5a2bd8febd843d0a13e5844e962d75')
 assert theorem['all_passed'] and theorem['eta_values_unread']
 # OriginalID matching and full ancestor closure, independently of the
 # alpha-specific or eta-specific encoder.
 original=pin(ROOT/'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63/compiled-original164-DAG-and-provided-pure-ell-inputs.json',
  'cec7274d71c0f361e1d580e6001396c4374f8ffe70359665bd9789882f440399')
 orig={r[0]:r for r in original['operator_records_originalID_op_args']}
 rec=c['operator_records_originalID_op_args'];meta=c['root_pairs_physical4']
 roots=[r[k] for r in meta for k in ['aA_node','aB_node']]
 assert [r['physical_row'] for r in meta]==[19,28,29,72]
 assert roots==[20117,40895,20132,40910,20135,40913,20159,40937]
 assert all(orig[r[0]]==r for r in rec)
 wanted=set();stack=list(roots)
 while stack:
  i=stack.pop()
  if i in wanted:continue
  wanted.add(i);r=orig[i]
  if r[1]not in ['c','input']:stack.extend(r[2:])
 assert sorted(wanted)==[r[0] for r in rec] and len(rec)==7945
 prior_inputs=original['provided105_originalID_pure_ell_coefficients'];eqs=original['equation105_order']
 own_source_peak=nc(original)+nc(c)+nc(prof)+nc(s)+nc(custody)+nc(cw)+len(orig)+len(wanted)+len(roots)+8000
 assert own_source_peak<150000
 del original,orig,wanted,stack
 joint=pin(ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
  '1909547a8d80ad59905518282c8d266d2c92f791c2c9c120b01abc3cd0860d44')
 jp=pin(ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW/report.json',
  '513e63044f8048da7c6761f1b02c63c0bcd6f643beb65339a2d3eb0de18c1401')
 assert jp['all_passed'] and joint['all105_equations_zero_through_degree3']
 assert joint['all105_equation_original_node_order']==eqs
 assert joint['monomials']==[[0,0],[1,0],[0,1],[2,0],[1,1],[0,2],[3,0],[2,1],[1,2],[0,3]]
 old=pin(ROOT/'reports/AMBIENT_B18_C3_ALPHA4_CACHED_SOURCE_20261105_v63/compiled-eight-even-roots-and-provided-joint105-plane3-cache.json',
  '83169df1e4a195a6cfef1f78024edf2cff17d64b38dd89c17a3a3e99feaabafa')
 assert joint['all105_unknown_jets']==old['provided105_originalID_r_s_coefficients10']
 assert rec==old['operator_records_originalID_op_args'] and meta==old['root_pairs_physical4']
 provided=[];lower=dict(prior_inputs)
 for i,vals in joint['all105_unknown_jets']:
  v=[Q(vals[k])/5**n for n,k in enumerate([0,2,5,9])]
  assert v[:3]==list(map(Q,lower[i]));provided.append([i,list(map(str,v))])
 assert len(provided)==105 and provided==c['provided105_originalID_pure_ell4']
 assert prior_inputs==c['provided105_originalID_pure_ell3_prior']
 ell=pin(ROOT/'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63/all164-cached-pure-ell-readbacks-and-physical12-second-values.json',
  'f213f27f92b77d734f864dbc48f49d964e47980038c3139a958608ad3647b408')
 expected=dict(ell['all164_actual_pure_ell_root_values3'])
 assert c['expected8_roots_prior_pure_ell3']==[[i,expected[i]] for i in roots]
 own_phase2=nc(c)+nc(prof)+nc(joint)+nc(old)+nc(ell)+nc(provided)+nc(lower)+16000
 own_source_peak=max(own_source_peak,own_phase2);assert own_source_peak<150000
 del joint,old,ell,lower
 # New ell-only ring; do NOT substitute the earlier r*ell² majorants.
 phi={i:[(Q(v[0]),1)]+[mark((ceiling(Q(x)),Q(x).denominator)) for x in v[1:]] for i,v in provided}
 uses=Counter()
 for r in rec:
  if r[1]not in ['c','input']:uses.update(r[2:])
 for i in roots:uses[i]+=1
 cache={};peak=0;digest=hashlib.sha256();inverses=[];constants=[]
 for r in rec:
  i,op=r[:2]
  if op=='input':v=phi[i]
  elif op=='c':v=zseries();v[0]=mark((Q(r[2],r[3]),1));constants.append([i,r[2],r[3]])
  elif op=='a':v=[sum_bound(x,y) for x,y in zip(cache[r[2]],cache[r[3]])]
  elif op=='n':v=list(cache[r[2]]);v[0]=mark((-v[0][0],1))
  elif op=='m':v=product_series(cache[r[2]],cache[r[3]])
  elif op=='i':
   v=inverse_bounds(cache[r[2]]);product_series(cache[r[2]],v);product_series(v,cache[r[2]])
   inverses.append([i,r[2],'1'])
  else:raise AssertionError(op)
  digest.update((enc([i,[[str(a),str(b)] for a,b in v]])+'\n').encode())
  cache[i]=v;peak=max(peak,len(cache))
  if op not in ['c','input']:
   for k in r[2:]:
    uses[k]-=1
    if not uses[k]:del cache[k]
  if not uses[i]:del cache[i]
 deltabounds=[]
 for row in meta:
  a,b=cache[row['aA_node']],cache[row['aB_node']]
  v=[mark((a[0][0]-b[0][0],1))]+[sum_bound(x,y) for x,y in zip(a[1:],b[1:])]
  deltabounds.append([row['physical_row'],[[str(x),str(y)] for x,y in v]])
 assert peak==184 and inverses==c['unit_inverse_centers']
 assert digest.hexdigest()==prof['every_node_majorant_sha256']
 assert deltabounds==prof['root_delta_final_majorants']
 assert constants==prof['literal_all_n_d_constants'] and len(constants)==45
 assert sum(d!=1 for _,n,d in constants)==32
 assert height==35==prof['max_stored_or_diagnostic_partial_height_majorant']
 assert 2*height+1==71==prof['max_operational_before_add_mul_bound']
 # Literal unrun worker: exact rationals n/d, guarded every QQ primitive,
 # unit inverse recurrence, both products and last-use freeing.
 tree=ast.parse(worker);fs={n.name:n for n in tree.body if isinstance(n,ast.FunctionDef)}
 seg={k:ast.get_source_segment(worker,n) for k,n in fs.items()}
 assert 'parse(str(r[2])+\'/\'+str(r[3]))' in seg['main']
 assert 'return note(a+b)' in seg['add'] and 'check(' in seg['add']
 assert 'return note(a*b)' in seg['times'] and 'check(' in seg['times']
 assert seg['add'].index('check(')<seg['add'].index('return note(a+b)')
 assert seg['times'].index('check(')<seg['times'].index('return note(a*b)')
 assert 'assert x[0]==1' in seg['inverse']
 assert 'add(v[n],neg(times(x[k],v[n-k])))' in seg['inverse']
 assert 'conv(cache[r[2]],v)==one==conv(v,cache[r[2]])' in seg['main']
 assert 'MAX_PRECHECK<=profile[' in seg['main'] and "assert delta[:3]==[Q(0)]*3" in seg['main']
 assert "if not uses[k]:del cache[k]" in seg['main']
 assert "--source-peer-report" in seg['main']
 # Worker count with this report below5000 means its literal whole bound
 # remains below75000, including every decoded metadata field and QQ buffer.
 future_base=nc(s)+nc(c)+nc(prof)+2*len(rec)+525+4*184+16+9000
 future_with_peer=future_base+5000
 assert future_with_peer<=75000 and prof['future_owned_numeric_upper']<=75000
 assert prof['source_owned_numeric_peak_bound']<150000
 assert rawbytes<2<<20
 assert not (ROOT/'reports/AMBIENT_B18_C_FULL_ETA4_CACHE_VALUES_20261105_v63').exists()
 proof={'fresh_eta03_not_old_alpha_rell2':True,
  'same_original7945_records_and_exact_ancestor_closure':True,
  'same_actual105_unknown_and_equation_namespace_unit_Jacobian_peer':True,
  'r0_s_ell_div5_ringhom_all315_lower_input_overlaps':True,
  'all24_prior_root_overlap_values_rebound_future_worker_mandatory':True,
  'all45_literal_n_d_constants_32_nonintegral_lossless':True,
  'independent_closed_positive_inverse_bounds_and_reverse_convolution':True,
  'every7945_node_majorant_hash':digest.hexdigest(),
  'all_unit_inverse_diagnostic_and_delta_partial_bounds_included':True,
  'new_max_majorant_height':height,'all_future_operation_bound':2*height+1,
  'lastuse_peak_nodes':peak,'whole_own_SOURCE_upper':own_source_peak,
  'whole_future_VALUE_upper_with5000_peer':future_with_peer,
  'prepared_worker_literal_guards_lifetimes_namespace_and_stop_reviewed':True,
  'no_requested_positive_root_or_eta_value_evaluated':True,
  'conditional_parent_ONE_job_requirements_fulfilled':True,
  'actual_scope':'only fullC HF3=5 or6 by accepted complete155; no e/H0/nil inference',
  'hard_value_slots':100000,'hard_global_slots':150000,'hard_bits':2048,
  'hard_RAW_bytes':32<<20,'hard_seconds':60,'hard_memory_bytes':2<<30}
 ep=out/'independent-ETA4-original-DAG-fresh-majorants-unrun-worker-and-whole-caps.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,
  'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},
  'source_frame_code_custody_and_lifetime_PASS':True,'requested_eta_values_evaluated':False,
  'parent_conditional_ONE_cache_job_ready':True,'elapsed_seconds':time.monotonic()-start}
 assert nc(report)<5000
 save(out/'report.json',report);signal.alarm(0)
 print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()

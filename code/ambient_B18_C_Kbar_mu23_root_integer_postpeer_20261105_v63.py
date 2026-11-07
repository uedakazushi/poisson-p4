"""Independent common-denominator integer replay of all37+14 cached dots."""
from pathlib import Path
import hashlib,json,math,resource,signal,time
R=Path(__file__).resolve().parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
out=R/'reports/AMBIENT_B18_C_KBAR_MU23_ROOT_INTEGER_POSTPEER_20261105_v63'
out.mkdir(exist_ok=False);start=time.monotonic()
resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
pins={};state={'phase':'SOURCE','cursor':{},'pending':{},'lower':[],'t37':[],'actions14':[],'direct14':[]}
maxbits=0;ops=0
def load(p,h):
 q=R/p;assert sha(q)==h;pins[p]=h;return json.loads(q.read_text())
def parse(q):
 a=str(q).split('/');assert len(a)<=2
 n=int(a[0]);d=int(a[1])if len(a)==2 else 1;assert d>0
 assert max(abs(n).bit_length(),d.bit_length())<=136
 g=math.gcd(n,d);return n//g,d//g
def text(q):return str(q[0])if q[1]==1 else str(q[0])+'/'+str(q[1])
def mul(a,b):
 global maxbits,ops
 state['pending']={'operation':'integer_multiply','a':a,'b':b}
 assert abs(a).bit_length()+abs(b).bit_length()<=273
 z=a*b;ops+=1;maxbits=max(maxbits,abs(z).bit_length());state['pending']={}
 return z
def plus(a,b):
 global maxbits
 state['pending']={'operation':'integer_add','a':a,'b':b}
 assert max(abs(a).bit_length(),abs(b).bit_length())+1<=273
 z=a+b;maxbits=max(maxbits,abs(z).bit_length());state['pending']={}
 return z
def dot(weights,values):
 terms=[]
 for j,(key,w)in enumerate(reversed(weights)):
  state['cursor']['term']=j;wn,wd=parse(w);vn,vd=values[int(key)]
  n=mul(wn,vn);d=mul(wd,vd);g=math.gcd(n,d);terms.append((n//g,d//g))
 L=1
 for n,d in terms:L=mul(L//math.gcd(L,d),d)
 total=0;state['row_common_denominator']=L;state['row_integer_partial']=0
 for n,d in terms:total=plus(total,mul(n,L//d));state['row_integer_partial']=total
 g=math.gcd(total,L);q=(total//g,L//g)
 assert max(abs(q[0]).bit_length(),q[1].bit_length())<=136
 return q
try:
 cp='reports/AMBIENT_B18_C_KBAR_MU23_CACHE_SOURCE_20261105_v63/complete-old-Pg-second-cache-source-and-fourteen-action-positive-gate.json'
 D=load(cp,'317b2c5e491e7f8369e73a7578982731415ac45d77444219308dcc5cc34b0127')
 pfx='reports/AMBIENT_B18_C_KBAR_MU23_CACHED_VALUES_20261105_v63/'
 pr=load(pfx+'report.json','0e46469e8ad8aee87e487e43d8794251b7992368c20e2531081c9756864126e2')
 assert pr['all_passed'] and pr['actual_mu_R0_K']==3
 one='phase1-all37-actual-Pg-second-values-and-lower-readbacks.json'
 two='phase2-all14-actual-actions-two-compositions-and-generator-count.json'
 v1=load(pfx+one,pr['evidence_sha256'][one]);v2=load(pfx+two,pr['evidence_sha256'][two])
 (out/'phase.json').write_text(json.dumps(dict(SOURCE_prior=True,input_sha256=pins,method='two-pass common row LCD and integer scaled sums, reversed source term order; no producer Fraction accumulator',committed136=True,temporary273=True,static_numeric14076=True,RAW800k=True))+'\n')
 C={i:{int(k):parse(c)for k,c in row}for i,row in D['all21_given_scalar_cache_rows_not_selected_rr']}
 assert len(C)==21 and D['actual_cache_monomials_original8'][19]==[7,7]
 for j,w in enumerate(D['all37_Pg_source_weights']):
  state.update(phase='independent_lower',cursor={'row':j})
  a=dot(w,{i:p.get(0,(0,1))for i,p in C.items()})
  b=dot(w,{i:p.get(5,(0,1))for i,p in C.items()})
  assert a==(0,1) and b==(0,1);state['lower'].append([j,text(a),text(b)])
 rr={i:p.get(19,(0,1))for i,p in C.items()}
 for j,w in enumerate(D['all37_Pg_source_weights']):
  state.update(phase='independent_t37',cursor={'row':j})
  q=dot(w,rr);assert text(q)==v1['all37_actual_t_head4r2'][j];state['t37'].append(q)
 t=dict(enumerate(state['t37']))
 for j,(w,dr)in enumerate(zip(D['all14_same37_weights'],D['all14_direct_Pg_source_weights'])):
  state.update(phase='independent_actions14',cursor={'row':j})
  q=dot([(i,c)for i,c in enumerate(w)],t);z=dot(dr,rr)
  assert q==z==(0,1) and text(q)==v2['all14_actual_q2_mod_rKbar'][j]
  state['actions14'].append(q);state['direct14'].append(z)
 assert len(state['lower'])==37 and len(state['t37'])==37 and len(state['actions14'])==14
 ep=out/'independent-integer-LCD-all37-and14-certificate.json'
 ep.write_text(json.dumps(state,indent=2)+'\n')
 r=dict(id=out.name,all_passed=True,independent_signed_common_denominator_integer_postpeer=True,input_sha256=pins,
  all37_origin_r_readbacks_zero=True,all37_rr_values_match=True,all14_quotient_actions_zero_in_two_compositions=True,
  original_mu_R0_K_exact=3,source_exhaustive14_action_and_ziq2_domain_proof_reused=True,
  D_module_minimal_generators_unchanged=7,original_H0_and_nil_unchanged=True,
  integer_multiplications=ops,observed_integer_temporary_bits=maxbits,hard_integer_temporary_bits=273,
  committed_prior=136,static_numeric_bound=14076,RAW_cap=800000,
  code_sha256=sha(Path(__file__)),evidence_sha256={ep.name:sha(ep),'phase.json':sha(out/'phase.json')},seconds=time.monotonic()-start)
 p=out/'report.json';p.write_text(json.dumps(r,indent=2)+'\n')
 assert sum((R/k).stat().st_size for k in pins)+sum(p.stat().st_size for p in out.iterdir())<800000
 print(sha(p),ops,maxbits)
except BaseException as e:
 (out/'failure-scientific-prefix.json').write_text(json.dumps(dict(error=repr(e),state=state,no_retry=True),indent=2)+'\n')
 raise

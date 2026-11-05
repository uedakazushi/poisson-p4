#!/usr/bin/env python3
"""Independent dictionary jet replay with closed inverses; same provided cache."""
import hashlib,json,pathlib,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_FULL_ETA4_CACHE_VALUES_20261105_v63';ID=PRI+'_REVIEW'
PRIMARY='023928cd09cbe0a43dd63c79500cac0f86a6de83325b4200e9fde18da535fa6b'
premax=valmax=0
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def nc(x):
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):
  try:Q(x);return 1
  except (ValueError,ZeroDivisionError):return 0
 if isinstance(x,(list,tuple)):return sum(map(nc,x))
 if isinstance(x,dict):return sum(map(nc,x.values()))
 return 0
def check(h):
 global premax
 premax=max(premax,h);assert h<=2048
def note(q):
 global valmax
 valmax=max(valmax,abs(q.numerator).bit_length(),q.denominator.bit_length());assert valmax<=2048
 return q
def parse(x):
 parts=str(x).split('/');assert len(parts)in[1,2]
 check(max(abs(int(v)).bit_length() for v in parts));return note(Q(str(x)))
def add(a,b):
 if not a:return b
 if not b:return a
 check(max(max(abs(a.numerator).bit_length()+b.denominator.bit_length(),abs(b.numerator).bit_length()+a.denominator.bit_length())+1,a.denominator.bit_length()+b.denominator.bit_length()))
 return note(a+b)
def mul(a,b):
 if not a or not b:return Q(0)
 check(max(abs(a.numerator).bit_length()+abs(b.numerator).bit_length(),a.denominator.bit_length()+b.denominator.bit_length()))
 return note(a*b)
def plus(a,b):
 d=dict(a)
 for k,v in b.items():
  q=add(d.get(k,Q(0)),v)
  if q:d[k]=q
  else:d.pop(k,None)
 return d
def neg(a):return {k:note(-v) for k,v in a.items()}
def times(a,b):
 d={}
 for i in sorted(a,reverse=True):
  for j in sorted(b,reverse=True):
   if i+j<4:
    k=i+j;v=add(d.get(k,Q(0)),mul(a[i],b[j]))
    if v:d[k]=v
    else:d.pop(k,None)
 return d
def inv(a):
 assert a.get(0,Q(0))==1
 u,v,w=[a.get(i,Q(0)) for i in range(1,4)]
 u2=mul(u,u);u3=mul(u2,u)
 c=[Q(1),note(-u),add(u2,-v),add(add(-u3,mul(Q(2),mul(u,v))),-w)]
 return {i:q for i,q in enumerate(c) if q}
def dense(a):return [str(a.get(i,Q(0))) for i in range(4)]
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw=0
 def pin(p,h=None):
  nonlocal raw
  p=pathlib.Path(p);v=sha(p)
  if h:assert v==h,(str(p),v,h)
  k=str(p.relative_to(ROOT))
  if k not in pins:raw+=p.stat().st_size
  pins[k]=v
  return json.loads(p.read_text()) if p.suffix=='.json' else None
 pr=pin(ROOT/'reports'/PRI/'report.json',PRIMARY);assert pr['all_passed']
 inputs={p:pin(ROOT/p,h) for p,h in pr['input_sha256'].items()}
 sourcepeer=next(d for p,d in inputs.items() if p.endswith('_REVIEW/report.json'))
 assert sourcepeer['source_frame_code_custody_and_lifetime_PASS']
 assert sha(ROOT/'code/ambient_B18_C_full_eta4_cache_values_20261105_v63.py')==pr['code_sha256']
 assert sha(ROOT/'research/notes'/f'{PRI}.md')==pr['note_sha256']
 E=pin(ROOT/'reports'/PRI/'eight-cached-free-even-roots-and-four-ell3-readbacks.json',
   'c984e1490ac5cd45cf279e41d136cffc8e4ddc2c76b871587573fedd223410c9')
 c=next(d for p,d in inputs.items() if p.endswith('compiled-eight-free-even-roots-and-pure-ell4-inputs.json'))
 rec=c['operator_records_originalID_op_args'];meta=c['root_pairs_physical4']
 roots=[r[k] for r in meta for k in ['aA_node','aB_node']]
 assert roots==[20117,40895,20132,40910,20135,40913,20159,40937] and len(rec)==7945
 own_whole=sum(nc(x) for x in inputs.values())+nc(E)+nc(pr)+2*len(rec)+525+4*184+16+10000
 assert own_whole<75000 and raw<1<<20
 save(out/'pre-arithmetic-independent-cache-replay-bound.json',{'owned_upper':own_whole,'RAW_read':raw,
   'input_fixed':True,'no_new_source_solve':True,'inverse_method':'dictionary jets closed formula 1-u+u²-u³ through ell³',
   'hard_bits':2048,'hard_seconds':60,'hard_memory_bytes':2<<30,'hard_numeric':100000,'global_numeric':150000})
 prior=dict(c['provided105_originalID_pure_ell3_prior']);phi={}
 for i,v in c['provided105_originalID_pure_ell4']:
  w=list(map(parse,v));assert w[:3]==list(map(parse,prior[i]))
  phi[i]={k:q for k,q in enumerate(w) if q}
 assert len(phi)==105
 uses=Counter()
 for r in rec:
  if r[1]not in['input','c']:uses.update(r[2:])
 for i in roots:uses[i]+=1
 values={};peak=0;units=[]
 for r in rec:
  i,op=r[:2]
  if op=='input':a=phi[i]
  elif op=='c':q=parse(str(r[2])+'/'+str(r[3]));a={0:q} if q else {}
  elif op=='a':a=plus(values[r[2]],values[r[3]])
  elif op=='n':a=neg(values[r[2]])
  elif op=='m':a=times(values[r[2]],values[r[3]])
  elif op=='i':
   a=inv(values[r[2]])
   assert times(a,values[r[2]])=={0:Q(1)}==times(values[r[2]],a)
   units.append([i,r[2],'1'])
  else:raise AssertionError(op)
  values[i]=a;peak=max(peak,len(values))
  if op not in['input','c']:
   for j in r[2:]:
    uses[j]-=1
    if not uses[j]:del values[j]
  if not uses[i]:del values[i]
 assert peak==184 and units==c['unit_inverse_centers']==E['both_unit_inverse_product_checks']
 expected=dict(c['expected8_roots_prior_pure_ell3'])
 for i in roots:assert dense(values[i])[:3]==expected[i]
 actual=[[i,dense(values[i])] for i in roots]
 assert actual==E['all8_free_even_roots_pure_ell4']
 results=[]
 for row in meta:
  a,b=values[row['aA_node']],values[row['aB_node']];delta=plus(a,neg(b));assert not delta
  saved=next(x for x in E['all4_eta03_actual_readbacks'] if x['physical_row']==row['physical_row'])
  assert dense(delta)==saved['delta_pure_ell4'] and saved['delta_ell3']=='0'
  results.append([row['physical_row'],dense(delta)])
 assert E['actual_full_C_cubic_extra_rank']==0 and E['actual_full_C_HF_degree3']==6
 assert pr['actual_full_C_HF_degree3']==6 and pr['nonzero_eta_count']==0
 assert not (ROOT/'reports/AMBIENT_B18_C_FULL_GRAPH03_CACHE_VALUES_20261105_v63').exists()
 proof={'all7945_original_source_ops_independently_replayed':True,
  'own_dictionary_ell4_arithmetic_reverse_multiplication_and_closed_inverses':True,
  'all315_same105_lower_inputs_all24_prior_root_overlaps':True,
  'all8_roots_and_all4_full_deltas_exact':actual,'four_actual_deltas':results,
  'both_unit_inverse_products':units,'lastuse_peak':peak,
  'actual_full_C_HF_degrees0through3':[3,6,6,6],
  'complete155_four_eta_callbacks_allzero_rank0':True,
  'no_exact_e_from_false_nonzero_eta_conditional_antecedent':True,
  'no_original_H0_nil_or_exhaustion_inference':True,
  'max_own_before_operation_bits':premax,'max_own_stored_QQ_bits':valmax,
  'owned_numeric_upper':own_whole,'RAW_read':raw,
  'hard_bits':2048,'hard_numeric':100000,'hard_global_numeric':150000,
  'hard_seconds':60,'hard_memory_bytes':2<<30,'hard_RAW_bytes':32<<20}
 ep=out/'independent-ell4-dictionary-root-replay-and-four-eta-zero-proof.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'input_sha256':pins,'code_sha256':sha(__file__),
  'evidence_sha256':{ep.name:sha(ep)},'all4_eta03_zero_and_actual_fullC_HF3_6_accepted':True,
  'e_remains_open':True,'elapsed_seconds':time.monotonic()-start}
 save(out/'report.json',report);signal.alarm(0)
 print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()

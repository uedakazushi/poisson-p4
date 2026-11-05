#!/usr/bin/env python3
"""One gated existing-cache graph ell^3 scan in QQ[ell]/ell^4."""
import argparse,hashlib,json,re,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_FULL_GRAPH03_CACHE_VALUES_20261105_v63'
SOURCE='reports/AMBIENT_B18_C_FULL_GRAPH03_CACHE_SOURCE_20261105_v63/report.json'
SOURCE_SHA='bc904617d476ff6bca01deef26691926d38decd7bc66176516d95d1af7a3af52'
COMPILED='reports/AMBIENT_B18_C_FULL_GRAPH03_CACHE_SOURCE_20261105_v63/compiled-ten-graph-even-roots-and-pure-ell4-inputs.json'
PROFILE='reports/AMBIENT_B18_C_FULL_GRAPH03_CACHE_SOURCE_20261105_v63/full-graph03-source-frame-and-whole-operational-profile.json'
MAX_PRECHECK=MAX_VALUE=0

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.values()))
 return 0
def check(n):
 global MAX_PRECHECK
 MAX_PRECHECK=max(MAX_PRECHECK,n);assert n<=2048
def note(q):
 global MAX_VALUE
 MAX_VALUE=max(MAX_VALUE,abs(q.numerator).bit_length(),q.denominator.bit_length());assert MAX_VALUE<=2048
 return q
def parse(x):
 p=str(x).split('/');assert len(p)in[1,2]
 check(max(abs(int(v)).bit_length()for v in p));return note(Q(str(x)))
def add(a,b):
 if not a:return b
 if not b:return a
 check(max(max(abs(a.numerator).bit_length()+b.denominator.bit_length(),abs(b.numerator).bit_length()+a.denominator.bit_length())+1,a.denominator.bit_length()+b.denominator.bit_length()))
 return note(a+b)
def times(a,b):
 if not a or not b:return Q(0)
 check(max(abs(a.numerator).bit_length()+abs(b.numerator).bit_length(),a.denominator.bit_length()+b.denominator.bit_length()))
 return note(a*b)
def neg(a):return note(-a)
def conv(a,b):
 v=[Q(0)]*4
 for n in range(4):
  for k in range(n+1):v[n]=add(v[n],times(a[k],b[n-k]))
 return v
def inverse(x):
 assert x[0]==1
 v=[Q(1),Q(0),Q(0),Q(0)]
 for n in range(1,4):
  for k in range(1,n+1):v[n]=add(v[n],neg(times(x[k],v[n-k])))
 return v

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source-peer-report',required=True);args=parser.parse_args()
 signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2147483648,2147483648));started=time.monotonic()
 assert sha(ROOT/SOURCE)==SOURCE_SHA
 source=json.loads((ROOT/SOURCE).read_text());peerpath=Path(args.source_peer_report)
 if not peerpath.is_absolute():peerpath=ROOT/peerpath
 peer=json.loads(peerpath.read_text());assert peer['all_passed']and SOURCE_SHA in peer['input_sha256'].values()
 for p in[COMPILED,PROFILE]:assert sha(ROOT/p)==source['evidence_sha256'][Path(p).name]
 c=json.loads((ROOT/COMPILED).read_text());profile=json.loads((ROOT/PROFILE).read_text())
 records=c['operator_records_originalID_op_args'];meta=c['root_pairs_physical5']
 assert len(records)==9823 and[r['physical_row']for r in meta]==[63,64,68,73,74]
 roots=[r[n]for r in meta for n in['aA_node','aB_node']]
 assert roots==[20144,40922,20147,40925,20156,40934,20162,40940,20165,40943]
 owned=count(source)+count(peer)+count(c)+count(profile)+2*len(records)+525+4*198+16+9000
 assert owned<=75000<100000<150000
 raw=sum(p.stat().st_size for p in[ROOT/SOURCE,ROOT/COMPILED,ROOT/PROFILE,peerpath]);assert raw+150000<=1000000
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 (out/'pre-arithmetic-source-and-whole-budget-custody.json').write_text(json.dumps({
  'source_sha256':SOURCE_SHA,'source_peer_sha256':sha(peerpath),'code_sha256':sha(Path(__file__)),
  'owned_numeric_upper':owned,'RAW_input_bytes':raw,'hard_bits':2048,'hard_seconds':60,
  'requested_positive_root_values_not_yet_evaluated':True},indent=2)+'\n')
 phi={i:list(map(parse,v))for i,v in c['provided105_originalID_pure_ell4']};prior=dict(c['provided105_originalID_pure_ell3_prior'])
 assert len(phi)==105 and all(v[:3]==list(map(parse,prior[i]))for i,v in phi.items())
 uses=Counter()
 for r in records:
  if r[1]not in['c','input']:uses.update(r[2:])
 for i in roots:uses[i]+=1
 cache={};peak=0;inverses=[]
 try:
  for r in records:
   i,op=r[:2]
   if op=='input':v=phi[i]
   elif op=='c':v=[parse(str(r[2])+'/'+str(r[3]))]+[Q(0)]*3
   elif op=='a':v=[add(a,b)for a,b in zip(cache[r[2]],cache[r[3]])]
   elif op=='n':v=[neg(a)for a in cache[r[2]]]
   elif op=='m':v=conv(cache[r[2]],cache[r[3]])
   elif op=='i':
    v=inverse(cache[r[2]]);one=[Q(1),Q(0),Q(0),Q(0)]
    assert conv(cache[r[2]],v)==one==conv(v,cache[r[2]]);inverses.append([i,r[2],'1'])
   else:raise AssertionError(op)
   cache[i]=v;peak=max(peak,len(cache))
   if op not in['c','input']:
    for k in r[2:]:
     uses[k]-=1
     if not uses[k]:del cache[k]
   if not uses[i]:del cache[i]
  assert peak==198 and inverses==c['unit_inverse_centers']
  expected=dict(c['expected10_roots_prior_pure_ell3'])
  assert all(cache[i][:3]==list(map(parse,expected[i]))for i in roots)
  results=[]
  for row in meta:
   a,b=cache[row['aA_node']],cache[row['aB_node']];delta=[add(x,neg(y))for x,y in zip(a,b)]
   assert delta[:3]==[Q(0)]*3
   results.append({'physical_row':row['physical_row'],'aA_node':row['aA_node'],'aB_node':row['aB_node'],
    'aA_pure_ell4':list(map(str,a)),'aB_pure_ell4':list(map(str,b)),
    'delta_pure_ell4':list(map(str,delta)),'delta_ell3':str(delta[3])})
 except BaseException as err:
  (out/'ROOT_ONLY_STOP_CHECKPOINT.json').write_text(json.dumps({'status':'STOP_no_cap_raise','error':repr(err),'completed_original_node':i,
   'max_value_bits':MAX_VALUE,'max_precheck_bits':MAX_PRECHECK,'already_available_root_values':[[j,list(map(str,cache[j]))]for j in roots if j in cache]},indent=2)+'\n');raise
 assert MAX_PRECHECK<=profile['max_operational_before_add_mul_bound']
 nonzero=[r['physical_row']for r in results if parse(r['delta_ell3'])]
 data={'all10_graph_even_roots_pure_ell4':[[i,list(map(str,cache[i]))]for i in roots],
  'all5_graph03_actual_readbacks':results,'nonzero_rows':nonzero,
  'all315_same105_lower_input_overlaps':True,'all30_prior_root_coefficients':True,'both_unit_inverse_product_checks':inverses,
  'lastuse_peak_nodes':peak,'max_actual_QQ_bits':MAX_VALUE,'max_before_operation_bits':MAX_PRECHECK,
  'source_before_operation_majorant':profile['max_operational_before_add_mul_bound'],
  'owned_numeric_upper':owned,'RAW_input_bytes':raw,'new_physical_IFT_gauge_or_jet_coefficients':0,
  'full_C_exponent3_theorem_triggered':bool(nonzero),
  'full_C_DVR_torsion_exponents_if_nonzero':[1,1,2,3]if nonzero else'OPEN_after_this_subset',
  'full_C_HF_if_nonzero':[3,6,6,'4 at every degree>=3']if nonzero else'UNDETERMINED_by_this_subset',
  'comparison_C_torsion_length_if_nonzero':7 if nonzero else None,
  'no_original_H0_length_nilcube_or_exact_nilindex_transfer':True,
  'allzero_scope':'Five graph03 witnesses negative only; no full C HF/exponent conclusion.'}
 ep=out/'ten-cached-graph-even-roots-and-five-ell3-readbacks.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'NONZERO_GRAPH03_FULL_C_EXPONENT3'if nonzero else'FIVE_GRAPH03_WITNESSES_ZERO_ONLY',
  'nonzero_graph_count':len(nonzero),'full_C_exact_torsion1123_proved':bool(nonzero),
  'input_sha256':{SOURCE:SOURCE_SHA,COMPILED:sha(ROOT/COMPILED),PROFILE:sha(ROOT/PROFILE),str(peerpath.relative_to(ROOT)):sha(peerpath)},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 assert raw+sum(p.stat().st_size for p in out.iterdir())<=1000000
 print(json.dumps(report));signal.alarm(0)

if __name__=='__main__':main()

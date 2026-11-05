#!/usr/bin/env python3
"""Gated ten-root composition of the provided joint105 cache, no new solve."""
import argparse
import hashlib
import json
import re
import resource
import signal
import time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C3_GRAPH5_CACHED_VALUES_20261105_v63'
SOURCE='reports/AMBIENT_B18_C3_GRAPH5_CACHED_SOURCE_20261105_v63/report.json'
SOURCE_SHA='84de2029149eee7ddc187004d2f356ab3675f044e4a2198c3f75b3e56672a8df'
COMPILED='reports/AMBIENT_B18_C3_GRAPH5_CACHED_SOURCE_20261105_v63/compiled-ten-even-roots-and-provided-joint105-plane3-cache.json'
PROFILE='reports/AMBIENT_B18_C3_GRAPH5_CACHED_SOURCE_20261105_v63/graph5-source-carry-cache-binding-and-operational-height-profile.json'
Q20='reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63/phase1-common105-and-114-source-equations.json'
MON=[(0,0),(1,0),(0,1),(1,1),(0,2),(1,2)]
MI={m:i for i,m in enumerate(MON)}
PAIRS=[[(a,b)for a,u in enumerate(MON)for b,v in enumerate(MON)
        if(u[0]+v[0],u[1]+v[1])==m]for m in MON]
MAX_PRECHECK=0
MAX_VALUE=0


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.values()))
 return 0
def field(text,key):
 token=json.dumps(key);assert text.count(token)==1;i=text.index(':',text.index(token))+1
 while text[i].isspace():i+=1
 return json.JSONDecoder().raw_decode(text,i)[0]
def checked_height(n):
 global MAX_PRECHECK
 MAX_PRECHECK=max(MAX_PRECHECK,n);assert n<=2048
def note(q):
 global MAX_VALUE
 MAX_VALUE=max(MAX_VALUE,abs(q.numerator).bit_length(),q.denominator.bit_length())
 assert MAX_VALUE<=2048
 return q
def parse(s):
 p=str(s).split('/');checked_height(max(abs(int(v)).bit_length()for v in p))
 return note(Q(s))
def add(a,b):
 if not a:return b
 if not b:return a
 checked_height(max(max(abs(a.numerator).bit_length()+b.denominator.bit_length(),
                         abs(b.numerator).bit_length()+a.denominator.bit_length())+1,
                    a.denominator.bit_length()+b.denominator.bit_length()))
 return note(a+b)
def times(a,b):
 if not a or not b:return Q(0)
 checked_height(max(abs(a.numerator).bit_length()+abs(b.numerator).bit_length(),
                    a.denominator.bit_length()+b.denominator.bit_length()))
 return note(a*b)
def neg(a):return note(-a)
def conv(a,b):
 v=[Q(0)]*6
 for m,pairs in enumerate(PAIRS):
  for i,j in pairs:v[m]=add(v[m],times(a[i],b[j]))
 return v
def inverse(x):
 assert x[0]==1
 v=[Q(1)]+[Q(0)]*5
 for m in range(1,6):
  for a,b in PAIRS[m]:
   if a:v[m]=add(v[m],neg(times(x[a],v[b])))
 return v


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--source-peer-report',required=True)
 args=parser.parse_args();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2147483648,2147483648))
 started=time.monotonic();assert sha(ROOT/SOURCE)==SOURCE_SHA
 source=json.loads((ROOT/SOURCE).read_text())
 peer_path=Path(args.source_peer_report)
 if not peer_path.is_absolute():peer_path=ROOT/peer_path
 peer=json.loads(peer_path.read_text());assert peer['all_passed']
 assert SOURCE_SHA in set(peer.get('input_sha256',{}).values())
 for p in[COMPILED,PROFILE]:assert sha(ROOT/p)==source['evidence_sha256'][Path(p).name]
 assert sha(ROOT/Q20)==source['input_sha256'][Q20]
 c=json.loads((ROOT/COMPILED).read_text());profile=json.loads((ROOT/PROFILE).read_text())
 text=(ROOT/Q20).read_text();qm=field(text,'base27_monomials');qvals=dict(field(text,'all105_cached_composite_values'));del text
 qi=[qm.index(m)for m in[[],[7],[8],[7,7],[7,8],[8,8]]]
 records=c['operator_records_originalID_op_args'];assert len(records)==9823
 input_count=count(source)+count(peer)+count(c)+count(profile)
 # Includes the largest temporary selected Q20 field (2940), converted
 # inputs, dependency counts, every last-use buffer, output copies and
 # serialization/controller metadata. Dense None controls are not scalars.
 owned_bound=(input_count+2940+735+9823+1200+20000)
 assert owned_bound<=75000<100000
 raw_inputs=sum(p.stat().st_size for p in[ROOT/SOURCE,ROOT/COMPILED,ROOT/PROFILE,ROOT/Q20,peer_path])
 assert raw_inputs+150000<=1000000
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 (out/'pre-arithmetic-source-and-budget-custody.json').write_text(json.dumps({
     'source_report_sha256':SOURCE_SHA,'worker_sha256':sha(Path(__file__)),
     'source_peer_sha256':sha(peer_path),'owned_numeric_bound':owned_bound,
     'raw_input_bytes':raw_inputs,'hard_bits':2048,'hard_seconds':60,
     'root_values_not_yet_evaluated':True},indent=2)+'\n')
 phi={};overlap_count=0
 for i,v in c['provided105_originalID_r_s_coefficients10']:
  a=list(map(parse,v));lower=[a[0],add(a[1],a[2]),times(a[2],parse('1/5')),
      add(add(a[3],a[4]),a[5]),times(add(a[4],times(parse('2'),a[5])),parse('1/5')),
      times(a[5],parse('1/25'))]
  assert lower==[parse(qvals[i][k])for k in qi];overlap_count+=6
  # Ordinary Taylor coefficients under s=r+ell/5; B6 discards r^2.
  phi[i]=[lower[0],lower[1],lower[2],lower[4],lower[5],
          times(add(a[8],times(parse('3'),a[9])),parse('1/25'))]
 assert overlap_count==630 and set(phi)==set(qvals)
 del qvals,qm
 max_id=max(r[0]for r in records)
 uses=[None]*(max_id+1);cache=[None]*(max_id+1)
 for r in records:uses[r[0]]=0
 roots=[r[n]for r in c['root_pairs_physical5']for n in['aA_node','aB_node']]
 for r in records:
  if r[1]not in['input','c']:
   for j in r[2:]:assert uses[j]is not None;uses[j]+=1
 for i in roots:uses[i]+=1
 live=peak=0;unit_checks=[]
 try:
  for r in records:
   i,op=r[:2]
   if op=='input':v=phi[i]
   elif op=='c':v=[parse(r[2])]+[Q(0)]*5
   elif op=='a':v=[add(a,b)for a,b in zip(cache[r[2]],cache[r[3]])]
   elif op=='n':v=[neg(a)for a in cache[r[2]]]
   elif op=='m':v=conv(cache[r[2]],cache[r[3]])
   elif op=='i':
    v=inverse(cache[r[2]]);one=[Q(1)]+[Q(0)]*5
    assert conv(v,cache[r[2]])==one==conv(cache[r[2]],v)
    unit_checks.append([i,r[2],True])
   else:raise AssertionError(op)
   cache[i]=v;live+=1;peak=max(peak,live)
   if op not in['input','c']:
    for j in r[2:]:
     uses[j]-=1
     if not uses[j]:cache[j]=None;live-=1
   if not uses[i]:cache[i]=None;live-=1
 except BaseException as err:
  checkpoint={'status':'STOP_no_cap_raise','error':repr(err),'completed_original_node':i,
      'known_root_values':[[n,list(map(str,cache[n]))]for n in roots if cache[n]is not None],
      'max_precheck_bits':MAX_PRECHECK,'max_value_bits':MAX_VALUE}
  (out/'ROOT_ONLY_STOP_CHECKPOINT.json').write_text(json.dumps(checkpoint,indent=2)+'\n')
  raise
 assert peak==198 and unit_checks==[[a,b,True]for a,b,z in c['both_original_inverse_centers']]
 expected=dict(c['expected10_lower_B6_coefficients_order_1_r_ell_rell_ell2'])
 for i in roots:assert cache[i][:5]==list(map(parse,expected[i]))
 actual=[]
 for r in c['root_pairs_physical5']:
  aa,ab=cache[r['aA_node']],cache[r['aB_node']]
  delta=[add(x,neg(y))for x,y in zip(aa,ab)]
  assert delta[:5]==[Q(0)]*5
  actual.append({'physical_row':r['physical_row'],'aA_node':r['aA_node'],'aB_node':r['aB_node'],
      'aA_B6':list(map(str,aa)),'aB_B6':list(map(str,ab)),
      'delta_B6':list(map(str,delta)),'delta_r_ell2':str(delta[5])})
 nonzero=[r['physical_row']for r in actual if parse(r['delta_r_ell2'])]
 data={'all10_even_roots_B6':[[i,list(map(str,cache[i]))]for i in roots],
       'all5_actual_graph_delta_B6':actual,
       'nonzero_rows':nonzero,'all630_same105_lower_cache_overlaps_verified':True,
       'all50_root_prior_lower_overlaps_verified':True,'both_unit_inverse_products':unit_checks,
       'source105_coefficients_provided_not_solved':True,'lastuse_peak_nodes':peak,
       'max_actual_value_bits':MAX_VALUE,'max_operational_precheck_bits':MAX_PRECHECK,
       'source_majorant_operation_bound':profile['max_future_QQ_temporary_precheck_bound'],
       'owned_numeric_bound_including_input_and_output_metadata':owned_bound,
       'input_RAW_bytes':raw_inputs,'new_physical_gauge_IFT_or_jet_coefficients':0,
       'known_C3_upper19_from_nonzero_cubic_and_quartic_multiple':bool(nonzero),
       'actual_C3_length18_or_full_torsion1122_determined':False,
       'allzero_scope':'Only this five-graph cubic subset; no full cubic/quartic source vanishing inferred.'}
 assert MAX_PRECHECK<=profile['max_future_QQ_temporary_precheck_bound']
 ep=out/'ten-cached-even-root-values-and-five-graph-cubic-readbacks.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'POSITIVE_C3_UPPER19'if nonzero else'ALLZERO_FIVE_GRAPH_CUBICS_SUBSET_ONLY',
         'nonzero_count':len(nonzero),'actual_C3_upper19_proved':bool(nonzero),
         'actual_C3_length_or_all_torsion_exponents_determined':False,
         'input_sha256':{SOURCE:SOURCE_SHA,COMPILED:sha(ROOT/COMPILED),PROFILE:sha(ROOT/PROFILE),Q20:sha(ROOT/Q20),str(peer_path.relative_to(ROOT)):sha(peer_path)},
         'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
         'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 assert raw_inputs+sum(p.stat().st_size for p in out.iterdir())<=1000000
 print(json.dumps(report));signal.alarm(0)


if __name__=='__main__':main()

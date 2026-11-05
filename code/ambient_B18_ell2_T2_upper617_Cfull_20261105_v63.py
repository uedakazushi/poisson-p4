#!/usr/bin/env python3
"""Prepared ONE cache-only224-Cfull construction. No K2 or upper images."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_T2_UPPER617_CFULL_20261105_v63'
GATE='reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_INPUT_GATE_20261105_v63/'
GATE_SHA='7ac51097df48c4416f377fa94380f3abb910f6f10bae77945f96e93ad76b30b4'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(v):return {int(i):F(x) for i,x in v if F(x)}
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except (ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))

def main(out,peer,custody):
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(exist_ok=False,parents=True);pins={};cursor={};maximum=0;C=[];XP=[]
 def read(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def dump(name,v):
  (out/name).write_text(json.dumps(v,separators=(',',':'))+'\n')
  assert sum(p.stat().st_size for p in out.iterdir())<1_000_000
 g=read(GATE+'report.json');assert g['all_passed'] and pins[GATE+'report.json']==GATE_SHA
 ip=GATE+'pruned-known-Cfull224-input-and-positive-source-envelopes.json';d=read(ip)
 assert pins[ip]==g['evidence_sha256'][Path(ip).name]
 p=read(peer);c=read(custody)
 assert p['all_passed'] and GATE_SHA in json.dumps(p)
 assert 'b4d00c49577561f7fddb4e928ca905adbc40699703a96dbfd9e204071cec9edb' in json.dumps(p)
 assert c['all_passed'] and c['safe_help_literal_PASS'] and c['worker_sha256']==sha(Path(__file__))
 for path,h in c['literal_input_sha256'].items():assert sha(ROOT/path)==h
 assert count([g,p,c,pins])<1000 and count(d)==g['pruned_numeric_records']==2401
 prior=g['source_positive_prior_bits'];cap=g['future_Cfull_numeric_bound']
 assert prior==147 and cap==11644 and c['derived_committed_bits']==prior
 assert c['derived_temporary_bits']==3*prior+1==442
 assert d['native_pivot_ell_is_SUM_L0Se_LeS0'] and d['full_coordinate_not_pureU']
 G=list(map(parse,d['Gamma_two']));L0=list(map(parse,d['native_L0']));Le=list(map(parse,d['native_Le']));S0=list(map(parse,d['native_S0']));Se=list(map(parse,d['native_Se']))
 P=list(map(parse,d['P27_columns56']));Pi=list(map(parse,d['Pi3_columns56']))
 Xr=list(map(parse,d['pivot27_r_columns']));Xe=[]
 delta=[(parse(v['r_row56']),parse(v['ell_row56'])) for v in d['DeltaC_known_point_rows']]
 alpha=list(map(F,d['alpha_three']));beta=list(map(F,d['beta_three']))
 assert len(P)==len(Pi)==56 and len(G)==2 and len(Xr)==27
 # Support unions cost no rational response arithmetic. Charge every future
 # Xe entry twice(index,value), all108 XP and224 C values, their332 serialized
 # copies,600 native-product scratch records and1000 metadata/frontier records.
 support_bound=0
 for kind,arg in d['pivot27_ell_recipes']:
  if kind=='literal':support_bound+=len(arg)
  else:
   assert kind=='native'
   support_bound+=len({i for j in Se[arg] for i in L0[j]}|{i for j in S0[arg] for i in Le[j]})
 before=2*count(d)+count([G,L0,Le,S0,Se,P,Pi,Xr,delta,alpha,beta])+2*support_bound+108+224+332+600+1000
 assert before<cap
 dump('phase.json',{'phase':'literal_source_peer_custody_and_preoperation_bounds_PASS_before_signed_C','gate_sha256':GATE_SHA,'source_committed_bits':147,'integer_temporary_bits':442,'hard_guards':[384,1153],'whole_numeric_cap':cap,'literal_source_numeric':count(d),'prearith_live_upper':before,'RAW_output_cap':1000000,'no_K_X2_X3_or_relation_query':True})
 def stop(reason,x):
  state={'reason':reason,'cursor':cursor,'uncommitted_fraction':str(x),'completed_C_prefix':[[[str(x) for x in row] for row in ch] for ch in C],'input_sha256':pins,'no_retry_or_cap_raise':True}
  assert 2*count(d)+count([state,C,XP,G,L0,Le,S0,Se,P,Pi,Xr,Xe,delta,alpha,beta])+1000<cap
  dump('guard-frontier.json',state);raise RuntimeError(reason)
 def ck(x):
  nonlocal maximum
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length());maximum=max(maximum,b)
  if b>prior:stop('source147_before_commit',x)
  return x
 def add(a,b):return ck(a+b)
 def mul(a,b):return ck(a*b)
 def mv(cols,x):
  y={}
  for j,a in x.items():
   for i,b in cols[j].items():y[i]=add(y.get(i,F(0)),mul(a,b))
  return {i:x for i,x in y.items() if x}
 def pluscols(a,b):
  y=dict(a)
  for i,x in b.items():y[i]=add(y.get(i,F(0)),x)
  return {i:x for i,x in y.items() if x}
 def dot(row,col):
  q=F(0)
  for i,x in row.items():q=add(q,mul(x,col.get(i,F(0))))
  return q
 for n,(kind,arg) in enumerate(d['pivot27_ell_recipes']):
  cursor={'phase':'individual_actual_pivot_ell_SUM','pivot27':n}
  if kind=='literal':Xe.append(parse(arg))
  else:
   assert kind=='native';Xe.append(pluscols(mv(L0,Se[arg]),mv(Le,S0[arg])))
 XP=[[[dot(q,x) for x in Xr],[dot(q,x) for x in Xe]] for q in G]
 for t in range(2):
  ch=[];C.append(ch)
  for k in range(2):
   row=[];ch.append(row)
   for j in range(56):
    cursor={'phase':'known_Cfull224_signed_cache_construction','channel':t,'r_or_ell':k,'column56':j}
    q=delta[t][k].get(j,F(0))
    if t==0:
     weights=alpha if k==0 else beta
     for i,x in Pi[j].items():q=add(q,mul(weights[i],x))
    for i,x in P[j].items():q=add(q,-mul(XP[t][k][i],x))
    D,U=map(int,d['positive_Cfull_two_by2_by56'][t][k][j])
    assert D>0 and abs(q.numerator)*D<=U*q.denominator and D%q.denominator==0
    row.append(ck(q))
 ep='actual224-Cfull-coefficients-with-native-pivot-SUM-and-Psi.json'
 evidence={'Cfull_two_by2_by56':[[[str(x) for x in row] for row in ch] for ch in C],'Gamma_Xp1_two_by2_by27':[[[str(x) for x in row] for row in ch] for ch in XP],'all224_positive_source_envelopes_verified':True,'full_Psi_bottom_correction_included':True,'individual_native_ell_includes_L0Se_PLUS_LeS0':True,'source_formula':d['source_formula'],'coherent_existence_frame_only_not_canonical_complementary_M1_values':True,'no_actual_K_Xstar2_Xc3_or_upper_relation_response':True}
 owned=2*count(d)+count([G,L0,Le,S0,Se,P,Pi,Xr,Xe,delta,alpha,beta,XP,C,evidence,g,p,c,pins,cursor])+1000
 assert owned<cap
 dump(ep,evidence)
 report={'id':STEM,'all_passed':True,'known_Cfull224_cache_construction_only':True,'all224_positive_source_bounds_checked':True,'actual_Cfull_maximum_bits':maximum,'derived_committed_bits':147,'derived_temporary_bits':442,'hard_guards':[384,1153],'whole_numeric_bound':cap,'actual_owned_numeric_upper':owned,'RAW_cap':1000000,'no_K_upper_HF_or_T2length_value':True,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 dump('report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','actual_Cfull_maximum_bits','actual_owned_numeric_upper','seconds']}))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output,a.source_peer,a.custody)

#!/usr/bin/env python3
"""Prepared ONE existing-cache396 Kfull values; no b1/upper responses."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_K99_20261105_v63'
GATE='reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_CACHE_INPUT_GATE_V2_20261105_v63/'
GATE_SHA='d3bf63533e2ab2a2882a51d09639ccba98d9d966d00edffd73a3f80e54c26bc1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except (ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def parse(v):return {int(i):F(x) for i,x in v if F(x)}
def rec(v):return [[i,str(x)] for i,x in sorted(v.items()) if x]
def main(out,peer,custody):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True)
 pins={};cursor={};K=[];pieces=[];maximum=0;jets={};oldg={};R=[];P=[];Q=[]
 def read(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def dump(name,v):
  (out/name).write_text(json.dumps(v,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in out.iterdir())<2_000_000
 g=read(GATE+'report.json');assert g['all_passed'] and pins[GATE+'report.json']==GATE_SHA
 ip=GATE+'pruned99-Kfull-P5-native23-and-prior396-input.json';d=read(ip);assert pins[ip]==g['evidence_sha256'][Path(ip).name]
 p=read(peer);c=read(custody)
 assert p['all_passed'] and GATE_SHA in json.dumps(p) and c['all_passed'] and c['safe_help_literal_PASS'] and c['worker_sha256']==sha(Path(__file__))
 for path,h in c['literal_input_sha256'].items():assert sha(ROOT/path)==h
 assert count([g,p,c,pins])<1000 and count(d)==4458
 prior=g['positive_committed_prior_bits'];cap=g['future_value_numeric_bound'];assert prior==147 and 3*prior+1==g['temporary_prior_bits']==442 and cap==23302
 assert g['all128_B11_S0_Sr_matches'] and g['all_Lr_Le_C16_matches'] and d['all99_weighted_nativeS0zero_source_order_peer_pinned']
 G=list(map(parse,d['Gamma_two']));C=[[[F(x) for x in row] for row in ch] for ch in d['Cfull_two_by2_by56']]
 L0=list(map(parse,d['native_L0']));Lr=list(map(parse,d['native_Lr']));Le=list(map(parse,d['native_Le']))
 Sr={j:parse(v) for j,v in d['native23_Sr']};Se={j:parse(v) for j,v in d['native23_Se']};Srr={j:parse(v) for j,v in d['native23_S2rr']};Sre={j:parse(v) for j,v in d['native23_S2rell']}
 oldr={j:parse(v) for j,v in d['old80_r_first_columns']};olde={j:parse(v) for j,v in d['old80_ell_first_columns']};leaf={i:list(map(F,v)) for i,v in d['physical24_plane5_values']}
 source_parsed=count([G,C,L0,Lr,Le,Sr,Se,Srr,Sre,oldr,olde,leaf])
 # Explicit worst case:288 shared-row(index,value)entries twice=1152;
 # 50 jets+IDs300;80 old-column
 # scalar arrays+IDs400;192 weighted native entries charged index+value384;
 # currentold vectors200;all396values+1980pieces twice4752;metadata1000.
 liveprior=2*count(d)+source_parsed+1152+300+400+384+200+4752+1000;assert liveprior<cap
 dump('phase.json',{'phase':'literal_source_peer_custody_and_livebounds_PASS_before_any_signed_K','gate_sha256':GATE_SHA,'explicit_live_prior':liveprior,'source_input_numeric4458_twice':8916,'source_committed_bits':147,'integer_temporary_bits':442,'hard_guards':[384,1153],'whole_numeric_cap':cap,'RAW_output_cap':2000000,'no_b1_full918_or_upper_response_query':True})
 def stop(reason,x):
  state={'reason':reason,'cursor':cursor,'uncommitted_fraction':str(x),'completed_K_prefix':K,'completed_piece_prefix':pieces,'input_sha256':pins,'no_retry_or_cap_raise':True};assert 2*count(d)+2*count(state)+count([G,C,L0,Lr,Le,Sr,Se,Srr,Sre,oldr,olde,leaf,jets,oldg,R,P,Q])+1000<cap;dump('guard-frontier.json',state);raise RuntimeError(reason)
 def ck(x):
  nonlocal maximum
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length());maximum=max(maximum,b)
  if b>prior:stop('source147_before_commit',x)
  return x
 def add(a,b):return ck(a+b)
 def mul(a,b):return ck(a*b)
 def dot(row,col):
  q=F(0)
  for i,x in row.items():q=add(q,mul(x,col.get(i,F(0))))
  return q
 def left(row,cols):return {j:dot(row,col) for j,col in enumerate(cols) if dot(row,col)}
 def plus(a,b):
  y=dict(a)
  for i,x in b.items():y[i]=add(y.get(i,F(0)),x)
  return {i:x for i,x in y.items() if x}
 for t in range(2):
  cursor={'phase':'known_shared_native_rows','channel':t}
  cr={i:x for i,x in enumerate(C[t][0]) if x};ce={i:x for i,x in enumerate(C[t][1]) if x}
  R.append(left(G[t],L0));P.append(plus(left(G[t],Lr),left(cr,L0)));Q.append(plus(left(G[t],Le),left(ce,L0)))
 mon=[(0,0),(1,0),(0,1),(2,0),(1,1)];idx={m:i for i,m in enumerate(mon)};table={(i,j):idx[(a[0]+b[0],a[1]+b[1])] for i,a in enumerate(mon) for j,b in enumerate(mon) if (a[0]+b[0],a[1]+b[1]) in idx}
 def jadd(a,b):return [add(x,y) for x,y in zip(a,b)]
 def jscale(a,c):return [mul(x,c) for x in a]
 def jmul(a,b):
  q=[F(0)]*5
  for (i,j),k in table.items():q[k]=add(q[k],mul(a[i],b[j]))
  return q
 skip=set(d['source_order_skip_seven_products']);honest=d['honest327554_opaque_typed_leaf'];gb=d['gB_plane_order3_typed_leaf'];assert len(skip)==7
 for i,n in d['typed_program50']:
  cursor={'phase':'OLD50_typed_P5_cache','typed_node':i}
  if i==honest:jets[i]=None
  elif i in skip or i==gb:jets[i]=[F(0)]*5
  elif n[0]=='actual_C0_scalar_function':assert n[1] in leaf;jets[i]=leaf[n[1]]
  elif n[0]=='sum':
   v=[F(0)]*5
   for k,a,b in n[1]:assert jets[k] is not None;v=jadd(v,jscale(jets[k],F(a,b)))
   jets[i]=v
  elif n[0]=='product':
   v=[F(n[1],n[2]),F(0),F(0),F(0),F(0)]
   for k in n[3]:assert jets[k] is not None;v=jmul(v,jets[k])
   jets[i]=v
  else:raise AssertionError(n)
 for j,col in d['old80_root_columns']:
  q=[[F(0),F(0)] for t in range(2)]
  for row,(kind,root) in col:
   if kind=='C0':
    if root==0:continue
    assert root in leaf;jet=leaf[root]
   else:assert kind=='typed' and root in jets and jets[root] is not None;jet=jets[root]
   for t in range(2):
    for k,mi in enumerate([3,4]):q[t][k]=add(q[t][k],mul(G[t].get(row,F(0)),jet[mi]))
  oldg[j]=q
 def weighted(cols,c):
  v={}
  for j,x in c.items():
   for i,y in cols[j].items():v[i]=add(v.get(i,F(0)),mul(x,y))
  return {i:x for i,x in v.items() if x}
 for n,c0 in enumerate(d['selected99_original768_constant_combinations']):
  co={j:F(x) for j,x in c0 if j<640};cn={j-640:F(x) for j,x in c0 if j>=640};sr=weighted(Sr,cn);se=weighted(Se,cn);srr=weighted(Srr,cn);sre=weighted(Sre,cn);ar=weighted(oldr,co);ae=weighted(olde,co);ks=[];ps=[]
  for t in range(2):
   cr={i:x for i,x in enumerate(C[t][0]) if x};ce={i:x for i,x in enumerate(C[t][1]) if x};kv=[];pv=[]
   for k in range(2):
    cursor={'phase':'396_actual_Kfull_cache_coefficient','selected99':n,'channel':t,'rr_or_rell':k};a=F(0)
    for j,x in co.items():a=add(a,mul(x,oldg[j][t][k]))
    if k==0:terms=[a,dot(cr,ar),dot(R[t],srr),dot(P[t],sr),F(0)]
    else:terms=[a,add(dot(cr,ae),dot(ce,ar)),dot(R[t],sre),dot(P[t],se),dot(Q[t],sr)]
    q=F(0)
    for x in terms:q=add(q,x)
    D,U=map(int,d['positive396_Kfull_envelopes'][n][t][k]);assert D>0 and D%q.denominator==0 and abs(q.numerator)*D<=U*q.denominator
    kv.append(str(q));pv.append(list(map(str,terms)))
   ks.append(kv);ps.append(pv)
  K.append(ks);pieces.append(ps)
 ep='actual396-Kfull-P5-and-native23-existing-cache-values.json'
 ev={'K99_two_channels_rr_rell':K,'all99_two_by2_five_piece_values':pieces,'five_piece_order':['OLD_Gamma_second','Cfull_OLDfirst_sum','R_intrinsic_second','P_intrinsic_first','Q_intrinsic_rfirst_or0'],'native_shared_R_two48':[rec(v) for v in R],'native_shared_P_two48':[rec(v) for v in P],'native_shared_Q_two48':[rec(v) for v in Q],'all396_source_positive_envelopes_checked':True,'fullPsi_Cfull_and_native_three_term_sum_retained':True,'honest327554_never_evaluated_or_zero_filled':True,'seven_h_products_zero_ONLY_by_source_order':True,'gBzero7_order3_not_rawf63':True,'scope':'ONLYcachedKfull396coefficients.Noall918adjoint/sourceZ2values/b1pairing/upper617/HF2orlength.'}
 owned=2*count(d)+count([G,C,L0,Lr,Le,Sr,Se,Srr,Sre,oldr,olde,leaf,jets,oldg,R,P,Q,K,pieces,ev,g,p,c,pins,cursor])+1000;assert owned<cap;dump(ep,ev)
 report={'id':STEM,'all_passed':True,'existing_cache_Kfull396_only':True,'all396_positive_bounds_checked':True,'source_committed_prior147':147,'temporary_prior442':442,'hard_guards':[384,1153],'actual_maximum_bits':maximum,'whole_numeric_cap':cap,'actual_owned_numeric_upper':owned,'RAW_output_cap':2000000,'no_b1_upper_HF_or_length_values':True,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};dump('report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','actual_maximum_bits','actual_owned_numeric_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output,a.source_peer,a.custody)

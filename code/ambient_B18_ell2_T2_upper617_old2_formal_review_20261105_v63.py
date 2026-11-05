#!/usr/bin/env python3
"""Own all617 formal maps from original source; no function responses."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63';STEM=A+'_REVIEW'
S='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_SOURCE_20261105_v63/'
G='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_INPUT_GATE_20261105_v63/'
P='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_INPUT_GATE_20261105_v63_REVIEW_V2/'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except (ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};maximum=0;positive_max=0;cursor={}
 def read(p):q=ROOT/p;pins[p]=sha(q);return json.loads(q.read_text())
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 r=read(f'reports/{A}/report.json');assert pins[f'reports/{A}/report.json']=='643f72906d5db1b74cbf100ed401160336d8c94619f3fd35f5b7d4e8872415e0' and r['all_passed']
 ev=read(f'reports/{A}/actual617-formal-combined-OLD2-second-symbol-coefficients.json');assert pins[f'reports/{A}/actual617-formal-combined-OLD2-second-symbol-coefficients.json']==r['evidence_sha256']['actual617-formal-combined-OLD2-second-symbol-coefficients.json']
 p=read(P+'report.json');assert p['all_passed'] and pins[P+'report.json']=='9fd284e851f7966fefe220fb83a61cf5991b72fa42b9f887e3962282f0f74187'
 g=read(G+'report.json');d=read(G+'617-combined-OLD2-independent-symbol-positive-input.json');assert pins[G+'617-combined-OLD2-independent-symbol-positive-input.json']==g['evidence_sha256']['617-combined-OLD2-independent-symbol-positive-input.json']
 sr=read(S+'report.json');s=read(S+'combined-OLD2-D-and-four-pivot-source-factor-inventory.json');assert pins[S+'combined-OLD2-D-and-four-pivot-source-factor-inventory.json']==sr['evidence_sha256']['combined-OLD2-D-and-four-pivot-source-factor-inventory.json']
 assert count([r,ev,p,g,d,sr,s,pins])+50000<150000
 # Independent positive prebound for the alternate grouping: precontract each
 # original source column with D before source/basis subtraction. Producer
 # contracts individual root weights first. No signed gradients yet.
 Z=(1,0)
 def posck(n):
  nonlocal positive_max
  positive_max=max(positive_max,abs(n).bit_length());assert abs(n).bit_length()<=120
  return n
 def ex(x):
  q=F(x);return (posck(q.denominator),posck(abs(q.numerator))) if q else Z
 def pa(a,b):
  if not a[1]:return b
  if not b[1]:return a
  den=posck(math.lcm(a[0],b[0]));return den,posck(posck(a[1]*(den//a[0]))+posck(b[1]*(den//b[0])))
 def pm(a,b):return (posck(a[0]*b[0]),posck(a[1]*b[1])) if a[1] and b[1] else Z
 pg={i:{f:(int(a),int(b)) for f,a,b in v} for i,v in d['positive_linear_node_gradients']}
 positive_cols={}
 for j,col in s['D_old_root_columns']:
  for t in range(2):
   z={}
   for row,node in col:
    for pair in range(10):
     if F(s['D_two_by10_by56'][t][pair][row]):
      for f,q in pg.get(node,{}).items():key=(pair,f);z[key]=pa(z.get(key,Z),pm(ex(s['D_two_by10_by56'][t][pair][row]),q))
   if z:positive_cols[j,t]=z
 assert positive_max<=120
 dump('positive-before-signed-checkpoint.json',{'alternate_original_column_positive_height':positive_max,
  'producer_original617_envelopes_pinned_by_9fd284':True,'signed_formal_or_response_values_not_started':True})
 def ck(q):
  nonlocal maximum
  maximum=max(maximum,abs(q.numerator).bit_length(),q.denominator.bit_length())
  if maximum>120:dump('guard-frontier.json',{'cursor':cursor,'offender':str(q),'maximum':maximum,'no_response_values':True});raise RuntimeError('independent120 STOP')
  return q
 def madd(dst,key,value):
  q=ck(dst.get(key,F(0))+value)
  if q:dst[key]=q
  else:dst.pop(key,None)
 def axpy(dst,src,scale):
  if not scale:return
  for k,v in src.items():madd(dst,k,ck(scale*v))
 # Full original757 node program, not the producer's pruned104 list.
 grad={}
 leaves={i for i,k,q in d['literal_function_origins71']}
 for i,n in s['typed_program']:
  cursor={'stage':'own_full757_point_differential','node':i}
  if n[0]=='actual_C0_scalar_function':assert n[1] in leaves;z={n[1]:F(1)}
  elif n[0]=='sum':
   z={}
   for j,a,b in n[1]:axpy(z,grad[j],F(a,b))
  else:
   assert n[0]=='product' and n[3];z={}
   if len(n[3])==1:axpy(z,grad[n[3][0]],F(n[1],n[2]))
  grad[i]=z
 oldcols={}
 for j,col in s['D_old_root_columns']:
  channels=[]
  for t in range(2):
   z={}
   for row,node in col:
    for pair in range(10):
     scale=F(s['D_two_by10_by56'][t][pair][row])
     if scale:
      for f,v in grad[node].items():madd(z,(pair,f),ck(scale*v))
   env=positive_cols.get((j,t),{});assert set(z)<=set(env)
   for key,q in z.items():den,u=env[key];assert den%q.denominator==0 and abs(q.numerator)*(den//q.denominator)<=u
   channels.append(z)
  if any(channels):oldcols[j]=channels
 # Literal negative OLD-minus-NEW edges are evaluated independently of typed
 # derivatives; all71 base function symbols stay independent.
 outer=s['outer_four_pivot_source'];cg={i:{i:F(1)} for i in leaves}
 for i,n in outer['literal_C0_program']:
  if n[0]=='g':assert i in cg
  elif n[0]=='n':cg[i]={f:ck(-v) for f,v in cg[n[1]].items()}
  else:assert n[0]=='a';cg[i]=dict(cg[n[1]]);axpy(cg[i],cg[n[2]],F(1))
 Gm=[dict(v) for v in outer['Gamma_two']];piv=dict(outer['four_source_root_columns']);old4=outer['original4_old_pivots']
 pivot_maps={}
 for ip,j in enumerate(old4):
  for var,pair in [(8,8),(7,9)]:
   for t in range(2):
    z={}
    for row,(kind,node) in piv[j]:
     if row not in Gm[t]:continue
     gg=grad[node] if kind=='typed' else cg[node]
     for f,v in gg.items():madd(z,(pair,f),ck(-F(Gm[t][row])*v))
    pivot_maps[ip,var,t]=z
 c619=[{j:F(x) for j,x in v if j<640} for v in s['original619_c0_source_combinations']]
 a619={j:{i:F(x) for i,x in v} for j,v in s['raw619_old4_plane_a1_nonempty28coefficients']}
 expected={j:[{(p,f):F(q) for p,f,q in v} for v in channels] for j,channels in ev['all617_combined_OLD2_independent_symbol_coefficients']}
 bound={(j,t):{(p,f):(int(a),int(b)) for p,f,a,b in v} for j,t,v in d['positive617_combined_symbol_coefficients']}
 owned0=count([r,ev,p,g,d,sr,s,pins,pg,positive_cols,grad,oldcols,cg,Gm,piv,pivot_maps,c619,a619,expected,bound])+10000
 assert owned0<150000
 records=[];nonempty=0;slots=0;symbols=set()
 for rec in s['kernel617_source_recipes']:
  j=rec['free_case619'];cursor={'stage':'independent_source_before_D_contraction','case':j};c=dict(c619[j]);a=dict(a619.get(j,{}))
  for b,v in zip([563,609],rec['subtract_actual_basis563609'],strict=True):axpy(c,c619[b],-F(v));axpy(a,a619.get(b,{}),-F(v))
  channels=[]
  for t in range(2):
   z={}
   for col,v in c.items():
    if col in oldcols:axpy(z,oldcols[col][t],v)
   for k,v in a.items():ip,var=divmod(k,9);assert ip<4 and var in [7,8];axpy(z,pivot_maps[ip,var,t],v)
   assert z==expected[j][t],(j,t,'full signed formal map mismatch')
   env=bound.get((j,t),{});assert set(z)<=set(env)
   for key,q in z.items():den,u=env[key];assert den%q.denominator==0 and abs(q.numerator)*(den//q.denominator)<=u
   channels.append([[p,f,str(v)] for (p,f),v in sorted(z.items())]);nonempty+=bool(z);slots+=len(z);symbols.update(f for pair,f in z)
  records.append([j,channels])
 assert len(records)==617 and nonempty==106 and slots==1909 and len(symbols)==42
 assert sorted(symbols)==r['actual_symbol_function_ids']==ev['actual_symbol_function_ids']
 fingerprint=hashlib.sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest()
 owned=count([r,ev,p,g,d,sr,s,pins,pg,positive_cols,grad,oldcols,cg,Gm,piv,pivot_maps,c619,a619,expected,bound,records])+12000
 assert owned<150000
 ep='independent-all617-original-column-grouped-formal-OLD2-identity.json'
 evidence={'all617_by2_maps_including_omitted_zeros_verified':True,'all1909_signed_symbol_slots_and106_nonemptychannels_verified':True,
  'full757_program_and_all71_independent_symbols_used':True,'outer_negative_edges_and_all28_a1_carries_retained':True,
  'grouping':'First originalOLDcolumn atom maps/D; then source563609 subtraction; separately signed original outer4 maps.',
  'all_own_formal_records_sha256':fingerprint,'surviving42_literal_function_ids':sorted(symbols),
  'positive_original_column_prebound_bits':positive_max,'actual_independent_formal_bits':maximum,
  'no_actual_second_firstcross_native_or_upper_response_values':True,'whole_owned_numeric_upper':owned}
 dump(ep,evidence);report={'id':STEM,'all_passed':True,'formal_identity_postpeer_only':True,'input_sha256':pins,
  'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
  'source_prior120_unchanged':True,'actual_formal_bits':maximum,'owned_numeric_upper':owned,'seconds':time.monotonic()-start}
 dump('report.json',report);dump('phase.json',{'phase':'FINAL_FORMAL_POSTPEER','report_sha256':sha(out/'report.json')});print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'owned':owned,'bits':maximum,'seconds':report['seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

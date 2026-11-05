#!/usr/bin/env python3
"""Prepared ONE92 Gamma-Xp2 worker; no execution without root/source peer.

Full Res4(N^t) and intrinsic S coefficients retain monic trace carry and every
individual Lrr*S0/Lrell*S0 term. No physical/gauge solve or honest evaluation.
"""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NATIVE_PIVOT23_SECOND_VALUES_20261105_v63'
GATE='reports/AMBIENT_B18_ELL2_NATIVE_PIVOT23_SECOND_SOURCE_20261105_v63_V2/'
GATE_SHA='ab2e8d901a32da8e2f0e410f86e5bf46512b74fcc3342ecc12edbb52ccb316a0'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except (ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def main(out,peer,custody):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 pins={};cursor={};N=[];G=[];nu=[];Gamma=[];A=[];S=[];targets=[];pieces=[];maximum=0
 def read(p):q=ROOT/p;pins[p]=sha(q);return json.loads(q.read_text())
 def dump(n,v):
  (out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
  assert sum(p.stat().st_size for p in out.iterdir())<1_000_000
 r=read(GATE+'report.json');assert r['all_passed'] and pins[GATE+'report.json']==GATE_SHA
 ip=GATE+'factored-native23-Gamma-pivot-plane2-source-and-full-L2S0.json';d=read(ip);assert pins[ip]==r['evidence_sha256'][Path(ip).name]
 p=read(peer);c=read(custody)
 assert p['all_passed'] and GATE_SHA in json.dumps(p) and c['all_passed'] and c['safe_help_literal_PASS']
 assert c['worker_sha256']==sha(Path(__file__)) and c['gate_sha256']==GATE_SHA
 for path,h in c['literal_input_sha256'].items():assert sha(ROOT/path)==h
 assert count([r,p,c,pins])<1000
 assert r['source_positive_committed_bits']==d['coarse_committed_bits']==280 and d['coarse_temporary_bits']==841
 assert r['prospective_value_numeric_bound']==d['prospective_value_whole_numeric_bound']==21034
 assert c['explicit_success_and_guard_numeric_upper']<21034 and c['whole_RAW_upper']<32*1024*1024
 cap=21034
 global_den=math.prod(int(v['D']) for v in d['source_positive_classes_D_U'].values())
 global_abs=math.prod(int(v['U']) for v in d['source_positive_classes_D_U'].values())*d['source_path_count_upper']
 assert max(global_den.bit_length(),global_abs.bit_length())==280
 def owned(extra=()):return 2*count(d)+count([r,p,c,pins,N,G,nu,Gamma,A,S,targets,pieces,extra])+1000
 dump('phase-before-values.json',{'phase':'literal_source_peer_custody_and_whole_owner_PASS','source280_temporary841':True,'cap21034':True,'fullL2S0_retained':True,'no_new_physical_gauge_or_honest_quotient_query':True})
 def stop(q,reason):
  state={'reason':reason,'cursor':cursor,'uncommitted_fraction':str(q),'Gamma_L_prefix':[[[[i,str(v)] for i,v in row.items()] for row in t] for t in A],
   'current_intrinsic_S':[[[i,str(v)] for i,v in row.items()] for row in S],
   'completed_target_prefix':targets,'completed_piece_prefix':pieces,'input_sha256':pins,'no_cap_raise_or_retry':True}
  assert owned(state)<cap;dump('guard-frontier.json',state);raise RuntimeError(reason)
 def ck(q):
  nonlocal maximum
  maximum=max(maximum,abs(q.numerator).bit_length(),q.denominator.bit_length())
  if maximum>280:stop(q,'source280 STOP before commit')
  return q
 def add(a,b):return ck(a+b)
 def mul(a,b):return ck(a*b)
 def madd(dst,key,v):
  q=add(dst.get(key,F(0)),v)
  if q:dst[key]=q
  else:dst.pop(key,None)
 def dot(a,b):
  q=F(0)
  for i,v in a.items():
   if i in b:q=add(q,mul(v,b[i]))
  return q
 N=[[[[ck(F(q)) for q in mask] for mask in a] for a in row] for row in d['cached_N_12_by12_plane5_free4_UNEVALUATED_IN_PIVOT_RESPONSE']]
 G=[[[ck(F(q)) for q in mask] for mask in a] for a in d['canonical_G8_plane5_free4_UNEVALUATED_IN_PIVOT_RESPONSE']]
 nu=[[ck(F(q)) for q in a] for a in d['canonical_nu_two_plane5']]
 # Full Gamma56 is source-pinned; native X's bottom8 rows are exactly zero.
 Gamma=[{i:ck(F(q)) for i,q in row if i<48 and F(q)} for row in d['Gamma_two']]
 mon=[(0,0),(1,0),(0,1),(2,0),(1,1)];idx={v:i for i,v in enumerate(mon)}
 assert all(v[:3]==[F(0)]*3 for v in nu)
 # Routing a scalar-C coefficient against a trace input basis. L2 includes all
 # N0*nu2 collisions. Never transpose mu4, and never set individual S0 to zero.
 def trace_routes(base,mask,inputmask):
  collision=mask&inputmask;outmask=mask^inputmask
  if not collision:return [(base,outmask,F(1))]
  if collision==3:return [] # nuA*nuB starts degree4.
  ans=[]
  for b,v in enumerate(nu[collision-1]):
   pair=(mon[base][0]+mon[b][0],mon[base][1]+mon[b][1])
   if v and pair in idx:ans.append((idx[pair],outmask,v))
  return ans
 A=[[{} for base in range(5)] for t in range(2)]
 for source in range(12):
  for inputmask in range(4):
   col=4*source+inputmask
   for physical in range(12):
    for base in range(5):
     for mask,v in enumerate(N[source][physical][base]):
      if not v:continue
      for b,om,carry in trace_routes(base,mask,inputmask):
       row=4*physical+om
       for t in range(2):
        if row not in Gamma[t]:continue
        cursor={'phase':'Gamma_times_full_Res4_Ntranspose_coefficient','source12':source,'inputmask':inputmask,'physical12':physical,'base':b,'channel':t}
        madd(A[t][b],col,mul(mul(Gamma[t][row],v),carry))
 assert owned()<cap
 for recipe in d['native23_pivot_recipes']:
  S=[{} for b in range(5)];im=recipe['input_trace_mask']
  for entry in recipe['S_source_entries']:
   if entry[0]=='G8':_,row,g,sign=entry;f=G[g]
   else:
    _,row,c0=entry;sign=1;f=[[F(c0),F(0),F(0),F(0)]]+[[F(0)]*4 for _ in range(4)]
   for base in range(5):
    for mask,v in enumerate(f[base]):
     if not v:continue
     for b,om,carry in trace_routes(base,mask,im):
      cursor={'phase':'full_intrinsic_Res4_S_coefficient','original768':recipe['original768_source'],'base':b,'source12':row}
      madd(S[b],4*row+om,mul(mul(F(sign),v),carry))
  tv=[];pv=[]
  for t in range(2):
   cursor={'phase':'individual_Gamma_Xp2_full_product','original768':recipe['original768_source'],'channel':t}
   rr=[dot(A[t][0],S[3]),dot(A[t][1],S[1]),dot(A[t][3],S[0])]
   re=[dot(A[t][0],S[4]),dot(A[t][1],S[2]),dot(A[t][2],S[1]),dot(A[t][4],S[0])]
   vr=F(0);ve=F(0)
   for v in rr:vr=add(vr,v)
   for v in re:ve=add(ve,v)
   for q in [vr,ve]:assert global_den%q.denominator==0 and abs(q.numerator)*global_den<=global_abs*q.denominator
   tv.append([str(vr),str(ve)]);pv.append([list(map(str,rr)),list(map(str,re))])
  targets.append(tv);pieces.append(pv);assert owned()<cap
 ep='actual92-Gamma-nativepivot23-plane2-and322-full-product-pieces.json'
 ev={'original768_native23_sources':[v['original768_source'] for v in d['native23_pivot_recipes']],
  'Gamma_Xp2_native23_twochannels_rr_rell':targets,'all23_two_by2_full_product_piece_values':pieces,
  'rr_piece_order':['GammaL0_Srr','GammaLr_Sr','GammaLrr_S0'],
  'rell_piece_order':['GammaL0_Srell','GammaLr_Sell','GammaLell_Sr','GammaLrell_S0'],
  'all_individual_L2S0_and_varying_monic_carries_retained':True,
  'all92_common_denominator_and_positive_source_envelopes_verified':True,
  'source_Gamma56_current27_and_original_Q48_physical_C0_frame_retained':True,
  'no_honest_quotient_new_source_gauge_N3_or_617_upper_HF_length_query':True}
 owner=owned(ev);assert owner<cap;dump(ep,ev)
 rep={'id':STEM,'all_passed':True,'only92_Gamma_nativeXp2_coefficients_and_full_product_pieces':True,
  'actual_max_bits':maximum,'source_committed280_temporary841':True,'hard_guards':[384,1153],
  'whole_numeric_cap':cap,'actual_owned_upper':owner,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 dump('report.json',rep);print(json.dumps({k:rep[k] for k in ['all_passed','actual_max_bits','actual_owned_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output,a.source_peer,a.custody)

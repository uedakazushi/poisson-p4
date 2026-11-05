#!/usr/bin/env python3
"""ONE554 weighted native/outer positive input, no signed response pairing."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NATIVE_ONE554_POSITIVE_INPUT_20261105_v63'
SP='reports/AMBIENT_B18_NATIVE_ONE554_FRAME_SOURCE_20261105_v63/'
BP='reports/AMBIENT_B18_NATIVE_ONE549_POSITIVE_INPUT_20261105_v63/'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0;bits=0;cursor={};S=[];gparts=[];dparts=[];outer=[]
 def read(p):
  nonlocal raw
  q=ROOT/p;pins[p]=sha(q);raw+=2*q.stat().st_size;return json.loads(q.read_text())
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 sr=read(SP+'report.json');assert sr['all_passed'] and pins[SP+'report.json']=='cb03c3d5cfbc9cba86f9c15280b303a032dafa08fcddaad9ae782168171b13f5'
 sn='ONE554-four-native-weighted-S0-currentframe-nonzero644-outer-source.json';s=read(SP+sn);assert pins[SP+sn]==sr['evidence_sha256'][sn]
 br=read(BP+'report.json');assert br['all_passed'] and pins[BP+'report.json']=='62d66bc4406244a863719b5e0ec9cdb2b780c5b2d61475c39ecc6796512f1b4f'
 bn='ONE549-positive-native649-GammaD-input-and-known648zero-coupling.json';b=read(BP+bn);assert pins[BP+bn]==br['evidence_sha256'][bn]
 source=b['native_primitive_input_from_frozen_source'];assert source['Gamma_two56']==s['current_Gamma_two56'] and source['D_two_by10_by56']==s['current_D_two_by10_by56']
 assert all(v['full_point_S0_column48']==[] for v in s['parameterized_four_native_source_records'])
 assert s['outer_pivot644_mask0_FULL_S0_column_nonzero']==[[11,'-1/2']]
 assert any(q!='0' for row in s['outer_full_GammaXp644_rr_rell_actual_cached_twochannels'] for q in row)
 initial=count([sr,s,br,b,source,pins])+40000;assert initial<150000 and raw+4_000_000<32*1024*1024
 dump('phase-before-positive.json',{'phase':'same_J18_N_G_nu_GammaD_four_masked_source_columns_and_true644outer_PASS','source_owner_upper':initial,'no_signed_response':True})
 Z=(1,0)
 def ck(n):
  nonlocal bits
  bits=max(bits,abs(n).bit_length())
  if abs(n).bit_length()>384:
   dump('positive-guard-frontier.json',{'cursor':cursor,'uncommitted_integer':str(n),'S_positive_prefix':S,'Gamma_part_prefix':gparts,'D_part_prefix':dparts,'outer_prefix':outer,'inputs':pins,'no_signed_values_or_capraise_retry':True});raise RuntimeError('SOURCE384 STOP')
  return n
 def ex(x):q=F(x);return (ck(q.denominator),ck(abs(q.numerator))) if q else Z
 def add(a,b):
  if not a[1]:return b
  if not b[1]:return a
  d=ck(math.lcm(a[0],b[0]));return d,ck(ck(a[1]*(d//a[0]))+ck(b[1]*(d//b[0])))
 def mul(a,b):return (ck(a[0]*b[0]),ck(a[1]*b[1])) if a[1] and b[1] else Z
 def madd(v,k,q):
  if q[1]:v[k]=add(v.get(k,Z),q)
 basis=[tuple(v) for v in source['basis18']];bi={m:i for i,m in enumerate(basis)}
 G=[[ex(q) for q in v] for v in source['G8_C72']];nu=[[ex(q) for q in v] for v in source['nu_two_J18']]
 S=[{} for _ in basis]
 for rec in s['parameterized_four_native_source_records']:
  weight=ex(rec['constant_source_coefficient']);im=rec['input_trace_mask']
  for kind,row,g,sign in rec['intrinsic_source_entries']:
   assert kind=='G8' and sign in [-1,1]
   for base in range(18):
    for mask in range(4):
     v=G[g][4*base+mask]
     if not v[1]:continue
     collision=mask&im;om=mask^im
     if not collision:madd(S[base],4*row+om,mul(weight,v))
     elif collision in [1,2]:
      for nb,q in enumerate(nu[collision-1]):
       m=tuple(a+b for a,b in zip(basis[base],basis[nb]))
       if q[1] and m in bi:madd(S[bi[m]],4*row+om,mul(mul(weight,v),q))
     else:assert collision==3 # inherited literal νAνB supportzero in J18.
 assert S[0]=={} # independently source-masked point columns, not a postulated zero.
 rows={tuple(k):{i:(int(d),int(u)) for i,d,u in v} for k,v in b['positive_covectorL_row_cache']}
 # These68 positive operator rows have identical Γ/D/N/ν source inputs here.
 # They are reused, not recomputed. A fresh554 audit must check numeric bound
 # applicability beyond the Boolean zero-support proof used for case549.
 def R(kind,t,p,k):return rows[kind,t,p,k]
 def dot(a,c):
  z=Z
  for i,v in a.items():
   if i in c:z=add(z,mul(v,c[i]))
  return z
 vf={0:1,1:2,4:3,5:4,7:5,8:6};cubic=[(0,17),(5,16),(6,15),(15,6),(16,5)];pairs=s['D_ordinary_pair_order10']
 total=[]
 for t in range(2):
  gg=[dot(R(0,t,-1,k),S[j]) for k,j in cubic];gparts.append(gg);dd=[]
  for p,(a,c) in enumerate(pairs):
   vals=[dot(R(1,t,p,0),S[7+p])]
   for x,y in ([(a,c)] if a==c else [(a,c),(c,a)]):vals.append(dot(R(1,t,p,vf[x]),S[vf[y]]))
   dd.append(vals)
  dparts.append(dd)
  # Keep the genuinely nonzero644 prior. Only ABS products are formed here;
  # signed outer pairing belongs to the future ONE554 whole value job.
  qrr=ex('-377/60');qre=ex('-175/3');rr,re=s['outer_full_GammaXp644_rr_rell_actual_cached_twochannels'][t]
  oo=[mul(qrr,ex(rr)),mul(qre,ex(re))];outer.append(oo);z=Z
  for q in gg:z=add(z,q)
  for v in dd:
   for q in v:z=add(z,q)
  for q in oo:z=add(z,q)
  total.append(z)
 assert bits<=384 and 3*bits+1<=1153
 enc=lambda q:[str(q[0]),str(q[1])]
 data={'selected_case619':554,'six_source_and_all29_a1_native_frame_report':SP+'report.json',
  'full_primitive_input_J18_N_G_nu_GammaD':source,'parameterized_four_native_source_records':s['parameterized_four_native_source_records'],
  'all29_a1_kept_UNEVALUATED':s['all29_a1_kept_UNEVALUATED'],
  'positive_weighted_intrinsic_S18':[[[i,*enc(q)] for i,q in sorted(v.items())] for v in S],
  'reused68_positive_covectorL_rows':b['positive_covectorL_row_cache'],
  'positive_Gamma_weighted_cubic_two_by5':[[enc(q) for q in v] for v in gparts],
  'positive_D_weighted_second_two_by10_unique_first_decomp':[[[enc(q) for q in v] for v in ch] for ch in dparts],
  'actual644_GammaXp_rr_rell_twochannels_UNPAIRED':s['outer_full_GammaXp644_rr_rell_actual_cached_twochannels'],
  'plane_a1_pivot5_original644':s['only_plane_a1_pivot5_original644'],
  'positive_outer644_two_by2':[[enc(q) for q in v] for v in outer],
  'positive_native_plus_absouter_two':list(map(enc,total)),
  'signed_future_formula':'nativeΓcubic+nativeDsecond SUBTRACT[ΓXp644_rr*(-377/60)+ΓXp644_re*(-175/3)]; addOLDpartandSUBTRACTallFIRST57 inwholebundle.',
  'full_varying_nu_unknownN2N3_kernel_fullL2S0_distinctions_preserved':True,
  'old648_zero_template_NOT_used':True,'same68_operator_inputs_checked_row_numeric_prior_fresh_audit_required':True,
  'external_OLD_FIRST57_lane':s['external_OLD_and57_FIRST_lane'],
  'positive_committed_bits':bits,'positive_temporary_bits':3*bits+1,
  'whole_joint_VALUE_gate_and_worker_custody_STILL_REQUIRED':True,'no_signed_native_outer_or_upper_pairing':True}
 n=count(data);future=2*(n+1000)+count([source['N_known12_by12_by9_by4'],source['G8_C72'],source['nu_two_J18']])+count([rows,S])+8000
 data['native_plus_outer_only_future_numeric_bound']=future
 owner=count([sr,s,br,b,source,pins,G,nu,S,rows,gparts,dparts,outer,total,data])+20000
 assert owner<150000 and future<150000 and raw+4_000_000<32*1024*1024
 ep='ONE554-positive-four-native-and-nonzero644-outer-input.json';dump(ep,data)
 rep={'id':STEM,'all_passed':True,'positive_only_no_signed_native_outer_or_upper':True,'prior_bits':bits,'temporary_bits':3*bits+1,
  'source_owned_numeric_upper':owner,'native_outer_only_future_numeric_bound':future,'whole_VALUE_gate_closed':False,
  'fresh_numeric68_row_envelope_audit_required_for_nonzero554_scope':True,'whole_SOURCE_RAW_upper':raw+4_000_000,
  'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 dump('report.json',rep);print(json.dumps({k:rep[k] for k in ['all_passed','prior_bits','temporary_bits','source_owned_numeric_upper','native_outer_only_future_numeric_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

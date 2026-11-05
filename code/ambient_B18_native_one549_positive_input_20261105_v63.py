#!/usr/bin/env python3
"""Native ONE549 positive-only coefficient templates, not signed responses."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NATIVE_ONE549_POSITIVE_INPUT_20261105_v63'
SP='reports/AMBIENT_B18_NATIVE_ONE549_FRAME_SOURCE_20261105_v63/'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple,set)):return sum(map(count,x))
 if isinstance(x,str):
  try:F(x);return 1
  except(ValueError,ZeroDivisionError):return 0
 return int(isinstance(x,(int,float,F)))
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0;bits=0;cursor={};S=[];rows={};gparts=[];dparts=[];totals=[]
 def read(p):
  nonlocal raw
  q=ROOT/p;pins[p]=sha(q);raw+=2*q.stat().st_size;return json.loads(q.read_text())
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 r=read(SP+'report.json');assert r['all_passed'] and pins[SP+'report.json']=='6fe41836acc116ca03338d59ca91515351ce61f61fd9d66b27bb68538892f66b'
 name='ONE549-J18-with-r2ell-native-full-S0-column-and-formal-N2N3-frontier.json';s=read(SP+name);assert pins[SP+name]==r['evidence_sha256'][name]
 vp='reports/AMBIENT_B18_NATIVE_PIVOT23_SECOND_VALUES_20261105_v63/actual92-Gamma-nativepivot23-plane2-and322-full-product-pieces.json';v92=read(vp)
 dp='reports/AMBIENT_B18_NATIVE92_ROOT_DIRECT_REVIEW_20261105_v63/report.json';direct=read(dp);assert direct['all_passed'] and pins[dp]=='10844d8720465054d6b9f7a7aa408dd95dd613a2e857eb52e50cd5a71a847e92'
 pos=v92['original768_native23_sources'].index(648);assert pos==2 and v92['Gamma_Xp2_native23_twochannels_rr_rell'][pos]==[['0','0'],['0','0']]
 external=s['external_OLD298_quarter_and_full32_knownFIRSTcarry_lane'];oldr=read(external['source_report']);assert oldr['all_passed'] and pins[external['source_report']]==external['source_sha256']
 assert s['full_intrinsic_Res4_S0_column9_exactempty']==[] and s['comparison_pivot648_mask0_S0_column8_nonzero']==[[15,'-1/2']]
 initial=count([r,s,v92,direct,oldr,pins]);assert initial+70000<150000 and raw+4_000_000<32*1024*1024
 dump('phase-before-positive.json',{'phase':'literal_source_known92pivotzero_and_external_lane_PASS','initial_numeric':initial,'source_reserve70000':True,'no_signed_native_responses':True})
 Z=(1,0)
 def ck(n):
  nonlocal bits
  bits=max(bits,abs(n).bit_length())
  if abs(n).bit_length()>384:
   dump('positive-guard-frontier.json',{'cursor':cursor,'uncommitted_integer':str(n),'intrinsic_S_positive_prefix':S,'row_cache_positive_prefix':[[list(k),[[i,str(q[0]),str(q[1])] for i,q in x.items()]] for k,x in rows.items()],
    'Gamma_part_prefix':gparts,'D_part_prefix':dparts,'inputs':pins,'no_actual_native_response_or_retry_capraise':True})
   raise RuntimeError('SOURCE positive384 STOP')
  return n
 def ex(q):
  q=F(q);return (ck(q.denominator),ck(abs(q.numerator))) if q else Z
 def add(a,b):
  if not a[1]:return b
  if not b[1]:return a
  d=ck(math.lcm(a[0],b[0]));return d,ck(ck(a[1]*(d//a[0]))+ck(b[1]*(d//b[0])))
 def mul(a,b):return (ck(a[0]*b[0]),ck(a[1]*b[1])) if a[1] and b[1] else Z
 def madd(x,k,q):
  if q[1]:x[k]=add(x.get(k,Z),q)
 basis=[tuple(v) for v in s['ordinary_basis_J18']];bi={m:i for i,m in enumerate(basis)};known=s['N_known_J18_positions'];Nk=s['N_known12_by12_by9_by4_cached_operands_ONLY']
 G=[[ex(q) for q in v] for v in s['already_accepted_G8_J18_C72_images_COPIED_NOT_REEVALUATED']]
 nu=[[ex(q) for q in v] for v in s['already_accepted_nu_two_J18_images']]
 N=[[[[ex(q) for q in v] for v in a] for a in row] for row in Nk]
 # Intrinsic source2/inputmask1: G3 inrow0, -G0 inrow3. Absolute bounds retain
 # both rows; nuA collision is inserted from the SAME canonical monics.
 S=[{} for _ in range(18)]
 for row,g in [(0,3),(3,0)]:
  for base in range(18):
   for mask in range(4):
    x=G[g][4*base+mask]
    if not x[1]:continue
    if not(mask&1):madd(S[base],4*row+(mask^1),x)
    else:
     for nb,y in enumerate(nu[0]):
      m=tuple(a+b for a,b in zip(basis[base],basis[nb]))
      if y[1] and m in bi:madd(S[bi[m]],4*row+(mask^1),mul(x,y))
 assert S[0]=={} # FULL Res4 input column9, not entrywise G8=0.
 Gamma=[{i:ex(q) for i,q in v if i<48} for v in s['Gamma_two_current56']]
 D=[[{i:ex(q) for i,q in enumerate(v) if i<48 and q!='0'} for v in ch] for ch in s['D_two_by10_by56']]
 def R(kind,t,pair,target):
  key=(kind,t,pair,target)
  if key in rows:return rows[key]
  cov=Gamma[t] if kind==0 else D[t][pair];ans={};cursor.update({'phase':'positive_covector_L_known_coefficient','rowkey':key})
  for source in range(12):
   for im in range(4):
    col=4*source+im
    for physical in range(12):
     for kpos,nb in enumerate(known):
      for mask,x in enumerate(N[source][physical][kpos]):
       om=mask^im;outrow=4*physical+om
       if not x[1] or outrow not in cov:continue
       collision=mask&im
       if not collision:
        if nb==target:madd(ans,col,mul(cov[outrow],x))
       elif collision in [1,2]:
        for vb,y in enumerate(nu[collision-1]):
         if y[1] and tuple(a+b for a,b in zip(basis[nb],basis[vb]))==basis[target]:madd(ans,col,mul(mul(cov[outrow],x),y))
       else:assert collision==3 # nuA*nuB=0 support certificate in SOURCE6fe.
  rows[key]=ans;return ans
 def dot(a,b):
  q=Z
  for i,x in a.items():
   if i in b:q=add(q,mul(x,b[i]))
  return q
 variable_first={0:1,1:2,4:3,5:4,7:5,8:6}
 cubic=[(0,17),(5,16),(6,15),(15,6),(16,5)]
 pairs=[tuple(v) for v in s['D_ordinary_pair_order10']]
 for t in range(2):
  cursor={'phase':'positive_native649_Gamma_cubic','channel':t};gg=[dot(R(0,t,-1,k),S[j]) for k,j in cubic];gparts.append(gg)
  dd=[]
  for p,(a,b) in enumerate(pairs):
   cursor={'phase':'positive_native649_D_second','channel':t,'pair':p};q=[dot(R(1,t,p,0),S[7+p])]
   # Ordinary convolution has rr exactly ONCE; distinct pairs have2 terms.
   decomp=[(a,b)] if a==b else [(a,b),(b,a)]
   for x,y in decomp:q.append(dot(R(1,t,p,variable_first[x]),S[variable_first[y]]))
   dd.append(q)
  dparts.append(dd);total=Z
  for q in gg:total=add(total,q)
  for qq in dd:
   for q in qq:total=add(total,q)
  totals.append(total)
 assert bits<=384 and 3*bits+1<=1153
 enc=lambda q:[str(q[0]),str(q[1])]
 payload={'selected_case619':549,'native_original768_source649_coefficient':'1','J18_with_r2ell_source_report':SP+'report.json',
  'native_primitive_input_from_frozen_source':{'basis18':s['ordinary_basis_J18'],'N_known_positions':known,'N_known12_by12_by9_by4':Nk,
   'G8_C72':s['already_accepted_G8_J18_C72_images_COPIED_NOT_REEVALUATED'],'nu_two_J18':s['already_accepted_nu_two_J18_images'],
   'Gamma_two56':s['Gamma_two_current56'],'D_two_by10_by56':s['D_two_by10_by56']},
  'all32_a1_source_entries_kept_UNEVALUATED':s['full_a1_27times9_32_entries_UNEVALUATED'],
  'positive_intrinsic_source9_columns18':[[[i,*enc(q)] for i,q in sorted(v.items())] for v in S],
  'positive_covectorL_row_cache':[[list(k),[[i,*enc(q)] for i,q in sorted(v.items())]] for k,v in sorted(rows.items())],
  'native_Gamma_cubic5_term_order':cubic,'native_Gamma_cubic_two_by5_D_U':[[enc(q) for q in v] for v in gparts],
  'native_D_second10_unique_first_decomposition_D_U':[[[enc(q) for q in v] for v in ch] for ch in dparts],
  'positive_native_Gamma_plus_D_two_D_U':list(map(enc,totals)),
  'outer_GammaXp648_plane_a1_exact0_by_independent92cache':{'data':vp,'peer':dp,'four_Gamma_rr_rell_targets':v92['Gamma_Xp2_native23_twochannels_rr_rell'][pos]},
  'OLD298quarter_and_all55_FIRST_carries_external_lane':external,
  'external_positive_payload_and_joint_execution_custody_STILL_REQUIRED':True,
  'unknown_N2mixed_and_N3_never_set_to_actual0':True,'native_S0_column9_zero_is_the_only_whole_kernel_cancellation':True,
  'full_varying_nu_Res4_first_and_planeL2_carries_retained':True,'ordinary_rr_first_product_occurs_once_not_twice':True,
  'no_signed_native_functional_or_whole_upper_response_computed':True,
  'positive_committed_bits':bits,'positive_temporary_bits':3*bits+1}
 # Native-only future count, with raw+parsed input copies, scalarN/G/nu,
 # allpotentialactualrow/support copies and complete guard-prefix duplicate.
 pcount=count(payload);operators=count(payload['positive_covectorL_row_cache'])+count(payload['positive_intrinsic_source9_columns18'])
 future=2*(pcount+1000)+count([Nk,s['already_accepted_G8_J18_C72_images_COPIED_NOT_REEVALUATED'],s['already_accepted_nu_two_J18_images']])+operators+7000
 payload['prospective_native_only_numeric_bound']=future
 payload['NO_claim_this_is_complete_joint_whole_VALUE_bound']=True
 owner=count([r,s,v92,direct,oldr,pins,N,G,nu,S,Gamma,D,rows,gparts,dparts,totals,payload])+20000
 assert owner<150000 and future<150000 and raw+4_000_000<32*1024*1024
 ep='ONE549-positive-native649-GammaD-input-and-known648zero-coupling.json';dump(ep,payload)
 rep={'id':STEM,'all_passed':True,'positive_native_templates_only_no_signed_response':True,'positive_committed_bits':bits,'positive_temporary_bits':3*bits+1,
  'source_owned_numeric_upper':owner,'native_only_future_numeric_bound':future,'whole_joint_VALUE_gate_closed':False,
  'external_OLD_and55_FIRST_positive_payload_needed_before_any_values':True,'whole_SOURCE_RAW_upper':raw+4_000_000,
  'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 dump('report.json',rep);print(json.dumps({k:rep[k] for k in ['all_passed','positive_committed_bits','positive_temporary_bits','source_owned_numeric_upper','native_only_future_numeric_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

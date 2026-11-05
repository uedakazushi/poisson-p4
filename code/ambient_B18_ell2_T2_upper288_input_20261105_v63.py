#!/usr/bin/env python3
"""Positive-only ONE288 cache input: three physical functions and FIRST carry."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER288_INPUT_GATE_20261105_v63';Z=(1,0)
def main(out):
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(exist_ok=False,parents=True);pins={};raw=0;bits=0;cursor={}
 def read(p):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024
  return json.loads((ROOT/p).read_text())
 def field(p,k):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024
  return H.selected(p,k)
 def ck(n):
  nonlocal bits
  bits=max(bits,abs(n).bit_length())
  if abs(n).bit_length()>384:
   H.dump(out/'BEFORE_RESPONSE_POSITIVE_GUARD.json',{'cursor':cursor,'uncommitted_positive_integer':str(n),'no_signed_response_values':True});raise RuntimeError('positive384STOP')
  return n
 def ex(q):q=F(q);return (ck(q.denominator),ck(abs(q.numerator))) if q else Z
 def add(a,b):
  if not a[1]:return b
  if not b[1]:return a
  den=ck(math.lcm(a[0],b[0]));u=ck(a[1]*(den//a[0]));v=ck(b[1]*(den//b[0]));return den,ck(u+v)
 def mul(a,b):return (ck(a[0]*b[0]),ck(a[1]*b[1])) if a[1] and b[1] else Z
 def dot(a,v):
  z=Z
  for i,q in a.items():z=add(z,mul(ex(q),v.get(i,Z)))
  return z
 def mv(cols,v):
  z={}
  for j,q in v.items():
   for i,x in cols[j].items():z[i]=add(z.get(i,Z),mul(ex(x),ex(q)))
  return z
 def va(a,b):
  z=dict(a)
  for i,q in b.items():z[i]=add(z.get(i,Z),q)
  return z
 sp='reports/AMBIENT_B18_ELL2_T2_UPPER288_SOURCE_20261105_v63/'
 sr=read(sp+'report.json');assert sr['all_passed'] and pins[sp+'report.json']=='91fff8d1815310c4b5e197382b99667b6c8e303e62e5e7bac5b69c7a3f35ea0c'
 s=read(sp+'ONE288-opaque-free-source-and-eight-known-first-pivot-carries.json');assert s['function_leaf_ids']==[20117,27948,40895] and s['actual_a1_r_and_ell_both_exact0'] and s['native_c0_exact_empty']
 assert len(s['full_typed_program'])==5 and s['literal_C0_program']==[[0,['c',0,1]]]
 rawp='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json';allphys=field(rawp,'all_requested164_root_values');ph={i:dict(v) for i,v in allphys if i in s['function_leaf_ids']};del allphys
 fpath='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json';f=field(fpath,'first640');fst={i:dict(v) for i,v in f['all79_atom_values_full10'] if i in ph};del f
 phys=[]
 for i in s['function_leaf_ids']:
  v=[ph[i].get(k,'0') for k in s['raw33_to_B21_selector21']];assert list(map(F,v[:10]))==[F(fst[i].get(k,'0')) for k in range(10)];phys.append([i,v])
 D=s['D_two_by10_by56'];Dph=[[[v[4],v[6]] for v in ch] for ch in D];Dcarry=[[ch[2],ch[3]] for ch in D];g4=[dict(v)[4] for v in s['Gamma_two']]
 L0={j:dict(v) for j,v in enumerate(s['native_L0'])};Le={j:dict(v) for j,v in enumerate(s['native_Le'])};S0=dict(enumerate(s['native_S0']));Se=dict(enumerate(s['native_Se']))
 piv=[];rbound=[];ebound=[];useL0=set();useLe=set()
 for rec in s['known_first_carry_pivots8']:
  cursor={'stage':'positive_actual_individual_first_pivot','pivot27':rec['pivot27']}
  r={i:ex(q) for i,q in rec['Xp_r_literal56']};kind,arg=rec['Xp_ell_recipe'];v=dict(rec)
  if kind=='literal':e={i:ex(q) for i,q in arg}
  else:
   assert kind=='native';p0=dict(S0[arg]);pe=dict(Se[arg]);useL0.update(pe);useLe.update(p0)
   # Actual SUM; neither source point S0 nor the Le*S0 term is discarded.
   e=va(mv(L0,pe),mv(Le,p0));v['native_Se_column']=Se[arg];v['native_S0_column']=S0[arg]
  piv.append(v);rbound.append(r);ebound.append(e)
 source_pre=count([s,ph,fst,phys,Dph,Dcarry,L0,Le,piv,rbound,ebound,pins])+20000;assert source_pre<60000
 pe={i:[ex(q) for q in v] for i,v in phys};delta=[add(a,b) for a,b in zip(pe[20117],pe[40895],strict=True)];z1=pe[27948]
 direct=[];second=[];carry=[];upper=[]
 for t in range(2):
  a=mul(ex(g4[t]),delta[20]);b=Z;v=Z
  for p in range(10):b=add(b,mul(ex(Dph[t][p][0]),delta[10+p]));b=add(b,mul(ex(Dph[t][p][1]),z1[10+p]))
  for j,rec in enumerate(piv):
   cursor={'stage':'positive_known_FIRST_carry_SUM','channel':t,'pivot27':rec['pivot27']}
   r=dot({i:q for i,q in enumerate(Dcarry[t][0]) if q!='0'},rbound[j]);e=dot({i:q for i,q in enumerate(Dcarry[t][1]) if q!='0'},ebound[j]);v=add(v,mul(ex(rec['a1_z1']),add(r,e)))
  direct.append(a);second.append(b);carry.append(v);upper.append(add(add(a,b),v))
 qpath='reports/AMBIENT_B18_ELL2_T2_UPPER617_FULL_COORDINATE_SOURCE_20261105_v63/known-Psi-point-correction-and-full-b1-source-annihilation.json';Q=field(qpath,'Q_two6');six=[]
 for k in range(6):
  v=Z
  for t in range(2):v=add(v,mul(ex(dict(Q[t]).get(k,'0')),upper[t]))
  six.append(v)
 enc=lambda q:[str(q[0]),str(q[1])]
 pv=lambda v:[[i,*enc(q)] for i,q in sorted(v.items())]
 data={'selected_case619':288,'physical_B21_functions3':phys,'Gamma_row4_two':g4,'D_physical_rows4_6_two_by10_by2':Dph,'D_z1r_and_z1ell_two_by2_by56':Dcarry,'known_FIRST_carry_pivots8':piv,'native_L0_selected_columns':[[j,s['native_L0'][j]] for j in sorted(useL0)],'native_Le_selected_columns':[[j,s['native_Le'][j]] for j in sorted(useLe)],'fixed_Q_embedding_two6':Q,'full_literal_typed_program5':s['full_typed_program'],'literal_C0_zero':s['literal_C0_program'],'positive_Xp_r_eight':list(map(pv,rbound)),'positive_Xp_ell_eight':list(map(pv,ebound)),'positive_Gamma_cubic_two':list(map(enc,direct)),'positive_D_physical_second_two':list(map(enc,second)),'positive_known_FIRST_carry_two':list(map(enc,carry)),'positive_full_upper_two':list(map(enc,upper)),'positive_full_six':list(map(enc,six)),'source_formula':'u_t=Gamma4*(phi20117-phi40895)[20]+sum_pD[t,p,4]*(phi20117-phi40895)[10+p]+D[t,p,6]*phi27948[10+p] -sum_a1z1*(D[t,z1r]Xp_r+D[t,z1ell]Xp_ell).','native_Xpell_is_SUM_L0Se_plus_LeS0_not_union':True,'a1_only_z1_so_pure_r2ell_GammaXp2a1_exact0':True,'native_c0_empty_so_no_nativeX3':True,'all_firsts_true_and_physical_higher_actual_raw33_ring_images':True,'no_honest_second_third_or_conductor_order_assumption':True,'no_signed_upper_response_evaluated':True}
 n=count(data);scratch=count([rbound,ebound,L0,Le,ph,delta,z1])+4000;cap=2*(n+1000)+scratch+3000;rawcap=2_000_000+cap*(2*math.ceil(bits*math.log10(2))+50)
 assert cap<60000 and rawcap<32*1024*1024 and bits<=384 and 3*bits+1<=1153
 data.update({'source_committed_bits':bits,'source_temporary_bits':3*bits+1,'input_numeric_before_finalmeta':n,'future_scratch_upper':scratch,'future_numeric_cap':cap,'future_RAW_cap':rawcap,'metadata_reserve':1000})
 ep=out/'ONE288-three-physical-functions-first-carry-positive-input.json';H.dump(ep,data)
 for p in ['reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_UPPER617_D_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json']:
  assert read(p)['all_passed']
 owned=count([s,phys,ph,fst,L0,Le,piv,rbound,ebound,data,Q,pins])+10000;assert owned<60000
 report={'id':STEM,'all_passed':True,'positive_known_cache_only':True,'no_new_function_or_actual_upper_response':True,'only_three_literal164_functions_and_eight_FIRST_carries':True,'all_individual_Xpell_SUMs_preserved':True,'no_native_second_third_or_honest_higher_needed':True,'source_prearith_bound':source_pre,'source_owned_bound':owned,'committed_prior_bits':bits,'temporary_prior_bits':3*bits+1,'future_numeric_cap':cap,'future_RAW_cap':rawcap,'whole_SOURCE_RAW_upper':raw,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','committed_prior_bits','temporary_prior_bits','source_owned_bound','future_numeric_cap','future_RAW_cap','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

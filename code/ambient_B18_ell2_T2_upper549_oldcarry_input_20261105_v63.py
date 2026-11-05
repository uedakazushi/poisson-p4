#!/usr/bin/env python3
"""Positive OLDquarter and all55 known-FIRST carries, no ONE549 responses."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER549_OLD_CARRY_INPUT_GATE_20261105_v63';Z=(1,0)
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
   H.dump(out/'BEFORE_RESPONSE_POSITIVE_GUARD.json',{'cursor':cursor,'uncommitted_positive_integer':str(n),'no_signed_response':True});raise RuntimeError('positive384STOP')
  return n
 def ex(q):q=F(q);return (ck(q.denominator),ck(abs(q.numerator))) if q else Z
 def add(a,b):
  if not a[1]:return b
  if not b[1]:return a
  den=ck(math.lcm(a[0],b[0]));u=ck(a[1]*(den//a[0]));v=ck(b[1]*(den//b[0]));return den,ck(u+v)
 def mul(a,b):return (ck(a[0]*b[0]),ck(a[1]*b[1])) if a[1] and b[1] else Z
 def dot(row,v):
  z=Z
  for i,q in enumerate(row):
   if q!='0':z=add(z,mul(ex(q),v.get(i,Z)))
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
 sp='reports/AMBIENT_B18_ELL2_T2_UPPER549_OLD_CARRY_SOURCE_20261105_v63/'
 sr=read(sp+'report.json');assert sr['all_passed'] and pins[sp+'report.json']=='6c5f77ba80584c564f22a99e6dc4e5d55f262fbd591a50fc06e8e13a5d642db2'
 s=read(sp+'ONE549-OLD298-and-full32-a1-known-FIRST-carry-factor-source.json')
 assert len(s['all_known_FIRSTcarry_factor_recipes'])==55 and len(s['all_needed_Xp_first_column_recipes'])==22
 physall=field('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json','all_requested164_root_values');ph={i:dict(v) for i,v in physall if i in s['OLD_physical_function_leaf_ids3']};del physall
 first=field('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json','first640');fst={i:dict(v) for i,v in first['all79_atom_values_full10'] if i in ph};del first
 phys=[]
 for i in sorted(ph):
  v=[ph[i].get(k,'0') for k in s['raw33_to_B21_selector21']];assert list(map(F,v[:10]))==[F(fst[i].get(k,'0')) for k in range(10)];phys.append([i,v])
 assert [i for i,v in phys]==[7454,20132,40910]
 L0={j:dict(v) for j,v in enumerate(s['native_L0'])};Le={j:dict(v) for j,v in enumerate(s['native_Le'])};use0=set();usee=set();env={};recipes=s['all_needed_Xp_first_column_recipes']
 for rec in recipes:
  p,var=rec['pivot27'],rec['variable9'];r=rec['first_column_recipe'];cursor={'stage':'positive_individual_cached_Xp_FIRST','pivot27':p,'variable9':var}
  if r['kind']=='literal56':v={i:ex(q) for i,q in r['column']}
  else:
   assert r['kind']=='native_ell_SUM' and var==8;a=dict(r['Se']);b=dict(r['S0']);use0.update(a);usee.update(b)
   v=va(mv(L0,a),mv(Le,b))
  env[p,var]=v
 pe={i:[ex(q) for q in v] for i,v in phys};delta=[add(a,b) for a,b in zip(pe[20132],pe[40910],strict=True)];coord=pe[7454]
 D=s['D_two_by10_by56'];G=[dict(v) for v in s['Gamma_two']];old3=[];old2=[];carry=[];factor_bounds=[];partial=[]
 for t in range(2):
  a=mul(ex('1/4'),mul(ex(G[t][1]),coord[20]));b=Z;v=Z;ff=[]
  for p in range(10):
   b=add(b,mul(ex(D[t][p][0]),delta[10+p]));b=add(b,mul(ex(D[t][p][1]),coord[10+p]))
  b=mul(ex('1/4'),b)
  for k,(p,pair,var,q) in enumerate(s['all_known_FIRSTcarry_factor_recipes']):
   cursor={'stage':'positive_all55_FIRST_carry_factors','channel':t,'factor':k}
   z=mul(ex(q),dot(D[t][pair],env[p,var]));ff.append(z);v=add(v,z)
  old3.append(a);old2.append(b);carry.append(v);factor_bounds.append(ff);partial.append(add(add(a,b),v))
 enc=lambda q:[str(q[0]),str(q[1])]
 pv=lambda v:[[i,*enc(q)] for i,q in sorted(v.items())]
 data={'selected_case619':549,'OLD_original768_source298_scale':'1/4','physical_B21_functions3':phys,'full_a1_32':s['actual_full_a1_27times9'],'Gamma_two':s['Gamma_two'],'D_two_by10_by56':D,'OLD_full_typed_program5':s['OLD_full_typed_program'],'OLD_Gamma_source_roots':s['OLD_Gamma_source_roots'],'OLD_D_source_roots':s['OLD_D_source_roots'],'known_FIRSTcarry_factor_recipes55':s['all_known_FIRSTcarry_factor_recipes'],'Xp_FIRST_column_recipes22':recipes,'native_L0_selected_columns':[[j,s['native_L0'][j]] for j in sorted(use0)],'native_Le_selected_columns':[[j,s['native_Le'][j]] for j in sorted(usee)],'positive_Xp_FIRST_columns22':[[p,var,pv(v)] for (p,var),v in sorted(env.items())],'positive_OLD_quarter_Gamma_cubic_two':list(map(enc,old3)),'positive_OLD_quarter_D_second_two':list(map(enc,old2)),'positive_FIRST55_factor_bounds_two_by55':[[enc(q) for q in v] for v in factor_bounds],'positive_FIRSTcarry_total_two':list(map(enc,carry)),'positive_OLD_and_FIRST_partial_two':list(map(enc,partial)),'future_signed_component_formula':'OLDpart=(Gamma X2983 +D X2982)/4; knownFIRSTcarry=SUM55(a1coef*Dpair*Xpfirst); fullONE549 addsnativeGamma3+Dnative2 andSUBTRACTSknownFIRSTcarry andacceptedGammaXp2a1.','native_individual_Xpell_is_actual_SUM_L0Se_LeS0':True,'all_first_sources_originalcurrent27_frame':True,'no_nativeS2_S3_or_outerXp2_values_here':True,'no_signedOLD_FIRST_or_upper_response_computed':True,'this_is_component_only_not_a_separate_VALUE_permission':True}
 n=count(data);scratch=count([env,L0,Le,ph,delta,coord])+6000;cap=2*(n+1000)+scratch+3000;rawcap=2_000_000+cap*(2*math.ceil(bits*math.log10(2))+50)
 assert cap<60000 and rawcap<32*1024*1024 and bits<=384 and 3*bits+1<=1153
 data.update({'source_committed_bits':bits,'source_temporary_bits':3*bits+1,'input_numeric_before_finalmeta':n,'future_scratch_upper':scratch,'future_numeric_cap':cap,'future_RAW_cap':rawcap,'metadata_reserve':1000})
 ep=out/'ONE549-OLD-quarter-and55-FIRST-positive-component-input.json';H.dump(ep,data)
 owned=count([s,phys,ph,fst,L0,Le,env,data,old3,old2,carry,partial,pins])+10000;assert owned<60000
 report={'id':STEM,'all_passed':True,'positiveOLD_and_FIRST_component_only':True,'no_signed_response_new_source_or_NATIVE_values':True,'all55_FIRST_factors_and22_pointrecipes_sourcebound':True,'all_individual_ell_SUMs_retained':True,'committed_prior_bits':bits,'temporary_prior_bits':3*bits+1,'source_owned_bound':owned,'future_component_numeric_bound':cap,'future_component_RAW_bound':rawcap,'whole_SOURCE_RAW_upper':raw,'combined_native_whole_gate_still_required':True,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','committed_prior_bits','temporary_prior_bits','source_owned_bound','future_component_numeric_bound','future_component_RAW_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

#!/usr/bin/env python3
"""Positive-only OLD two columns and all57 known-FIRST carries for ONE554."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER554_OLD_CARRY_INPUT_GATE_20261105_v63';Z=(1,0)
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
 def va(a,b,coefficient=None):
  z=dict(a)
  for i,q in b.items():z[i]=add(z.get(i,Z),q if coefficient is None else mul(q,coefficient))
  return z
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
 sp='reports/AMBIENT_B18_ELL2_T2_UPPER554_OLD_CARRY_SOURCE_20261105_v63/'
 sr=read(sp+'report.json');assert sr['all_passed'] and pins[sp+'report.json']=='1a5b1726d6592dbd67fe647cfbd9e65bd36997e615aa8f66696767d6a43ccc15'
 s=read(sp+'ONE554-six-c0-OLD-and-full29-a1-known-FIRST-source.json')
 assert len(s['all_known_FIRSTcarry_factor_recipes'])==57 and len(s['all_needed_Xp_FIRST_column_recipes'])==28
 allphi=field('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json','all_requested164_root_values');ph={i:dict(v) for i,v in allphi if i in s['OLD_physical_function_leaf_ids7']};del allphi
 first=field('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json','first640');fst={i:dict(v) for i,v in first['all79_atom_values_full10'] if i in ph};del first
 phys=[]
 for i in sorted(ph):
  v=[ph[i].get(k,'0') for k in s['raw33_to_B21_selector21']];assert list(map(F,v[:10]))==[F(fst[i].get(k,'0')) for k in range(10)];phys.append([i,v])
 assert len(phys)==7
 b=[tuple(v) for v in s['ordinary_basis_B21']];bi={v:i for i,v in enumerate(b)};table=[[bi.get(tuple(x+y for x,y in zip(a,c,strict=True))) for c in b] for a in b]
 def vm(a,b):
  z={}
  for i,x in a.items():
   for j,y in b.items():
    k=table[i][j]
    if k is not None:z[k]=add(z.get(k,Z),mul(x,y))
  return z
 pe={i:{k:ex(q) for k,q in enumerate(v) if F(q)} for i,v in phys};to={}
 for i,n in s['OLD_full_typed_program']:
  cursor={'stage':'positive_full11_OLD_program','node':i}
  if n[0]=='actual_C0_scalar_function':v=dict(pe[n[1]])
  elif n[0]=='sum':
   v={}
   for j,a,c in n[1]:v=va(v,to[j],ex(F(a,c)))
  else:
   assert n[0]=='product';v={0:ex(F(n[1],n[2]))} if n[1] else {}
   for j in n[3]:v=vm(v,to[j])
  to[i]=v
 L0={j:dict(v) for j,v in enumerate(s['native_L0'])};Le={j:dict(v) for j,v in enumerate(s['native_Le'])};env={};use0=set();usee=set();recipes=s['all_needed_Xp_FIRST_column_recipes']
 for rec in recipes:
  p,var=rec['pivot27'],rec['variable9'];r=rec['first_column_recipe'];cursor={'stage':'positive_individual_FIRST_column','pivot27':p,'variable9':var}
  if r['kind']=='literal56':v={i:ex(q) for i,q in r['column']}
  else:
   assert r['kind']=='native_ell_SUM' and var==8;a=dict(r['Se']);c=dict(r['S0']);use0.update(a);usee.update(c);v=va(mv(L0,a),mv(Le,c))
  env[p,var]=v
 source_pre=count([s,phys,ph,fst,L0,Le,env,to,table,pins])+20000;assert source_pre<60000
 D=s['D_two_by10_by56'];G=[dict(v) for v in s['Gamma_two']];old3=[];old2=[];carry=[];ff=[];partial=[]
 def C(i):return {} if i==0 else pe[i]
 for t in range(2):
  a=Z;c=Z;v=Z;factors=[]
  for old,scale,roots in s['OLD_Gamma_source_columns']:
   for row,(kind,i) in roots:
    q=to[i].get(20,Z) if kind=='typed' else C(i).get(20,Z);a=add(a,mul(ex(scale),mul(ex(G[t][row]),q)))
  for old,scale,roots in s['OLD_D_source_columns']:
   for row,i in roots:
    for pair in range(10):
     if D[t][pair][row]!='0':c=add(c,mul(ex(scale),mul(ex(D[t][pair][row]),to[i].get(10+pair,Z))))
  for k,(p,pair,var,q) in enumerate(s['all_known_FIRSTcarry_factor_recipes']):
   cursor={'stage':'positive_all57_FIRST_factors','channel':t,'factor':k};z=mul(ex(q),dot(D[t][pair],env[p,var]));factors.append(z);v=add(v,z)
  old3.append(a);old2.append(c);carry.append(v);ff.append(factors);partial.append(add(add(a,c),v))
 enc=lambda q:[str(q[0]),str(q[1])];pv=lambda v:[[i,*enc(q)] for i,q in sorted(v.items())]
 data={'selected_case619':554,'OLD_source2':s['OLD_source2'],'physical_B21_functions7':phys,'full_a1_29':s['actual_full_a1_27times9'],'Gamma_two':s['Gamma_two'],'D_two_by10_by56':D,'OLD_full_typed_program11':s['OLD_full_typed_program'],'OLD_Gamma_source_columns':s['OLD_Gamma_source_columns'],'OLD_D_source_columns':s['OLD_D_source_columns'],'known_FIRSTcarry_factor_recipes57':s['all_known_FIRSTcarry_factor_recipes'],'Xp_FIRST_column_recipes28':recipes,'native_L0_selected_columns':[[j,s['native_L0'][j]] for j in sorted(use0)],'native_Le_selected_columns':[[j,s['native_Le'][j]] for j in sorted(usee)],'positive_OLD_typed_B21':[[i,pv(v)] for i,v in sorted(to.items())],'positive_Xp_FIRST_columns28':[[p,var,pv(v)] for (p,var),v in sorted(env.items())],'positive_OLD_Gamma_cubic_two':list(map(enc,old3)),'positive_OLD_D_second_two':list(map(enc,old2)),'positive_FIRST57_factor_bounds_two_by57':[[enc(q) for q in v] for v in ff],'positive_FIRSTcarry_total_two':list(map(enc,carry)),'positive_OLD_and_FIRST_partial_two':list(map(enc,partial)),'future_signed_component_formula':'OLD=sum[(Gamma Xj3 +D Xj2)*OLDcoeff]; full addsnativeGamma3+Dnative2 thenSUBTRACTSknownFIRSTcarry andnonzeroacceptedGammaXp6442*a1plane.','all_individual_native_Xpell_SUM_L0Se_LeS0_retained':True,'all_physical_first_overlap_literal_current9_frame':True,'no_honest_higher_or_conductor_order_assumption_needed_in_OLDcomponent':True,'nonzero_p644_outer_target_owned_by_NATIVEpayload_not_dropped':True,'no_signed_OLD_FIRST_native_or_upper_response_computed':True,'component_only_fresh_whole_owner_gate_and_parentGO_required':True}
 n=count(data);scratch=count([env,to,L0,Le,ph,table])+6000;cap=2*(n+1000)+scratch+3000;rawcap=2_000_000+cap*(2*math.ceil(bits*math.log10(2))+50)
 assert cap<60000 and rawcap<32*1024*1024 and bits<=384 and 3*bits+1<=1153
 data.update({'source_committed_bits':bits,'source_temporary_bits':3*bits+1,'input_numeric_before_finalmeta':n,'future_scratch_upper':scratch,'future_numeric_cap':cap,'future_RAW_cap':rawcap,'metadata_reserve':1000})
 ep=out/'ONE554-OLD-two-columns-and57-FIRST-positive-component-input.json';H.dump(ep,data)
 owned=count([s,phys,ph,fst,L0,Le,env,to,table,data,old3,old2,carry,partial,pins])+10000;assert owned<60000
 report={'id':STEM,'all_passed':True,'positive_component_only_no_signed_response':True,'allOLD2_sources_and57_FIRST_factors_retained':True,'all28_individual_first_recipes_SUMbound':True,'source_prearith_bound':source_pre,'source_owned_bound':owned,'committed_prior_bits':bits,'temporary_prior_bits':3*bits+1,'future_component_numeric_bound':cap,'future_component_RAW_bound':rawcap,'whole_SOURCE_RAW_upper':raw,'whole_native_outer_coupling_owner_gate_pending':True,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','committed_prior_bits','temporary_prior_bits','source_owned_bound','future_component_numeric_bound','future_component_RAW_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

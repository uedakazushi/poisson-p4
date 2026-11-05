#!/usr/bin/env python3
"""Positive/cache input for one complete opaque-free upper508 response."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER508_INPUT_GATE_20261105_v63';Z=(1,0)
def main(out):
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);pins={};raw=0;bits=0;cursor={}
 def read(p):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024
  return json.loads((ROOT/p).read_text())
 def field(p,k):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024
  return H.selected(p,k)
 def ck(v):
  nonlocal bits
  bits=max(bits,abs(v).bit_length())
  if abs(v).bit_length()>384:
   H.dump(out/'BEFORE_RESPONSE_POSITIVE_GUARD.json',{'cursor':cursor,'uncommitted_positive_integer':str(v),'no_actual_response_values':True});raise RuntimeError('positive384 STOP')
  return v
 def ex(q):q=F(q);return (ck(q.denominator),ck(abs(q.numerator))) if q else Z
 def add(a,b):
  if not a[1]:return b
  if not b[1]:return a
  den=ck(math.lcm(a[0],b[0]));u=ck(a[1]*(den//a[0]));v=ck(b[1]*(den//b[0]));return den,ck(u+v)
 def mul(a,b):return (ck(a[0]*b[0]),ck(a[1]*b[1])) if a[1] and b[1] else Z
 def va(a,b,coef=None):
  z=dict(a)
  for i,q in b.items():z[i]=add(z.get(i,Z),q if coef is None else mul(q,coef))
  return z
 sp='reports/AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63/'
 sr=read(sp+'report.json');assert sr['all_passed'] and pins[sp+'report.json']=='7b2bd2588c1d7f1cb30a8c6c7822a80d1a770417cfd114d324531f4925a93c22'
 s=read(sp+'ONE508-source-closure-B21-and-exact-carry-zeros.json');assert s['a1_all27times9_exact_empty'] and s['native_c0_exact_empty']
 assert read('reports/AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63_REVIEW/report.json')['all_passed']
 jsrc='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_FORMAL_20261105_v63/617-formal-direct-old-third-and-second-symbol-sensitivities.json'
 jr=next(q for q in read(jsrc)['all617_formal_directOLD_J0_Jr_Jell'] if q['free_case619']==508)
 pr='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json';ph={i:dict(v) for i,v in field(pr,'all_requested164_root_values')}
 fp='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json';fst={i:dict(v) for i,v in field(fp,'first640')['all79_atom_values_full10']}
 leaves=s['function_leaf39_ids'];sel=s['raw33_to_B21_selector21'];physical=[];arbitrary=[];leafvalues={};catalog=[]
 for i in leaves:
  assert i in fst
  first=[fst[i].get(k,'0') for k in range(10)];assert first[0]=='0'
  if i in ph:
   v=[ph[i].get(k,'0') for k in sel];assert list(map(F,v[:10]))==list(map(F,first));physical.append(i);kind='original_physical164_actual_B21'
  else:
   # This is a whole-functional extension, not an actual higher-jet claim.
   v=first+['0']*11;arbitrary.append(i);kind='whole_functional_higher0_extension'
  leafvalues[i]=v;catalog.append([i,kind,v])
 assert len(physical)==23 and len(arbitrary)==16 and {299650,505938}<=set(arbitrary)
 for i in [299650,505938]:assert all(F(q)==0 for q in leafvalues[i][:10])
 for ch in jr['two_formal_J0_Jr_Jell_symbol_rows']:
  for i,v in ch:
   if i in arbitrary:assert i in [299650,505938] and F(v[0])==0
 dform=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63/actual617-formal-combined-OLD2-second-symbol-coefficients.json');dc=dict(dform['all617_combined_OLD2_independent_symbol_coefficients'])[508]
 for ch in dc:
  for p,i,q in ch:
   if i in arbitrary:assert i in [299650,505938] and p in [8,9]
 for p in ['reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_FORMAL_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_OLD_CASE508_ORDER_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json']:
  assert read(p)['all_passed']
 qpath='reports/AMBIENT_B18_ELL2_T2_UPPER617_FULL_COORDINATE_SOURCE_20261105_v63/known-Psi-point-correction-and-full-b1-source-annihilation.json';Q=field(qpath,'Q_two6')
 assert len(Q)==2 and all(len(dict(v)) in [3,6] for v in Q)
 source_pre=count([s,jr,ph,fst,leafvalues,catalog,dform,Q,pins])+15000;assert source_pre<150000
 del ph,fst,dform
 b=[tuple(v) for v in s['ordinary_basis_B21']];bi={m:i for i,m in enumerate(b)}
 table=[[bi.get(tuple(x+y for x,y in zip(a,c,strict=True))) for c in b] for a in b]
 def vm(a,b):
  z={}
  for i,x in a.items():
   for j,y in b.items():
    k=table[i][j]
    if k is not None:z[k]=add(z.get(k,Z),mul(x,y))
  return z
 leafenv={i:{j:ex(q) for j,q in enumerate(v) if F(q)} for i,v in leafvalues.items()};co={};to={}
 def C(i):return co[i] if i in co else leafenv[i]
 for i,n in s['literal_C0_program38']:
  cursor={'stage':'positive_C0_B21','node':i}
  if n[0]=='c':v={0:ex(F(n[1],n[2]))} if n[1] else {}
  elif n[0]=='n':v=dict(C(n[1]))
  elif n[0]=='a':v=va(C(n[1]),C(n[2]))
  else:assert n[0]=='m';v=vm(C(n[1]),C(n[2]))
  co[i]=v
 for i,n in s['full_typed_program86']:
  cursor={'stage':'positive_full_typed_B21','node':i}
  if n[0]=='actual_C0_scalar_function':v=dict(C(n[1]))
  elif n[0]=='sum':
   v={}
   for j,a,c in n[1]:v=va(v,to[j],ex(F(a,c)))
  else:
   assert n[0]=='product';v={0:ex(F(n[1],n[2]))} if n[1] else {}
   for j in n[3]:v=vm(v,to[j])
  to[i]=v
 G=[dict(v) for v in s['Gamma_two']];D=s['D_two_by10_by56'];gparts=[];dparts=[];total=[]
 for ch in range(2):
  a=Z;u=Z
  for row,(kind,i) in s['Gamma_old595_source_roots']:
   a=add(a,mul(ex(G[ch][row]),(to[i] if kind=='typed' else C(i)).get(20,Z)))
  for row,i in s['D_old595_source_roots']:
   for pair in range(10):
    if D[ch][pair][row]!='0':u=add(u,mul(ex(D[ch][pair][row]),to[i].get(10+pair,Z)))
  gparts.append(a);dparts.append(u);total.append(add(a,u))
 six=[]
 for k in range(6):
  v=Z
  for ch in range(2):v=add(v,mul(ex(dict(Q[ch]).get(k,'0')),total[ch]))
  six.append(v)
 encode=lambda q:[str(q[0]),str(q[1])]
 data={'selected_case619':508,'original768_source':[[595,'1']],'exact_a1_and_native_zero_source':True,'ordinary_basis_B21':s['ordinary_basis_B21'],'physical164_actual_B21_and_only79_truefirst_extension_catalog39':catalog,'physical_function_ids23':physical,'arbitrary_higher_extension_ids16':arbitrary,'actual_conductor_second_scope_only_rr_rell_zero':True,'both_conductor_third_symbol_coefficients0_IN_THIS_CASE':True,'all_other_honest_second_third_symbol_coefficients0_IN_THIS_CASE':True,'original_case508_J1_certificate':jr,'original_case508_OLD2_certificate':dc,'literal_C0_program38':s['literal_C0_program38'],'full_typed_program86':s['full_typed_program86'],'Gamma_two':s['Gamma_two'],'D_two_by10_by56':D,'Gamma_old595_source_roots':s['Gamma_old595_source_roots'],'D_old595_source_roots':s['D_old595_source_roots'],'fixed_Q_embedding_two6':Q,'positive_C0_B21':[[i,[[k,*encode(q)] for k,q in sorted(v.items())]] for i,v in co.items()],'positive_typed_B21':[[i,[[k,*encode(q)] for k,q in sorted(v.items())]] for i,v in to.items()],'positive_Gamma_direct_cubic_two':list(map(encode,gparts)),'positive_D_full_firstcross_second_two':list(map(encode,dparts)),'positive_actual_upper_two':list(map(encode,total)),'positive_full_six':list(map(encode,six)),'no_actual_response_value_computed':True,'source_formula':'upper_t=Gamma_t X5953[r²ell]+D_t X5952; fullC/Psi alreadyinD/Gamma; noXp*a1/native becausebothFULLsourcevectors0.','extension_scope':'Only COMPLETEfunctional usesarbitraryhigher0 for16honest/conductorleaves, aftertwopeeredformalcertificates. Conductorr/resecondsareactual0byplaneorder3; nootheractualhigherzeroassertion.'}
 n=count(data);scratch=count([leafenv,co,to,table,G,gparts,dparts,total,six])+5000;success=2*(n+1000)+scratch+3000;rawupper=2_000_000+success*(2*math.ceil(bits*math.log10(2))+50)
 assert success<60000 and rawupper<32*1024*1024 and bits<=384 and 3*bits+1<=1153
 data.update({'source_committed_bits':bits,'source_temporary_bits':3*bits+1,'input_numeric_before_finalmeta':n,'future_scratch_upper':scratch,'future_numeric_cap':success,'future_RAW_cap':rawupper,'metadata_reserve':1000})
 ep=out/'ONE508-positive-full-cache-input-and-honest-extension-certificate.json';H.dump(ep,data)
 owned=count([s,jr,leafvalues,catalog,co,to,leafenv,data,Q,pins])+10000;assert owned<150000
 report={'id':STEM,'all_passed':True,'only_positive_known_cache_and_certificate_inputs':True,'no_actual_conductor_or_other_response':True,'all_full39firsts_and23physical_overlap_verified':True,'source_prearith_bound':source_pre,'source_owned_bound':owned,'committed_prior_bits':bits,'temporary_prior_bits':3*bits+1,'future_numeric_cap':success,'future_RAW_cap':rawupper,'whole_SOURCE_RAW_reads_upper':raw,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','committed_prior_bits','temporary_prior_bits','source_owned_bound','future_numeric_cap','future_RAW_cap','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

#!/usr/bin/env python3
"""SOURCE positive shared OLD/FIRST/outer/native-fine/Q bounds, no values."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER605_POSITIVE_WHOLE_20261105_v63';Z=(1,0)
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 pins={};raw=0;cursor={};pending={};maximum=0;C={};T={};first={};table_first={};rawbounds={};unions={}
 def read(p):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024;return json.loads((ROOT/p).read_text())
 def ck(n):
  nonlocal maximum
  maximum=max(maximum,abs(n).bit_length())
  if abs(n).bit_length()>384:cursor.update(offending_uncommitted_integer=str(n),bits=abs(n).bit_length());raise RuntimeError('POSITIVE384_STOP')
  return n
 def binary(a,b,kind):
  nonlocal pending
  pending={'kind':kind,'a':str(a),'b':str(b),'cursor':dict(cursor)}
  z=ck(a*b if kind=='mul' else a+b if kind=='add' else math.lcm(a,b));pending={};return z
 def ex(v):
  q=F(v);return (ck(q.denominator),ck(abs(q.numerator))) if q else Z
 def add(a,b):
  if not a[1]:return b
  if not b[1]:return a
  d=binary(a[0],b[0],'lcm');u=binary(a[1],d//a[0],'mul');v=binary(b[1],d//b[0],'mul');return d,binary(u,v,'add')
 def union(a,b):
  if not a[1]:return b
  if not b[1]:return a
  d=binary(a[0],b[0],'lcm');u=binary(a[1],d//a[0],'mul');v=binary(b[1],d//b[0],'mul');return d,max(u,v)
 def mul(a,b):return (binary(a[0],b[0],'mul'),binary(a[1],b[1],'mul')) if a[1] and b[1] else Z
 def va(a,b,c=(1,1)):
  z=dict(a)
  for k,q in b.items():z[k]=add(z.get(k,Z),mul(c,q))
  return z
 def dot(a,b):
  z=Z
  for i,q in a.items():
   if i in b:z=add(z,mul(q,b[i]))
  return z
 def save_union(name,t,q):unions[name,t]=union(unions.get((name,t),Z),q)
 try:
  sp='reports/AMBIENT_B18_ELL2_T2_UPPER605_OLD_FIRST_CACHE_SOURCE_20261105_v63/';gr=read(sp+'report.json');assert gr['all_passed'];s=read(sp+'pruned605-OLD757-C0532-B21-FIRST162-and-native92-source-input.json')
  bp='reports/AMBIENT_B18_ELL2_T2_UPPER605_BOOLEAN_PROFILE_20261105_v63/';br=read(bp+'report.json');assert br['all_passed'];bo=read(bp+'Boolean-shared-program-supports-FIRST324-and-literal-LCM-classes.json')
  np='reports/AMBIENT_B18_NATIVE_SHARED72_FINE_POSITIVE_20261105_v63/';nr=read(np+'report.json');assert nr['all_passed'];nf=read(np+'native72-fine-Gamma-D-and605-constantmap-positive-envelopes.json')
  native=read('reports/AMBIENT_B18_NATIVE_SHARED72_INPUT_SOURCE_20261105_v63/native-shared75-formal-to72-actual-to605-source-input.json')
  npeer=read('reports/AMBIENT_B18_NATIVE_SHARED72_SOURCE_FINE_20261105_v63_REVIEW/report.json');assert npeer['all_passed']
  qpath='reports/AMBIENT_B18_ELL2_T2_UPPER617_FULL_COORDINATE_SOURCE_20261105_v63/known-Psi-point-correction-and-full-b1-source-annihilation.json';pins[qpath]=H.sha(ROOT/qpath);raw+=2*(ROOT/qpath).stat().st_size;Q=H.selected(qpath,'Q_two6')
  assert count([s,bo,nf,native,gr,br,nr,npeer,Q])+30000<100000
  b=[tuple(x) for x in s['ordinary_B21_basis']];bi={x:i for i,x in enumerate(b)};table=[[bi.get(tuple(x+y for x,y in zip(a,c))) for c in b] for a in b]
  def vm(a,b):
   z={}
   for i,x in a.items():
    for j,y in b.items():
     k=table[i][j]
     if k is not None:z[k]=add(z.get(k,Z),mul(x,y))
   return z
  leaves={fn:{i:ex(q) for i,q in enumerate(v) if F(q)} for fn,mode,v in s['actual_or_formal_B21_boundary96']};assert sum(map(len,leaves.values()))==152
  C=dict(leaves)
  for i,d in s['shared_literal_C0_program532']:
   cursor={'stage':'positive_C0532','node':i}
   if d[0]=='c':v={0:ex(F(d[1],d[2]))} if d[1] else {}
   elif d[0]=='n':v=dict(C[d[1]])
   elif d[0]=='a':v=va(C[d[1]],C[d[2]])
   else:assert d[0]=='m';v=vm(C[d[1]],C[d[2]])
   C[i]=v
  assert sum(len(C[i]) for i,d in s['shared_literal_C0_program532'])==119
  G=[{i:ex(q) for i,q in row} for row in s['Gamma_two']];D=[[{i:ex(q) for i,q in enumerate(row) if F(q)} for row in ch] for ch in s['D_two_by10_by56']]
  gc=dict(s['OLD_Gamma_source_columns618']);dc=dict(s['OLD_D_source_columns618']);old3={};old2={}
  for old,roots in gc.items():
   v=[Z,Z]
   for row,(kind,node) in roots:
    if kind!='C0':continue
    for t in range(2):
     if row in G[t]:v[t]=add(v[t],mul(G[t][row],C[node].get(20,Z)))
   old3[old]=v
  Ccells=sum(len(C[i]) for i,d in s['shared_literal_C0_program532']);C.clear() # phase-disjoint lifetime
  for i,d in s['shared_full_typed_program757']:
   cursor={'stage':'positive_typed757','node':i}
   if d[0]=='actual_C0_scalar_function':v=dict(leaves[d[1]])
   elif d[0]=='sum':
    v={}
    for k,a,den in d[1]:v=va(v,T[k],ex(F(a,den)))
   else:
    assert d[0]=='product';v={0:ex(F(d[1],d[2]))} if d[1] else {}
    for k in d[3]:v=vm(v,T[k])
   T[i]=v
  assert sum(map(len,T.values()))==324
  for old,roots in gc.items():
   for row,(kind,node) in roots:
    if kind!='typed':continue
    for t in range(2):
     if row in G[t]:old3[old][t]=add(old3[old][t],mul(G[t][row],T[node].get(20,Z)))
  for old,roots in dc.items():
   v=[Z,Z]
   for row,node in roots:
    for t in range(2):
     for pair in range(10):
      if row in D[t][pair]:v[t]=add(v[t],mul(D[t][pair][row],T[node].get(10+pair,Z)))
   old2[old]=v
  T.clear()
  L0={j:{i:ex(q) for i,q in row} for j,row in s['native_L0_selected_columns']};Le={j:{i:ex(q) for i,q in row} for j,row in s['native_Le_selected_columns']}
  for rec in s['individual_FIRST_column_recipes162']:
   p,var=rec['pivot27'],rec['variable9'];d=rec['first_column_recipe'];cursor={'stage':'positive_FIRST162','pivot':p,'var':var}
   if d['kind']=='literal56':v={i:ex(q) for i,q in d['column'] if F(q)}
   else:
    assert d['kind']=='native_ell_SUM';v={}
    for name,mat in [('Se',L0),('S0',Le)]:
     for k,q in d[name]:v=va(v,mat[k],ex(q))
   first[p,var]=v
  assert sum(map(len,first.values()))==1295
  pairs=s['fixed_pair_order10']
  for t in range(2):
   for p in range(27):
    for var in [0,1,4,5,7,8]:
     cursor={'stage':'positive_FIRST324','channel':t,'pivot':p,'var':var};z=Z
     for pair,(a,b) in enumerate(pairs):
      if var==a:z=add(z,dot(D[t][pair],first[p,b]))
      if var==b and a!=b:z=add(z,dot(D[t][pair],first[p,a]))
     table_first[t,p,var]=z
  source92=s['accepted_native92_original_sources23'];i92={orig:p for p,orig in enumerate(source92)};xp92=s['accepted_native92_Gamma_plane_targets23_twochannels_rr_rell'];piv=s['point27_source_indices']
  for case,(c0,a1) in enumerate(zip(s['raw619_original768_constant_sources'],s['raw619_actual_a1_27times9'])):
   rawbounds[case]=[]
   for t in range(2):
    cursor={'stage':'positive_raw619_FORMALextension_partial','case':case,'channel':t};v=[Z,Z,Z,Z]
    for col,q in c0:
     if col<640:v[0]=add(v[0],mul(ex(q),old3[col][t]));v[1]=add(v[1],mul(ex(q),old2[col][t]))
    for slot,q in a1:
     p,var=divmod(slot,9)
     if var in [0,1,4,5,7,8]:v[2]=add(v[2],mul(ex(q),table_first[t,p,var]))
     else:assert var in [2,3,6]
     if p>=4 and var in [7,8]:
      direction=0 if var==8 else 1;v[3]=add(v[3],mul(ex(q),ex(xp92[i92[piv[p]]][t][direction])))
    rawbounds[case].append(v)
    for k,name in enumerate(['raw_OLD3','raw_OLD2','raw_FIRST','raw_nativeOuter']):save_union(name,t,v[k])
  ncase={i:[tuple(map(int,q)) for q in row] for i,row in nf['all605_native_sum_fine_twochannels']}
  assert set(ncase)==set(s['cases605']);case_positive_nonzero=[]
  for case,(a,beta) in s['corrected605_recipes_basis563609']:
   nonzero=False
   for t in range(2):
    cursor={'stage':'positive_corrected605_full_assembly','case':case,'channel':t};v=list(rawbounds[case][t])
    for base,q in [(563,a),(609,beta)]:
     if F(q):v=[add(x,mul(ex(q),y)) for x,y in zip(v,rawbounds[base][t])]
    for k,name in enumerate(['OLD3','OLD2','FIRST','nativeOuter']):save_union(name,t,v[k])
    save_union('native',t,ncase[case][t]);z=ncase[case][t]
    for x in v:z=add(z,x)
    save_union('upper',t,z);nonzero|=bool(z[1])
   if nonzero:case_positive_nonzero.append(case)
  qdense=[dict(row) for row in Q]
  for i in range(6):
   z=Z
   for t in range(2):z=add(z,mul(ex(qdense[t].get(i,'0')),unions['upper',t]))
   save_union('six',i,z)
  # All numeric coefficient operations above are positive envelopes only.
  # The runtime must compare actual native output to its FINE bound before
  # using the upper/Q assembly, never multiply intermediate354 byQ.
  union_prior=max(maximum,354);assert union_prior<=384;temporary=3*union_prior+1;assert temporary<=1153
  enc=lambda q:[str(q[0]),str(q[1])]
  payload={'SOURCE_positive_ONLY_no_signed_response':True,'OLDcache_source_report':sp+'report.json','native_cache_source_report':'reports/AMBIENT_B18_NATIVE_SHARED72_INPUT_SOURCE_20261105_v63/report.json','native_FINE_source_report':np+'report.json','all605_case_coverage_checked':s['cases605'],'positive_common_envelopes':[[name,t,enc(q)] for (name,t),q in sorted(unions.items())],'native72_fine_bounds':nf['native_kernel72_Gamma_D_sum_fine'],'structural_positive_nonzero_case_labels':case_positive_nonzero,'positive_OLD_compiler_cell_counts':[Ccells,324],'positive_FIRST_table_cells324':True,'native354_intermediate_and_FINE_before_Q_rule':True,'positive_whole_operation_bits':maximum,'future_value_committed_prior':union_prior,'future_value_temporary_prior':temporary,'unchanged_hard_guards':[384,1153],'raw619_partial_values_are_FORMAL_extensions_until_corrected605':True,'allOLD4_opaque24_keptformal_and_exactlyannihilated_ONLYafter605correction':True,'future_native_raw75_FORMAL_kernel72_ACTUAL_scope_retained':True,'whole_immutableworker_lifetime_codec_and_freshpeer_ROOTGO_STILLrequired':True}
  ep=out/'positive-shared605-OLD-FIRST-nativeOuter-nativeFINE-and-Q-common-bounds.json';H.dump(ep,payload)
  owned=2*count([s,bo,nf,native,gr,br,nr,npeer,Q])+count([leaves,G,D,L0,Le,first,table_first,old3,old2,rawbounds,unions,payload,table])+15000;assert owned<150000
  report={'id':STEM,'all_passed':True,'only_nonnegative_SOURCEbounds_no_actual605_values':True,'all605_complete_source_bounds':True,'source_owned_numeric_upper':owned,'RAW_read_upper':raw,'positive_whole_operation_bits':maximum,'future_value_committed_prior':union_prior,'future_value_temporary_prior':temporary,'native_FINE_checked_before_upperQ':True,'whole_VALUE_immutablecode_custody_peer_and_ROOTGO_pending':True,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','positive_whole_operation_bits','future_value_committed_prior','future_value_temporary_prior','source_owned_numeric_upper','seconds']}))
 except BaseException as error:
  H.dump(out/'BEFORE_RESPONSE_POSITIVE_STOP.json',{'error':repr(error),'cursor':cursor,'pending_integer_operation':pending,'C0_current':C,'typed_current':T,'FIRST_current':first,'FIRST_table':[[list(k),list(v)] for k,v in table_first.items()],'raw_positive_bound_prefix':rawbounds,'whole_positive_unions':[[list(k),list(v)] for k,v in unions.items()],'no_signed_response_values':True,'no_auto_retry_or_capraise':True});raise
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())

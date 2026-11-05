#!/usr/bin/env python3
"""Prepared single20 cache bundle; fresh SOURCE peer/custody/ROOT GO required."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER20_VALUES_20261105_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def count(v):
 if isinstance(v,(F,int)):return 1
 if isinstance(v,str):
  try:F(v);return 1
  except (ValueError,ZeroDivisionError):return 0
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple,set)):return sum(map(count,v))
 return 0
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x) for k,x in v.items()}
 if isinstance(v,(tuple,list,set)):return list(map(enc,v))
 return v
def dump(p,v):p.write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
def main(out,gate,peer,custody):
 def timeout(signum,frame):raise TimeoutError('WHOLE60_STOP')
 start=time.monotonic();signal.signal(signal.SIGALRM,timeout);signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 pins={};raw=0;cursor={};pending={};work={};maximum=0;tempmax=0;opindex=0
 C={};T={};first={};aux={};points={};kernels={};rawparts={};FIRSTfactors={};old3={};old2={};classes={};results=[];integer=[];factors=[];kernelcert=[];numerators=[];selected=[];rows=[];minor=0;rank=None
 def pin(p):
  nonlocal raw
  pp=ROOT/p;raw+=pp.stat().st_size;pins[p]=sha(pp);return pp
 def read(p):return json.loads(pin(p).read_text())
 def ck(x,temporary=False):
  nonlocal maximum,tempmax
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length()) if isinstance(x,F) else abs(x).bit_length()
  tempmax=max(tempmax,b)
  if b>(1063 if temporary else 354):cursor.update(offending_uncommitted_value=str(x),bits=b);raise RuntimeError('TEMP1063_STOP' if temporary else 'COMMIT354_STOP')
  if not temporary:maximum=max(maximum,b)
  return x
 def parse(q):return ck(F(q))
 def binary(a,b,kind,temporary=False):
  nonlocal pending,opindex,tempmax
  opindex+=1;cursor['binary_index']=opindex;pending={'operation':kind,'a':str(a),'b':str(b),'temporary_only':temporary,'cursor':dict(cursor)}
  bit=lambda x:max(abs(x.numerator).bit_length(),x.denominator.bit_length()) if isinstance(x,F) else abs(x).bit_length()
  operand_bound=bit(a)+bit(b)+(kind=='add') if isinstance(a,F) or isinstance(b,F) or kind=='mul' else max(bit(a),bit(b))+1
  tempmax=max(tempmax,operand_bound);assert operand_bound<=1063
  z=ck(a+b if kind=='add' else a*b,temporary);pending={};return z
 def add(a,b):return binary(a,b,'add')
 def mul(a,b):return binary(a,b,'mul')
 def neg(x):return ck(-x)
 def put(v,k,q):
  z=add(v.get(k,F(0)),q)
  if z:v[k]=z
  else:v.pop(k,None)
 def dot(a,b):
  z=F(0);work['dot_partial']=z
  for i,x in a.items():
   cursor['dot_row']=i
   if i in b:z=add(z,mul(x,b[i]));work['dot_partial']=z
  return z
 def env(q,b):
  D,U=map(int,b);assert D%q.denominator==0;assert ck(abs(q.numerator)*(D//q.denominator),True)<=U
 try:
  cu=read(custody);pr=read(peer);g=read(gate);assert cu['all_passed'] and cu['safe_help_literal_PASS'] and pr['all_passed'] and g['all_passed']
  assert sha(Path(__file__))==cu['worker_sha256'] and pins[gate]==cu['source_gate_sha256'] and pins[gate] in json.dumps(pr)
  for p,h in cu['worker_and_source_sha256'].items():assert sha(pin(p))==h
  ip=g['value_input_path'];d=read(ip);assert pins[ip]==g['value_input_sha256']
  assert d['committed_prior']==354 and d['temporary_prior']==1063 and d['hard_guards']==[384,1153]
  cap=d['whole_numeric_cap'];rawcap=d['whole_RAW_cap'];assert cap<60000 and rawcap<32*1024*1024
  assert count([cu,pr,g,pins])<d['metadata_reserve']
  o=d['OLD'];n=d['NATIVE'];E={(name,t):b for name,t,b in d['positive_global_envelopes']};case20=o['exceptional_cases20'];assert len(case20)==20
  dump(out/'phase.json',{'phase':'ROOT_authorized_ONE20_cache_bundle','source_gate':pins[gate],'source_peer':pins[peer],'custody':pins[custody],'605coverage_sourcezero584_plus_accepted554_plus20':True,'no_new_physical_gauge_IFT_jets':True})
  Gamma=[{i:parse(q) for i,q in row} for row in o['Gamma_two']];D=[[{i:parse(q) for i,q in enumerate(row) if F(q)} for row in ch] for ch in o['D_two_by10_by56']]
  assert all(i<48 for ch in D for row in ch for i in row)
  # All leaf jets here are either physical cache or a certified arbitrary
  # extension. No individual only79 higher jet is called an actual value.
  leaves={fn:{i:parse(q) for i,q in enumerate(v) if F(q)} for fn,mode,v in o['boundary_cache']};b=[tuple(x) for x in o['ordinary_B21_basis']];bi={x:i for i,x in enumerate(b)}
  table=[[bi.get(tuple(x+y for x,y in zip(a,c))) for c in b] for a in b]
  def va(v,w,c):
   work['polynomial_partial']=v
   for i,x in w.items():cursor['polynomial_add_index']=i;put(v,i,mul(c,x))
  def vm(a,b):
   z={};work['polynomial_product_partial']=z
   for i,x in a.items():
    for j,y in b.items():
     k=table[i][j];cursor['polynomial_product_indices']=[i,j]
     if k is not None:put(z,k,mul(x,y))
   return z
  C=dict(leaves)
  for i,desc in o['literal_C0_program']:
   cursor={'stage':'shared_C0_FORMALextension','node':i}
   if desc[0]=='c':v={0:parse(F(desc[1],desc[2]))} if desc[1] else {}
   elif desc[0]=='n':v={k:neg(q) for k,q in C[desc[1]].items()}
   elif desc[0]=='a':v=dict(C[desc[1]]);va(v,C[desc[2]],F(1))
   else:assert desc[0]=='m';v=vm(C[desc[1]],C[desc[2]])
   C[i]=v
  gc=dict(o['OLD_Gamma_columns']);dc=dict(o['OLD_D_columns']);old3={j:[F(0),F(0)] for j in gc};old2={j:[F(0),F(0)] for j in dc}
  for j,roots in gc.items():
   for row,(kind,node) in roots:
    if kind=='C0':
     for t in range(2):
      if row in Gamma[t]:cursor={'stage':'OLD_C0_cubic_FORMALcolumn','source':j,'channel':t,'root':node};old3[j][t]=add(old3[j][t],mul(Gamma[t][row],C[node].get(20,F(0))))
  C.clear();work.clear()
  for i,desc in o['typed_program']:
   cursor={'stage':'shared_typed_FORMALextension','node':i}
   if desc[0]=='actual_C0_scalar_function':v=dict(leaves[desc[1]])
   elif desc[0]=='sum':
    v={}
    for k,num,den in desc[1]:va(v,T[k],parse(F(num,den)))
   else:
    assert desc[0]=='product';v={0:parse(F(desc[1],desc[2]))} if desc[1] else {}
    for k in desc[3]:v=vm(v,T[k])
   T[i]=v
  for j,roots in gc.items():
   for row,(kind,node) in roots:
    if kind=='typed':
     for t in range(2):
      if row in Gamma[t]:cursor={'stage':'OLD_typed_cubic_FORMALcolumn','source':j,'channel':t,'root':node};old3[j][t]=add(old3[j][t],mul(Gamma[t][row],T[node].get(20,F(0))))
  for j,roots in dc.items():
   for row,node in roots:
    for t in range(2):
     for pair in range(10):
      if row in D[t][pair]:cursor={'stage':'OLD_second_FORMALcolumn','source':j,'channel':t,'pair':pair,'root':node};old2[j][t]=add(old2[j][t],mul(D[t][pair][row],T[node].get(10+pair,F(0))))
  T.clear();work.clear()
  L0={j:{i:parse(q) for i,q in row} for j,row in o['native_L0_selected_columns']};Le={j:{i:parse(q) for i,q in row} for j,row in o['native_Le_selected_columns']}
  for rec in o['FIRST_column_recipes']:
   key=(rec['pivot27'],rec['variable9']);r=rec['first_column_recipe'];v={};first[key]=v;cursor={'stage':'actual_individual_FIRST','pivot_var':key}
   if r['kind']=='literal56':v={i:parse(q) for i,q in r['column']};first[key]=v
   else:
    assert r['kind']=='native_ell_SUM'
    for name,mat in [('Se',L0),('S0',Le)]:
     for j,q in r[name]:
      for i,x in mat[j].items():cursor.update(SUM_branch=name,input=j,row=i);put(v,i,mul(parse(q),x))
  pairs=o['pair_order10'];targets92={orig:v for orig,v in o['accepted_native_Gamma_plane_targets']};piv=o['point27_source_indices'];c0=dict(o['original_source_allc0_raw21']);a1=dict(o['original_full_a1_raw21']);FIRSTfactors={}
  for case in o['required_raw_cases21']:
   rawparts[case]=[];FIRSTfactors[case]=[]
   for t in range(2):
    cursor={'stage':'raw21_FORMAL_partial_NOT_actualbasisupper','case':case,'channel':t};v=[F(0)]*4;rawparts[case].append(v);fv=[];FIRSTfactors[case].append(fv)
    for j,q in c0[case]:
     if j<640:v[0]=add(v[0],mul(parse(q),old3[j][t]));v[1]=add(v[1],mul(parse(q),old2[j][t]))
    for slot,q in a1[case]:
     p,var=divmod(slot,9)
     for pair,(a,beta) in enumerate(pairs):
      for other in ([beta] if var==a else [])+([a] if var==beta and a!=beta else []):
       cursor.update(a1_slot=slot,pair=pair,firstvar=other);z=mul(parse(q),dot(D[t][pair],first[p,other]));fv.append(z);v[2]=add(v[2],z)
     if p>=4 and var in [7,8]:
      direction=0 if var==8 else 1;cursor.update(a1_slot=slot,native92_direction=direction);v[3]=add(v[3],mul(parse(q),parse(targets92[piv[p]][t][direction])))
    for k,name in enumerate(['raw_OLD3','raw_OLD2','raw_FIRST','raw_nativeOuter']):env(v[k],E[name,t])
  first.clear();work.clear()
  # Native raw auxiliary outputs are FORMAL partial-N extensions. Only
  # fixed kernel restrictions below are genuine native functionals.
  p=n['primitive'];basis=[tuple(x) for x in p['basis18']];idx={x:i for i,x in enumerate(basis)};known=p['N_known_positions'];assert known==[0,1,2,3,4,5,6,15,16]
  N=[[[[parse(q) for q in cell] for cell in row] for row in src] for src in p['N_known_first8_by12_by9_by4']];G=[[parse(q) for q in row] for row in p['G8_C72']];nu=[[parse(q) for q in row] for row in p['nu_two_J18']]
  assert all(tuple(a+b for a,b in zip(basis[i],basis[j])) not in idx for x in nu for y in nu for i,a in enumerate(x) if a for j,b in enumerate(y) if b)
  needed=sorted({i//4 for row in Gamma for i in row if i<48}|{i//4 for ch in D for row in ch for i in row});assert len(needed)==7
  def dest(a,b):return idx.get(tuple(x+y for x,y in zip(basis[a],basis[b])))
  def route(v,base,mask,im,x,y):
   collision=mask&im;om=mask^im
   if collision==3:return
   if not collision:put(v,4*base+om,mul(x,y));return
   for nb,z in enumerate(nu[collision-1]):
    if z:
     target=dest(base,nb)
     if target is not None:cursor['monicbase']=nb;put(v,4*target+om,mul(mul(x,y),z))
  for rec in n['raw_program']:
   j=rec['native128'];Sc={};NX={};work.update(current_raw_native=j,Sc=Sc,NX=NX);im=rec['inputmask']
   for entry,(row,gen,sign) in enumerate(rec['Koszul_entries']):
    sv=Sc.setdefault(row,{})
    for flat,x in enumerate(G[gen]):
     if x:cursor={'stage':'FORMAL_raw_intrinsic_S','raw':j,'entry':entry,'Gflat':flat};route(sv,flat//4,flat%4,im,x,F(sign))
   points[j]={4*row+k:x for row,sv in Sc.items() for k,x in sv.items() if k<4}
   assert points[j]=={i:parse(q) for i,q in dict(n['raw_point_S0_columns48'])[j]}
   for physical in needed:
    v=NX.setdefault(physical,{})
    for row,sv in Sc.items():
     for kp,nb in enumerate(known):
      for mask,x in enumerate(N[row][physical][kp]):
       if not x:continue
       for flat,y in sv.items():
        target=dest(nb,flat//4)
        if target is not None:cursor={'stage':'FORMAL_partialN_raw_extension','raw':j,'physical':physical,'source':row,'Nbase':nb,'Nmask':mask,'Sflat':flat};route(v,target,mask,flat%4,x,y)
   aux[j]=[]
   for t in range(2):
    v=[F(0),F(0)];aux[j].append(v)
    for i,x in Gamma[t].items():
     if i<48:cursor={'stage':'FORMAL_aux_Gamma','raw':j,'channel':t,'row':i};v[0]=add(v[0],mul(x,NX[i//4].get(4*17+i%4,F(0))))
    for pair,rr in enumerate(D[t]):
     for i,x in rr.items():cursor={'stage':'FORMAL_aux_D','raw':j,'channel':t,'pair':pair,'row':i};v[1]=add(v[1],mul(x,NX[i//4].get(4*(7+pair)+i%4,F(0))))
  del Sc,NX
  work.clear();fines={r['kernel105']:r['two_channels_Gamma_D_sum_DU'] for r in n['kernel_fine']}
  for rec in n['restrictions']:
   k=rec['kernel105'];cursor={'stage':'actualkernel_point_restriction','kernel':k};point={};work['kernel_point_partial']=point;kernels[k]=[]
   for j,q in rec['original_native_source_combination']:
    for row,x in points[j].items():cursor.update(raw=j,trace_row=row);put(point,row,mul(parse(q),x))
   assert not point
   for t in range(2):
    v=[F(0),F(0)];work['kernel_component_partial']=v
    for j,q in rec['original_native_source_combination']:
     for component in range(2):cursor={'stage':'ACTUAL_fixedkernel_restriction','kernel':k,'channel':t,'component':component,'raw':j};v[component]=add(v[component],mul(parse(q),aux[j][t][component]))
    z=add(*v)
    for x,bound in zip(v+[z],fines[k][t]):env(x,bound)
    kernels[k].append(z)
  for k in n['source_zero_kernel_ids']:
   assert all(int(bound[1])==0 for row in fines[k] for bound in row);kernels[k]=[F(0),F(0)]
  coord=dict(n['case_coordinates']);Q=[{i:parse(q) for i,q in row} for row in o['Q_two6']]
  for case,weights in o['corrected20_source_recipes']:
   cursor={'stage':'actual20_corrected_whole','case':case};piece=[];up=[];work['case_upper_partial']=up
   for t in range(2):
    cursor.update(channel=t);v=list(rawparts[case][t]);work['case_component_partial']=v
    for base,q in zip([563,609],weights):
     if F(q):
      for component in range(4):v[component]=add(v[component],neg(mul(parse(q),rawparts[base][t][component])))
    names=['OLD3','OLD2','FIRST','nativeOuter']
    for name,x in zip(names,v):env(x,E[name,t])
    native=F(0);work['case_native_partial']=native
    for k,q in coord[case]:cursor.update(native_kernel=k);native=add(native,mul(parse(q),kernels[k][t]));work['case_native_partial']=native
    env(native,E['native',t]);z=add(v[0],v[1]);work['whole_scalar_partial']=z;z=add(z,native);work['whole_scalar_partial']=z;z=add(z,neg(v[2]));work['whole_scalar_partial']=z;z=add(z,neg(v[3]));work['whole_scalar_partial']=z;env(z,E['upper',t]);up.append(z);piece.append({'OLD_Gamma':v[0],'OLD_D':v[1],'native_Gamma_plus_D':native,'FIRST_to_subtract':v[2],'native92_outer_to_subtract':v[3],'OLD4_outer_exact_coherent_ZERO_not_true_opaque_values0':True})
   full=[];work['case_six_partial']=full
   for i in range(6):
    z=F(0);work['Q_scalar_partial']=z
    for t in range(2):cursor={'stage':'fixed_Q_after_full20case','case':case,'coordinate6':i,'channel':t};z=add(z,mul(Q[t].get(i,F(0)),up[t]));work['Q_scalar_partial']=z
    env(z,E['six',i]);full.append(z)
   classes[case]=up;results.append({'case619':case,'two_source_component_pieces':piece,'actual_upper_two':up,'actual_full_six':full})
  ep=out/'actual20-full-components-upper2-six-and-shared-tables.json';dump(ep,{'all20':results,'raw21_partial_components_FORMAL_until605correction':rawparts,'all_raw21_FIRST_factor_values':FIRSTfactors,'FIRST_factor_order':'original a1 order, then pair10 order, then (other=b if var=a) followed by (other=a if var=b and a!=b); ordinary rr once','native_raw_auxiliaries_FORMAL_partialN_extensions':aux,'native_actual_kernel_restrictions':kernels,'all_evaluated_actualkernel_fullS0_zero_checked':True,'sourcezero_kernel_ids':n['source_zero_kernel_ids'],'coverage605':o['all605_coverage'],'only605_closed_after_peer_remaining12opaque_SECOND_open':True})
  # Fixed rank in a row-wise common integer lattice. Only comparisons use
  #501/502-bit TEMP products; stored minors/factors are <=334<354 bits.
  ds=[int(E['upper',t][0]) for t in range(2)];us=[int(E['upper',t][1]) for t in range(2)]
  for case in case20:
   v=[]
   for t in range(2):cursor={'stage':'integer_row_clear','case':case,'row':t};q=classes[case][t];assert ds[t]%q.denominator==0;x=binary(q.numerator,ck(ds[t]//q.denominator),'mul');assert abs(x)<=us[t];v.append(x)
   integer.append(v)
  selected=[];rows=[];minor=0
  for j,v in enumerate(integer):
   cursor={'stage':'rank_first_pivot_candidate','column20':j,'case':case20[j]}
   if any(v):selected=[j];rows=[next(i for i,x in enumerate(v) if x)];minor=v[rows[0]];break
  if selected:
   a=integer[selected[0]]
   for j,v in enumerate(integer):
    cursor={'stage':'rank_second_pivot_candidate','column20':j,'case':case20[j]}
    det=ck(binary(binary(a[0],v[1],'mul',True),-binary(a[1],v[0],'mul',True),'add',True))
    if det:selected.append(j);rows=[0,1];minor=det;break
  rank=len(selected);numerators=[]
  for j,v in enumerate(integer):
   cursor={'stage':'fixed_integer_rank_factor','column20':j,'rank':rank}
   if not rank:assert not any(v);factor_numerators=[];factors.append([])
   elif rank==1:
    factor_numerators=[v[rows[0]]];work['rank_numerators_current']=factor_numerators;a=integer[selected[0]]
    for i in range(2):cursor['row']=i;assert binary(minor,v[i],'mul',True)==binary(a[i],factor_numerators[0],'mul',True)
    factors.append([ck(F(factor_numerators[0],minor))])
   else:
    a,beta=[integer[i] for i in selected];factor_numerators=[ck(binary(binary(beta[1],v[0],'mul',True),-binary(beta[0],v[1],'mul',True),'add',True)),ck(binary(binary(a[0],v[1],'mul',True),-binary(a[1],v[0],'mul',True),'add',True))];work['rank_numerators_current']=factor_numerators
    for i in range(2):cursor['row']=i;assert binary(minor,v[i],'mul',True)==binary(binary(a[i],factor_numerators[0],'mul',True),binary(beta[i],factor_numerators[1],'mul',True),'add',True)
    factors.append([ck(F(x,minor)) for x in factor_numerators])
   numerators.append(factor_numerators)
  for j in range(20):
   if j in selected:continue
   cursor={'stage':'cleared_kernel_free_column','column20':j,'case':case20[j]}
   w={j:minor if rank else 1}
   for k,x in enumerate(numerators[j]):
    if x:w[selected[k]]=ck(-x)
   for row in range(2):
    cursor['row']=row
    z=0
    for k,x in w.items():cursor['source_column20']=k;z=binary(z,binary(integer[k][row],x,'mul',True),'add',True)
    assert z==0
   kernelcert.append(sorted(w.items()))
  rp=out/'actual20-fixed-integer-rank-all-factors-and-kernel.json';dump(rp,{'case619_order20':case20,'actual_integer_matrix2_by20':integer,'row_common_denominators':ds,'rank':rank,'selected_original20_columns':selected,'fixed_minor_rows':rows,'fixed_minor':minor,'all20_factor_numerators':numerators,'all20_rational_factors':factors,'cleared_kernel_vectors':kernelcert,'kernel_free_block':'minorI (I ifrank0)','all20_factor_and_kernel_identities_checked':True,'complete605_upper_rank_equals20rank_by584SOURCEzero_and554acceptedzero':True,'remaining12opaque_cases_open':True,'HF2_upper_only_without_fresh_parent_corollary':4-rank,'no_total_T2_H0_torsion_or_nilcube_inference':True})
  live=2*count(d)+count([N,G,nu,Gamma,D,leaves,C,T,first,L0,Le,aux,points,kernels,rawparts,FIRSTfactors,classes,results,integer,factors,kernelcert,numerators,old3,old2,Q])+count([cu,pr,g,pins])+4000;assert live<cap
  assert raw+sum(p.stat().st_size for p in out.iterdir())<rawcap
  report={'id':STEM,'all_passed':True,'actual_upper_rank20_and_complete605':rank,'selected_cases619':[case20[i] for i in selected],'all20_factor_and_kernel_identities':True,'remaining12opaque_not_evaluated':True,'source_committed_prior':354,'source_temporary_prior':1063,'actual_committed_bits':maximum,'actual_temporary_bits':tempmax,'actual_live_numeric_upper':live,'whole_numeric_cap':cap,'whole_RAW_cap':rawcap,'no_new_physical_gauge_IFT_jet_values':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'evidence_sha256':{ep.name:sha(ep),rp.name:sha(rp)},'seconds':time.monotonic()-start};dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','actual_upper_rank20_and_complete605','selected_cases619','actual_committed_bits','actual_live_numeric_upper','seconds']}))
 except BaseException as error:
  dump(out/'FAILURE_PREFIX.json',{'error':repr(error),'cursor':cursor,'pending_uncommitted_operation':pending,'current_work':work,'all_NEW_C0':C,'all_NEW_typed':T,'all_NEW_FIRST':first,'FORMAL_OLD_Gamma_table':old3,'FORMAL_OLD_D_table':old2,'all_FORMAL_aux':aux,'all_point_certificates':points,'all_actual_kernel_prefix':kernels,'all_raw_partial_prefix':rawparts,'all_FIRSTfactor_prefix':FIRSTfactors,'all_class_prefix':classes,'all_result_prefix':results,'integer_prefix':integer,'rank_state':{'selected':selected,'rows':rows,'minor':minor,'rank':rank,'numerators':numerators},'factors':factors,'kernel_prefix':kernelcert,'known_inputs_by_literal_pins':pins,'no_retry_or_capraise':True});raise
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output.resolve(),a.gate,a.source_peer,a.custody)

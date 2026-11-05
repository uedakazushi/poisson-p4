#!/usr/bin/env python3
"""Prepared ONE554 cache value bundle; requires fresh peer, custody and ROOT GO."""
import argparse, hashlib, json, resource, signal, time
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_ELL2_T2_UPPER554_VALUES_V2_20261105_v63'
GATE = 'reports/AMBIENT_B18_ELL2_T2_UPPER554_WHOLE_INPUT_GATE_20261105_v63/report.json'
GSHA = 'e250de02066f1df5abdbc512c4a268a9e3f83ca1722e0f34ba97986c01d9cdfa'
INPUT = 'ONE554-whole-positive-input-fine-assembly-and-native339-height.json'

def sha(p):
 h = hashlib.sha256()
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

def encoded(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):encoded(x) for k,x in v.items()}
 if isinstance(v,(list,tuple)):return list(map(encoded,v))
 return v

def dump(p,v):p.write_text(json.dumps(encoded(v),ensure_ascii=False,separators=(',',':'))+'\n')

def main(out,source_peer,custody):
 def timeout(signum,frame):raise TimeoutError('WHOLE60s_STOP')
 start=time.monotonic();signal.signal(signal.SIGALRM,timeout);signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False)
 pins={};raw=0;cursor={};pending={};maximum=0;temporary=0
 Sc={};NX={};oldjets={};firstcache={};pieces={};factor_values=[];upper=[];six=[];working={};operation_index=0
 def pin(p):
  nonlocal raw
  p=Path(p);p=p if p.is_absolute() else ROOT/p
  k=str(p.relative_to(ROOT));raw+=p.stat().st_size
  assert raw<15_767_720;pins[k]=sha(p);return p
 def read(p):return json.loads(pin(p).read_text())
 def ck(q):
  nonlocal maximum
  b=max(abs(q.numerator).bit_length(),q.denominator.bit_length());maximum=max(maximum,b)
  if b>339:cursor.update(offending_uncommitted_value=str(q),bits=b);raise RuntimeError('SOURCE339_STOP')
  return q
 def op(a,b,kind):
  nonlocal pending,temporary,operation_index
  operation_index+=1;cursor['binary_operation_index']=operation_index
  pending={'operation':kind,'operand_a':str(a),'operand_b':str(b),'cursor':dict(cursor)}
  ab=max(abs(a.numerator).bit_length(),a.denominator.bit_length());bb=max(abs(b.numerator).bit_length(),b.denominator.bit_length())
  temporary=max(temporary,ab+bb+(kind=='add'));assert temporary<=1018
  z=ck(a+b if kind=='add' else a*b);pending={};return z
 def add(a,b):return op(a,b,'add')
 def mul(a,b):return op(a,b,'mul')
 def neg(a):return ck(-a)
 def parse(q):return ck(F(q))
 def env(q,v):
  den,u=map(int,v);assert den%q.denominator==0 and abs(q.numerator)*(den//q.denominator)<=u
 def put(v,k,q):
  z=add(v.get(k,F(0)),q)
  if z:v[k]=z
  else:v.pop(k,None)
 def dot(a,b):
  z=F(0)
  working['dot_partial']=z
  for k,x in a.items():
   cursor['dot_row_index']=k
   if k in b:z=add(z,mul(x,b[k]));working['dot_partial']=z
  return z
 try:
  cu=read(custody);pr=read(source_peer);g=read(GATE)
  assert cu['all_passed'] and cu['safe_help_literal_PASS'] and pr['all_passed'] and g['all_passed']
  assert pins[GATE]==GSHA==cu['source_gate_sha256'] and GSHA in json.dumps(pr)
  assert sha(Path(__file__))==cu['worker_sha256'] and cu['all_actual_values_UNRUN_at_custody']
  for p,h in cu['worker_and_source_sha256'].items():assert sha(pin(p))==h
  for p,h in g['input_sha256'].items():assert sha(pin(p))==h
  ip=str((Path(GATE).parent/INPUT));d=read(ip);assert pins[ip]==g['evidence_sha256'][INPUT]
  lower=cu['case554_lower_zero_literal_record'];assert lower['case619']==554 and lower['actual_lower4_values']==['0']*4
  # Custody binds the original prefix record and its independent prefix14 peer.
  assert cu['accepted_lower_zero_peer_sha256']=='868ba1c858dcc03946180b2cbf952437cb8f0f7323bc1dc3e48363b43911bb3f'
  assert count([g,pr,cu])<2000 and 2*count(d)+9905+16000<55515
  assert d['selected_case619']==554 and d['global_value_committed_bits']==339 and d['global_value_temporary_bits']==1018
  o=d['OLD_FIRST_component'];n=d['NATIVE_outer_component'];p=n['full_primitive_input_J18_N_G_nu_GammaD']
  assert p['Gamma_two56']==o['Gamma_two'] and p['D_two_by10_by56']==o['D_two_by10_by56']
  assert n['all29_a1_kept_UNEVALUATED']==o['full_a1_29'] and n['plane_a1_pivot5_original644']==[[52,'-175/3'],[53,'-377/60']]
  dump(out/'phase.json',{'phase':'authorized_ONE554_cache_values','source_gate':GSHA,'source_peer':pins[str(source_peer)],'custody':pins[str(custody)],'accepted_lower4_zero':lower,'no_new_source_solve':True})
  Gamma=[{i:parse(q) for i,q in row} for row in p['Gamma_two56']]
  D=[[{i:parse(q) for i,q in enumerate(row) if q!='0'} for row in ch] for ch in p['D_two_by10_by56']]
  assert all(i<48 for ch in D for row in ch for i in row)
  needed=sorted({i//4 for row in Gamma for i in row if i<48}|{i//4 for ch in D for row in ch for i in row});assert len(needed)==7
  basis=[tuple(v) for v in p['basis18']];bi={m:i for i,m in enumerate(basis)};assert len(basis)==18
  known=p['N_known_positions'];assert known==[0,1,2,3,4,5,6,15,16]
  N=[[[[parse(q) for q in cell] for cell in row] for row in src] for src in p['N_known12_by12_by9_by4'][:4]]
  G=[[parse(q) for q in row] for row in p['G8_C72']];nu=[[parse(q) for q in row] for row in p['nu_two_J18']]
  assert all(tuple(a+b for a,b in zip(basis[i],basis[j])) not in bi for x in nu for y in nu for i,a in enumerate(x) if a for j,b in enumerate(y) if b)
  def target(a,b):return bi.get(tuple(x+y for x,y in zip(basis[a],basis[b])))
  def route(outv,base,mask,im,x,y):
   # Every unsupported base route is rejected BEFORE scalar products.
   collision=mask&im;om=mask^im
   if collision==3:return
   if not collision:put(outv,4*base+om,mul(x,y));return
   for vb,z in enumerate(nu[collision-1]):
    cursor['monic_base_index']=vb
    if not z:continue
    dest=target(base,vb)
    if dest is not None:put(outv,4*dest+om,mul(mul(x,y),z))
  for rec_index,rec in enumerate(n['parameterized_four_native_source_records']):
   cursor={'stage':'weighted_intrinsic_S','source':rec['original768'],'record_index':rec_index};im=rec['input_trace_mask'];w=parse(rec['constant_source_coefficient'])
   for entry_index,(kind,row,gen,sign) in enumerate(rec['intrinsic_source_entries']):
    cursor.update(intrinsic_entry_index=entry_index,G_generator=gen)
    assert kind=='G8' and row<4 and sign in [-1,1];v=Sc.setdefault(row,{})
    for flat,x in enumerate(G[gen]):
     cursor['G_flat_index']=flat
     if x:route(v,flat//4,flat%4,im,x,neg(w) if sign==-1 else w)
  assert all(not any(k<4 for k in row) for row in Sc.values())
  expectedS=[{i:(den,u) for i,den,u in row} for row in n['positive_weighted_intrinsic_S18']]
  for row,v in Sc.items():
   for flat,x in v.items():env(x,expectedS[flat//4].get(4*row+flat%4,(1,0)))
  assert not expectedS[0]
  for physical in needed:
   v=NX.setdefault(physical,{})
   for source,s in Sc.items():
    for kp,nb in enumerate(known):
     for mask,x in enumerate(N[source][physical][kp]):
      if not x:continue
      for flat,y in s.items():
       dest=target(nb,flat//4)
       if dest is None:continue
       cursor={'stage':'known_partial_N_transpose_S','physical':physical,'source':source,'Nbase':nb,'Nmask':mask,'Sflat':flat}
       route(v,dest,mask,flat%4,x,y)
  for t in range(2):
   ng=F(0);nd=F(0);working.update(native_Gamma=ng,native_D=nd)
   for i,a in Gamma[t].items():
    if i<48:
     cursor={'stage':'native_Gamma_contraction','channel':t,'row':i};ng=add(ng,mul(a,NX[i//4].get(4*17+i%4,F(0))));working['native_Gamma']=ng
   for pair,row in enumerate(D[t]):
    for i,a in row.items():
     cursor={'stage':'native_D_contraction','channel':t,'pair':pair,'row':i};nd=add(nd,mul(a,NX[i//4].get(4*(7+pair)+i%4,F(0))));working['native_D']=nd
   env(ng,d['positive_fine_native_Gamma_two'][t]);env(nd,d['positive_fine_native_D_two'][t])
   pieces.setdefault(t,{}).update(native_Gamma=ng,native_D=nd)
  # Native contraction is now back under the fine bounds; all following
  # assembly uses the fine214 proof, never multiplies a coarse339 envelope.
  zero9=(0,)*9;pairs=[(0,7),(0,8),(1,7),(1,8),(4,7),(4,8),(5,7),(5,8),(7,7),(7,8)]
  b21=[zero9]+[tuple(int(j==i) for j in range(9)) for i in range(9)]+[tuple(int(j==a)+int(j==b) for j in range(9)) for a,b in pairs]+[(0,0,0,0,0,0,0,2,1)]
  bidx={m:i for i,m in enumerate(b21)};assert len(bidx)==21
  physical={k:{i:parse(q) for i,q in enumerate(v) if q!='0'} for k,v in o['physical_B21_functions7']}
  def polyadd(v,w,c):
   working['current_OLD_polynomial']=v
   for k,x in w.items():cursor['OLD_polynomial_add_index']=k;put(v,k,mul(c,x))
  def polymul(v,w):
   ans={};working['current_OLD_polynomial_product']=ans
   for i,x in v.items():
    for j,y in w.items():
     cursor['OLD_polynomial_product_indices']=[i,j];k=bidx.get(tuple(a+b for a,b in zip(b21[i],b21[j])))
     if k is not None:put(ans,k,mul(x,y))
   return ans
  for node,desc in o['OLD_full_typed_program11']:
   cursor={'stage':'actual_OLD_B21','typednode':node}
   if desc[0]=='actual_C0_scalar_function':v=dict(physical[desc[1]])
   elif desc[0]=='sum':
    v={}
    for dep,num,den in desc[1]:polyadd(v,oldjets[dep],parse(F(num,den)))
   else:
    assert desc[0]=='product';v={0:parse(F(desc[1],desc[2]))}
    for dep in desc[3]:v=polymul(v,oldjets[dep])
   oldjets[node]=v
  for t in range(2):
   og=F(0);od=F(0);working.update(OLD_Gamma=og,OLD_D=od)
   for source,w,roots in o['OLD_Gamma_source_columns']:
    for row,tag in roots:
     if row not in Gamma[t]:continue
     v=oldjets[tag[1]] if tag[0]=='typed' else physical.get(tag[1],{})
     cursor={'stage':'OLD_Gamma_contraction','channel':t,'source':source,'row':row};og=add(og,mul(mul(parse(w),Gamma[t][row]),v.get(20,F(0))));working['OLD_Gamma']=og
   for source,w,roots in o['OLD_D_source_columns']:
    for row,node in roots:
     for pair in range(10):
      if row in D[t][pair]:
       cursor={'stage':'OLD_D_contraction','channel':t,'source':source,'row':row,'node':node,'pair':pair};od=add(od,mul(mul(parse(w),D[t][pair][row]),oldjets[node].get(10+pair,F(0))));working['OLD_D']=od
   env(og,o['positive_OLD_Gamma_cubic_two'][t]);env(od,o['positive_OLD_D_second_two'][t]);pieces[t].update(OLD_Gamma=og,OLD_D=od)
  L0={j:{i:parse(x) for i,x in v} for j,v in o['native_L0_selected_columns']};Le={j:{i:parse(x) for i,x in v} for j,v in o['native_Le_selected_columns']}
  for rec in o['Xp_FIRST_column_recipes28']:
   key=(rec['pivot27'],rec['variable9']);r=rec['first_column_recipe'];v={};firstcache[key]=v;cursor={'stage':'individual_Xp_FIRST','pivot_variable':key}
   if r['kind']=='literal56':v={i:parse(q) for i,q in r['column']}
   else:
    assert r['kind']=='native_ell_SUM'
    for label,mat in [('Se',L0),('S0',Le)]:
     for j,q in r[label]:
      assert j in mat
      for i,x in mat[j].items():cursor.update(SUM_branch=label,S_input=j,output_row=i);put(v,i,mul(parse(q),x))
   firstcache[key]=v
  for pvt,var,bound in o['positive_Xp_FIRST_columns28']:
   v=firstcache[pvt,var];bb={i:(den,u) for i,den,u in bound}
   for i,x in v.items():env(x,bb.get(i,(1,0)))
  for t in range(2):
   carry=F(0);fv=[];factor_values.append(fv);working['FIRST_carry']=carry
   for j,(pvt,pair,var,a) in enumerate(o['known_FIRSTcarry_factor_recipes57']):
    cursor={'stage':'FIRST57_contraction','channel':t,'factor':j};z=mul(parse(a),dot(D[t][pair],firstcache[pvt,var]));env(z,o['positive_FIRST57_factor_bounds_two_by57'][t][j]);fv.append(z);carry=add(carry,z);working['FIRST_carry']=carry
   env(carry,o['positive_FIRSTcarry_total_two'][t])
   xp=n['actual644_GammaXp_rr_rell_twochannels_UNPAIRED'][t]
   cursor={'stage':'actual_outer644_pairing','channel':t};outer_terms=[];working['outer_terms']=outer_terms
   for direction,coef in enumerate(['-377/60','-175/3']):cursor['rr_or_rell']=direction;outer_terms.append(mul(parse(xp[direction]),parse(coef)))
   for z,b in zip(outer_terms,n['positive_outer644_two_by2'][t]):env(z,b)
   outer=add(*outer_terms);env(outer,d['positive_outer_total_two'][t]);pieces[t].update(FIRST_to_subtract=carry,outer644_terms=outer_terms,outer644_to_subtract=outer)
   z=F(0);working['upper_assembly_partial']=z;cursor={'stage':'whole_upper_assembly','channel':t}
   for k in ['OLD_Gamma','OLD_D','native_Gamma','native_D']:cursor['source_component']=k;z=add(z,pieces[t][k]);working['upper_assembly_partial']=z
   z=add(add(z,neg(carry)),neg(outer));env(z,d['positive_WHOLE_upper_two'][t]);upper.append(z)
  Q=[{i:parse(q) for i,q in row} for row in d['fixed_Q_embedding_two6']]
  for i in range(6):
   z=F(0);working['Q_partial']=z
   for t in range(2):
    cursor={'stage':'fixed_Q_embedding','coordinate6':i,'channel':t};z=add(z,mul(Q[t].get(i,F(0)),upper[t]));working['Q_partial']=z
   env(z,d['positive_WHOLE_full_six'][i]);six.append(z)
  state=[N,G,nu,Gamma,D,Sc,NX,physical,oldjets,L0,Le,firstcache,pieces,factor_values,upper,six,Q,working]
  actual_live=2*count(d)+count(state)+count([g,cu,pr,pins])+4000;assert actual_live<=55515
  ep=out/'actual-ONE554-complete-source-pieces-and-quadratic-class.json'
  dump(ep,{'case619':554,'original768_source_all6':d['literal_source_all6'],'all29_a1':o['full_a1_29'],'accepted_lower4_zero':lower,'source_components_two_channels':pieces,'all114_FIRSTfactor_values':factor_values,'actual_upper_two':upper,'actual_full_six':six,'nonzero_class':any(six),'known_partial_N_extension_independence_source_not_actual_unknownNzero':True,'all_fine_component_and_whole_envelopes_checked':True,'no_new_physical_gauge_implicit_or_jet_coefficients':True,'only_ONE554_no_all605_or_HF2_equality_total_length_H0_claim':True})
  assert raw+sum(p.stat().st_size for p in out.iterdir())<15_767_720
  report={'id':STEM,'all_passed':True,'case619':554,'actual_upper_two':list(map(str,upper)),'actual_full_six':list(map(str,six)),'nonzero_class':any(six),'all_actual_response_values_before_this_job_UNRUN':True,'no_new_physical_gauge_implicit_or_jet_coefficients':True,'source_committed_bits':339,'source_temporary_bits':1018,'actual_committed_bits':maximum,'actual_binary_temporary_bits':temporary,'actual_live_numeric_upper':actual_live,'whole_numeric_cap':55515,'whole_RAW_cap':15767720,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start}
  dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','actual_upper_two','nonzero_class','actual_committed_bits','actual_live_numeric_upper','seconds']}))
 except BaseException as error:
  dump(out/'FAILURE_PREFIX.json',{'error':repr(error),'cursor':cursor,'pending_uncommitted_binary_operation':pending,'current_local_working_accumulators':working,'all_NEW_Sc':Sc,'all_NEW_NSc':NX,'all_NEW_OLDjets':oldjets,'all_NEW_FIRSTcache':firstcache,'completed_pieces':pieces,'completed_FIRSTfactors':factor_values,'completed_upper':upper,'completed_six':six,'immutable_known_inputs_by_pins':pins,'no_retry_or_cap_raise':True})
  raise

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True)
 a=p.parse_args();main(a.output.resolve(),a.source_peer,a.custody)

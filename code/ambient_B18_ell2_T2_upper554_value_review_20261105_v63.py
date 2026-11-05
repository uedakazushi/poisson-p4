#!/usr/bin/env python3
"""Fixed ONE554 signed postcheck by covector-L and dense OLD, no new cases."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER554_VALUES_V2_20261105_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple,set)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0;bits=0;cursor={};R={};S=[];jets={};first={};own=[]
 def read(path):
  nonlocal raw
  p=ROOT/path;pins[path]=sha(p);raw+=2*p.stat().st_size;assert raw<32*1024*1024;return json.loads(p.read_text())
 def dump(name,v):(out/name).write_text(json.dumps(v,separators=(',',':'))+'\n')
 def encode(v):
  if isinstance(v,F):return str(v)
  if isinstance(v,dict):return {str(k):encode(x) for k,x in v.items()}
  if isinstance(v,(list,tuple)):return list(map(encode,v))
  return v
 def ck(q):
  nonlocal bits
  bits=max(bits,abs(q.numerator).bit_length(),q.denominator.bit_length());assert bits<=339;return q
 def add(a,b):return ck(a+b)
 def mul(a,b):return ck(a*b)
 def env(q,z):
  d,u=map(int,z);assert d%q.denominator==0 and abs(q.numerator)*(d//q.denominator)<=u
 try:
  gp='reports/AMBIENT_B18_ELL2_T2_UPPER554_WHOLE_INPUT_GATE_20261105_v63/';g=read(gp+'report.json');w=read(gp+'ONE554-whole-positive-input-fine-assembly-and-native339-height.json')
  peer=read('reports/AMBIENT_B18_ELL2_T2_UPPER554_WHOLE_GATE_20261105_v63_REVIEW/report.json');assert peer['all_passed'] and pins['reports/AMBIENT_B18_ELL2_T2_UPPER554_WHOLE_GATE_20261105_v63_REVIEW/report.json']=='e2c34f0b5c399d60781d5e01c2df79208449905897f7dbc490374456215a1e48'
  bp='reports/AMBIENT_B18_ELL2_T2_UPPER554_VALUES_V2_20261105_v63/';pr=read(bp+'report.json');assert pr['all_passed'] and pins[bp+'report.json']=='b339566cba6c4b3e297d05561a80e033c1d463346d156cc79cf8da7887191e2f'
  pn='actual-ONE554-complete-source-pieces-and-quadratic-class.json';saved=read(bp+pn);assert pins[bp+pn]==pr['evidence_sha256'][pn]
  n=w['NATIVE_outer_component'];o=w['OLD_FIRST_component'];p=n['full_primitive_input_J18_N_G_nu_GammaD'];assert saved['case619']==554 and saved['original768_source_all6']==w['literal_source_all6'] and saved['all29_a1']==o['full_a1_29']
  assert w['global_value_committed_bits']==339 and w['global_value_temporary_bits']==1018
  before=count([g,w,peer,pr,saved,pins])+2*13820+20000;assert before<150000 and raw+2_000_000<32*1024*1024
  dump('phase-before-signed-fixed-replay.json',{'onlyONE554_fixed_inputs':True,'source_prior339_temporary1018':True,'own_source_before_arithmetic_upper':before,'native_method':'covectorL68 thenweightedS/fiveGamma+tenDordinaryterms, not direct NSc','OLD_method':'dense21 signed polynomial with reverse sums','no_new_cases_or_callback_functions':True})
  N=[[[list(map(lambda q:ck(F(q)),cell)) for cell in row] for row in src] for src in p['N_known12_by12_by9_by4']]
  G=[[ck(F(q)) for q in row] for row in p['G8_C72']];nu=[[ck(F(q)) for q in row] for row in p['nu_two_J18']]
  Gamma=[{i:ck(F(q)) for i,q in row} for row in p['Gamma_two56']];D=[[list(map(lambda q:ck(F(q)),row)) for row in ch] for ch in p['D_two_by10_by56']]
  basis=[tuple(v) for v in p['basis18']];bi={m:i for i,m in enumerate(basis)};known=p['N_known_positions']
  table=[[bi.get(tuple(x+y for x,y in zip(a,b))) for b in basis] for a in basis]
  assert all(table[i][j] is None for a in nu for b in nu for i,x in enumerate(a) if x for j,y in enumerate(b) if y)
  bounds={tuple(k):{i:(d,u) for i,d,u in row} for k,row in n['reused68_positive_covectorL_rows']}
  # Build signed operator rows first, unrelated to producer's NSc route.
  for key in sorted(bounds,reverse=True):
   kind,t,pair,target=key;cov=Gamma[t] if kind==0 else {i:q for i,q in enumerate(D[t][pair]) if q};row=[F(0)]*48
   for src in range(11,-1,-1):
    for im in range(3,-1,-1):
     z=F(0);cursor={'stage':'signed_covectorL','key':key,'source':src,'inputmask':im}
     for physical in range(11,-1,-1):
      for kp in range(8,-1,-1):
       nb=known[kp]
       for mask in range(3,-1,-1):
        x=N[src][physical][kp][mask];a=cov.get(4*physical+(mask^im),F(0))
        if not x or not a:continue
        collision=mask&im
        if collision==0:
         if nb==target:z=add(z,mul(a,x))
        elif collision in [1,2]:
         for vb,y in enumerate(nu[collision-1]):
          if y and table[nb][vb]==target:z=add(z,mul(mul(a,x),y))
        else:assert collision==3
     row[4*src+im]=z;env(z,bounds[key].get(4*src+im,(1,0)))
   R[key]=row
  assert len(R)==68
  S=[[F(0)]*48 for _ in basis]
  for rec in reversed(n['parameterized_four_native_source_records']):
   weight=F(rec['constant_source_coefficient']);im=rec['input_trace_mask']
   for kind,row,gen,sign in reversed(rec['intrinsic_source_entries']):
    for flat in range(71,-1,-1):
     x=G[gen][flat]
     if not x:continue
     gb,mask=divmod(flat,4);collision=mask&im;col=4*row+(mask^im);q=mul(weight*sign,x)
     if not collision:S[gb][col]=add(S[gb][col],q)
     elif collision in [1,2]:
      for vb,y in enumerate(nu[collision-1]):
       dest=table[gb][vb]
       if y and dest is not None:S[dest][col]=add(S[dest][col],mul(q,y))
  assert S[0]==[F(0)]*48
  for b,eb in enumerate(n['positive_weighted_intrinsic_S18']):
   bound={i:(d,u) for i,d,u in eb}
   for i,q in enumerate(S[b]):env(q,bound.get(i,(1,0)))
  def dot(a,b):
   z=F(0)
   for x,y in reversed(list(zip(a,b))):
    if x and y:z=add(z,mul(x,y))
   return z
  vf={0:1,1:2,4:3,5:4,7:5,8:6};pairs=[(0,7),(0,8),(1,7),(1,8),(4,7),(4,8),(5,7),(5,8),(7,7),(7,8)];cubic=[(0,17),(5,16),(6,15),(15,6),(16,5)];nativeparts=[]
  for t in range(2):
   cg=[];ds=[]
   for j,(a,b) in enumerate(cubic):q=dot(R[0,t,-1,a],S[b]);env(q,n['positive_Gamma_weighted_cubic_two_by5'][t][j]);cg.append(q)
   for pair,(a,b) in enumerate(pairs):
    q=[dot(R[1,t,pair,0],S[7+pair])]
    for x,y in ([(a,b)] if a==b else [(a,b),(b,a)]):q.append(dot(R[1,t,pair,vf[x]],S[vf[y]]))
    for x,bound in zip(q,n['positive_D_weighted_second_two_by10_unique_first_decomp'][t][pair]):env(x,bound)
    ds.append(q)
   ng=F(0);nd=F(0)
   for q in reversed(cg):ng=add(ng,q)
   for q in reversed([x for row in ds for x in row]):nd=add(nd,q)
   own.append({'native_Gamma':ng,'native_D':nd});nativeparts.append({'Gamma_five':cg,'D_ten_ordinary_pieces':ds})
  # Dense B21 calculation preserves all physical first-cross terms.
  zero=(0,)*9;b21=[zero]+[tuple(int(i==j) for i in range(9)) for j in range(9)]+[tuple(int(i==a)+int(i==b) for i in range(9)) for a,b in pairs]+[(0,0,0,0,0,0,0,2,1)];bidx={m:i for i,m in enumerate(b21)}
  ph={i:[ck(F(q)) for q in row] for i,row in o['physical_B21_functions7']}
  def pm(a,b):
   z=[F(0)]*21
   for i in range(20,-1,-1):
    if a[i]:
     for j in range(20,-1,-1):
      k=bidx.get(tuple(x+y for x,y in zip(b21[i],b21[j])))
      if b[j] and k is not None:z[k]=add(z[k],mul(a[i],b[j]))
   return z
  for i,node in o['OLD_full_typed_program11']:
   cursor={'stage':'denseOLD','node':i}
   if node[0]=='actual_C0_scalar_function':v=ph[node[1]][:]
   elif node[0]=='sum':
    v=[F(0)]*21
    for j,a,b in reversed(node[1]):
     q=ck(F(a,b));v=[add(x,mul(q,y)) for x,y in zip(v,jets[j])]
   else:
    assert node[0]=='product';v=[F(0)]*21;v[0]=ck(F(node[1],node[2]))
    for j in reversed(node[3]):v=pm(v,jets[j])
   jets[i]=v
  for t in range(2):
   a=F(0);b=F(0)
   for source,c,roots in reversed(o['OLD_Gamma_source_columns']):
    for row,(kind,i) in reversed(roots):
     v=jets[i] if kind=='typed' else ph.get(i,[F(0)]*21);a=add(a,mul(mul(F(c),Gamma[t].get(row,F(0))),v[20]))
   for source,c,roots in reversed(o['OLD_D_source_columns']):
    for row,i in reversed(roots):
     for pair in range(9,-1,-1):b=add(b,mul(mul(F(c),D[t][pair][row]),jets[i][10+pair]))
   own[t].update(OLD_Gamma=a,OLD_D=b)
  L0={j:{i:F(q) for i,q in row} for j,row in o['native_L0_selected_columns']};Le={j:{i:F(q) for i,q in row} for j,row in o['native_Le_selected_columns']}
  for rec in o['Xp_FIRST_column_recipes28']:
   key=(rec['pivot27'],rec['variable9']);r=rec['first_column_recipe'];v=[F(0)]*56
   if r['kind']=='literal56':
    for i,q in r['column']:v[i]=ck(F(q))
   else:
    assert r['kind']=='native_ell_SUM'
    for label,mat in [('S0',Le),('Se',L0)]:
     for j,q in reversed(r[label]):
      for i,a in mat[j].items():v[i]=add(v[i],mul(F(q),a))
   first[key]=v
  allf=[]
  r92=read('reports/AMBIENT_B18_NATIVE_PIVOT23_SECOND_VALUES_20261105_v63/actual92-Gamma-nativepivot23-plane2-and322-full-product-pieces.json');idx=r92['original768_native23_sources'].index(644);assert idx==1
  for t in range(2):
   vals=[]
   for k,(pvt,pair,var,c) in enumerate(o['known_FIRSTcarry_factor_recipes57']):cursor={'stage':'FIRST114','channel':t,'factor':k};q=mul(F(c),dot(D[t][pair],first[pvt,var]));env(q,o['positive_FIRST57_factor_bounds_two_by57'][t][k]);vals.append(q)
   carry=F(0)
   for q in reversed(vals):carry=add(carry,q)
   assert list(map(F,saved['all114_FIRSTfactor_values'][t]))==vals;allf.append(vals)
   outer=[mul(F(r92['Gamma_Xp2_native23_twochannels_rr_rell'][idx][t][0]),F('-377/60')),mul(F(r92['Gamma_Xp2_native23_twochannels_rr_rell'][idx][t][1]),F('-175/3'))]
   own[t].update(FIRST_to_subtract=carry,outer644_terms=outer,outer644_to_subtract=add(outer[1],outer[0]))
   for key,value in own[t].items():
    old=saved['source_components_two_channels'][str(t)][key];assert value==(list(map(F,old)) if isinstance(value,list) else F(old)),(t,key)
  upper=[]
  for t in range(2):
   z=add(own[t]['native_Gamma'],own[t]['OLD_Gamma']);z=add(z,add(own[t]['native_D'],own[t]['OLD_D']));z=add(z,-add(own[t]['outer644_to_subtract'],own[t]['FIRST_to_subtract']));upper.append(z)
  Q=[dict(v) for v in w['fixed_Q_embedding_two6']];six=[add(mul(F(Q[1].get(i,'0')),upper[1]),mul(F(Q[0].get(i,'0')),upper[0])) for i in range(6)]
  assert upper==list(map(F,saved['actual_upper_two']))==[F(0)]*2 and six==list(map(F,saved['actual_full_six']))==[F(0)]*6
  for path,h in pr['input_sha256'].items():assert sha(ROOT/path)==h;pins[path]=h
  owned=count([g,w,peer,pr,saved,N,G,nu,Gamma,D,basis,bi,known,table,R,S,jets,ph,first,L0,Le,own,nativeparts,allf,r92,upper,six,Q,pins])+15000;assert owned<150000
  ep='independent-covectorL-native-denseOLD-FIRST114-and644-complete554-check.json';ev={'all68_signed_covectorL_rows_and864_intrinsic_S_coefficients_bounded':True,'native_method_independent_of_direct_partialNSc_worker':True,'all_native_fiveGamma_andtenDordinary_piece_values':encode(nativeparts),'all_OLD_native_FIRST_andouter_saved_components_verified':True,'all114_FIRST_exact':True,'outer644_fullL2S0_actual92_pairings_verified':True,'all_source_components_two_channels':encode(own),'actual_upper_two':encode(upper),'actual_full_six':encode(six),'one554_only_no_all605_HF_orlength':True,'no_new_physical_gauge_higher_functions_or_cases':True,'actual_QQ_bits':bits,'whole_owned_numeric_upper':owned}
  dump(ep,ev);assert raw+sum(q.stat().st_size for q in out.iterdir())<32*1024*1024
  report={'id':STEM,'all_passed':True,'independent_signed_FIXED_ONE554_postcheck':True,'actual_upper_two':['0','0'],'actual_full_six':['0']*6,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep),'phase-before-signed-fixed-replay.json':sha(out/'phase-before-signed-fixed-replay.json')},'source_prior339_temporary1018':True,'actual_bits':bits,'owned_numeric_upper':owned,'whole_RAW_upper':raw+sum(q.stat().st_size for q in out.iterdir()),'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
  dump('report.json',report);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'bits':bits,'owned':owned,'seconds':report['seconds']}))
 except BaseException as e:
  dump('FAILURE_PREFIX.json',encode({'error':repr(e),'cursor':cursor,'R':R,'S':S,'OLDjets':jets,'first':first,'components':own,'input_sha256':pins}));raise
if __name__=='__main__':
 pa=argparse.ArgumentParser();pa.add_argument('--output',required=True,type=Path);a=pa.parse_args();main(a.output)

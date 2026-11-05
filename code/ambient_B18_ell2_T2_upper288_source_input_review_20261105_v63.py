#!/usr/bin/env python3
"""Independent ONE288 source and positive-zero proof; no signed response."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER288_SOURCE_INPUT_20261105_v63_REVIEW'
SP='reports/AMBIENT_B18_ELL2_T2_UPPER288_SOURCE_20261105_v63/'
IP='reports/AMBIENT_B18_ELL2_T2_UPPER288_INPUT_GATE_20261105_v63/'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple,set)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def fields(p,names):
 txt=p.read_text();dec=json.JSONDecoder();i=txt.index('{')+1;d={}
 while True:
  while txt[i].isspace() or txt[i]==',':i+=1
  if txt[i]=='}':break
  k,i=dec.raw_decode(txt,i)
  while txt[i].isspace() or txt[i]==':':i+=1
  v,i=dec.raw_decode(txt,i)
  if k in names:d[k]=v
  del v
 assert set(d)==set(names);return d
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0;maximum=0
 def read(p,names=None):
  nonlocal raw
  q=ROOT/p;pins[p]=sha(q);raw+=2*q.stat().st_size
  return json.loads(q.read_text()) if names is None else fields(q,names)
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 sr=read(SP+'report.json');assert sr['all_passed'] and pins[SP+'report.json']=='91fff8d1815310c4b5e197382b99667b6c8e303e62e5e7bac5b69c7a3f35ea0c'
 se='ONE288-opaque-free-source-and-eight-known-first-pivot-carries.json';s=read(SP+se);assert pins[SP+se]==sr['evidence_sha256'][se]
 ir=read(IP+'report.json');assert ir['all_passed'] and pins[IP+'report.json']=='d9286144743b0ef44192357ed8339c787bcc76c74a0fc8852774854fbf3bd7bd'
 ie='ONE288-three-physical-functions-first-carry-positive-input.json';d=read(IP+ie);assert pins[IP+ie]==ir['evidence_sha256'][ie]
 op='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_SOURCE_20261105_v63/combined-OLD2-D-and-four-pivot-source-factor-inventory.json'
 orig=read(op,['original619_c0_source_combinations','kernel617_source_recipes','typed_program','D_old_root_columns','D_two_by10_by56'])
 fr='reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json';a=read(fr,['all_zero_bottom_and_Y0_source_combinations'])['all_zero_bottom_and_Y0_source_combinations'][288]
 assert a['original768_combination']==orig['original619_c0_source_combinations'][288]==s['original768_source']==[[352,'1']]
 assert a['actual_a1_27_times9']==s['actual_full27times9_a1'] and len(a['actual_a1_27_times9'])==8 and all(i%9==1 for i,q in a['actual_a1_27_times9'])
 recipe=next(v for v in orig['kernel617_source_recipes'] if v['free_case619']==288);assert recipe==s['same617_source_recipe'] and recipe['subtract_actual_basis563609']==['0','0']
 hp='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_INPUT_GATE_20261105_v63/617-direct-old-known-first-formal-Jacobian-positive-input.json';h=read(hp,['typed_source_program','literal_C0_program','pruned_original_old_columns_roots','Gamma_two'])
 ops=dict(h['typed_source_program']);u=dict(orig['typed_program']);assert all(ops[i]==n for i,n in u.items() if i in ops);ops.update(u)
 assert dict(orig['D_old_root_columns'])[352]==s['D_old352_source_roots'] and dict(h['pruned_original_old_columns_roots'])[352]==s['Gamma_old352_source_roots']
 assert h['Gamma_two']==s['Gamma_two'] and orig['D_two_by10_by56']==s['D_two_by10_by56']
 roots={i for row,i in s['D_old352_source_roots']}|{i for row,(kind,i) in s['Gamma_old352_source_roots'] if kind=='typed'};todo=list(roots);keep=set();leaves=set()
 while todo:
  i=todo.pop()
  if i in keep:continue
  keep.add(i);n=ops[i]
  if n[0]=='actual_C0_scalar_function':leaves.add(n[1])
  elif n[0]=='sum':todo.extend(j for j,a,b in n[1])
  else:assert n[0]=='product';todo.extend(n[3])
 assert [[i,ops[i]] for i in sorted(keep)]==s['full_typed_program']==d['full_literal_typed_program5']
 assert leaves=={20117,27948,40895} and s['literal_C0_program']==d['literal_C0_zero']==[[0,['c',0,1]]]
 # Exact original scalar graph is delta=phi20117-phi40895, plus zB1. No
 # multiproduct/unknown response is discarded in obtaining this expression.
 assert ops[101]==['actual_C0_scalar_function',27948] and ops[186]==['actual_C0_scalar_function',20117] and ops[187]==['actual_C0_scalar_function',40895]
 assert ops[237]==['product',-1,1,[187]] and ops[238]==['sum',[[186,1,1],[237,1,1]]]
 ccpath='reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_INPUT_GATE_20261105_v63/pruned-known-Cfull224-input-and-positive-source-envelopes.json'
 cc=read(ccpath);assert cc['Gamma_two']==s['Gamma_two']
 assert [(v['pivot27'],v['a1_z1']) for v in s['known_first_carry_pivots8']]==[(i//9,q) for i,q in a['actual_a1_27_times9']]
 for rec in s['known_first_carry_pivots8']:
  j=rec['pivot27'];assert rec['original768']==cc['pivot27_source_indices'][j] and rec['Xp_r_literal56']==cc['pivot27_r_columns'][j] and rec['Xp_ell_recipe']==cc['pivot27_ell_recipes'][j]
 for key in ['native_L0','native_Le','native_S0','native_Se']:assert s[key]==cc[key]
 # Same B21 and fullPsi current frame already independently source-proved.
 bp='reports/AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63/ONE508-source-closure-B21-and-exact-carry-zeros.json';bsource=read(bp,['ordinary_basis_B21','raw33_to_B21_selector21']);assert s['ordinary_basis_B21']==bsource['ordinary_basis_B21'] and s['raw33_to_B21_selector21']==bsource['raw33_to_B21_selector21']
 assert read('reports/AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63_REVIEW/report.json')['all_passed']
 j1=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_FORMAL_20261105_v63/617-formal-direct-old-third-and-second-symbol-sensitivities.json',['all617_formal_directOLD_J0_Jr_Jell'])['all617_formal_directOLD_J0_Jr_Jell'];j=next(v for v in j1 if v['free_case619']==288);assert j==s['original_case288_J1_certificate'] and {i for ch in j['two_formal_J0_Jr_Jell_symbol_rows'] for i,q in ch}=={20117,40895}
 o2=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63/actual617-formal-combined-OLD2-second-symbol-coefficients.json',['all617_combined_OLD2_independent_symbol_coefficients'])['all617_combined_OLD2_independent_symbol_coefficients'];o=next(ch for i,ch in o2 if i==288);assert o==s['original_case288_OLD2_certificate'] and {i for ch in o for pair,i,q in ch}==leaves
 pv=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',['all_requested164_root_values'])['all_requested164_root_values'];p3={i:dict(v) for i,v in pv if i in leaves};del pv
 fv=read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json',['first640'])['first640']['all79_atom_values_full10'];f3={i:dict(v) for i,v in fv if i in leaves};del fv
 for i,v in d['physical_B21_functions3']:
  assert list(map(F,v))==[F(p3[i].get(k,'0')) for k in s['raw33_to_B21_selector21']]
  assert list(map(F,v[:10]))==[F(f3[i].get(k,'0')) for k in range(10)]
 Q=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_FULL_COORDINATE_SOURCE_20261105_v63/known-Psi-point-correction-and-full-b1-source-annihilation.json',['Q_two6'])['Q_two6'];assert Q==d['fixed_Q_embedding_two6']
 owner0=count([r for r in [sr,s,ir,d,orig,a,h,ops,u,cc,j1,o2,p3,f3,Q,pins]])+15000;assert owner0<150000 and raw+1_000_000<32*1024*1024
 dump('phase-before-positive.json',{'phase':'original288_source_z1carry_frame_cache_and_owner_PASS','source90_temporary271':True,'owner_upper':owner0,'no_signed_response':True})
 def ck(v):
  nonlocal maximum
  assert v>=0;maximum=max(maximum,v.numerator.bit_length(),v.denominator.bit_length());assert maximum<=90;return v
 def plus(a,b):return ck(a+b)
 def mul(a,b):return ck(a*b)
 def av(v):return {i:ck(abs(F(x))) for i,x in v if F(x)}
 def pos_mv(cols,v):
  z=[F(0)]*56
  for j,a in v.items():
   for i,b in cols[j].items():z[i]=plus(z[i],mul(a,b))
  return z
 def dot(row,v):
  q=F(0)
  for a,b in zip(row,v):q=plus(q,mul(ck(abs(F(a))),b))
  return q
 L0=[av(v) for v in s['native_L0']];Le=[av(v) for v in s['native_Le']];S0=[av(v) for v in s['native_S0']];Se=[av(v) for v in s['native_Se']]
 rpos=[];epos=[]
 for n,rec in enumerate(s['known_first_carry_pivots8']):
  rr=av(rec['Xp_r_literal56']);rv=[rr.get(i,F(0)) for i in range(56)];kind,arg=rec['Xp_ell_recipe']
  if kind=='literal':ee=av(arg);ev=[ee.get(i,F(0)) for i in range(56)]
  else:
   assert kind=='native';x=pos_mv(L0,Se[arg]);y=pos_mv(Le,S0[arg]);ev=[plus(a,b) for a,b in zip(x,y)]
   assert d['known_FIRST_carry_pivots8'][n]['native_Se_column']==s['native_Se'][arg] and d['known_FIRST_carry_pivots8'][n]['native_S0_column']==s['native_S0'][arg]
  for key,v in [('positive_Xp_r_eight',rv),('positive_Xp_ell_eight',ev)]:
   expected={i:F(int(u),int(den)) for i,den,u in d[key][n]};assert v==[expected.get(i,F(0)) for i in range(56)]
  rpos.append(rv);epos.append(ev)
 phys={i:[ck(abs(F(x))) for x in v] for i,v in d['physical_B21_functions3']};delta=[plus(a,b) for a,b in zip(phys[20117],phys[40895])];z1=phys[27948]
 D=s['D_two_by10_by56'];Gamma=[dict(v) for v in s['Gamma_two']];gd=[];ds=[];cs=[]
 for t in range(2):
  gd.append(mul(ck(abs(F(Gamma[t][4]))),delta[20]));dv=F(0);cv=F(0)
  for pair in range(10):dv=plus(dv,mul(ck(abs(F(D[t][pair][4]))),delta[10+pair]));dv=plus(dv,mul(ck(abs(F(D[t][pair][6]))),z1[10+pair]))
  for rec,rp,ep0 in zip(s['known_first_carry_pivots8'],rpos,epos):cv=plus(cv,mul(ck(abs(F(rec['a1_z1']))),plus(dot(D[t][2],rp),dot(D[t][3],ep0))))
  ds.append(dv);cs.append(cv)
 assert gd==ds==cs==[F(0),F(0)]
 total=[plus(plus(a,b),c) for a,b,c in zip(gd,ds,cs)];six=[]
 for k in range(6):
  v=F(0)
  for t in range(2):v=plus(v,mul(ck(abs(F(dict(Q[t]).get(k,'0')))),total[t]))
  six.append(v)
 for key,own in [('positive_Gamma_cubic_two',gd),('positive_D_physical_second_two',ds),('positive_known_FIRST_carry_two',cs),('positive_full_upper_two',total),('positive_full_six',six)]:assert [F(int(u),int(den)) for den,u in d[key]]==own==[F(0)]*len(own)
 for rep in [sr,ir]:
  for path,ha in rep['input_sha256'].items():
   if path not in pins:raw+=(ROOT/path).stat().st_size;pins[path]=sha(ROOT/path)
   assert pins[path]==ha,path
 owner=count([sr,s,ir,d,orig,a,h,ops,u,cc,j1,o2,p3,f3,Q,pins,L0,Le,S0,Se,rpos,epos,phys,delta,z1,gd,ds,cs,total,six])+15000
 assert owner<150000 and raw+1_000_000<32*1024*1024
 ep='independent-ONE288-source-eight-FIRST-SUM-carries-positive-zero.json';ev={'exact_single352_nativeempty_fulla1eight_onlyz1_verified':True,
  'pureplane_GammaXp2_times_z1_factor_zero_before_any_N2_response':True,'fullPsi_D_carry_uses_z1r_and_z1ell_ordinary_coefficients':True,
  'all_native_Xpell_SUM_L0Se_plus_LeS0_eight_first_envelopes_verified':True,'all_three_physicalB21_andtruefull9firsts_verified':True,
  'no_honest_or_conductor_order_hypothesis_needed':True,'Gamma_D_and_eight_FIRST_carry_positive_envelopes_separately0':True,
  'actual_ONE288_upper_class_zero_by_nonnegative_bound_only':True,'all_two_upper_six_embedding_envelopes0':True,
  'no_signed_response_or_zero_only_VALUEjob':True,'only_case288_no_other605_HF_or_length_claim':True,'own_positive_bits':maximum,'whole_owned_upper':owner}
 dump(ep,ev);r={'id':STEM,'all_passed':True,'source_and_positive_zero_proof_only':True,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'whole_numeric_upper':owner,'whole_RAW_upper':raw+1_000_000,'seconds':time.monotonic()-start}
 dump('report.json',r);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'owned':owner,'bits':maximum,'seconds':r['seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

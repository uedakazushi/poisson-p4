#!/usr/bin/env python3
"""Own OLDtwo/57-FIRST source and nonnegative-zero proof, no responses."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER554_OLD_FIRST_20261105_v63_REVIEW'
SP='reports/AMBIENT_B18_ELL2_T2_UPPER554_OLD_CARRY_SOURCE_20261105_v63/';IP='reports/AMBIENT_B18_ELL2_T2_UPPER554_OLD_CARRY_INPUT_GATE_20261105_v63/'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple,set)):return sum(map(count,x))
 if isinstance(x,str):
  try:F(x);return 1
  except(ValueError,ZeroDivisionError):return 0
 return int(isinstance(x,(int,float,F)))
def fields(p,names):
 s=p.read_text();d=json.JSONDecoder();i=s.index('{')+1;z={}
 while True:
  while s[i].isspace() or s[i]==',':i+=1
  if s[i]=='}':break
  k,i=d.raw_decode(s,i)
  while s[i].isspace() or s[i]==':':i+=1
  v,i=d.raw_decode(s,i)
  if k in names:z[k]=v
  del v
 assert set(z)==set(names);return z
def main(out):
 t0=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0;maxbits=0
 def read(p,names=None):
  nonlocal raw
  q=ROOT/p;pins[p]=sha(q);raw+=2*q.stat().st_size
  return json.loads(q.read_text()) if names is None else fields(q,names)
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 sr=read(SP+'report.json');assert sr['all_passed'] and pins[SP+'report.json']=='1a5b1726d6592dbd67fe647cfbd9e65bd36997e615aa8f66696767d6a43ccc15'
 se='ONE554-six-c0-OLD-and-full29-a1-known-FIRST-source.json';s=read(SP+se);assert pins[SP+se]==sr['evidence_sha256'][se]
 ir=read(IP+'report.json');assert ir['all_passed'] and pins[IP+'report.json']=='8373bcb85e3fb3f2bff7751e6f451d69db5850fd2747537102b82af69bc29d5d'
 ie='ONE554-OLD-two-columns-and57-FIRST-positive-component-input.json';d=read(IP+ie);assert pins[IP+ie]==ir['evidence_sha256'][ie]
 src='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_SOURCE_20261105_v63/combined-OLD2-D-and-four-pivot-source-factor-inventory.json';o=read(src,['original619_c0_source_combinations','kernel617_source_recipes','typed_program','D_old_root_columns','D_two_by10_by56','pair_order10'])
 aa=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json',['all_zero_bottom_and_Y0_source_combinations'])['all_zero_bottom_and_Y0_source_combinations'][554]
 assert aa['original768_combination']==s['original768_source_all6']==o['original619_c0_source_combinations'][554]==[[289,'25/4'],[379,'-1'],[643,'29/4'],[645,'25'],[651,'-25/2'],[670,'1']]
 assert aa['actual_a1_27_times9']==s['actual_full_a1_27times9']==d['full_a1_29'] and len(d['full_a1_29'])==29
 recipe=next(v for v in o['kernel617_source_recipes'] if v['free_case619']==554);assert recipe==s['same617_source_recipe'] and recipe['subtract_actual_basis563609']==['0','0']
 factor=[];uses=set();omitted=[]
 for slot,q in aa['actual_a1_27_times9']:
  pi,var=divmod(slot,9);hit=False
  for p,(i,j) in enumerate(o['pair_order10']):
   if i==var:factor.append([pi,p,j,q]);uses.add((pi,j));hit=True
   if j==var and i!=j:factor.append([pi,p,i,q]);uses.add((pi,i));hit=True
  if not hit:assert var in [2,3,6];omitted.append([slot,q])
 assert factor==s['all_known_FIRSTcarry_factor_recipes']==d['known_FIRSTcarry_factor_recipes57'] and len(factor)==57
 assert omitted==s['omitted_a1_slots_only_inactive_Dpair_zero']
 assert sorted(uses)==[(v['pivot27'],v['variable9']) for v in d['Xp_FIRST_column_recipes28']]
 cc=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_INPUT_GATE_20261105_v63/pruned-known-Cfull224-input-and-positive-source-envelopes.json')
 assert s['Gamma_two']==d['Gamma_two']==cc['Gamma_two'] and s['D_two_by10_by56']==d['D_two_by10_by56']==o['D_two_by10_by56']
 old=read('reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',['source_columns640_rows56_sparse_dual9'])['source_columns640_rows56_sparse_dual9']
 native=read('reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json',['source_columns128_rows48_sparse_dual9'])['source_columns128_rows48_sparse_dual9']
 assert d['Xp_FIRST_column_recipes28']==s['all_needed_Xp_FIRST_column_recipes']
 need0=set();neede=set()
 for rec in d['Xp_FIRST_column_recipes28']:
  p,var,orig=rec['pivot27'],rec['variable9'],rec['original768'];assert orig==cc['pivot27_source_indices'][p];f=rec['first_column_recipe']
  if var==8:
   kind,arg=cc['pivot27_ell_recipes'][p]
   if kind=='literal':assert f['kind']=='literal56' and f['column']==arg
   else:
    assert f['kind']=='native_ell_SUM' and f['native128_source']==arg==orig-640 and f['Se']==cc['native_Se'][arg] and f['S0']==cc['native_S0'][arg]
    need0.update(i for i,q in f['Se']);neede.update(i for i,q in f['S0'])
  else:
   dual=8 if var==7 else var+1;cols=old[orig] if orig<640 else native[orig-640]
   actual=[[row,dict(v)[dual]] for row,v in cols if dual in dict(v) and dict(v)[dual]!='0'];assert f['kind']=='literal56' and f['column']==actual
 del old,native
 assert d['native_L0_selected_columns']==[[j,cc['native_L0'][j]] for j in sorted(need0)] and d['native_Le_selected_columns']==[[j,cc['native_Le'][j]] for j in sorted(neede)]
 assert s['native_L0']==cc['native_L0'] and s['native_Le']==cc['native_Le']
 # Reconstruct the complete two-column graph and direct C0 roots independently.
 hj=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_INPUT_GATE_20261105_v63/617-direct-old-known-first-formal-Jacobian-positive-input.json',['typed_source_program','pruned_original_old_columns_roots','Gamma_two'])
 ops=dict(hj['typed_source_program']);ops.update(dict(o['typed_program']))
 dcols=dict(o['D_old_root_columns']);gcols=dict(hj['pruned_original_old_columns_roots']);oldsrc=[[289,'25/4'],[379,'-1']]
 todo=[i for j,q in oldsrc for row,i in dcols[j]];direct={0};seen=set()
 for j,q in oldsrc:
  for row,(kind,i) in gcols[j]:
   if kind=='typed':todo.append(i)
   else:assert kind=='C0';direct.add(i)
 while todo:
  i=todo.pop()
  if i in seen:continue
  seen.add(i);node=ops[i]
  if node[0]=='sum':todo.extend(j for j,a,b in node[1])
  elif node[0]=='product':todo.extend(node[3])
  else:assert node[0]=='actual_C0_scalar_function'
 assert [[i,ops[i]] for i in sorted(seen)]==s['OLD_full_typed_program']==d['OLD_full_typed_program11'] and len(seen)==11
 assert s['OLD_D_source_columns']==d['OLD_D_source_columns']==[[j,q,dcols[j]] for j,q in oldsrc]
 assert s['OLD_Gamma_source_columns']==d['OLD_Gamma_source_columns']==[[j,q,gcols[j]] for j,q in oldsrc]
 leaves={ops[i][1] for i in seen if ops[i][0]=='actual_C0_scalar_function'}|(direct-{0})
 assert leaves==set(s['OLD_physical_function_leaf_ids7'])=={7199,20117,20138,27897,27948,40895,40916}
 pp='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json';pv=read(pp,['all_requested164_root_values'])['all_requested164_root_values'];ph={i:dict(v) for i,v in pv if i in leaves};del pv
 fp='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json';fv=read(fp,['first640'])['first640']['all79_atom_values_full10'];ff={i:dict(v) for i,v in fv if i in ph};del fv
 for i,v in d['physical_B21_functions7']:
  assert list(map(F,v))==[F(ph[i].get(k,'0')) for k in s['raw33_to_B21_selector21']]
  assert list(map(F,v[:10]))==[F(ff[i].get(k,'0')) for k in range(10)]
 initial=count([sr,s,ir,d,o,aa,cc,ops,seen,ph,ff,pins])+20000;assert initial<150000 and raw+1_000_000<32*1024*1024
 dump('phase-before-positive.json',{'phase':'all29a1_57uniqueordinaryfactors_28originalFIRSTrecipes_andsevenphysical_inputs_PASS','source122_temporary367':True,'owner':initial})
 def ck(v):
  nonlocal maxbits
  assert v>=0;maxbits=max(maxbits,v.numerator.bit_length(),v.denominator.bit_length());assert maxbits<=122;return v
 def av(v):return {i:ck(abs(F(q))) for i,q in v if F(q)}
 def plus(a,b):return ck(a+b)
 def mul(a,b):return ck(a*b)
 def pmv(cols,v):
  z=[F(0)]*56
  for j,a in v.items():
   for i,b in cols[j].items():z[i]=plus(z[i],mul(a,b))
  return z
 L0={j:av(cc['native_L0'][j]) for j in need0};Le={j:av(cc['native_Le'][j]) for j in neede};env={}
 for rec in d['Xp_FIRST_column_recipes28']:
  f=rec['first_column_recipe']
  if f['kind']=='literal56':a=av(f['column']);v=[a.get(i,F(0)) for i in range(56)]
  else:
   a=pmv(L0,av(f['Se']));b=pmv(Le,av(f['S0']));v=[plus(x,y) for x,y in zip(a,b)]
  env[rec['pivot27'],rec['variable9']]=v
 for p,var,v in d['positive_Xp_FIRST_columns28']:
  bound={i:F(int(u),int(den)) for i,den,u in v};assert env[p,var]==[bound.get(i,F(0)) for i in range(56)]
 physical={i:[ck(abs(F(q))) for q in v] for i,v in d['physical_B21_functions7']}
 basis=[tuple(v) for v in s['ordinary_basis_B21']];idx={v:i for i,v in enumerate(basis)}
 table=[[idx.get(tuple(x+y for x,y in zip(a,b))) for b in basis] for a in basis]
 def vplus(a,b,c=F(1)):return [plus(x,mul(c,y)) for x,y in zip(a,b)]
 def vprod(a,b):
  out=[F(0)]*21
  for i,x in enumerate(a):
   if x:
    for j,y in enumerate(b):
     k=table[i][j]
     if y and k is not None:out[k]=plus(out[k],mul(x,y))
  return out
 jets={}
 for i in sorted(seen):
  node=ops[i]
  if node[0]=='actual_C0_scalar_function':v=physical[node[1]][:]
  elif node[0]=='sum':
   v=[F(0)]*21
   for j,a,b in node[1]:v=vplus(v,jets[j],ck(abs(F(a,b))))
  else:
   assert node[0]=='product';v=[F(0)]*21;v[0]=ck(abs(F(node[1],node[2])))
   for j in node[3]:v=vprod(v,jets[j])
  jets[i]=v
 for i,v in d['positive_OLD_typed_B21']:
  saved={k:F(int(u),int(den)) for k,den,u in v};assert jets[i]==[saved.get(k,F(0)) for k in range(21)]
 G=[dict(v) for v in s['Gamma_two']];D=s['D_two_by10_by56'];old3=[];old2=[];factorbounds=[];carry=[];partial=[]
 for t in range(2):
  q3=F(0);q2=F(0)
  for old,c,roots in d['OLD_Gamma_source_columns']:
   for row,(kind,i) in roots:
    v=jets[i] if kind=='typed' else ([F(0)]*21 if i==0 else physical[i])
    q3=plus(q3,mul(ck(abs(F(c))),mul(ck(abs(F(G[t][row]))),v[20])))
  for old,c,roots in d['OLD_D_source_columns']:
   for row,i in roots:
    for p in range(10):q2=plus(q2,mul(ck(abs(F(c))),mul(ck(abs(F(D[t][p][row]))),jets[i][10+p])))
  old3.append(q3);old2.append(q2);fb=[];total=F(0)
  for p,pair,var,c in factor:
   term=F(0)
   for i,x in enumerate(env[p,var]):
    if x:term=plus(term,mul(ck(abs(F(D[t][pair][i]))),x))
   term=mul(ck(abs(F(c))),term);fb.append(term);total=plus(total,term)
  factorbounds.append(fb);carry.append(total);partial.append(plus(plus(q3,q2),total))
 assert all(q>0 for q in old3+old2+carry)
 for k,own in [('positive_OLD_Gamma_cubic_two',old3),('positive_OLD_D_second_two',old2),('positive_FIRSTcarry_total_two',carry),('positive_OLD_and_FIRST_partial_two',partial)]:assert [F(int(u),int(den)) for den,u in d[k]]==own
 assert [[F(int(u),int(den)) for den,u in row] for row in d['positive_FIRST57_factor_bounds_two_by57']]==factorbounds
 for rep in [sr,ir]:
  for path,h in rep['input_sha256'].items():
   if path not in pins:raw+=(ROOT/path).stat().st_size;pins[path]=sha(ROOT/path)
   assert pins[path]==h,path
 owner=count([sr,s,ir,d,o,aa,cc,ops,seen,ph,ff,pins,env,L0,Le,physical,jets,table,old3,old2,factorbounds,carry,partial,hj])+15000
 assert owner<150000 and raw+1_000_000<32*1024*1024
 ep='independent-OLDtwo-and-all114-FIRST-positive-component-checks.json';e={'all29a1_and57ordinaryfactors_no_rr_doublecount_verified':True,
  'all28_originalFIRSTrecipes_verified_in_current27_frame':True,'all_native_ell_FIRST_actualSUM_L0Se_LeS0_retained':True,
  'all_seven_actual_physicalB21_andfull9first_inputs_verified':True,'OLDtwo_cubic_Dsecond_andall114_FIRSTfactor_positive_envelopes_verified':True,'OLD_Gamma_D_and_FIRSTtotals_nonzero':True,
  'only_OLD_FIRST_component_no_native_or_whole554_response_inference':True,'no_signed_response_or_new_function_values':True,'own_positive_bits':maxbits,'whole_owned':owner}
 dump(ep,e);rep={'id':STEM,'all_passed':True,'source_and_positive_component_proof_only':True,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'whole_numeric_upper':owner,'whole_RAW_upper':raw+1_000_000,'seconds':time.monotonic()-t0}
 dump('report.json',rep);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'owned':owner,'bits':maxbits,'seconds':rep['seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

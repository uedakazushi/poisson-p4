#!/usr/bin/env python3
"""Independent third-class source audit; no new inverse or source jet.

Only QQ[r]/r^3 is needed for the third functional.  Opaque regular
quotients retain distinct symbolic r^2 coefficients throughout the check.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
import hashlib,json,resource,signal,time

ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63'
STEM=AUTHOR+'_REVIEW'
ZERO=(Q(0),Q(0),Q(0));ONE=(Q(1),Q(0),Q(0))

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def load(p):return json.loads((ROOT/p).read_text())
def sparse(p):return {int(i):Q(q) for i,q in p if Q(q)}
def accum(a,i,q):
 v=a.get(i,Q(0))+q
 if v:a[i]=v
 else:a.pop(i,None)
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def scale(a,q):return tuple(q*x for x in a)
def mul(a,b):return tuple(sum((a[i]*b[n-i] for i in range(n+1)),Q(0)) for n in range(3))
def inv(a):
 assert a[0],('nonunit inverse',a)
 b=(1/a[0],-a[1]/a[0]**2,a[1]**2/a[0]**3-a[2]/a[0]**2)
 assert mul(a,b)==ONE
 return b
def enc(a):return [str(q) for q in a]
def encsp(a):return [[i,str(q)] for i,q in sorted(a.items()) if q]

class Fingerprints:
 """Canonical local expressions, with variables distinguished by node ID."""
 def __init__(self):self.pool={};self.keys=[];self.zero=self.const(0);self.one=self.const(1)
 def put(self,key):
  if key not in self.pool:self.pool[key]=len(self.keys);self.keys.append(key)
  return self.pool[key]
 def const(self,q):return self.put(('constant',Q(q)))
 def neg(self,a):
  k=self.keys[a]
  if k[0]=='constant':return self.const(-k[1])
  if k[0]=='negative':return k[1]
  return self.put(('negative',a))
 def plus(self,a,b):
  if a==self.zero:return b
  if b==self.zero:return a
  if self.keys[a]==('negative',b) or self.keys[b]==('negative',a):return self.zero
  if self.keys[a][0]==self.keys[b][0]=='constant':return self.const(self.keys[a][1]+self.keys[b][1])
  return self.put(('sum',min(a,b),max(a,b)))
 def times(self,a,b):
  if a==self.zero or b==self.zero:return self.zero
  if a==self.one:return b
  if b==self.one:return a
  if self.keys[a][0]==self.keys[b][0]=='constant':return self.const(self.keys[a][1]*self.keys[b][1])
  return self.put(('product',min(a,b),max(a,b)))
 def inverse(self,a):
  if self.keys[a][0]=='constant':return self.const(1/self.keys[a][1])
  return self.put(('inverse',a))
 def evaluate(self,nodes,needed,replacements=None):
  vals={};replacements=replacements or {}
  for i in sorted(needed):
   if i in replacements:vals[i]=replacements[i];continue
   n=nodes[i];op=n[0]
   if op=='c':v=self.const(Q(n[1],n[2]))
   elif op=='v':v=self.put(('variable_node',i))
   elif op=='a':v=self.plus(vals[n[1]],vals[n[2]])
   elif op=='m':v=self.times(vals[n[1]],vals[n[2]])
   elif op=='n':v=self.neg(vals[n[1]])
   elif op=='i':v=self.inverse(vals[n[1]])
   else:raise AssertionError(('not a scalar source operation',i,n))
   vals[i]=v
  return vals

def closure(nodes,roots,stops=(),allow_variables=False):
 seen=set();todo=list(roots);stops=set(stops)
 while todo:
  i=todo.pop()
  if i in seen:continue
  seen.add(i)
  if i in stops:continue
  n=nodes[i];op=n[0]
  if op in ('a','m'):todo.extend(n[1:3])
  elif op in ('n','i'):todo.append(n[1])
  elif op=='c' or (allow_variables and op=='v'):pass
  else:raise AssertionError(('uncached actual leaf',i,n))
 return seen

def evaluate(nodes,roots,cache):
 needed=closure(nodes,roots,cache);values={};inverse_nodes=[]
 for i in sorted(needed):
  n=nodes[i];op=n[0]
  if i in cache:v=cache[i]
  elif op=='c':v=(Q(n[1],n[2]),Q(0),Q(0))
  elif op=='a':v=add(values[n[1]],values[n[2]])
  elif op=='m':v=mul(values[n[1]],values[n[2]])
  elif op=='n':v=scale(values[n[1]],-1)
  elif op=='i':v=inv(values[n[1]]);inverse_nodes.append(i)
  else:raise AssertionError((i,n))
  values[i]=v
 return {i:values[i] for i in roots},len(needed),inverse_nodes

def main():
 began=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_CPU,(60,60))
 out=ROOT/'reports'/STEM;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 work=ROOT/'work/B18-next3-third-independent-review';work.mkdir(parents=True,exist_ok=True)
 hashes={}
 def get(path):hashes[path]=sha(path);return load(path)
 def data(stem,name):return get('reports/'+stem+'_20261104_v63/'+name)
 def save(name,x):
  p=out/name;p.write_text(json.dumps(x,indent=2)+'\n');return p
 def stage(x):(work/'stage.json').write_text(json.dumps(x,indent=2)+'\n')
 R=get('reports/'+AUTHOR+'/report.json');assert R['all_passed']
 assert hashes['reports/'+AUTHOR+'/report.json']=='5233957a7ab9dac28f16d9f51085fe375ffaac1d293b82d065f0381c0062bf92'
 for p,h in R['input_sha256'].items():assert sha(p)==h;hashes[p]=h
 for name,h in R['evidence_sha256'].items():
  p='reports/'+AUTHOR+'/'+name;assert sha(p)==h;hashes[p]=h
 p='code/ambient_B18_T_next3_existing_cache_transport_20261104_v63.py';assert sha(p)==R['code_sha256'];hashes[p]=sha(p)
 p='research/notes/'+AUTHOR+'.md';assert sha(p)==R['note_sha256'];hashes[p]=sha(p)
 accepted={
  'AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63_REVIEW':'bb6e7e7c5327528f2be31b1d29db5891461d12adb319daa184b6c6addd63d693',
  'AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW':'ca127c4dc0f3c12a4d4d80dfc1ca316f42dee8b90ecc74cb3e5e423752c906b7',
  'AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63_REVIEW':'c933016430a1779bed947b0c81484b0aae23afbee23e4df15a386ea959bf0dd9',
  'AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63_REVIEW':'0000a1aa2790891a329f7ff708ef6dcfcaea8b04571de07b79ee5b503caabcb1',
  'AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63_REVIEW':'baaec4550a00c25ee24b2ab5d03e9eebaf67943caa216a01abe624bf19f3dd17'}
 for st,h in accepted.items():
  p='reports/'+st+'/report.json';d=get(p);assert d['all_passed'] and hashes[p]==h,(st,hashes[p])
 C=data('AMBIENT_B18_T0_HONEST_C0','honest-C0-D0-D4-DAG.json')
 Full=data('AMBIENT_B18_FULL_RELATIVE_INPUT','full-relative186-DAG.json')
 S=data('AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT','actual-source-coefficients.json')
 J=data('AMBIENT_B18_T3_FULL9_FIRSTJET','full9-firstjet-certificate.json')
 Art=data('AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION','actual-source105-Artin20-values.json')
 H=data('AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20','actual-special-fiber20-coefficients.json')
 AF=data('AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR','phase1-actual-atom-firsts.json')
 M=data('AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT','actual-projected-three-class-X2-operator-preflight.json')
 T=data('AMBIENT_B18_HONEST_T_EULER_OPERATOR','actual-honest-T116-936-Euler-operator.json')
 P2=get('reports/'+AUTHOR+'/phase2-existing-cache-B11-copy-values.json')
 P3=get('reports/'+AUTHOR+'/phase3-composite-second-source-functional.json')
 nodes=C['DAG_scalar_nodes'];original=S['DAG_scalar_nodes'];fullnodes=Full['DAG_scalar_nodes']
 assert nodes[:len(original)]==original and nodes[:len(fullnodes)]==fullnodes
 assert len(fullnodes)==339669
 inst=Full['all16_pointed_unit57_inverse_instances'][14];bind=dict(inst['bindings']);cb=J['both57_source2']['CB']
 assert inst['retained_Y_indices']==[7,8] and set(bind)==set(cb['unknown57']) and len(bind)==57
 assert cb['equations57']==Full['original_unit48_equations_CB']+Full['CB_coordinate_Y9_nodes']
 assert original[11]==['v','CA_trace_p'] and S['families']['CA']['p_quadratic_relation']['s']==11
 assert 11 in Art['joint105_unknown_order']
 # Independent expression fingerprints avoid the producer's index-replay
 # table and remain insensitive to later duplicate negative-node aliases.
 fp=Fingerprints();targets=inst['all57_pointed_equations'];assert len(targets)==57
 actualneed=closure(nodes,targets+list(bind.values())+cb['F9_nodes'][7:9],allow_variables=True)
 actual=fp.evaluate(nodes,actualneed)
 oldneed=closure(original,cb['equations57'],bind,allow_variables=True)
 rewired=fp.evaluate(original,oldneed,{i:actual[j] for i,j in bind.items()})
 expected=[rewired[i] for i in cb['equations57'][:48]]
 expected += [fp.plus(rewired[i],fp.neg(actual[i] if k in (7,8) else fp.zero)) for k,i in enumerate(cb['F9_nodes'])]
 assert expected==[actual[i] for i in targets]
 stage({'phase':'57 actual copy equations independently fingerprinted','source_hashes':hashes})
 # The third class selects row55 and r^2.  All first-order correction
 # coefficients contributing to that monomial are literally zero.
 assert M['Pi3_180_dual_rows'][2]==[[179,'1']]
 assert M['pi_free32_rows'][-1]==[[31,'1']] and M['free_kernel_rows32'][31]==55
 assert M['quadratic_monomials36'][35]==[7,7]
 assert all(sparse(p).get(39,Q(0))==0 for p in M['Ystar1_40_by84_columns'])
 assert all(sparse(p).get(55*8+7,Q(0))==0 for p in M['actual_pivot_X1_by27_sparse56_times8'])
 X=T['actual_X_source']['source_X_sparse_columns'];assert len(X)==768
 assert all(all(i<48 for i,_ in col) for col in X[640:])
 op=T['scalar_and_vector_operator_nodes'];folded=[]
 for combo in M['620_actual_source_combinations_in_original768']:
  roots={}
  for j,q in sparse(combo).items():
   if j>=640:continue
   root=C['honest_D0_relation_columns'][j][8]
   # Verify the root's actual X row, not only the C0 display table.
   entry=dict(X[j]).get(55)
   if entry is None:assert nodes[root][0]=='c' and nodes[root][1]==0;continue
   assert op[entry]==['actual_C0_scalar_function_restrict_ell0',root]
   if nodes[root][0]=='c':continue
   accum(roots,root,q)
  folded.append(roots)
 active=sorted(set().union(*(set(x) for x in folded)))
 assert len(folded)==620 and len(active)==90
 stage({'phase':'literal third functional and 620 source combinations folded','active90':active})
 # H20 -> QQ[r]/r^3 is checked on every one of its 160 action identities.
 hb=H['quotient_basis_monomials'];assert hb==Art['A20_basis']
 images=[]
 for m in hb:images.append(tuple(Q(int(not any(m[:7]) and m[7]==k)) for k in range(3)))
 def phi(raw):
  v=ZERO
  for i,q in raw:v=add(v,scale(images[i],Q(q)))
  return v
 vars8=[ZERO]*7+[(Q(0),Q(1),Q(0))]
 for k,act in enumerate(H['multiplication8_column_sparse']):
  for j,p in enumerate(act):assert mul(vars8[k],images[j])==phi(p)
 cached={i:phi(p) for i,p in Art['all_requested_original_source_root_values']}
 for i,p in zip(Art['joint105_unknown_order'],Art['all105_unknown_final_values']):
  v=phi(p)
  if i in cached:assert cached[i]==v
  cached[i]=v
 yroots=cb['F9_nodes'][7:9];yb,ycount,_=evaluate(nodes,yroots,cached)
 args=[ZERO]*7+[yb[i] for i in yroots];copyvalues={}
 for i,coefs in cb['all57_source2jets']:
  v=ZERO
  for mon,c in zip(J['source2_monomials'],coefs):
   if not Q(c):continue
   term=scale(ONE,Q(c))
   for k in mon:term=mul(term,args[k])
   v=add(v,term)
  copyvalues[bind[i]]=v
 cached.update(copyvalues)
 eqvals,eqcount,eqinverses=evaluate(nodes,targets,cached);assert all(v==ZERO for v in eqvals.values())
 for i,p in P2['all57_copied_unknown_values_through2']:
  raw=sparse(p);assert copyvalues[i]==tuple(raw.get(k,Q(0)) for k in range(3))
 # First-order values of opaque regular quotients are accepted source
 # caches; their r^2 coefficients remain independent symbols theta_i.
 first={}
 for key in ('all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA'):
  for i,p in J[key].items():first[int(i)]=(Q(p[0]),Q(p[8])+Q(p[9]))
 for item in AF['all35_ordered_CA_prefix_firsts']:
  p=sparse(item['first9_bar']);first[item['node']]=(p.get(0,Q(0)),p.get(8,Q(0)))
 for i,p in AF['all8_g_quotient_overrides']:
  p=sparse(p);first[i]=(p.get(0,Q(0)),p.get(8,Q(0)))
 # Actual C0 gA has no terms of degree <=2, by the independently accepted
 # cubic cache; it is not set to zero merely because it vanishes at P0.
 ga=C['cached_actual_CA_gA_all_degree3_coefficients'];assert len(ga)==10 and all(Q(q)==0 for q in ga[:6])
 assert [Q(q) for q in ga[6:]]==[Q(3,8),Q(2),Q(13,6),Q(2,3)]
 for i in C['actual_C0_basis8_scalar_nodes']:cached[i]=ZERO
 need=closure(nodes,active,set(cached)|set(first));opaque=sorted(set(need)&set(first)-set(cached))
 deps={}
 for i in sorted(need):
  if i in cached:deps[i]=set()
  elif i in first:deps[i]={i}
  else:
   n=nodes[i];ch=n[1:3] if n[0] in ('a','m') else n[1:] if n[0] in ('n','i') else []
   deps[i]=set().union(*(deps[j] for j in ch))
 forecast=3*len(need)+sum(map(len,deps.values()))+15000
 assert forecast<150000,forecast
 save('preflight.json',{'active_actual_P_roots':90,'closure_nodes':len(need),'opaque_regular_quotients':len(opaque),'support_upper':forecast,'cap':150000,'no_new_source_solve_or_jet':True})
 vals={};grads={};unitinverses=[]
 for i in sorted(need):
  n=nodes[i];o=n[0]
  if i in cached:v=cached[i];g={}
  elif i in first:v=(first[i][0],first[i][1],Q(0));g={i:Q(1)}
  elif o=='c':v=(Q(n[1],n[2]),Q(0),Q(0));g={}
  elif o=='a':
   v=add(vals[n[1]],vals[n[2]]);g=dict(grads[n[1]])
   for j,q in grads[n[2]].items():accum(g,j,q)
  elif o=='n':v=scale(vals[n[1]],-1);g={j:-q for j,q in grads[n[1]].items()}
  elif o=='m':
   a,b=vals[n[1]],vals[n[2]];v=mul(a,b);g={}
   for j,q in grads[n[1]].items():accum(g,j,q*b[0])
   for j,q in grads[n[2]].items():accum(g,j,q*a[0])
  elif o=='i':
   a=vals[n[1]];v=inv(a);g={j:-q/a[0]**2 for j,q in grads[n[1]].items()};unitinverses.append([i,str(a[0])])
  else:raise AssertionError(('uncached actual regular function',i,n))
  assert set(g)<=deps[i];vals[i]=v;grads[i]=g
 # Compare only these independently selected 90 roots, not the producer's
 # first two projected functionals or its old root-adjoint recipe.
 primaryroots={i:(sparse(v),sparse(g)) for i,v,g in P3['all108_composite_values_and_opaque_second_gradients']}
 for i in active:
  p,g=primaryroots[i];assert vals[i]==tuple(p.get(k,Q(0)) for k in range(3)) and grads[i]==g
 nonzero={i:v[2] for i,v in vals.items() if i in active and v[2]}
 assert nonzero=={cb['F9_nodes'][0]:Q(-445,72),cb['F9_nodes'][4]:Q(-6625,288),cb['F9_nodes'][5]:Q(-145,72)}
 checks=[]
 for k,combo in enumerate(folded):
  known=Q(0);symbols={}
  for i,q in combo.items():
   known+=q*vals[i][2]
   for j,c in grads[i].items():accum(symbols,j,q*c)
  assert known==0 and not symbols,(k,known,symbols)
  old=P3['620_projected_source_functionals'][k]
  assert old['index']==k and old['three_unknown_division_second_functionals'][2]==[]
  assert Q(old['three_known_oldB_P_and_pivot_correction_values'][2])==0 and old['three_native_X2_functionals'][2]==[]
  checks.append([k,str(known),encsp(symbols)])
 total=sum(3+len(grads[i]) for i in need);assert total<150000
 evidence={'third_Pi3':[[179,'1']],'third_raw_source_row':55,'third_monomial':'r_CA^2',
  'third_Euler_and_first_pivot_correction_coefficients_literal_zero':True,
  'all128_native_bottom8_source_rows_literal_zero':True,
  'all620_original_source_combinations_folded_without_root_adjoint':[[[i,str(q)] for i,q in sorted(c.items())] for c in folded],
  'active90_actual_P_roots':active,'all90_known_univariate_parts_and_symbolic_second_sensitivities':[[i,enc(vals[i]),encsp(grads[i])] for i in active],
  'three_nonzero_known_rr_parts':encsp(nonzero),'opaque_regular_quotient_nodes':opaque,
  'all620_known_and_symbolic_rr_totals':checks,
  'exact57_copy_bindings':inst['bindings'],'copied57_all_equations_zero_in_QQ_r_mod_r3':True,
  'copied57_values':[[i,enc(copyvalues[i])] for i in sorted(copyvalues)],
  'actual_YB7_8':[[i,enc(yb[i])] for i in yroots],
  'H20_to_QQ_r_mod_r3_action_identities':160,'CA_trace_p_node11_is_actual_saved105_root':True,
  'actual_CA_gA_source_coefficients_through3':ga,'unit_inverse_origins':unitinverses,
  'actual_source_functional_equality_scope':'Actual minimal5 Schur quadratic class modulo known rank177, using independently accepted M1 cancellation; no literal unprojected equality asserted.',
  'all_opaque_rr_values_kept_symbolic_not_assigned_zero':True,'source_slots':total}
 save('third-class-source-certificate.json',evidence)
 for p in ['code/ambient_B18_T_next3_existing_cache_transport_review_20261104_v63.py','research/notes/'+STEM+'.md']:hashes[p]=sha(p)
 report={'id':STEM,'all_passed':True,
  'actual_third_projected_quadratic_class_vanishes_on_all620':True,
  'exact57_copy_source_and90_P_rr_transport_verified':True,
  'third_functional_is_nonzero_and_annihilates_all_known464_inputs':True,
  'known464_annihilation_and_M1_projected_scope_reused_from_independent_B11_review':True,
  'all620_opaque_second_sensitivities_cancel_symbolically':True,
  'no_source_solve_newjet_point_prime_or_GB':True,
  'source_scope':'Third class only; first two classes/full projected rank/total T length/full nil index are not audited or inferred here.',
  'depth_consequence_requires_complete_initial_space_proof':'Nonzero third functional annihilates all actual quadratic initial relations once the separately audited complete5x825/linear-kernel proof is applied; then lengthT>=9 and accepted comparison criterion gives cone depth19/projective18.',
  'input_sha256':hashes,'evidence_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()},
  'seconds':time.monotonic()-began,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 save('report.json',report);stage({'phase':'FINAL','report_sha256':sha('reports/'+STEM+'/report.json')})
 assert sum(p.stat().st_size for p in out.iterdir())<32<<20
 print(json.dumps({k:v for k,v in report.items() if k not in ('input_sha256','evidence_sha256')}),flush=True)

if __name__=='__main__':main()

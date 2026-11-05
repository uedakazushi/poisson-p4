#!/usr/bin/env python3
"""Independent native70 SOURCE/positive review. Never compute signed images."""
import argparse,hashlib,json,math,os,resource,signal,subprocess,sys,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_INPUT_GATE_20261105_v63_REVIEW'
BASE='reports/AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_INPUT_GATE_20261105_v63/'
PAIR=[(7,7),(7,8),(0,7),(0,8),(4,7),(4,8),(5,7),(5,8)]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sp(a):return {int(i):Q(x) for i,x in a if Q(x)}
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 return int(isinstance(x,(int,float,Q)) or isinstance(x,str) and x.lstrip('-').replace('/','').isdigit())
def img(cols,v):
 z={}
 for j,c in v.items():
  for i,x in cols[j].items():
   q=z.get(i,Q(0))+c*x
   if q:z[i]=q
   elif i in z:del z[i]
 return z
def main(out,custody,live):
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={}
 def read(p,keys=None):pins[p]=sha(ROOT/p);return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 r=read(BASE+'report.json');assert pins[BASE+'report.json']=='c3c2e9e2105ac60a26c8d969689946ab61df5e70ee4f0fbe1bc89edcda5f0eef' and r['all_passed']
 d=read(BASE+'complete70-cached-input-source-and-positive-whole-image-bounds.json');assert pins[BASE+'complete70-cached-input-source-and-positive-whole-image-bounds.json']==r['evidence_sha256']['complete70-cached-input-source-and-positive-whole-image-bounds.json']
 for p,z in r['input_sha256'].items():assert sha(ROOT/p)==z['sha256']
 fbase='reports/AMBIENT_B18_ELL2_T2_NATIVE70_C72_FRAME_SOURCE_20261105_v63/'
 frame=read(fbase+'actual-C72-cache-frame-and-nine-direction-Nfirst-source-bundle.json')
 kb='reports/AMBIENT_B18_ELL2_T2_NATIVE70_A18_SOURCE_20261105_v63/'
 kernel=read(kb+'all70-C4-kernels-and-full-raw33-A18-ring-map-source.json')
 nb='reports/AMBIENT_B18_ELL2_T2_NATIVE70_C72_POSITIVE_SOURCE_20261105_v63/'
 nd=read(nb+'positive-C72-program-and-native70-partial-coefficient-envelopes.json')
 mons=list(map(tuple,d['A18_basis_monomials']));mi={m:i for i,m in enumerate(mons)};assert mons==[()]+[(i,) for i in range(9)]+PAIR
 rec=d['source70_records'];assert rec==kernel['all70_post27_native_point_kernel_records'] and [x['case619'] for x in rec]==list(range(549,619))
 original=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json',['all_zero_bottom_and_Y0_source_combinations'])['all_zero_bottom_and_Y0_source_combinations']
 for a in rec:
  b=original[a['case619']];assert a['original768_constant_combination']==b['original768_combination'] and a['actual_a1_27_times9']==b['actual_a1_27_times9']
 del original
 sourcegamma=read('reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_X2_SOURCE_20261105_v63/actual-four-Gamma-X2-source-rows-and-lazy619-recipe.json',['Gamma4_on_source56_times8','ordinary_second_pair_order8'])
 assert sourcegamma['Gamma4_on_source56_times8']==d['Gamma4_source56_times8'] and list(map(tuple,sourcegamma['ordinary_second_pair_order8']))==PAIR
 S0=list(map(sp,read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['native_S0_columns'])['native_S0_columns']))
 for a in rec:
  cn={i-640:x for i,x in sp(a['original768_constant_combination']).items() if i>=640};assert cn==sp(a['native128_constant_projection']) and max(cn)<112 and not img(S0,cn)
 # Rebuild raw33 quotient from actual monomials, including every product.
 raw=read('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',['actual33_basis','coordinate_actions9_sparse_columns'])
 def word(t):s,r,e=t;return tuple(sorted(([s-1] if s else [])+[7]*r+[8]*e))
 Pi=[{mi[word(t)]:Q(1)} if word(t) in mi else {} for t in raw['actual33_basis']]
 assert Pi==list(map(sp,kernel['raw33_to_A18_columns']))
 aa=[[sp(v) for v in m] for m in raw['coordinate_actions9_sparse_columns']]
 def mul18(v,w):
  z={}
  for i,a in v.items():
   for j,b in w.items():
    k=mi.get(tuple(sorted(mons[i]+mons[j])))
    if k is not None:z[k]=z.get(k,Q(0))+a*b
  return {i:x for i,x in z.items() if x}
 for i,t in enumerate(raw['actual33_basis']):
  for j in range(33):
   v={j:Q(1)}
   for k in word(t):v=img(aa[k],v)
   assert img(Pi,v)==mul18(Pi[i],Pi[j])
 for k,cols in enumerate(aa):
  for j,v in enumerate(cols):assert img(Pi,v)==mul18({mi[(k,)]:Q(1)},Pi[j])
 copy=[next(i for i,p in enumerate(Pi) if p=={k:Q(1)}) for k in range(18)];assert copy==frame['raw33_to_A18_coeff_copy_source_indices18']
 # Independently restrict actual canonical33 coefficients, not a name match.
 ms='reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63/'
 actual=read(ms+'phase3-actual-kappa9-and-normalized-nilframe.json',['actual14_C132','kappa9_C132'])
 for key,dest in [('actual14_C132','actual14_C72'),('kappa9_C132','actual_kappa9_C72')]:
  assert [[str(Q(v[4*i+k])) for i in copy for k in range(4)] for v in actual[key]]==frame[dest]
 oldproper=read('reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase4-native-source-and-R12.json',['actual14','kappa9'])
 ellproper=read('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_INPUT_GATE_20261105_v63/actual-native-Rell4-cache-source-input-and-adapter-gate.json',['actual14_C16','actual_kappa9_C16','actual_nu_A_B_Rell4'])
 for fk,ok,ek in [('actual14_C72','actual14','actual14_C16'),('actual_kappa9_C72','kappa9','actual_kappa9_C16')]:
  for a,b,c in zip(frame[fk],oldproper[ok],ellproper[ek],strict=True):
   assert list(map(Q,a[:36]))==list(map(Q,b))
   assert [Q(a[4*i+k]) for i in [0,8,9,11] for k in range(4)]==list(map(Q,c))
 nu33=read(ms+'phase1-source150-readback-fingerprints.json',['nu_A_B33'])['nu_A_B33'];nu=[[Q(v[i]) for i in copy] for v in nu33];assert nu==[list(map(Q,v)) for v in d['actual_nu_A_B_A18']]
 assert d['actual_kappa9_C72']==frame['actual_kappa9_C72'] and d['actual_G8_point_and_first40']==frame['actual_G8_point_and_first40']
 for a,b in zip(nu,ellproper['actual_nu_A_B_Rell4'],strict=True):assert [a[i] for i in [0,8,9,11]]==list(map(Q,b))
 # Literal N point/first blocks from both accepted source matrices.
 n36=read('reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase5-native-adapter-firsts.json',['N','G8'])
 n16=read('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase2-actual-T-N-G8-C16.json',['N','actualG8'])
 for i,j,v in d['actual_needed_N_point_and_first40']:
  a,b=list(map(Q,n36['N'][i][j])),list(map(Q,n16['N'][i][j]));assert a[:4]==b[:4] and a[32:36]==b[4:8] and list(map(Q,v))==a+b[8:12]
 for a,b,c in zip(d['actual_G8_point_and_first40'],n36['G8'],n16['actualG8'],strict=True):assert list(map(Q,a))==list(map(Q,b))+list(map(Q,c[8:12]))
 # Bind the regular (not radical) source program to the original literals.
 roots=read('reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63/actual4variable7core-polynomials.json',['source9_variable_order','selected5_pivot_order','five_actual_rational_roots'])
 ci=read('reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/preparation.json',['generators3'])['generators3'];prog=d['regular_G8_source_program']
 assert prog==frame['regular_G8_source_program'] and prog['core_kappa_indices4']==[0,1,2,6] and prog['five_actual_rational_roots']==roots['five_actual_rational_roots'] and prog['regular_CI_generators3']==ci
 assert prog['five_pivot_order']==roots['selected5_pivot_order'] and prog['five_pivot_kappa_indices']==[roots['source9_variable_order'].index(x) for x in roots['selected5_pivot_order']]
 # Rebuild exact trace multiplication using full actual monics.
 tab={};tr=[]
 for i in range(72):
  for j in range(72):
   p=mi.get(tuple(sorted(mons[i//4]+mons[j//4])))
   if p is None:continue
   common=(i%4)&(j%4);mask=(i%4)^(j%4);v={4*p+mask:Q(1)}
   for family in range(2):
    if common&(1<<family):
     w={}
     for slot,c in v.items():
      for k,x in enumerate(nu[family]):
       q=mi.get(tuple(sorted(mons[slot//4]+mons[k])))
       if x and q is not None:w[4*q+slot%4]=w.get(4*q+slot%4,Q(0))+c*x
     v={k:x for k,x in w.items() if x}
   if v:tab[i,j]=v;tr.append([i,j,[[k,str(x)] for k,x in sorted(v.items())]])
 assert tr==d['C72_multiplication_records']
 # Independent positive sparse arithmetic; NEVER evaluate signed G8 images.
 zero=(1,0);maxbits=0
 def mark(x):
  nonlocal maxbits
  maxbits=max(maxbits,x[0].bit_length(),x[1].bit_length());assert maxbits<=384;return x
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator))) if x else zero
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  q=math.lcm(a[0],b[0]);return mark((q,a[1]*(q//a[0])+b[1]*(q//b[0])))
 def times(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else zero
 def add(a,b):
  z=dict(a)
  for i,x in b.items():z[i]=plus(z.get(i,zero),x)
  return {i:x for i,x in z.items() if x[1]}
 def scale(a,c):return {i:times(x,ex(c)) for i,x in a.items() if c}
 def mult(a,b):
  z={}
  for i,x in a.items():
   for j,y in b.items():
    for k,c in tab.get((i,j),{}).items():z[k]=plus(z.get(k,zero),times(times(x,y),ex(c)))
  return z
 one={0:ex(1)};kap=[{i:ex(x) for i,x in enumerate(v) if Q(x)} for v in d['actual_kappa9_C72']];powers=[]
 for idx in prog['core_kappa_indices4']:
  p=[one,kap[idx]]
  for _ in range(2,7):p.append(mult(p[-1],kap[idx]))
  assert not mult(p[-1],kap[idx]);powers.append(p)
 def poly(terms):
  z={}
  for exp,n,den in terms:
   if sum(exp)>=7:continue
   v=one
   for i,e in enumerate(exp):
    if e:v=mult(v,powers[i][e])
   z=add(z,scale(v,Q(n,den)))
  return z
 def inv(k,p):
  if not k:return one
  z={}
  for e in range(4):z=add(z,scale(powers[1][2*e],Q(math.comb(k+e-1,e),5**(k+e)) if p else math.comb(k+e-1,e)))
  return z
 G=[]
 for name,idx in zip(prog['five_pivot_order'],prog['five_pivot_kappa_indices'],strict=True):
  q=prog['five_actual_rational_roots'][name];G.append(add(kap[idx],mult(mult(poly(q['numerator']),inv(q['D_power'],False)),inv(q['P_power'],True))))
 G.extend(poly(v) for v in ci)
 assert [[[str(x) for x in g.get(i,zero)] for i in range(72)] for g in G]==nd['regular_G8_raw_positive_bounds8_by72']
 for i,row in enumerate(d['actual_G8_point_and_first40']):
  for j,x in enumerate(row):
   if Q(x):G[i][j]=ex(x)
   else:G[i].pop(j,None)
 sop=read('reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',['native_S12_by32'])['native_S12_by32'];sdefs={j:[] for j in range(28)}
 for row,col,g,sign in sop['koszul28_in_first8_rows']:sdefs[col].append((row,g,sign))
 N={(i,j):{k:ex(x) for k,x in enumerate(v) if Q(x)} for i,j,v in d['actual_needed_N_point_and_first40']};gamma=list(map(sp,d['Gamma4_source56_times8']))
 native_bounds={}
 for rr in rec:
  ss=[{} for _ in range(8)]
  for j,c in sp(rr['native128_constant_projection']).items():
   for row,g,_ in sdefs[j//4]:ss[row]=add(ss[row],mult({j%4:ex(c)},G[g]))
  for v in ss:
   for k in range(4):v.pop(k,None)
  xx={j:{} for j in [0,4,5,7]}
  for j in xx:
   for i in range(8):xx[j]=add(xx[j],mult(N[i,j],ss[i]))
  vals=[]
  for row in gamma:
   z=zero
   for slot,c in row.items():
    source,pair=divmod(slot,8);g,mask=divmod(source,4);z=plus(z,times(ex(c),xx[g].get(4*(10+pair)+mask,zero)))
   vals.append(z)
  native_bounds[rr['case619']]=vals
 assert [[[str(x) for x in z] for z in native_bounds[a['case619']]] for a in rec]==[a['four_native_partial_positive_bounds'] for a in nd['all70_native_partial_positive_bounds']]
 # Literal old first/second and carry sources are independently cache-bound.
 allfirst=read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json',['first640'])['first640']['all79_atom_values_full10'];firstmap={i:sp(v) for i,v in allfirst};first={i:sp(v) for i,v in d['actual_C0_leaf_first10']}
 assert all(firstmap[i]==v for i,v in first.items()) and all(not v.get(0,Q(0)) for v in first.values())
 rawvals=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',['all_requested164_root_values'])['all_requested164_root_values'];rawvals={i:sp(v) for i,v in rawvals};secs={i:list(map(Q,v)) for i,v in d['actual_13_honest_second8']}
 primitive=[27897,27948,28203,28254,28437,28692,28794]
 for i in primitive:
  vv=[rawvals[i].get(k,Q(0)) for k in copy];assert secs[i]==vv[10:]
  assert vv[:10]==[firstmap[i].get(k,Q(0)) for k in range(10)]
 assert all(not any(secs[i]) for i in set(secs)-set(primitive))
 assert read('reports/AMBIENT_B18_ELL2_T2_OLD_CASE508_ORDER_20261105_v63_REVIEW/report.json')['all_passed']
 basis=read('reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',['actual_C0_basis8_scalar_nodes'])['actual_C0_basis8_scalar_nodes']
 assert {7199,7454,7505,8045}<=set(basis[:7]) and basis[7]==505938
 # Rebuild the exact pruned old source DAG, not merely its cached values.
 oldsrc=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['old640_source_columns','old640_typed_closure_nodes','native_S1ell_columns','native_S0_columns','L0_Res4_Ntranspose_columns'])
 rows=sorted({k//8 for row in gamma for k in row});origB=list(map(dict,oldsrc['old640_source_columns']));origops=dict(oldsrc['old640_typed_closure_nodes']);closure=set()
 def visit(k):
  if k in closure:return
  closure.add(k);n=origops[k]
  for j in ([j for j,a,b in n[1]] if n[0]=='sum' else n[3] if n[0]=='product' else []):visit(j)
 oldids=sorted({j for a in rec for j,x in a['original768_constant_combination'] if j<640})
 want={j:{i:k for i,k in origB[j].items() if i in rows} for j in oldids}
 assert {j:dict(v) for j,v in d['oldB_source_rows_pruned']}==want
 for row in want.values():
  for k in row.values():visit(k)
 assert dict(d['typed_closure_nodes'])=={k:origops[k] for k in closure}
 del origB,origops
 # Recompute every positive old jet and every individual ell SUM.
 ops=dict(d['typed_closure_nodes']);jets={}
 def jm(a,b):
  z=[times(a[0],b[0])]+[plus(times(a[0],b[k]),times(a[k],b[0])) for k in range(1,10)]
  for k,(u,v) in enumerate(PAIR):
   t=plus(times(a[0],b[10+k]),times(a[10+k],b[0]));t=plus(t,times(a[u+1],b[v+1]))
   if u!=v:t=plus(t,times(a[v+1],b[u+1]))
   z.append(t)
  return z
 for k,n in sorted(ops.items()):
  if n[0]=='actual_C0_scalar_function':z=[ex(first[n[1]].get(i,Q(0))) for i in range(10)]+[ex(x) for x in secs.get(n[1],[Q(0)]*8)]
  elif n[0]=='c':z=[ex(Q(n[1],n[2]))]+[zero]*17
  elif n[0]=='sum':
   z=[zero]*18
   for child,a,b in n[1]:z=[plus(x,times(ex(Q(a,b)),y)) for x,y in zip(z,jets[child],strict=True)]
  else:
   assert n[0]=='product';z=[ex(Q(n[1],n[2]))]+[zero]*17
   for child in n[3]:z=jm(z,jets[child])
  jets[k]=z
 L0=list(map(sp,d['native_L0_rows_pruned']));Le=list(map(sp,d['native_L1ell_rows_pruned']));S0p={i:sp(v) for i,v in d['native_pivot_S0_columns']};Sep={i:sp(v) for i,v in d['native_pivot_S1ell_columns']}
 Xp={i:[sp(v) for v in a] for i,a in d['actual_pivot_old8_first_pruned']};ell={i:{k:ex(x) for k,x in sp(v).items()} for i,v in d['actual_old_pivot_ell_first_pruned']}
 assert L0==[{i:x for i,x in sp(v).items() if i in rows} for v in oldsrc['L0_Res4_Ntranspose_columns']]
 oldLe=read('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase3-actual-native-S2rell-L1ell-and619-coherent-parts.json',['actual_L1ell_columns'])['actual_L1ell_columns']
 assert Le==[{i:x for i,x in sp(v).items() if i in rows} for v in oldLe]
 piv=d['point_pivot27_original_source_ids'];assert piv==read('reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json',['original27_X_source_columns'])['original27_X_source_columns']
 for p,j in d['native_pivot_local128_ids']:assert j==piv[p]-640 and S0p[p]==sp(oldsrc['native_S0_columns'][j]) and Sep[p]==sp(oldsrc['native_S1ell_columns'][j])
 for offset,path,key in [(0,'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json','source_columns640_rows56_sparse_dual9'),(640,'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json','source_columns128_rows48_sparse_dual9')]:
  arr=read(path,[key])[key]
  for p in d['needed_pivot_indices']:
   if offset<=piv[p]<offset+len(arr):
    src={i:sp(v) for i,v in arr[piv[p]-offset]}
    for var in range(8):assert Xp[p][var]=={i:v[var+1] for i,v in src.items() if i in rows and v.get(var+1)}
  del arr
 oldell=read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json',['actual_first640_ell_columns'])['actual_first640_ell_columns']
 assert {p:sp(v) for p,v in d['actual_old_pivot_ell_first_pruned']}=={p:{i:x for i,x in sp(oldell[piv[p]]).items() if i in rows} for p in d['needed_pivot_indices'] if piv[p]<640}
 del oldsrc,oldLe,oldell
 for p in S0p:
  z={}
  for mat,v in [(L0,Sep[p]),(Le,S0p[p])]:
   for j,c in v.items():
    for k,x in mat[j].items():z[k]=plus(z.get(k,zero),times(ex(c),ex(x)))
  ell[p]=z
 B={j:dict(v) for j,v in d['oldB_source_rows_pruned']};bounds=[]
 for rr in rec:
  old={i:x for i,x in sp(rr['original768_constant_combination']).items() if i<640};a=sp(rr['actual_a1_27_times9']);four=[]
  for chan,row in enumerate(gamma):
   ob=cb=zero
   for slot,c in row.items():
    source,pair=divmod(slot,8);u,v=PAIR[pair]
    for j,c0 in old.items():
     node=B[j].get(source)
     if node is not None:ob=plus(ob,times(ex(c),times(ex(c0),jets[node][10+pair])))
    for p in d['needed_pivot_indices']:
     au,av=a.get(9*p+u,Q(0)),a.get(9*p+v,Q(0))
     if au:cb=plus(cb,times(ex(c),times(ex(au),ex(Xp[p][v].get(source,Q(0))) if v<8 else ell[p].get(source,zero))))
     if u!=v and av:cb=plus(cb,times(ex(c),times(ex(av),ex(Xp[p][u].get(source,Q(0))) if u<8 else ell[p].get(source,zero))))
   nbnd=native_bounds[rr['case619']][chan];total=plus(plus(ob,nbnd),cb);four.append({k:list(map(str,val)) for k,val in [('old',ob),('native',nbnd),('individual_pivot_carry',cb),('combined',total)]})
  bounds.append({'case619':rr['case619'],'four_component_positive_envelopes':four})
 assert bounds==d['all70_four_component_positive_envelopes'] and maxbits<=128 and 3*maxbits+1<=385
 # Immutable prepared worker / live addendum are mandatory final premises.
 cp=read(custody);lp=read(live);assert cp['all_passed'] and lp['all_passed']
 assert cp['worker_sha256']=='cfcba5133a81026bc782ba85e8b40ccf68e58a38908e4b0a3cae4e152b000a47' and cp['source_gate_sha256']==pins[BASE+'report.json'] and cp['safe_help_literal_PASS']
 for p,h in cp['worker_and_helper_sha256'].items():assert sha(ROOT/p)==h
 for p,h in cp['literal_external_source_sha256'].items():assert sha(ROOT/p)==h
 help_run=subprocess.run([sys.executable,str(ROOT/cp['worker']),'--help'],cwd=ROOT,capture_output=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},timeout=10)
 assert help_run.returncode==cp['safe_help_exit']==0 and not help_run.stderr and hashlib.sha256(help_run.stdout).hexdigest()==cp['safe_help_stdout_sha256']
 lev=read(str(Path(live).parent/'explicit-final-input-and-worker-overlap-count.json'))
 assert count(d)==lev['frozen_input_numeric_including_all_final_envelope_fields']==15878
 assert lev['two_decoded_input_copies_conservative']==2*count(d)==31756
 assert sum(lev['parsed_helper_categories'].values())==lev['helper_actual_count_upper']<=lev['helper_reserved']==6000
 assert sum(lev['signed_vector_working_categories'].values())==8078
 assert 31756+6000+8078+lev['all_output_and_serialization_reserve']+lev['gate_peer_custody_metadata_reserve']==lev['actual_worker_conservative_owned_bound']==49734<lev['existing_cap_unchanged']==51263
 assert lp['source_gate_sha256']==pins[BASE+'report.json'] and lp['source_input_sha256']==pins[BASE+'complete70-cached-input-source-and-positive-whole-image-bounds.json']
 assert cp['live_count_addendum_sha256']==pins[live] and cp['explicit_worker_owned_bound']==49734 and cp['unchanged_whole_numeric_cap']==51263
 # The separately retained C16 overlap array fits even if it is conservatively
 # added to every listed helper category:5509+128=5637<6000.
 assert lev['helper_actual_count_upper']+8*16<lev['helper_reserved']
 owned=count([d,frame,kernel,nd,raw,aa,Pi,S0,kap,powers,G,N,first,rawvals,jets,Xp,ell,bounds])+14000
 assert owned<150000
 ep=out/'independent70-source-frame-positive-bound-and-worker-proof.json'
 evidence={'all70_fullC4_S0_kernels':True,'all1089_raw33_A18_products_and297_intertwiners':True,'actual_canonical_kappa_monic_and_Nfirst_cache_bindings':True,'regular_G8_literal_program_not_radical9':True,'all480_C72_products':True,'all70_native_positive_envelopes':True,'all70_combined_positive_envelopes':True,'individual_ell_pivot_is_SUM':True,'no_signed_G8_second_or_native_or_combined_images_evaluated':True,'actual_positive_max_bits':maxbits,'source_owned_bound':owned,'custody':{'path':custody,'sha256':pins[custody]},'live':{'path':live,'sha256':pins[live]}}
 ep.write_text(json.dumps(evidence,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'source_only_no70_response_values':True,'source_positive_bound_verified':[128,385],'input_sha256':pins,'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start,'source_owned_bound':owned}
 (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'seconds':report['seconds'],'owned':owned,'positive_bits':maxbits}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--custody',required=True);p.add_argument('--live',required=True);a=p.parse_args();main(a.output,a.custody,a.live)

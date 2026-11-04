#!/usr/bin/env python3
"""Independent fixed-source adjoint/jet proof for the actual T quadratic layer.
Standard-library Fraction only; no producer imports, source solve or rank search.
"""
import json,hashlib,itertools,signal,resource,time,sys
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63_REVIEW'
OUT=ROOT/'reports'/STEM
PATHS={
 'V56':'AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
 'Minimal':'AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/actual-minimal5-by825-Schur-operator.json',
 'M1':'AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/actual-projected-three-class-X2-operator-preflight.json',
 'W':'AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63/actual-linear-cycle-muT-certificate.json',
 'Q':'AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/actual-linear-quadratic-T-certificate.json',
 'L1':'AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
 'L2':'AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase2-certified-quadratic-input.json',
 'Dual':'AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/three-quadratic-duals.json',
 'Comb':'AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/zero-linear620-source-combinations.json',
 'Sens':'AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/projected-source-X2-sensitivity.json',
 'B':'AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'AF':'AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
 'C0':'AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'Full':'AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json',
 'Source':'AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'J':'AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'Cache':'AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63/actual-source105-Artin20-values.json',
 'H20':'AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'Fold':'AMBIENT_B18_T_NEXT3_ROOT_ATOM_ADJOINT_20261104_v63/phase1-folded-root-functionals.json',
 'Adj':'AMBIENT_B18_T_NEXT3_ROOT_ATOM_ADJOINT_20261104_v63/phase2-atom-second-adjoint.json',
 'AD':'AMBIENT_B18_T_NEXT3_ROOT_ATOM_ADJOINT_20261104_v63/oldB-AD-fixed-node-values-and-gradients.json',
 'Copy':'AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase2-existing-cache-B11-copy-values.json',
 'Comp':'AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
 'Native':'AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json',
 'Result':'AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63/actual-three-by620-quadratic-image-certificate.json',
 'Xfirst':'AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
 'Nfirst':'AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json',
 'X0':'AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json'}
GATES={
 'AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63_REVIEW':'998bcf929242837dee44de326f05950b1b295a8828f9cfc8785c082503804293',
 'AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW':'ca127c4dc0f3c12a4d4d80dfc1ca316f42dee8b90ecc74cb3e5e423752c906b7',
 'AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63_REVIEW':'c933016430a1779bed947b0c81484b0aae23afbee23e4df15a386ea959bf0dd9',
 'AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63_REVIEW':'bb6e7e7c5327528f2be31b1d29db5891461d12adb319daa184b6c6addd63d693'}
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def put(d,k,c):
 c=F(c)+d.get(k,F(0))
 if c:d[k]=c
 else:d.pop(k,None)
def add(d,e,c=F(1)):
 for k,v in e.items():put(d,k,c*v)
 return d
def dec(xs):return {tuple(k) if isinstance(k,list) else k:F(v) for k,v in xs if F(v)}
def enc(d):return [[list(k) if isinstance(k,tuple) else k,str(v)] for k,v in sorted(d.items())]
def dump(name,x):(OUT/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def load(p):return json.loads((ROOT/'reports'/p).read_text())
def dot(a,b):return sum((v*b.get(k,F(0)) for k,v in a.items()),F(0))
def closure(nodes,roots,stops=(),typed=False):
 seen=set();todo=list(roots);stops=set(stops)
 while todo:
  i=todo.pop()
  if i in seen:continue
  seen.add(i)
  if i in stops:continue
  n=nodes[i];o=n[0]
  if typed:
   if o=='sum':todo.extend(j for j,a,b in n[1])
   elif o=='product':todo.extend(n[3])
   else:assert o in ['c','actual_C0_scalar_function'],(i,n)
  elif o in ['a','m']:todo.extend(n[1:3])
  elif o in ['n','i']:todo.append(n[1])
  else:assert o=='c',(i,n)
 return sorted(seen)
class Quotient:
 def __init__(self,mons):
  self.mons=mons;index={m:i for i,m in enumerate(mons)}
  self.mi={(i,j):index.get(tuple(a+b for a,b in zip(m,n))) for i,m in enumerate(mons) for j,n in enumerate(mons)}
 def mul(self,a,b):
  out={}
  for i,c in a.items():
   for j,d in b.items():
    k=self.mi[i,j]
    if k is not None:put(out,k,c*d)
  return out
 def inv(self,a):
  c=a.get(0,F(0));assert c
  h={i:v for i,v in a.items() if i};out={0:1/c};add(out,h,-1/c**2);add(out,self.mul(h,h),1/c**3)
  assert self.mul(a,out)=={0:F(1)}
  return out
 def evaluate(self,nodes,roots,cache,first=None):
  need=closure(nodes,roots,set(cache)|set(first or {}));values={};grads={};inverses=0
  for i in need:
   n=nodes[i];o=n[0]
   if i in cache:v=cache[i];g={}
   elif first is not None and i in first:v=first[i];g={i:F(1)}
   elif o=='c':v={0:F(n[1],n[2])} if n[1] else {};g={}
   elif o=='a':v=add(dict(values[n[1]]),values[n[2]]);g=add(dict(grads[n[1]]),grads[n[2]])
   elif o=='n':v={j:-c for j,c in values[n[1]].items()};g={j:-c for j,c in grads[n[1]].items()}
   elif o=='m':
    v=self.mul(values[n[1]],values[n[2]]);g={};add(g,grads[n[1]],values[n[2]].get(0,F(0)));add(g,grads[n[2]],values[n[1]].get(0,F(0)))
   elif o=='i':v=self.inv(values[n[1]]);g={j:-c/values[n[1]][0]**2 for j,c in grads[n[1]].items()};inverses+=1
   else:raise AssertionError((i,n))
   values[i]=v;grads[i]=g
  assert sum(len(values[i])+len(grads[i]) for i in need)<150000
  return values,grads,len(need),inverses

def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));sys.set_int_max_str_digits(0)
 OUT.mkdir(exist_ok=True);hashes={str(Path('reports')/p):sha(ROOT/'reports'/p) for p in PATHS.values()}
 dump('phase.json',{'stage':'standalone_input_hashes_before_arithmetic','input_sha256':hashes})
 D={k:load(p) for k,p in PATHS.items()}
 gates={}
 for stem,want in GATES.items():
  p=ROOT/'reports'/stem/'report.json';assert sha(p)==want,(stem,sha(p));assert json.loads(p.read_text())['all_passed'];gates[str(p.relative_to(ROOT))]=want
 M,W,Q,L1,L2=D['M1'],D['W'],D['Q'],D['L1'],D['L2'];mons=list(itertools.combinations_with_replacement(range(8),2));mi={m:i for i,m in enumerate(mons)}
 # Backward dualization of fixed normalized columns (no basis search).
 qb=list(map(dec,Q['quadratic_certificate']['basis_normalized_columns']));piv=[min(b) for b in qb];assert len(qb)==177 and len(set(piv))==177
 reps=[i for i in range(180) if i not in piv];assert reps==[143,176,179]
 dual=[]
 for r in reps:
  x={r:F(1)}
  for b,p in sorted(zip(qb,piv),key=lambda z:z[1],reverse=True):put(x,p,-sum((v*x.get(k,F(0)) for k,v in b.items() if k!=p),F(0)))
  assert all(dot(x,b)==0 for b in qb);dual.append(x)
 assert dual==list(map(dec,D['Dual']['projection3_by180_rows_sparse']))==list(map(dec,M['Pi3_180_dual_rows']))
 assert all(dot(x,dec(c))==0 for x in dual for c in L2['all464_certified_quadratic_columns_sparse'])
 # Adjoint Euler construction of the full3x56x36 functional.
 wb=list(map(dec,W['complete_rank_certificate']['basis_normalized_columns']));wp=[min(b) for b in wb];assert len(wb)==84
 coordinates=[]
 for j in range(84):
  x={j:F(1)};c={}
  for p in sorted(range(84),key=lambda p:wp[p]):
   b=wb[p];a=x.get(wp[p],F(0))
   if a:c[p]=a;add(x,b,-a)
  assert not x;coordinates.append(c)
 qnorm=list(map(dec,L1['normalized84_unitcycle_Q1_columns_sparse']));assert len(qnorm)==84
 pi=[{M['free_kernel_rows32'][i]:v for i,v in dec(row).items()} for row in M['pi_free32_rows']];m0=[[F(v) for v in row] for row in D['V56']['all24_by56_m0']]
 assert m0==[[F(v) for v in row] for row in M['m0_actual24_by56']] and M['pi_free32_rows']==D['Minimal']['constant_free32_quotient_projection5_rows']
 assert D['Minimal']['actual_minimal_generator_presentation_shape']==[5,825]
 sens=[]
 for t,dt in enumerate(dual):
  beta=[]
  for q in qnorm:
   beta.append({k:-sum((b*dt.get(36*(i//8)+mi[tuple(sorted((i%8,k)))],F(0)) for i,b in q.items()),F(0)) for k in range(8)})
  s={}
  for target,row in enumerate(pi):
   for src,c in row.items():
    for m in range(36):put(s,36*src+m,c*dt.get(36*target+m,F(0)))
  for w,(sector,(i,j)) in enumerate(W['literal_wedge_basis84']):
   for v in range(8):
    c=sum((a*beta[p][v] for p,a in coordinates[w].items()),F(0))
    if not c:continue
    terms=[(i,j,1),(j,i,-1)] if j<7 else [(i,j,1)]
    for derivative,outrow,sign in terms:
     k=int(derivative<7)+int(v<7);den=k+1 if j<7 else k
     assert den
     factor=c*sign*F(2 if derivative==v else 1,den);m=mi[tuple(sorted((derivative,v)))]
     for src,a in enumerate(m0[8*sector+outrow]):
      if a:put(s,36*src+m,factor*a)
  sens.append(s)
 assert sens==list(map(dec,D['Sens']['three_rows_sparse']));assert list(map(len,sens))==[2,28,1]
 dump('phase.json',{'stage':'own_adjoint31_sensitivity_and_known464_passed'})
 # Independently reconstruct source combinations, constant/first cancellation and27 pivots.
 def fixed_zero_combinations(cert):
  normalized=[];zeros=[]
  for j,op in enumerate(cert['fixed_elimination_receipts']):
   assert j==op['column'];c={j:F(1)}
   for k,a in op['subtract']:add(c,normalized[k],-F(a))
   if 'pivot_row' in op:normalized.append({i:a/F(op['divide']) for i,a in c.items()})
   else:assert op['zero'];zeros.append(c)
  return zeros
 def compose(c,columns):
  out={}
  for i,a in c.items():add(out,columns[i],a)
  return out
 wz=fixed_zero_combinations(W['complete_rank_certificate']);lz=fixed_zero_combinations(Q['complete_linear_rank_certificate'])
 base=list(map(dec,W['actual_constant_source_combinations']))
 comb=D['Comb']['all620'];assert len(comb)==620
 for row,c657 in zip(comb,lz):
  c741=compose(c657,wz);c768=compose(c741,base)
  assert c657==dec(row['combination_in657']) and c741==dec(row['combination_in741']) and c768==dec(row['combination_in_actual768'])
 assert len(wz)==657 and len(lz)==620
 assert [x['combination_in_actual768'] for x in comb]==M['620_actual_source_combinations_in_original768']
 assert [x['linear_Xp27_correction_row8j_plus_variable'] for x in comb]==M['620_actual_linear_pivot_coefficients_a1_by27_times8']
 first=D['Xfirst']['source_columns640_rows56_sparse_dual9']+D['Nfirst']['source_columns128_rows48_sparse_dual9']
 first=[{int(r):dec(v) for r,v in col} for col in first]
 xp1=list(map(dec,M['actual_pivot_X1_by27_sparse56_times8']));xp0=list(map(dec,M['actual_pivot_X0_columns27_sparse56']))
 for j,c0 in enumerate(D['Comb']['normalized_Xp27_original768_combinations']):
  a={};b={}
  for src,c in dec(c0).items():
   for r,p in first[src].items():
    put(a,r,c*p.get(0,F(0)))
    for k,v in p.items():
     if k:put(b,8*r+k-1,c*v)
  assert a==xp0[j] and b==xp1[j]
 for row in comb:
  origin={};linear={};correct={}
  for src,c in dec(row['combination_in_actual768']).items():
   for r,p in first[src].items():
    put(origin,r,c*p.get(0,F(0)))
    for k,v in p.items():
     if k:put(linear,8*r+k-1,c*v)
  for a,c in dec(row['linear_Xp27_correction_row8j_plus_variable']).items():
   p,k=divmod(a,8)
   for r,v in xp0[p].items():put(correct,8*r+k,c*v)
  assert not origin and linear==correct
 # Fold actual source roots and retain every P component and native effective callback.
 B,C0,AF=D['B'],D['C0'],D['AF'];ops=B['operator_scalar_nodes'];cnodes=C0['DAG_scalar_nodes'];bcols=list(map(dict,B['B48_by640_sparse_columns']))
 pairs=sorted({mons[i%36] for s in sens for i in s});assert pairs==[(0,7),(1,7),(4,7),(5,7),(7,7)];pairid={v:i for i,v in enumerate(pairs)}
 folded=[];native=[];known=[];roots=set()
 for num,row in enumerate(comb):
  cs=dec(row['combination_in_actual768']);alpha=dec(row['linear_Xp27_correction_row8j_plus_variable']);fs=[];ns=[];ks=[]
  for s in sens:
   f={};n={};k=F(0)
   for slot,a in s.items():
    src,m=divmod(slot,36);u,v=mons[m];p=pairid[u,v]
    for j,x in enumerate(xp1):
     k-=a*alpha.get(8*j+u,F(0))*x.get(8*src+v,F(0))
     if u!=v:k-=a*alpha.get(8*j+v,F(0))*x.get(8*src+u,F(0))
    for j,b in cs.items():
     if j>=640:
      if src<48:put(n,(j-640,src,p),a*b)
     elif src<48:
      r=bcols[j].get(src)
      if r is None or ops[r][0]=='c':continue
      if ops[r][0]=='actual_C0_scalar_function':key=(1,ops[r][1],p)
      else:key=(0,r,p);roots.add(r)
      put(f,key,a*b)
     else:
      r=C0['honest_D0_relation_columns'][j][src-47]
      if cnodes[r][0]!='c':put(f,(1,r,p),a*b)
   fs.append(f);ns.append(n);ks.append(k)
  saved=D['Fold']['620_rows'][num]
  assert fs==list(map(dec,saved['three_oldB_P_root_functionals'])) and ns==list(map(dec,saved['three_native_X2_functionals'])) and ks==list(map(F,saved['three_known_Xp1alpha1_corrections']))
  folded.append(fs);native.append(ns);known.append(ks)
 # A14 quotient multiplication simultaneously obtains origin,firsts,selectedcross terms.
 zero=(0,)*8;units=[tuple(int(j==i) for j in range(8)) for i in range(8)];qmon=[tuple(int(i==u)+int(i==v) for i in range(8)) for u,v in pairs];A14=Quotient([zero]+units+qmon)
 atom={i:dec(p) for i,p in AF['all79_atom_values']};needed=closure(ops,roots,typed=True);vals={};grads={}
 for i in needed:
  n=ops[i];o=n[0]
  if o=='c':v={0:F(n[1],n[2])} if n[1] else {};g={}
  elif o=='actual_C0_scalar_function':v=atom[n[1]];g={n[1]:F(1)}
  elif o=='sum':
   v={};g={}
   for j,a,b in n[1]:add(v,vals[j],F(a,b));add(g,grads[j],F(a,b))
  elif o=='product':
   v={0:F(n[1],n[2])};g={}
   for j in n[3]:
    gg={};add(gg,g,vals[j].get(0,F(0)));add(gg,grads[j],v.get(0,F(0)));v=A14.mul(v,vals[j]);g=gg
  else:raise AssertionError(n)
  vals[i]=v;grads[i]=g
 assert len(needed)==762
 for i,c,g,f,q in D['AD']['node_data']:
  assert vals[i].get(0,F(0))==F(c) and grads[i]==dec(g)
  assert {k-1:v for k,v in vals[i].items() if 1<=k<=8}==dec(f)
  assert {k-9:v for k,v in vals[i].items() if k>=9}==dec(q)
 atomrows=[];constant=[]
 for num,fs in enumerate(folded):
  ars=[];ks=[]
  for t,f in enumerate(fs):
   a={};k=known[num][t]
   for (kind,i,p),c in f.items():
    if kind:put(a,(i,p),c)
    else:
     for j,v in grads[i].items():put(a,(j,p),c*v)
     k+=c*vals[i].get(9+p,F(0))
   ars.append(a);ks.append(k)
  saved=D['Adj']['620_rows'][num]
  assert ars==list(map(dec,saved['three_atom_or_opaque_P_second_functionals'])) and ks==list(map(F,saved['known_first_product_and_Xp_correction_three_vector']))
  atomrows.append(ars);constant.append(ks)
 dump('phase.json',{'stage':'all620_root_fold_and_762_node_independent_jet_AD_passed'})
 # Full B11 coupling. Copy syntax gate is independently pinned; all values/equations are fresh here.
 Full,Source,J,Cache,H=D['Full'],D['Source'],D['J'],D['Cache'],D['H20'];inst=Full['all16_pointed_unit57_inverse_instances'][14];CB=J['both57_source2']['CB'];bind=dict(inst['bindings']);eq=inst['all57_pointed_equations']
 assert inst['retained_Y_indices']==[7,8] and set(bind)==set(CB['unknown57'])
 assert CB['equations57']==Full['original_unit48_equations_CB']+Full['CB_coordinate_Y9_nodes']
 assert cnodes[:len(Full['DAG_scalar_nodes'])]==Full['DAG_scalar_nodes'] and cnodes[:len(Source['DAG_scalar_nodes'])]==Source['DAG_scalar_nodes']
 assert cnodes[11]==['v','CA_trace_p'] and Source['families']['CA']['p_quadratic_relation']['s']==11 and 11 in Cache['joint105_unknown_order']
 bmons=[zero,units[7],tuple(2*int(i==7) for i in range(8))]+[units[k] for k in [0,1,4,5]]+qmon[:4];B11=Quotient(bmons);bi={m:i for i,m in enumerate(bmons)};variables=[{bi[m]:F(1)} if m in bi else {} for m in units]
 images=[{bi[tuple(m)]:F(1)} if tuple(m) in bi else {} for m in H['quotient_basis_monomials']]
 def phi(p):
  v={}
  for i,a in p.items():add(v,images[i],a)
  return v
 for k,action in enumerate(H['multiplication8_column_sparse']):
  for j,p in enumerate(action):assert B11.mul(variables[k],images[j])==phi(dec(p))
 assert Cache['A20_basis']==H['quotient_basis_monomials']
 cached={i:phi(dec(p)) for i,p in Cache['all_requested_original_source_root_values']}
 cached.update({i:phi(dec(p)) for i,p in zip(Cache['joint105_unknown_order'],Cache['all105_unknown_final_values'])})
 yroots=CB['F9_nodes'][7:9];yb,_,yn,_=B11.evaluate(cnodes,yroots,cached);args=[{}]*7+[yb[i] for i in yroots];copy={}
 for original,p in CB['all57_source2jets']:
  v={}
  for m,c in zip(J['source2_monomials'],p):
   if not F(c):continue
   term={0:F(c)}
   for k in m:term=B11.mul(term,args[k])
   add(v,term)
  copy[bind[original]]=v
 assert copy=={i:dec(p) for i,p in D['Copy']['all57_copied_unknown_values_through2']}
 eqvals,_,en,_=B11.evaluate(cnodes,eq,cached|copy);assert all(not eqvals[i] for i in eq)
 unknown={}
 for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
  for i,p in J[key].items():unknown[int(i)]=list(map(F,p[:8]))+[F(p[8])+F(p[9])]
 for row in AF['all35_ordered_CA_prefix_firsts']:
  p=dec(row['first9_bar']);unknown[row['node']]=[p.get(k,F(0)) for k in range(9)]
 for i,p in AF['all8_g_quotient_overrides']:
  p=dec(p);unknown[i]=[p.get(k,F(0)) for k in range(9)]
 cache=cached|copy
 assert C0['cached_actual_CA_gA_all_degree3_coefficients'][:6]==['0']*6
 for k,i in enumerate(C0['actual_C0_basis8_scalar_nodes']):cache[i]=variables[k] if k<7 else {}
 first={}
 for i,p in unknown.items():
  v={0:p[0]} if p[0] else {}
  for k in range(8):add(v,variables[k],p[k+1])
  first[i]=v
 active=sorted({i for ars in atomrows for a in ars for i,p in a});assert len(active)==108
 cv,cg,cl,iv=B11.evaluate(cnodes,active,cache,first)
 for i,v,g in D['Comp']['all108_composite_values_and_opaque_second_gradients']:assert cv[i]==dec(v) and cg[i]==dec(g)
 dump('phase.json',{'stage':'full57_copy_equations_H20_B11_and108_symbolic_composites_passed','composite_nodes':cl})
 # All opaque second symbols cancel BEFORE effective native second assembly.
 E=D['Native'];eff=[{r:dec(p) for r,p in col} for col in E['source_columns128_rows48_sparse_B11_effective_second']];raw=[{r:dec(p) for r,p in col} for col in E['actual_S_source_columns128_rows48_sparse_B11']];ptob=[7,8,9,10,2];actual=[];splits=[]
 assert E['individual_actual_native_X2_claimed'] is False and E['valid_for_actual_zero_origin_source_combinations_only_by_L2_cancellation']
 for number,ars in enumerate(atomrows):
  raw0={}
  for j,c in dec(comb[number]['combination_in_actual768']).items():
   if j>=640:
    for r,p in raw[j-640].items():put(raw0,r,c*p.get(0,F(0)))
  assert not raw0
  old=[];nv=[]
  for t,a in enumerate(ars):
   value=constant[number][t];opaque={}
   for (i,p),c in a.items():
    value+=c*cv[i].get(ptob[p],F(0))
    for j,v in cg[i].items():put(opaque,(j,p),c*v)
   assert not opaque,('UNKNOWN_SECOND_NOT_CANCELLED',number,t,opaque)
   old.append(value);nv.append(sum((c*eff[j].get(r,{}).get(ptob[p],F(0)) for (j,r,p),c in native[number][t].items()),F(0)))
  v=[a+b for a,b in zip(old,nv)];s=D['Result']['all620_source_value_splits'][number]
  assert old==list(map(F,s['old_actual_and_pivot_correction'])) and nv==list(map(F,s['native_effective_second_contribution'])) and v==list(map(F,s['actual_projected_quadratic_vector']))
  assert {i:x for i,x in enumerate(v) if x}==dec(D['Result']['actual620_three_class_columns'][number])
  actual.append(v);splits.append({'index':number,'old':list(map(str,old)),'native':list(map(str,nv)),'full':list(map(str,v))})
 nz=[i for i,v in enumerate(actual) if any(v)];assert len(nz)==6 and actual[563]==[F(-250,3),F(-1470875,5478),F(0)]
 assert all(v[1]==F(11767,3652)*v[0] and not v[2] for v in actual)
 dump('independent-fixed-source-certificate.json',{'own_dual3':list(map(enc,dual)),'own_sensitivity31':list(map(enc,sens)),'source_combinations':620,'copy57_equations_full_B11_zero':True,'all160_H20_actions':True,'old_AD_nodes':len(needed),'composite_nodes':cl,'unit_inverse_checks':iv,'all620_opaque_second_maps_zero':True,'all620_native_source_origin_maps_zero':True,'all620_three_class_splits':splits,'nonzero_column_indices':nz,'fixed_selected_column563':list(map(str,actual[563])),'common_normalized_direction':['1','11767/3652','0'],'rank':1,'HF_through2':[5,3,2],'length_mod_m3':10})
 report={'id':STEM,'all_passed':True,'independent_standard_library_Fraction_only':True,'producer_imported_or_rerun':False,'new_source_jet_solve_or_rank_search':False,'source31_adjoint_independently_reconstructed':True,'fixed620_old_AD_and_full_B11_couplings':True,'all_opaque_second_coefficients_cancel_before_numeric_assembly':True,'native_L2_cancel_checked_on_all620_raw_source_origin_maps':True,'actual_projected_quadratic_rank':1,'actual_HF_through_degree2':[5,3,2],'actual_length_T_mod_m3':10,'total_T_length_or_full_nil_cube_not_computed':True,'input_sha256':hashes,'independent_dependency_reports_sha256':gates,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{'independent-fixed-source-certificate.json':sha(OUT/'independent-fixed-source-certificate.json')},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 assert sum(p.stat().st_size for p in OUT.iterdir())<33554432
 dump('report.json',report);dump('phase.json',{'stage':'FINAL','all_passed':True});print(json.dumps(report))
if __name__=='__main__':
 try:main()
 except BaseException as e:
  dump('failure.json',{'exception':repr(e),'code_sha256':sha(Path(__file__))});raise

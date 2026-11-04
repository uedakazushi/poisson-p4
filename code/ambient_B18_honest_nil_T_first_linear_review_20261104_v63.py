#!/usr/bin/env python3
"""Independent cached source2/monics/m1 and all849 Tnil first relations.
No producer imports, source solve, missing callback assignment, or rank search.
"""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63';ID=BASE+'_REVIEW';OUT=ROOT/'reports'/ID
M=[()]+[(i,)for i in range(8)]+list(itertools.combinations_with_replacement(range(8),2));INDEX={m:i for i,m in enumerate(M)}
EXP=[]
for mon in M:
 e=[0]*8
 for j in mon:e[j]+=1
 EXP.append(e)
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(tuple,list)):return [enc(v)for v in x]
 return x
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={};cache={}
 def load(st,name='report.json'):
  p=ROOT/'reports'/st/name;k=str(p.relative_to(ROOT));pins[k]=sha(p)
  if k not in cache:cache[k]=json.loads(p.read_text())
  return cache[k]
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),indent=2)+'\n')
 report=load(BASE);assert pins[f'reports/{BASE}/report.json']=='8bbfd553d06781f8f3277865c2518fc9709622bf8a6cdfb58eaab49c605829b8'
 for k,h in report['input_sha256'].items():assert sha(ROOT/k)==h
 for k,h in report['evidence_sha256'].items():assert sha(ROOT/'reports'/BASE/k)==h
 for st in ['AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES','AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR','AMBIENT_B18_T3_ZERO_STRUCTURAL','AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20']:
  assert load(st+'_20261104_v63_REVIEW')['all_passed']
 p1=load(BASE,'phase1-cached-source-monics-and-physical-adapter.json');p2=load(BASE,'phase2-row20-surrogates-and-true-m1row7.json');p3=load(BASE,'phase3-all849-minimal-first-relations-and-rank.json')
 completion=load('AMBIENT_B18_HONEST_NIL_T_LINEAR_CACHED_COMPLETION_PREFLIGHT_20261104_v63','actual-final-row20-cancellations-and-eleven-node-monic-source.json')
 whole=load('AMBIENT_B18_HONEST_NIL_T_LINEAR_SOURCE_WHOLE_PREFLIGHT_20261104_v63','actual-final681-source-sensitivities-and-cached-directions.json')
 minimal=load('AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63','actual-minimal-Tnil5-849-and-d3-linear-operator.json')
 m0data=load('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase1-actual-m0-and-section.json');X=load('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase2-actual-radical-X-firsts.json');W=load('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase3-actual-W84-rank-and-mu.json')
 oldB=load('AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63','actual-first640-X-mod-m2.json')
 G=load('AMBIENT_B18_NATIVE128_FIRST_20261104_v63','phase5-native-adapter-firsts.json');proper=load('AMBIENT_B18_NATIVE128_FIRST_20261104_v63','phase4-native-source-and-R12.json');nfirst=load('AMBIENT_B18_NATIVE128_FIRST_20261104_v63','phase1-source-firsts.json')
 B11=load('AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63','phase3-actual-proper-G8-B11.json');nb=load('AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63','phase1-source-B11-values.json');bmons=load('AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63','actual-B11-source-and-gauge-preflight.json')['B11_basis_monomials8']
 source2=load('AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63','full9-firstjet-certificate.json');adapter=load('AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63','actual-special-fiber20-coefficients.json');R=load('AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63','actual-radical-normalizer-comparison-operator.json');GL=load('AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63','actual-E4-nil-two-generators-and-relations.json')
 save('phase.json',dict(phase='pinned source caches and independent baseline peers before arithmetic',no_new_source=True))
 # Independent commutative degree2 substitution, with physical ell=s-r set0.
 def bmul(p,q):
  v={}
  for i,a in p.items():
   for j,b in q.items():
    mon=tuple(sorted(M[i]+M[j]))
    if len(mon)<=2:Q.put(v,INDEX[mon],a*b)
  return v
 def fold_mons(v,mons):
  p={}
  for mon,c in zip(mons,v):
   if F(c):Q.put(p,INDEX[tuple(sorted(7 if i==8 else i for i in mon))],F(c))
  return p
 def arr(p,n):return [p.get(i,F(0))for i in range(n)]
 Y=[]
 for lin,quad in zip(adapter['YB_in_YA_linear'],adapter['YB_in_YA_quadratic']):
  v={}
  for i,c in enumerate(lin):
   if F(c):Q.put(v,INDEX[(7 if i==8 else i,)],F(c))
  for e,c in quad:
   if F(c):Q.put(v,INDEX[tuple(sorted(7 if i==8 else i for i,n in enumerate(e)for _ in range(n)))],F(c))
  Y.append(v)
 own={f:dict(source2['both57_source2'][f]['all57_source2jets'])for f in ['CA','CB']}
 def ca(i):return fold_mons(own['CA'][i],source2['source2_monomials'])
 def cb(i):
  v={}
  for mon,c in zip(source2['source2_monomials'],own['CB'][i]):
   if not F(c):continue
   q={0:F(c)}
   for j in mon:q=bmul(q,Y[j])
   v=Q.plus(v,q)
  return v
 assert [arr(v,45)for v in Y]==[[F(c)for c in v]for v in p1['CB_F9_in_CA_F9_ell0_degree2']]
 shared=[]
 for i in range(2,11):assert ca(i)==cb(i);shared.append(ca(i))
 assert [arr(v,45)for v in shared]==[[F(c)for c in v]for v in p1['all9_shared_physical_degree2_polynomials']]
 trace=[ca(11),ca(12),cb(20537),cb(20538)];nu=[Q.plus(trace[1],Q.scaled(bmul(trace[0],trace[0]),F(1,4))),Q.plus(trace[3],Q.scaled(bmul(trace[2],trace[2]),F(1,4)))]
 assert [arr(v,45)for v in trace]==[[F(c)for c in v]for v in p1['four_literal_source_s_t_degree2']]
 assert [arr(v,45)for v in nu]==[[F(c)for c in v]for v in p1['actual_nu_A_B_degree2']]
 assert [[v.get(i,F(0))for i in range(9)]for v in nu]==[[F(c)for c in v]for v in nfirst['nu_first9']]
 BMON=[(),(7,),(7,7),(0,),(1,),(4,),(5,),(0,7),(1,7),(4,7),(5,7)]
 assert [[v.get(INDEX[m],F(0))for m in BMON]for v in nu]==[[F(c)for c in v]for v in nb['nu_A_B_B11']]
 src=load('AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63','actual-source-coefficients.json');nodes=src['DAG_scalar_nodes'];roots=[src['families'][f]['actual_tracefree_nu']for f in ['CA','CB']];todo=roots[:];seen=set()
 while todo:
  i=todo.pop()
  if i not in seen:
   seen.add(i)
   if nodes[i][0]in['a','m','n','i']:todo.extend(nodes[i][1:])
 assert len(seen)==11 and sorted(i for i in seen if nodes[i][0]=='v')==[11,12,20537,20538]
 assert [[i,nodes[i]]for i in sorted(seen)]==completion['actual_monic2_cached_source_DAG']
 del src,nodes
 C180=Q.FiniteC(EXP,[arr(v,45)for v in nu]);C36=Q.FiniteC(EXP[:9],nfirst['nu_first9']);C44=Q.FiniteC(bmons,nb['nu_A_B_B11'])
 def variables(C,kappa,u):
  core=[Q.sparse(kappa[i])for i in[0,1,2,6]]
  return [Q.sparse(v)for v in u[:5]]+[C.polynomial(GL['forward_GL4'][s],core)for s in ['h','c','r','l']]
 va=variables(C180,proper['kappa9'],G['G8']);vb=variables(C44,B11['kappa9_C_B44'],B11['actualG8_C_B44'])
 native=[]
 for col in R['native_tensor_Koszul6_HB2_resolution'][1]['columns']:
  pol=next((p for row,p in col if row==5),[]);assert all(not e[5]and sum(e)<=2 for e,c in pol)
  v=C180.polynomial(pol,va)
  for k in range(4):native.append(C180.component(C180.mul(v,{k:F(1)}),0))
 assert len(native)==140
 d=C44.polynomial(R['d_polynomial'],vb);m1=[]
 for pol in R['R9_generators_polynomials9']:
  v=C44.mul(d,C44.polynomial(pol,vb))
  for k in range(4):m1.append(C44.component(C44.mul(v,{k:F(1)}),1).get(2,F(0)))
 m1 += [F(0)]*8
 assert m1==list(map(F,p2['true_actual_m1_row7_r_coefficients44']))
 AF=[]
 for row in X['new_A12_by9_C36']:
  assert not F(row[5][0])
  for k in range(4):
   v=C36.component(C36.mul(Q.sparse(row[5]),{k:F(1)}),0);assert not v.get(0,F(0));AF.append(v)
 assert [arr(v,9)for v in AF]==[[F(c)for c in v]for v in p2['all48_actual_adapter_row20_first9']]
 old=[];bold0=[]
 for col in oldB['source_columns640_rows56_sparse_dual9']:
  v={};b0={}
  for row,sp in col:
   if row<48:
    z=Q.pairs(sp);v=Q.plus(v,bmul(AF[row],z))
    if z.get(0,F(0)):b0[row]=z[0]
  old.append(v);bold0.append(b0)
 raw=old+native;bold0 += [{}for _ in range(140)];needed=[tuple(m)for m in whole['all15_needed_row20_X2_monomials']];assert len(needed)==15
 raw15=[[v.get(INDEX[m],F(0))for m in needed]for v in raw];assert raw15==[[F(c)for c in v]for v in p2['all780_row20_second_surrogates15']]
 # Exact source sensitivities of UNKNOWN proper u5 even second coefficients.
 sensitiv=[{}for _ in range(780)]
 for j,col in enumerate(R['native_tensor_Koszul6_HB2_resolution'][1]['columns']):
  for row,pol in col:
   if row==5:
    for e,c in pol:
     if sum(e)==1 and any(e[:5]):Q.put(sensitiv[640+4*j],e.index(1),F(c))
 assert sensitiv==[Q.pairs(v)for v in completion['unknown_u5_even2_raw_source_coefficients780_by5_sparse']]
 # C/m_R C exterior trace algebra proves all OTHER unknown core2 effects
 # have zero even component, not that those source functions vanish.
 C4=Q.FiniteC([EXP[0]],[[0],[0]])
 for i in [0,1,2,6]:
  origin={j:F(c)for j,c in enumerate(proper['kappa9'][i][:4])if F(c)};assert not origin.get(0,F(0))and not origin.get(3,F(0))
  for a,b in itertools.product(range(4),repeat=2):assert not C4.mul(C4.mul(origin,{a:F(1)}),{b:F(1)}).get(0,F(0))
 m0=[Q.sparse(v)for v in m0data['m0_24_by44']];assert m0[7]=={20:F(5,2)}
 Xcols=[{int(i):Q.pairs(sp)for i,sp in col}for col in X['all780_source_columns44_sparse_dual9']];xp=X['X0_pivot_columns'];xc=[Q.sparse(c)for c in X['all780_X0_factorization_coefficients']]
 O=W['residual_original_X_columns'];wp=W['W_pivot_residual_columns'];wr=W['W_pivot_rows'];Wi=[Q.sparse(v)for v in W['W_fixed_minor_inverse']];wcols=[Q.pairs(v)for v in W['all_W84_columns_sparse']]
 C=[i for i in range(44)if i not in m0data['selected_source_columns24']];I=[C[i]for i in minimal['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=[Q.sparse(v)for v in minimal['fixed_minor_data']['inverse']]
 hi=[[F(5,2)*Xcols[j].get(20,{}).get(i+1,F(0))for j in xp]for i in[0,4,5]];assert hi==[[F(c)for c in row]for row in p3['all54_actual_pivot_h0_readback']]
 def cres(j):return Q.plus({j:F(1)},*({xp[k]:-a}for k,a in xc[j].items()))
 def comb(vals,c):return Q.plus(*(Q.scaled(vals[k],a)for k,a in c.items()))
 columns=[];labels=[];row20values=[];checks=0
 for j,w in enumerate(wcols):
  if j in wp:continue
  alpha=Q.mvec(Wi,{i:w[r]for i,r in enumerate(wr)if w.get(r,F(0))});co=Q.plus(cres(O[j]),*(Q.scaled(cres(O[wp[k]]),-a)for k,a in alpha.items()))
  assert not comb(bold0,co)and not comb(sensitiv,co)
  x1=[{}for _ in range(44)];x0={}
  for k,a in co.items():
   for i,sp in Xcols[k].items():
    if sp.get(0,F(0)):Q.put(x0,i,a*sp[0])
    for v,b in sp.items():
     if v:Q.put(x1[i],v-1,a*b)
  assert not x0
  for row in m0:
   assert not Q.plus(*(Q.scaled(x1[i],a)for i,a in row.items()));checks+=8
  row20=[sum((a*raw15[k][i]for k,a in co.items()),F(0))for i in range(15)];row20values.append(row20);labels.append(O[j])
  a1=[Q.plus(*(Q.scaled(x1[I[k]],a)for k,a in row.items()))for row in Ji]
  final=[]
  for t in [42,43]:final.append(Q.plus(x1[t],*(Q.scaled(a1[k],-Xcols[xp[k]].get(t,{}).get(0,F(0)))for k in range(18))))
  for hindex,i in enumerate([0,4,5]):
   v={}
   for q in range(8):
    z=F(5,2)*row20[needed.index(tuple(sorted((i,q))))]if q<=i or q==7 else F(0)
    if q==7:z+=sum((m1[k]*x1[k].get(i,F(0))for k in range(44)),F(0))
    z-=sum((hi[hindex][k]*a1[k].get(q,F(0))for k in range(18)),F(0))
    if z:v[q]=z
   final.append(v)
  columns.append({8*t+i:a for t,v in enumerate(final)for i,a in v.items()})
 assert len(columns)==681 and labels==minimal['retained_original_relation_columns849'][:681]==p3['actual_final_source_labels681']
 assert row20values==[[F(c)for c in v]for v in p3['all681_actual_final_row20_second_coefficients15']]
 assert checks==681*192==p3['all681_m0_finalX1_zero_coefficients_checked']
 columns += [Q.pairs(v)for v in minimal['d3_only_linear']['all168_known_minimal_linear_columns40_sparse']]
 assert len(columns)==849 and columns==[Q.pairs(v)for v in p3['all849_actual_minimal_first_columns40_sparse']]
 cp=p3['rank_minor_columns'];rr=p3['rank_minor_rows'];inv=[Q.sparse(v)for v in p3['rank_minor_inverse']];assert len(cp)==len(rr)==35
 minor=[{j:columns[c][r]for j,c in enumerate(cp)if columns[c].get(r,F(0))}for r in rr];assert Q.mprod(minor,inv)==Q.mprod(inv,minor)==Q.eye(35)and Q.det_fixed(minor)==F(p3['rank_minor_determinant'])
 hs=hashlib.sha256()
 for col in columns:
  c=Q.mvec(inv,{i:col[r]for i,r in enumerate(rr)if col.get(r,F(0))});assert Q.plus(*(Q.scaled(columns[cp[k]],a)for k,a in c.items()))==col
  hs.update((json.dumps([[i,str(a)]for i,a in sorted(c.items())],separators=(',',':'))+'\n').encode())
 assert hs.hexdigest()==p3['rank_factorization_stream_sha256']
 assert [5,40-35]==report['actual_HF_Tnil_degree0_1']==[5,5]and report['actual_length_mod_m_squared']==10
 save('independent-source-monic-m1-and-all849-proof.json',{'own_sparse_monomial_C180_C36_C44_algebra':True,'cached_both57_source2_full_CB_F9_adapter_and_all405_shared_physical_coefficients':True,'actual11node_monic_DAG_and_all_nu45_firsts_B11_match':True,'actual_ResAT_row20_first48_and_zero_origins':True,'all780_targeted_surrogates15_verified_not_claimed_actual_individually':True,'all681_true_constant_Bold0_and_u5_even2_sensitivities_zero':True,'unknown_other_core2_even_projection_zero_C4_proof':True,'actual_m1_row7_r44_from_C44_retained':True,'all681_times192_first_Koszul_source_checks':checks,'all54_pivot_h0_and18_origin_source_cleanup':True,'ordered_Hadamard_not_symmetrized_or_zero_guessed':True,'all681_actual_row20_final_values_and849_minimal_first_columns':True,'fixed35_both_inverse_products_determinant_and_all849_factorizations':True,'actual_HF':[5,5],'actual_length_mod_m_squared':10,'no_new_source_jet_gauge_or_rank_search':True,'not_total_length_nilgens_cube_or_index':True})
 for p in[Path(__file__),Path(Q.__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 result=dict(id=ID,all_passed=True,frozen=True,actual_first_relation_rank=35,actual_HF_Tnil_degree0_1=[5,5],actual_length_mod_m_squared=10,actual_source_monic_m1_and_all849_columns_verified=True,no_producer_import_or_newsource_jet_gauge_solve_ranksearch=True,total_length_nil_generators_cube_or_exact_index_NOT_decided=True,input_sha256=pins,evidence_sha256={'independent-source-monic-m1-and-all849-proof.json':sha(OUT/'independent-source-monic-m1-and-all849-proof.json')},seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 save('report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()

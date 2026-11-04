#!/usr/bin/env python3
"""ONE existing-cache evaluation of ALL849 minimal T_nil first relations.
No source/gauge solve or new source jet. Unknown second callbacks cancel only
after the proven exact constant source combinations; raw X2 is a surrogate.
"""
import collections,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_honest_nil_T_first_values_20261104_v63 import pivots,invq,mv,polynomial,enc,nqq,sparse
from ambient_B18_native128_first_20261104_v63 import Algebra
from ambient_B18_T_native_S2_B11_20261104_v63 import Finite
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63'
PATHS={
 'completion':'reports/AMBIENT_B18_HONEST_NIL_T_LINEAR_CACHED_COMPLETION_PREFLIGHT_20261104_v63/actual-final-row20-cancellations-and-eleven-node-monic-source.json',
 'completion_report':'reports/AMBIENT_B18_HONEST_NIL_T_LINEAR_CACHED_COMPLETION_PREFLIGHT_20261104_v63/report.json',
 'whole':'reports/AMBIENT_B18_HONEST_NIL_T_LINEAR_SOURCE_WHOLE_PREFLIGHT_20261104_v63/actual-final681-source-sensitivities-and-cached-directions.json',
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'minimal_peer':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63_REVIEW/report.json',
 'm0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json',
 'X1':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'W0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase3-actual-W84-rank-and-mu.json',
 'oldB1':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
 'native1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase5-native-adapter-firsts.json',
 'proper1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase4-native-source-and-R12.json',
 'monic1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase1-source-firsts.json',
 'B11':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase3-actual-proper-G8-B11.json',
 'monicB11':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase1-source-B11-values.json',
 'source2':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'adapter2':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'nilGL':'reports/AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63/actual-E4-nil-two-generators-and-relations.json'}
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
MON=[()]+[(i,)for i in range(8)]+list(itertools.combinations_with_replacement(range(8),2));MI={m:i for i,m in enumerate(MON)}
Z=(Q(0),)*45
def add(*ps):return tuple(sum(v,Q(0))for v in zip(*ps))if ps else Z
def scale(p,c):return tuple(x*c for x in p)
def mul(p,q):
 if p==Z or q==Z:return Z
 v=[Q(0)]*45
 for i,a in enumerate(p):
  if a:
   for j,b in enumerate(q):
    if b and len(MON[i])+len(MON[j])<=2:v[MI[tuple(sorted(MON[i]+MON[j]))]]+=a*b
 return tuple(v)
def const(c):return (Q(c),)+(Q(0),)*44
def fold9(p,mons):
 v=[Q(0)]*45
 for mon,c in zip(mons,p):
  if Q(c):v[MI[tuple(sorted(7 if i==8 else i for i in mon))]]+=Q(c)
 return tuple(v)
class C180(Algebra):
 def __init__(self,nus):
  self.n=180;self.Z=(Q(0),)*180;self.ONE=self.c(1);self.nilcap=9;self.inverse_checks=0;self.table={}
  ns=[tuple(map(Q,v))for v in nus]
  for i in range(180):
   a,u=divmod(i,4)
   for j in range(180):
    b,v=divmod(j,4);mon=tuple(sorted(MON[a]+MON[b]))
    if len(mon)>2:continue
    carry=u&v;mask=u^v;coef=const(0);z=[Q(0)]*45;z[MI[mon]]=Q(1);coef=tuple(z)
    for k in range(2):
     if carry&(1<<k):coef=mul(coef,ns[k])
    nz=[(4*k+mask,c)for k,c in enumerate(coef)if c]
    if nz:self.table[i,j]=nz
def basis(a,j):v=[Q(0)]*a.n;v[j]=Q(1);return tuple(v)
def lift36(v):return tuple(map(Q,v))+(Q(0),)*144
def component(v,i):return tuple(v[4*k+i]for k in range(len(v)//4))
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_FSIZE,(32<<20,32<<20))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False);saved=0;bytesout=0
 def save(name,data):
  nonlocal saved,bytesout
  saved+=nqq(data);assert saved<150000
  blob=json.dumps(enc(data),separators=(',',':'))+'\n';bytesout+=len(blob.encode());assert bytesout<32<<20
  (out/name).write_text(blob);(out/'phase.json').write_text(json.dumps({'phase':name,'seconds':time.monotonic()-start,'QQ_saved':saved})+'\n')
  print(json.dumps({'phase':name,'seconds':time.monotonic()-start}),flush=True)
 D={k:read(ROOT/p)for k,p in PATHS.items()}
 assert D['completion_report']['all_passed']and D['minimal_peer']['all_passed']
 for p,h in D['completion_report']['input_sha256'].items():assert sha(ROOT/p)==h
 J=D['source2'];H=D['adapter2'];mons=J['source2_monomials']
 assert len(mons)==55
 YB=[]
 for lin,quad in zip(H['YB_in_YA_linear'],H['YB_in_YA_quadratic']):
  v=[Q(0)]*45
  for i,c in enumerate(lin):v[MI[(7 if i==8 else i,)]]+=Q(c)
  for mon,c in quad:
   ids=[7 if i==8 else i for i,e in enumerate(mon)for _ in range(e)];v[MI[tuple(sorted(ids))]]+=Q(c)
  YB.append(tuple(v))
 own={f:dict(J['both57_source2'][f]['all57_source2jets'])for f in ['CA','CB']}
 def ca(i):return fold9(own['CA'][i],mons)
 def cb(i):
  p=own['CB'][i];v=const(p[0])
  for mon,c in zip(mons[1:],p[1:]):
   if Q(c):
    q=const(1)
    for k in mon:q=mul(q,YB[k])
    v=add(v,scale(q,Q(c)))
  return v
 shared=[]
 for i in range(2,11):
  av,bv=ca(i),cb(i);assert av==bv,('PHYSICAL_SOURCE_ADAPTER_SECOND',i);shared.append(av)
 trace=[ca(11),ca(12),cb(20537),cb(20538)]
 nus=[add(trace[1],scale(mul(trace[0],trace[0]),Q(1,4))),add(trace[3],scale(mul(trace[2],trace[2]),Q(1,4)))]
 assert [list(v[:9])for v in nus]==[[Q(c)for c in v]for v in D['monic1']['nu_first9']]
 BMON=[(),(7,),(7,7),(0,),(1,),(4,),(5,),(0,7),(1,7),(4,7),(5,7)]
 assert [[v[MI[m]]for m in BMON]for v in nus]==[[Q(c)for c in v]for v in D['monicB11']['nu_A_B_B11']]
 save('phase1-cached-source-monics-and-physical-adapter.json',{'base_monomials45':MON,'CB_F9_in_CA_F9_ell0_degree2':YB,
  'four_literal_source_s_t_degree2':trace,'actual_nu_A_B_degree2':nus,'all9_shared_physical_degree2_polynomials':shared,
  'all_shared45_coefficients_verified':405,'all_nu_firsts_and_B11_projections_match':True,
  'new_source_gauge_equation_solve':False})
 a=C180(nus);a36=Algebra(8,2,[[Q(c)for c in v]for v in D['monic1']['nu_first9']])
 P,G,R,GL=(D[k]for k in ['proper1','native1','radical','nilGL'])
 core=[lift36(P['kappa9'][i])for i in [0,1,2,6]]
 cv=[lift36(v)for v in G['G8'][:5]]+[polynomial(a,GL['forward_GL4'][s],core)for s in ['h','c','r','l']]
 # All raw native row20 outputs are evaluated with unknown proper2 omitted.
 # Their exact final-functional insensitivity is a separately pinned premise.
 native=[]
 for col in R['native_tensor_Koszul6_HB2_resolution'][1]['columns']:
  pol=next((p for row,p in col if row==5),[]);v=polynomial(a,pol,cv)
  for k in range(4):native.append(component(a.mul(v,basis(a,k)),0))
 assert len(native)==140
 # Saved C44 provides the TRUE m1 row7, rather than assigning it zero.
 b=Finite(2,[[Q(c)for c in v]for v in D['monicB11']['nu_A_B_B11']])
 coreb=[tuple(map(Q,D['B11']['kappa9_C_B44'][i]))for i in [0,1,2,6]]
 vb=[tuple(map(Q,v))for v in D['B11']['actualG8_C_B44'][:5]]+[polynomial(b,GL['forward_GL4'][s],coreb)for s in ['h','c','r','l']]
 db=polynomial(b,R['d_polynomial'],vb);m1=[]
 for pol in R['R9_generators_polynomials9']:
  q=b.mul(db,polynomial(b,pol,vb))
  for k in range(4):m1.append(component(b.mul(q,basis(b,k)),1)[2])
 m1 += [Q(0)]*8;assert len(m1)==44
 # Actual source adapter firsts, with varying monics, not naive coefficient differentiation.
 A=D['X1']['new_A12_by9_C36'];AF=[]
 for i in range(12):
  assert not Q(A[i][5][0])
  for k in range(4):
   v=component(a36.mul(tuple(map(Q,A[i][5])),basis(a36,k)),0);assert not v[0];AF.append(v)
 old=[]
 for col in D['oldB1']['source_columns640_rows56_sparse_dual9']:
  v=Z
  for row,sp in col:
   if row<48:
    q=[Q(0)]*45
    for k,c in sp:q[k]=Q(c)
    aa=tuple(AF[row])+(Q(0),)*36
    v=add(v,mul(aa,tuple(q)))
  old.append(v)
 raw=old+native;assert len(raw)==780
 needed=[tuple(m)for m in D['whole']['all15_needed_row20_X2_monomials']]
 raw15=[[v[MI[m]]for m in needed]for v in raw]
 save('phase2-row20-surrogates-and-true-m1row7.json',{'monomials15':needed,'all780_row20_second_surrogates15':raw15,
  'raw_surrogates_are_NOT_claimed_individual_actual_X2':True,
  'exact_final681_functional_insensitivity_bound_to_completion_preflight':True,
  'true_actual_m1_row7_r_coefficients44':m1,'m1_row7_z_coefficients_zero_by_definition':True,
  'all48_actual_adapter_row20_first9':AF,'unknown_Bold2_adapter2_proper2_never_evaluated':True})
 M,X,W,Min=(D[k]for k in ['m0','X1','W0','minimal'])
 m0=[[Q(c)for c in v]for v in M['m0_24_by44']];msp=[[(i,c)for i,c in enumerate(v)if c]for v in m0]
 xp=X['X0_pivot_columns'];xco=[[Q(c)for c in v]for v in X['all780_X0_factorization_coefficients']]
 Xs=[]
 for col in X['all780_source_columns44_sparse_dual9']:
  v={}
  for row,sp in col:
   q=[Q(0)]*9
   for k,c in sp:q[k]=Q(c)
   v[row]=q
  Xs.append(v)
 zero9=[Q(0)]*9
 O=W['residual_original_X_columns'];wp=W['W_pivot_residual_columns'];wr=W['W_pivot_rows'];Wi=[[Q(c)for c in v]for v in W['W_fixed_minor_inverse']]
 # Retain only the15 targeted coefficients in the constant source transport.
 res15=[[raw15[j][i]-sum((c*raw15[xp[k]][i]for k,c in enumerate(xco[j])if c),Q(0))for i in range(15)]for j in O]
 C=[i for i in range(44)if i not in M['selected_source_columns24']]
 I=[C[i]for i in Min['fixed_minor_data']['compatible_X0_top18_rows_in_C20']]
 Ji=[[Q(c)for c in v]for v in Min['fixed_minor_data']['inverse']];Jsp=[[(i,c)for i,c in enumerate(v)if c]for v in Ji]
 out2=[v['index']for v in Min['minimal_target_labels'][:2]];assert out2==[42,43]
 hi0=[[Q(5,2)*Xs[j].get(20,zero9)[i+1]for j in xp]for i in [0,4,5]]
 def residual1(j):
  cof=[(k,c)for k,c in enumerate(xco[j])if c]
  return [[Xs[j].get(i,zero9)[q+1]-sum((c*Xs[xp[k]].get(i,zero9)[q+1]for k,c in cof),Q(0))for q in range(8)]for i in range(44)]
 pivres=[residual1(O[k])for k in wp]
 columns=[];realrow20=[];labels=[];source_checks=0
 for j,sp in enumerate(W['all_W84_columns_sparse']):
  if j in wp:continue
  w={i:Q(c)for i,c in sp};ww=[(t,w[r])for t,r in enumerate(wr)if r in w]
  aa=[sum((row[t]*c for t,c in ww),Q(0))for row in Wi];an=[(k,c)for k,c in enumerate(aa)if c];rx=residual1(O[j])
  x1=[[rx[i][q]-sum((c*pivres[k][i][q]for k,c in an),Q(0))for q in range(8)]for i in range(44)]
  assert all(sum((c*x1[i][q]for i,c in row),Q(0))==0 for row in msp for q in range(8));source_checks+=192
  row20=[res15[j][t]-sum((c*res15[wp[k]][t]for k,c in an),Q(0))for t in range(15)]
  realrow20.append(row20);labels.append(O[j])
  a1=[[sum((c*x1[I[k]][q]for k,c in row),Q(0))for q in range(8)]for row in Jsp]
  top=[[x1[t][q]-sum((Xs[xp[k]].get(t,zero9)[0]*a1[k][q]for k in range(18)),Q(0))for q in range(8)]for t in out2]
  bottom=[]
  for n,i in enumerate([0,4,5]):
   v=[]
   for q in range(8):
    mon=tuple(sorted((i,q)));value=Q(5,2)*row20[needed.index(mon)]if q<=i or q==7 else Q(0)
    if q==7:value+=sum((m1[k]*x1[k][i]for k in range(44)),Q(0))
    value-=sum((hi0[n][k]*a1[k][q]for k in range(18)),Q(0));v.append(value)
   bottom.append(v)
  columns.append(sum(top+bottom,[]))
 assert len(columns)==681 and labels==Min['retained_original_relation_columns849'][:681]
 for col in Min['d3_only_linear']['all168_known_minimal_linear_columns40_sparse']:
  v=[Q(0)]*40
  for i,c in col:v[i]=Q(c)
  columns.append(v)
 assert len(columns)==849
 cp,cr=pivots(columns);rank=len(cp);Ai,det=invq([[columns[j][i]for j in cp]for i in cr])if cp else([],Q(1))
 factorhash=hashlib.sha256()
 for c in columns:
  v=mv(Ai,[c[i]for i in cr])if cp else[]
  assert all(sum((columns[cp[k]][i]*v[k]for k in range(rank)),Q(0))==c[i]for i in range(40))
  factorhash.update((json.dumps(enc(sparse(v)),separators=(',',':'))+'\n').encode())
 hf1=40-rank
 save('phase3-all849-minimal-first-relations-and-rank.json',{'target_order':Min['minimal_target_labels'],'variable_order':D['whole']['same_S_frame'],
  'all681_actual_final_row20_second_coefficients15':realrow20,'actual_final_source_labels681':labels,
  'all849_actual_minimal_first_columns40_sparse':[sparse(v)for v in columns],
  'all681_m0_finalX1_zero_coefficients_checked':source_checks,'all54_actual_pivot_h0_readback':hi0,
  'actual_first_relation_rank':rank,'rank_minor_columns':cp,'rank_minor_rows':cr,'rank_minor_determinant':det,'rank_minor_inverse':Ai,
  'all849_rank_factorizations_verified':True,'rank_factorization_stream_sha256':factorhash.hexdigest(),
  'actual_HF_Tnil_degree0_1':[5,hf1],'actual_length_Tnil_mod_m_squared':5+hf1,
  'no_source_gauge_solve_or_fullX2_m1_evaluation':True})
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'status':'PROVED_actual_complete_first_Loewy_layer_of_Tnil_from_existing_caches',
  'minimal_targets':5,'actual_complete_linear_relation_columns':849,'actual_linear_rank':rank,
  'actual_HF_Tnil_degree0_1':[5,hf1],'actual_length_mod_m_squared':5+hf1,
  'all_unknown_second_callbacks_cancel_only_on_actual_constant_sources':True,
  'true_m1_row7_retained_and_nu2_same_physical_frame':True,'new_source_gauge_IFT_jet_or_GB':False,
  'Tnil_total_length_or_Nil_generators_cube_exact_index_decided':False,
  'saved_QQ_slots':saved,'bytes':bytesout,'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{p.name:sha(p)for p in out.glob('phase[123]-*.json')},
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()

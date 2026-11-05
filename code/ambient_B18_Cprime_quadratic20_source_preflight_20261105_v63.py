#!/usr/bin/env python3
"""Complete4x20 quadratic source gate and single cached ell² gauge profile."""
import collections,hashlib,itertools,json,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_Cprime_ell2_nilimage_source_preflight_20261105_v63 import PATHS as OLD
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_CPRIME_QUADRATIC20_SOURCE_PREFLIGHT_20261105_v63'
PATHS={k:OLD[k]for k in ['source','interface','proper','second5','preparation','J30','firstpreflight','table33','source33']}
PATHS.update({
 'Vfirst':'reports/AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63/phase1-actual-fullnine-first-V-columns.json',
 'Vpeer':'reports/AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63_REVIEW/report.json',
 'q4':'reports/AMBIENT_B18_CPRIME_PRODUCT_ANN_QUOTIENT_SOURCE_PREFLIGHT_20261105_v63/report.json',
 'annpeer':'reports/AMBIENT_B18_CPRIME_PRODUCT_ANN_SOURCE_20261105_v63_REVIEW/report.json',
 'gauge33':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63/phase2-canonical30-exact-values.json',
 'proper33':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63/phase3-actual-kappa9-and-normalized-nilframe.json',
 'gauge33peer':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63_REVIEW/report.json',
 'source2':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'adapter2':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json'})
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def main():
 t=time.monotonic();out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 for k in ['Vpeer','annpeer','gauge33peer']:assert D[k]['all_passed']
 L=D['Vfirst'];cols=[{i:Q(c)for i,c in v}for v in L['all44_actual_linear_columns27_sparse']];cp=L['fixed_pivot_columns'];rr=L['fixed_pivot_rows'];iv=[[Q(c)for c in r]for r in L['fixed_minor_inverse']];assert len(cp)==24;free=[i for i in range(44)if i not in cp];kernel=[]
 for j in free:
  co=[sum((a*cols[j].get(rr[k],Q(0))for k,a in enumerate(row)),Q(0))for row in iv];v={j:Q(1)}
  for p,c in zip(cp,co):
   if c:v[p]=-c
  assert all(sum((c*cols[k].get(i,Q(0))for k,c in v.items()),Q(0))==0 for i in range(27));kernel.append(v)
 assert len(kernel)==20 and sum(map(len,kernel))==49
 mon2=list(itertools.combinations_with_replacement(range(9),2));pi=[{i:Q(c)for i,c in v}for v in D['q4']['literal_fixed_pi4_rows_by135_sparse']];support=sorted({mon2[i%45]for row in pi for i in row});assert support==[(7,7),(7,8),(8,8)]
 MON=[()]+[(i,)for i in range(9)]+[m for m in mon2 if not(m[0]<7 and m[1]<7)];assert len(MON)==27 and MON[-1]==(8,8);M26=MON[:-1];idx26={m:i for i,m in enumerate(M26)};idx27={m:i for i,m in enumerate(MON)};plane=[(),(7,),(8,),(7,7),(7,8),(8,8)];assert all(m in MON for m in plane)
 B=D['table33'];maps=[]
 for z,r,e in B['basis33']:
  mon=tuple(([z-1]if z else[])+[7]*r+[8]*e);maps.append({idx26[mon]:Q(1)}if mon in idx26 else{})
 def multiply(v,w):
  z={}
  for i,c in v.items():
   for j,a in w.items():
    k=idx26.get(tuple(sorted(M26[i]+M26[j])))
    if k is not None:z[k]=z.get(k,Q(0))+c*a
  return {i:c for i,c in z.items()if c}
 def project(v):
  z={}
  for i,c in v:
   for j,a in maps[i].items():z[j]=z.get(j,Q(0))+Q(c)*a
  return {i:c for i,c in z.items()if c}
 assert all(project(B['table_columns_for_each_left_basis'][i][j])==multiply(maps[i],maps[j])for i in range(33)for j in range(33))
 ids=[{i for i,v in D['source2']['both57_source2'][f]['all57_source2jets']}for f in['CA','CB']];assert len(ids[0]&ids[1])==9 and len(ids[0]|ids[1])==105
 S=D['source'];nodes=S['DAG_scalar_nodes'];roots=[i for f in S['families'].values()for p in f['complete_normalized_physical_w_coefficients'][0]for i in p]+[f['actual_tracefree_nu']for f in S['families'].values()]+[i for f in['CA','CB']for i in D['source2']['both57_source2'][f]['equations57']]+[7790];seen=set();todo=list(roots);bind=ids[0]|ids[1];bad=[]
 while todo:
  i=todo.pop()
  if i in seen:continue
  seen.add(i)
  if i in bind:continue
  n=nodes[i]
  if n[0]in ['a','m','n','i']:todo.extend(n[1:])
  elif n[0]!='c':bad.append([i,n])
 assert not bad
 uses=collections.Counter(roots)
 for i in seen:
  if i not in bind and nodes[i][0]in ['a','m','n','i']:uses.update(nodes[i][1:])
 live=peak=0
 for i in sorted(seen):
  live+=27;peak=max(peak,live)
  if i not in bind and nodes[i][0]in ['a','m','n','i']:
   for j in nodes[i][1:]:
    uses[j]-=1
    if not uses[j]:live-=27
 budget={'stream_actual_physical_DAG_peak':peak,'cached_full9_source2_used':6270,'common105_base27':2835,'physical150_pairs_and_nus':16254,'gauge30_base27':810,'J30_and_RHS':930,'actual14_kappa9_C108':2484,'finite_GL_source_work':43200,'fixedVfirst_q4_and_kernel_metadata':8000,'actual44_projection20_and_cubic_stream':6000,'new_outputs_and_overlap_receipts':18000};whole=sum(budget.values());assert whole<150000
 data={'actual_full9_order':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA','ell'],'base27_monomials':MON,'base26_monomials':M26,'new_scalar_monomial_only':(8,8),'required_plane6_monomials':plane,
  'formal_base27_definition':'R0/(m³, all z_i*z_j), keeping allnine firsts and17seconds. Its26quotient additionallykillsell².',
  'all1089_existing33_to_base26_products_verified':True,'actual33_to26_columns':[[[i,c]for i,c in sorted(v.items())]for v in maps],
  'full9_V1_kernel_free_source_labels':free,'all20_constant_kernel_columns44_sparse':[[[i,c]for i,c in sorted(v.items())]for v in kernel],'all540_first_kernel_coefficients_zero':True,'free_tail_I20':True,
  'complete_quadratic_formula':'NF4 * V2 * K20; everydegree2initial is V2*a0+V1*a1 with V1*a0=0, hencea0=K20*c and NF4 killsV1*a1. ExactHF2(Cprime)=4-rank(actual4x20).',
  'actual_NF4_rows_by135_sparse':D['q4']['literal_fixed_pi4_rows_by135_sparse'],'required_second_monomials_only':support,'all_nonplane_highest_second_callbacks_zero_by_literal_support':True,
  'physical_cached_transport':'Compose ownCA57 full9source2 with YA=(z7,r,r+ell/5). Compute YB using ALLnine full linear+quadratic adapter2 rows, retaining CBzero7curvature. Compose ownCB57 withthatYB, verifyall9 sharedphysicalvariables andall114originalequations beforegauge.',
  'all_old16_second_cache_overlap_gate':'Project recomposed105 andrequested164roots toB26 andcompareliteralactual33projections. Forremainingphysical150roots, equalityfollowsfromsameoriginalDAG andthe105coefficients/algebra-mapidentity; recordallrootsandverifytheirsourceevaluations.',
  'canonical30_recipe':'Projectactual33gauge toB26, appendzero ONLYasinitialell² seed, evaluatefull30constraints inB27 and requireallold26coefficients0. Correct30newcoefficients using SAMEfixedJ30^-1; verifyall30residuals exactly0 andallB26 gauge/properoverlaps. Neverassignell² actualvalue0.',
  'new_derived_gauge_scalar_count':30,'new_physical_source_solve_or_source3':False,'unittrace_carry':'Keep actualnu through2, including contributions fromnonzerooddgroupcenter; evaluatefullproperpolynomial/rational formula BEFOREtrace extraction.',
  'V2_scope':'OnlypureCA/CB normalizerunit functions/proper30; no oldr/g divisions, honestX, m2 orunknownquotient gradients occur inVfull. R0-degree2functions depend on source/gauge through2 by regular algebra pullback, regardless nilpotenttau constants.',
  'optional_cubic_source_count_bound':2046,'optional_cubic_gate':'Ifh2<=3, canonicalNFrepresentatives ofselectedactualquadraticclasses are trueinitialrelations modknownlinearproducts. Stream44*S2 + actualrank<=4quadrics*S1 +3ell³ +allconsistency. Positivecum9/zero3 proveslongexponent2; otherwisefrontieronly.',
  'source_DAG_ancestor_count':len(seen),'budget':budget,'whole_saved_live_QQ_upper':whole,'hard_seconds':60,'hard_memory_bytes':2<<30,'hard_total_bytes':32<<20,'derived_bit_guard':2048,
  'actual_gauge_or_quadratic_values_not_computed':True,'prepared_worker_hash_required_before_value_GO':True}
 p=out/'complete-four-by20-and-ell-square-cache-source-gate.json';p.write_text(json.dumps(enc(data),separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_PROFILE_complete4x20_and_cached_single_gauge_pair','complete_actual_quadratic_shape':[4,20],'source_kernel_nonzero_coefficients':49,'new_derived_gauge_coefficients':30,'new_physical_source_values':0,'all114_residuals_required':True,'whole_saved_live_QQ_upper':whole,'actual_quadratic_rank_or_long_exponent_computed':False,'prepared_worker_hash_gate_pending':True,'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()

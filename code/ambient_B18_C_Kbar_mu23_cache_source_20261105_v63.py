#!/usr/bin/env python3
"""SOURCE/cache profile only: old M2 lastrow reduces to21 given scalar roots."""
import hashlib,json,math,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
import ambient_B18_C_Kbar_r2_initial_source_20261105_v63 as I
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KBAR_MU23_CACHE_SOURCE_20261105_v63'
PATHS={
 'inventory':'reports/AMBIENT_B18_C_KBAR_MU23_SOURCE_INVENTORY_20261105_v63/fourteen-z-actions-common37-M2-source-functionals-and-minimal-generator-criterion.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'C0_peer':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json',
 'source37':'reports/AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63/actual119-cubic-X2-only-source-operator-preflight.json',
 'unit':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
 'unit_peer':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63_REVIEW/report.json',
 'schur':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/actual-minimal5-by825-Schur-operator.json',
 'schur_peer':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW/report.json',
 'first640':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
 'first640_peer':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63_REVIEW/report.json',
 'firstQ':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
 'm1':'reports/AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/actual-projected-three-class-X2-operator-preflight.json',
 'cache':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD/phase3-oldB-common-source2-values-and-symbols.json',
 'cache_report':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD/report.json',
 'cache_peer':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def fold(dst,k,c):
 if c:
  v=I.addq(dst.get(k,Q(0)),c)
  if v:dst[k]=v
  elif k in dst:del dst[k]
def env(xs):
 xs=[I.parse(x)for x in xs];den=math.lcm(*(x.denominator for x in xs));num=max(abs(x.numerator)*(den//x.denominator)for x in xs);return num,den
def hb(pair):return max(abs(pair[0]).bit_length(),pair[1].bit_length())
def main():
 t=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={
 'inventory':['all14_affine_rows_constant_then37_unknown_t','all14_z_U_S_full_source_and_coherent_M2_functionals','known_rU_rS_firsts','actual_pi_q2_degree2_sparse'],
 'C0':['honest_D0_relation_columns'],
 'source37':['actual_selected37_original768_constant_combinations','all37_m0_Xj1_zero_and_pi_Xj1_equal_actual_Lj'],
 'unit':['typed_n56_by24','typed_regular_function_operators'],
 'schur':['actual_free_five_source_rows','constant_free32_quotient_projection5_rows','first_schur_formula','second_schur_formula','actual_relation_order'],
 'first640':['source_columns640_rows56_sparse_dual9'],
 'firstQ':['normalized84_unitcycle_Q1_columns_sparse'],
 'm1':['actual_pivot_X1_by27_sparse56_times8'],
 'cache':['base20_monomials_original8','C0_composite_values_and_opaque_second_gradients'],
 }
 d={}
 for k,p in PATHS.items():
  d[k]=S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k.endswith('_peer'):assert d[k]['all_passed']
 sourceprior=S.count(d)+30000;assert sourceprior<150000
 assert d['unit']['typed_n56_by24'][55]==[0]*24 and d['unit']['typed_regular_function_operators'][0]==['zero']
 assert d['schur']['actual_free_five_source_rows']==[22,23,29,30,55]
 assert d['schur']['constant_free32_quotient_projection5_rows'][4]==[[31,'1']]
 assert d['schur']['actual_relation_order']==['168 literal d3 images','657 remaining actual source cycles']
 assert d['inventory']['actual_pi_q2_degree2_sparse']==[[179,'1']]
 first=d['first640']['source_columns640_rows56_sparse_dual9'];assert len(first)==640
 for col in first:
  for row,values in col:
   if row==55:assert not Q(dict(values).get(8,'0'))
 assert all(not Q(dict(col).get(447,'0'))for col in d['m1']['actual_pivot_X1_by27_sparse56_times8'])
 assert all(not Q(dict(col).get(39,'0'))for col in d['firstQ']['normalized84_unitcycle_Q1_columns_sparse'])
 assert d['source37']['all37_m0_Xj1_zero_and_pi_Xj1_equal_actual_Lj']
 P=d['C0']['honest_D0_relation_columns'];assert len(P)==640 and all(len(col)==9 for col in P)
 combos=d['source37']['actual_selected37_original768_constant_combinations'];assert len(combos)==37
 roots37=[]
 for col in combos:
  row={}
  for k,c in col:
   assert 0<=k<768
   if k<640 and P[k][8]!=0:fold(row,P[k][8],I.parse(c))
  roots37.append(row)
 roots=sorted({j for row in roots37 for j in row});assert len(roots)==21 and sum(map(len,roots37))==98
 cache={i:(v,g)for i,v,g in d['cache']['C0_composite_values_and_opaque_second_gradients']}
 assert all(i in cache and cache[i][1]==[]for i in roots)
 MON=d['cache']['base20_monomials_original8'];assert len(MON)==20 and MON[19]==[7,7]
 forms=[[I.parse(c)for c in row]for row in d['inventory']['all14_affine_rows_constant_then37_unknown_t']];assert len(forms)==14 and all(row[0]==0 for row in forms)
 roots14=[]
 for f in forms:
  row={}
  for a,r in zip(f[1:],roots37):
   for j,c in r.items():fold(row,j,I.mulq(a,c))
  roots14.append(row)
 # All bounds use complete given polynomial entries, no requested rr
 # coefficient is selected or paired during SOURCE.
 Ecache=env(c for i in roots for m,c in cache[i][0]);E37=env(c for row in roots37 for c in row.values());E14=env(c for row in roots14 for c in row.values());Ef=env(c for row in forms for c in row[1:])
 def bound(weights,values,n):return(n*weights[0]*values[0],weights[1]*values[1])
 Et=bound(E37,Ecache,21);Eaction=bound(E14,Ecache,21);Etwo=bound(Ef,Et,37)
 lam=I.parse(d['inventory']['known_rU_rS_firsts']['lambdaS']);Elam=(abs(lam.numerator),lam.denominator)
 Eminor=(Eaction[0]*Elam[0]**2,Eaction[1]*Elam[1]**2)
 H=max(map(hb,[Ecache,E37,E14,Ef,Et,Eaction,Etwo,Eminor]));assert H<=384 and 2*H+1<=1153
 evidence={
 'scope':'FreshSOURCE/cache gate only. No rootrr callbackpairedwithweights, no37M2or14actionvalues.',
 'input_sha256':pins,
 'exact_old_lastrow_schur_proof':{'same_OLD5x825':True,'typed_nrow55_all24_zero':True,'pi5_lastrow_literal_raw55':True,'all640_raw55_r_firsts_zero':True,'all27_pivot_raw55_r_firsts_zero':True,'all84_Ystar_head4_r_firsts_zero':True,'selected37_zero_origin_constant_source_combinations':True,'formula':'R37_2(head4,r²)=Xraw37_2(row55,r²). B1_r*C1_r andYstar1_r*Wrest1_r arezero; lowerC0/Wrest0zero. NoM2/homotopyvaluequeried.'},
 'all37_Pg_source_weights':[[[j,str(c)]for j,c in sorted(row.items())]for row in roots37],
 'all14_direct_Pg_source_weights':[[[j,str(c)]for j,c in sorted(row.items())]for row in roots14],
 'all14_same37_weights':[list(map(str,row[1:]))for row in forms],
 'all14_z_action_labels':[v['source']for v in d['inventory']['all14_z_U_S_full_source_and_coherent_M2_functionals']],
 'all21_given_scalar_cache_rows_not_selected_rr':[[i,cache[i][0]]for i in roots],
 'all21_regular_second_opaque_gradients_empty':True,
 'actual_cache_monomials_original8':MON,'future_rr_index':19,'known_lambdaS':str(lam),
 'source_namespace':'LiteralP9lastrowfromactualhonestC0relation640, native128lowerP0, sameOLD768sourceweights. ProvidednonlinearCBtransport/common105cache. NoC0gA/f63namespace substitution.',
 'future_actual_worker_checks':['all37constant_and_r_first_readback_zero BEFORErr','all37rawtvalues','all14two_stage_equals_direct21rootfoldtermwise','mu2ifanyactionnonzero_elsemu3;unpairedunitminor lambdaS²*nonzerowitness'],
 'positive_whole_future_height':{'cache':Ecache,'weights37':E37,'weights14':E14,'weights14_by37':Ef,'all37t':Et,'all14direct':Eaction,'all14two_stage':Etwo,'witness_unit_minor':Eminor,'committed':H,'temporary':2*H+1,'hard_committed':384,'hard_temporary':1153},
 'future_value_scope':'Oneexistingcache37+14dotjob ONLYafterfreshcode/custody/independentSOURCEpeer+ROOTGO. No newIFT/gauge/physical2/3, nofullactions/lifts/H0/nil.',
 'all_actual_t_actions_UNCOMPUTED':True,'selective_receipts':receipts,
 }
 numeric=S.count(evidence);future_static=2*numeric+12000;assert future_static<60000
 freshraw=sum((ROOT/p).stat().st_size for p in PATHS.values());assert freshraw+800000<33554432
 evidence['future_numeric_profile']={'compiled_input_numeric_withlabels':numeric,'static_owned_upper':future_static,'hard_owned':60000,'global_cap':150000,'RAW_compiled_plus_output_upper':800000,'RAW_SOURCE_allinputs_plusreserve':freshraw+800000,'seconds':60,'memory_bytes':2000000000}
 owned=max(sourceprior,S.count(d)+S.count(roots37)+S.count(roots14)+S.count(evidence)+25000);assert owned<150000;evidence['whole_SOURCE_owned_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);p=out/'complete-old-Pg-second-cache-source-and-fourteen-action-positive-gate.json';p.write_text(json.dumps(evidence,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_ALL_MU23_INPUTS_ALREADY_CACHED_NO_ACTIONVALUES','all37_actual_second_callbacks_evaluated':False,'all14_actual_actions_evaluated':False,'roots':len(roots),'weights37':sum(map(len,roots37)),'whole_SOURCE_owned_upper':owned,'future_static_owned':future_static,'future_height':H,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()

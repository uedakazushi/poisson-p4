#!/usr/bin/env python3
"""Typed Told/Tnil native overlap and conditional third-source inventory only."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_JOINT_NATIVE_THIRD_SOURCE_BINDING_PREFLIGHT_20261104_v63'
PATHS={
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'native':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'joint3':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
 'source2':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63/phase1-actual105-source2-and-copy.json',
 'source2_peer':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD_REVIEW/report.json',
 'gauge2':'reports/AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63/phase2-canonical270-zz-gauge.json',
 'proper2':'reports/AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63/phase3-actual-proper-G8-base20.json',
 'gauge2_peer':'reports/AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63_REVIEW/report.json',
 'nu2':'reports/AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63/phase1-actual-physical-source-base20.json',
 'high':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_HIGH_ORDER_ADJOINT_20261104_v63/preflight.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def poly(p):return {tuple(m):Q(c) for m,c in p if Q(c)}
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()};dump(out,'phase.json',{'phase':'joint_binding_inputs_saved','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert D['source2_peer']['all_passed'] and D['gauge2_peer']['all_passed']
 R=D['radical'];N=D['native'];assert R['C_basis']==N['basis']==['1','tauA','tauB','tauA*tauB']
 assert R['accepted_native_C_coefficient_interface']['sha256']==sha(ROOT/PATHS['native'])
 assert D['source']['physical_zero_base_rows9']==[27,108,109,112,113,117,118,119,61]
 rules=R['actual_C_variable_rules'];assert rules['GL4_variable_input_order']==['a','d','x','t']
 want={
  'h':{(1,0,0,0):Q(1),(0,1,0,0):Q(-1),(0,0,1,0):Q(1),(0,0,0,1):Q(-2)},
  'c':{(0,1,0,0):Q(1,2),(0,0,1,0):Q(1,2)},
  'r':{(1,0,0,0):Q(-1,3),(0,1,0,0):Q(1,6),(0,0,1,0):Q(-1,6)},
  'l':{(1,0,0,0):Q(1,3),(0,1,0,0):Q(1,3),(0,0,1,0):Q(-1,3)}}
 assert {k:poly(v) for k,v in rules['actual_core_GL4'].items()}==want
 d2=R['native_tensor_Koszul6_HB2_resolution'][1];assert d2['shape']==[9,35]
 row={j:poly(p) for j,col in enumerate(d2['columns']) for i,p in col if i==5}
 def unit(i):return tuple(int(j==i) for j in range(9))
 expected={j:{unit(i):Q(1)} for j,i in zip([4,8,11,13,14],range(5))}
 for col,(i,j) in zip([30,31,32],[(6,7),(6,8),(7,8)]):expected[col]={tuple(int(k==i)+int(k==j) for k in range(9)):Q(-1)}
 assert row==expected
 told=[(0,7,7),(1,7,7),(4,7,7),(5,7,7),(7,7,7)];tnil=[(0,7,7),(4,7,7),(5,7,7)]
 assert D['high']['high_order_X3_monomials']==list(map(list,told))
 old=[tuple(m) for m in D['source2']['base20_monomials_original8']];base=old+told
 assert len(set(base))==25 and all(tuple(m[:k]+m[k+1:]) in base for m in base if m for k in range(len(m)))
 assert all(tuple(sorted(m+(v,))) not in base for m in told for v in range(8))
 source105={i for i,p in D['source2']['all105_actual_source2_values']};joint105={i for i,p in D['joint3']['all105_unknown_jets']}
 assert source105==joint105 and len(source105)==105
 assert all(not Q(p[0]) for p in D['nu2']['nu_A_B_base20'])
 ev={'source_manifest':inputs,'same_Cfree4_basis_and_source_interface_checked':True,
  'Told_native_operator':'Res4(N^t*S(G8)), scalar physical output rows0..47',
  'Tnil_native_operator':'Res4(Srad35), scalar row20=generatorh index5 with trace mask0; no N^t on this native block',
  'Tnil_native_h_row_nonzero_columns':[[j,[[list(m),str(c)] for m,c in p.items()]] for j,p in sorted(row.items())],
  'shared_actual_proper_inputs':'same kappa9, core=kappa[0,1,2,6], u5=G8 first5, canonicalgroup30, tauA/tauB monics and external CA F9 cut variables',
  'Told_target_triples':told,'Tnil_target_triples':tnil,'native_overlap_triples':sorted(set(told)&set(tnil)),
  'joint_base25_monomials':base,'base25_is_downward_closed_monomial_quotient':True,'monic_cover_rank4_dimension':100,
  'new_cubic_base_ideal_annihilated_by_maximal_base_ideal':True,
  'already_cached_pure_r3':'all105 source0..3 on originalr/s plane, restricted r=s; same105 IDs as current nonlinear source2 cache',
  'candidate_uncomputed_mixed_physical_thirds':{'triples':told[:4],'coefficients_per105':4,'candidate_total':420,
   'interpretation':'four contracted common source directions, not a full third jet; sourceCA57 plus commonCB48 block count 57*4+48*4',
   'necessity_or_execution_authorization':False},
  'candidate_joint_canonical_gauge_thirds':{'five_triples_times30':150,'fixed_unit_J30_would_suffice_if_RHS_is_bound':True,'necessity_or_execution_authorization':False},
  'Tnil_only_candidate_counts':{'three_mixed_physical105_coefficients':315,'three_gauge30_coefficients':90,'not_separate_jobs':True},
  'native_primitive_coverage':'u5 plus actual corea,d,x,t and monicnuA/nuB third callbacks determine both native blocks by finite polynomial/rational evaluation; adjoints should select a smaller subset before any allocation',
  'uncovered_old_side_gates':['Told rawoldP4/Xp3 ancestor and honest quotient sensitivities','Tnil oldBrad=A^t*Bold with A=N^-1*D8x9 requires its own higher cancellation proof'],
  'row20_native_blocks_not_identified_with_each_other':True,'no_source3_gauge3_source6_or_new_coefficient_values_computed':True}
 ep=out/'joint-native-third-source-binding-and-conditional-counts.json';dump(out,ep.name,ev)
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'exact_native_basis_interface_and_h_row_binding_checked':True,
  'conditional_joint_base25_cover100_inventory_proved':True,'new_source_or_gauge_coefficients_authorized_or_computed':False,
  'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL','report_sha256':sha(out/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()

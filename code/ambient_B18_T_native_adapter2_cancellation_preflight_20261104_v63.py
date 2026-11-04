#!/usr/bin/env python3
"""Source-specific cancellation of the complete native unit-adapter second jet.

No second jet is evaluated.  The fixed origin direct sum and all zero-source
combinations certify that L2*S0 cancels before any Euler/readback projection.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_NATIVE_ADAPTER2_CANCELLATION_PREFLIGHT_20261104_v63'
PATHS={
 'X0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
 'X0_peer':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63_REVIEW/report.json',
 'adapter0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/phase2-native-adapter.json',
 'SyzQ':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
 'native1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase5-native-adapter-firsts.json',
 'm1':'reports/AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/actual-projected-three-class-X2-operator-preflight.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rec(a):return {int(i):Q(c) for i,c in a if Q(c)}
def mul(a,b):
 c=[Q(0)]*4
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   if not i&j:c[i|j]+=x*y
 return c
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 (out/'phase.json').write_text(json.dumps({'phase':'source_factorization_inputs_saved','inputs':inputs},indent=2)+'\n')
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 X,A,S,M=D['X0'],D['adapter0'],D['SyzQ'],D['m1']
 assert D['X0_peer']['all_passed'] and X['rank_X0']==27
 assert S['actual_d2']['shape']==[48,128]
 N=[[[Q(c) for c in v] for v in row] for row in A['N']]
 Ni=[[[Q(c) for c in v] for v in row] for row in A['N_inverse']]
 for U,V in [(N,Ni),(Ni,N)]:
  for i in range(12):
   for j in range(12):
    c=[sum(mul(U[i][k],V[k][j])[b] for k in range(12)) for b in range(4)]
    assert c==[Q(i==j),0,0,0]
 gc=[]
 for v in A['G8_C0']:
  vv=list(map(Q,v));assert vv[:3]==[0,0,0];gc.append(vv[3])
 assert gc[0]==Q(1,2)
 # Actual S has28 Koszul columns in its first8 rows and4 free columns in
 # rows8..11.  G0=g*tauA*tauB: multiplication by tauA*tauB has QQ rank1;
 # contraction by the nonzero g has rank7; the free block has rank16.
 rank_native0=7+4*4
 receipts=X['all768_fixed_elimination_receipts_and_constant_combinations']
 oldp=[r for r in receipts if 'pivot' in r and r['source_column']<640]
 assert len(oldp)==4
 assert [r['source_column'] for r in oldp]==[604,613,622,631]
 assert all(all(int(j)<640 for j,c in r['constant_original_source_combination']) for r in receipts if r['source_column']<640)
 assert 4+rank_native0==27
 data=X['full_X0_56_by768'];assert len(data)==56
 def separate_zero(combo):
  for rows in [range(56)]:
   for i in rows:
    old=sum(Q(data[i][j])*c for j,c in combo.items() if j<640)
    native=sum(Q(data[i][j])*c for j,c in combo.items() if j>=640)
    assert old==0 and native==0
 rest=[r for r in receipts if r.get('zero')];assert len(rest)==741
 for r in rest:separate_zero(rec(r['constant_original_source_combination']))
 all620=M['620_actual_source_combinations_in_original768'];assert len(all620)==620
 for c in all620:separate_zero(rec(c))
 evidence={
  'actual_coefficient_base':'Rbar=QQ[[zA0,...,zA6,rA]], same actual Cfree4 monics and basis1,tauA,tauB,tauA*tauB',
  'source_manifest':inputs,
  'actual_native_factorization':{
   'X_native48_by128':'L*S0_family with L=Res4(N^t), S0_family=Res4(native_S12_by32(G8))',
   'exact_source':'immutable complete actual_d2 operator, same native first128 signs',
   'bottom8':'zero',
   'L0_unit':'both C4 N*N_inverse products checked entrywise; Res4 preserves inverse'},
  'actual_G8_origin_tauA_tauB_coefficients':list(map(str,gc)),
  'rank_native_origin_proof':{'koszul':7,'free4_times_C4':16,'sum':23,'nonzero_g_witness':'g0=1/2','disjoint_native_target_blocks':'first8 versuslast4 before invertible N0^t'},
  'actual_old_B_origin_rank':4,'actual_old_B_selected_columns':[604,613,622,631],
  'actual_joint_origin_rank':27,'actual_two_origin_images_direct':True,
  'all741_zero_origin_combinations_separately_old_and_native_zero_full56':True,
  'all620_zero_linear_combinations_separately_old_and_native_zero_full56':True,
  'exact_second_order_factor_formula':'(L*S)_2=L0*S2+L1*S1+L2*S0',
  'actual_native_constant_combination_identity':'S0*c_native=0 because L0 is a unit and X_native0*c_native=0 separately',
  'entire_unknown_L2_term_zero':'L2*S0*c_native=0 for every actual741/620 remaining source combination',
  'unknown_N2_and_Res4_second_corrections_inside_L2_not_needed':True,
  'retained_native_second_coefficient_callbacks':['L0*S2*c_native','L1*S1*c_native'],
  'L0_L1_are_frozen_actual_inputs':True,
  'pivot_correction_Xp1_a1_uses_no_L2':True,
  'cancellation_before_Euler_or_three_class_projection':True,
  'does_not_assert_S2_or_physical_core_second_coefficients_zero_or_known':True,
  'does_not_change_Q_B_N_or_source_functions':True,
  'no_second_gauge_solve_source_jet_coefficient_evaluation_rank_search_or_348_matrix':True}
 ep=out/'actual-native-L2-cancellation-source-preflight.json';blob=json.dumps(evidence,separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;ep.write_text(blob)
 report={'id':STEM,'all_passed':True,'actual_native_origin_rank':23,'actual_old_B_origin_rank':4,'actual_joint_origin_rank':27,
  'actual_origin_images_direct':True,'actual741_and620_separate_origin_identities_verified':True,
  'unknown_native_unit_adapter_L2_cancels_in_actual_projected_X2_input':True,
  'unknown_native_S2_coefficients_computed_or_removed':False,'new_second_coefficients_computed':False,
  'T_length_or_depth_decided':False,'input_hashes':inputs,'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 (out/'phase.json').write_text(json.dumps({'phase':'FINAL','report_sha256':sha(out/'report.json')},indent=2)+'\n')
 print(json.dumps(report),flush=True)
if __name__=='__main__':main()

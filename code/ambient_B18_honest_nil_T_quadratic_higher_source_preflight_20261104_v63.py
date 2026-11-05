#!/usr/bin/env python3
"""Full661 cached-first correction and exact highest native/m2 source gate.
No actual X3, m2, higher source/gauge coefficients are evaluated.
"""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_honest_nil_T_first_values_20261104_v63 import mv,enc,nqq,sparse,polynomial
from ambient_B18_native128_first_20261104_v63 import Algebra
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_QUADRATIC_HIGHER_SOURCE_PREFLIGHT_20261104_v63'
PATHS={
 'operator661':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_661_OPERATOR_20261104_v63/actual-complete-six-by661-quadratic-source-operator.json',
 'X1':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'W0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase3-actual-W84-rank-and-mu.json',
 'first2':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase2-row20-surrogates-and-true-m1row7.json',
 'first3':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json',
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'm0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json',
 'proper1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase4-native-source-and-R12.json',
 'G1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase5-native-adapter-firsts.json',
 'GL':'reports/AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63/actual-E4-nil-two-generators-and-relations.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'u5_gate':'reports/AMBIENT_B18_HONEST_NIL_T_LINEAR_CACHED_COMPLETION_PREFLIGHT_20261104_v63/actual-final-row20-cancellations-and-eleven-node-monic-source.json',
 'u5_report':'reports/AMBIENT_B18_HONEST_NIL_T_LINEAR_CACHED_COMPLETION_PREFLIGHT_20261104_v63/report.json',
 'N2B11':'reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63/phase2-actual-N-B11.json',
 'N2B11_peer':'reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63_REVIEW/report.json',
 'new270':'reports/AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63/phase3-actual-proper-G8-base20.json',
 'new270_peer':'reports/AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63_REVIEW/report.json'}
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 D={k:read(ROOT/p)for k,p in PATHS.items()}
 assert D['N2B11_peer']['all_passed']and D['new270_peer']['all_passed']
 assert D['u5_report']['all_passed']and D['u5_gate']['all681_unknown_u5_even2_sensitivities_zero']
 (out/'phase.json').write_text(json.dumps({'phase':'all_frozen_inputs_before_cached_first_arithmetic','source_hashes':{p:sha(ROOT/p)for p in PATHS.values()}},indent=2)+'\n')
 X,W,Min,M,O=(D[k]for k in ['X1','W0','minimal','m0','operator661'])
 raw=[]
 for col in X['all780_source_columns44_sparse_dual9']:
  a={}
  for row,sp in col:
   v=[Q(0)]*9
   for k,c in sp:v[k]=Q(c)
   a[row]=v
  raw.append(a)
 Z=[Q(0)]*9;xp=X['X0_pivot_columns'];xc=[[Q(c)for c in v]for v in X['all780_X0_factorization_coefficients']]
 order=W['residual_original_X_columns'];wp=W['W_pivot_residual_columns'];wr=W['W_pivot_rows'];Wi=[[Q(c)for c in v]for v in W['W_fixed_minor_inverse']]
 final=[j for j in range(762)if j not in wp]
 assert [order[j]for j in final]==D['first3']['actual_final_source_labels681']
 def first_after18(j):
  nz=[(k,c)for k,c in enumerate(xc[j])if c]
  return [[raw[j].get(i,Z)[q+1]-sum((c*raw[xp[k]].get(i,Z)[q+1]for k,c in nz),Q(0))for q in range(8)]for i in range(44)]
 pivfirst=[first_after18(order[k])for k in wp]
 wfac=[]
 for sp in W['all_W84_columns_sparse']:
  a={i:Q(c)for i,c in sp};nz=[(j,a[r])for j,r in enumerate(wr)if r in a]
  wfac.append([sum((row[j]*c for j,c in nz),Q(0))for row in Wi])
 def first_minimal(j):
  k=final[j];v=first_after18(order[k]);nz=[(p,c)for p,c in enumerate(wfac[k])if c]
  return [[v[i][q]-sum((c*pivfirst[p][i][q]for p,c in nz),Q(0))for q in range(8)]for i in range(44)]
 # All source661 columns have a unique optional free-X entry and old35 basis.
 lbase=sorted({j for c in O['all661_source_columns_in_original849']for j,a in c if j<681 and j not in O['all646_Xonly_zero_linear_source_labels']})
 assert len(lbase)==35
 fixed={j:first_minimal(j)for j in lbase}
 C=[i for i in range(44)if i not in M['selected_source_columns24']]
 I=[C[j]for j in Min['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=[[Q(c)for c in v]for v in Min['fixed_minor_data']['inverse']]
 a1all=[];fullchecks=0;mu2checks=0
 for source in O['all661_source_columns_in_original849']:
  z=[[Q(0)]*8 for _ in range(44)]
  for j,a in source:
   if j>=681:continue # mixed d3 source retained separately, X-part is zero
   v=fixed[j]if j in fixed else first_minimal(j);a=Q(a)
   for i in range(44):
    for q in range(8):
     if v[i][q]:z[i][q]+=a*v[i][q]
  a1=[[sum((row[k]*z[I[k]][q]for k in range(18)if row[k]),Q(0))for q in range(8)]for row in Ji]
  for i in range(44):
   for q in range(8):
    assert z[i][q]==sum((raw[p].get(i,Z)[0]*a1[k][q]for k,p in enumerate(xp)if raw[p].get(i,Z)[0]),Q(0))
    fullchecks+=1
  # Formal sensitivity coefficients of all44 independent mu2 row entries.
  for zi in [0,4,5]:
   for i in range(44):
    assert z[i][zi]-sum((raw[p].get(i,Z)[0]*a1[k][zi]for k,p in enumerate(xp)if raw[p].get(i,Z)[0]),Q(0))==0;mu2checks+=1
  a1all.append(sparse(sum(a1,[])))
 assert sum(len(v)for v in a1all)<40000
 mu1=[Q(c)for c in D['first2']['true_actual_m1_row7_r_coefficients44']]
 # Actual restricted cycle identity (mX)_7|z0=0, checked through known order1.
 for c in raw:assert sum((mu1[i]*v[0]for i,v in c.items()),Q(0))+Q(5,2)*c.get(20,Z)[8]==0
 m2contracts=[{'pivot_original_source':p,'known_first_offset':-sum((mu1[i]*v[8]for i,v in raw[p].items()),Q(0)),
  'remaining_actual_callback':['X2_scalar_row20',p,'rA^2'],'callback_coefficient':Q(-5,2)}for p in xp]
 # Source-highest native responses: u5_even,n and both monic_n for ANY n>=2.
 # They are degree-n base perturbations, so coefficient response is at C4 origin.
 R=D['radical'];proper=D['proper1'];G=D['G1'];nativecols=R['native_tensor_Koszul6_HB2_resolution'][1]['columns']
 def cv0(alg):
  core=[tuple(map(Q,proper['kappa9'][i][:4]))+(Q(0),)*4 for i in [0,1,2,6]]
  return [tuple(map(Q,v[:4]))+(Q(0),)*4 for v in G['G8'][:5]]+[polynomial(alg,D['GL']['forward_GL4'][s],core)for s in ['h','c','r','l']]
 native_responses=[]
 for monic in [0,1]:
  nu=[[Q(0),Q(int(k==monic))]for k in range(2)];a=Algebra(1,2,nu);vars=cv0(a);v=[]
  for col in nativecols:
   p=next((p for row,p in col if row==5),[]);q=polynomial(a,p,vars)
   for mask in range(4):
    b=[Q(0)]*8;b[mask]=Q(1);v.append(a.mul(q,tuple(b))[4])
  native_responses.append([Q(0)]*640+v)
 highest=[]
 for v in native_responses:
  res=[v[j]-sum((c*v[xp[k]]for k,c in enumerate(xc[j])if c),Q(0))for j in order]
  ans=[res[k]-sum((c*res[wp[p]]for p,c in enumerate(wfac[k])if c),Q(0))for k in final]
  highest.append(ans)
 native_allzero=all(not c for row in highest for c in row)
 data={'actual_source_split':[646,15],'all661_full44_by8_X1u_equals_Xp0_a1_verified':True,
  'all661_a1_18_by8_sparse':a1all,'all232672_source_first_coefficient_identities':fullchecks,
  'all87252_formal_mu2_sensitivity_identity_coefficients':mu2checks,
  'factored_actual_X_source661':'originalX780 * (injection_nonXp−Xp_factor)*(injection_nonW−W_factor)*C661_Xblock; mixed d3 block is kept separately and has no X component',
  'canonical_m7_structure':'mu0=(5/2)e20, mu1=r*saved44, mu2=r²*unknown44; no z derivatives at any order',
  'all780_known_pure_r_first_cycle_coefficients_zero':True,
  'exact18_mu2_Xp0_replacement_identities':m2contracts,
  'replacement_proof':'Restrict d1*m*Xp=0 to z0: r*(mXp)_7=0 in QQ[[r]], so (mXp)_7=0. Its r² coefficient is mu2Xp0+mu1Xp1_r+(5/2)Xp2_row20,r²=0. Since X1u=Xp0a1, every retained mu2X1u_zi term equals these18 contractions times a1_zi. No mu2 value is set0.',
  'm2_separate_higher_source_input_eliminated':True,
  'genuine_second_unit_inverse_derivatives_not_needed':'After constant unit normalization, M2 on a zero-M1 source is D2−C1U0^-1B1; U2 cannot occur because both offdiagonal constant blocks vanish.',
  'pivot_first_carries_NOT_declared_cancelled':True,
  'native_highest_monic_responses_raw2_by780':native_responses,
  'native_highest_monic_responses_after_final681_2_by681':highest,
  'all_native_monic_highest_responses_cancel':native_allzero,
  'native_u5_even_highest_responses':'same order-independent constant coefficient map as proven CACHED_COMPLETION u5_even2; all final681 responses0, hence all661 too',
  'other_unknown_proper_highest_core_responses_zero':'positive-tau core origins make scalar response of a core quadratic to an independent degree-n coefficient zero at nu0=0; n>=2',
  'highest_source_gate_NOT_a_numeric_X3_evaluation':True,
  'shared_cached_source_binding':{'new270_base':'QQ[z0,z1,z4,z5,r]/(m³,z1²), original8 CA-F9, sA=rA; NOT conductorH20',
   'new270_scope':'all needed proper/gauge second pairs in active set; no source3',
   'N2_B11':'literal actual N_C_B44, exact T/R12/N source equality; supplies r² and zir for i=0,1,4,5',
   'row20_X2_pivot_r_squared':'obtainable from these pinned actual N2_B11/D and cachedA1B1/nativeB11, no new source/gauge',
   'old619_shared_union':'r³,z0r²,z1r²,z4r²,z5r²; my candidate three mixed directions are subset, but no new job authorized here'},
  'remaining_first640_X3_input_identity':'A0_row20=0 and Bold0*u=0 remove Bold3 and A3; keep A1*Bold2*u+A2*Bold1*u before source/pivot correction',
  'remaining_old_Bold2_callbacks':'only r²,z0r,z4r,z5r contracted with actual A1; do not request whole old coefficient matrices',
  'remaining_pivot_carry_gate':'retain the projected81 first-carry adjoints from the literal minimal Schur operator, including15 mixed d3 sources; no full cancellation claimed',
  'no_new_X3_m2_source_gauge_jet_values_computed':True}
 assert nqq(data)<150000
 p=out/'actual661-first-corrections-m2-replacement-and-native-highest-sensitivities.json';blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;p.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'SOURCE_GATE_full661_m2_input_replaced_and_native_highest_sensitivities_tested',
  'all661_full44_source_first_corrections_verified':True,'separate_m2_input_eliminated_by18_pivot_second_callbacks':True,
  'native_highest_monic_sensitivities_all_zero':native_allzero,
  'pivot_first_higher_carries_or_all661_quadratic_values_closed':False,
  'no_new_source_gauge_jet_X3_m2_coefficients':True,'visible_QQ_support':nqq(data),'bytes':p.stat().st_size,
  'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()

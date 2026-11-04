#!/usr/bin/env python3
"""Cache/source/support preflight only; no new coefficient/rank worker."""
import hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_PREFLIGHT_20261104_v63'
PATHS={
 'honest':'reports/AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63/actual-honest-Tnil128-972-formal-operator.json',
 'honest_peer':'reports/AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63_REVIEW/report.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'radical_peer':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63_REVIEW/report.json',
 'native1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase5-native-adapter-firsts.json',
 'proper1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase4-native-source-and-R12.json',
 'native1_report':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/report.json',
 'native1_peer':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63_REVIEW/report.json',
 'properB11':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase3-actual-proper-G8-B11.json',
 'monicB11':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase1-source-B11-values.json',
 'B11_basis':'reports/AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63/actual-B11-source-and-gauge-preflight.json',
 'B11_peer':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63_REVIEW/report.json',
 'Bfirst':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
 'Bfirst_peer':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63_REVIEW/report.json',
 'nil':'reports/AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63/actual-E4-nil-two-generators-and-relations.json',
 'Koszul':'reports/AMBIENT_B18_REDUCED_CBAR_KOSZUL_INPUT_20261104_v63/actual-reduced-Cbar-three-sector-Koszul.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 ds={k:read(ROOT/p)for k,p in PATHS.items()}
 for k in ['honest_peer','radical_peer','native1_peer','B11_peer','Bfirst_peer']:assert ds[k]['all_passed']
 A,G,B,R=ds['proper1'],ds['native1'],ds['properB11'],ds['radical']
 assert len(A['kappa9'])==9 and all(len(v)==36 for v in A['kappa9'])
 assert len(G['G8'])==8 and all(len(v)==36 for v in G['G8'])
 assert len(G['N'])==12 and all(len(r)==12 and all(len(v)==36 for v in r)for r in G['N'])
 assert len(B['kappa9_C_B44'])==9 and all(len(v)==44 for v in B['kappa9_C_B44'])
 assert len(B['actualG8_C_B44'])==8 and all(len(v)==44 for v in B['actualG8_C_B44'])
 assert B['all_constants_and_firsts_match_frozen_native128']
 coords=['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA']
 old=ds['Bfirst'];assert old['coordinate_order']==['constant']+coords
 col=old['source_columns640_rows56_sparse_dual9'];assert len(col)==640
 stats={'oldB_entries':0,'oldB_constants':0,'oldB_firsts':0,'P_entries':0,'P_constants':0,'P_firsts':0}
 for c in col:
  for i,v in c:
   p='oldB'if i<48 else'P';stats[p+'_entries']+=1
   stats[p+'_constants']+=sum(k==0 for k,q in v);stats[p+'_firsts']+=sum(k>0 for k,q in v)
 assert stats=={'oldB_entries':155,'oldB_constants':4,'oldB_firsts':245,'P_entries':105,'P_constants':0,'P_firsts':160}
 native=R['native_tensor_Koszul6_HB2_resolution'];assert native[1]['shape']==[9,35]
 syz_support=sum(len(c)for c in native[1]['columns']);assert syz_support==70
 Dterms=sum(len(p)for row in R['G8_polynomial_factor_D8_by9']for p in row);assert Dterms==54
 # Bound sparse coefficients before computing a single actual scalar/rank.
 xbound=36*(9*stats['oldB_constants']+stats['oldB_firsts'])+stats['P_firsts']+syz_support*4*4*9
 assert xbound==20356
 bounds={'new_X0_and_X1_nonzero_output_QQ_upper':xbound,'m0_24_by44_dense':1056,
  'm0_selected_minor_and_inverse_twice':1152,'X0_pivot_vectors44_by20':880,
  'X0_constant_pivot_coefficients20_by760':15200,'W0_84_by780_upper':65520,
  'W0_fixed_pivot_minor84_squared':7056,'simultaneous_C36_N_Ninverse_A_slots':14256,
  'C36_C44_proper_variables_and_dR_values':3000,'native_Koszul_sparse_linear_slots':336}
 assert sum(bounds.values())==128812 and sum(bounds.values())<150000
 plan={'actual_coefficient_frame':coords,'ell0':'same actual CA-F9 sA=rA; not physical rho approximation',
  'source_coefficient_availability':{
   'all9_proper_kappa_full8_firsts':'proper1.kappa9, C36 interleaved four trace components',
   'all5_actual_u_full8_firsts':'native1.G8[0:5], not rational pivot roots on the intersection',
   'native_N12_full8_firsts':'native1.N; invert the actual unit matrix in C36',
   'third_gamma_sector_r_squared':'properB11.kappa9_C_B44 and actualG8_C_B44 plus actual monics/B11 multiplication',
   'old640_B_and_P_all8_firsts':'Bfirst source_columns640_rows56_sparse_dual9',
   'radical_native35_source':'literal Koszul6 tensor HilbertBurch2; 70 nonzero polynomial entries',
   'missing_source_directions':[]},
  'coefficient_construction':[
   'Evaluate actual dR9 and D8x9 in C36 using u5=G8[:5], core=(kappa[0],kappa[1],kappa[2],kappa[6]) and pinned GL4.',
   'Preserve tauA^2=nuA_1 and tauB^2=nuB_1 in C36; do not freeze both squares to zero.',
   'Evaluate dR9 in C44 to extract only third gamma rA^2; B11 retains every required monomial.',
   'Vrad=[Qrad non-even3, -beta8*tauA non-even3], beta=zA7,gA; gA has order3 and hence zero firsts/rA^2.',
   'm0 consists of gamma-sector z firsts and r firsts for sectors0,1; sector2 has rA^2 coefficient.',
   'A12x9=N12_inverse*injection12x8*D8x9; Res4(A^t) transforms old Bfirst, with unchanged P8.',
   'Last140 X columns are the exact restricted native Syz(dR9), evaluated in C36.'],
  'conditional_numeric_minimalization':[
   'First check rank(m0)=24. Select literal24 columns E, put n=E*(mE)^-1. Then mn=I and Vgamma*n=d1.',
   'Compute r=rank(X0), select pivot source columns; residual constant combinations c satisfy X0*c=0.',
   'For each residual, (mX*c)_1=m0*(X1*c), so unknown m1 and regular-quotient first derivatives cancel.',
   'Extract the unique84 constant Koszul lift W0 from coefficients zi in component j (including pairs zi,rA^2). Verify all24x8 derivative equations.',
   'mu_S(Tnil)=104-r-rank(W0). A constant rank is a minimal-generator result, not total length or nil-generator coefficients.'],
  'bounds':bounds,'total_saved_plus_live_QQ_bound':sum(bounds.values()),
  'C36_matrix_inverse_C_product_bound':3*12**3,'C36_A_matrix_product_bound':12*8*9,
  'worker_cpu_seconds':60,'worker_address_bytes':2<<30,'total_output_bytes':32<<20,
  'unknown_quotient_derivatives_assumed_zero':False,'new_source_jet_or_value_worker_run':False}
 ev=out/'actual-cache-source-and-minimalization-preflight.json';ev.write_text(json.dumps(plan,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'PREFLIGHT_ONLY_all_required_cached_source_directions_present',
  'actual_m0_X0_X1_or_ranks_computed':False,'missing_source_directions':[],
  'existing_sparse_Bfirst_support':stats,'source_visible_QQ_bound':sum(bounds.values()),
  'new_source_solve_or_jet_or_rank_worker':False,'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ev.name:sha(ev)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Exact H20->B11 source quotient and bounded native-S2 source preflight.

Only saved finite algebra values are mapped.  No gauge correction or source
jet is solved; future scalarB gauge30 corrections need150QQ coefficients.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63'
PATHS={
 'H20':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'H20_peer':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63_REVIEW/report.json',
 'source105':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63/actual-source105-Artin20-values.json',
 'source105_peer':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63_REVIEW/report.json',
 'gauge1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase3-gauge30-firsts.json',
 'J30':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase2-fixed30-Jacobian.json',
 'first_preflight':'reports/AMBIENT_B18_NATIVE128_FIRST_SOURCE_PREFLIGHT_20261104_v63/actual-native128-first-input-and-mixed-index-profile.json',
 'L2_cancel':'reports/AMBIENT_B18_T_NATIVE_ADAPTER2_CANCELLATION_PREFLIGHT_20261104_v63/report.json'}
SURVIVE_Z=[0,1,4,5]
BASIS=[(0,)*8,tuple([0]*7+[1]),tuple([0]*7+[2])]
for i in SURVIVE_Z:BASIS.append(tuple(int(j==i) for j in range(8)))
for i in SURVIVE_Z:BASIS.append(tuple(int(j==i or j==7) for j in range(8)))
INDEX={m:i for i,m in enumerate(BASIS)}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rec(a):return {int(i):Q(c) for i,c in a if Q(c)}
def enc(a):return [[i,str(c)] for i,c in sorted(a.items()) if c]
def add(a,b,c=Q(1)):
 a=dict(a)
 for i,v in b.items():
  z=a.get(i,Q(0))+c*v
  if z:a[i]=z
  else:a.pop(i,None)
 return a
def bm(a,b):
 z={}
 for i,c in a.items():
  for j,d in b.items():
   m=tuple(x+y for x,y in zip(BASIS[i],BASIS[j]));k=INDEX.get(m)
   if k is not None:
    v=z.get(k,Q(0))+c*d
    if v:z[k]=v
    else:z.pop(k,None)
 return z
def bpower(a,n):
 v={0:Q(1)}
 for _ in range(n):v=bm(v,a)
 return v
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 (out/'phase.json').write_text(json.dumps({'phase':'source_ring_preflight_inputs_persisted','inputs':inputs},indent=2)+'\n')
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 H,S,G,J,P=D['H20'],D['source105'],D['gauge1'],D['J30'],D['first_preflight']
 assert D['H20_peer']['all_passed'] and D['source105_peer']['all_passed'] and D['L2_cancel']['all_passed']
 hb=H['quotient_basis_monomials'];assert hb==S['A20_basis'] and len(hb)==20
 mp=[{INDEX[tuple(m)]:Q(1)} if tuple(m) in INDEX else {} for m in hb]
 lifts=[hb.index(list(m)) for m in BASIS];assert lifts==[0,8,16,1,2,5,6,9,10,13,14]
 actions=[[rec(c) for c in v] for v in H['multiplication8_column_sparse']]
 def haction(i,v):
  z={}
  for j,c in v.items():z=add(z,actions[i][j],c)
  return z
 def hprod(i,j):
  v={i:Q(1)}
  for k,n in enumerate(hb[j]):
   for _ in range(n):v=haction(k,v)
  return v
 def project(v):
  z={}
  for i,c in v.items():z=add(z,mp[i],c)
  return z
 for i in range(20):
  for j in range(20):assert project(hprod(i,j))==bm(mp[i],mp[j])
 coord=[mp[hb.index([int(j==i) for j in range(8)])] for i in range(8)]
 assert all(not bm(coord[i],coord[j]) for i in range(7) for j in range(7))
 assert all(not bm(coord[i],bpower(coord[7],2)) for i in range(7))
 assert not bpower(coord[7],3)
 unknowns=[[int(i),enc(project(rec(v)))] for i,v in zip(S['joint105_unknown_order'],S['all105_unknown_final_values'])]
 requested=[[int(i),enc(project(rec(v)))] for i,v in S['all_requested_original_source_root_values']]
 nus=[project(rec(v)) for v in S['nu_CA_CB_values']]
 assert all(not v.get(0,Q(0)) for v in nus)
 assert all(not bm(a,b) for a in nus for b in nus)
 # Build the actual C_B multiplication table with both true source monics.
 # No unknown second adapter is evaluated here.  A matrix product can reuse
 # this sparse table in the separately authorized prospective worker.
 def cm_basis(i,j):
  bi,ti=divmod(i,4);bj,tj=divmod(j,4);v=bm({bi:Q(1)},{bj:Q(1)});mask=ti^tj;carry=ti&tj
  for k in range(2):
   if carry&(1<<k):v=bm(v,nus[k])
  return {4*b+mask:c for b,c in v.items()}
 table=[[cm_basis(i,j) for j in range(44)] for i in range(44)]
 table_support=sum(len(v) for row in table for v in row)
 def cm(a,b):
  z={}
  for i,c in a.items():
   for j,d in b.items():z=add(z,table[i][j],c*d)
  return z
 # All44^3 associativity products are finite table checks, not source jets.
 for i in range(44):
  for j in range(44):
   for k in range(44):assert cm(table[i][j],{k:Q(1)})==cm({i:Q(1)},table[j][k])
 source_evaluation={'closure_count':len(P['dual9_source_evaluation']['closure_order']),
  'same_full150_tracepair_roots_as_first_worker':True,
  'base_B11_values':'map all105 saved H20 roots; no implicit source solve',
  'maximum_live_base_coefficients':11*3600,
  'mapped105_saved_coefficients_upper':105*11,
  'mapped164_saved_coefficients_upper':164*11,
  'physical300_trace_coefficients_upper':300*11,
  'C_B44_physical150_tensor_upper':150*44,
  'C_B44_matrix5_by5_upper':25*44,
  'simultaneous_matrix_tensor_scratch_upper':10*25*44+3*150*44,
  'gauge30_scalar_B_unknown_quadratic_coefficients':30*5,
  'redundant_C4_baseextended_gauge30_unknown_quadratic_coefficients_upper':30*5*4,
  'fixed_J30_inverse_existing_QQ_slots':30*30,
  'new_Jacobian_or_Jacobian_second_coefficient_required':False,
  'new_physical_source_or_105_solve_required':False,
  'native_R12_T_N_second_calculation_required':False}
 qqguard=(11*3600+105*11+164*11+300*11+10*25*44+3*150*44+30*5*4+900+table_support)
 assert qqguard<150000
 ev={
  'B11_variable_order8':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA'],
  'B11_exact_ideal':'(rA^3, all z_i*z_j, all z_i*rA^2, zA2,zA3,zA6)',
  'B11_basis_monomials8':BASIS,'B11_basis_names':['1','r','r2','z0','z1','z4','z5','z0r','z1r','z4r','z5r'],
  'H20_basis_to_B11_columns':list(map(enc,mp)),'B11_basis_lifts_in_H20':lifts,
  'all400_H20_products_intertwined':True,'actual_conductor_product_I_maps_zero':True,
  'actual_I_zero_source_proof':'Use accepted all-order conductor pair IA=(z7,r^3), IB=(z7+c*r^2,r^3), with actual unit normalizations. Every crossproduct maps0; all higher plane terms are killed by r^3 and all z*r^2.',
  'all105_actual_B11_root_values':unknowns,'all164_requested_actual_B11_source_values':requested,
  'nu_A_B_actual_B11':list(map(enc,nus)),
  'nu_A_square_nu_A_nu_B_nu_B_square_all_zero':True,
  'C_B44_basis':'B11_basis in given order, each multipliedby(1,tauA,tauB,tauA*tauB)',
  'C_B44_multiplication_basis_products':[[enc(v) for v in row] for row in table],
  'C_B44_table_support':table_support,'all85184_C_B44_associativity_checks':True,
  'C_B44_unit_inverse_safe_geometric_nilpotence_guard':9,
  'nilpotence_proof':'m_B^3=0; tauA^2,tauB^2 in m_B. Any product of9 maximal-ideal generators contains at least3 base factors after square reduction, hence0. The bound is conservative.',
  'gauge_source_convention':{
   'actual_unknowns':'15even scalarB coefficients and15odd scalarB coefficients of AUXILIARY etaA',
   'actual_equations':'15 actual CA gauge/source constraints in B11[etaA]/etaA^2=nuA, split into30 scalarB equations',
   'old_gauge_first9':'same frozen phase3 even15_first9/odd15_first9, mapped toB11',
   'unknowns150':'only five quadratic base monomials for each of30 scalar functions',
   'aux_eta_vs_physical_tau':'after solving scalar30 system, embed same coefficients into C_B44 then evaluate etaA=physical tauA',
   'fixed_J30':'same frozen origin30_by30 inverse; correction lies m_B², which is killed bym_B, so all higher baseJacobian terms givezero',
   'redundant600':'baseextend the scalar30 system byC4 ifdesired; valid but notnecessary',
   'source_first_constraint_gate':'all known constant/linear residuals must be checked zero before the quadratic correction; do not zero them by convention'},
  'prospective_source_evaluation_profile':source_evaluation,'combined_live_and_saved_QQ_guard':qqguard,
  'prospective_worker_goal':'only gauge second correction, actual proper kappa/G8 and nativeS2, then L0S2+L1S1; no N2/T2/R12 second adapter',
  'prospective_worker_not_run':True,'no_new_physical_jet_source_solve_gauge_correction_or_X2_evaluation':True,
  'source_manifest':inputs}
 ep=out/'actual-B11-source-and-gauge-preflight.json';blob=json.dumps(ev,separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;ep.write_text(blob)
 report={'id':STEM,'all_passed':True,'actual_H20_surjects_B11':True,'all400_products_verified':True,
  'all105_saved_roots_mapped_without_solve':True,'actual_C_B_dimension':44,'C_B_table_associativity_checks':44**3,
  'actual_scalar_gauge_new_QQ_unknowns':150,'redundant_bound600_optional':True,'combined_QQ_guard':qqguard,
  'gauge_or_source_second_worker_executed':False,'X2_T_length_depth_computed':False,
  'input_hashes':inputs,'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()

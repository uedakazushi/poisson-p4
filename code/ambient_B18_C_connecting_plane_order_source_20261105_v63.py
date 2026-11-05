#!/usr/bin/env python3
"""Actual 56-column plane-order/jet routing; preserve the surviving sigma2 adjoint."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_CONNECTING_PLANE_ORDER_SOURCE_20261105_v63'
PATHS={
 'five':'reports/AMBIENT_B18_C_CONNECTING_FIVE_FORMAL_SOURCE_20261105_v63/actual-five-head-formal-RegDivell-operator-and-source-binding.json',
 'five_report':'reports/AMBIENT_B18_C_CONNECTING_FIVE_FORMAL_SOURCE_20261105_v63/report.json',
 'q7_peer':'reports/AMBIENT_B18_C_GENUINE_PROJECTION7_SOURCE_20261105_v63_REVIEW/report.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'C2_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
 'ell':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63/all164-cached-pure-ell-readbacks-and-physical12-second-values.json',
 'ell_peer':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63_REVIEW/report.json',
 'frame':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63/actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json',
 'frame_peer':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63_REVIEW/report.json',
 'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'chain_peer':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_20261104_v63_REVIEW/report.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'C0_peer':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json',
 'CA':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 started=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={
  'five':['full_V56_reference_pairs_columns56_targets3','required_raw_numerator_coefficient_functionals','old_five_heads_free56_rows'],
  'C2':['physical12_coefficients'],'raw33':['actual33_basis'],
  'ell':['all12_actual_physical_coefficients_pure_ell3'],
  'frame':['all12_verified_source_records','both_nu_constant_and_pureplane_first_zero'],
  'chain':['operator_scalar_nodes','old_actual_CA_Q0_multipliers9_C4','exact_existing_U9_C4','accepted_actual_regular_function_identities'],
  'C0':['actual_C0_basis8_scalar_nodes','all35_regular_coordinate_divisions','cached_actual_CA_gA_all_degree3_coefficients'],
  'CA':['families'],
 }
 # Only the first ten chain scalar nodes are needed for C0 U9. Do not
 # decode its entire scalar array or the unrelated 640 boundary columns.
 chain_fields=fields.pop('chain');d={};family_parse_count=0
 for k,p in PATHS.items():
  if k=='chain':continue
  d[k]=S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k=='CA':
   family_parse_count=S.count(d[k]);assert family_parse_count<20000
   d[k]={'actual_zero_odd_I0_generators':d[k]['families']['CA']['actual_zero_odd_I0_generators']}
  if k.endswith('_peer'):assert d[k]['all_passed']
 import ambient_B18_C_connecting_five_formal_source_20261105_v63 as F
 c=S.selective(ROOT/PATHS['chain'],chain_fields[1:],receipts)
 prefix=F.scalar_prefix(ROOT/PATHS['chain'],9,receipts)
 basis=d['C0']['actual_C0_basis8_scalar_nodes'];assert basis==[7148,7199,7454,7505,7688,7943,8045,505938]
 assert[prefix[i]for i in range(2,10)]==[['actual_C0_scalar_function',j]for j in basis]
 assert c['exact_existing_U9_C4']==[[1,0,0,0]]+[[0,j+2,0,0]for j in range(8)]
 I0=d['CA']['actual_zero_odd_I0_generators'];assert[i for row,i in I0[:7]]==basis[:7]
 assert I0[7]==[63,7790]and basis[-1]!=I0[7][1]
 tail=d['C0']['all35_regular_coordinate_divisions'][0]
 assert tail['original_pure_CA_source_function']==7790 and tail['restriction_r_s_node']==505938
 tail_identity=next(r for r in c['accepted_actual_regular_function_identities']if r['kind']=='accepted_actual_CA_f63_tail')
 assert tail_identity['function']==7790 and tail_identity['gA']==505938
 assert tail_identity['coordinate_quotients']==[r['node']for r in tail['all7_regular_coordinate_divisions']]
 gA=list(map(Q,d['C0']['cached_actual_CA_gA_all_degree3_coefficients']));assert len(gA)==10 and not any(gA[:6])and any(gA[6:])
 old_alias={r['physical_odd_coefficient']:r['odd_generator']for r in c['accepted_actual_regular_function_identities']if r['kind']=='exact_QQ_linear_source_alias'}
 assert d['frame']['both_nu_constant_and_pureplane_first_zero']
 raw=d['raw33']['actual33_basis'];C2=d['C2']['physical12_coefficients'];ell={r['physical_row']:r for r in d['ell']['all12_actual_physical_coefficients_pure_ell3']}
 assert len(C2)==12
 ri,ei=raw.index([0,1,0]),raw.index([0,0,1]);rr,re=raw.index([0,2,0]),raw.index([0,1,1])
 B2=[[Q(0)]*56 for _ in range(3)];H2=[[Q(0)]*56 for _ in range(3)];checks=[]
 for j,rec in enumerate(C2):
  row=rec['physical_row'];data={k:{int(i):Q(q)for i,q in rec[k]}for k in['delta','bA','bB']}
  for k in data:assert all(not data[k].get(i,Q(0))for i in[0,ri,ei])
  for k,ek in [('delta','delta_jet3'),('bA','bA_jet3')]:
   assert not data[k].get(rr,Q(0))and not data[k].get(re,Q(0))and Q(ell[row][ek][2])==0
  bb=[data['bB'].get(rr,Q(0)),data['bB'].get(re,Q(0)),Q(ell[row]['bB_jet3'][2])]
  for k in range(3):B2[k][4*j]=-bb[k];H2[k][4*j+1]=-bb[k]
  if j<7:assert old_alias[rec['roots']['bA_node']]==basis[j]
  checks.append({'physical_row':row,'all_delta_bA_bB_plane_orders0_1_zero':True,'delta_bA_plane2_zero':True,'bB_plane2':list(map(str,bb))})
 support=sorted({j for row in B2 for j,q in enumerate(row)if q});assert support==[0,16,20]
 C=[B2[2][j]-Q(96,275)*B2[1][j]for j in range(56)]
 assert any(C)and all(not q for j,q in enumerate(C)if j not in support)
 # List routing in the exact selected q7 functionals. sigma lives in Rbar,
 # so its exponent at ell is ZERO; V1 with a positive ell factor is zero.
 routed=[];required_sigma=set();dropped=[]
 for req in d['five']['required_raw_numerator_coefficient_functionals']:
  M=tuple(req['numerator_monomial']);t=req['raw_target'];paths=[]
  for sm in S.mons(sum(M)):
   if sm[8]or any(a>b for a,b in zip(sm,M)):continue
   vm=tuple(b-a for a,b in zip(sm,M));order=sum(vm)
   if not order or order==1:
    dropped.append([t,list(M),list(sm),'V0=0 or V1_ell=0']);continue
   if order==2 and not any(vm[:7])and t==0:
    dropped.append([t,list(M),list(sm),'V2_A,plane=0']);continue
   required_sigma.add(sm);paths.append({'sigma_monomial':list(sm),'V_monomial':list(vm),'V_order':order})
  routed.append({'raw_target':t,'numerator_monomial':list(M),'retained_product_paths':paths})
 expected={S.Z,S.mon(r=1),S.mon(z0=1),S.mon(z4=1),S.mon(z5=1),S.mon(r=2)}
 assert required_sigma==expected,(required_sigma,expected)
 assert not any(sum(m)>=3 for m in required_sigma)
 evidence={
  'scope':'Source-only actual full56 plane-order and exact required sigma-jet routing. No connecting coefficient or newly requested cached response is evaluated.',
  'namespace_correction_to_frozen_five_formal_prose':{
   'superseded_field':'old_eighth_generator in FIVE_FORMAL evidence; corresponding note sentence is overbroad',
   'correct_actual_U9_C0_basis':'U9 honest C0 free9 uses zA0..6,gA505938, exactly as its 168 reference pairs already encode.',
   'correct_OLD640_multiplier':'OLD Q0 relations use physical f63=7790, preserved by the actual identity f63=gA505938+Σe_i*zAi. These are distinct functions; the chain old_multiplier_last_is_f63 flag refers only to OLD640 multipliers.',
   'actual_full_V56_refs_signs_sigma_RegDiv_no_change':True,
   'actual_C0_basis8_original_nodes':basis,'physical_f63_original_node':7790,'tail_gA_original_node':505938,
   'exact_tail_identity_record':tail_identity,
  },
  'C0_U9_origin_plane1_plane2_zero_proof':'First7 Uodd are the exact accepted CA physical odd source aliases and have plane order≥3 by the literal bA2/first/constant checks. Eighth Uodd is the actual zero7 tail gA with its first6 bivariate coefficients zero; it has plane order3. No f63 versus gA substitution is made in the actual U9 matrix.',
  'all12_literal_low_order_checks':checks,
  'universal_C4_odd_components_per_Fj_basis_multiplier':[
    ['bA','-bB','0'],['delta','0','-bB'],['0','delta','bA'],['-nuB*bB','nuA*bA','delta']],
  'universal_derivation':'Fj=delta+bA*tauA−bB*tauB; tauA²=nuA,tauB²=nuB. AllF/U coefficients vanish at0, their plane1 vanish, nu0/plane1 vanish, and delta2/bA2_plane vanish. Therefore ALL56 V0=0 and V1_plane=0; ALL56 V2_A,plane=0. V2_B andV2_H are the literal matrices below. Full Q48 frame is the independently checked actual source, not an equivalent substitute.',
  'actual_plane2_B_three_monomials_rows3_columns56':[[str(q)for q in r]for r in B2],
  'actual_plane2_H_three_monomials_rows3_columns56':[[str(q)for q in r]for r in H2],
  'plane2_monomial_order':['r²','r*ell','ell²'],
  'sigma2_surviving_B_high_coordinate_adjoint':[[j,str(q)]for j,q in enumerate(C)if q],
  'sigma2_scope':'The B-high functional N_(r²ell²)−b1*N_(r³ell) retains (V2_B,ell²−b1*V2_B,r ell)*sigma2_rr. It is nonzero on source rows0,16,20 and is NOT dropped. The A-high has no sigma2 term. Jordan coordinates may combine these, but no further cancellation is assumed.',
  'all_selected_q7_product_paths':routed,'formally_dropped_paths_with_literal_reason':dropped,
  'only_required_sigma_jets':['sigma0','sigma1_r','sigma1_z0','sigma1_z4','sigma1_z5','sigma2_rr'],
  'sigma3_drops_exactly':True,
  'highest_physical_frontier':'Projected V4_A/B*(sigma0) at r³ell,r²ell², together with coherent V3*sigma1, V2 mixed*sigma1 and the retained B2*sigma2_rr. Old r/g source functions stay regular/opaque. No blanket source4 or sigma2 cache availability is asserted.',
  'no_new_IFT_gauge_physical_jet_edge_or_H0_value':True,
  'selective_input_receipts':receipts,
 }
 owned=S.count(d)+S.count(c)+S.count(prefix)+S.count(evidence)+20000+family_parse_count;assert owned<150000
 evidence['whole_source_owned_numeric_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'all56-plane-orders-surviving-sigma2-adjoint-and-U9-namespace-correction.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_ALL56_PLANE_SIGMA3_ZERO_SIGMA2_ONE_ADJOINT_RETAINED',
  'actual_U9_and_old_Q0_eighth_namespace_distinguished':True,'sigma3_needed':False,'sigma2_rr_adjoint_retained':True,
  'actual_connecting_values_computed':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{ep.name:sha(ep)},
  'whole_numeric_upper':owned,'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Independent four-graph/RHS/six-adjoint input audit, without solve or M1 query."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,itertools,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale,mv
from ambient_B18_ell2_T2_mu_review_20261105_v63 import full_G2
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_LINEAR_GRAPH4_SOURCE_20261105_v63';STEM=A+'_REVIEW'
P={'primary':f'reports/{A}/report.json','gate':f'reports/{A}/actual23-linear-graph-four-common-contractions-and-six-adjoint-input.json',
'transport':'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63/actual-active-unit-bases-and26-frame-transport.json',
'active':'reports/AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_PREFLIGHT_20261105_v63/active-five-d3-source-split-and-budget.json',
'resolution':'reports/AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63/actual-C2-resolution-and11-unit-source-cancellation.json',
'G2':'reports/AMBIENT_B18_ELL2_C2_G2_LINEAR_SOURCE_PROFILE_20261105_v63/actual-G2-linear-fixed-source-profile.json',
'G2_peer':'reports/AMBIENT_B18_ELL2_C2_G2_LINEAR_SOURCE_PROFILE_20261105_v63_REVIEW/report.json',
'source_peer':'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63_REVIEW/report.json',
'known':'reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json',
'mini':'reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json',
'adjoint':'reports/AMBIENT_B18_ELL2_T2_FIRST_ADJOINT52_20261105_v63/actual-two-full918-source-adjoint-certificates-or-obstructions.json',
'adjoint_peer':'reports/AMBIENT_B18_ELL2_T2_FIRST_ADJOINT52_20261105_v63_REVIEW/report.json',
'gauge':'reports/AMBIENT_B18_ELL2_T2_FIRST_GRAPH_COMPARISON_GAUGE_PREFLIGHT_20261105_v63/report.json',
'HF66_peer':'reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_RANK_20261105_v63_REVIEW/report.json',
'primitive':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/phase2-all105-equations-and-primitive56-source-columns.json',
'primitive_peer':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'GRAPH_SOURCE_AND_ADJOINT_INPUT_ONLY','inputs':inputs});D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert D['primary']['all_passed'] and all(D[k]['all_passed'] for k in D if k.endswith('peer'))
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 E=D['gate'];T=D['transport'];M=list(map(sp,T['actual_Mnew0_columns']));J=T['same_selected26_V_columns'];factors=list(map(sp,T['all56_new_column_factorizations_in_same_selectedJ_basis']))
 W=[add({j:Q(1)},scale({J[i]:c for i,c in factors[j].items()},-1)) for j in D['mini']['remaining_free3_originalV_rows']];v={0:Q(-3375,11),1:Q(-50),2:Q(25,2)};w=mv(W,v);assert w=={23:Q(-3375,11),29:Q(-50),30:Q(25,2)} and not mv(M,w)
 assert list(map(enc,W))==E['H0_free3_columns'] and enc(w)==E['actual_H0v_source56']
 bottom=list(map(sp,D['known']['all642_bottom27_columns']));selected=D['known']['selected_bottom_columns'];patterns=[]
 for k,(j,record) in enumerate(zip(selected,E['all23_linear_graph_source_records'])):
  assert record['bottom_basis_column']==k and record['source642_column']==j and record['original768_combination']==D['mini']['all642_zero_constant_combinations'][j]['original768_combination'] and enc(bottom[j])==record['known_b1_3_times9']
  br={i:bottom[j].get(i*9+7,Q(0)) for i in range(3) if bottom[j].get(i*9+7,Q(0))};be={i:bottom[j].get(i*9+8,Q(0)) for i in range(3) if bottom[j].get(i*9+8,Q(0))};cr=br.get(0,Q(0))/v[0];ce=be.get(0,Q(0))/v[0];assert br==scale(v,cr) and be==scale(v,ce)
  assert Q(record['four_unknown_mu_coefficients_per_channel']['M_r_H0v'])==ce and Q(record['four_unknown_mu_coefficients_per_channel']['M_ell_H0v'])==cr
  if cr or ce:patterns.append([k,j,str(cr),str(ce)])
 assert patterns==[[19,582,'1','14/125'],[21,634,'8/3','32/75']]
 # New source inclusion has all constants already peered; differentiate
 # only the canonical variable part of periodic/Koszul d1.
 eta={i:j for j,i in enumerate(D['resolution']['physical_target3_basis_positions'])};G1=[]
 for col in map(sp,T['new_full_I1_origin26_columns']):
  p={}
  for k,c in col.items():
   block,row=divmod(k,14)
   if row in eta:p=add(p,{eta[row]*9+(8 if block==0 else block-1):c})
  G1.append(p)
 assert len(G1)==26 and list(map(enc,G1))==E['actual_G1linear3_by26_columns']
 G2=full_G2(T,D['active'],D['resolution'],D['G2']);mons=[p for p in itertools.combinations_with_replacement(range(9),2) if p!=(8,8)];mi={p:i for i,p in enumerate(mons)}
 for col in G2:
  p={}
  for index,c in col.items():
   row,var=divmod(index,9)
   for k,a in G1[row].items():
    target,u=divmod(k,9);m=tuple(sorted((u,var)))
    if m in mi:p=add(p,{target*44+mi[m]:a*c})
  assert not p
 lam=[{i//2:Q(c) for i,c in r['lambda52']} for r in D['adjoint']['results']];assert all(i%2 for r in D['adjoint']['results'] for i,c in r['lambda52']);invariants=[{},{}]
 for i,c in lam[0].items():invariants[0]=add(invariants[0],{i*9+7:c});invariants[1]=add(invariants[1],{i*9+8:c})
 for i,c in lam[1].items():invariants[1]=add(invariants[1],{i*9+7:-c})
 assert list(map(enc,invariants))==E['two_invariant_functionals_on_all234_F1_coefficients']
 pairs=[(7,7),(7,8)];adjcols=[{} for _ in range(6)]
 for j,col in enumerate(G1):
  for index,c in col.items():
   row,u=divmod(index,9)
   for d in range(9):
    p=tuple(sorted((u,d)))
    if p in pairs:adjcols[2*row+pairs.index(p)]=add(adjcols[2*row+pairs.index(p)],{j*9+d:c})
 assert list(map(enc,adjcols))==E['six_rr_rell_source_adjoint_columns_in234_domain']
 primitive=D['primitive'];basis=list(map(tuple,primitive['basis56']));assert len([m for m in basis if sum(m)==2])==44
 assert primitive['all105_equations_exact'] and primitive['exact_allorder_primitive_frame_binding'] and len(primitive['all56_primitive_V_columns'])==56
 assert E['primitive_V2_input_provenance']['path']==P['primitive'] and E['primitive_V2_input_provenance']['source_combination']==enc(w)
 assert E['whole_six_adjoint_slot_bound']==10817<150000 and E['six_adjoint_solve_not_started']
 ep=o/'independent-four-scalar-graph-and-full234-adjoint-input.json';dump(o,ep.name,{'inputs':inputs,'all23_original_source_graph_columns_and_two_common_patterns_verified':True,'actual_Wv_fixed_source23_29_30_and_M0Wv0_verified':True,'actual_G1linear26_from_peered_source_inclusion_and_all102_G1G2_zero_verified':True,'two_invariant_covectors_and_all_six234_coefficient_adjoint_columns_verified':True,
'exact_RHS_proof':'G1(0)=0 and M0*Wv=0. Expanding the fullformalidentityG1*M=V on the FIXED constant56vectorWv throughdegree2 givesG1_1*(M1Wv)=(VWv)_2. TermsG1_0*M2 andG1_2*M0Wv vanish separately; no change from rA toell trace is used.',
'adjoint_required_before_values':'Only an explicit rho withrho*G1linear=invariant onALL234 F1slots permitsinvariant=rho*(VWv)_2. Checkingonlysixdesiredoutputslotsisnot sufficient.',
'coherent_gauge_scope':'Changes(0,u,u,v) removeonlytwogaugecoordinates. ThecanonicalfourM1valuesremainuncomputed; source/targetunitnormalizationisnotfourarbitraryzeroassignments.',
'primitive_RHS_all44_quadratic_coefficients_existing_source_bound':True,'new_physical_or_M1_or_adjoint_values_computed':False,
'coarse_closure_limit':'Discardingall23graphrelationsleavesanabstractbottom(B/I²)^3 ofdimensionatleast99. ThisisnotanactualT2lowerbound. No unstatedidentificationofnonhomogeneousI²withitsinitialidealisusedhere.',
'whole_six_adjoint_slot_bound':10817})
 report={'id':STEM,'all_passed':True,'graph4_source_and_full234_six_adjoint_gate_accepted':True,'actual_M1_or_invariant_or_adjoint_values_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json'),'seconds':report['seconds']}))
if __name__=='__main__':main()

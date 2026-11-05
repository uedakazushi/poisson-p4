#!/usr/bin/env python3
"""Next coherent HEAD29 source, retaining exactly three plane4 callbacks."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
from ambient_B18_C_connecting_head23_source_20261105_v63 import add,sc,mul,const,enc,sum_poly
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_CONNECTING_HEAD29_SOURCE_20261105_v63'
PATHS={
 'q7':'reports/AMBIENT_B18_C_GENUINE_PROJECTION7_SOURCE_20261105_v63/sparse-actual-Dlinear-projector7-source-and-five-head-order-frontier.json',
 'q7_peer':'reports/AMBIENT_B18_C_GENUINE_PROJECTION7_SOURCE_20261105_v63_REVIEW/report.json',
 'unit':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
 'unit_peer':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63_REVIEW/report.json',
 'head23':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_SOURCE_20261105_v63/full-coherent-head23-two-socle-cache-polynomials-and-low-two-row-source.json',
 'head23_peer':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_SOURCE_20261105_v63_REVIEW/report.json',
 'plane':'reports/AMBIENT_B18_C_CONNECTING_PLANE_ORDER_SOURCE_20261105_v63/all56-plane-orders-surviving-sigma2-adjoint-and-U9-namespace-correction.json',
 'adjoint':'reports/AMBIENT_B18_C_CONNECTING_SIGMA2_ADJOINT_SOURCE_20261105_v63/actual-three-row-sigma2-source-adjoint-and-coherent-fourth-frontier.json',
 'adjoint_peer':'reports/AMBIENT_B18_C_CONNECTING_SIGMA2_ADJOINT_SOURCE_20261105_v63_REVIEW/report.json',
 'five':'reports/AMBIENT_B18_C_CONNECTING_FIVE_FORMAL_SOURCE_20261105_v63/actual-five-head-formal-RegDivell-operator-and-source-binding.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'lambda':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63/phase1-six-root-existing-cache-readbacks.json',
 'new_lambda':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63/phase1-single-root-existing-cache-readback.json',
 'value_peer':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 st=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={
  'q7':['kappa_k_A_B_UNPAIRED','projector7_by660_sparse_columns','beta2_plane_after_actual_s_equals_r_plus_ell_over5'],
  'unit':['fixed_selected24_columns','fixed24_minor_inverse'],
  'head23':['sigma1_z_complete_leftadjoints_Rz_n0'],
  'plane':['actual_plane2_B_three_monomials_rows3_columns56','actual_plane2_H_three_monomials_rows3_columns56'],
  'adjoint':['actual_sigma0_56_by5_sparse','fixed_L_Cn0_row24_sparse'],
  'five':['full_V56_reference_pairs_columns56_targets3'],
  'C2':['physical12_coefficients'],'lambda':['all3_lambda_readbacks'],
  'new_lambda':['requested_bB63_r_ell2'],'raw33':['actual33_basis'],
 }
 d={}
 for k,p in PATHS.items():
  d[k]=S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k.endswith('_peer'):assert d[k]['all_passed']
 assert d['adjoint']['actual_sigma0_56_by5_sparse'][2]==[[29,'1']]
 assert d['five']['full_V56_reference_pairs_columns56_targets3'][29][1]==[0,1]
 inputs={};loc={}
 def var(k,value=None,locator=None):
  k=tuple(k)
  if value is not None:
   if k in inputs:assert inputs[k]==Q(value)
   inputs[k]=Q(value);loc[k]=locator
  return{(k,):Q(1)}
 phys={r['physical_row']:r for r in d['C2']['physical12_coefficients']};raw=d['raw33']['actual33_basis']
 def delta(row,ix):return var(('cache',row,ix),dict(phys[row]['delta']).get(ix,'0'),['C2',row,'delta',ix])
 def bb(row,ix):return var(('bB',row,ix),dict(phys[row]['bB']).get(ix,'0'),['C2',row,'bB',ix])
 lam={r['physical_row']:r['delta_r_ell2']for r in d['lambda']['all3_lambda_readbacks']}
 def dl(row):return var(('lambda12',row),lam[row],['lambda',row,'delta_r_ell2'])
 bblambda=var(('bB63_lambda12',),d['new_lambda']['requested_bB63_r_ell2'],['new_lambda','requested_bB63_r_ell2'])
 missing={k:var(('NEW_delta63_plane4',k))for k in['r4','r3ell','r2ell2']}
 bb03=var(('NEW_CACHED_bB63_ell3',))
 def affine(expr):
  return sum_poly(sc(delta(k[1],k[2]),Q(q))for k,q in expr)
 kap=[[affine(expr)for expr in row]for row in d['q7']['kappa_k_A_B_UNPAIRED']]
 qi={}
 for col in d['q7']['projector7_by660_sparse_columns']:
  e=col['monomial']
  if col['raw_target']==2 and sum(e)==1 and any(e[:7]):qi[e.index(1)]=affine(col['target7_affine'][0][1])
 u=d['unit'];n0=[[Q(0)]*24 for _ in range(56)]
 for i,j in enumerate(u['fixed_selected24_columns']):n0[j]=list(map(Q,u['fixed24_minor_inverse'][i]))
 w=[{}for _ in range(24)];w[7]=missing['r4'];w[23]=sc(bb(63,29),Q(-1))
 for i in range(7):w[i]=delta(63,raw.index([i+1,1,0]));w[16+i]=sc(bb(63,raw.index([i+1,1,0])),Q(-1))
 v=[sum_poly(sc(w[k],-n0[j][k])for k in range(24))for j in range(56)]
 assert all(not v[j]for j in[0,4,8,12,16,20,24])
 H2=[[Q(q)for q in row]for row in d['plane']['actual_plane2_H_three_monomials_rows3_columns56']]
 h11=sum_poly(sc(v[j],H2[1][j])for j in range(56));h02=sum_poly(sc(v[j],H2[2][j])for j in range(56))
 a1,b1=Q(8,25),Q(96,275);L=dict((i,Q(q))for i,q in d['adjoint']['fixed_L_Cn0_row24_sparse'])
 beta=[[Q(q)for q in row]for row in d['q7']['beta2_plane_after_actual_s_equals_r_plus_ell_over5']]
 hh=[sc(bb(63,30),Q(-1)),sc(bblambda,Q(-1)),sc(bb03,Q(-1))]
 Ah=add(add(missing['r2ell2'],sc(missing['r3ell'],-a1)),sum_poly(mul(kap[0][k],p)for k,p in enumerate(hh)))
 Bh=sum_poly(mul(kap[1][k],p)for k,p in enumerate(hh))
 RAv={};RBv={}
 for row,j in[(18,1),(60,17),(67,21)]:RAv=add(RAv,mul(add(dl(row),sc(delta(row,30),-a1)),v[j]))
 RAv=add(RAv,sc(v[28],Q(1,6)-a1*Q(5,3)))
 for row,j in[(18,2),(60,18)]:RBv=add(RBv,mul(add(dl(row),sc(delta(row,30),-b1)),v[j]))
 RBv=add(RBv,mul(add(sc(bblambda,Q(-1)),sc(bb(63,30),b1)),v[28]))
 Ah=add(Ah,add(RAv,add(mul(kap[0][0],h11),mul(kap[0][1],h02))))
 Bh=add(Bh,add(RBv,add(mul(kap[1][0],h11),mul(kap[1][1],h02))))
 for i in[0,4,5]:
  Ai=b1*beta[i][0]-beta[i][1]
  mixe=sum_poly(mul(p,v[j])for p,j in[(sc(bb(63,raw.index([i+1,0,1])),Q(-1)),28),(delta(18,raw.index([i+1,0,1])),2),(delta(60,raw.index([i+1,0,1])),18)])
  mixr=sum_poly(mul(p,v[j])for p,j in[(sc(bb(63,raw.index([i+1,1,0])),Q(-1)),28),(delta(18,raw.index([i+1,1,0])),2),(delta(60,raw.index([i+1,1,0])),18)])
  Bh=add(Bh,add(sc(mixe,Ai),sc(mixr,-L[8+i])))
  Bh=add(Bh,mul(qi[i],sc(bb(63,raw.index([i+1,0,1])),Q(-1))))
 canonical=[sc(Ah,Q(5,3)),add(Bh,sc(Ah,Q(-121,48)))];polys=[Ah,Bh]+canonical
 used={label for p in polys for mon in p for label in mon}
 required=used-set(inputs);assert required=={('NEW_delta63_plane4',k)for k in['r4','r3ell','r2ell2']}|{('NEW_CACHED_bB63_ell3',)}
 # Extract sensitivity OPERATORS, not their unknown responsevalues. Every
 # missingcallbackislinear; there is no productoftwo missingcallbacks.
 sensitivities=[]
 for p in polys:
  rows=[]
  for label in sorted(required):
   op={}
   for mon,q in p.items():
    if label in mon:
     assert mon.count(label)==1 and sum(x in required for x in mon)==1
     tail=tuple(x for x in mon if x!=label);op[tail]=op.get(tail,Q(0))+q
   rows.append([list(label),enc(op)])
  sensitivities.append(rows)
 outdata={
  'scope':'NEXT SOURCE-only coherentoldhead2=free29=physical63*tauA. NO actualfourthfunctional,mixedthird,rootell³compositionorhead29edgevaluecomputed.',
  'actual_sigma0':'e29','literal_full_B_component_column29':'0 exactly in the accepted originalQ48 frame',
  'sigma1_r_w24':[[i,enc(p)]for i,p in enumerate(w)if p],
  'sigma1_r_v56':[[i,enc(p)]for i,p in enumerate(v)if p],
  'sigma1_r_source_proof':'Asectorcurvequotientneedsdelta63_r4 (retainedNEW); Bsectoridentically0; Hsectorunshiftedr²quotientneeds−bB63_r³. CoordinateHadamardfirstsarecachedδ63zi*r and−bB63zi*r. NoprefixBshift/mixedsource3isneededforTHIShead.',
  'sigma1_z_and_m2_cancellation':'Ahighleftsigma1_z coefficients0. BhighleftRz*n0 hasonlyBzrows8/12/13, alreadyverifiedindependently; allthosem1_z*e29=0sinceentireBcomponent0. LikewiseL*m2*e29=0. ThisdoesNOTdrop−L*m1_r*sigma1_r,whichisexplicitlyretainedinbothmixedcarrypolynomialterms.',
  'native_fourth_and_Hthird':'A4sigma0=delta63plane4, B4sigma0=0,H3sigma0=−bB63plane3. AllgenuineκH2/qHzi andmixedB V2sigma1/cyclicprefixtermsareincluded.',
  'known_cached_inputs':[[list(k),str(q),loc[k]]for k,q in sorted(inputs.items())],
  'four_genuine_and_canonical_high_polynomials_UNPAIRED':[enc(p)for p in polys],
  'retained_missing_callbacks':[list(k)for k in sorted(required)],
  'linear_missing_response_sensitivity_operators_UNPAIRED':sensitivities,
  'new_cache_composition_scope':'bB63ell³ is existingprovidedjoint105degree3/root38648 underpureell4; B6head23jobannihilatesell³soDO_NOTinferitfromthe2039/360r ell²readback.',
  'selected_IFT4_theoretical_gate':'For g=delta63(root20144−40922), sourceF105(Φ,Y)=0 has SAMEfixedunitJ105. Forα∈{r4,r³ell,r²ell²}, unknownΦ4 enterslinearly: g4_α=g_lower4_α−(gradΦg|0 J105^-1)*RHS4_α. RHS4andg_lower4useONLYalreadyprovidedFULLplaneΦ≤3 andlinearYA=(z7,r,r+ell/5). Onefixed105rowρservesall3monomials. Thisisaderivation/formulaONLY: actualρ/sourceclosure/positiveheight/owned/RAWpreparedcodepeerandROOTselectedfunctionalONEGOarenotyetclosed.',
  'no_higher_value_or_worker_GO_from_this_packet':True,
  'all_other_uncomputed_head_highs_retained':True,
  'no_originalKbar_type_H0_nilcube_claim':True,
  'selective_source_receipts':receipts,
 }
 owned=S.count(d)+S.count(outdata)+S.count(v)+S.count(w)+20000;assert owned<100000;outdata['whole_source_owned_numeric_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'coherent-head29-high-operator-and-three-selected-plane4-source-gate.json';ep.write_text(json.dumps(outdata,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_HEAD29_COHERENT_THREE_PLANE4_AND_ONE_CACHE03_CALLBACKS_RETAINED',
  'new_mixedthird_callbacks_needed_for_head29':False,'three_plane_fourth_callbacks_retained':True,
  'actual_head29_or_newfourth_values_computed':False,'selectedIFT4_code_budget_or_GO_closed':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'polynomial_helper_sha256':sha(ROOT/'code/ambient_B18_C_connecting_head23_source_20261105_v63.py'),
  'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{ep.name:sha(ep)},'whole_numeric_upper':owned,'seconds':time.monotonic()-st}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()

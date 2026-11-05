#!/usr/bin/env python3
"""SOURCE only: full coherent head30, with safe B63 callbacks unpaired."""
import hashlib,json,math,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
from ambient_B18_C_connecting_head23_source_20261105_v63 import add,sc,mul,const,enc,sum_poly
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_CONNECTING_HEAD30_SOURCE_20261105_v63'
PATHS={
 'q7':'reports/AMBIENT_B18_C_GENUINE_PROJECTION7_SOURCE_20261105_v63/sparse-actual-Dlinear-projector7-source-and-five-head-order-frontier.json',
 'q7_peer':'reports/AMBIENT_B18_C_GENUINE_PROJECTION7_SOURCE_20261105_v63_REVIEW/report.json',
 'unit':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
 'unit_peer':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63_REVIEW/report.json',
 'head23':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_SOURCE_20261105_v63/full-coherent-head23-two-socle-cache-polynomials-and-low-two-row-source.json',
 'head23_peer':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_SOURCE_20261105_v63_REVIEW/report.json',
 'plane':'reports/AMBIENT_B18_C_CONNECTING_PLANE_ORDER_SOURCE_20261105_v63/all56-plane-orders-surviving-sigma2-adjoint-and-U9-namespace-correction.json',
 'plane_peer':'reports/AMBIENT_B18_C_CONNECTING_PLANE_ORDER_SOURCE_20261105_v63_REVIEW/report.json',
 'adjoint':'reports/AMBIENT_B18_C_CONNECTING_SIGMA2_ADJOINT_SOURCE_20261105_v63/actual-three-row-sigma2-source-adjoint-and-coherent-fourth-frontier.json',
 'adjoint_peer':'reports/AMBIENT_B18_C_CONNECTING_SIGMA2_ADJOINT_SOURCE_20261105_v63_REVIEW/report.json',
 'five':'reports/AMBIENT_B18_C_CONNECTING_FIVE_FORMAL_SOURCE_20261105_v63/actual-five-head-formal-RegDivell-operator-and-source-binding.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'lambda':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63/phase1-six-root-existing-cache-readbacks.json',
 'lambda_peer':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
 'beta':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'beta_peer':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63_REVIEW/report.json',
 'alpha_peer':'reports/AMBIENT_B18_C3_ALPHA4_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'graph_peer':'reports/AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63_REVIEW/report.json',
 'B63':'reports/AMBIENT_B18_C_CONNECTING_B63_SAFE_OUTPUT_SOURCE_20261105_v63/finite63-socle-safe-output-theorem-and-actual-lower-input-overlaps.json',
 'B63_peer':'reports/AMBIENT_B18_C_CONNECTING_B63_SAFE_OUTPUT_SOURCE_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={'q7':['kappa_k_A_B_UNPAIRED','projector7_by660_sparse_columns','beta2_plane_after_actual_s_equals_r_plus_ell_over5','beta3_plane_after_actual_s_equals_r_plus_ell_over5'],
  'unit':['fixed_selected24_columns','fixed24_minor_inverse'],
  'head23':['sigma1_z_complete_leftadjoints_Rz_n0'],
  'plane':['actual_plane2_B_three_monomials_rows3_columns56','actual_plane2_H_three_monomials_rows3_columns56'],
  'adjoint':['actual_sigma0_56_by5_sparse','fixed_L_Cn0_row24_sparse','actual_CB_shift_c7'],
  'five':['full_V56_reference_pairs_columns56_targets3'],
  'C2':['physical12_coefficients'],'lambda':['all3_lambda_readbacks'],
  'raw33':['actual33_basis'],'beta':['gA_plane3'],
  'B63':['basis63_monomials','socle32_indices']}
 d={}
 for k,p in PATHS.items():
  if k in fields:d[k]=S.selective(ROOT/p,fields[k],receipts)
  elif k.endswith('_peer'):
   d[k]=json.loads((ROOT/p).read_text());assert d[k]['all_passed']
 assert d['adjoint']['actual_sigma0_56_by5_sparse'][3]==[[30,'1']]
 assert d['five']['full_V56_reference_pairs_columns56_targets3'][30][0]==[0,1]
 inputs={};loc={}
 def var(k,value=None,locator=None):
  k=tuple(k)
  if value is not None:
   if k in inputs:assert inputs[k]==Q(value)
   inputs[k]=Q(value);loc[k]=locator
  return{(k,):Q(1)}
 phys={r['physical_row']:r for r in d['C2']['physical12_coefficients']};raw=d['raw33']['actual33_basis']
 def cached(side,row,ix):return var((side,row,ix),dict(phys[row][side]).get(ix,'0'),['C2',row,side,ix])
 def delta(row,ix):return cached('delta',row,ix)
 def odd(side,row,ix):return cached(side,row,ix)
 safe_basis={tuple(m)for m in d['B63']['basis63_monomials']}
 def safe(m):
  m=tuple(sorted(m));assert m in safe_basis
  return var(('SAFE_delta63',)+m)
 def olambda(side,row):return var(('TRUSTED_plane3',side,row,1,2))
 def dl(row):
  v={r['physical_row']:r['delta_r_ell2']for r in d['lambda']['all3_lambda_readbacks']}
  if row in v:return var(('delta_lambda12',row),v[row],['lambda',row,'delta_r_ell2'])
  assert row in[19,28,29,72]
  return var(('delta_lambda12',row),'0',['alpha_peer_and_graph_peer',row,'actual_delta_r_ell2_zero'])
 ga=[Q(0)]*4
 for (nr,ns),q in d['beta']['gA_plane3']:
  for j in range(ns+1):ga[j]+=Q(q)*Q(math.comb(ns,j),5**j)
 assert ga==list(map(Q,['125/24','5/3','1/6','2/375']))
 gA=[var(('gA_plane3',i),q,['beta','gA_plane3','s=r+ell/5',i])for i,q in enumerate(ga)]
 def affine(expr):return sum_poly(sc(delta(k[1],k[2]),Q(q))for k,q in expr)
 kap=[[affine(expr)for expr in row]for row in d['q7']['kappa_k_A_B_UNPAIRED']]
 qi={}
 for col in d['q7']['projector7_by660_sparse_columns']:
  e=col['monomial']
  if col['raw_target']==2 and sum(e)==1 and any(e[:7]):qi[e.index(1)]=affine(col['target7_affine'][0][1])
 u=d['unit'];n0=[[Q(0)]*24 for _ in range(56)]
 for i,j in enumerate(u['fixed_selected24_columns']):n0[j]=list(map(Q,u['fixed24_minor_inverse'][i]))
 shifts=list(map(Q,d['adjoint']['actual_CB_shift_c7']));active=[i for i,q in enumerate(shifts)if q];assert active==[0,4,5]
 # Adapted y=z+c*r² curve is used only for B's r³ quotient.
 shifted4=add(safe((7,7,7,7)),sum_poly(sc(safe((i,7,7)),-shifts[i])for i in active))
 shifted4=add(shifted4,sum_poly(sc(safe((i,j)),shifts[i]*shifts[j])for i in active for j in active if j>=i))
 w=[{}for _ in range(24)];w[15]=shifted4;w[23]=gA[0]
 for i in range(7):
  w[8+i]=delta(63,raw.index([i+1,1,0]));w[16+i]=odd('bA',63,raw.index([i+1,1,0]))
 v=[sum_poly(sc(w[k],-n0[j][k])for k in range(24))for j in range(56)]
 H2=[[Q(q)for q in row]for row in d['plane']['actual_plane2_H_three_monomials_rows3_columns56']]
 h11=sum_poly(sc(v[j],H2[1][j])for j in range(56));h02=sum_poly(sc(v[j],H2[2][j])for j in range(56))
 a1,b1=Q(8,25),Q(96,275);L={i:Q(q)for i,q in d['adjoint']['fixed_L_Cn0_row24_sparse']}
 beta=[list(map(Q,row))for row in d['q7']['beta2_plane_after_actual_s_equals_r_plus_ell_over5']]
 beta3=[list(map(Q,row))for row in d['q7']['beta3_plane_after_actual_s_equals_r_plus_ell_over5']]
 Ah=sum_poly(mul(kap[0][k],gA[k+1])for k in range(3))
 Bh=add(add(safe((7,7,8,8)),sc(safe((7,7,7,8)),-b1)),sum_poly(mul(kap[1][k],gA[k+1])for k in range(3)))
 # Full native H mixed and B mixed are preserved before any sigma terms.
 for i in active:
  Ai=b1*beta[i][0]-beta[i][1];Bi=-beta[i][0]
  Bh=add(Bh,mul(qi[i],odd('bA',63,raw.index([i+1,0,1]))))
  Bh=add(Bh,add(sc(safe((i,7,8)),Ai),sc(safe((i,8,8)),Bi)))
  gamma=b1*beta3[i][0]-beta3[i][1]
  Bh=add(Bh,sc(delta(63,raw.index([i+1,0,1])),gamma))
 rows=[18,19,28,29,60,67,72,63]
 rA={};rB={};mixr={i:{}for i in active};mixe={i:{}for i in active}
 for p,row in enumerate(rows):
  j=4*p
  # Mask0 has (bA,-bB,0), masks1/2 have delta in A/B.
  ar=add(gA[2],sc(gA[1],-a1))if row==63 else add(olambda('bA',row),sc(odd('bA',row,30),-a1))
  br=add(olambda('bB',row),sc(odd('bB',row,30),-b1))
  rA=add(rA,mul(ar,v[j]));rB=add(rB,mul(sc(br,Q(-1)),v[j]))
  if p<7:
   rA=add(rA,mul(add(dl(row),sc(delta(row,30),-a1)),v[j+1]))
   rB=add(rB,mul(add(dl(row),sc(delta(row,30),-b1)),v[j+2]))
  else:assert not v[j+1]and not v[j+2]and not v[j+3]
  for i in active:
   ir,ie=raw.index([i+1,1,0]),raw.index([i+1,0,1])
   mixr[i]=add(mixr[i],mul(sc(odd('bB',row,ir),Q(-1)),v[j]))
   mixe[i]=add(mixe[i],mul(sc(odd('bB',row,ie),Q(-1)),v[j]))
   if p<7:
    mixr[i]=add(mixr[i],mul(delta(row,ir),v[j+2]));mixe[i]=add(mixe[i],mul(delta(row,ie),v[j+2]))
 # The first seven old C0 odd entries use the same bA plane functions.
 for i,row in enumerate(rows[:7]):
  ar=add(olambda('bA',row),sc(odd('bA',row,30),-a1))
  rA=add(rA,mul(sc(ar,Q(-1)),v[48+i]))
 assert not v[55]
 Ah=add(Ah,add(rA,add(mul(kap[0][0],h11),mul(kap[0][1],h02))))
 Bh=add(Bh,add(rB,add(mul(kap[1][0],h11),mul(kap[1][1],h02))))
 rzmap={i:{k:Q(q)for k,q in r}for i,r in d['head23']['sigma1_z_complete_leftadjoints_Rz_n0']}
 sigma1z_terms=[];m2_terms=[]
 for i in active:
  Ai=b1*beta[i][0]-beta[i][1]
  Bh=add(Bh,add(sc(mixe[i],Ai),sc(mixr[i],-L[8+i])))
  # Ordered Hadamard: m1_zi rows8+k, k<=i, read delta_(zk*zi).
  zz=sum_poly(sc(safe((k,i)),-q)for k,q in[(k-8,q)for k,q in rzmap[i].items()]if k<=i)
  Bh=add(Bh,zz);sigma1z_terms.append([i,enc(zz)])
  # ORIGINAL-z second, including its adapted-coordinate prefix correction.
  hi=add(safe((i,7,7)),sc(safe((i,i)),-shifts[i]))
  hi=add(hi,sum_poly(sc(safe((i,j)),-shifts[j])for j in active if j>i))
  Bh=add(Bh,sc(hi,-L[8+i]));m2_terms.append([8+i,enc(hi)])
 polys=[Ah,Bh,sc(Ah,Q(5,3)),add(Bh,sc(Ah,Q(-121,48)))]
 labels={k for p in polys for mon in p for k in mon};required=labels-set(inputs)
 assert all(k[0]in['SAFE_delta63','TRUSTED_plane3']for k in required)
 data={
  'scope':'SOURCE full coherent head3/free30/physical63*tauB, not values. Entire delta callbacks are safe B63 outputs; other callback outputs are restricted to trusted provided plane3.',
  'actual_cycle_origin':'e30; native FULL V=(0,delta63,bA63), low two0. Fixed old section gives sigma1=-n0*m1*e30, with its complete n1 recurrence retained.',
  'sigma1_r_source24':[[i,enc(p)]for i,p in enumerate(w)if p],
  'sigma1_r56':[[i,enc(p)]for i,p in enumerate(v)if p],
  'adapted_curve_fourth_proof':'B r³ quotient only: [r4]delta(z=-c r²,r)=delta_r4-sum c_i delta_zi*r²+sum_(i<=j)c_i c_j delta_zi*zj. H r² quotient remains ORIGINAL z=0 and reads bA_r³. All shift-active coefficients are kept, including diagonal zz without an extra factor2.',
  'sigma1_z_projected_terms':sigma1z_terms,
  'sigma1_z_proof':'High-A left adjoint0. High-B Rz*n0 supported Bz rows8/12/13. For head30 B=delta63, so original ordered-prefix first in zi is delta_zk*zi in Hadamard row8+k iff k<=i. H bA_zz terms are killed by the accepted leftadjoint; no unknown zz input is zero-filled.',
  'm2_original_z_source_terms':m2_terms,
  'all_native_H_and_B_mixed_carries_retained':True,
  'full_sigma1_r_primitive_B_mixed_and_minus_L_m1sigma1_retained':True,
  'all_sigma2_contribution':'-L*m2*e30-L*m1_r*sigma1r, exactly the accepted recurrence; all3 original-z Hadamard rows retained.',
  'Fj_mask0_and_old_C0_scope':'Use actual mask0(bA,-bB,0), mask1(delta,0,-bB), mask2(0,delta,bA). Plane3 mask3 is0 by accepted nu/b odd low orders. Never infer first7 bA plane3zero: their trusted root lambda callbacks are included explicitly and enter old -U entries with the same sign/basis. v55=0 is literal selected-section support.',
  'known_cached_inputs':[[list(k),str(q),loc[k]]for k,q in sorted(inputs.items())],
  'four_genuine_and_canonical_high_polynomials_UNPAIRED':[enc(p)for p in polys],
  'required_safe_or_trusted_callbacks':[list(k)for k in sorted(required)],
  'maximum_polynomial_degree':max(len(m)for p in polys for m in p),
  'no_unknown_Phi_coefficients_dropped':'Safe entire delta63 follows B63 socle theorem and independent 105-gradient0. Other functions retain every required plane3 callback; full unknown Phi3mixed/Phi2zz/Phi4 never asserted0.',
  'no_head_values_or_new_source_solve':True,'no_original_H0_or_nil_claim':True,
  'selective_source_receipts':receipts,
 }
 owned=S.count(d)+S.count(inputs)+S.count(v)+S.count(w)+S.count(polys)+S.count(data)+20000;assert owned<150000
 data['whole_source_owned_numeric_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'full-coherent-head30-and-safe63-cached-callbacks-source.json';ep.write_text(json.dumps(data,separators=(',',':'))+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_HEAD30_FULL_COHERENT_B63_SAFE_AND_TRUSTED_CALLBACKS_UNPAIRED',
  'actual_head30_or_new_root_values_evaluated':False,'source_solve_needed':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{ep.name:sha(ep)},'whole_source_numeric_upper':owned,'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()

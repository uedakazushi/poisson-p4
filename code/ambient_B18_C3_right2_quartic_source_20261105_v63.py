#!/usr/bin/env python3
"""Coherent two Hilbert--Burch quartic sources; all coefficient leaves formal."""
import hashlib,json,math,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63'
PATHS={
 'frame':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63/actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json',
 'kernel':'reports/AMBIENT_B18_C3_TAUA_CACHED9_SOURCE_20261105_v63/nine-complete-tauA-source-kernel-and-unpaired-cubic-cache-adjoints.json',
 'cover':'reports/AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63/actual-cover-NF-and-five-triangular-cubic-quartic-contradictions.json',
 'B4':'reports/AMBIENT_B18_C3_TAUA_B4_SOURCE_V2_20261105_v63/four-actual-cubic-sources-two-cached-adjoints-and-rank2-completion.json',
 'B4_peer':'reports/AMBIENT_B18_C3_TAUA_B4_SOURCE_V2_20261105_v63_REVIEW/report.json',
 'alpha':'reports/AMBIENT_B18_C3_ALPHA4_CACHED_VALUES_20261105_v63/eight-cached-even-root-values-and-four-complete-alpha-readbacks.json',
 'alpha_peer':'reports/AMBIENT_B18_C3_ALPHA4_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'complete_peer':'reports/AMBIENT_B18_C3_COMPLETE_CUBIC_LINE_SOURCE_20261105_v63_REVIEW/report.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(a,b,scale=Q(1)):
 c=a.copy()
 for k,v in b.items():
  c[k]=c.get(k,Q(0))+scale*v
  if not c[k]:del c[k]
 return c
def scal(a,q):return{k:v*q for k,v in a.items()if v*q}
def combination(polys,weights):
 r={}
 for p,w in zip(polys,weights):r=add(r,p,w)
 return r
def mm(a,b):return[[sum((x*y for x,y in zip(row,col)),Q(0))for col in zip(*b)]for row in a]
def encpoly(p):return[[list(k),str(v)]for k,v in sorted(p.items())]


def main():
 started=time.monotonic();signal.alarm(60);out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 pins={p:sha(ROOT/p)for p in PATHS.values()};d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 assert d['alpha']['complete_cubic_extra_rank']==0 and d['alpha']['actual_C3_HF_degree3']==4
 records=d['frame']['all12_verified_source_records'];piv=[0,4,5];rows=[18,60,67]
 assert[records[i]['physical_row']for i in piv]==rows
 B=[[Q(d['kernel']['bB_plane2_rows3_columns12'][m][i])for i in piv]for m in range(3)]
 T=[[-Q(v)for v in row]for row in d['kernel']['fixed_inverse']]
 minusB=[[-v for v in row]for row in B];eye=[[Q(i==j)for j in range(3)]for i in range(3)]
 assert mm(minusB,T)==eye==mm(T,minusB)
 H=[[Q(records[p]['bB_first_z7'][i])for i in piv]for p in piv]
 for p in piv:assert all(Q(records[p]['bB_first_z7'][i])==0 for i in[1,2,3])
 HT=mm([list(col)for col in zip(*H)],T)
 beta=[[Q(v)for v in row]for row in d['B4']['beta2_plane_r2_rell_ell2']]
 active=d['B4']['active_beta_indices'];assert active==[0,4,5]
 basis=d['raw33']['actual33_basis'];ir3=basis.index([0,3,0]);ir2l=basis.index([0,2,1])
 izr=[basis.index([i+1,1,0])for i in range(7)];izl=[basis.index([i+1,0,1])for i in range(7)]
 c2={r['physical_row']:r for r in d['C2']['physical12_coefficients']}
 # Drop row72 only by the independently proved full set of needed
 # low/pure cubic values, never by a tangent or field coincidence.
 assert dict((r['physical_row'],r['delta_r_ell2'])for r in d['alpha']['all4_actual_alpha_delta_B6'])[72]=='0'
 v72=dict(c2[72]['delta'])
 assert all(Q(v72.get(i,'0'))==0 for i in[ir3,ir2l]+[i for z in active for i in[izr[z],izl[z]]])
 curves=[list(map(Q,r['curve3_r3_r2ell_rell2_ell3']))for r in d['cover']['two_block_source_curve_NF_receipts']]
 a1,a2=curves[0][1]/curves[0][0],curves[0][2]/curves[0][0]
 b1,b2=curves[1][1]/curves[1][0],curves[1][2]/curves[1][0]
 # D_i3 and E_i3=D_i3-gamma_i3 are linear forms in24 old C2 leaves,
 # three prospective lambda readbacks and three unused ell^3 leaves.
 D=[];E=[]
 for physical in rows:
  ds=[{('cache',physical,ir3):Q(1)},{('cache',physical,ir2l):Q(1)},
      {('lambda',physical,-1):Q(1)},{('opaque_ell3',physical,-1):Q(1)}]
  gamma=[{}for _ in range(4)]
  for z in active:
   for ld in range(3):
    gamma[ld]=add(gamma[ld],{('cache',physical,izr[z]):beta[z][ld]})
    gamma[ld+1]=add(gamma[ld+1],{('cache',physical,izl[z]):beta[z][ld]})
  D.append(ds);E.append([add(x,y,Q(-1))for x,y in zip(ds,gamma)])
 def L(poly,c1):return add(poly[1],poly[0],-c1)
 def R(poly,c1,c2):return combination([poly[2],poly[1],poly[0]],[Q(1),-c1,c1*c1-c2])
 LA=[L(p,a1)for p in D];RA=[R(p,a1,a2)for p in D]
 LB=[L(p,b1)for p in E];RB=[R(p,b1,b2)for p in E]
 operator=[[{}for j in range(2)]for i in range(2)]
 for j in range(2):
  operator[0][j]=add(combination(LA,[T[i][j]for i in range(3)]),combination(RA,[T[i][j+1]for i in range(3)]),Q(-1))
  operator[1][j]=add(combination(LB,[HT[i][j]for i in range(3)]),combination(RB,[HT[i][j+1]for i in range(3)]),Q(-1))
 assert all(not any(k[0]=='opaque_ell3'for k in p)for r in operator for p in r)
 sensitivity=[[[operator[t][j].get(('lambda',p,-1),Q(0))for p in rows]for j in range(2)]for t in range(2)]
 dual=[add(operator[1][j],operator[0][j],Q(1,4))for j in range(2)]
 cached_combo=add(dual[1],dual[0],Q(10))
 assert not any(k[0]=='lambda'for k in cached_combo)
 data={
  'actual_scope':'Two known quartic initial relations of original C3 after the exact cubic extra rank0 result. A positive2x2 rank2 decides C3length19; this is not a complete quartic presentation ifrank<2.',
  'pivot_physical_rows':rows,'pivot_original_roots':[c2[p]['roots']for p in rows],
  'q_p_actual_raw116_columns':[[[80+3*p,'1']]+[[81+3*i,str(Q(records[p]['bB_first_z7'][i]))]for i in range(7)if Q(records[p]['bB_first_z7'][i])]for p in piv],
  'q_p_exact_raw_sector_vectors':['delta_p','sum_i h_pi*delta_i','-bB_p+sum_i h_pi*bA_i'],
  'source_first7_h_rows':[[str(Q(v))for v in records[p]['bB_first_z7']]for p in piv],
  'plane2_matrix_rows_r2_rell_ell2_columns_qp':[[str(-v)for v in row]for row in B],
  'constant_source_normalization_T':[[str(v)for v in row]for row in T],
  'both_inverse_products_identity':True,
  'P_k_actual_allorder_recipe':'q_p haszeroM1. Apply canonical ordered z-division to its FULL actual Taylor2 part and subtract linearfirst-framecolumns toleaveonly -bB_p2*eH. ThenmultiplysourcecolumnsbyQQunitT. ThusPk has FULL second term (r2,r*ell,ell2)eH; unknownzzcoefficientsareactual inputfunctionals, notassignedzero.',
  'source_normalized_plane3_A':'sum_p T_pk*delta_p3',
  'source_normalized_plane3_B':'sum_i (H^t*T)_ik*(delta_i3-sum_z(delta_i2_zi*r*r+delta_i2_zi*ell*ell)*beta_z2)',
  'negative_RB_beta_carry_retained':True,
  'other_h_row72_removed_exactly':'alpha72=0modell3 plusoldr3/r2ell0 andallbeta-active mixedleaves0, soitsBgraph-evaluatedcubic0modell3. Itsunknownell3coefficientisannihilatedbythequarticNF, notsetzero.',
  'HB_sources':[{'label':'HB0','formula':'ell*P0-r*P1'},{'label':'HB1','formula':'ell*P1-r*P2'}],
  'HB_full_order3_cancellation':'The FULL Pk second terms are exactly the displayed plane H quadratics. The two HB combinations therefore havezeroorder3 in every rawsector, includingz; theystartat4.',
  'quartic_projection_carry_proof':'Their fourthterms are ell/r times Pk thirdterms. Everyz-containingthird retainsz and vanishesintheknownfirstinitialquotient; degree3 firstframeclearsintroducecarryonlyat5. H planequarticsareclearedbydegree2 timesactualHquadraticpivots,whoseA/Bthirdtermscontributeonlyat5. Pk fourthandhighertermstimeslinearcoeffstartat5. ThusONLYplane3A/B coefficients aboveenterquarticNF.',
  'highest_unavailable_sensitivities_zero_by_degree':{'q_p_second_zz_coefficients':252,'q_p_second_H_mixed_coefficients':42,'q_p_third_H_plane_coefficients':12,'q_p_third_z_containing_coefficients':1449,'q_p_fourth_all_coefficients':4455},
  'no_missing_mixed_source3_delta4_beta3_or_gauge_input':True,
  'curve_ratios':{'a1':str(a1),'a2':str(a2),'b1':str(b1),'b2':str(b2)},
  'matrix2_rows_A4_B4_columns_HB0_HB1_linear_leaf_operators':[[encpoly(p)for p in row]for row in operator],
  'lambda_sensitivity2_by2_by3':[[[str(v)for v in x]for x in row]for row in sensitivity],
  'lambda_callbacks3':[{'physical_row':p,'functional':'coefficient r*ell2 ofsame-R0 delta_p=aA_p-aB_p','original_roots':c2[p]['roots'],'status':'UNREAD_provided_joint105_plane3_callback'}for p in rows],
  'old_C2_leaf_scope':'Onlyindices29/30andmixedzi*r/zi*ell foractivebetaindices0/4/5, onphysical18/60/67. NoindividualuncachedTaylorleafsetzero.',
  'dual_B4_plus_A4_over4_operators':[encpoly(p)for p in dual],
  'HB1_plus10HB0_dual_all_lambda_sensitivities_zero':True,
  'all_cached_dual_HB1_plus10HB0_operator_NOT_PAIRED':encpoly(cached_combo),
  'positive_rank2_completion':'Two independent actual quarticinitialclasses removebothknowncoverquartics. ActualC3HF36640andlength19,matchingacceptedlower19. DVRtypeatcut3is(1,1,2,e>=3); thelongereandgenuineC[ell]remainopen.',
  'rank0_or1_scope':'Knownsource subsetonly; doNOTinferfullquarticrankorlength20.',
  'all_leaf_values_and_rank_uncomputed':True,
  'future_cache_source_scope':'Onlythreeexistingjoint105plane3even-delta readbacks; noIFT/gauge/newphysicalsource3/4 orpoint. Mustfreeze their literalDAGclosure/height/fullbudgetandobtainindependentSOURCEcodepeerbeforevalues.'
 }
 ep=out/'two-coherent-HB-quartic-sources-and-unpaired-affine-operator.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_ONLY_RIGHT2_quartic_operator_no_lambda_or_value_pairing',
   'right2_values_or_rank_evaluated':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),
   'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)


if __name__=='__main__':main()

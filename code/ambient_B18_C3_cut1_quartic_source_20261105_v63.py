#!/usr/bin/env python3
"""The extra ell*P2-ell^3*eH quartic source; cache adjoints left unpaired."""
import hashlib,json,math,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C3_CUT1_QUARTIC_SOURCE_20261105_v63'
PATHS={
 'right_source':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63/two-coherent-HB-quartic-sources-and-unpaired-affine-operator.json',
 'right_source_peer':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63_REVIEW/report.json',
 'right_value':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63/phase2-actual-RIGHT2-quartics-two-factors-and-rank.json',
 'right_value_peer':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'B4':'reports/AMBIENT_B18_C3_TAUA_B4_SOURCE_V2_20261105_v63/four-actual-cubic-sources-two-cached-adjoints-and-rank2-completion.json',
 'cover':'reports/AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63/actual-cover-NF-and-five-triangular-cubic-quartic-contradictions.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'C2_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(a,b,scale=Q(1)):
 c=a.copy()
 for key,v in b.items():
  c[key]=c.get(key,Q(0))+scale*v
  if not c[key]:del c[key]
 return c
def linear(polys,weights):
 p={}
 for q,w in zip(polys,weights):p=add(p,q,w)
 return p
def en(p):return[[list(k),str(v)]for k,v in sorted(p.items())]


def main():
 started=time.monotonic();signal.alarm(60);out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 pins={p:sha(ROOT/p)for p in PATHS.values()};d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 assert d['right_value']['actual_RIGHT2_rank']==1
 v=list(map(Q,d['right_value']['rank_certificate']['pivot_column']));assert v==[Q(-3,5),Q(-121,80)]
 slope=Q(121,48);assert v[1]-slope*v[0]==0
 source=d['right_source'];T=[list(map(Q,row))for row in source['constant_source_normalization_T']]
 H=[[Q(row[i])for i in[0,4,5]]for row in source['source_first7_h_rows']]
 HT=[[sum((H[j][i]*T[j][k]for j in range(3)),Q(0))for k in range(3)]for i in range(3)]
 beta=[list(map(Q,row))for row in d['B4']['beta2_plane_r2_rell_ell2']]
 active=[i for i,b in enumerate(beta)if any(b)];assert active==[0,4,5]
 curves=[list(map(Q,r['curve3_r3_r2ell_rell2_ell3']))for r in d['cover']['two_block_source_curve_NF_receipts']]
 a1,b1=curves[0][1]/curves[0][0],curves[1][1]/curves[1][0]
 basis=d['raw33']['actual33_basis'];ir3=basis.index([0,3,0]);ir2l=basis.index([0,2,1])
 izr=[basis.index([i+1,1,0])for i in range(7)];izl=[basis.index([i+1,0,1])for i in range(7)]
 A={};B={};rows=[18,60,67]
 for i,physical in enumerate(rows):
  A=add(A,{('cache',physical,ir2l):Q(1),('cache',physical,ir3):-a1},T[i][2])
  E={('cache',physical,ir2l):Q(1),('cache',physical,ir3):-b1}
  for z in active:
   E=add(E,{('cache',physical,izr[z]):b1*beta[z][0]-beta[z][1]})
   E=add(E,{('cache',physical,izl[z]):-beta[z][0]})
  B=add(B,E,HT[i][2])
 residual=add(B,A,-slope)
 assert all(key[0]=='cache'for p in[A,B,residual]for key in p)
 # Literal degree4 ell multiplier on each actual cubic NF kills the
 # lambda and opaque ell3 channels for all possible input coefficients.
 for sector in range(2):
  n4=list(map(Q,d['cover']['two_block_source_curve_NF_receipts'][sector]['NF4_to_r2ell2']))
  n3=list(map(Q,d['cover']['two_block_source_curve_NF_receipts'][sector]['NF3_to_r2ell_rell2'][0]))
  assert n4[1:]==n3 and n4[3:]==[Q(0),Q(0)]
 c2={r['physical_row']:dict(r['delta'])for r in d['C2']['physical12_coefficients']}
 bounds=[]
 for name,p in[('A4',A),('B4',B),('residual',residual)]:
  U=0;D=1;product=0
  for(_,physical,index),w in p.items():
   leaf=Q(c2[physical].get(index,'0'))
   U+=((abs(w.numerator)+w.denominator-1)//w.denominator)*((abs(leaf.numerator)+leaf.denominator-1)//leaf.denominator)
   D=math.lcm(D,w.denominator*leaf.denominator)
   product=max(product,max(abs(w.numerator).bit_length(),w.denominator.bit_length())+max(abs(leaf.numerator).bit_length(),leaf.denominator.bit_length()))
  h=max((U*D).bit_length(),D.bit_length());bounds.append({'output':name,'U':str(U),'D':str(D),'height':h,'direct_product_precheck':product,'terms':len(p)})
 h=max(b['height']for b in bounds);post=12*h+3;assert post<=2048
 data={
  'actual_source':'S_cut=ell*P2 - ell^3*eH in the FULL expanded R0-module presentation of C3=C/ell^3C',
  'literal_cut_scope':'ell^3 eH is an exact additional quotient column with no higher tail. It is not replaced byell*P2; theirdifferenceISthenewquarticsource.',
  'same_actual_P2':'P=q*T from the reviewed RIGHT2 operator, physical pivots18/60/67, after full actual quadraticzclears. FULLTaylor2(P2)=ell²eH inallrawsectors, notmerelyitsplaneprojection.',
  'full_cubic_cancellation':'ell*P2 leadingcubicisexactlyell³eH; subtractingliteralcut cancelsitineverysector. Noformalell-divisionorArtinkernelargument.',
  'quartic_full_carry':'Fourthpartisell*P2,3. Everythirdztermretainszandisclearedbydegree3firstframe, whosehighercarrybegins5. Hquarticsclearbydegree2Hquadraticpivots, whoseA/Bcarrybegins5. P2 fourthterms multipliedbyellstart5. Thus quarticA/B NFisell*theFULLplane3P2A/BfromRIGHT2.',
  'negative_RB_beta_carry_retained':'A=ΣT_i2(delta_i21-a1 delta_i30). B=Σ(HT)_i2[(delta_i21-gamma_i21)-b1(delta_i30-gamma_i30)], gamma_i3=Σδ_i2_zmixed*beta2. Fullr³gamma30isretained.',
  'row72_and_other_source_callbacks':'Same exactrow72dropasRIGHT2. Alllambda_r*ell2andopaqueell3cubechannelsannihilatedbyell/NF, notsetzero. Noν,delta4,mixedΦ3,beta3,newλorgaugeinput.',
  'cubic_to_quartic_NF_identity_bothblocks_verified':True,
  'A4_cached_linear_operator':en(A),'B4_cached_linear_operator':en(B),'residual_B_minus_121A_over48_operator':en(residual),
  'known_quartic_vector':list(map(str,v)),'residual_dual_slope':str(slope),
  'constant_source_normalization_T':source['constant_source_normalization_T'],
  'output_majorants_no_pairing':bounds,'max_output_height_bound':h,'whole_fold_minor_inverse_rank_operation_bound':post,
  'all_values_uncomputed':True,
  'positive_completion':'Anyresidualnonzero givesaquarticindependentofknownHBvector, determinant=(-3/5)*residual. ThusC3HF36640,length19. Longerexponentremains>=3unknown; nooriginalH0/nilclaim.',
  'zero_scope':'Onlytheextra literalcutquarticsourceisdependentonknownHBline; nofullquarticrank/length20inference.',
  'future_worker_budget':{'static_owned_numeric':18000,'hard_owned_numeric':30000,'RAW_input_output_bytes':400000,'hard_seconds':60,'hard_memory_bytes':2147483648,'hard_QQ_bits':2048,'new_physical_source_gauge_IFT_or_jet':0}
 }
 ep=out/'literal-cut-quartic-source-and-three-unpaired-cache-operators.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_ONLY_extra_literal_cut_quartic_NOT_PAIRED','input_sha256':pins,
    'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)


if __name__=='__main__':main()

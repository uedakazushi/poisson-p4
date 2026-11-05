#!/usr/bin/env python3
"""Two cached cubic adjoints on each of four actual tauA sources; no dots."""
import hashlib,json,math,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C3_TAUA_B4_SOURCE_V2_20261105_v63'
PATHS={
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
 'pair':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'frame':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63/actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json',
 'frame_peer':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63_REVIEW/report.json',
 'kernel':'reports/AMBIENT_B18_C3_TAUA_CACHED9_SOURCE_20261105_v63/nine-complete-tauA-source-kernel-and-unpaired-cubic-cache-adjoints.json',
 'kernel_peer':'reports/AMBIENT_B18_C3_TAUA_CACHED9_SOURCE_20261105_v63_REVIEW/report.json',
 'cover':'reports/AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63/actual-cover-NF-and-five-triangular-cubic-quartic-contradictions.json',
 'cover_peer':'reports/AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ceilabs(q):return(abs(q.numerator)+q.denominator-1)//q.denominator
def main():
 started=time.monotonic();signal.alarm(60);out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 pins={p:sha(ROOT/p)for p in PATHS.values()};d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 basis=d['raw33']['actual33_basis'];i30=basis.index([0,3,0]);i21=basis.index([0,2,1]);assert[i30,i21]==[29,30]
 izr=[basis.index([i+1,1,0])for i in range(7)];izl=[basis.index([i+1,0,1])for i in range(7)]
 # Exact nonzz Taylor duals: no unknown zz scalar coefficient can enter.
 duals=izr+izl+[i30,i21]
 for m,nf in zip(d['raw33']['basis74'],d['raw33']['all74_normal_forms']):
  nf=dict(nf)
  assert all(Q(nf.get(i,'0'))==Q(m==basis[i])for i in duals)
 beta=[]
 for row in d['pair']['beta2_before_ell']:
  co=[Q(0)]*3
  for mon,s in row:
   assert not any(mon[:7])and sum(mon)==2
   nr,ns=mon[7:]
   for ld in range(ns+1):co[ld]+=Q(s)*math.comb(ns,ld)/5**ld
  beta.append(co)
 active=[i for i,v in enumerate(beta)if any(v)];assert active==[0,4,5]
 a,b=[list(map(Q,r['curve3_r3_r2ell_rell2_ell3']))for r in d['cover']['two_block_source_curve_NF_receipts']]
 a2,b1,b2=a[2]/a[0],b[1]/b[0],b[2]/b[0]
 rows=[19,28,29,72];physical=d['kernel']['physical_row_order'];free=d['kernel']['free_columns9'];kernel=[[Q(x)for x in row]for row in d['kernel']['complete_constant_kernel12_by9']]
 frame={r['physical_row']:r for r in d['frame']['all12_verified_source_records']};c2={r['physical_row']:r for r in d['C2']['physical12_coefficients']}
 ops=[];bounds=[];source_rows=[];generic_nf_checks=[];literal_zero_leaves=[]
 for physical_row in rows:
  i=physical.index(physical_row);ki=free.index(i)
  assert[kernel[j][ki]for j in range(12)]==[Q(j==i)for j in range(12)]
  expected=[Q(0)]*7;expected[i]=Q(-1,4)
  assert list(map(Q,frame[physical_row]['bB_first_z7']))==expected
  assert list(map(Q,frame[physical_row]['bA_first_z7']))==[Q(j==i)for j in range(7)]
  source_rows.append({'physical_row':physical_row,'first7_index':i,'literal_kernel_column':ki,
      'actual_source_column':'F_j*(tauA-tauB/4)',
      'literal_RH_binding':'RH_i is the same full physical F_j*tauB, not just a first-jet replacement',
      'original_physical_coefficient_roots':c2[physical_row]['roots'],
      'full_C4_source_before_R0_unit_elimination':['nuA*bA+nuB*bB/4','delta','-delta/4','-bB-bA/4'],
      'all_bB_plane2_zero':all(Q(v)==0 for v in frame[physical_row]['bB_second_plane'])})
  assert source_rows[-1]['all_bB_plane2_zero']
  psi={i21:Q(-1,4),i30:b1/4};chi={i30:(b2-a2)/4}
  for z in active:
   psi[izr[z]]=(beta[z][1]-b1*beta[z][0])/4
   psi[izl[z]]=beta[z][0]/4
   chi[izr[z]]=(beta[z][2]-b2*beta[z][0])/4
   chi[izl[z]]=beta[z][1]/4
  # Independently derive the adjoints from generic full cubic vectors,
  # including the r^3 component of EVERY beta2 mixed carry.
  # Leaves: D30,D21,D12,D03, then seven (delta_zi*r,delta_zi*ell).
  D=[[Q(j==k)for j in range(18)]for k in range(4)]
  carry=[[Q(0)]*18 for _ in range(4)]
  for z in range(7):
   for ell_degree in range(3):
    carry[ell_degree][4+2*z]+=beta[z][ell_degree]/4
    carry[ell_degree+1][5+2*z]+=beta[z][ell_degree]/4
  B=[[(-D[k][j]/4+carry[k][j])for j in range(18)]for k in range(4)]
  generic_psi=[B[1][j]-b1*B[0][j]for j in range(18)]
  generic_chi=[B[2][j]-b2*B[0][j]+(D[2][j]-a2*D[0][j])/4 for j in range(18)]
  leaves=[i30,i21,None,None]+[idx for z in range(7)for idx in[izr[z],izl[z]]]
  for label,actual,weights in[('psi',generic_psi,psi),('chi',generic_chi,chi)]:
   expected_coeff=[weights.get(idx,Q(0))if idx is not None else Q(0)for idx in leaves]
   assert actual==expected_coeff
   generic_nf_checks.append([physical_row,label,list(map(str,actual))])
  pair=[];values=dict(c2[physical_row]['delta'])
  for label,weights in[('psi_B_r2ell',psi),('chi_B_rell2_plus_A_rell2_over4',chi)]:
   terms=[[physical_row,k,str(v)]for k,v in weights.items()if v]
   pair.append({'label':label,'terms':terms})
   # Integer abs/denominator majorants, never signed summand products/dots.
   den=1;absolute_ceil=0;wheight=lheight=0
   for _,idx,s in terms:
    w=Q(s);leaf=Q(values.get(idx,'0'))
    assert leaf==0
    literal_zero_leaves.append([physical_row,label,idx,s,str(leaf)])
    den=math.lcm(den,w.denominator*leaf.denominator)
    absolute_ceil+=ceilabs(w)*ceilabs(leaf)
    wheight=max(wheight,abs(w.numerator).bit_length(),w.denominator.bit_length())
    lheight=max(lheight,abs(leaf.numerator).bit_length(),leaf.denominator.bit_length())
   h=max((absolute_ceil*den).bit_length(),den.bit_length())
   bounds.append({'physical_row':physical_row,'label':label,'denominator_divisor':str(den),
       'integer_absolute_bound':absolute_ceil,'coefficient_height_bound':h,
       'product_precheck_bound':wheight+lheight,'summand_count':len(terms)})
  ops.append({'physical_row':physical_row,'adjoints2':pair})
 h=max(v['coefficient_height_bound']for v in bounds)
 operation=max(12*h+3,max(v['product_precheck_bound']for v in bounds))
 assert operation<2048
 # Both linear functionals are literal rows on the accepted whole cover.
 projection=[[Q(0),Q(0),Q(1),Q(0)],[Q(0),Q(1,4),Q(0),Q(1)]]
 assert projection[0][2]==1 and projection[1][3]==1
 nf4=d['cover']['whole_NF4_rows2_columns15'];nf3=d['cover']['whole_NF3_rows4_columns12']
 assert list(map(Q,nf4[1][6:10]))==list(map(Q,nf3[2][4:8]))
 data={
  'scope':'Actual original C3 module. Four actual cubic sources and two cached cover functionals each; NOT all cubic coefficients or an algebra splitting.',
  'source_rows4':source_rows,
  'beta2_plane_r2_rell_ell2':[[str(v)for v in row]for row in beta],
  'active_beta_indices':active,'curve_ratios':{'a2':str(a2),'b1':str(b1),'b2':str(b2)},
  'raw33_nonzz_degree2_and_plane_degree3_Taylor_duals_verified':len(duals)*len(d['raw33']['basis74']),
  'complete_order2_clearing':'q_j=F_j tauA-RH_j/4 has first0 and plane2zero. q2_A=delta_j2, q2_B=-delta_j2/4. Clear mixed A using exact RA (no higher B); clear mixed B using RB=z+beta, giving +delta_j2mixed*beta2/4. H mixed clears use r/ell RH with eB plane2=actualdelta_i2_plane0. Every zz clearing carries z and no plane3. No beta or H mixed term is omitted.',
  'full_cubic_A_B_projection':'A=delta_j3; B=-delta_j3/4+sum_i(delta_j2_zi*r*r+delta_j2_zi*ell*ell)*beta_i2/4. H is retained then killed by known plane2 H.',
  'psi_formula':'-delta_(r2ell)/4+b1*delta_(r3)/4+sum_i[delta_(zi*r)*(beta_i1-b1*beta_i0)+delta_(zi*ell)*beta_i0]/4',
  'chi_formula':'(b2-a2)*delta_(r3)/4+sum_i[delta_(zi*r)*(beta_i2-b2*beta_i0)+delta_(zi*ell)*beta_i1]/4',
  'unknown_third_cancellation':'chi=B_rell2+A_rell2/4 cancels the SAME actual delta_j3_rell2 exactly. psi assigns that leaf coefficient0 by the literal NF. Neither sets it zero. The other channels and full cubic vectors remain formal.',
  'all144_generic_full_NF_functional_coefficients_verified':generic_nf_checks,
  'all60_nonzero_weight_cached_leaf_accesses_are_literal_zero':literal_zero_leaves,
  'all_eight_corrected_adjoints_zero_without_dot_worker':True,
  'remaining_formal_cubic_vectors':'alpha_j*(A_rell2-B_rell2/4); alpha_j is unknown, not setzero.',
  'cover_functionals_rows2_columns4':[[str(v)for v in row]for row in projection],
  'four_two_cached_adjoints_NOT_PAIRED':ops,
  'each_dot_integer_height_majorant':bounds,
  'max_entry_height_bound':h,'all_dots_minor_inverse_factor_QQ_operation_bound':operation,
  'rank2_completion_theorem':'Rank2 of the actual 2x4 projection certifies at least2 independent cubic initial classes in the known cover. Its psi row is nonzero, so ell times a corresponding source has B quartic NF=psi!=0. Thus known length21 loses>=2 cubic+>=1 quartic dimensions: C3<=18. Accepted lower18 forces equality, HF36621 and original C torsion exponents(1,1,2,2). No unknown full third vector is needed.',
  'rank1_scope':'A nonzero projected class gives at least one cubic and one nonzero r/ell multiple, since the joint r/ell multiplication map on the two CI(3,3) degree3 blocks is injective. Hence only C3 upper19, without equality or both exponents2.',
  'rank0_scope':'Only these eight functionals vanish. Cubic classes may remain unknown multiples of A_rell2-B_rell2/4; quartic sources remain open.',
  'future_worker_profile':{'cache_leaves_upper':32,'dot_products_upper':60,'actual_dot_terms':sum(len(f['terms'])for r in ops for f in r['adjoints2']),
      'all_six_2x2_minors':True,'static_owned_numeric':22000,'hard_owned_numeric':40000,
      'fresh_RAW_total_bytes_bound':500000,'hard_bits':2048,'hard_seconds':60,
      'hard_memory_bytes':2147483648,'new_physical_gauge_IFT_or_jet_coefficients':0},
  'eight_values_and_rank_not_evaluated':True,
 }
 ep=out/'four-actual-cubic-sources-two-cached-adjoints-and-rank2-completion.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'PURE_SOURCE_eight_corrected_projections_zero_by_literal_leaves',
      'eight_values_evaluated':False,'rank_evaluated':False,'eight_projections_proved_zero_without_dot_job':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)


if __name__=='__main__':main()

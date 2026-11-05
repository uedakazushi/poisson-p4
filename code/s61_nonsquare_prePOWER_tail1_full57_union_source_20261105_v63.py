#!/usr/bin/env python3
"""Combine accepted tail1 nonzero/source gates; compare stored supports only."""
import argparse,hashlib,json,math,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_TAIL1_FULL57_UNION_SOURCE_20261105_v63'
MA,ME=1048576,1024
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 t=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={}
 def load(st,name='report.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def accepted(st):
  r=load(st);assert r['all_passed'];return r
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 one=accepted('S61_NONSQUARE_PREPOWER_TAIL1_CACHED_KAPPA_ENTRY_20261105_v63')
 peer=accepted('S61_NONSQUARE_PREPOWER_TAIL1_AXIS_REVIEW_20261105_v63')
 assert peer['T1_polynomial_nonzero_independently_proved']and peer['only_axis_subset_of_all7595_coefficients_peered']
 accepted('S61_NONSQUARE_PREPOWER_DIRECT_M12_FULL15_SOURCE_20261105_v63_REVIEW')
 accepted('S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63_REVIEW')
 accepted('S61_NONSQUARE_PREPOWER_M12_AXIS_ONE_COEFFICIENT_20261105_v63_REVIEW')
 accepted('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63_REVIEW')
 family=load('S61_NONSQUARE_PREPOWER_TAIL8_COMPLEMENT_SOURCE_20261105_v63','typed-family-and-source-plans.json')
 assert family['cofactor_sign']=='kappa_i=(-1)^i det(A_without_row_i); det([A,v])=-kappa*v'
 plan=family['source_plans'][0];assert plan['target']==1 and plan['common_B_power']==1
 assert plan['source_current15_rows']==[2,4,5,8]
 rows=[0,1,3,6,2,11,12,9,7,13,10,4,5,8,14]
 J11=[22,11,13,25,15,19,14,18,28,30,29]
 targets=J11+[27,31,1,17];assert len(set(targets))==15
 assert (-1)**sum(rows[i]>rows[j]for i in range(15)for j in range(i+1,15))==-1
 dump('literal-D15j-source-and-sign.json',dict(
  source15_order_before_unit_operations=rows,
  source_row_operations=['Q9','G18/a','G39/a','G68','G20','G115/c10','(G120-G99)/c11','G92/c8','G85/c6','a*G129/c12','G99/c9','current4','current5','current8','G134'],
  full_source_precondition='normalized row1 <- normalized row1 +2*normalized row2, on EVERY target column including newtail1 and17; determinant1',
  source_permutation_sign=-1,source_GL15_determinant='-1/(a*c6*c8*c9*c10*c11*c12)',
  target15_order=targets,column_units=dict(H27='1/L',H31='1/(L*gamma)',H3_28_30_29='1/L',tail1='B',otherpure='1'),
  total_target_unit='B/(L^5*gamma)',L='a*F*C^2',gamma='-4*F*B*C^2/625',
  source_constants_reference='M12_AXIS_WORKER_GATE/recipe.json#/constants, allnonzeroQQ',
  old_generic16var_identity='DIRECT_M12_FULL15_SOURCE_REVIEW det(A*V-c*u)=A^2 detA4; replace only the arbitrary fourth column by B*M[:,1]',
  permutation_anchor_first=[2,0,1,3],permutation_anchor_first_sign=1,
  pair_H3_literal_zero=[28,30,29],pair_fixed9_raw_pure_rank1=[22,11,13,25,15,19,14,18,17],
  T1='B*kappa*M[current15 rows2,4,5,8;original32 column1]',
  four_by_four_det='det([H27/L,H31/(L*gamma),anchor22,B*tail1])=-T1',
  exact_domain_identity='A_anchor*D15_1=-T1*D12',
  proof='multiply three pair rows by A and subtract their anchor coefficients times G20; Laplace sign positive; detHcorr=A^2*(-T1); cancel A^2 in the source domain',
  no_A_kappa3_Delta_Saux_N6_N7_P_inverse=True,all7_other_tail_effects_retained_in_global15x22_source=True))
 p=ROOT/'reports/S61_NONSQUARE_PREPOWER_TAIL1_CACHED_KAPPA_ENTRY_20261105_v63/integer-column4-and-result.txt';pins[str(p.relative_to(ROOT))]=sha(p)
 with p.open()as f:
  assert f.readline().strip()=='J1KV001'
  for j in range(5):
   node,D=map(int,f.readline().split());n=int(f.readline());S=set()
   for _ in range(n):
    k,c=map(int,f.readline().split());assert c
    if j==4:S.add(k-3*MA)
  D1=D;S1=S;assert n==7595
 two=load('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63','root239.json')['polynomial']
 S2={m[0]*MA+m[1]*ME+m[2]-4*MA for m,c in two};assert len(S2)==7208
 def decode(k):
  a,r=divmod(k,MA);e,d=divmod(r,ME);return(a,e,d)
 def info(S):
  m=[decode(k)for k in S]
  return dict(box=[[min(v[i]for v in m),max(v[i]for v in m)]for i in range(3)],total_degree=max(sum(v)for v in m))
 i1,i2=info(S1),info(S2)
 assert i1==dict(box=[[0,20],[0,41],[0,18]],total_degree=58)
 assert i2==dict(box=[[0,20],[0,40],[0,18]],total_degree=57)
 outside=sorted(S1-(S2|{k+ME for k in S2}));assert len(outside)==78
 witness=decode(outside[0]);assert witness==(11,0,0)
 dump('stored-support-comparison.json',dict(
  U1='T1/a^3',U2='T2/a^4; T2=existingT',supportU1=i1,supportU2=i2,
  if_polynomial_division='deg_a quotient=0,deg_d quotient=0,deg_e quotient=1, so quotient=alpha*e+beta',
  support_obstruction_count=78,first_support_obstruction=list(witness),
  both_possible_preimages_absent_in_U2=True,ordinary_QQ_polynomial_divisibility_excluded=True,
  NOT_proved=['localized divisibility after named-unit factors','radical coprimality','unit ideal','strict FULLU cover enlargement','empty intersection'],
  no_new_coefficients_products_gcd_factor_or_parameter_evaluation=True))
 old=load('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_PREFLIGHT_20261104_v63','operator-and-cost.json')
 D2=int(old['nodes'][old['output']]['QQ_denominator']);assert D2.bit_length()==197
 assert all(D2%Q(c).denominator==0 for m,c in two)
 # Preflight only: C=100(a-a0), a0=(-18e^2+9e+2)/100.
 # Bound all Hasse orders0..20 after ONLY the already visible a shifts.
 maxa=20;maxe=41+2*maxa;maxd=18
 accumulator=(maxe+1)*(maxd+1);powers=(maxa+1)**2
 owned=2*(7595+powers+accumulator)+15000+1000
 coeffbits=363+20+(100**20).bit_length()+2
 products=(7595+7208)*(maxa+1)**2
 RAW=4_000_000+2*(maxa+1)*accumulator*(math.ceil(coeffbits*math.log10(2))+80)
 assert owned<150000 and RAW<32*1024*1024 and coeffbits<2048
 dump('NEXT-known-C-valuation-PREFLIGHT-ONLY.json',dict(
  method='fixed knownFULLU unit C Taylor/Hasse evaluation, no polynomial gcd/factorization or 3var quotient',
  a0='(-18e^2+9e+2)/100',inputs=['T1/a^3','T2/a^4'],
  literal_operator='D*100^20*[z^k]U(a0+z,e,d)=sum(c_integer*binom(i,k)*100^(20-i+k)*(-18e^2+9e+2)^(i-k)*e^j*d^l)',
  orders='0..20, stop at FIRST nonzero coefficient for each input; exact lowest order=v_C',
  source_coefficient_l1_bits_upper363=363,root_numerator_l1=29,root_denominator100=100,
  height_proof='29<100 and binom(i,k)<=2^20, so all integer partial sums <=L1(input)*2^20*100^20',
  integer_coefficient_bits_upper=coeffbits,one_accumulator_support_upper=accumulator,
  cached_univariate_root_power_coefficients_upper=powers,total_mul_pairs_upper=products,
  owned_coefficient_key_metadata_upper=owned,all_output_RAW_upper=RAW,
  no_3var_quotient_or_actual_unit_valuation_evaluated=True,new_parent_budget_REQUIRED=True,
  purpose='certify exact known-unit valuations to shrink T1/T2 before any subsequent multivariate gcd/elimination; by itself removes no base locus'))
 note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,
  exact_j1_AD15_source_units_and_signs_pinned=True,
  full57_sufficient_union='FULLU D(D12) intersect (D(T) union D(T1))',
  T1_nonzero_peer_scope='ONEaxis subset and four column axes, not all7595coefficientidentity',
  no_strict_cover_enlargement_or_bad_ideal_equality=True,
  ordinary_polynomial_nondivisibility_SOURCE_only=True,
  next_C_unit_valuation_ONLY_preflight=True,no_new_parameter_coefficient_job=True,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-t)
 dump('report.json',report);print(json.dumps({k:v for k,v in report.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

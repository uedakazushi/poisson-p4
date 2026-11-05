#!/usr/bin/env python3
"""SOURCE only: common-D two-tail operators and ONE named-J case cost, no new coefficient values."""
from pathlib import Path
import argparse,hashlib,json,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_TWO_TAIL_CANONICAL_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):return json.loads(pin(ROOT/'reports'/st/name).read_text())
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 p0=load('S61_NONSQUARE_PREPOWER_T0_C4_ROOT_LONGDIVISION_REVIEW_20261105_v63');assert p0['all_passed']
 p1=load('S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63');assert p1['all_passed']
 f1=load('S61_NONSQUARE_PREPOWER_Q1_F_ROOT_HORNER_REVIEW_20261105_v63');assert f1['all_passed']
 c0=load('S61_NONSQUARE_PREPOWER_T0_C_ROOT_JET_HORNER_REVIEW_20261105_v63');assert c0['all_passed']
 c1=load('S61_NONSQUARE_PREPOWER_T1_C_HASSE_HORNER_REVIEW_20261105_v63');assert c1['all_passed']
 union=load('S61_NONSQUARE_PREPOWER_TAIL1_FULL57_UNION_SOURCE_20261105_v63');assert union['all_passed']
 literal=load('S61_NONSQUARE_PREPOWER_TAIL1_FULL57_UNION_SOURCE_20261105_v63','literal-D15j-source-and-sign.json')
 assert literal['exact_domain_identity']=='A_anchor*D15_1=-T1*D12'
 globalpeer=load('S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63_REVIEW');assert globalpeer['all_passed']
 records=[]
 for j,st,name,magic,n in [(0,'S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_ENTRY_20261105_v63','quotient-integer-common-D.txt','T0C4QV001',4500),(1,'S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63','quotient-original-common-D.txt','T1C4QN001',4811)]:
  p=pin(ROOT/'reports'/st/name);mins=[10**9]*3;maxs=[0]*3;total=0;l1=0;mb=0;previous=None
  with p.open()as f:
   assert f.readline().strip()==magic;D=int(f.readline());assert D>0 and D.bit_length()==197;assert int(f.readline())==n
   count=0
   for line in f:
    key,c=map(int,line.split());assert c and(previous is None or previous<key);previous=key
    a,x=divmod(key,1048576);e,d=divmod(x,1024)
    mins=[min(u,v)for u,v in zip(mins,(a,e,d))];maxs=[max(u,v)for u,v in zip(maxs,(a,e,d))]
    total=max(total,a+e+d);l1+=abs(c);mb=max(mb,abs(c).bit_length());count+=1
   assert count==n
  records.append(dict(j=j,terms=n,common_D197=D,min_exponents=mins,max_exponents=maxs,total_degree=total,max_integer_bits=mb,L1_bits=l1.bit_length(),literal_sha256=sha(p)))
 assert records[0]['common_D197']==records[1]['common_D197']
 assert [r['min_exponents']for r in records]==[[0,0,0],[0,0,0]]
 assert [r['max_exponents']for r in records]==[[16,32,18],[16,33,18]]
 # Tiny named-unit identity only; no products of J0/J1 or new source values.
 B={(1,0,1):50,(0,2,1):18,(0,1,1):-9,(0,0,1):-2,(0,1,0):60,(0,0,0):-40}
 J={(1,0,1):5,(0,1,0):-6,(0,0,0):-4}
 C={(1,0,0):100,(0,2,0):18,(0,1,0):-9,(0,0,0):-2}
 check=dict(B)
 for k,c in J.items():check[k]=check.get(k,0)+10*c
 for(a,e,d),c in C.items():k=(a,e,d+1);check[k]=check.get(k,0)-c
 check[(0,0,0)]=check.get((0,0,0),0)+80
 assert not any(check.values())
 canonical=dict(SOURCE_ONLY_NO_NEW_COEFFICIENT_VALUES=True,
  ring='R=QQ[a,e,d,1/(a*d*F*B*C*J)]',common_D197=records[0]['common_D197'],
  independently_verified_polynomials=records,
  exact_quotient_identities=['T0=a4*C4*J0/D','T1=a3*C4*J1/D'],
  localized_scheme_ideal_equality='(T0,T1)R=(J0,J1)R; only named FULLU units and nonzero QQ D are inverted',
  sufficient_cover='D(D12) intersect (D(J0) union D(J1)) inside FULLU',
  retained_unknown_complement='V(D12) union V(J0,J1) in FULLU; merely complement of this sufficient cover, NOT actual bad Fitting equality',
  ordinary_monomial_valuations_a_d_zero='both canonical numerators have min_a=min_d=0 by literal nonzero support; no new evaluation',
  named_common_factors_excluded=['a,d by savedsupport','C by exactvC(T0)=vC(T1)=4 and knowna coprimeC','F by independentvF(J1)=0'],
  remaining_named_unit_factors_NOT_tested=['B','J'],
  named_unit_identity='B+10*J-C*d=-80',no_new_unknown_base_chart=True)
 dump('common-D-and-localized-source-identities.json',canonical)
 operators=dict(accepted_pair_rows_current15=[2,4,5,8],same_cached_kappa=True,
  old_exact_source=['T0=B*kappa*Mpair[:,2]','T1=B*kappa*Mpair[:,1]'],
  normalized_covector='Phi=D*B/(a^3*C^4)*kappa',
  normalized_source_columns=['v1=Mpair[:,1]','v0=a^-1*Mpair[:,2]'],
  exact_contractions=['Phi*v1=J1','Phi*v0=J0'],
  source_column_change='on a two-tail source enlargement, diag(1,a^-1) is invertible over FULLU; no kappa3,anchor,Delta,Saux,N6,N7 orP inverse',
  linear_combinations='Phi*(alpha*v1+beta*v0)=alpha*J1+beta*J0 for every parameter-only alpha,beta; no selected new combination value is computed',
  common_D15_source_order=literal['source15_order_before_unit_operations'],
  common_source_GL15=literal['source_GL15_determinant'],
  common_target_prefix13=literal['target15_order'][:13],last_common_target17=17,
  selected_D15='replace ONLY targettail by B*(alpha*M[:,1]+beta*a^-1*M[:,2]); rowprefix/other14 columns unchanged',
  normalized_selected_minors=['E1=D/(a^3*C^4)*D15_1','E0=D/(a^4*C^4)*D15_0'],
  domain_identities=['A*E1=-J1*D12','A*E0=-J0*D12','A*(alpha*E1+beta*E0)=-(alpha*J1+beta*J0)*D12'],
  no_A_inverse=True,all_full57_330source_and_other7tails_retained=True,
  NOT_all_tail_generation_or_actual_bad_ideal_claim=True)
 dump('two-column-finite-source-operator.json',operators)
 # Proposed ONE known-unit case only. No Hasse value, factor, worker or gcd.
 counts=[sum(a+e<=68-k for a in range(35)for e in range(52-k))for k in range(19)]
 maxlayer=max(counts);owned=2*(4811+190+maxlayer)+16000
 assert maxlayer==1667 and owned==29336
 proposal=dict(SOURCE_ONLY_NO_J_VALUES_OR_WORKER=True,
  minimal_case='ONE named-J valuation of J1, not arbitrary GCD; vJ(J1)=0 would exclude commonnamedJ, v>0 would give a knownunit reduction before two-tail complement algebra',
  J='5*a*d-6*e-4',root_d='(6e+4)/(5a)',root_LC='5a is already a FULLU unit; no b1/unknown pivot chart',
  primitive_over_QQa_e='gcd(5a,6e+4)=1, hence ordinaryJ valuation agrees with linear-d generic valuation',
  input='canonical4811 integers at originalD197, exponentsa<=16/e<=33/d<=18,total<=50',
  finite_layer='P_k=sum_{l>=k} c_ijl*binom(l,k)*(6e+4)^(l-k)*5^(18-l+k)*a^(i+18-l+k)*e^j',
  identity='[z^k]J1(a,e,(6e+4)/(5a)+z)=P_k/(D*5^18*a^18)',
  scope='k0..18 firstnonzero only; no evaluation at a numerical point, prime, dEuclid or quotient',
  layer_support_bounds=dict(a_max34=34,e_max51_minus_k=True,total_degree_max68_minus_k=True,all19_counts=counts,max1667=maxlayer),
  root_cache190=190,whole_owned29336=owned,proposed_owned34376=34376,
  height='every L1 layer/partial <=inputL1*15^18; coefficientwise multinomial bound; inputL1<2^301 gives372bits',
  source_integer_bits372=372,proposed_hard518=518,
  all_product_upper1200000=1200000,
  detailed_upper='4811*190=914090 rootproducts,<=2*4811*19=182818 scalings,rootcache/powers/binom<2000; all<1200000',
  no_extra_unknown_inversions=True,RAW_proposed8000000=8000000,soft45_hard60_RAM2GB=True,
  no_Jvaluation_GO_and_fresh_source_code_audit_REQUIRED=True)
 dump('ONE-known-J-case-PREFLIGHT-ONLY.json',proposal)
 pin(ROOT/f'research/notes/{ID}.md');pin(Path(__file__))
 r=dict(id=ID,all_passed=True,frozen=True,source_only=True,canonical_commonD_scheme_and_source_identity=True,
  no_new_polynomial_coefficients_factor_GCD_tail_values=True,
  proposed_J1_case_ONLY_unexecuted=True,input_sha256=pins,
  evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump('report.json',r);print(json.dumps({k:v for k,v in r.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

#!/usr/bin/env python3
"""SOURCE only: last named-B case with exact coefficient-coprimality, finite cache and budget; no B values."""
from pathlib import Path
import argparse,hashlib,json,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_Q1_B_HASSE_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 t=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):return json.loads(pin(ROOT/'reports'/st/name).read_text())
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 tw=load('S61_NONSQUARE_PREPOWER_TWO_TAIL_CANONICAL_SOURCE_20261105_v63');assert tw['all_passed']
 jp=load('S61_NONSQUARE_PREPOWER_Q1_J_ROOT_DHORNER_REVIEW_20261105_v63');assert jp['all_passed']
 load('S61_NONSQUARE_PREPOWER_Q1_J_ROOT_DHORNER_REVIEW_PUBLICATION_20261105_v63')
 fp=load('S61_NONSQUARE_PREPOWER_Q1_F_ROOT_HORNER_REVIEW_20261105_v63');assert fp['all_passed']
 norm=load('S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63');assert norm['all_passed']
 lit=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63/quotient-original-common-D.txt')
 assert sha(lit)==norm['evidence_sha256'][lit.name]
 hist=load('S61_NONSQUARE_PREPOWER_Q1_J_HASSE_WORKER_CUSTODY_20261105_v63','fixed-J-formula-height-lifetime-and-cursors.json')['source_d_histogram']
 assert sum(hist)==4811 and len(hist)==19
 # ONLY a tiny defining-factor identity, never a product involving Q1.
 b1={(1,0):50,(0,2):18,(0,1):-9,(0,0):-2};R={(0,1):-60,(0,0):40};w={(0,1):6,(0,0):1}
 identity={k:20*c for k,c in b1.items()}
 for x,c in R.items():
  for y,v in w.items():k=(x[0]+y[0],x[1]+y[1]);identity[k]=identity.get(k,0)+c*v
 identity[(1,0)]-=1000;assert not any(identity.values())
 lemma=dict(SOURCE_ONLY_NO_BVALUES=True,
  ring='ordinaryQQ[a,e,d]; actualsourceR_FULLU=QQ[a,e,d,1/(a*d*F*B*C*J)]',
  B='(50*a+18*e^2-9*e-2)*d+60*e-40=b1*d-R',
  b1='50*a+(3e-2)*(6e+1)',R='40-60e=-20*(3e-2)',
  leading_b1_is_NOT_assumed_FULLU_unit=True,
  primitive_QQa_e='onlypossible commonfactor with R is3e-2; b1 at e=2/3 equals50a notzero; hence coefficientgcd(b1,R)=1',
  primitive_ZZ='integercoefficients include-9 and-2; -1*(-9)+4*(-2)=1',
  ordinary_valuation='primitiveB linearind is irreducible overQQ[a,e,d]; firstnonzero Taylorcoefficient at d0=R/b1 in QQ(a,e)[d] gives exactvB in ordinaryring',
  no_b1_chart_discard='pseudo-Hasse numerators are fullpolynomial identities. Never infer pointwise solution coverage by settingb1 nonzero or invertb1 in actualFULLU adapter.',
  stronger_unit_identity='20*b1+R*(6e+1)=1000*a',
  integer_GL2=[['b1','-R'],['6e+1','20']],GL2det='1000a, alreadyFULLUunit',
  adjugate_substitution=['x=20X+R*Y','y=-(6e+1)X+b1*Y'],
  homogeneous_B_identity='b1*x-R*y=1000a*X',
  source_use='global homogeneous factor normalization overR[a^-1] available withoutinvertingb1; doesNOTprovebad-complementempty',
  valuation_upper16_by_a_degree='deg_a I=16 anddeg_a B=1, so ordinaryvB<=16; retainproposedorders0..18 as full d-degree recipe')
 dump('ordinary-B-valuation-and-no-unknown-chart.json',lemma)
 sizes=[19*(19-h)for h in range(19)];assert sum(sizes)==3610
 layers=[sum(a+e<=86-k for a in range(35)for e in range(70))for k in range(19)]
 assert max(layers)==2297
 root_R=2*sum(range(1,19));root_P=4*sum(j*j for j in range(1,19))
 root_S=sum((19-h)**2*(h+1)for h in range(19));assert(root_R,root_P,root_S)==(342,8436,13300)
 H=[sum(n*19*(19-l+k)for l,n in enumerate(hist)if l>=k)for k in range(19)]
 N=[sum(n for l,n in enumerate(hist)if l>=k)for k in range(19)]
 schedule=[root_R+root_P+root_S+sum(H[:k+1])+2*H[k]+sum(N[:k+1])+2*N[k]for k in range(19)]
 assert max(schedule)==12768337
 construct=2*(4811+3610+685+190)+16000
 value=2*(4811+3610+2297)+16000
 assert(construct,value)==(34592,37436)
 proof=dict(SOURCE_ONLY_NO_WORKER_OR_BVALUES=True,ONE_canonicalQ1_namedB_ONLY=True,
  input='4811 literalinteger numeratorI1 andoriginalD197, unchangedbody; NO newfactorvalue or sourcequery',
  fixed_root_cache='S_h(a,e)=R(e)^h*b1(a,e)^(18-h),h0..18',
  cache_support='S_h: a<=18-h,e<=36-h-2a; 19*(19-h) slots upper',all19cache_bounds=sizes,total_cache3610=3610,
  finite_formula='P_k=sum_{l>=k} I_ijl*choose(l,k)*S_(l-k)*a^i*e^j',
  identity='[z^k]Q1(a,e,R/b1+z)=P_k/(D197*b1^18)',
  denominator='D197 times SYMBOLIC b1^18, NOT a numericalvalue or assumption thatb1isFULLUunit',
  orders='0..18 FIRSTnonzero only; reversefirstlayer and reproduce in ONEaccumulator',
  support='a<=34,e<=69,totaldegree<=86-k; doNOTreplacee bound by69-k for this pseudo-Hasseformula',
  all19_layer_bounds=layers,maxlayer2297=2297,
  construction='R^h onevariablecache190; onecurrentP_t=b1^t andone nextP (atmost324+361); for t0..18 putS_(18-t)=P_t*R^(18-t) directlyintofinalrootmap. No fullroot-producttemp.',
  last_use='releaseALLRcache andcurrentP workspace afterALL19S ready BEFORE allocatingHasseaccumulator',
  input4811=4811,cache3610=3610,formula_owned37436=value,construction_owned34592=construct,
  owned='two cells percoefficient+packedkey, ALLmetadata/scratch16000; owned=2*(4811+3610+2297)+16000',
  metadata16000='choose19x19,all exponents/source/rootkeys/cursors/counts/l1/D,b1/Rdefinitions,numerictemporaryZs andprefix/outputscratch',
  height='L1I<2^301; for l<=18 sum_k choose(l,k)100^(l-k)79^(18-l+k)=79^(18-l)*179^l<=179^18<2^135. Allsourceweightedpartials<2^436.',
  root_cache_height120=120,source_main_height436=436,proposed_hard518=518,
  actual_index_schedule_including_threepasses=schedule,onepass_inner_upper_by_k=H,
  fixed_root_construction_products22078=root_R+root_P+root_S,
  ALL_product_upper12768337=max(schedule),proposed_product_guard13000000=13000000,
  runtime='13m finite integer multiply-add upper; CPUcompletion not measured. One hard45ssoft60shard guardedjob only afterfreshcode/rootaudit, saveexactprefixifSTOP; no retry promise.',
  prefix_FULLstate_bound4642060=(10718+16000)*170+100000,
  prefix_or_output_proposal5000000=5000000,prepared_inputs_codes_must_fit2800000=2800000,
  all_future_receipts200000=200000,RAW_proposed8000000=8000000,soft45_hard60_RAM2GB=True,
  no_Bvalues_GCD_tail_or_badlocus_and_new_budget_REQUIRED=True)
 dump('finite-cache-formula-lifetime-height-and-cost.json',proof)
 conditional=dict(SOURCE_ONLY_IF_futurevB_positive_NO_quotientcomputed=True,
  integer_Gauss='BprimitiveinZZ[a,e,d]; if I1 divisibleB^s rationally, quotientI1/B^s isinteger, so originalD197 canbe retained',
  no_unknown_LC_inverse='use ordinary multivariatelex a>e>d: LM(B)=50*a*d. Leadingcoefficient50 is a constant, not b1.',
  possible_algorithm='repeatedfixed6-termB integerlongdivision withleadingLC50; every integerdivision andallcoefficientmultiplyback mustbeverified; no b1 inverse or arbitrarypolyGCD',
  true_quotient_degrees_ifmultiplicity_s='a<=16-s,e<=33-2s,d<=18-s,totaldegree<=50-3s',
  no_ready_quotient_budget='actuals REQUIRED before derivingheight/live/output/worker gate; thislemma isnotpermission tocomputequotient',
  complement_scope='normalizingB changesnoFULLUzero set sinceBnamedunit; stillnotunitideal or actualbadFitting equality')
 dump('conditional-known-B-removal-method-ONLY.json',conditional)
 pin(ROOT/f'research/notes/{ID}.md');pin(Path(__file__))
 r=dict(id=ID,all_passed=True,frozen=True,source_only=True,
  ordinary_B_valuation_no_b1_chart_assumption=True,fixedcache_method=True,
  source_owned37436=37436,source436_hard518=True,ALL_products_upper12768337=12768337,
  no_B_Hasse_values_worker_quotient_GCD_or_tail=True,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-t)
 dump('report.json',r);print(json.dumps({k:v for k,v in r.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

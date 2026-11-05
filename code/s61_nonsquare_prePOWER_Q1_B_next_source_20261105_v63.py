#!/usr/bin/env python3
"""Publish Bvaluation4 primary; SOURCE only for distinct fivejet peer and fixed-B lex quotient."""
from pathlib import Path
from fractions import Fraction
import argparse,hashlib,json,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_Q1_B_NEXT_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):return json.loads(pin(ROOT/'reports'/st/name).read_text())
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 actual=load('S61_NONSQUARE_PREPOWER_Q1_B_HASSE_ENTRY_20261105_v63');assert actual['all_passed']and actual['result']['B_valuation']==4
 assert actual['result']['first_nonzero_terms']==1107
 zero=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_B_HASSE_ENTRY_20261105_v63/completed-zero-orders.txt');assert zero.read_text().splitlines()==[f'{i} ZERO'for i in range(4)]
 value=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_B_HASSE_ENTRY_20261105_v63/first-nonzero-B-Hasse.txt')
 assert sha(value)==actual['evidence_sha256'][value.name]
 source=load('S61_NONSQUARE_PREPOWER_Q1_B_HASSE_SOURCE_20261105_v63');assert source['all_passed']
 cp=load('S61_NONSQUARE_PREPOWER_Q1_B_HASSE_WORKER_CUSTODY_20261105_v63');assert cp['all_passed']
 hist=load('S61_NONSQUARE_PREPOWER_Q1_J_HASSE_WORKER_CUSTODY_20261105_v63','fixed-J-formula-height-lifetime-and-cursors.json')['source_d_histogram']
 norm=load('S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63');assert norm['all_passed']
 lit=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63/quotient-original-common-D.txt');assert sha(lit)==norm['evidence_sha256'][lit.name]
 def count(l,r):
  if r>18-l:return 0
  return sum(a+e<=86-3*l-r for a in range(35-l)for e in range(70-2*l))
 stages=[[count(l,r)for r in range(5)]for l in range(19)]
 assert sum(stages[0])==11295 and max(map(max,stages))==2297
 whole=2*(4811+11295+2297+685)+16000;assert whole==54176
 products=8436
 for l in range(17,-1,-1):
  products+=2*count(l+1,0)+hist[l]*(19-l)**2
  for r in range(4,0,-1):products+=2*count(l+1,r)+4*count(l+1,r-1)
 peer=dict(SOURCE_ONLY_NO_PEER_VALUES=True,original_input='canonical4811 integers/D197, notproducerS-cache; allnonnegativea/e andd0..18',
  independent_algorithm='H0=I18,H1..4=0; for l17..0, r4..1 descending: newHr=R*Hr+b1*H(r-1); newH0=R*H0+b1^(18-l)*I_l',
  proof='H_r is [z^r]sum_{j>=l}I_j*(R+b1*z)^(j-l)*b1^(18-j); at l0 ALLfivejets shareD*b1^18',
  R='40-60e',b1='50a+18e2-9e-2',
  no_binomial_Scache_or_nativeengine=True,
  required='H0,H1,H2,H3 exactlyZERO; H4 all1107 integercoefficients match afterinputfree; D197 andsymbolicb1^18 literalmatch; doNOTinvertb1 or discardb1zero',
  support='stage(l,r): a<=34-l,e<=69-2l,total<=86-3l-r; Hr=0 ifr>18-l',all19_fivejet_support_bounds=stages,
  jet_total_max11295=11295,one_temp2297=2297,Pcurrent_next685=685,
  lifetime='input4811+allfivejets+ONEtemp+current/nextb1power only; descendingr preserves oldlowerjet; no completefivejetduplicate/expected1107map',
  whole_owned54176=whole,proposed_owned60000=60000,metadata16000_included=True,
  height='weighted binomialsum79^(18-j)*179^(j-l) bounds allfivejets by L1I*179^18<2^436; conservativebeforecancel437,hard518',
  source_stored436=436,conservative_preop437=437,hard518=518,
  exact_stage_product_upper=products,proposed_all_product_upper3000000=3000000,
  RAW_proposed8000000=8000000,soft45_hard60_RAM2GB=True,
  no_Bquotient_GCD_otherfactor_tail_or_badlocus=True)
 assert products<3000000
 dump('independent-fivejet-dHorner-SOURCE.json',peer)
 Ns=[];degree=[]
 for j in range(5):
  a,e,d,t=16-j,33-2*j,18-j,50-3*j
  nn=sum(x+y+z<=t for x in range(a+1)for y in range(e+1)for z in range(d+1));Ns.append(nn);degree.append([a,e,d,t])
 assert Ns==[10013,8400,6970,5712,4615]
 L=(1<<301)-1;bounds=[]
 for j in range(4):
  m=16-j;factor=sum((Fraction(129,50)**h)/50 for h in range(m));q=Fraction(L)*factor
  Q=(q.numerator+q.denominator-1)//q.denominator;RR=L+179*Q
  bounds.append(dict(stage=j+1,input_L1_upper=str(L),factor_num=str(factor.numerator),factor_den=str(factor.denominator),quotient_L1_upper=str(Q),quotient_bits=Q.bit_length(),all_partial_remainder_bits=RR.bit_length(),scalarproduct_bits=(60*Q).bit_length()));L=Q
 verify=((1<<301)-1)+179**4*L
 assert verify.bit_length()==386
 sourceproducts=6*sum(Ns[1:])+2*sum(Ns[1:])+225*Ns[-1]+5400
 assert sourceproducts==1249351
 quot=dict(SOURCE_ONLY_NO_QUOTIENT_VALUES=True,independent_Bvaluation_peer_REQUIRED=True,
  target='J2=I1/B4 integer polynomial; normalizedtailQ1/B4=J2/originalD197',
  primitive_Gauss='B integerprimitive: -1*(-9)+4*(-2)=1, so rationalB4divisibility impliesintegerJ2',
  fixed_B_sixterms='50ad+18e2d-9ed-2d+60e-40',
  order='lex a>e>d via keya*1048576+e*1024+d overprovedbounds; LM(B)=50ad, integerLC50',
  algorithm='FOUR standard sparseintegerlongdivisions by fixedB; require leadinga,d>=1 andleadingcoef divisible50; subtract all6terms; eachstage FULLremainderZERO; swapcompletedquotientintoremainder fornextstage',
  no_b1_inverse_or_polyGCD=True,
  true_stage_degree_boxes_a_e_d_total=degree,stage_support_bounds=Ns,
  justification='truequotient degrees subtractB-degree(a1,e2,d1,total3); everypartialR=stageinput-B*(removedtruequotient terms), so allpartialR stayinoldstagebox',
  L1_recurrence='B=50ad+U(e,d),||U||1=129,||B||1=179. quotientbound <=inputL1/50 *sum_{h=0..m-1}(129/50)^h, m=deg_a stageinput; divisionbyd is monomialreindexing only.',
  exact_L1_bounds=bounds,final_full_B4multiplyback_partial_bits386=386,hard_allinteger518_UNCHANGED=True,
  source_tmp='largeststageR363bits,finalintegerB4multiplyback386; quotient*50reversechecks andall6termproducts covered',
  stage_owned=[2*(Ns[j]+Ns[j+1])+16000 for j in range(4)],whole_owned52826=52826,
  final_verify_owned45706=2*(Ns[0]+Ns[-1]+225)+16000,
  lifecycle='readoriginal4811 onceintoremainder, no secondinputvector; R/Q maps only. Stagehashes streamed notfullintermediatefiles. After4stages buildB4cache<=225 withQfinal, releaseoldfactorpower before rereading originalI intoverifyremainder.',
  final_identity='ALLcoeff B4*J2=I1, remainderZERO; no discardednegativeexponent/denominator orfieldcomponent',
  all_product_upper1249351=sourceproducts,proposed_product_guard1300000=1300000,
  output_terms_upper4615=4615,one_final_output_ONLY=True,
  full_prefix_text_bound5950210=(18413+16000)*170+100000,prefix_reserve6000000=6000000,
  future_prepared_all_inputs_history_codes_mustfit1800000=1800000,receipt_reserve200000=200000,RAW_proposed8000000=8000000,
  source_goal='onecompleteexplicitknownunitB4removal usefulforcommonkappa/twotailcomparison; noordinaryfactorGCD/unitideal/complementclaim',
  new_code_source_audit_and_ONE_GO_REQUIRED=True,soft45_hard60_RAM2GB=True)
 dump('fixed-primitive-B4-lex-quotient-SOURCE-ONLY.json',quot)
 dump('primary-publication-binding.json',dict(primary_Bvaluation4=True,firstlayer1107=True,
  primary_sha256=sha(ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_B_HASSE_ENTRY_20261105_v63/report.json'),
  independent_peer_NOT_run_by_this_source=True,no_Bquotient_otherfactor_or_GCD=True))
 pin(ROOT/f'research/notes/{ID}.md');pin(Path(__file__))
 result=dict(id=ID,all_passed=True,frozen=True,source_only=True,actual_primary_publication=True,
  independent_dHorner_NOT_run=True,fixed_B4_quotient_NOT_computed=True,
  peer_owned54176=54176,quotient_owned52826=52826,quotient_hard518=True,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump('report.json',result);print(json.dumps({k:v for k,v in result.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

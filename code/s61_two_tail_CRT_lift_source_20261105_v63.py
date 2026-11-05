#!/usr/bin/env python3
"""SOURCE ONLY for one fixed40-prime two-tail CRT candidate and conditional exact verification."""
from pathlib import Path
from math import gcd,isqrt,prod
import argparse,hashlib,json,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_TWO_TAIL_CRT_LIFT_SOURCE_20261105_v63'
PRIMES=[32003,32009,32027,32029,32051,32057,32059,32063,32069,32077,32083,32089,32099,32117,32119,32141,32143,32159,32173,32183,32189,32191,32203,32213,32233,32237,32251,32257,32261,32297,32299,32303,32309,32321,32323,32327,32341,32353,32359,32363]
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):return json.loads(pin(ROOT/'reports'/st/name).read_text())
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 diag=load('S61_TWO_TAIL_FIXED_KRONECKER_ENTRY_20261105_v63');assert diag['all_passed']and diag['result']['gcd_degree']==3899
 peer=load('S61_TWO_TAIL_FIXED_KRONECKER_QUOTIENT_REVIEW_V2_20261105_v63');assert peer['all_passed']
 oldLM=load('S61_NONSQUARE_PREPOWER_TWO_TAIL_LM_ROOT_REVIEW_20261105_v63');assert oldLM['all_passed']
 seed=pin(ROOT/'reports/S61_TWO_TAIL_FIXED_KRONECKER_ENTRY_20261105_v63/monic-gcd-mod32003.txt');assert sha(seed)==diag['evidence_sha256'][seed.name]
 support=[]
 with seed.open()as f:
  assert int(f.readline())==3899
  for line in f:
   E,c=map(int,line.split());assert 0<c<32003;a,r=divmod(E,646);e,d=divmod(r,19);support.append((a,e,d))
 assert len(support)==94 and max(support)==(6,1,4)
 assert[max(v[i]for v in support)for i in range(3)]==[6,9,4]and max(map(sum,support))==11
 assert max(2*a+e for a,e,d in support)==13
 assert[v for v in support if 2*v[0]+v[1]==13]==[(6,1,4)]
 spbytes=''.join(f'{a} {e} {d}\n'for a,e,d in sorted(support)).encode()
 inputs=[]
 for st,magic,n in [('S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_ENTRY_20261105_v63','T0C4QV001',4500),('S61_NONSQUARE_PREPOWER_Q1_B4_QUOTIENT_ENTRY_20261105_v63','Q1B4QV001',2317)]:
  r=load(st);assert r['all_passed'];p=pin(ROOT/'reports'/st/'quotient-integer-common-D.txt');assert sha(p)==r['evidence_sha256'][p.name]
  L=0;degree=-1;lc=0;maxW=-1;count=0;maxs=[0,0,0];total=0
  with p.open()as f:
   assert f.readline().strip()==magic;D=int(f.readline());assert D>0 and D.bit_length()==197;assert int(f.readline())==n
   for line in f:
    key,c=map(int,line.split());assert c;a,x=divmod(key,1048576);e,d=divmod(x,1024)
    E=646*a+19*e+d
    if E>degree:degree=E;lc=c
    L+=abs(c);count+=1;maxW=max(maxW,2*a+e);maxs=[max(s,t)for s,t in zip(maxs,[a,e,d])];total=max(total,a+e+d)
  assert count==n
  inputs.append(dict(path=str(p.relative_to(ROOT)),magic=magic,terms=n,D197=D,Kronecker_degree=degree,Kronecker_LC=lc,LC_bits=abs(lc).bit_length(),L1_bits=L.bit_length(),maxW=maxW,degrees=maxs,total_degree=total))
 assert[x['Kronecker_degree']for x in inputs]==[10544,7975]and[x['maxW']for x in inputs]==[42,35]
 gamma=gcd(inputs[0]['Kronecker_LC'],inputs[1]['Kronecker_LC']);assert gamma>0 and gamma.bit_length()==226
 assert len(PRIMES)==40 and len(set(PRIMES))==40 and PRIMES[0]==32003 and max(PRIMES)<32768
 primeledger=[]
 for p in PRIMES:
  assert all(p%d for d in range(2,isqrt(p)+1))
  residues=[x['Kronecker_LC']%p for x in inputs]+[gamma%p,inputs[0]['D197']%p];assert all(residues)
  primeledger.append(dict(p=p,certification='allintegerdivisors2..isqrt testednonzero',sqrt_bound=isqrt(p),LC0_LC2_gamma_D_residues=residues,reuse_saved_seed=(p==32003)))
 modulus=prod(PRIMES);assert modulus.bit_length()==599
 bound=dict(SOURCE_ONLY_NO_NEW_GCD_CRT_OR_FACTOR_VALUES=True,
  verified_modular_seed_degree3899=True,seed_review_sha256=sha(ROOT/'reports/S61_TWO_TAIL_FIXED_KRONECKER_QUOTIENT_REVIEW_V2_20261105_v63/report.json'),
  upper_bound_lemma='Any primitivechar0commonpolynomialH has e-degree<=25,d-degree<=14, soKroneckerimage isinjective/nonconstant; Gaussgivesintegerquotients and bothmappedinputLCs survive32003, henceHmappedLCsurvives. Hbar dividesverifiedgcd3899, soKrondegreeH<=3899.',
  consequence='deg_a H<=6 because7*646=4522>3899; ifdeg_a6 then19e+d<=23 onleadinga6terms, e<=1 andd<=4ife1.',
  no_common_factor_existence_irreducibility_or_maximality_claim=True,
  fixed_seed_support94=dict(sha256=hashlib.sha256(spbytes).hexdigest(),box=[6,9,4],total11=11,unique_maxW13=[6,1,4]),
  named_candidate_inventory='supportprofile differs fromeachof6literalnamedunitfactors; ifexactfactorofN2 passes, acceptedN2coprimeall6 ensuresHnonunitinBASE6. No generalcandidatefactorsearch or unitpolynomialproducts computed.')
 dump('verified-seed-and-commonfactor-degree-upper-bound.json',bound)
 schedule=dict(SOURCE_ONLY=True,inputs=inputs,gamma=gamma,gamma_bits226=226,
  gamma_purpose='gamma=gcd(integerLCN0,integerLCN2)>0 scalesmonicunivariategcd tointegercandidate; removeCONTENT before anyintegerfactor division',
  fixed40_primes=primeledger,max_modulus_bits599=599,
  existing32003_seed_reused_NO_rerun=True,new_prime_jobs_atmost39=True,
  degrees='inputweightedLCs mustremainnonzero anddegrees10544/7975 remain; no parameterpoint evaluation',
  degree_filter='target3899 fixed. gcddegree>3899:markbadandSKIP withinthisfixedlist; degree<3899:STOPfrontier/noresetornewcandidate; degree3899:monicnormalizeandCRT. InputLC/gamma/Dpolebad:skiponlythatprelistedprime.',
  no_additional_prime_reset_farm_or_perprime_candidate_retry=True,
  one_reconstruction_afterfixed_list_ONLY=True,
  crt='store3900implicitdegreecoefficientsC in[0,M). r_i=gamma modp*monicgcd_i; t_i=(r_i-C_i modp)*(M^-1 modp); C_i+=M*t_i. MchangesONLYafterALL3900coefficientscommit.',
  leading_reconstruction_gate='requirefinalM>2gamma, thenbalancedleadingcoeffmustEXACTgamma; ifnotSTOPcandidatefrontier',
  complete_bound_not_claimed='599bitmodulus isNOT a proven coefficientheight bound fortrueintegergcd. Fixed40 isoneattempt; correctnesscomesONLYfromfinalexactmultivariatedivisions/multiplybacks.',
  all_modular_products_cap6000000000=6000000000,per_newprime_cap150000000=150000000,
  standard_Euclid_source_upper='same2vectorinplaceEuclid bound112.5M perprime fromverifiedproducerdegrees;39*150M=5.85B, plusCRT/inverse/controlmargin<6B',
  CRT_bigint_products_upper153000=153000,
  scalar_bounds='p<32768,modularproducts<2^30,decimalparser<2^19; M/C/scaledCRTcoefficients<=600bits,absoluteintegercontents<=600bits',
  runtime='NO deterministiccompletionpromise; hard45ssoft60shard2GBsharedclock; exactfullprefixonSTOP; noautoretry')
 dump('fixed-primes-degreefilter-and-CRT-normalization-SOURCE.json',schedule)
 Hshape=[(a,e,d)for a in range(7)for e in range(10)for d in range(5)if a+e+d<=11 and 2*a+e<=13 and(2*a+e<13 or(a,e,d)==(6,1,4))]
 assert len(Hshape)==207
 R0=19*sum(43-2*a for a in range(17));Q0=15*sum(30-2*a for a in range(11))
 R2=15*sum(36-2*a for a in range(13));Q2=11*sum(23-2*a for a in range(7))
 assert(R0,Q0,R2,Q2)==(8721,3300,4680,1309)
 height=dict(SOURCE_ONLY_ADAPTIVE_CANDIDATE_GATE=True,
  first_normalize='balancedCRTcoefficients; g=positiveintegergcd ofALL3900coefficients; exactdivideallbyg andbackcheck; choosepositiveleadingcoefficient. This isINTEGERcontent, NOT multivariatefactorsearch.',
  mandatory='NEVER tryintegerdivisionbygamma-scaledgcd beforeprimitivecontentremoval',
  selected_candidate_only='decodeall3900exponentsa646+e19+d. AcceptONLYoneprimitivecandidate whoseactualsupportliesin fixedshapea<=6/e<=9/d<=4,total<=11,ALLotherW<13; uniqueW13term(6,1,4)nonzero;actualdegreeprofile[6,9,4,total11]. OtherwiseSTOPfrontier,nostartnewfamily.',
  shape_terms_upper207=207,may_include_restored_zero_mod32003_coefficients=True,
  LM_proof='candidatelexLM(6,1,4)checkedfromactualselectedcoefficients, NOT inferredchar0frommodularseed; alsouniqueleadingforweightedorder(W=2a+e,lexa/e/d).',
  key='integerdivision usesW*2^26+legacy(a*2^20+e*2^10+d), notKronecker646/19/1; e canexceed33inpartialremainders. legacykey<2^26andallcoordinateaddsno1024carry.',
  partial_false_candidate_invariants=dict(R0=['a<=16','d<=18','W<=42','e<=42'],Q0=['a<=10','d<=14','W<=29'],R2=['a<=12','d<=14','W<=35','e<=35'],Q2=['a<=6','d<=10','W<=22']),
  why='HtailsstrictlydropW andhavea<=6,d<=4. Acceptedleadingdivisionplusmonomialshiftnevernigrowsa/d,andallremainderproductsdropW; appliesEVENFALSEcandidate. Neverassumeordinarye/totaldegreedropfromsuccessfulfactor.',
  support_bounds=dict(R0=R0,Q0=Q0,R2=R2,Q2=Q2),
  input_L1_bits=[inputs[0]['L1_bits'],inputs[1]['L1_bits']],quotient_weight_depths=[29,22],
  actual_height_test='c=absLC(H)>0,T=L1(H)-c,BH=bitlengthL1(H). r=minnonnegativeintegerwithT<=c*2^r. ForinputBf andm29/22: Bq=max(0,Bf-floorlog2c+ceilLog2(m+1)+m*r), BR=max(Bf,Bq+BH)+1. Computeusingboundedinteger/logmetadata, nohugepowerallocation.',
  complete_condition='proceedintegerworkONLYifBOTH BR<=8192; elseSTOPsavedcandidate. ThegeometricweightedL1boundcoversALLpartialR/Qandfreshmultiplyback,evenwrongcandidateuntilfirstdivisionfailure.',
  useful_sufficient_gate='candidatecoeff<=600bits/207terms=>BH<=608; r<=250 givesN0BR<=8162<8192 usingBf298,soafeasiblepositivegate withoutassumingitwillpass',
  new_algorithm_cap_NOT_old518_capraise=True,new_integer_coefficient_hard8192=8192,
  division='HprimitiveZZ; GaussmeansQQdivisibilityimpliesintegerquotient. RequireLMHdividescurrentLMandLCintegerexactdivision; ANYfailureSTOPcandidatefrontier/nofactorclaim. OneR/Qmap,no productpolynomial.',
  proof_complete='candidateNONCONSTANT + bothintegercofactors + ALLcoeffH*Q0=N0andH*Q2=N2 matchedtoORIGINAL3variableinputs. UnivariategcdorunivariatedivisionisNEVERenough.',
  final_scope='propercommonfactorofselectedtwo tailsONLY; noirreducibility/gcdmaximality/badFitting/full57/239/exhaustionclaim',
  bigint_products_upper2100000=2100000,
  bigint_count='2*207*(3300+1309)+2*(3300+1309)+3900+153000=2074244<2.1M,plusmargin; coefficientproductcostseparatefrom6Bmodularuintproducts')
 dump('primitive-candidate-and-false-candidate-safe-exact-verification.json',height)
 field_owned=18521+3900+2*6817+16000
 integer_owned=2*(8721+3300+207+6817)+16000
 assert(field_owned,integer_owned)==(52055,54090)
 int_prefix=(8721+3300)*1040+207*96+6817*48+65536
 field_prefix=18521*4+3900*80+6817*48+65536
 cofs=(3300+1309)*1040
 assert int_prefix==12914464 and field_prefix==778836 and cofs==4793360
 raw=3*2**20+16*2**20+5*2**20+2**20+2**20;assert raw==26*2**20<32*2**20
 budget=dict(SOURCE_ONLY_NO_WORKER_OR_VALUES=True,
  old_private_jobs_UNCHANGED=True,new_fixed_pair_factor_algorithm=True,
  immutable_input_lifetime='store6817sparserecords(legacykey+integercoefficient)ONCE,sourcecoeff<=300bits. Reuseforprimefieldsandtwomultivariableidentities,no sourceoracle/newpoint.',
  field_lifetime='2uintvectorscapacity18521+CRT3900+inputpairs13634+metadata16000; normalizecurrentgcdINPLACE, CRTscaleonthefly,nothirdmonicvector',
  field_owned52055=field_owned,
  candidate_lifetime='inplacebalancedCRT/contentnormalization;extractonlyactualacceptedshapeH<=207;publishCRTsnapshotthenreleaseBOTHfieldvectorsanddenseCRTBEFOREintegerR/Q',
  integer_lifetime='input6817pairs+R8721+Q3300+H207pairs+metadata16000; releaseverifiedcofactor0mapafterbinarypublicationbeforeseconddivision',
  integer_owned54090=integer_owned,proposed_owned60000=60000,globalcap150000_UNCHANGED=True,
  binary_records='coefficient8192bits:key64+sign8+len16+1024bytes<=1040;input301bits<=48bytes;candidate600bits<=96bytes;denseCRT<=80bytes/coeff',
  typed_prefix_metadata='<=64KiB:40primeledger*10uint64,40primes,32cursors,16bigscalars<=8192bits,sourceSHApathsanderrors; NOT16000arbitrary8kbigints',
  integer_binary_prefix_upper12914464=int_prefix,field_prefix_upper778836=field_prefix,prefix_reserve16MiB=16*2**20,
  required_prefix_cursor='mode/primeindex/degreefilter,Euclidpass/kshift/j_next/currentfactor/inverse; CRTi_next/oldM/newM/acceptedledger;contenti/g/sign;candidatefactor andadaptiveheight;whichinput/weightedR/Q/activequotient/Hterm_next;allopsandfield/statearrays',
  CRT_interrupted_update='MnotcommitteduntilALLcoefficientsdone; prefixhasoldM,newM,nexti andmixedold/newC; neverrepeatfinishedcoefficients',
  allstate='onebinaryprefixonlyonSTOP,novalueclaim/automaticretry; noatomicduplicate;IOfailurecannottriggersecondfullprefix',
  cofactor_outputs_upper4793360=cofs,cofactor_reserve5MiB=5*2**20,
  CRT_candidate_normalization_reserve1MiB=2**20,allreports_logs_reserve1MiB=2**20,
  all_prepared_history_inputs_codes_mustfit3MiB=3*2**20,
  RAW_union_upper26MiB=raw,RAW_hard32MiB=32*2**20,
  source_modular_product_cap6000000000=6000000000,source_bigint_product_cap2100000=2100000,
  soft45_hard60_RAM2GB=True,
  SOURCE_feasible_conditional_NOT_completion_guarantee=True,fresh_worker_custody_ROOT_source_code_audit_and_ONE_GO_REQUIRED=True)
 dump('whole-ownership-binary-prefix-and-RAW-budget.json',budget)
 pin(ROOT/f'research/notes/{ID}.md');pin(Path(__file__))
 result=dict(id=ID,all_passed=True,frozen=True,source_only=True,no_newprime_gcd_CRT_candidate_or_exactdivision_run=True,
  fixed40_primes_seed32003_reused=True,independent_seed_peer_pinned=True,
  conditional_8192_verification_height=True,proposed_owned60000=60000,RAW_union26MiB=raw,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump('report.json',result);print(json.dumps({k:v for k,v in result.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

#!/usr/bin/env python3
"""Complete degree-two source coverage, 617 kernel, and r-tail corollary."""
import json,math,resource,signal
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_LOWER_RANK2_HF_TAIL_SOURCE_20261105_v63'
def image(A,b):
 z={}
 for j,c in b.items():
  for i,x in A[j].items():
   v=z.get(i,Q(0))+c*x
   if v:z[i]=v
   else:z.pop(i,None)
 return z
def main():
 signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));pins={}
 def read(p,keys=None):pins[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 def passed(p):v=read(p);assert v['all_passed'];return v
 # Three disjoint value checks are all needed for the exact70 statement.
 for p in ['reports/AMBIENT_B18_ELL2_T2_NATIVE70_PREFIX14_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_NATIVE70_WITNESS563_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_NATIVE70_CONTINUATION55_V2_20261105_v63_REVIEW/report.json']:passed(p)
 nb='reports/AMBIENT_B18_ELL2_T2_NATIVE70_CONTINUATION55_V2_20261105_v63';pr=passed(nb+'/report.json');assert pr['actual_lower4_by70_rank']==2
 all70=read(nb+'/actual-native70-complete-lower-images-from-cached15-and-new55.json')['computed_prefix'];assert [r['case619'] for r in all70]==list(range(549,619))
 val={r['case619']:list(map(Q,r['actual_lower4_values'])) for r in all70};w,q=val[563],val[609];assert w==list(map(Q,['-1000/3','-212/3','-19000','-3968']))
 # ALL old549 are covered: numeric548 plus the independently proven order508.
 passed('reports/AMBIENT_B18_ELL2_T2_OLD548_FIRST_ONLY_20261105_v63_REVIEW/report.json');passed('reports/AMBIENT_B18_ELL2_T2_OLD_CASE508_ORDER_20261105_v63_REVIEW/report.json')
 old=read('reports/AMBIENT_B18_ELL2_T2_OLD548_FIRST_ONLY_20261105_v63/actual548-first-only-lower-prefix-and-witness-if-any.json');assert old['all548_zero'] and {r['case619'] for r in old['computed_source_prefix']}==set(range(549))-{508}
 for r in old['computed_source_prefix']:assert not any(map(Q,r['actual_lower4_values']));val[r['case619']]=[Q(0)]*4
 val[508]=[Q(0)]*4;assert set(val)==set(range(619))
 # The complete G3 kernel393 has been checked in disjoint active/nonactive families.
 for p in ['reports/AMBIENT_B18_ELL2_T2_G3_PLANE_QUADRATIC_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_NONACTIVE171_SOURCE_PREFLIGHT_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_VALUES_20261105_v63_REVIEW/report.json']:passed(p)
 kb='reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63';A=read(kb+'/phase1-known-G3linear-fixed-source-rank.json');B=read(kb+'/actual-known-bottom-first-rank-and-corrected-source-kernel.json');passed(kb+'_REVIEW/report.json')
 GA=list(map(sp,A['actual_projected_G3linear_columns']));selA=A['selected_columns'];facA=list(map(sp,A['all418_original_basis_factors']));BB=list(map(sp,B['all642_bottom27_columns']));selB=B['selected_bottom_columns'];facB=list(map(sp,B['all642_bottom_original_basis_factors']));src=B['all_zero_bottom_and_Y0_source_combinations']
 assert len(GA)==418 and len(selA)==25 and sum(j<238 for j in selA)==16 and len(BB)==642 and len(selB)==23 and len(src)==619
 # Fixed span receipts, not a new rank search. Source rank minors are already independently checked.
 basA=[GA[j] for j in selA];basB=[BB[j] for j in selB]
 assert all(image(basA,c)==v for c,v in zip(facA,GA,strict=True));assert all(image(basB,c)==v for c,v in zip(facB,BB,strict=True))
 freeA=[j for j in range(418) if j not in selA];freeB=[j for j in range(642) if j not in selB];assert len(freeA)==393 and len(freeB)==619
 for j in freeA:
  c={j:Q(1)}
  for k,x in facA[j].items():c[selA[k]]=c.get(selA[k],Q(0))-x
  assert not image(GA,{i:x for i,x in c.items() if x})
 for j,r in zip(freeB,src,strict=True):
  c={j:Q(1)}
  for k,x in facB[j].items():c[selB[k]]=c.get(selB[k],Q(0))-x
  assert sp(r['source642_combination'])=={i:x for i,x in c.items() if x} and r['full56_Xgamma1_equals_Xp0a1'] and not image(BB,sp(r['source642_combination']))
 hf=passed('reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_RANK_20261105_v63/report.json');assert hf['actual_HF_T2_through1']==[6,6];passed('reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_RANK_20261105_v63_REVIEW/report.json');passed('reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63_REVIEW/report.json')
 assert 418+642==1060 and 25+23==48 and 393+619==1060-48==1012
 lin=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json',['degree2_linear_standard_module_monomials','allowed_variables_per_target']);assert lin['allowed_variables_per_target']==[[],[],[],[7,8],[7,8],[7,8]]
 passed('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_GB_COVER_20261105_v63_REVIEW/report.json');passed('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_GB_COVER_20261105_v63/report.json')
 Q2=list(map(sp,read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE171_SOURCE_PREFLIGHT_20261105_v63/exact171-decomposition-162-zero-9-lower-first-source-input.json',['Q_two_upper_quadratic_columns'])['Q_two_upper_quadratic_columns']));assert len(Q2)==2
 # Positive bounds before computing the new constant617 source basis.
 positive=0
 def ex(x):x=Q(x);return (x.denominator,abs(x.numerator))
 def pos(a,b):D=math.lcm(a[0],b[0]);return D,a[1]*(D//a[0])+b[1]*(D//b[0])
 def product(a,b):return a[0]*b[0],a[1]*b[1]
 for x in val.values():
  for coeff in [[Q(-59,23200),Q(-1,464)],[Q(-53,20000),Q(1,80)]]:
   D,U=pos(product(ex(coeff[0]),ex(x[0])),product(ex(coeff[1]),ex(x[1])));positive=max(positive,D.bit_length(),U.bit_length())
 assert positive<128
 kernels=[]
 for j in range(619):
  a=Q(-59,23200)*val[j][0]-val[j][1]/464;b=Q(-53,20000)*val[j][0]+val[j][1]/80
  assert val[j]==[a*x+b*y for x,y in zip(w,q,strict=True)]
  if j not in [563,609]:kernels.append({'free_case619':j,'subtract_actual_basis563609':[str(a),str(b)],'constant_source_recipe':'corrected_relation[j]-a*corrected_relation[563]-b*corrected_relation[609] in the SAME coherent full27/99 frame; not basis-upper zero filling'})
 assert len(kernels)==617
 u=[Q(-59,23200)*x-Q(53,20000)*y for x,y in zip(w,q,strict=True)];v=[-x/464+y/80 for x,y in zip(w,q,strict=True)];assert u==[Q(1),Q(0),Q(57),Q(-9,50)] and v==[Q(0),Q(1),Q(0),Q(57)]
 bounds=passed('reports/AMBIENT_B18_ELL2_HONEST_T2_LENGTH_BOUNDS_20261105_v63/report.json');assert bounds['actual_T2_length_interval']==[17,20];passed('reports/AMBIENT_B18_ELL2_HONEST_T2_LENGTH_BOUNDS_20261105_v63_REVIEW/report.json')
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'complete1012-kernel617-and-filtered-r-tail-source-proof.json'
 proof={'degree2_completeness':'For any formal combination of the1060 minimal relations R, itsdegree1zero condition is a0 inker(R1). Itsdegree2 class modulo L48*m is thena0R2; a1R1vanishes in theknownNF. The block rank25+23 decomposition gives anindependentcomplete constantkernel393+619. The619 basis includesactualfixed25G3corrections, evenwhen theircoefficientvalues aretypedunknown; triangularprojectiontofree642coordinatesprovesindependence. All393 G3quadraticclasseszero, so all newdegree2classescome from619.',
 'HF2':'The complete lowerprojection of the619 image hasrank2 andkernel617. ItsupperkernelQ2 hasdimension2. Thereforedimgr2=6-(2+rank(upper2 restrictedtokerLOWER617))=4-rank(upper2|617), lyingin{2,3,4}. No newuppervalue isassigned.',
 'filtered_r_tail':'The completeknownhomogeneousGBhasstandardsbottom3*(r^d,r^(d-1)ell) ford>=1. Hence multiplicationbyrfromcoverdegree2todegree3isiso andsurjectsaftertheactualgradedquotient. Thusm3T=rm2T+m4T. Nakayamaonm3T/rm2Tgivesm3T=rm2T; inductiongivesm^dT=r^(d-2)m2Tford>=2. Thefinitek[[rA]]tailm2T haslengthL-12in[5,8], minimalgeneratorcountdimgr2in[2,4], andnonincreasingHilbertlayers. SinceL>=17anddimgr2<=4, m3Tcannotvanish; NakayamaalsogivesHF3>=1. ThisrAisnottheellDVRparameter.',
 'source617_scope':'The617 displayedconstantcombinations formthefullkerLOWER basis. Allactual27/99/G3 corrections, b1,A2,M1 andnative/honesthighertermsmustbecombinedwiththesameconstantsbeforecoefficientqueries. Basis563/609 upper2 valuesremainunknown; neverzero-fillthem.'}
 ev={'complete_constant_initial_kernel_dimensions':[393,619,1012],'canonical_lower_basis2':list(map(lambda a:list(map(str,a)),[u,v])),'all617_kernel_lower_constant_source_recipes':kernels,'kernel_weight_positive_prior_bits':positive,'actual_HF2_interval':[2,4],'actual_HF3_interval':[1,4],'tail_length_interval':[5,8],'tail_parameter':'rA, not ell','exact_filtered_identity':'m^d T2=rA^(d-2)*m^2 T2 for d>=2','HF2_exact_formula':'4-rank(upper2_on_kerLOWER617)','proof':proof,'no_upper617_values_or_total_T2_or_H0_decision':True};ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
 rr={'id':STEM,'all_passed':True,'full_degree2_constant_kernel_coverage':1012,'full_lower619_rank':2,'remaining_upper_source_kernel_count':617,'actual_HF2_interval':[2,4],'actual_HF3_interval':[1,4],'actual_tail_length_interval':[5,8],'actual_T2_length_or_H0_exponents_not_decided':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)}};(out/'report.json').write_text(json.dumps(rr,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'kernel617':len(kernels),'weight_prior':positive}))
if __name__=='__main__':main()

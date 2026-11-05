#!/usr/bin/env python3
"""Pure source preflight for four lower quadratic mixed619 channels."""
import hashlib,json,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_SOURCE_PREFLIGHT_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sp(v):return {int(i):Q(c) for i,c in v}
def main():
 t=time.monotonic();inputs={}
 def read(p):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return json.loads((ROOT/p).read_text())
 def select(p,k):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,k)
 kernel=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json');records=kernel['all_zero_bottom_and_Y0_source_combinations'];assert len(records)==619
 mini=read('reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json')
 assert mini['remaining_free3_originalV_rows']==[23,29,30] and mini['remaining_free3_Y_rows']==[1,36,38]
 qdata=read('reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63/whole-cache-eight-psi-four-quadratic-classes-source-input.json');cols=list(map(sp,qdata['Q_two6_columns']))
 duals=[{2:Q(1),4:Q(4)},{3:Q(1),5:Q(4)},{0:Q(11),4:Q(270)},{1:Q(11),5:Q(270)}]
 assert all(sum((x*c.get(i,Q(0)) for i,x in a.items()),Q(0))==0 for a in duals for c in cols)
 # The fixed4minor on columns(2,3,0,1) has determinant121.
 assert all(duals[i].get(j,0)==(1 if i==k and i<2 else 11 if i==k else 0) for i in range(4) for k,j in enumerate([2,3,0,1]))
 full='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json';f=select(full,['PiY3_from102_columns','actual_full_G3linear102_by418_columns']);PiY=list(map(sp,f['PiY3_from102_columns']));assert [PiY[i] for i in (1,36,38)]==[{0:Q(1)},{1:Q(1)},{2:Q(1)}]
 assert kernel['all_kernel_full56_Xgamma1_Xp0a1_identities_verified']
 total_source=sum(len(r['original768_combination']) for r in records);total_a1=sum(len(r['actual_a1_27_times9']) for r in records)
 assert all(r['full56_Xgamma1_equals_Xp0a1'] for r in records)
 for r in records:
  assert all(0<=j<768 for j,x in r['original768_combination']) and all(0<=j<27*9 for j,x in r['actual_a1_27_times9'])
 rank=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json');assert len(rank['selected_columns'])==25
 lin=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json');standards=lin['degree2_linear_standard_module_monomials'];assert [r[0] for r in standards]==[3,3,4,4,5,5]
 for stem in ['AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63_REVIEW','AMBIENT_B18_ELL2_T2_FIRST_HILBERT_RANK_20261105_v63_REVIEW','AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63_REVIEW','AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_GB_COVER_20261105_v63_REVIEW']:
  r=read(f'reports/{stem}/report.json');assert r['all_passed']
 proof={
 'Pi3':'last3 coordinates of fixedbottom30inverse composed with restriction to free30Vindices. Pi3Xp0=0, and it gives actualpost27bottom3. Never replace by three rawsource rows withoutthismap.',
 'B99':'selected99 constantY0basis withfreeE(1,36,38) complement. Apply fixed99receipts to eachunitcolumn; residualEfreePiY0 remains, sourcecoefficients give B99. Verify B99Ystar0=I99 and B99Efree=0. No newGaussian.',
 'Hstar1':'Pi3 times eachselected99 originalsource first vector; subtractionofXp0a1 hasPi3zero. Cachedfull9source firsts suffice.',
 'post27':'For619, Xgamma0=Ygamma0=0 andXgamma1=Xp0a1. Z=Xgamma-Xp*a1 startsdegree2. Hence G2linear*Ybar1=M0*(Xgamma2-Xp1*a1): M1 cancels andG2quadY0=0. Unknowna2drops M0Xp0 andPi3Xp0.',
 'G3clearing':'Let J be full102x9 firstcolumns ofthefixed25G3basis, F25itsfixedleftmap onprojected27. Foractual619, PiY0Ybar1 liesinthatspan byacceptedall619firstzeroquotient. Define C=I-JF25PiY0 onlinearcoefficientspace. Itsuse isessential; Ybar1 isnotsetzero.',
 'lower_formula':'Correctedlowerquadratic=Pi3*(Xgamma2-Xp1*a1)-Hstar1*B99*C(Ybar1), withordinaryrr/rell convolution. The secondtermalsoincludes theactualquadraticlowerpartofthe25G3clearings. NounknownupperPiY1 entersbecauseallupperquadraticNF isannihilatedbythefourduals.',
 'source_order':'FourlowerchannelscanbeseekedfromX2only via adjoints ofG2linear. Top2channelsremainseparateandmayneedX3. NativeS0zero mustberecheckedforeachnewfold; oldhonestsecondunknownsstayformal untilsensitivitydescent. No newphysical/gauge/sourcecoefficients.',
 'next_small_input':'Firstbuildfourknownlinear covectors onYbar1 usingPi3/B99/Hstar1/J/F25; retainall918 linearYcoeff slots. Try52 rr/rell outputadjoints ofG2linear against204 planeYslots onlyifall714zslots vanish; otherwise freezeobstructionorexpandonlywithanewsourcebudget. Noadjoint solvehere.',
 'scope':'Anyverifiednonzero correctedlowerclass isatruequadraticinitialrelationanduppercoverinput. Allzeroonlyclosesfourlowerchannelsfor619,notremainingtop2orcompleteHF2.'}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'mixed619-lower4-source-operator-proof-and-next-input-schema.json';ep.write_text(json.dumps({'source619_count':619,'total_original768_coefficients':total_source,'total_saved_a1_coefficients':total_a1,'known_lower_duals4_by6':[[[i,str(x)] for i,x in sorted(a.items())] for a in duals],'all_duals_annihilate_Q2':True,'fixed4minor':121,'primitive_source_operator_recipe':proof,'weighted_covectors_adjoint_matrix_or_values_constructed':False,'next_prior_height_and_owned_gate_required':True},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'pure_source_X2_only_four_lower_channel_preflight':True,'unknown_weighted_functions_or_values_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'source619':619,'source_coeffs':total_source,'a1_coeffs':total_a1}))
if __name__=='__main__':main()

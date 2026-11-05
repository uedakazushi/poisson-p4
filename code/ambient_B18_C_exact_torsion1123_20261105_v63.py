#!/usr/bin/env python3
"""Source-bound exact comparison-C theorem after independent tauAB witness."""
import hashlib,json,signal,time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_EXACT_TORSION1123_20261105_v63'
PATHS={
 'value':'reports/AMBIENT_B18_C_FULL_TAUAB12_VALUES_20261105_v63/report.json',
 'value_data':'reports/AMBIENT_B18_C_FULL_TAUAB12_VALUES_20261105_v63/all-full-tauAB12-values-and-fixed-known-rank-C4-certificates.json',
 'value_peer':'reports/AMBIENT_B18_C_FULL_TAUAB12_VALUES_20261105_v63_REVIEW/report.json',
 'source_peer':'reports/AMBIENT_B18_C_FULL_TAUAB12_SOURCE_20261105_v63_REVIEW/report.json',
 'eta':'reports/AMBIENT_B18_C_FULL_ETA4_CACHE_VALUES_20261105_v63/report.json',
 'eta_peer':'reports/AMBIENT_B18_C_FULL_ETA4_CACHE_VALUES_20261105_v63_REVIEW/report.json',
 'reduction':'reports/AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63/report.json',
 'reduction_peer':'reports/AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63_REVIEW/report.json',
 'quartic':'reports/AMBIENT_B18_C_FULL_QUARTIC_BOUND_SOURCE_20261105_v63/report.json',
 'quartic_peer':'reports/AMBIENT_B18_C_FULL_QUARTIC_BOUND_SOURCE_20261105_v63_REVIEW/report.json',
 'C4':'reports/AMBIENT_B18_C4_HB_RESIDUAL_SOURCE_20261105_v63/report.json',
 'C4_peer':'reports/AMBIENT_B18_C4_HB_RESIDUAL_SOURCE_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 started=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()}
 d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()};assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 assert d['value']['exact_comparison_C_e3_proved']and d['value_data']['exact_comparison_C_e3_proved']
 assert d['value_data']['C4_rank4_certificate_if_positive']['universal_block_minor']=='-1331/15360'
 assert d['eta']['actual_full_C_HF_degree3']==6
 assert d['reduction']['unconditional_m3C_equals_ell_m2C_proved']
 assert d['C4']['unconditional_C4_m5_zero_proved']
 pattern=[1,1,2,3]
 dlengths=[[n,4*n+sum(min(n,a)for a in pattern)]for n in range(1,9)]
 assert dlengths[:4]==[[1,8],[2,14],[3,19],[4,23]]
 data={
  'scope':'ACTUAL original comparison module C=coker(D0→E4), over samephysicalregular9 R0 andD=QQ[[ell]]. Not originalH0(D0) or K torsion.',
  'exact_noncanonical_D_module_type':{'free_rank':4,'torsion_exponents':[1,1,2,3]},
  'comparison_C_DVR_torsion_length':7,
  'source_witness':'ActualtauAB18 quarticresidual1331/9216 andfixed2/C4fourminor−1331/15360, checkedindependentlyfrom39-symbolpolynomial. This givesC4upper23, matchingcurrentDVRlower23, hencee3.',
  'exact_full_m_adic_HF':[3,6,6,6,'4 at everydegree>=4'],
  'full_Hilbert_series':'(3+3t-2t^4)/(1-t)',
  'full_HF_proof':'ExactfullHF3=6 fromcomplete155andeta4zero. KnownactualHB+tauABrank2 removes2fromfullquarticcover6, soHF4<=4. Reductionm5=ellm4andrank4 givesHF4>=4; m4C isD-free4. Allhigherm-powersareellmultiplesofthisfree4, soallhigherHF4.',
  'm_powers_D_types':{'m2C':'D4 ⊕ (D/ell²)^2','m3C':'D4 ⊕ k²','m4C':'Dfree4','m_nC_n_ge4':'ell^(n-4)*m4C'},
  'comparison_C_torsion_annihilators':['ell³*Tor_D(C)=0','m⁴*Tor_D(C)=0'],
  'torsion_annihilator_scope':'No m³-annihilation or originalH0annihilator isclaimed.',
  'ell_power_quotient_lengths_first8':dlengths,
  'ell_power_quotient_length_formula':'4n+2+min(n,2)+min(n,3); =4n+7 forn>=3',
  'm_adic_quotient_lengths':'lengthC/m^n=3,9,15 forn1,2,3; =4n+5 forn>=4',
  'actual_C3_HF':[3,6,6,4,0],'actual_C3_length':19,
  'actual_C4_HF':[3,6,6,6,2,0],'actual_C4_length':23,
  'canonical_genuine_torsion_target':{
   'module':'C[ell]≅ker(E:C4→C4)/im(E³:C4→C4)',
   'E_ranks':[23,15,9,4,0],'kerE_dimension':8,'imE3_dimension':4,'genuine_dimension':4,
   'proof':'InD4+tors1123,thefreeD/ell4blockscontributeequal4-dim kerEandimE³; torsionblocksareallkilledbyell³andcontributetheir4socles. Intrinsically, ifellx=ell4y thenx−ell³y isgenuineell-torsion; itsclassmodimE³ isindependentofchoices. No arbitraryfinitekerE identification.',
   'source_value_gate':'Actual23actions/basis/genuine4representativesandconnectingvalues areNOTcomputedbythiscorollary.'},
  'no_original_H0_length_K_torsion_exponents_nilcube_or_exact_nilindex_transfer':True,
  'no_fresh_coefficient_job_or_source_replay':True,
  'next_source_goal':'ConstructactualfiniteC4/genuine4coefficientinterfaceandoriginalcomparisonconnectingmapwithsameR0sourcebinding. DoNOTcontinueprimitiveδ4foralready-solvedC-exponentgoal.',
 }
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'exact-comparison-C-torsion-full-HF-and-genuine-target-corollary.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'EXACT_COMPARISON_C_TORSION1123_FULL_HF3666444_GENUINE_TARGET_SOURCE_COROLLARY',
  'actual_comparison_C_torsion_exponents':[1,1,2,3],'actual_comparison_C_torsion_length':7,'actual_C4_length':23,
  'original_H0_or_nilcube_determined':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)

if __name__=='__main__':main()

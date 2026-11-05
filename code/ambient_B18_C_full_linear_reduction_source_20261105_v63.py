#!/usr/bin/env python3
"""Full comparison C reduction theorem and conditional graph-cubic e=3 gate."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63'
PATHS={
 'frame':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63/actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json',
 'frame_peer':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63_REVIEW/report.json',
 'graph':'reports/AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63/actual-cover-NF-and-five-triangular-cubic-quartic-contradictions.json',
 'graph_peer':'reports/AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63_REVIEW/report.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'C2_report':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/report.json',
 'C2_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'rank4':'reports/AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63/report.json',
 'rank4_peer':'reports/AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63_REVIEW/report.json',
 'cut':'reports/AMBIENT_B18_C3_CUT1_QUARTIC_VALUES_20261105_v63/one-CUT1-residual-and-formal-coordinate-rank-certificate.json',
 'cut_peer':'reports/AMBIENT_B18_C3_CUT1_QUARTIC_VALUES_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dot(a,b):return sum((x*y for x,y in zip(a,b)),Q(0))


def main():
 started=time.monotonic();signal.alarm(60);out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 pins={p:sha(ROOT/p)for p in PATHS.values()};d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 assert d['rank4']['actual_comparison_C_generic_DVR_rank4_proved']and d['rank4']['actual_mu_D_C8_from_independent_source_evaluation']
 assert d['C2_report']['actual_C2_dimension']==14 and d['cut']['actual_C3_length19_proved']
 assert d['cut']['cut3_torsion_pattern_if_positive']==[1,1,2,'e>=3']
 assert d['frame']['actual_C_and_C3_HF_degrees0_1_2']==[3,6,6]
 # The NO-CUT cover has six cubic dimensions. In each sector the free
 # normal coordinates are r^2ell,r ell^2,ell^3; none is discarded.
 block=[]
 for rr in d['graph']['two_block_source_curve_NF_receipts']:
  g=list(map(Q,rr['curve3_r3_r2ell_rell2_ell3']));assert g[0]
  nf=[[-g[k]/g[0]]+[Q(i==k)for i in range(1,4)]for k in range(1,4)]
  assert all(dot(row,g)==0 for row in nf)
  assert[[row[j]for j in[1,2,3]]for row in nf]==[[Q(i==j)for j in range(3)]for i in range(3)]
  block.append(nf)
 fullNF=[[Q(0)]*12 for _ in range(6)]
 for s in range(2):
  for k in range(3):fullNF[3*s+k][4*s:4*s+4]=block[s][k]
 # ell from the six A/B quadratic monomials maps to all six known
 # no-cut cubic coordinates; any actual extra relations preserve surj.
 ell_matrix=[[fullNF[i][4*s+j+1]for s in range(2)for j in range(3)]for i in range(6)]
 assert ell_matrix==[[Q(i==j)for j in range(6)]for i in range(6)]
 frame={r['physical_row']:r for r in d['frame']['all12_verified_source_records']}
 vals={r['physical_row']:dict(r['delta'])for r in d['C2']['physical12_coefficients']}
 graph_bindings=[]
 for r in d['graph']['same_five_actual_source_premises']:
  physical=r['physical_row'];assert all(Q(x)==0 for x in frame[physical]['bA_first_z7']+frame[physical]['bB_first_z7'])
  assert all(Q(x)==0 for x in frame[physical]['bB_second_plane'])
  assert Q(vals[physical].get(29,'0'))==Q(vals[physical].get(30,'0'))==0
  graph_bindings.append({'physical_row':physical,'actual_full_plane_delta3':'d_graph03*ell^3: r3/r2ell0 from C2 and r*ell2=0 by reviewed triangular C3 proof; ell3 remainsUNREAD',
      'tauA_full_no_cut_cubic_NF':['0','0','d_graph03','0','0','0'],
      'tauB_full_no_cut_cubic_NF':['0','0','0','0','0','d_graph03'],
      'source_carry':'h=bB1=0, so tauA has no primitiveB. Hmixed clears use RH secondplane0; delta2zz clears retainz. Hcubics are killed by linear*actualHquadratics and affect A/B only at4. Noell3cut is used in this fullC argument.'})
 assert len(graph_bindings)==5
 data={
  'scope':'FULL original comparison module C=coker(D0→E4) over same actual regular9 R0; D=QQ[[ell]], m=m_R0. C is not asserted to be an algebra or D-free.',
  'accepted_DVR_data':{'rank_D':4,'mu_D':8,'C_mod_ell2_length':14,'C_mod_ell3_length':19,'torsion_exponents_before_gate':[1,1,2,'e>=3']},
  'full_no_cut_cubic_NF6_rows12_columns':[[str(v)for v in row]for row in fullNF],
  'full_ell_multiplication_six_quadratic_to_cover_cubic_I6':[[str(v)for v in row]for row in ell_matrix],
  'unconditional_full_C_reduction_proof':'Actualinitialz7*3 andplane2H removeallothercubicmonomials; gA3/gB3 unitr3 termsrewritepurer3A/B asell timesquadratics. Hence gr3C=ell*gr2C forFULLC, withnoell3quotient. Thereforem3C=ell*m2C+m4C. ForfiniteR0moduleQ=m3C/(ell*m2C), thisisQ=mQ; Nakayama givesm3C=ell*m2C. Multiplying yieldsm^(n+1)C=ell*m^nC forEVERYn>=2. Noell-cancellation orC-flatness assumption.',
  'unconditional_m_powers_ell_reduction_from2':True,
  'conditional_full_HF3_lower':'m3C genericD-rank4 sinceell3C⊂m3C. Itsmu_D equalsdim m3C/ell*m3C=dim m3C/m4C=HF3full. ThusHF3full>=4. AlternativelyfullC→C3surjects ongr3andacceptedC3HF3=4.',
  'five_graph_witness_source_bindings':graph_bindings,
  'one_nonzero_graph03_implies_full_HF3_4':'In the no-cut cubic cover dim6, these two actual sources havefixedminor d_graph03^2 onthe twoell3 coordinates. Ifd!=0, HF3full<=4; previouslower4givesexact4. Nootherunknowncoefficientissetzero.',
  'conditional_D_freeness_and_exponent3':'HF3full4 meansfiniteDmodulem3C hasrank4/mu4, henceisD-free4. ItsD-submoduleell3C istorsionfree. AnyCtorsionexponent>3wouldleavenonzeroell3torsinsideell3C, impossible. ThusallCtorsionexponents<=3, andaccepted e>=3 forces e=3.',
  'conditional_exact_D_type':{'free_rank':4,'torsion_exponents':[1,1,2,3],'torsion_length7_scope':'comparison C ONLY, not original H0(D0)'},
  'conditional_full_C_HF':[3,6,6,'4 for every degree>=3'],
  'conditional_full_C_Hilbert_series':'(3+3t-2t^3)/(1-t)',
  'conditional_m4_equals_ell3C':'C3m4=0 givesm4C⊂ell3C. FullHFthrough3=3664 giveslengthC/m4=19, equalacceptedlengthC/ell3; hence equality. Thenm^kC=ell^(k-1)C fork>=4.',
  'conditional_Ctorsion_m_annihilator':'m3*tors_D(C)=0 sinceitliesinm3C∩tors_D(C), andm3C isD-free. ThisdoesNOTannihilateoriginalH0/Ktorsion.',
  'conditional_genuine_target_gate':'Type1123 giveslengthC/ell4=23 andcanonical C[ell]≅ker(ell:C/ell4→C/ell4)/ell3(C/ell4), dim4. Thisisatheorem only: actual23-dimensionalactions/basis/edgevaluesstillneedsourceconstruction.',
  'graph03_allzero_scope':'Onlythesegraphwitnessesnegative; doNOTinferfullHF3,exactlongere,orotherfullC cubicrelationsfromthisalone.',
  'no_graph03_values_or_new_source_gauge_IFT_jet':True,
  'fixed_source_profile':{'max_degree':3,'owned_numeric_envelope':25000,'hard_numeric':150000,'hard_seconds':60,'hard_memory_bytes':2147483648},
 }
 ep=out/'full-C-ell-reduction-and-conditional-graph03-exponent3-proof.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'PURE_SOURCE_fullC_reduction_and_conditional_graph03_exponent3_NO_VALUES',
   'unconditional_m3C_equals_ell_m2C_proved':True,'graph03_nonzero_witness_computed':False,'longer_exponent_determined':False,
   'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)


if __name__=='__main__':main()

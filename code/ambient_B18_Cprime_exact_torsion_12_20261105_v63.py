#!/usr/bin/env python3
"""Dependency-bound DVR/Hilbert corollary; no source or value computation."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_CPRIME_EXACT_TORSION_12_20261105_v63'
INPUTS=[
 'reports/AMBIENT_B18_REDUCED_COMPARISON_GENERIC_RANK2_20261104_v63_REVIEW/report.json',
 'reports/AMBIENT_B18_REDUCED_COMPARISON_CBAR_VALUES_20261104_v63_REVIEW/report.json',
 'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63_REVIEW/report.json',
 'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63_REVIEW/report.json',
 'reports/AMBIENT_B18_CPRIME_RIGHT85_VALUES_20261105_v63/report.json',
 'reports/AMBIENT_B18_CPRIME_RIGHT85_VALUES_20261105_v63_REVIEW/report.json',
 'reports/AMBIENT_B18_CPRIME_CUBIC_ELL3_VALUES_20261105_v63_REVIEW/report.json']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);D={p:json.loads((ROOT/p).read_text())for p in INPUTS};assert all(d['all_passed']for d in D.values())
 value=D[INPUTS[4]];assert value['actual_right85_rank']==1 and value['positive_Cprime_ell3_length9_and_long_exponent2']
 lengths=[2*n+min(n,1)+min(n,2)for n in range(1,13)];assert lengths[:3]==[4,7,9]
 hf=[3,3,3]+[2]*12;numerator=[hf[0]]+[hf[j]-hf[j-1]for j in range(1,len(hf))];assert numerator[:4]==[3,0,0,-1]and not any(numerator[4:])
 proof={'same_DVR':'D=QQ[[ell]], ell=5*(sA-rA); rA itself is not the DVR parameter',
  'actual_Cprime_D_module_type':'D^2 + D/(ell) + D/(ell^2), noncanonical D-module isomorphism only',
  'actual_ell_power_quotient_lengths_n1_to12':lengths,'actual_Cprime_torsion_D_length':3,
  'm3_equals_ell3_proof':'Knowninitial+genuineV2u1 hasfull27cubicrank afteraddingell³, so m³(Cprime/ell³)=0. The automaticell³Cprime⊂m³Cprime then gives equality.',
  'reduction_proof':'ell³Cprime⊂ell*m²Cprime⊂m³Cprime=ell³Cprime, hence m³=ell*m². Multiplying bym gives m^{j+1}=ell*m^j for everyj>=2; combiningm³=ell³ givesm^j=ell^j forj>=3.',
  'full_actual_m_adic_Hilbert_function':'(3,3,3,2,2,...)','Hilbert_series_numerator':[3,0,0,-1],'Hilbert_series_denominator':'1-t','Hilbert_Samuel_multiplicity':2,
  'genuine_torsion_extraction':'Cprime[ell] canonically equals ker(ell:Cprime/ell³→Cprime/ell³)/ell²(Cprime/ell³); dimension2. It does not equal the full finite kernel before quotientingell².',
  'new_actual_target_action_or_connecting_values':False,'original_H0_D0_length_nilcube_or_generators_claimed':False,
  'Cprime_torsion_not_original_H0_D0':True,'Cprime_m3_not_zero_globally':True}
 ep=out/'exact-exponents-full-Hilbert-and-genuine-target-corollary.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'PROVED_dependency_corollary_Cprime_exact_DVR_exponents_and_full_HF','Cprime_DVR_torsion_exponents':[1,2],'Cprime_torsion_length':3,'Cprime_ell3_quotient_length':9,'Cprime_full_HF_prefix':[3,3,3,2,2,2],'genuine_Cprime_ell_torsion_dimension':2,
  'input_sha256':{p:sha(ROOT/p)for p in INPUTS},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'new_source_gauge_jet_or_scalar_job':False}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()

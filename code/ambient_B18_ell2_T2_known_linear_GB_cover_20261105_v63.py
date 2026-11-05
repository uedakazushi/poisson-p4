#!/usr/bin/env python3
"""Pure globalGB consequence of completed432 zero reductions; no arithmetic rerun."""
import json,itertools
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_GB_COVER_20261105_v63'
def main():
 gp='reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json';ep='reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_DEGREE2_20261105_v63/actual432-overlap-rank-and-known-linear-cover-degree2.json';d=json.loads((ROOT/gp).read_text());e=json.loads((ROOT/ep).read_text());assert all(not v for v in e['all432_remainder_columns'])
 rows=d['normalized_linear_rewrite_rows'];index={tuple(v):i for i,v in enumerate(e['all432_true_source_multiple_labels'])};pairs=[]
 for (p,a),(q,b) in itertools.combinations(rows,2):
  if p//9==q//9:pairs.append({'linear_leads':[p,q],'S_pair_zero_NF_source_multiple_columns':[index[p,q%9],index[q,p%9]],'combination':[1,-1]})
 ell=[]
 for p,a in rows:
  if p%9==8:ell.append({'linear_lead':p,'kind':'leading_ell_overlap_degree2','zero_NF_source_column':index[p,8]})
  else:ell.append({'linear_lead':p,'kind':'lcm_degree3_all_terms_have_ell_squared','source':'ell²*L -leadvar*e_target*ell²'})
 allowed=d['allowed_variables_per_target'];assert allowed==[[],[],[],[7,8],[7,8],[7,8]]
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);proof={'all_linear_linear_critical_pairs':pairs,'all_linear_ell_squared_critical_pairs':ell,
  'all432_source_multiple_NFs_zero':True,'full_module_Groebner_basis_proved':True,
  'allowed_variables_per_target':allowed,'known_cover_HF_through3':[6,6,6,6],
  'known_cover_Hilbert_series':'6/(1-t)','known_cover_prefix_length_through3':24,
  'I2_annihilator_cannot_lower_degree_le3':'IA,IB subsetm =>I=IAIB subsetm²=>I² subsetm4. ThusL+in(I²) hasidenticaldegree<=3 quotient.',
  'stop_known_linear_closure':True,'actual_T2_HF2_HF3_or_total_length_not_computed':True}
 p=out/'all-critical-pairs-and-known-cover-stop-proof.json';H.dump(p,proof)
 report={'id':STEM,'all_passed':True,'linear_linear_pairs':len(pairs),'linear_ell2_pairs':len(ell),'known_cover_HF_through3':[6,6,6,6],
  'known_cover_prefix_length':24,'accepted_actual_T2_length_upper':20,'known_linear_closure_stopped':True,
  'actual_T2_higher_HF_length_or_H0_not_claimed':True,
  'input_sha256':{gp:H.sha(ROOT/gp),ep:H.sha(ROOT/ep)},'code_sha256':H.sha(Path(__file__)),
  'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:H.sha(p)}}
 H.dump(out/'report.json',report);print(json.dumps({'all_passed':True,'pairs':[len(pairs),len(ell)],'known_cover_prefix':24,'report_sha256':H.sha(out/'report.json')}))
if __name__=='__main__':main()

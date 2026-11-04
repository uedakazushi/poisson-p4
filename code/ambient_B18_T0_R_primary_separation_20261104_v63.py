#!/usr/bin/env python3
"""Proof-gate receipt for canonical removal of the R-primary CM quotient."""
import hashlib,json
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T0_R_PRIMARY_SEPARATION_20261104_v63'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def main():
    out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
    gates=[
      ROOT/'reports/AMBIENT_B18_R_ZERO_SLICE_IMAGE_20261104_v63/report.json',
      ROOT/'reports/AMBIENT_B18_R_N_PROJECTION_NORM_20261104_v63/report.json',
      ROOT/'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_REDUCED_INTERSECTION_DVR_GLUING_20261103_v63/report.json',
      ROOT/'reports/AMBIENT_CA_CB_MINIMAL_LENGTHS_AND_R_CM_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_T0_COKERNEL_Z_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_T1_T2_QUADRATIC_PHI_20261104_v63_REVIEW/report.json']
    data=[json.loads(p.read_text()) for p in gates]
    assert all(r['all_passed'] for r in data)
    assert data[0]['actual_R0_uniformizer_row']==112
    assert data[1]['projection_norm_first_nonzero_degree']==4
    assert F(data[1]['projection_norm_coefficient_T4'])==-81600
    report={'id':STEM,'all_passed':True,'proof_type':'Noetherian_saturation_and_DVR_module_argument_no_new_CAS',
      'coefficient_ring':'actual_R0_regular9',
      'actual_scalar_F':'v^2-nu*z1^2 from accepted N same-R0 binding',
      'actual_scalar_U':'centered physical rho112',
      'actual_R0_R_image':'eight regular graph equations of accepted unit20 map',
      'canonical_F_torsion_submodule':'S_F=union_n(0:_{D0}F^n)',
      'R_primary_quotient':'B_R=D0/S_F',
      'R_primary_quotient_is_1D_CM_and_free_rank2_over_QQ_U':True,
      'R_primary_quotient_abstract_algebra':'QQ[[U,epsilon]]/(epsilon^2)',
      'R0_embedding_of_epsilon_or_saturation_exponent_computed':False,
      'T0_R_quotient':'T0/(0:_{T0}F^infinity)',
      'T0_R_quotient_is_unspecified_finite_colength_ideal_in_B_R_and_free_rank2_over_QQ_U':True,
      'socle_exact_reduction':'H0_m(D0)=H0_m(T0)=H0_m(S_F)=H0_m(T0 intersect S_F)',
      'residual_S_F_finite_curve_supported_or_socle_free_claimed':False,
      'exact_T0_annihilator_saturation_exponent_or_D_depth_claimed':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in gates+[Path(__file__),ROOT/'research/notes'/f'{STEM}.md']}}
    p=out/'report.json';p.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'report_sha256':sha(p)}))
if __name__=='__main__':main()

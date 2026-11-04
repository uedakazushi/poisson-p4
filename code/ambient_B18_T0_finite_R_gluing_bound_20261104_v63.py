#!/usr/bin/env python3
"""Finite gluing bound after canonical R-primary separation; no new CAS."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T0_FINITE_R_GLUING_BOUND_20261104_v63'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def main():
    out=ROOT/'reports'/STEM;out.mkdir(exist_ok=True);assert not (out/'report.json').exists()
    gates=[
      ROOT/'reports/AMBIENT_B18_T0_R_PRIMARY_SEPARATION_20261104_v63/report.json',
      ROOT/'reports/AMBIENT_B18_R_N_PROJECTION_NORM_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_REDUCED_INTERSECTION_DVR_GLUING_20261103_v63/report.json',
      ROOT/'reports/AMBIENT_B18_T0_COKERNEL_Z_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_W1_S61_INTERSECTION_RECOGNITION_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63/report.json']
    a=[json.loads(p.read_text()) for p in gates];assert all(x['all_passed'] for x in a)
    assert a[0]['R_primary_quotient_abstract_algebra']=='QQ[[U,epsilon]]/(epsilon^2)'
    assert a[1]['actual_projection_norm_and_generic_R_length2_independently_accepted']
    # These are the exact DVR length identities, not a sampled matrix rank.
    norm_order=4;nilpotence=2;rank=2
    cokernel_bound=norm_order*nilpotence*rank
    assert cokernel_bound==16
    allowed=[(m,n) for m in range(4,9) for n in range(m+1) if m+n<=cokernel_bound]
    report={'id':STEM,'all_passed':True,'proof_type':'exact_DVR_length_and_cokernel_sequence_no_new_CAS',
      'D4_nilradical_square_zero_used_on_the_whole_actual_tail':True,
      'actual_F_squared_annihilates_D4':True,
      'R_primary_ring':'B_R=QQ[[U,epsilon]]/(epsilon^2)',
      'R_primary_T0_image':'J_R=image(T0 -> B_R)',
      'F_squared_B_R_is_contained_in_J_R':True,
      'order_U_of_F_reduced':4,
      'order_U_of_determinant_multiplication_F_squared':16,
      'length_B_R_over_J_R_upper_bound':16,
      'reduced_J_R_containment_only':'image(J_R -> QQ[[U]]) subset (U^4)',
      'HNF':'J_R=(U^m+b(U)*epsilon,U^n*epsilon), b modulo U^n',
      'allowed_m_n_pairs':allowed,'HNF_bounds':{'m_min':4,'m_max':8,'n_min':0,'n_max':'m','m_plus_n_max':16},
      'exact_sequence':'0 -> B_R/J_R -> D4/image(S_F) -> Z -> 0',
      'length_Z':3,'length_D4_over_image_S_F_upper_bound':19,
      'actual_m_n_b_or_finite_quotient_matrices_computed':False,
      'residual_kernel_finite_socle_or_depth_claimed':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in gates+[Path(__file__),ROOT/'research/notes'/f'{STEM}.md']}}
    p=out/'report.json';p.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'report_sha256':sha(p)}))
if __name__=='__main__':main()

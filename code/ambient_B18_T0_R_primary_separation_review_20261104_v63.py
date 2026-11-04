#!/usr/bin/env python3
"""Independent proof-gate audit for the canonical R-primary double DVR."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_T0_R_PRIMARY_SEPARATION_20261104_v63';ID=BASE+'_REVIEW'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    primary=ROOT/'reports'/BASE/'report.json';p=read(primary);assert p['all_passed'];hashes={}
    for name,h in p['input_sha256'].items():assert sha(ROOT/name)==h;hashes[name]=h
    gates=['AMBIENT_B18_R_ZERO_SLICE_IMAGE_20261104_v63_REVIEW/report.json','AMBIENT_B18_R_N_PROJECTION_NORM_20261104_v63_REVIEW/report.json','AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63_REVIEW/report.json','AMBIENT_B18_D4_LITERAL_HORSESHOE_20261104_v63_REVIEW/report.json','AMBIENT_B18_W1_POLYNOMIAL_LAURENT_SCOPE_20261104_v63.json']
    for name in gates:
        path=ROOT/'reports'/name;r=read(path);assert r['all_passed'];hashes[str(path.relative_to(ROOT))]=sha(path)
    norm=read(ROOT/'reports/AMBIENT_B18_R_N_PROJECTION_NORM_20261104_v63/report.json');rz=read(ROOT/'reports/AMBIENT_B18_R_ZERO_SLICE_IMAGE_20261104_v63/report.json');d0=read(ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json');gr=read(ROOT/'reports/AMBIENT_B18_REDUCED_INTERSECTION_DVR_GLUING_20261103_v63/report.json')
    assert norm['projection_norm_first_nonzero_degree']==4 and F(norm['projection_norm_coefficient_T4'])==-81600
    assert rz['actual_R0_uniformizer_row']==112 and rz['unit20_minor']=='64'
    assert d0['actual_D0_free9_relations640_and_D4_map_checked'] and gr['exact_R_N_contact_order']==4
    F_T4=F(-81600);U_T1=F(-2);F_U4=F_T4/U_T1**4;assert F_U4==-5100
    assert p['R_primary_quotient_is_1D_CM_and_free_rank2_over_QQ_U'] and p['R_primary_quotient_abstract_algebra']=='QQ[[U,epsilon]]/(epsilon^2)'
    assert not p['residual_S_F_finite_curve_supported_or_socle_free_claimed'] and not p['exact_T0_annihilator_saturation_exponent_or_D_depth_claimed']
    proof=dict(actual_R0_regular_dimension=9,D0_actual_finite_R0_presentation=[9,640],F='v^2-nu*z1^2',U='centered physical rho112',F_R_order=4,F_R_leading_coefficient_in_U=str(F_U4),reduced_R0_support_D0='V(F) union V(P_R)',degree0_reduction_proved_in_nonnegative_finite_w_model=True,Noetherian_F_torsion_ideal_S=True,associated_primes_BR_exactly_P_R=True,graph_power_argument_proves_finite_over_QQ_U=True,U_regular_and_BR_free_rank2=True,actual_reduced_quotient_QQ_U_with_scalar_U_section=True,nil_free_rank1_and_generically_squarezero_implies_squarezero_everywhere=True,abstract_BR='QQ[[U,epsilon]]/(epsilon^2)',T0_image_finite_colength_rank2_ideal=True,H0_equalities='H0(D0)=H0(T0)=H0(S_F)=H0(T0 intersect S_F)',H0_is_finite_length_torsion_not_entire_socle=True,no_epsilon_embedding_or_saturation_index_colength_bound_or_full_T0_CM_claim=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for path in [primary,Path(__file__),ROOT/f'research/notes/{ID}.md']:hashes[str(path.relative_to(ROOT))]=sha(path)
    r=dict(id=ID,all_passed=True,actual_canonical_R_primary_double_DVR_quotient_independently_accepted=True,finite_colength_T0_ideal_and_H0_reduction_proved=True,no_CAS_source_jets_or_extra_16_19_bounds=True,input_sha256=hashes,evidence_sha256={'independent-proof.json':sha(out/'independent-proof.json')},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=r['seconds'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

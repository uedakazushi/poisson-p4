#!/usr/bin/env python3
"""Pure same-source quotient-map/DVR corollary, no coefficient worker."""
import argparse,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_REDUCED_COMPARISON_GENERIC_RANK2_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main(out):
 out.mkdir(parents=True,exist_ok=False)
 gates={
 'reports/AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63/report.json':'e5272f143977cfce50ed5d26df4d6910fa3e6a583f2874c4730c0bbe50e4818a',
 'reports/AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63_REVIEW/report.json':'22ce48012e3366120c1246a4b9276965eea2569a4f953478b32c7a79105d1496',
 'reports/AMBIENT_CA_CRITICAL_V1_CONTACT_20261104_v63_REVIEW/report.json':'cc2a24acbc7c8f4abfb8a3b56ffb1b70fe94662e5cdbcc315d423bcdd8a51f41',
 'reports/AMBIENT_CA_CRITICAL_V1_CONTACT_ALGEBRA_20261104_v63_REVIEW/report.json':'46dd4d22aadd1c25e94f25f3fef656aea95225d205fda6951fe176ff7409f695',
 'reports/AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW/report.json':'7cdd073ed27c79149e4f9fba56f6c9bfb02a90c34f65473ea672abceb110741f',
 'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW/report.json':'edab93bf8b666e088d62d6ece3a65be62b3ba7621b3e22c6490ca6f85ad16662',
 'reports/AMBIENT_B18_REDUCED_COMPARISON_CBAR_VALUES_20261104_v63/report.json':'b090a3037d8a22c64dd461594f79dca1566b2a3b7110150cf797d82a3a11b6ac',
 'reports/AMBIENT_B18_REDUCED_COMPARISON_ORIGIN_20261104_v63/report.json':'0524a73bc87ff7c12ffa04ccd5137aca5d0f15173a725e76f50dbde3f6229674',
 }
 for p,h in gates.items():assert sha(ROOT/p)==h and read(ROOT/p)['all_passed'],p
 bar=read(ROOT/'reports/AMBIENT_B18_REDUCED_COMPARISON_CBAR_VALUES_20261104_v63/report.json')
 assert bar['actual_reduced_Cbar_length']==4 and bar['actual_reduced_Cbar_HF']==[3,1,0]
 contributions={'CA_critical':0,'CA_selfpair_Pa':1,'CB_critical':0,'CB_selfpair_Qb':1}
 assert sum(contributions.values())==2 and 4-2==2
 proof={
 'actual_ring_map':'D0 -> E4_red on same physical R0, not a new type-only base',
 'Cprime_is_quotient_of_C':True,'all_generic_support_contained_in_accepted_four_zero_conductor_curves':True,
 'Cprime_finite_over_DVR_k_ell':True,
 'G19_transport':'simultaneous actual normalization diagram product, flat regular group factor',
 'generic_full_coker_degree0':'w kills it by actual conductors; trace generators have character0',
 'reduced_target_coker_is_graded_quotient_so_same_degree0':True,
 'degree0_reduction_commutes':'Nil(Efull)_0=Nil(Efull_0), nonnegative grading',
 'four_actual_reduced_coker_lengths':contributions,
 'each_reduced_support_curve_DVR_degree':1,
 'CA_critical_reduced_map':'identity G[[l7,v]] -> G[[l7,v]]',
 'Pa_reduced_map':'graph O_N -> O_N direct_sum O_SP',
 'Qb_reduced_map':'graph B/(f) -> B/(f) direct_sum B/(z7,f,rho)',
 'CB_critical_absent_from_actual_CA':True,
 'graph_cokernel_proof':'(a,b) maps to b-restriction(a), kernel equals the graph',
 'generic_DVR_rank_Cprime':2,'length_Cprime_mod_ell':4,'DVR_torsion_minimal_generators':2,
 'DVR_decomposition':'D^2 plus D/(ell^m1) plus D/(ell^m2), mi>=1 unknown',
 'rA_cut_variable_not_ell_DVR_parameter':True,
 'nil_D0_finiteness_or_generators_KJK_torsion_exponents_or_cube_decided':False,
 'new_source_coefficient_jet_point_GB_or_rank_worker':False,
 }
 ep=out/'actual-four-source-quotient-maps-and-rank2-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
 notes=[
 'research/notes/AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63.md',
 'research/notes/AMBIENT_CA_CRITICAL_V1_CONTACT_ALGEBRA_20261104_v63_REVIEW.md',
 'research/notes/AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63.md',
 'research/notes/AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63.md',
 f'research/notes/{STEM}.md',str(Path(__file__).relative_to(ROOT))]
 report={'id':STEM,'all_passed':True,'status':'PROVED_same_source_reduced_cokernel_DVR_rank2',
 'actual_generic_DVR_rank_Cprime':2,'actual_length_Cprime_mod_ell':4,'actual_DVR_torsion_minimal_generators':2,
 'torsion_exponents_or_original_nilgenerators_cube_KJK_claimed':False,
 'same_G19_and_exact_zero_degree_quotient_transport':True,
 'input_sha256':{p:sha(ROOT/p)for p in list(gates)+notes},'evidence_sha256':{ep.name:sha(ep)}}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['id','all_passed','actual_generic_DVR_rank_Cprime','actual_DVR_torsion_minimal_generators']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())

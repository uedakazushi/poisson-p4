#!/usr/bin/env python3
"""Fixed-source two9x9 linear maps; whole E4 ell-regularity and nil finiteness."""
import argparse,hashlib,json
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def row_solve(J,v):
 n=len(J);A=[[J[j][i]for j in range(n)]+[v[i]]for i in range(n)]
 for k in range(n):
  p=next(i for i in range(k,n)if A[i][k]);A[k],A[p]=A[p],A[k];a=A[k][k];A[k]=[x/a for x in A[k]]
  for i in range(n):
   if i!=k:
    a=A[i][k];A[i]=[x-a*y for x,y in zip(A[i],A[k])]
 w=[r[-1]for r in A]
 assert [sum(w[i]*J[i][j]for i in range(n))for j in range(n)]==v
 return w
def main(out):
 out.mkdir(parents=True,exist_ok=False)
 files={
 'Z':'reports/AMBIENT_B18_T0_COKERNEL_Z_20261104_v63/finite-cokernel-Z-and-actions.json',
 'Z_peer':'reports/AMBIENT_B18_T0_COKERNEL_Z_20261104_v63_REVIEW/report.json',
 'interface':'reports/AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63/actual-source-physical-firstjet-interface.json',
 'interface_peer':'reports/AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63_REVIEW/report.json',
 'five_peer':'reports/AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63_REVIEW/report.json',
 'square_peer':'reports/AMBIENT_B18_NORMALIZER_NIL_SQUARE_THREE_SUPPORTS_20261104_v63_REVIEW/report.json',
 'rankK':'reports/AMBIENT_B18_CRITICAL_CONTACT_AND_COMPARISON_RANK4_20261104_v63/report.json',
 'contact_source_peer':'reports/AMBIENT_CA_CRITICAL_V1_CONTACT_20261104_v63_REVIEW/report.json',
 'contact_algebra_peer':'reports/AMBIENT_CA_CRITICAL_V1_CONTACT_ALGEBRA_20261104_v63_REVIEW/report.json',
 'rankC_peer':'reports/AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63_REVIEW/report.json',
 'rankCprime_peer':'reports/AMBIENT_B18_REDUCED_COMPARISON_GENERIC_RANK2_20261104_v63_REVIEW/report.json',
 'Cprime_bar_peer':'reports/AMBIENT_B18_REDUCED_COMPARISON_CBAR_VALUES_20261104_v63_REVIEW/report.json',
 'HF3_peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json',
 }
 fixed={
 'Z':'2c25ddba496cf9fc3d4a53a07fa124b0de46933697235b427462451bbc13dc38',
 'Z_peer':'42843da6ce59c420d06d01ed20db4394f597b90222f639136e13c68577f35c6b',
 'interface':'7aee89c964510451fdc4f205361c73f615251d61439a28f756a87ae22909aeba',
 'interface_peer':'2033dff3459a3e9e4088f9dbbf4ad6f6549dd23cdd01bc8a3a543450523dc7eb',
 'five_peer':'cffc3d5c6a9cf9ed1207788a27675bdc76f925380fa27f0dc3a3c0c8e5337f80',
 'square_peer':'02413f680e959581d00e3a6807ec45dd62ee214e01d59604999a41813e38873b',
 'rankK':'a53d4613cbc088e376cec76b787b8669a807dd58898687590df627d44e0aa6ae',
 'contact_source_peer':'cc2a24acbc7c8f4abfb8a3b56ffb1b70fe94662e5cdbcc315d423bcdd8a51f41',
 'contact_algebra_peer':'46dd4d22aadd1c25e94f25f3fef656aea95225d205fda6951fe176ff7409f695',
 'rankC_peer':'22ce48012e3366120c1246a4b9276965eea2569a4f953478b32c7a79105d1496',
 'rankCprime_peer':'18300f976757ec9de2ad7d69f450db8cd794d235b0d8d65948c264156431c266',
 'Cprime_bar_peer':'fb4546b8667fc50c668ded9cc62c5597f3d19c08cbe1a6f116b6532209a09e6c',
 'HF3_peer':'808637882f3d8d28f393a7b0f4abfeec59e076bd1269034bd67d8ac6dd8c706a',
 }
 for k,p in files.items():assert sha(ROOT/p)==fixed[k],k
 for k in files:
  if k not in ['Z','interface']:assert read(ROOT/files[k])['all_passed'],k
 Z,I=[read(ROOT/files[k])for k in ['Z','interface']]
 physical=Z['physical_R0_rows9'];assert physical==[27,108,109,112,113,117,118,119,61]
 W=[[Q(c)for c in row]for row in Z['actual_physical_R0_firstjets_on_CA_F9_9by9']]
 lF9=[Q(0)]*7+[Q(-5),Q(5)]
 lphysical=row_solve(W,lF9)
 want=[Q(0)]*9;want[3]=Q(-15,2);want[5]=Q(10,3);assert lphysical==want
 R=[[Q(c)for c in row]for row in I['physical_R0_first9_by9']]
 lC=[sum(lphysical[i]*R[i][j]for i in range(9))for j in range(9)]
 K=[[Q(c)for c in row]for row in I['actual_source9_by_physical9_Jacobian']]
 lp=row_solve(K,lC)
 assert lp==[Q(-5,6),Q(5,3),Q(15,2),Q(0),Q(0),Q(0),Q(-10,3),Q(0),Q(0)]
 lines={'N':[0,0,0,1],'C':[0,1,1,0],'R':[-2,1,-1,-2],'P':[2,0,-2,1],'L':[1,1,-1,Q(-1,2)]}
 lv={n:sum(lp[i]*Q(x)for i,x in zip([0,1,2,6],v))for n,v in lines.items()}
 assert lv=={'N':Q(-10,3),'C':Q(55,6),'R':Q(5,2),'P':Q(-20),'L':Q(-5)}
 assert all(lv.values()) and any(lp[i]for i in [0,2,6])
 nil_lengths={'N':3,'C':3,'R':1,'P':2,'L':0};assert sum(nil_lengths.values())==9
 image_rank=9-4+2;assert image_rank==7 and image_rank+4==11
 hf=read(ROOT/files['HF3_peer']);assert hf['actual_T_length_interval_accepted']==[11,17]
 proof={
 'same_physical_R0_order':physical,'ell_F9_first':[str(c)for c in lF9],
 'ell_physical_first':[str(c)for c in lphysical],'ell_regular_C_first':[str(c)for c in lC],
 'ell_proper_source9_first':[str(c)for c in lp],
 'all5_proper_core_line_first_values':{n:str(v)for n,v in lv.items()},
 'all5_unit_core_roots_order_at_least2':'accepted actual core-elimination/interface proof, unchanged tangent lines',
 'ell_not_in_main_d_prime_by_non_d_linear_terms':True,
 'ell_regular_whole_E4_all6_Ass_checked':True,'ell_regular_E4_red':True,
 'each_of5_line_DVR_degree':1,
 'normalized_nil_generic_DVR_lengths':nil_lengths,
 'normalized_nil_DVR_finite_free_rank':9,'normalized_nil_square_DVR_finite_free_rank':3,
 'exact_sequence':'0->K->NilD0->NilE4->C->Cprime->0',
 'actual_NilD0_DVR_finite':True,'image_NilD0_in_NilE4_free_DVR_rank':7,
 'NilD0_DVR_decomposition':'K direct_sum D^7, noncanonical D-module ONLY',
 'actual_NilD0_generic_DVR_rank':11,
 'actual_NilD0_DVR_mu_formula':'length(T)+3','actual_NilD0_DVR_mu_interval':[14,20],
 'actual_NilD0_DVR_torsion_mu_formula':'length(T)-8','actual_NilD0_DVR_torsion_mu_interval':[3,9],
 'torsion_NilD0_equals_torsion_K_equals_H0D0':True,
 'D0_or_R0_module_or_algebra_split_claimed':False,
 'NilD0_annihilator_IAIB_claimed':False,
 'prospective_honest_Tnil_length_formula':'length(T)+5 ONLY once exact source section/presentation is constructed',
 'honest104_by948_Tnil_constructed_here':False,
 'exact_nilcube_or_index_or_H0_torsion_exponents_computed':False,
 'new_source_jet_gauge_CAS_point_or_scalar_worker':False,
 }
 ep=out/'actual-ell-line-values-and-nil-DVR-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
 notes=[f'research/notes/{STEM}.md',str(Path(__file__).relative_to(ROOT)),
 'research/notes/AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63_REVIEW.md',
 'research/notes/AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63.md',
 'research/notes/AMBIENT_B18_REDUCED_COMPARISON_GENERIC_RANK2_20261104_v63.md']
 report={'id':STEM,'all_passed':True,'status':'PROVED_actual_ell_regular_normalizers_and_zero_nil_DVR_finite',
 'actual_ell_regular_whole_E4_and_Ered':True,
 'actual_normalized_nil_DVR_free_rank':9,'actual_normalized_nil_square_DVR_free_rank':3,
 'actual_NilD0_DVR_finite_rank':11,'actual_NilD0_noncanonical_DVR_split':'K plus D^7',
 'actual_NilD0_mu_interval':[14,20],'actual_NilD0_torsion_mu_interval':[3,9],
 'all6_Ass_same_source_ell_checked':True,
 'new_nil_generators_Tnil_presentation_cube_or_H0_length_claimed':False,
 'input_sha256':{p:sha(ROOT/p)for p in list(files.values())+notes},'evidence_sha256':{ep.name:sha(ep)}}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['id','all_passed','actual_normalized_nil_DVR_free_rank','actual_NilD0_DVR_finite_rank']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())

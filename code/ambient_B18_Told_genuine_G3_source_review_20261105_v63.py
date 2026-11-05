#!/usr/bin/env python3
"""Independent same-base triangle and one pre-existing natural-column callback gate."""
import ast,hashlib,json,resource,signal,subprocess,sys,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_TOLD_GENUINE_G3_SOURCE_20261105_v63_REVIEW';PINS={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):PINS[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 rp='reports/AMBIENT_B18_TOLD_GENUINE_G3_SOURCE_20261105_v63/report.json';r=read(rp);assert r['all_passed']and not r['witness_evaluated']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p in r['input_sha256']:
  if p.endswith('/report.json'):assert read(p)['all_passed']
 S=read('reports/AMBIENT_B18_TOLD_GENUINE_G3_SOURCE_20261105_v63/actual-old-to-genuine-edge-one-cache-witness-and-diagram.json')
 N=read('reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63/actual-natural-Told-to-Tnil-source-map-and-origin.json')
 M=N['actual_origin_data']['actual_origin_map_columns5'];assert len(M)==5 and all(len(col)==5 and Q(col[0])==0 for col in M)
 assert M==S['all_five_natural_origin_columns5']and N['all32_native_full_polynomial_chain_identities_checked']
 assert N['source_C_ring_and_frame']=='same actual Cfree4 over R0, exact ell0 CA-F9; no abstract type substitution'
 assert N['old_minimal_generator_free56_rows']==[22,23,29,30,55]
 edge=read('reports/AMBIENT_B18_NBAR14_KERNEL_CACHED_VALUES_20261105_v63/actual-five-edge-functionals-and-genuine-Nbar14-cycles.json')
 assert edge['genuine_Nbar_dimension']==14
 vals=list(map(Q,edge['five_actual_cached_functional_values_a_b1_b2_b3_b4']));xi=Q(edge['xi_equals_minus11a_by16']);assert vals[0]and xi==Q(-11,16)*vals[0]
 W=S['one_cached_witness'];assert W['natural_origin_column5']==M[4]and W['old_generator_index']==4 and W['old_actual_free56_label']==55
 assert list(map(Q,W['genuine_edge_b1_b2_b3_b4']))==vals[1:]and Q(W['genuine_edge_rh0_coefficient_xi'])==xi
 old=read(W['natural_head0_r_first_callback_path']);J=Q(old[W['natural_head0_r_first_callback_key']]);assert J==Q(-43373,15129)
 ownJ=read('reports/AMBIENT_B18_TOLD_TNIL_HEAD4_Z5_20261105_v63_REVIEW/phase2-independent-physical-source2-map-and-complete-HF-proof.json');assert Q(ownJ['actual_Jr'])==J
 # This exact old first coefficient, including oldn1 and18/81 Schurs, is not
 # the head4*z5 quadratic class nor its R35 correction.
 finite=read('reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63/report.json');fproof=read('reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63/actual-ell-line-values-and-nil-DVR-proof.json')
 assert finite['actual_ell_regular_whole_E4_and_Ered']and finite['actual_normalized_nil_DVR_free_rank']==9 and finite['actual_NilD0_noncanonical_DVR_split']=='K plus D^7'
 assert fproof['image_NilD0_in_NilE4_free_DVR_rank']==7 and not fproof['D0_or_R0_module_or_algebra_split_claimed']
 oldL=read('reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/report.json');assert oldL['actual_T_length']==11 and oldL['actual_Tnil_length']==16
 G=read('reports/AMBIENT_B18_CPRIME_GENUINE_G3_VALUES_20261105_v63/actual-Cprime9-actions-and-genuine-G3-target.json');assert G['actual_genuine_G3_action_scalars9']==['0']*7+['-11/16','0']
 assert G['genuine_target_dimension']==2
 custody=read('reports/AMBIENT_B18_TOLD_GENUINE_G3_VALUE_CUSTODY_20261105_v63/report.json')
 for p,h in custody['input_sha256'].items():assert sha(ROOT/p)==h
 code='code/ambient_B18_Told_genuine_G3_values_20261105_v63.py';PINS[code]=sha(ROOT/code);assert sha(ROOT/code)==custody['prepared_code_sha256']=='4ca3a74aa578032aad055483a0df97c1509c28970c8120c0220f7e9a0829a631'
 tree=ast.parse((ROOT/code).read_text());assert not any(isinstance(n,ast.ImportFrom)and n.module not in ['fractions','pathlib']for n in ast.walk(tree))
 assert subprocess.run([sys.executable,str(ROOT/code),'--help'],capture_output=True,timeout=5).returncode==0
 assert not (ROOT/'reports/AMBIENT_B18_TOLD_GENUINE_G3_VALUES_20261105_v63').exists()
 data={'all5_natural_head0_origins_zero':True,'oldhead4_source_free56_label55':True,'actual_cached_Jr_not_quadratic_product':str(J),'prior_independent_source_derivative_reused':True,
 'rank_upper1_source_proof':'For arbitrary minimal coordinate v, epsilon(v)=a*v0(0)*K+(beta*v0(0)+sum bj*vj(0)+xi*dr(v0))*Q1. All5 natural origin head0 entries0. Hence their images lie socle kQ1; S-linearity bounds the entire oldimage dimension by1. Beta is uncomputed but multiplied by0; allz anddegree>=2 vanish inG3.',
 'genuine_source_diagram_proof':'The two accepted comparison cones use the same actual D0 free9/640 source and E4→Ered quotient, with actual old56→new44 transport and all32 full polynomial boundary identities. Short exact complexes0→NE[0]→[D0→E4]→[D0→Ered]→0, NE Dfree9, E4/Ered ellregular, give0→Told→Tnil→NEbar→Cbar→Cprimebar→0 on derived ellfibers. Tor1(D0) occurs in degree2 and is not discarded. Equivalently apply snake to the complete sameD0/ell maps and exact0NEbar→E4bar→Eredbar.',
 'compatibility_with_actual_nil_image':'ActualNbar→Tnil→NEbar is reduction of N→NE, by the same comparison chain and genuine kernel map. Its kernel equals Nbar∩imageTold. No abstract isomorphic target replacement is used.',
 'conditional_positive_dimensions':{'epsilon_old_rank':1,'intersection_Told_Nbar':10,'image_Nbar_in_NEbar':4,'I_Dfree_rank':7,'NE_Dfree_rank':9,'L_Dgeneric_rank':2,'L_torsion_mu':3},
 'conditional_Tor_argument':'Let I=im(N→NE), L=NE/I. Nbar surjects Ibar. Positiveoldedge rank1 givesintersection11-1=10, hence Nbarimage14-10=4. Tensor0D7→D9→L0 gives0L[ell]→k7→k9→Lbar0; therefore dimL[ell]=3, dimLbar=5 and rankDL2. These are torsion generator counts, not lengths or exponents.',
 'Q_not_Lbar':'Q=Tnil/Told is canonically a length5 submodule im(theta) ofNEbar; Lbar is a length5 quotient ofNEbar. Equal dimensions imply no identification, no S-module splitting, no extrapolation from the r-onlyJordan311 statement.',
 'five_term_witness_NOT_paired':True,'future_fixed_input117_static15000_hard30000':True,'source_bit163_under512_under2048':True,'prepared_no_new_derivatives_source_gauge_jet':True}
 ep=out/'independent-same-base-triangle-and-existing-Jr-witness-source.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 result={'id':STEM,'all_passed':True,'status':'SOURCE_PASS_actual_old_edge_one_cache_witness_and_normalizer_diagram','input_sha256':PINS,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'witness_paired':False,'prepared_code_sha256':sha(ROOT/code),'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));signal.alarm(0)
if __name__=='__main__':main()

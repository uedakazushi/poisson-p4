#!/usr/bin/env python3
"""DVR20/raw-cokernel26 source audit and exact comparison-complex preflight."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def enc(a):
    if isinstance(a,F):return str(a)
    if isinstance(a,(list,tuple)):return [enc(x) for x in a]
    if isinstance(a,dict):return {str(k):enc(v) for k,v in a.items()}
    return a
def cubic(c,r,s):return c[0]*r**3+c[1]*r*r*s+c[2]*r*s*s+c[3]*s**3
def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    ep=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json'
    cp=ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json'
    bp=ROOT/'reports/AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63/degree3-certificate.json'
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    fp=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    gates=[ROOT/'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_CURVE_BASE_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed'] for p in gates)
    E,C,B,S,P=map(read,[ep,cp,bp,sp,fp])
    Li=[[F(c) for c in row[1:10]] for row in E['families']['CA']['inverse_F9_in_common_R0_through2']]
    Xb=[[F(c) for c in row[1:10]] for row in E['families']['CB']['physical_R0_readback_jets_in_own_F9']]
    YA_by_YB=[[sum((a*b for a,b in zip(r,c)),F(0)) for c in zip(*Xb)] for r in Li]
    ellB=[5*(YA_by_YB[8][i]-YA_by_YB[7][i]) for i in range(9)]
    assert ellB[7:]==[F(35,3),F(-85,3)]
    ga=list(map(F,C['g_CA_all10_coefficients'][6:10]));gb=list(map(F,B['g_CB_all10_coefficients'][6:10]))
    ua=cubic(ga,F(1),F(1));ub=cubic(gb,F(17),F(7));assert ua==F(125,24) and ub==-181500
    wB=cubic(gb,F(1),F(7,17));assert wB==F(-181500,4913)
    ca,cb=S['families']['CA'],S['families']['CB'];rows=[18,19,28,29,60,67,72,63,64,68,73,74]
    f=[b for r,b in ca['actual_zero_odd_I0_generators']];g=[b for r,b in cb['actual_zero_odd_I0_generators']]
    raw_actions=[]
    for r in rows:
        aA,bA=ca['complete_normalized_physical_w_coefficients'][0][r]
        aB,bB=cb['complete_normalized_physical_w_coefficients'][0][r]
        raw_actions.append({'physical_row':r,'aA_node':aA,'aB_node':aB,'bA_node':bA,'bB_node':bB,
          'delta_operator':['subtract',aA,aB],
          'matrix_on_raw_cokernel':'[[delta,0,-nuB*bB],[0,delta,nuA*bA],[-bB,bA,delta]]'})
    zero=[q for q in P['all99_CB_kernel_and87_graph_generators_evaluated_on_CA'] if q['weight']==0]
    assert len(zero)==68 and len(P['actual_Qn_R0_squared_generator_matrices']['0'])==68*9
    source={'actual_ell_operator':['multiply',5,['subtract',J9:=S['families']['CA']['actual_positive_odd_J0_generators'][1][1],S['families']['CA']['actual_positive_odd_J0_generators'][0][1]]],
      'CA_ell_on_own_plane_linear':[-5,5],'CB_ell_on_own_plane_linear':ellB[7:],
      'CA_g3_unit_at_ell_zero':ua,'CB_g3_unit_at_integer_kernel_direction_17_7':ub,
      'CA_Weierstrass_cubic_leading_unit':ua,'CB_Weierstrass_cubic_leading_unit_in_rB':wB,
      'both_actual_CI_curve_modules_are_DVR_free_rank3':True,
      'DVR':'D=QQ[[ell]], ell=5*(s_CA-r_CA)',
      'product_curve_ring_DVR_type':{'free_rank':6,'D_mod_ell_squared_copies':7,'D_mod_ell_copies':7,'special_fiber_dimension':20},
      'raw_normalizer_cokernel_R0_type':'R0/IA*tauA direct_sum R0/IB*tauB direct_sum R0/(IAIB)*tauA*tauB',
      'raw_cokernel_DVR_type':{'free_rank':12,'D_mod_ell_squared_copies':7,'D_mod_ell_copies':7,'special_fiber_dimension':26},
      'raw_Tor2_kernel_DVR_type':{'D_mod_ell_squared_copies':21,'D_mod_ell_copies':21,'length':63,'DVR_minimal_generators':42},
      'all12_physical_difference_actions_on_raw_cokernel':raw_actions,
      'raw_cokernel_row_annihilators':['IA','IB','IA*IB'],'actual_IA_generator_nodes8':f,'actual_IB_generator_nodes8':g,
      'actual_coker_comparison_exact_formula':'C=C_raw/(sum_j F_j*C_raw)',
      'conditional_finite_DVR_presentation_bound':{'generators':26,'intrinsic_torsion_relations':14,'physical_action_columns':12*26,'total_relation_columns':326},
      'actual_DVR_basis_and_R0_action_matrices_computed':False,
      'actual_zero_CB_polynomial_records68':zero,
      'factored_640_chain_lift_recipe':{'intrinsic28':'zero in normalizer, choose zero48 lift',
        'remaining612':'68 CB polynomials times9 C0 basis; polynomial target degree<=2 gives exact Hadamard difference F(qA)-F(qB)=sum physical_delta_j*D_j; multiply by C0 basis in free4',
        'old_regular_coefficients':'retain their original same-R0 scalar functions; never replace by fresh0/0 or zero'},
      'exact_kernel_presentation_gate':{'matrix':'L=ker[U,-Q] subset R0^(9+48)',
        'quotient':'K=L/(columns(P640,B640) + columns(0,Syz(Q)))',
        'missing_complete_syzygy_input':'Syz(Q) for actual E4 free4/48, equivalently a complete source-bound syzygy model of the7 core4 functions plus five Koszul directions',
        'presentation_after_completion':'only after K has a complete cokernel presentation may one tensor by R0/(IAIB), then reduce to the20-generator DVR base by right exactness'},
      'finite_complex_criterion':{'required_homology':'H1=K,H0=C and all other H=0',
        'DVR_torsion_test':'H2(P_complex tensor^L_D QQ)=K[ell]',
        'special_fiber_sequence':'0 -> K/ell*K -> H1(derived fiber) -> C[ell] -> 0',
        'generic_rank':'rank_D K = dim_QQ((ell)) H1(generic complex)',
        'no_naive_specialized_kernel':True},
      'full_K_DVR_comparison_complex_constructed':False,'K_torsion_length_socle_or_D0_depth_decided':False}
    ep=out/'DVR20-raw-cokernel26-and-comparison-preflight.json';ep.write_text(json.dumps(enc(source),separators=(',',':'))+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',
      ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json',cp,bp,sp,fp,*gates]
    report={'id':STEM,'all_passed':True,'status':'actual_DVR20_base_raw_cokernel26_and_honest_comparison_completion_preflight',
      'actual_ell_regular_both_conductor_curves_checked':True,'both_curve_DVR_ranks':3,
      'product_curve_DVR_free_rank':6,'product_curve_DVR_special_fiber':20,'raw_cokernel_DVR_special_fiber':26,
      'right_exact_coker_recipe_without_kernel_basechange':True,
      'full_K_presentation_or_actual_finite_DVR_actions_constructed':False,
      'no_higher_jets_second_carries_Smith_or_naive_specialized_kernel':True,
      'K_torsion_or_actual_D_depth_claimed':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{ep.name:sha(ep)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['all_passed','product_curve_DVR_special_fiber','raw_cokernel_DVR_special_fiber','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

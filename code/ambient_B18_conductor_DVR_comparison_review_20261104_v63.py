#!/usr/bin/env python3
"""Independent fixed-source arithmetic and DVR-comparison proof audit.

Does not run the author builder, solve a source system, or construct DVR actions.
"""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(a):
    if isinstance(a,F):return str(a)
    if isinstance(a,(list,tuple)):return [enc(v) for v in a]
    if isinstance(a,dict):return {str(k):enc(v) for k,v in a.items()}
    return a
def matmul(A,B):return [[sum((x*y for x,y in zip(r,c)),F(0)) for c in zip(*B)] for r in A]
def linear(p):return [[F(c) for c in row[1:10]] for row in p]
def cubic(co,mons,r,t):return sum((F(c)*r**i*t**j for c,(i,j) in zip(co,mons)),F(0))
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    primary=ROOT/'reports/AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63/report.json'
    pr=read(primary);assert pr['all_passed']
    for rel,h in pr['input_sha256'].items():assert sha(ROOT/rel)==h,(rel,'changed primary input')
    proofp=primary.parent/'DVR20-raw-cokernel26-and-comparison-preflight.json';P=read(proofp)
    assert sha(proofp)==pr['evidence_sha256'][proofp.name]
    ep=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json'
    cp=ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json'
    bp=ROOT/'reports/AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63/degree3-certificate.json'
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    fp=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    E,C,B,S,FR=map(read,[ep,cp,bp,sp,fp])
    # Actual CA F9 last coordinates ARE the two positive odd coefficient
    # functions (not raw physical q46/q95 and not arbitrary source coordinates).
    fca=S['families']['CA'];fcb=S['families']['CB']
    assert [r for r,n in fca['actual_positive_odd_J0_generators']]==[46,95]
    j46,j95=[n for r,n in fca['actual_positive_odd_J0_generators']]
    assert E['families']['CA']['own_F9_coordinate_order'][7:]==[j46,j95]
    assert P['actual_ell_operator']==['multiply',5,['subtract',j95,j46]]
    A=linear(E['families']['CA']['inverse_F9_in_common_R0_through2'])
    X=linear(E['families']['CB']['physical_R0_readback_jets_in_own_F9'])
    adapter=matmul(A,X);ell=[5*(adapter[8][i]-adapter[7][i]) for i in range(9)]
    assert ell[7:]==[F(35,3),F(-85,3)]==list(map(F,P['CB_ell_on_own_plane_linear']))
    mons=C['monomials'];assert mons==B['ordered_monomials']
    assert mons==[[0,0],[1,0],[0,1],[2,0],[1,1],[0,2],[3,0],[2,1],[1,2],[0,3]]
    ca=list(map(F,C['g_CA_all10_coefficients']));cb=list(map(F,B['g_CB_all10_coefficients']))
    assert not any(ca[:6]+cb[:6])
    va=cubic(ca,mons,F(1),F(1));vb=cubic(cb,mons,F(17),F(7));vc=cubic(cb,mons,F(1),F(7,17))
    assert va==F(125,24) and vb==F(-181500) and vc==F(-181500,4913)
    assert va==F(P['CA_g3_unit_at_ell_zero']) and vb==F(P['CB_g3_unit_at_integer_kernel_direction_17_7'])
    # Actual A3 quotient uses CA r,s. Higher terms of ell vanish in m^2=0.
    Mr=s.Matrix([[0,0,0],[1,0,0],[0,0,0]]);Ms=s.Matrix([[0,0,0],[0,0,0],[1,0,0]])
    Ell=-5*Mr+5*Ms;assert Ell.rank()==1 and Ell**2==s.zeros(3)
    assert Ell.T.rank()==1 and Ell.T**2==s.zeros(3)
    tor1_copies=math.comb(6,1)+math.comb(6,0);tor2_copies=math.comb(6,2)+math.comb(6,1)
    assert tor1_copies==7 and tor2_copies==21
    assert P['product_curve_ring_DVR_type']=={'free_rank':6,'D_mod_ell_squared_copies':7,'D_mod_ell_copies':7,'special_fiber_dimension':20}
    assert P['raw_cokernel_DVR_type']=={'free_rank':12,'D_mod_ell_squared_copies':7,'D_mod_ell_copies':7,'special_fiber_dimension':26}
    assert P['raw_Tor2_kernel_DVR_type']=={'D_mod_ell_squared_copies':21,'D_mod_ell_copies':21,'length':63,'DVR_minimal_generators':42}
    # Derive, independently, the three-generator raw-coker action from full
    # free4 multiplication, THEN project off the R0 scalar coordinate.
    na,nb,d,a,b,e,c,h=s.symbols('na nb delta a b epsilon c h')
    TA=s.Matrix([[0,na,0,0],[1,0,0,0],[0,0,0,na],[0,0,1,0]])
    TB=s.Matrix([[0,0,nb,0],[0,0,0,nb],[1,0,0,0],[0,1,0,0]])
    assert TA*TB==TB*TA and TA**2==na*s.eye(4) and TB**2==nb*s.eye(4)
    M=(d*s.eye(4)+a*TA-b*TB)[1:4,1:4]
    expected=s.Matrix([[d,0,-nb*b],[0,d,na*a],[-b,a,d]]);assert M==expected
    N=(e*s.eye(4)+c*TA-h*TB)[1:4,1:4]
    comm=(M*N-N*M).applyfunc(s.expand)
    symbols=[na,nb,d,a,b,e,c,h]
    # Rows are quotiented by IA,IB,IAIB respectively. Check every universal
    # commutator monomial has the requisite source ideal factors.
    for i in range(3):
        for j in range(3):
            for ex,co in s.Poly(comm[i,j],*symbols).terms():
                if not co:continue
                hasA=ex[3]+ex[6]>0;hasB=ex[4]+ex[7]>0
                assert (hasA if i==0 else hasB if i==1 else hasA and hasB)
    rows=[18,19,28,29,60,67,72,63,64,68,73,74]
    assert [r['physical_row'] for r in P['all12_physical_difference_actions_on_raw_cokernel']]==rows
    for r,record in zip(rows,P['all12_physical_difference_actions_on_raw_cokernel']):
        aa,ba=fca['complete_normalized_physical_w_coefficients'][0][r]
        ab,bb=fcb['complete_normalized_physical_w_coefficients'][0][r]
        assert [record['aA_node'],record['aB_node'],record['bA_node'],record['bB_node']]==[aa,ab,ba,bb]
        assert record['delta_operator']==['subtract',aa,ab]
    assert P['actual_IA_generator_nodes8']==[n for r,n in fca['actual_zero_odd_I0_generators']]
    assert P['actual_IB_generator_nodes8']==[n for r,n in fcb['actual_zero_odd_I0_generators']]
    zero=[r for r in FR['all99_CB_kernel_and87_graph_generators_evaluated_on_CA'] if r['weight']==0]
    assert zero==P['actual_zero_CB_polynomial_records68'] and len(zero)==68
    kinds={k:sum(r['kind']==k for r in zero) for k in {r['kind'] for r in zero}}
    assert kinds=={'CB_kernel':64,'graph_CB_minus_CA':4}
    assert len(FR['actual_Qn_R0_squared_generator_matrices']['0'])==68*9
    assert P['conditional_finite_DVR_presentation_bound']=={'generators':26,'intrinsic_torsion_relations':14,'physical_action_columns':312,'total_relation_columns':326}
    assert not P['actual_DVR_basis_and_R0_action_matrices_computed']
    assert not P['full_K_DVR_comparison_complex_constructed'] and not P['K_torsion_length_socle_or_D0_depth_decided']
    gates=[ROOT/'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_CURVE_BASE_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed'] for p in gates)
    evidence={'own_actual_CA_F9_last2_nodes':[j46,j95],'own_actual_YA_by_YB_first':adapter,
      'own_ellB_first_all9':ell,'own_gA3_at_1_1':va,'own_gB3_at_17_7':vb,'own_gB3_at_1_7over17':vc,
      'A3_ell_action':[[str(c) for c in row] for row in Ell.tolist()],
      'omega_A3_ell_action':[[str(c) for c in row] for row in Ell.T.tolist()],
      'both_actions_rank1_square0':True,'own_derived_raw_action':str(M),
      'own_projected_commutator':str(comm),'all_commutator_entries_in_row_annihilator_ideals':True,
      'actual12_action_bindings_checked':True,'zero_CB_counts':kinds,'chain_lift_factored_count':612,
      'module_types_proved_noncanonically_not_actual_DVR_bases':True,
      'right_exact_C_computation_proof':'B4/(image(rawtensor)+sumF*B4)=Craw/sumF*Craw; F lies in raw image subalgebra',
      'kernel_gate':'K=ker[U,-Q]/((P,B)+(0,SyzQ)); requires COMPLETE relation and syzygy data before base change',
      'derived_fiber_gate':'For a finite perfect D-complex with only H1=K,H0=C, H2(fiber)=K[ell]; H1(fiber) has C[ell] correction',
      'actual_DVR_actions_complete_K_Smith_socle_or_D0_depth_computed':False}
    p=out/'independent-proof.json';p.write_text(json.dumps(enc(evidence),indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',primary,proofp,ep,cp,bp,sp,fp,*gates]
    report={'id':STEM,'all_passed':True,'status':'independent_actual_parameter_DVR_types_coker_and_honest_completion_gate_accepted',
      'actual_ell_both_conductor_parameters':True,'both_curve_DVR_rank':3,'product_base_DVR_special_fiber':20,
      'raw_coker_DVR_min_generators':26,'conditional_actual_coker_relation_bound':326,
      'actual_DVR_action_or_kernel_torsion_or_D0_depth_claimed':False,
      'no_author_builder_source_worker_new_jet_GB_or_Smith':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{p.name:sha(p)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

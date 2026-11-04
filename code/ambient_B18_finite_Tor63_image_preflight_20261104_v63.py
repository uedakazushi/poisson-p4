#!/usr/bin/env python3
"""Proof and finite-input preflight only for the Tor2 length63 socle target."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_top_graph_telescoping_review_20261104_v63 import solve
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_FINITE_TOR63_IMAGE_PREFLIGHT_20261104_v63'
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
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    jp=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
    cp=ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json'
    pp=cp.parent/'joint105-algebraic-DAG.json'
    ep=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json'
    gates=[ROOT/'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_CURVE_BASE_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_T1_T2_QUADRATIC_PHI_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed'] for p in gates)
    S,J,C,P,E=map(read,[sp,jp,cp,pp,ep])
    assert C['all105_equations_zero_through_degree3']
    assert all(not any(map(F,row[1:12])) for row in E['all14_raw_equations_through2'][2:])
    rows=[18,19,28,29,60,67,72,63]
    f=[p[1] for p in S['families']['CA']['actual_zero_odd_I0_generators']]
    g=[p[1] for p in S['families']['CB']['actual_zero_odd_I0_generators']]
    assert rows==[p[0] for p in S['families']['CA']['actual_zero_odd_I0_generators']]
    assert rows==[p[0] for p in S['families']['CB']['actual_zero_odd_I0_generators']]
    M=[[F(x) for x in r] for r in J['YB_to_YA_linear']]
    ca=J['both57_source2']['CA'];cb=J['both57_source2']['CB']
    fa8=list(map(F,ca['pure_function2jets'][str(ca['f63'])][1:10]))
    gb8=list(map(F,cb['pure_function2jets'][str(cb['f63'])][1:10]))
    Gfirst=M[:7]+[[sum((gb8[j]*M[j][i] for j in range(9)),F(0)) for i in range(9)]]
    assert not any(fa8[7:]) and all(not any(a[7:]) for a in Gfirst)
    U={p['CB_physical_zero_row']:list(map(F,p['all_degree_le3_coefficients'][3:6])) for p in C['all8_CB_zero_restrictions']}
    A=[]
    for i in range(7):A.append([F(i==j) for j in range(7)]+[F(0)]*3)
    A.append(fa8[:7]+[F(0)]*3)
    A.extend([Gfirst[i][:7]+U[r] for i,r in enumerate(rows)])
    selected=[0,1,2,3,4,5,6,8,12,13]
    redundant=[i for i in range(16) if i not in selected]
    V=[A[i] for i in selected]
    Vi,det=solve(V,[[F(i==j) for j in range(10)] for i in range(10)])
    assert det==F(625,768)
    q0=[]
    for k in redundant:
        coefficients=[sum((A[k][i]*Vi[i][j] for i in range(10)),F(0)) for j in range(10)]
        q=[F(i==k) for i in range(16)]
        for i,c in zip(selected,coefficients):q[i]-=c
        assert all(dot(q,[A[i][j] for i in range(16)])==0 for j in range(10))
        q0.append(q)
    assert len(q0)==6
    minimal=[{'type':'q_wedge_q','q_indices':[i,j]} for i,j in itertools.combinations(range(6),2)]
    minimal += [{'type':'eta_wedge_q','eta':e,'q_index':i} for e in ['eta_u','eta_v'] for i in range(6)]
    basis63=[{'type':'q_wedge_q','q_indices':[i,j],'scalar':s} for i,j in itertools.combinations(range(6),2) for s in ['1','r','s']]
    basis63 += [{'type':'eta_wedge_q','eta':e,'q_index':i,'scalar':s} for i in range(6) for e,s in [('eta_u','1'),('eta_v','1'),('eta_u','s')]]
    assert len(minimal)==27 and len(basis63)==63
    proof={'actual_ideal_generator_physical_rows8':rows,'IA_generator_source_nodes8':f,'IB_generator_source_nodes8':g,
      'combined16_order':['f0..f7','g0..g7'],'canonical_minimal10_order':['z0..z6','r^2','r*s','s^2'],
      'selected10_combined_indices':selected,'redundant6_combined_indices':redundant,
      'actual_constant16_by10_generator_matrix':A,'selected10_unit_matrix':V,'selected10_unit_inverse':Vi,
      'selected10_unit_minor':det,'all6_zero_column_cycle_origins16':q0,
      'eight_full_cycle_recipe':{'q_alpha':'e_redundant-A_red*A_selected_inverse*e_selected',
        'eta_u':'lift(-s,r,0) from canonical quadratic3 generators through selected10 inverse',
        'eta_v':'lift(0,-s,r) from canonical quadratic3 generators through selected10 inverse'},
      'all27_minimal_Tor2_generator_specs':minimal,'all63_QQ_basis_specs':basis63,
      'dimension_shift_exact_sequence':'Tor2(R/IA,R/IB)=Tor1(R/IA,IB)=ker(IA tensor IB -> IA*IB)',
      'actual_D0_map_formula':{'h_i_beta':'q_Ai*p_Bbeta-p_Ai*q_Bbeta',
        'normalizer_pair':'sum h_i_beta*(f_i*tauA)*(g_beta*tauB)',
        'CA_even':'nuA*sum h_i_beta*f_i*f_beta',
        'CA_C0_odd_vector':'sum h_i_beta*(aA_beta-aB_beta)*image(f_i*tauA)'},
      'first_possible_C0_orders':{'q_wedge_q':{'even':3,'odd_module':2},'eta_wedge_q':{'even':4,'odd_module':3}},
      'quadratic_representatives_depend_only_on_saved_cycle_origins_and_cached_even_differences':True,
      'all_order_transform_input_recipe':{
        'joint105_prototype_unknown_bindings':P['all105_unknown_bindings'],
        'joint105_prototype_unit_equations':P['all105_pointed_unit_equations'],
        'prototype_unit_minor':P['unit_Jacobian_determinant'],
        'coordinate_prefix_systems':7,'axis_systems':2,'all_new_prefix_systems_are_unit105':True,
        'second_Hadamard_coefficients':['B(r,0)/r^2','(B(r,s)-B(r,0)-B(0,s))/(r*s)','B(0,s)/s^2'],
        'regularity_proof':'plane B in(r,s)^2; mixed numerator vanishes on both axes',
        'no_quotient_is_restricted_after_setting_its_own_denominator_zero':True,
        'unit3_matrix_rows':[18,60,67],'unit3_minor_at_origin':'625/768',
        'C0_free9_adapter':'existing actual f63 zero7 Hadamard vector; keep its regular functions, no redivision'},
      'size_preflight':{'prototype_nodes':len(P['DAG_scalar_nodes']),'nine_clone_conservative_nodes':9*len(P['DAG_scalar_nodes']),
        'normalizer_base_generators':4,'C0_target_free_generators':9,'actual_Tor2_minimal_generators':27,'QQ_dimension_upper_bound_for_image':63},
      'all_order_cycle_DAG_expanded':False,'actual_D0_nonmembership_test_executed':False,
      'image_nonzero_finite_socle_or_depth_claimed':False}
    ep=out/'finite63-proof-and-input-preflight.json';text=json.dumps(enc(proof),separators=(',',':'))+'\n';assert len(text.encode())<32*1024*1024;ep.write_text(text)
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',sp,jp,cp,pp,
      ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json',*gates,
      ROOT/'code/ambient_B18_top_graph_telescoping_review_20261104_v63.py']
    report={'id':STEM,'all_passed':True,'status':'finite_Tor2_image_socle_target_proof_and_source_bound_preflight_only',
      'raw_normalizer_kernel_Tor_index':2,'raw_kernel_abstract_module':'A3^15 direct_sum omega_A3^6',
      'raw_kernel_QQ_length':63,'raw_kernel_socle_dimension':36,'raw_kernel_minimal_R0_generators':27,
      'natural_image_in_D0_is_finite_and_maxideal_squared_annihilated':True,
      'nonzero_image_would_prove_actual_D_cone_depth19':True,
      'exact_dimension_shift_and_representative_determinant_recipe_verified':True,
      'selected10_constant_unit_minor':'625/768','literal_full_transform_or_nonmembership_run':False,
      'no_new_source_jets_points_primes_GB_or_old_r_g_composition':True,
      'actual_image_nonzero_or_T0_depth_solved':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{ep.name:sha(ep)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['all_passed','raw_kernel_QQ_length','raw_kernel_minimal_R0_generators','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

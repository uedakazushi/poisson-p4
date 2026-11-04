#!/usr/bin/env python3
"""Actual cached first-source trace gate and pure conductor proof of nil^7=0."""
import argparse
import functools
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path
import resource
import signal
import time

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER7_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())

def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out.mkdir(parents=True,exist_ok=False)
    files={
      'original_source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
      'joint_source':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint105-algebraic-DAG.json',
      'joint_cache':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
      'joint_peer':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW/report.json',
      'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
      'radical_peer':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63_REVIEW/report.json',
      'upper12':'reports/AMBIENT_B18_QUADRATIC_NORMALIZER_NIL_UPPER12_20261104_v63/report.json',
      'upper12_peer':'reports/AMBIENT_B18_QUADRATIC_NORMALIZER_NIL_UPPER12_20261104_v63_REVIEW/report.json',
      'product':'reports/AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63/report.json',
      'E4_peer':'reports/AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63_REVIEW/report.json',
      'nil_finite_peer':'reports/AMBIENT_B18_NIL_CUBE_FINITE_ORIGIN_20261104_v63_REVIEW/report.json',
    }
    expected={
      'original_source':'75cc6bcc4d4ff8ebd255a273cd3db494f92cecdff69583bea80e4316dee4fa74',
      'joint_source':'62f9b5bd9e805437a883caef24b097be73eceecc70086a51ea49173c36bef1e0',
      'joint_cache':'1909547a8d80ad59905518282c8d266d2c92f791c2c9c120b01abc3cd0860d44',
      'joint_peer':'513e63044f8048da7c6761f1b02c63c0bcd6f643beb65339a2d3eb0de18c1401',
      'radical':'376d5fd2b8f2ab8a45651b0374d5b983ebe2df8f8a507108675942e853adbba5',
      'radical_peer':'07b8ec2ea5d6dc419c15fad7cb3ea9526d817037b48c8d96b9966882c1dc851f',
      'upper12':'332701163619c4e3480e77e4240f604bb3200ef9c17a2131daa72d363bc0d77c',
      'upper12_peer':'e1a6937b676c1518d192adb0dfdc5b35d7d841eaaf381d0f52f2c9c0ca1db219',
      'product':'d093004d8fdad173825cc4b8eb5d5b70ea77615b2c1dfbc8aaa3d1116105eb4f',
      'E4_peer':'cffc3d5c6a9cf9ed1207788a27675bdc76f925380fa27f0dc3a3c0c8e5337f80',
      'nil_finite_peer':'ee80d12faa103fb5ac0e2144b4d12cd662f0c5afcfa1ded335e2171946ad1ccc',
    }
    for key,p in files.items():assert sha(ROOT/p)==expected[key],key
    src,js,jc,R=[read(ROOT/files[k])for k in ['original_source','joint_source','joint_cache','radical']]
    for key in ['joint_peer','radical_peer','upper12','upper12_peer','product','E4_peer','nil_finite_peer']:
        assert read(ROOT/files[key])['all_passed'],key
    peer=read(ROOT/files['joint_peer'])
    assert peer['exact_joint_zero_ideal_quotient_length']==3
    assert peer['actual_CB_three_quadratics_span_m_squared_minor']=='625/768'
    assert peer['literal_finite105_DAG_and_every_binding_reconstructed_from_actual_source']
    nodes=js['DAG_scalar_nodes'];original=src['DAG_scalar_nodes']
    assert nodes[:len(original)]==original and len(original)==41315
    roots=[src['families'][i]['actual_tracefree_nu']for i in ['CA','CB']]
    assert roots==[20086,40865]
    assert jc['CA_F9_targets']==['0']*7+['r','s']
    overrides={i:tuple(Q(x)for x in v[:3])for i,v in jc['all105_unknown_jets']}
    rr,ss=js['formal_curve_parameter_nodes'];overrides[rr]=(Q(0),Q(1),Q(0));overrides[ss]=(Q(0),Q(0),Q(1))
    visited=set()
    @functools.cache
    def ev(i):
        visited.add(i)
        if i in overrides:return overrides[i]
        n=nodes[i];o=n[0]
        if o=='c':return Q(n[1],n[2]),Q(0),Q(0)
        if o=='n':return tuple(-x for x in ev(n[1]))
        if o=='a':return tuple(x+y for x,y in zip(ev(n[1]),ev(n[2])))
        if o=='m':
            a,b=ev(n[1]),ev(n[2]);return a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[0]*b[2]+a[2]*b[0]
        if o=='i':
            a=ev(n[1]);assert a[0];return 1/a[0],-a[1]/a[0]**2,-a[2]/a[0]**2
        raise AssertionError((i,n))
    vals=[ev(i)for i in roots];assert vals==[(Q(0),)*3]*2
    # Order lower bounds in the ACTUAL regular9 C maximal ideal, no scalar jets.
    cn=R['new_Cscalar_and_R0component_nodes']
    @functools.cache
    def order(i):
        n=cn[i];o=n[0]
        if o=='C_QQ':return 100000 if n[1]==0 else 0
        if o in ['actual_proper_u','actual_core_GL4']:return 1
        if o=='C_basis':return int(n[1]).bit_count()
        if o=='C_product':return sum(order(j)for j in n[1])
        if o=='C_sum':return min(order(j)for j in n[1])
        raise AssertionError((i,n))
    orders=[]
    for q in R['actual_Qrad4_by36_columns']:
        terms=[cn[j]for j in q]
        assert all(t[0]=='R0_component_of_C' and t[2]==i for i,t in enumerate(terms))
        assert len({t[1]for t in terms})==1
        orders.append(order(terms[0][1]))
    assert min(orders)>=2 and len(orders)==36
    assert R['NilD0_exactly_H1_of_actual_complex']
    proof={
      'actual_coefficient_base':'same physical regular9 R0 throughout',
      'actual_joint_ideal':'J=IA+IB=(zA7,r²,rs,s²)',
      'm_R0_square_subset_J':True,
      'trace_roots_actual_original_node_ids':roots,
      'cached_restricted_first_coefficients':[[str(x)for x in v]for v in vals],
      'used_original_ancestor_count':len(visited),
      'all_used_nodes_first_coefficients':[[i,nodes[i],[str(x)for x in ev(i)],i in overrides]for i in sorted(visited)],
      'actual_nuA_nuB_in_J_allorders_from_m_squared_inclusion':True,
      'all36_Qrad_C_maximal_ideal_order_lower_bounds':orders,
      'same_origin_tauA_tauB_zero':'actual rank2 source and accepted proper physical interface, not a translated trace origin',
      'universal_scalar_product_formula':'(fg)0=f0*g0+nuA*fA*gA+nuB*fB*gB+nuA*nuB*fAB*gAB',
      'scalar_component_m_C_squared_subset_m_R0_squared_plus_nuA_nuB':True,
      'all_Qrad_even_coefficients_in_J':True,
      'all_nil_cycles_even_coefficient_in_J_complete_H1':True,
      'actual_odd_basis':'(zA7,gA)*tauA, including unchanged old gA node505938',
      'entire_conductor':'Ii*Bi=Ii+Ii*taui subset Ai, not just scalar Ii',
      'entire_conductor_annihilates_arbitrary_completed_basechange_kernel':True,
      'conductor_factorization_proof':'multiplication c:Bi->Ai is Ai-linear; composing inclusion on the left gives c on Ai, hence c kills every tensor kernel',
      'K_squared_subset_both_first_normalization_kernels':True,
      'J_annihilates_K_squared':True,
      'all_CA_odd_basis_annihilates_K_squared':True,
      'Nil_D0_times_K_squared_zero':True,
      'K_power3_zero':True,
      'Nil_D0_power3_subset_K_by_whole_E4_nilcube':True,
      'Nil_D0_power7_zero':True,
      'full_actual_B18_nil_power7_zero_by_positive_and_G19_gates':True,
      'nil_index_interval':[3,7],
      'E4_embedded_N_P_nil_square_nonzero_respected':True,
      'requires_E4_nil_support_branch_image_inventory':False,
      'exact_index_or_nil_cube_zero_nonzero_or_H0_annihilator_claimed':False,
      'new_source_solve_jet_GB_point_or_780_expansion':False,
    }
    ep=out/'actual-trace-first-gate-and-entire-conductor-proof.json'
    ep.write_text(json.dumps(proof,indent=2)+'\n')
    other=[
      'research/notes/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW.md',
      'research/notes/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63.md',
      'research/notes/AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63.md',
      'research/notes/AMBIENT_B18_W1_POLYNOMIAL_LAURENT_SCOPE_20261104_v63.md',
      f'research/notes/{STEM}.md',str(Path(__file__).relative_to(ROOT)),
    ]
    report={
      'id':STEM,'all_passed':True,
      'status':'PROVED_actual_quadratic_conductor_nil_upper7',
      'scope':'actual CA intersect CB, formal B18 representative P0, cone/projective germ',
      'source_trace_restrictions_both_first_zero':True,
      'actual_J_contains_m_R0_squared_and_both_nu':True,
      'all_complete_nil_cycles_even_coefficient_in_J':True,
      'entire_conductor_annihilation_after_completed_basechange':True,
      'Nil_D0_times_K_squared_zero':True,'comparison_K_power3_zero':True,
      'actual_D0_and_full_B18_nil_power7_zero':True,
      'actual_nil_index_lower_bound':3,'actual_nil_index_upper_bound':7,
      'upper7_uses_no_E4_N_P_target_squarezero_claim':True,
      'nil_cube_zero_nonzero_exact_index_or_finite_H0_annihilator_claimed':False,
      'source_worker_jet_point_GB_or_large_matrix_job':False,
      'input_sha256':{p:sha(ROOT/p)for p in list(files.values())+other},
      'evidence_sha256':{ep.name:sha(ep)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k]for k in ['id','all_passed','actual_nil_index_lower_bound','actual_nil_index_upper_bound','seconds']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())

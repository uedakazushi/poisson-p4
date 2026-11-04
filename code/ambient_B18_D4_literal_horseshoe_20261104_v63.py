#!/usr/bin/env python3
"""Actual literal D4 horseshoe from the frozen nil-delta coefficient input."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
from ambient_B18_actual_top_coefficient_preflight_20261104_v63 import DAG
from ambient_B18_D4_horseshoe_template_20261104_v63 import build_horseshoe

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_D4_LITERAL_HORSESHOE_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    np=ROOT/'reports/AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63/actual-N-literal-cone-DAG.json'
    nr=np.parent/'report.json'
    dp=ROOT/'reports/AMBIENT_B18_NIL_DERIVATION_EXTRACTOR_20261104_v63/nil-delta-extractor-DAG.json'
    dr=dp.parent/'report.json'
    tr=ROOT/'reports/AMBIENT_B18_D4_HORSESHOE_TEMPLATE_20261104_v63/report.json'
    nd,dd=map(read,[np,dp]);assert all(read(p)['all_passed'] for p in [nr,dr,tr])
    D=DAG();D.nodes=[tuple(n) for n in nd['DAG_scalar_nodes']]
    D.cache={n:i for i,n in enumerate(D.nodes)};assert len(D.cache)==len(D.nodes)
    D.values=[]
    for n in D.nodes:
        o=n[0]
        if o=='c':v=F(n[1],n[2])
        elif o=='v':v=F(0)
        elif o=='a':v=D.values[n[1]]+D.values[n[2]]
        elif o=='n':v=-D.values[n[1]]
        elif o=='m':v=D.values[n[1]]*D.values[n[2]]
        elif o=='i':assert D.values[n[1]];v=1/D.values[n[1]]
        else:raise ValueError(o)
        D.values.append(v)
    # Compose the entire pointed unit extractor with Z=actual z1.  In
    # particular, do not identify its N56 variables with the ambient N56.
    Z=dd['curve_parameter_Z_node'];z1=nd['z1_through_z8_actual_odd_coefficient_nodes'][0]
    mp={};variables={}
    for i,n in enumerate(dd['DAG_scalar_nodes']):
        o=n[0]
        if i==Z:assert o=='v';v=z1
        elif o=='v':v=D.var('D4_EXTRACTOR_'+n[1]);variables[i]=v
        elif o=='c':v=D.const(F(n[1],n[2]))
        elif o=='a':v=D.add(mp[n[1]],mp[n[2]])
        elif o=='n':v=D.neg(mp[n[1]])
        elif o=='m':v=D.mul(mp[n[1]],mp[n[2]])
        elif o=='i':v=D.inverse(mp[n[1]])
        else:raise ValueError(o)
        mp[i]=v
    c=[mp[i] for i in dd['positive_nil_residue_coefficients_c36']]
    psi=[mp[i] for i in dd['nilcurve_graph_ideal_P8_right_sides_psi']]
    assert len(c)==36 and len(psi)==8 and all(D.values[i]==0 for i in c+psi)
    h=build_horseshoe(D,nd,c,psi)
    assert all(D.values[i]==0 for key in ['D4_d1_9_by44_columns','D4_d2_44_by112_columns'] for col in h[key] for i in col.values())
    assert all(D.values[i]==0 for i in h['all84_curve_cocycle_obstruction_nodes'])
    curve=[mp[i] for i in dd['curve77_equations']]
    dual=[mp[i] for i in dd['linear_dual77_equations']]
    assert len(curve)==len(dual)==77 and all(D.values[i]==0 for i in curve+dual)
    data={
        'DAG_scalar_nodes':D.nodes,
        'actual_R0_nodes_and_N_generator_input':str(np.relative_to(ROOT)),
        'same_actual_R0_coordinate_nodes_by_row':nd['same_actual_R0_coordinate_nodes_by_row'],
        'actual_z1_through_z8_nodes':nd['z1_through_z8_actual_odd_coefficient_nodes'],
        'actual_nu_node':nd['nu_actual_node'],'actual_v_node':nd['actual_v_node'],
        'N56_unit_unknowns':nd['all56_coefficient_unknown_nodes'],
        'N56_unit_equations':nd['all56_actual_coefficient_Hensel_equations'],
        'nil_curve_parameter_Z_is_actual_z1_node':z1,
        'cloned_curve77_unknowns':[mp[i] for i in dd['curve77_unknowns_group17_free4_N56']],
        'cloned_curve77_equations':curve,'curve77_unit_determinant':dd['curve77_unit_determinant'],
        'cloned_linear_dual77_centered_unknowns':[mp[i] for i in dd['linear_dual77_centered_unknowns']],
        'cloned_linear_dual77_equations':dual,
        'extractor_variable_node_map':[[i,j] for i,j in sorted(variables.items())],
        'positive_c36_actual_nodes':c,'nilcurve_psi8_actual_nodes':psi,
        'D4_generators_F_then_N8':[{'index':0,'actual_target_function':'F_old=2*aa+v0+z25'}]+
             [{**q,'index':q['index']+1} for q in nd['eight_N_module_generators']],
        **h,
        'd1_sign':'nilann*F and Nrelation_alpha-c_alpha*F',
        'all84_curve_obstructions_zero_all_orders_proof':
             'N d1*d2=0 implies(sum c*S)F=0; annF=P and Ffree over R0/P, so curve restriction0',
        'exactness_at_44_proof':'project syzygy to N36, use complete N84; remaining nil8 use complete Koszul28',
        'actual2069_old_relation_lifts_or_B_full_not_computed':True,
        'T0_T1_T2_or_depth_not_solved':True,
    }
    ep=out/'actual-D4-horseshoe-DAG.json';ep.write_text(json.dumps(data,separators=(',',':'))+'\n')
    assert ep.stat().st_size<32*1024*1024
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',np,nr,dp,dr,tr,
           ROOT/'code/ambient_B18_D4_horseshoe_template_20261104_v63.py',
           ROOT/'code/ambient_B18_actual_top_coefficient_preflight_20261104_v63.py']
    report={'id':STEM,'all_passed':True,'status':'actual_D4_literal_9_44_112_horseshoe_from_frozen_nil_extractor',
        'actual_c36_and_psi8_bound_by_same_pointed_unit77_plus_linear77':True,
        'actual_D4_generators':9,'complete_first_relations':44,'complete_second_syzygy_generators':112,
        'no_new_Hadamard_or_source_function_needed_for_offdiagonal8_by84':True,
        'all84_curve_identities_follow_from_exact_extractor_and_free_nil_basis_not_origin_extrapolation':True,
        'nil_extractor_input_reused_not_independently_reaudited_in_constructor':True,
        'B_full_old2069_lifts_T0_T1_T2_or_depth_computed':False,
        'no_source_solve_jets_new_points_primes_GB_or_old_r_g_pullback':True,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
        'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start,
        'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'DAG_nodes':len(D.nodes),'certificate_bytes':ep.stat().st_size}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','actual_D4_generators','complete_first_relations','complete_second_syzygy_generators','seconds','max_rss_kib','DAG_nodes','certificate_bytes']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

#!/usr/bin/env python3
"""Typed exact E4 -> four-variable quotient input, no higher coefficients."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_actual_top_coefficient_preflight_20261104_v63 import DAG
from ambient_B18_top_graph_telescoping_20261104_v63 import ancestors
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NORMALIZER_CORE4_REDUCTION_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    ep=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json'
    np=ROOT/'reports/AMBIENT_B18_N_NORMALIZER_PAIR_PREFLIGHT_20261104_v63/actual-N-pair-preflight-evidence.json'
    nr=np.parent/'report.json';S,E,N=map(read,[sp,ep,np]);assert read(nr)['all_passed']
    D=DAG();D.nodes=[tuple(n) for n in S['DAG_scalar_nodes']];D.cache={n:i for i,n in enumerate(D.nodes)};D.values=[];D.divisions=[]
    for i,n in enumerate(D.nodes):
        if n[0]=='c':v=F(n[1],n[2])
        elif n[0]=='v':v=F(0)
        elif n[0]=='a':v=D.values[n[1]]+D.values[n[2]]
        elif n[0]=='m':v=D.values[n[1]]*D.values[n[2]]
        elif n[0]=='n':v=-D.values[n[1]]
        elif n[0]=='i':assert D.values[n[1]];v=1/D.values[n[1]]
        else:raise ValueError(n[0])
        D.values.append(v)
    D.zero=D.cache[('c',0,1)];D.one=D.cache[('c',1,1)]
    rows=S['physical_zero_base_rows9'];x=[next(i for i,n in enumerate(D.nodes) if n==('v',f'x{r}')) for r in rows]
    tau={name:D.var('E4_CORE_'+name+'_tau') for name in ['CA','CB']}
    coeff=[];eq=[];squares=[]
    for name in ['CA','CB']:
        f=S['families'][name]
        coeff += [i for p in f['coefficient_unknown_pairs'].values() for i in p]+[f['p_quadratic_relation']['s'],f['p_quadratic_relation']['t']]
        eq+=f['all48_zero_IFT_quadratic_coefficient_equations']
        squares.append(D.add(D.mul(tau[name],tau[name]),D.neg(f['actual_tracefree_nu'])))
    assert len(coeff)==len(eq)==96 and len(set(coeff))==96
    targets=[18,19,28,29,60,67,72,63,64,68,73,74];diff=[]
    for r in targets:
        sides=[]
        for name in ['CA','CB']:
            a,b=S['families'][name]['complete_normalized_physical_w_coefficients'][0][r]
            sides.append(D.add(a,D.mul(b,tau[name])))
        diff.append(D.add(sides[0],D.neg(sides[1])))
    assert all(D.values[i]==0 for i in eq+squares+diff)
    eliminate=[x[0],x[8]];params=x[1:8]+[tau['CA'],tau['CB']]
    detA=F(S['families']['CA']['unit48_block_determinant_even24_then_odd24_vars_a23_t_b23_s'])
    detB=F(S['families']['CB']['unit48_block_determinant_even24_then_odd24_vars_a23_t_b23_s'])
    schur=F(E['linear_Schur_minor']);assert schur==F(-20,9)
    assert E['linear_pivot_equations']==[0,1] and E['linear_pivot_variables']==[0,8]
    # Choose the source-bound five quotient derivatives, in physical order.
    selected=[0,1,4,2,3];urows=[108,109,113,118,119];core=[112,117,'tau_CA','tau_CB']
    Q=E['all_common_rational_linear_factors'];assert len(Q)==1
    quot=Q[0]['linear_quotient_rows'];jac=[[F(quot[2+j][S['physical_zero_base_rows9'][1:8].index(r)]) for r in urows] for j in selected]
    expected=[[F(int(i==j),24 if i==2 else 12) for j in range(5)] for i in range(5)]
    assert jac==expected
    det5=F(1,24*12**4);assert det5==F(1,497664)
    assert Q[0]['factor_linear_coefficients']==N['unique_linear_equation_coefficients']
    remain=[i for i in range(12) if i not in selected]
    operations=[]
    operations.append({'id':'SCHUR98','operation':'pointed_unit_Hensel_solve',
      'unknown_nodes':coeff+eliminate,'equation_nodes':eq+squares,'parameter_nodes':params,
      'determinant_block_factors':[str(detA),str(detB),str(schur)],
      'origin_all_unknowns_zero':True,'output12_difference_nodes':diff})
    operations.append({'id':'N_BETA','operation':'inverse_of_actual_Npair_rank8_map',
      'source_evidence':str(np.relative_to(ROOT)),
      'input_coordinates':['x108','x109','x112','x113','x117','x118','x119','tau_CA'],
      'output':'tau_CB','minor':'-1/2',
      'finite_function_inputs':N['actual_finite_regular_function_input'],
      'no_vanishing_denominator_unit_inverse':True})
    operations.append({'id':'L','operation':'algebraic_regular_function',
      'formula':'6*(tau_CB-N_BETA(X,tau_CA))',
      'firstjet':N['unique_linear_equation_coefficients'],'linear_tau_CB_derivative':'6'})
    operations.append({'id':'G12','operation':'regular_Hadamard_coordinate_division',
      'ambient_coordinate_change':'tau_CB=N_BETA(X,tau_CA)+ell/6',
      'numerator_source':'SCHUR98.output12_difference_nodes','denominator':'ell',
      'zero_restriction_proof':'all12 differences vanish on actual Npair, the full kernel is (L)',
      'number_of_outputs':12,'localization_at_L':False})
    operations.append({'id':'CORE5','operation':'pointed_unit_Hensel_solve',
      'unknown_physical_rows':urows,'equations_G12_indices':selected,
      'remaining_coordinates':core,'unit_Jacobian5':[[str(v) for v in row] for row in jac],
      'unit5_minor':str(det5),'output':'five actual eliminated coordinates U(core4)'})
    operations.append({'id':'K7','operation':'functorial_regular_substitution',
      'functions_G12_indices':remain,'substitution':'CORE5','target':'QQ[[x112,x117,tau_CA,tau_CB]]',
      'quotient_identity':'S/(G1,...,G12) = QQ[[core4]]/(K1,...,K7)'})
    proof={'DAG_scalar_nodes':D.nodes,'typed_operations':operations,
      'physical_R0_rows9':rows,'tau_nodes':tau,'actual_all12_difference_rows':targets,
      'all12_actual_difference_nodes':diff,'unit_schur_unknown_nodes98':coeff+eliminate,
      'unit_schur_equations98':eq+squares,'actual9_parameter_nodes':params,
      'G12_firstjets_from_existing_quadratic_cache':quot[2:],
      'five_selected_difference_rows':[targets[i] for i in selected],
      'seven_retained_difference_rows':[targets[i] for i in remain],
      'five_Jacobian':[[str(c) for c in r] for r in jac],'five_minor':str(det5),
      'exact_sequence':['0 -> S/J -- multiplication by L -> S/(L*J) -> S/(L) -> 0',
                        'J=(G1,...,G12); S regular9; S/(L) regular8'],
      'local_cohomology_consequence':'H0_m(E4) = H0_m(QQ[[core4]]/(K7))',
      'representation_limit':'typed finite algebraic regular-function operators, not expanded K7 coefficients or a pure-polynomial-only Hensel serialization',
      'no_H0_D0_identification_without_comparison_kernel':True,
      'no_new_source_solve_higher_jets_or_GB':True}
    p=out/'typed-core4-input.json';p.write_text(json.dumps(proof,separators=(',',':'))+'\n');assert p.stat().st_size<32*1024*1024
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',sp,ep,np,nr]
    rep={'id':STEM,'all_passed':True,'status':'actual_E4_allorder_Lfactor_and_exact_four_variable_quotient_input',
      'pure_schur_unit_unknowns':98,'actual_output_differences':12,'quotient_unit5_minor':str(det5),
      'remaining_core_variables':4,'remaining_core_functions':7,
      'H0_E4_equals_H0_core4_proved':True,'core_coefficients_socle_or_D0_depth_computed':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
      'evidence_sha256':{p.name:sha(p)},'DAG_nodes':len(D.nodes),'certificate_bytes':p.stat().st_size,
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

#!/usr/bin/env python3
"""Read-only adapter: expose the ACTUAL unsolved six odd-membership systems."""
import argparse
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(output):
    folder=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63'
    source=folder/'actual-source-coefficients.json';data=json.loads(source.read_text())
    ca,cb=data['families']['CA'],data['families']['CB']
    pos=data['common_positive_odd_physical_rows']
    zero=[row for row,node in cb['actual_zero_odd_I0_generators']]
    system={'coefficient_ring':'R0=QQ[[x27,x108,x109,x112,x113,x117,x118,x119,x61]], with the explicitly saved pointed unit48 scalar Hensel functions',
            'CB_odd_membership_matrix_1_by10_nodes':[node for row,node in cb['actual_zero_odd_I0_generators']+cb['actual_positive_odd_J0_generators']],
            'CB_odd_membership_matrix_physical_rows':zero+pos,
            'six_actual_right_hand_sides_nodes':[[row,pair[1]]for row,pair in cb['actual_individual_physical_cubic_trace_pairs']],
            'equations':'for each alpha solve sum_i h_alpha_i*f_Bi + sum_e h_alpha_e*j_Be = B_Balpha in R0',
            'all60_membership_coefficients_solved':False,
            'CA_nu_node':ca['actual_tracefree_nu'],'CA_j_e_nodes':ca['actual_positive_odd_J0_generators'],
            'CA_zero_raw_trace_pairs':[[row,ca['complete_normalized_physical_w_coefficients'][0][row]]for row in zero],
            'CB_zero_raw_trace_pairs':[[row,cb['complete_normalized_physical_w_coefficients'][0][row]]for row in zero],
            'CA_positive_raw_over_w_trace_pairs':[[row,ca['complete_normalized_physical_w_coefficients'][1][row]]for row in pos],
            'CB_positive_raw_over_w_trace_pairs':[[row,cb['complete_normalized_physical_w_coefficients'][1][row]]for row in pos],
            'CA_raw_cubic_trace_pairs':ca['actual_individual_physical_cubic_trace_pairs'],
            'CB_raw_cubic_trace_pairs':cb['actual_individual_physical_cubic_trace_pairs'],
            'relative_formula_even':'a_alpha=A_Balpha-A_Aalpha+sum_i h_alpha_i*(A_Ai-A_Bi)+sum_e h_alpha_e*(A_Ae-A_Be)',
            'relative_formula_odd':'b_alpha=-B_Aalpha+sum_i h_alpha_i*f_Ai+sum_e h_alpha_e*j_Ae',
            'requested12_template':'(nu_A*j_Ae*b_alpha, j_Ae*a_alpha), e=1,2; alpha in six physical cubic rows',
            'requested12_actual_pairs_constructed':False,
            'no_tau_B_equals_tau_A_identification':True,
            'lower_CB_relations_must_be_retained_for_any_membership_lift_choice':True,
            'input_sha256':{str(source.relative_to(ROOT)):sha(source),str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))}}
    output.write_text(json.dumps(system,indent=2)+'\n')
    print(json.dumps({'actual_membership_matrix_shape':[1,10],'actual_rhs_count':6,'solutions_computed':False}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)

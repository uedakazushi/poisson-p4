#!/usr/bin/env python3
"""Exact deck/conductor scheme maps and fixed initial-direction audit only."""
import argparse,gzip,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_DECK_CONDUCTOR_BRANCH_PREFLIGHT_20261104_v63'
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    np=ROOT/'reports/AMBIENT_B18_N_NORMALIZER_PAIR_PREFLIGHT_20261104_v63/actual-N-pair-preflight-evidence.json'
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    gp=ROOT/'reports/AMBIENT_B18_CORE_NODAL_QUARTIC_20261104_v63/fixed-quadrics-nodal-quartic.json'
    qp=ROOT/'reports/AMBIENT_B18_CORE_CUBIC_INITIAL_20261104_v63/five-actual-core-cubics-and-binary12.json'
    meta=ROOT/'reports/AMBIENT_CB_ACTUAL_D_CUBIC_SYZYGY_20261103_v63/generic-setup.json.gz'
    N,S,G,Q,M=map(read,[np,sp,gp,qp,meta]);selected=[18,19,60,28,29]
    tauA=list(map(s.Rational,N['families']['CA']['tau_firstjet_in_N_source8']));tauB=list(map(s.Rational,N['families']['CB']['tau_firstjet_in_N_source8']))
    betaTau=tauB[-1]/tauA[-1];assert betaTau==-4
    Ua0=12*betaTau;assert Ua0==-48
    scheme=[]
    for family,expected in [('CA',16),('CB',-s.Rational(1,64))]:
        traces={r['physical_row']:r for r in N['families'][family]['all12_target_minus_even_and_conductor_coefficients']}
        A=s.Matrix([[s.Rational(c) for c in r['conductor_pullback_mixed2'][1:9]] for r in traces.values()]);assert A.rank()==5
        chosen=[0,1,3,5,6]
        mat=s.Matrix([[s.Rational(traces[r]['conductor_pullback_mixed2'][1+j]) for j in chosen] for r in selected])
        assert mat.det()==expected
        f=S['families'][family];odd=dict(f['actual_zero_odd_I0_generators'])
        scheme.append({'family':family,'source8_file':str(np.relative_to(ROOT)),
          'actual_source_function_recipe':N['actual_finite_regular_function_input'][family],
          'actual8_conductor_source_nodes':[[r,odd[r]] for r in [18,19,28,29,60,67,72,63]],
          'selected5_conductor_rows':selected,'selected5_first_matrix':[[str(c) for c in row] for row in mat.tolist()],
          'selected5_unit_minor':str(mat.det()),'source8_linear_rank':5,
          'remaining_three_source_functions_rows':[67,72,63],
          'remaining_three_source_parameters':['N_rho112','N_rho117','N_source_p'],
          'actual_map_to_core':'deckflip tau_CA only' if family=='CA' else 'deckflip tau_CB only',
          'pullback_J_equals_actual_conductor_ideal':True,
          'regular_root_restriction':'Restrict the already established regular root, never reevaluate its 0/0 defining ratio after quotient.'})
    u,w=s.symbols('u w');phi=s.Matrix([s.sympify(a) for a in G['normalization_P1_to_quartic_four_binary_forms']])
    forms={'L0':[5,-20,24,6],'LA':[5,-20,-24,6],'LB':[5,-20,24,-6],'LAB':[5,-20,-24,-6]}
    hyper={name:str(s.factor((s.Matrix([c])*phi)[0])) for name,c in forms.items()}
    factors=[u-w,9*u-44*w,19*u-124*w,43*u-148*w,611*u*u-4792*u*w+7856*w*w]
    allocation=[]
    for factor,ex in zip(factors,[3,2,2,1,2]):
        incident=[]
        for name,polynomial in hyper.items():
            q,r=s.div(s.Poly(s.sympify(polynomial),u,w),s.Poly(factor,u,w))
            if r.is_zero:incident.append(name)
        allocation.append({'binary_factor':str(factor),'exponent':ex,'degree':s.Poly(factor,u,w).total_degree(),
          'conjugate_hyperplanes_containing_direction':incident,'actual_branch_or_length_not_inferred_from_this':True})
    assert [r['conjugate_hyperplanes_containing_direction'] for r in allocation]==[['L0'],['L0','LA'],[],['LA'],['LB']]
    assert M['f_row']==18 and M['base29_rows'][M['rho_base_index']]==108
    def weight(row):
        pairs=[(i,j) for i in range(5) for j in range(i+1,5)]
        mons=[(i,j) for i in range(5) for j in range(i,5)]
        pair=pairs[row//15];mon=mons[row%15]
        return 1+mon.count(2)-pair.count(2)
    assert weight(18)==weight(108)==0 and weight(140)==1
    gates=[ROOT/'reports/AMBIENT_B18_N_NORMALIZER_PAIR_PREFLIGHT_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_NORMALIZER_CORE4_REDUCTION_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_CORE_CUBIC_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_R_COMPLEMENT_20261103_v63_REVIEW.json']
    assert all(read(p)['all_passed'] for p in gates)
    proof={'deck_action_on_actual_S':'R0 fixed; sigmaA(tauA)=-tauA, sigmaB(tauB)=-tauB; the unique monicSchur R0 solution is even in each root',
      'Npair_L_convention':'L=6*(tauB-beta(X,tauA))',
      'LA_on_N':'tauA*UA; UA=6*(beta(X,tauA)-beta(X,-tauA))/tauA, regular odd-Hadamard quotient',
      'UA_origin':str(Ua0),'LB_on_N':'-12*tauB exactly',
      'F_on_sigmaA_N':'-2*bA_j*tauA','F_on_sigmaB_N':'+2*bB_j*tauB',
      'G_on_sigmaA_N':'-2*bA_j/UA','G_on_sigmaB_N':'-bB_j/6',
      'cancellation_scope':'First cancel nonzero tau in normal DOMAIN B_N; then restrict to conductor quotient, even if tau vanishes there.',
      'two_actual_conductor_preimage_scheme_inputs':scheme,
      'frozen_binary12_expression':Q['binarydegree12_factored_expressions'][0],
      'four_conjugate_hyperplane_binary_sections':hyper,'initial_direction_incidence_only':allocation,
      'Qb_transverse_Z4_actual_target_tangent_generators':[18,108],
      'Qb_two_generators_torus_weights':[0,0],'positive_w140_weight':1,
      'known_generic_normalizer_pair_meeting_lengths_full_source':{'R':2,'Pa_extra_CA_pair':1,'Qb_extra_CB_pair':4},
      'proper_fixed_locus_transfer_argument':'Proper local normalization charts extend puncturedPa/Qb flag arcs; connectedGm fixes finite geometric flags; fixed-locus commutes with fiber products; afterG19 cut the Qb Z4 tangent generators are bothweight0 so its entire Artin algebra is weight0.',
      'literal_Pa_Qb_origin_core_arc_coefficients_or_tangent_matching_computed':False,
      'multiplicity_lower7_transfer_application_left_for_separate_acceptance':True,
      'missing_length5_sources':'mainNpair direction u=w(order3) and criticalCA direction9u=44w(order2) are ONLY tangent candidates; no actual branch/genericlength3or2 proved',
      'all_six_branches_e12_lower_bound_graph_redundancy_CM_or_socle_proved':False,
      'no_new_sourcejets_GB_beta_points_or_corehigher_coefficient_evaluation':True}
    p=out/'deck-conductor-source-schemes-and-direction-gates.json';p.write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',np,sp,gp,qp,meta,*gates]
    rep={'id':STEM,'all_passed':True,'status':'exact_two_deck_conductor_scheme_maps_and_scoped_branch_source_preflight',
      'exact_deck_conductor_scheme_identities':2,'two_source8_selected5_minors':['16','-1/64'],
      'conductor_preimage_inputs_three_variables_three_functions_each':True,
      'six_actual_branches_or_multiplicity12_lower_bound_proved':False,
      'literal_origin_arcs_generic_critical_lengths_CM_socle_or_depth_computed':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{p.name:sha(p)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

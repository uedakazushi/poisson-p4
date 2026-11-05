#!/usr/bin/env python3
"""SOURCE8 M1 gauge invariants and literal second-cache eligibility only."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_G3_SOURCE8_M1_CACHE_REDUCTION_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 data=read('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_INPUT_GATE_20261105_v63/actual-eight-source-root-folds-M1-operators-and-positive-AD-budget.json')
 ad=read('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_AD_20261105_v63/actual-eight-old-formal-rell-sensitivities-with-M1-native-retained.json')
 assert read('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_PSI_SOURCE_FACTOR_20261105_v63_REVIEW/report.json')['all_passed']
 gate='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json'
 assert not field(gate,'Mnew0_columns')[23]
 G2=list(map(H.parse,field(gate,'full_G2_linear234_by102_columns')));lam=list(map(H.parse,data['lambda26_two_channels']))
 def left(ch,var):return {j:sum((lam[ch].get(i//9,F(0))*c for i,c in col.items() if i%9==var),F(0)) for j,col in enumerate(G2) if sum((lam[ch].get(i//9,F(0))*c for i,c in col.items() if i%9==var),F(0))}
 assert not left(0,7);assert left(0,8)==left(1,7)
 by={}
 for r in data['eight_source_folds']:by[(r['extra_active_source'],r['contraction'])]=r
 kappas=[(77,F(3,55)),(101,F(0)),(117,F(15,44)),(125,F(3,55))];records=[]
 for extra,k in kappas:
  a=by[extra,'A'];b=by[extra,'B'];ar,ae=list(map(H.parse,a['M1_known_Xeff_r_two_channels'])),list(map(H.parse,a['M1_known_Xeff_ell_two_channels']));br,be=list(map(H.parse,b['M1_known_Xeff_r_two_channels'])),list(map(H.parse,b['M1_known_Xeff_ell_two_channels']))
  v={23:k} if k else {};nv={23:-k} if k else {}
  assert ar==[{},{}] and ae==[v,{}] and br==[nv,{}] and be==[{},v]
  assert not ar[0] and not add(ae[0],br[0]) and not be[0]
  records.append({'extra_active_source':extra,'kappa':str(k),'A_M1_response':'kappa*alpha23','B_M1_response':'-kappa*beta23',
   'W2_rr_full56_zero':True,'W2_rell_full56_zero':True,'W2_ell_squared_killed_in_base':True})
 old=read('reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_INPUT_GATE_20261105_v63/actual104-old-C0-ancestor-cache-bindings-and-Rell4-gate.json')
 oldlabels=set(old['formal_original_C0_second_labels104']);labels=ad['unknown_original_C0_second_label_list'];new=[i for i in labels if i not in oldlabels];assert new==[7148,7688,7943]
 basis=field(H.PATHS['C0'],'actual_C0_basis8_scalar_nodes');assert [basis[i] for i in [0,4,5]]==new
 # Exact coordinate identity is the accepted literal C0 basis construction:
 # basis=Y[:7]+[gA], Y[0..6] are the source CA-F9 first7 coordinates.
 assert old['literal_zero_node0']
 assert read('reports/AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63_REVIEW/report.json')['all_passed']
 assert read('reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63_REVIEW/report.json')['all_passed']
 rho=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63/actual3960-rho-factor-identities-or-obstructions.json')
 assert rho['combined_two_rho_factor_raw33'] and not rho['failures']
 ep=out/'actual-M1-alpha23-beta23-gauge-reduction-and14-cache-bindings.json';H.dump(ep,{
  'all_four_M1_and_gauge_reductions':records,
  'alpha23_definition':'lambda0*(partial_r M1)*e23','beta23_definition':'lambda0*(partial_ell M1)*e23-lambda1*(partial_r M1)*e23',
  'M0_e23_literal_zero':True,'lambda0_G2r_literal_zero':True,'lambda0_G2ell_equals_lambda1_G2r':True,
  'canonical_normalizer_mu_symbols_canceled_only_in_these_invariants':True,
  'W2_covariance_proof':'W2 pureplane rr=0,rell=kappa e23-kappa e23=0,ell²=0. Remainingz*linear hasNF0 forALL6targets, so ANYconstanttriangularnormalizer K*w2 hasNF0, as does knownbottom-w2. No canonicalgaugeparameter isset0.',
  'new_two_cache_dot_source':'TheSAMEtwo accepted raw33bar_rho covectors apply to actualprimitiveVcol23 becauseM0e23=0. Column23=mask3 ofrawphysical67: (-nuB*bB67,nuA*bA67,delta67). No actualdots evaluated.',
  'bar_rho_source_path':'reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63/actual3960-rho-factor-identities-or-obstructions.json',
  'all14_old_second_labels':labels,'eleven_literal_old104_bindings':[i for i in labels if i in oldlabels],
  'three_literal_coordinate_bindings':[[7148,0,'zA0'],[7688,4,'zA4'],[7943,5,'zA5']],
  'three_coordinate_rell_coefficients_literal_zero':True,'node0_constant_from_frozen_old_gate':True,
  'coordinate_source_proof':'acceptedC0producer basis=Y[:7]+[gA], withY theactualCA-F9coordinatevector; allorderprimitiveVfinal7are(-zAi,0,0) inacceptedconormal sourcepeer. NumericIDpositions checkedagainstimmutablebasis8.',
  'no_actual_second_invariant_native_G2quad_or_comparison_values_queried':True,
  'remaining_source_value_gates':['sourceheight/storagegate for two raw33dots andold/native/firstcorrections','contractedpureknownC2 G2quadratic_rell sourcecallback','knownrawactiveG3quadratic4 finiteunitcallback'],
  'actual_T2_HF2_or_length_not_computed':True})
 report={'id':STEM,'all_passed':True,'M1_gauge_invariant_scalar_count':2,'all_old14_second_callbacks_bound_to_frozen_cache_or_literal_zero':True,
  'actual_second_invariant_or_phi_values_computed':False,'new_physical_or_implicit_coefficients':False,
  'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','M1_gauge_invariant_scalar_count','all_old14_second_callbacks_bound_to_frozen_cache_or_literal_zero','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())

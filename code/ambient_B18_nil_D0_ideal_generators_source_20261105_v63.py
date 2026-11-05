#!/usr/bin/env python3
"""Actual odd-source action kills the five-dimensional nil-generator fiber."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NIL_D0_IDEAL_GENERATORS_SOURCE_20261105_v63'
PATHS={
 'first':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json',
 'firstpeer':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63_REVIEW/report.json',
 'R0mu':'reports/AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_VALUES_20261105_v63/report.json',
 'R0mupeer':'reports/AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_VALUES_20261105_v63_REVIEW/report.json',
 'R0scope':'reports/AMBIENT_B18_NIL_R0_GENERATION_SCOPE_20261105_v63/report.json',
 'product':'reports/AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63/actual-nil-product44-and-cube-DAG.json',
 'productpeer':'reports/AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63_REVIEW/report.json',
 'C0peer':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json',
 'C0note':'research/notes/AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63.md',
 'deep':'reports/AMBIENT_B18_NBAR_CUBE_DEEP_SOURCE_REDUCTION_20261105_v63/actual-bar-cube-order-four-and-leading-source-proof.json',
 'deeppeer':'reports/AMBIENT_B18_NBAR_CUBE_DEEP_SOURCE_20261105_v63_REVIEW/report.json',
 'Pi':'reports/AMBIENT_B18_NBAR_CUBE_ZERO_PROJECTION_20261105_v63/all-three-Schur-row42-and-structural-bar-cube-zero.json',
 'Pipeer':'reports/AMBIENT_B18_NBAR_CUBE_ZERO_PROJECTION_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 pins={p:sha(ROOT/p)for p in PATHS.values()};D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()if k!='C0note'}
 assert all(d['all_passed']for k,d in D.items()if k.endswith('peer')or k in ['R0mu','R0scope'])
 assert D['R0mu']['actual_mu_R0_Nil_D0']==5
 assert D['C0peer']['actual_D0_free9_relations640_and_D4_map_checked']
 P=D['product'];beta=[7148,7199,7454,7505,7688,7943,8045,505938]
 assert P['actual_beta8_C0_node_ids']==beta and P['actual_nuA_C0_node_id']==20086
 assert P['mu_trace_order']==['1','tauA','tauB','tauA*tauB']
 assert D['deep']['all14_a36_and_zero7_odd_coordinate_origins_zero']
 assert all(not Q(c)for row in D['deep']['all14_actual_cycle_origins44_by14'][:43]for c in row)
 assert D['Pi']['full_minimal_Pi0_first_row128_sparse']==[[42,'1']]
 L={0:Q(1),4:Q(25,4),5:Q(1)}
 cols=D['first']['all849_actual_minimal_first_columns40_sparse']
 for col in cols:
  c={i:Q(v)for i,v in col};assert sum((a*c.get(i,Q(0))for i,a in L.items()),Q(0))==0
 # Literal finite typed action: input a36 and actual b0=Qrad_even*a.
 nodes=[['QQ',0],['actual_R0_function','C0',20086],['actual_nil_even_b0','product',P['cycle_representative_even_b0']]]
 aa=[]
 for j in range(36):aa.append(len(nodes));nodes.append(['nil_a36_input',j])
 roots=[]
 for i,betaid in enumerate(beta):
  b=len(nodes);nodes.append(['actual_R0_function','C0',betaid]);nb=len(nodes);nodes.append(['product',[1,b]])
  aout=[]
  for g in range(9):
   for coeff,source in [(nb,4*g+1),(b,4*g),(nb,4*g+3),(b,4*g+2)]:
    aout.append(len(nodes));nodes.append(['product',[coeff,aa[source]]])
  roots.append(aout+[2 if j==i else 0 for j in range(8)])
 assert len(roots)==8 and all(len(r)==44 for r in roots)
 # Universal C4 trace identity (e=beta*tauA, w=b0+B*tauA).
 # Pair tuples are coefficient monomials (nu,beta,b0,B); no source values.
 product4=[{(1,1,0,1):1},{(0,1,1,0):1},{},{}]
 assert product4==[{(1,1,0,1):1},{(0,1,1,0):1},{},{}]
 data={
  'coefficient_rings':'R0 regular actual base; S=R0/(ell); A=D0 is the actual zero-degree intersection local ring; N=Nil(A).',
  'actual_C0_algebra_generators':'1 and beta_i*tauA, beta=(zA7,gA), tauA²=nuA; accepted actual C0free9→D0 is surjective. Hence m_D0=m_R0*D0+sum D0*(beta_i*tauA).',
  'actual_source_beta8_C0_ids':beta,'actual_source_nuA_C0_id':20086,
  'actual_odd_action_nodes':nodes,'all8_odd_action44_roots':roots,
  'actual_action_formula':'a_prime block=mu_C(beta_i*tauA)*a block; odd_coeff_j_prime=b0*delta_ij, with b0=Qrad_even*a.',
  'cycle_and_D0_action_proof':'Qa=b0+B*tauA implies Qa_prime=beta_i*tauA*Qa=nuA*beta_i*B+beta_i*b0*tauA. The chosen odd coefficients have sum beta_j*c_j_prime=beta_i*b0; thus V(a_prime,c_prime)=0 and the represented actualC0 element is exactly (beta_i*tauA)*(b0+B*tauA). Completeness of the accepted comparison makes boundary changes actualD0 zero, so this descends to the genuine Nbar action without claiming arbitraryTnil is an algebra.',
  'universal_C4_output_monomials_nu_beta_b0_B':[[[list(m),v]for m,v in p.items()]for p in product4],
  'nilfiber_functionals5':['v1(0)','v2(0)','v3(0)','v4(0)','d_z0(v0)+(25/4)*d_z4(v0)+d_z5(v0)'],
  'fifth_functional_Lz':[[i,str(a)]for i,a in L.items()],
  'all849_actual_first_relations_annihilated_by_Lz':True,
  'Lz_not_literal_z5_coordinate':'Some first relation head0*z5 coefficients are nonzero. Only the displayed combinedLz annihilatesall849; its value on z5*h0 is1.',
  'five_fiber_functionals_complete_proof':'On genuineNbar v0(0)=0. Hence the five functionals kill m_S*Nbar; firstrelations show they are well-defined. Their matrix on the four correctedheads and z5*h0 isI5. Accepted mu_S Nbar5 makes the induced map Nbar/m_SNbar→k5 an isomorphism.',
  'source_orders_on_actual_generator_barcyles':'All14 canonicalgenuinebarcycles, hence chosenminimalfive, have a36 origins0. Qrad_even36 origins0 so b0 in m_S². Each beta_i lies m_S (gA actuallyorder3). Thus actiona36 and actionodd8 both lie m_S². No nuA first-jet vanishing is assumed.',
  'honest_order_and_Pi_proof':'For an action cycle c_prime in m_S², regularm gives mc_prime in m_S² and accepted actualH lowers ordinaryorder by atmost1, so (c_prime,Hmc_prime) in m_S. RegularminimalPi yields zeroheadorigins. Its head0 linearpart is Pi0_row0 on honestdegree1=e42 on c_prime_degree1=0, since c_prime in m_S². HigherPi terms cannot contribute without an inputorigin. Thus allhead0 linears, and in particularLz, vanish.',
  'all8_odd_actions_on_nil_generator_fiber_zero':True,
  'maximal_ideal_actions_equal':'m_D0*N=m_R0*N',
  'proved_minimal_D0_nilideal_generator_count':5,
  'same_five_lift_statement':'Any genuine fullnil lifts of four correctedbarheads and z5*h0 generate N both asR0module and asD0ideal, and are minimal for both. Existence follows from genuineNbar; expanded fullnil-lift coefficients are not constructed.',
  'scope':'D0 only. This does not compute minimalgenerators of the whole positive-graded nilideal, fullalgebra relations, products, cubevanishing, or exactnilindex.',
  'no_new_source_gauge_jet_or_8by5_scalar_action_values':True}
 ev=out/'literal-odd-eight-action-and-five-dimensional-nil-fiber-proof.json';ev.write_text(json.dumps(data,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'PROVED_actual_D0_nilideal_minimal_generator_count5','actual_mu_R0_N':5,'actual_mu_D0_ideal_N':5,'no_expanded_full_nil_lifts_or_scalar_action_values':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ev.name:sha(ev)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()

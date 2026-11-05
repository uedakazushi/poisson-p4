#!/usr/bin/env python3
"""Pure selected-implicit-functional method; no new source coefficients."""
from pathlib import Path
import hashlib,itertools,json,time
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_JOINT135_SELECTED_FUNCTIONAL_METHOD_PREFLIGHT_20261104_v63'
PATHS={
 'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
 'joint_peer':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW/report.json',
 'binding':'reports/AMBIENT_B18_JOINT_NATIVE_THIRD_SOURCE_BINDING_PREFLIGHT_20261104_v63/joint-native-third-source-binding-and-conditional-counts.json',
 'gauge2':'reports/AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63/phase2-canonical270-zz-gauge.json',
 'gauge2_peer':'reports/AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63_REVIEW/report.json',
 'J30':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase2-fixed30-Jacobian.json',
 'J30_peer':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63_REVIEW/report.json',
 'source4_peer':'reports/AMBIENT_B18_T_SOURCE4_ADJOINT9_20261104_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 st=time.monotonic();out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()};dump(out,'phase.json',{'phase':'METHOD_ONLY_INPUTS_PINNED','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert all(D[k]['all_passed'] for k in ['joint_peer','gauge2_peer','J30_peer','source4_peer'])
 B=D['binding'];base=[tuple(m) for m in B['joint_base25_monomials']];old=[m for m in base if len(m)<3];high=[m for m in base if len(m)==3]
 assert len(base)==25 and len(old)==20 and high==[(0,7,7),(1,7,7),(4,7,7),(5,7,7),(7,7,7)]
 for m in base:
  for k in range(len(m)):assert m[:k]+m[k+1:] in base
 assert all(tuple(sorted(m+(i,))) not in base for m in high for i in range(8))
 assert all(tuple(sorted(m+n)) not in base for m in high for n in high)
 # Ordinary highest-degree increments form an ideal killed by all base
 # parameters. Hence their coefficients enter every rational unit expression
 # only through its center differential; lower response remains nonzero.
 rules=[]
 for a,b in itertools.product(range(4),repeat=2):
  carry=a&b;mask=a^b;powers=tuple(int(bool(carry&(1<<k))) for k in range(2))
  terms=[]
  for k,power in enumerate(powers):
   if power:
    rem=list(powers);rem[k]-=1;terms.append({'nu_powers':rem,'increment':'delta_nuA' if k==0 else 'delta_nuB','coefficient':power})
  assert len(terms)==carry.bit_count()
  rules.append({'left_mask':a,'right_mask':b,'output_mask':mask,'nu_powers':powers,'coefficient_differential_terms':terms})
 assert rules[5]['left_mask']==rules[5]['right_mask']==1 and rules[5]['coefficient_differential_terms']==[{'nu_powers':[0,0],'increment':'delta_nuA','coefficient':1}]
 assert rules[15]['nu_powers']==(1,1) and len(rules[15]['coefficient_differential_terms'])==2
 u=D['joint']['all105_unknown_original_node_order'];assert len(u)==len(set(u))==105
 ev=D['gauge2']['actual_even15_base20'];od=D['gauge2']['actual_odd15_base20'];assert len(ev)==len(od)==15 and all(len(p)==20 for p in ev+od)
 gauge_order=[{'parity':parity,'coefficient_position':i} for parity in ['even','odd'] for i in range(15)]
 evidence={'inputs':inputs,'proof_status':'EXACT METHOD, ACTUAL Jg_phi/RHS/SELECTED TARGETS STILL REQUIRE PRIMARY SOURCE VERIFICATION',
  'original_physical105_unknown_order':u,'canonical_scalar_gauge30_order':gauge_order,
  'canonical_even_odd_origins':[[str(p[0]) for p in ev],[str(p[0]) for p in od]],
  'downward_closed_base25':base,'highest_cubic_monomials':high,'highest_cubic_ideal_square_zero_and_base_maximal_ideal_annihilates_it':True,
  'monic_cover_basis':['1','tauA','tauB','tauA*tauB'],'all16_monic_coefficient_differential_rules':rules,
  'block_J135':[['actual_J105','0_105_by30'],['actual_J_gauge_physical30_by105','actual_J30']],
  'required_lower_block':'differentiate source-dependent monic nuA and nonzero odd gauge center before extracting 30 scalar constraints',
  'adjoint_conditions':['rho_gauge J30=v_gauge','rho_physical J105+rho_gauge J_gauge_physical=v_physical'],
  'exact_requested_formula':'target[m]=target_lower_response[m]-rho*constraint_residual_with_requested_unknown_cubic_coefficients_zero[m]',
  'proof':'J135*u_m+f_m=0; target_m=v*u_m+q_m; rho*J135=v. All coefficients ordinary Taylor coefficients, no factorial conventions.',
  'pure_r3_crosscheck_scope':'all105 physical source r3 already cached; not a complete previously cached canonical gauge30 r3 vector',
  'remaining_actual_gates':['source-equation and canonical-gauge exact input/rational-unit ancestry',
    'actual Jg_phi with monic differential carry', 'all requested rhoJ135=v identities',
    'actual source/gauge residuals from known lower input and direct target lower response',
    'native N3 and cofactor m2 cancellation in the selected coherent source combinations',
    'coherent copied57 third chain rule, without independent copied high coefficients',
    'total live/output/resource and requested-index guard'],
  'new_coefficients_computed':0,'actual_selected_third_values_computed':False,'HF4_or_total_length_claimed':False}
 ep=out/'selected135-implicit-functional-proof-and-required-source-gates.json';dump(out,ep.name,evidence)
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'scope':'pure method and actual scalar-coordinate binding only',
  'actual_J135_selected_rhs_or_values_verified':False,'new_source_or_gauge_coefficients':0,
  'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
  'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-st}
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL_METHOD_ONLY','report_sha256':sha(out/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Restricted118 pure quartic source formula; no source coefficient evaluation.

The scalar extension on the known cubic initial subspace is checked on all
1332 linear generators and all eight multiples of the canonical quadratic.
The actual118 source combinations are retained as factored finite operators.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,resource,signal,time

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_RESTRICTED118_QUARTIC_PREFLIGHT_20261104_v63'
PATHS={
 'X2only':'reports/AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63/actual119-cubic-X2-only-source-operator-preflight.json',
 'cubic_values':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_VALUES_20261104_v63/actual-linear-syzygy-cubic-source-values.json',
 'cubic_peer':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_VALUES_20261104_v63_REVIEW/report.json',
 'offset':'reports/AMBIENT_B18_T_QUADRATIC_L_SYZYGY_R4_20261104_v63/actual-quadratic-L-syzygy-scalar-offset-certificate.json',
 'zi_qB':'reports/AMBIENT_B18_T_ZI_QB_QUARTIC_EXTENSION_20261104_v63/actual-seven-zi-qB-scalar-offset-certificate.json',
 'rewrite':'reports/AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63/linear-module-rewrite-certificate.json',
 'known40':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63/known-initial40-and-length-upper17.json',
 'quadratic_values':'reports/AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63/actual-three-by620-quadratic-image-certificate.json',
 'Theta2':'reports/AMBIENT_B18_T_CONSTANT_QUADRATIC_KERNEL_R3_ADJOINT_20261104_v63/phase1-theta-and-pure-r3-second-adjoint.json',
 'd3':'reports/AMBIENT_B18_T_CUBIC_BOUNDARY_SOURCE3_PREFLIGHT_20261104_v63/phase1-corrected-d3-pure-r3-channel.json',
 'pivots':'reports/AMBIENT_B18_T_CUBIC_BOUNDARY_SOURCE3_PREFLIGHT_20261104_v63/phase2-known-Schur-second-r-gates.json',
 'HF3peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json'}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rec(v):return {int(i):Q(c) for i,c in v if Q(c)}
def enc(v):return [[i,str(c)] for i,c in sorted(v.items()) if c]
def dot(a,b):return sum((c*b.get(i,Q(0)) for i,c in a.items()),Q(0))
def dump(out,name,v):(out/name).write_text(json.dumps(v,separators=(',',':'))+'\n')

def main():
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 dump(out,'phase.json',{'phase':'source_and_scalar_extension_inputs_saved','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert D['cubic_peer']['all_passed'] and D['HF3peer']['all_passed']
 A,V,O,N,K=[D[k] for k in ['X2only','cubic_values','offset','rewrite','known40']]
 assert A['all119_full56_x2_equals_Xp0_a2'] and V['actual_new_cubic_relation_channel_rank']==1
 assert not O['nonzero_quartic_relation_witness']
 linear=list(map(rec,O['normalized_actual37_linear_relations']))
 h=list(map(Q,O['normalized_actual37_e5_r2_values']))
 mons2=list(itertools.combinations_with_replacement(range(8),2))
 mons3=list(itertools.combinations_with_replacement(range(8),3))
 m3={m:i for i,m in enumerate(mons3)}
 assert mons3==[tuple(m) for m in N['monomials3_120']]
 theta0={239:Q(-55495,20736),359:Q(27785,1536)}
 thetaB={359:Q(1,4),479:Q(1)}
 receipts=[]
 for p,L in enumerate(linear):
  for j,mon in enumerate(mons2):
   cubic={}
   for row,c in L.items():
    target,var=divmod(row,8);i=120*target+m3[tuple(sorted(mon+(var,)))]
    cubic[i]=cubic.get(i,Q(0))+c
   expect=h[p] if mon==(7,7) else Q(0)
   assert dot(theta0,cubic)==expect and dot(thetaB,cubic)==0
   receipts.append([36*p+j,str(expect)])
 assert len(receipts)==1332
 qb=K['homogeneous_module_generators40'][38]
 expected_qb={(3,(0,0,0,0,0,0,0,2)):Q(1),
              (4,(0,0,0,0,0,1,0,1)):Q(11767,3652)}
 assert {(t,tuple(m)):Q(c) for t,m,c in qb}==expected_qb
 qb_receipts=[]
 for var in range(8):
  cubic={}
  for target,exponents,c in qb:
   mon=tuple(i for i,n in enumerate(exponents) for _ in range(n))+(var,)
   i=120*target+m3[tuple(sorted(mon))];cubic[i]=Q(c)
  assert dot(theta0,cubic)==0 and dot(thetaB,cubic)==(1 if var==7 else 0)
  qb_receipts.append([var,'0',str(dot(thetaB,cubic))])
 vals=[list(map(Q,r['two_actual_cubic_relation_values'])) for r in V['all119_two_actual_cubic_relation_values']]
 assert len(vals)==119 and all(v[1]==0 for v in vals)
 pivot=71;scale=vals[pivot][0]
 assert scale==Q(-251099375,1051776)
 combinations=[]
 for j,v in enumerate(vals):
  if j==pivot:continue
  combo={j:Q(1)}
  if v[0]:combo[pivot]=-v[0]/scale
  assert all(sum((c*vals[i][k] for i,c in combo.items()),Q(0))==0 for k in range(2))
  combinations.append(enc(combo))
 assert len(combinations)==118
 qcert=D['quadratic_values']['complete_fixed_rank_certificate']
 assert qcert['selected_original_columns']==[563]
 assert Q(qcert['selected_minor_determinant'])==Q(-250,3)
 theta2=rec(D['Theta2']['theta_extension180'])
 assert theta2=={71:Q(-55495,20736),107:Q(27785,1536)}
 assert all(D['d3'][k] for k in ['actual_d3_linear_r_coefficient_identically_zero','all84_Ystar1_source55_r_zero'])
 assert not D['d3']['nonzero_boundary_channel_indices']
 # The r-multiple of Theta2 is exactly theta0.  This is a source-normalized,
 # coefficient-level identity, rather than an arbitrary extension on standards.
 for index,c in theta2.items():
  target,m=divmod(index,36)
  assert theta0[120*target+m3[tuple(sorted(mons2[m]+(7,)))]]==c
 source_support=sum(len(rec(x)) for x in A['actual_selected37_original768_constant_combinations'])
 source_support+=sum(len(rec(x)) for x in A['all119_actual_a2_sparse27_times36'])
 source_support+=sum(len(x) for x in combinations)+len(receipts)+sum(map(len,linear))
 assert source_support<150000
 dump(out,'phase.json',{'phase':'all1332_plus8_scalar_extension_identities_and118_kernel_checked',
  'QQ_support_bound':source_support,'actual_source_coefficient_callback_evaluated':False})
 fields={
  'actual_base':'QQ[[zA0,...,zA6,rA]]; last target=e4, source row55',
  'source_inputs':inputs,
  'actual_selected37_original768_constant_combinations':A['actual_selected37_original768_constant_combinations'],
  'actual_selected37_minimal657_indices':A['actual_selected37_minimal_relation_indices'],
  'all119_linear_source_syzygies296':A['all119_linear_syzygies_source296'],
  'all119_actual_a2_sparse27_times36':A['all119_actual_a2_sparse27_times36'],
  'actual118_constant_combinations_of119':combinations,
  'cubic_image_pivot119':pivot,'cubic_image_pivot_value':str(scale),
  'Theta3_on_cubic600':{'constant_part':enc(theta0),'canonical_qB_pure_cubic_coefficient_part':enc(thetaB)},
  'all1332_generator_checks_expected_values':receipts,
  'all8_qB_multiples_checks_Theta0_and_ThetaB':qb_receipts,
  'canonical_qB_initial':qb,
  'canonical_qB_source_normalization':{'actual620_index':563,'divide':str(Q(-250,3)),
   'subtract':'the chosen N177 lift using Theta2_extension180',
   'b3_formula':'F620_pure_third[563]/(-250/3)',
   'F620_formula':'rawP55,3-Xp55,2rr*a1_r-(Theta2+Ystar55,2rr*beta1_r)(Xc2-Xp1*a1)',
   'Theta2_vanishes_on_canonical_qB':True},
  'typed119_C3':'sum(linear_ell*actual_X2(selected37_source))-Xp1*actual_a2',
  'typed119_beta2_rr':'r² coefficient of U0^-1*kappa_degree3(m0*C3)',
  'typed119_Qproxy3':'pi*C3-Ystar1*U0^-1*kappa_degree3(m0*C3)',
  'typed119_raw_pure_quartic':'sum(ell_r*actual_raw_P55_source3)-sum(Xp55,2rr*a2_rr)-sum(Ystar55,2rr*beta2_rr)',
  'typed118_corrected_pure_quartic':'kernel118*(raw_pure_quartic-Theta3(Qproxy3))',
  'raw_P55_source3_rule':'sum original768 constant source weights; native columns640..767 contribute identically zero',
  'known_Xp_and_Ystar_source55_second':D['pivots'],
  'contraction':'wi8=(E+1)^-1 partial_i(v8); wij=(E+2)^-1(partial_i(vj)-partial_j(vi)); E on OUTPUT',
  'Euler_boundary_reason':'Both corrected wedge expressions start in m². Constant d3_linear on wedge³ constants is injective, so an exact d3-boundary has a preimage in m. Its pure quartic image is zero because actual minimal d3 pure degree2 and degree3 are zero.',
  'no_native_X3_N2_L2S1_or_commonCB_source4_in_this_restricted_formula':True,
  'no_Xp3_or_Ystar3_in_this_restricted_formula':'first27 correction a and second84 correction beta both start in degree2; Xp55 and Ystar55 start in degree2',
  'native_effective_S2_use':'Only source-selected zero-origin constant combinations; the accepted L2*S0 cancellation applies before multiplying by ell.',
  'remaining_actual_coefficient_gate':'Adjoint/evaluate raw P3, existing X2 pairs and the source-normalized scalar b3. Keep any uncached ordinary quotient third derivatives symbolic until aggregate cancellation.',
  'this_is_only_the_restricted118_family':True,
  'whole_quartic_initial_space_or_HF4_or_T_length_computed':False,
  'QQ_support_bound':source_support}
 ep=out/'restricted118-quartic-source-operator-preflight.json';dump(out,ep.name,fields)
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'source_kernel_combinations':118,
  'Theta3_generator_identities_checked':1340,'restricted_pure_quartic_formula_proved':True,
  'native_X3_N2_commonCB4_or_new_source_coefficients_evaluated':False,
  'quartic_relation_values_HF4_or_T_length_computed':False,
  'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
  'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start,
  'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL','report_sha256':sha(out/'report.json')})
 print(json.dumps(report),flush=True)

if __name__=='__main__':main()

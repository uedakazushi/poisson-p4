#!/usr/bin/env python3
"""Source-specific sensitivity only; no new source/gauge coefficient values."""
import collections,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_honest_nil_T_first_values_20261104_v63 import mv
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_LINEAR_SOURCE_WHOLE_PREFLIGHT_20261104_v63'
PATHS={
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'minimal_peer':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63_REVIEW/report.json',
 'm0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json',
 'X1':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'W0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase3-actual-W84-rank-and-mu.json',
 'oldB':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'atoms':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'source2':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'physical_adapter2':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'B11':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase3-actual-proper-G8-B11.json',
 'native1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase5-native-adapter-firsts.json',
 'proper1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase4-native-source-and-R12.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'source2_peer':'reports/AMBIENT_B18_T3_ZERO_STRUCTURAL_20261104_v63_REVIEW/report.json',
 'adapter_peer':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(x):
 if isinstance(x,Q):return str(x)
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 return x
def add(a,b,c=Q(1)):
 a=dict(a)
 for i,v in b.items():
  a[i]=a.get(i,Q(0))+c*v
  if not a[i]:del a[i]
 return a
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=True);assert not(out/'report.json').exists()
 (out/'phase.json').write_text(json.dumps({'phase':'source_manifest_before_sensitivity','input_sha256':{p:sha(ROOT/p)for p in PATHS.values()}},indent=2)+'\n')
 D={k:read(ROOT/p)for k,p in PATHS.items()};assert D['minimal_peer']['all_passed']
 M,X,W,Min=(D[k]for k in ['m0','X1','W0','minimal'])
 assert [(i,Q(c))for i,c in enumerate(M['m0_24_by44'][7])if Q(c)]==[(20,Q(5,2))]
 assert Min['d3_only_linear']['constant_projection3_by84']==[[str(int(i==j))for i in range(84)]for j in [6,24,26]]
 monsets={i:[tuple(sorted((i,j)))for j in list(range(i+1))+[7]]for i in [0,4,5]}
 monomials=list(dict.fromkeys(m for p in monsets.values()for m in p));assert len(monomials)==15
 cached=[m for m in monomials if m[1]==7];missing=[m for m in monomials if m[1]!=7]
 assert cached==[(0,7),(4,7),(5,7)]and len(missing)==12
 wr=W['W_pivot_rows'];wp=W['W_pivot_residual_columns'];Wi=[[Q(x)for x in row]for row in W['W_fixed_minor_inverse']]
 original=W['residual_original_X_columns'];coeff=[[Q(c)for c in row]for row in X['all780_X0_factorization_coefficients']];xp=X['X0_pivot_columns']
 # Unknown adapter2 contribution must survive both constant elimination stages,
 # not merely the first X0 elimination.  Prove all final681 Bold0 images zero.
 inter={j:{i:Q(c)for i,c in sp}for j,sp in Min['full_first_layer_missing_data_profile']['X2_data']['constant_residual_Bold0_nonzero_list']}
 Wfactor=[]
 for sp in W['all_W84_columns_sparse']:
  q={i:Q(c)for i,c in sp};Wfactor.append(mv(Wi,[q.get(i,Q(0))for i in wr]))
 finalids=[j for j in range(762)if j not in wp];assert len(finalids)==681
 for j in finalids:
  v=inter.get(original[j],{})
  for c,k in zip(Wfactor[j],wp):
   if c:v=add(v,inter.get(original[k],{}),-c)
  assert not v,('ACTUAL_ADAPTER2_CONSTANT_NOT_CANCELLED',j,v)
 # Origin derivatives of the old scalar operator give all unknown atom2
 # sensitivities. No second scalar value is assigned.
 A=D['atoms'];atoms={i:{k:Q(c)for k,c in v}.get(0,Q(0))for i,v in A['all79_atom_values']}
 B=D['oldB'];ops=B['operator_scalar_nodes'];columns=B['B48_by640_sparse_columns']
 AA=D['X1']['new_A12_by9_C36'];weights=[Q(row[5][0])for row in AA]
 roots=[v for col in columns for row,v in col if row%4==0 and weights[row//4]]
 needed=set();todo=list(roots)
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i);n=ops[i]
  if n[0]=='sum':todo.extend(j for j,a,b in n[1])
  elif n[0]=='product':todo.extend(n[3])
 val={};der={}
 for i in sorted(needed):
  n=ops[i];kind=n[0]
  if kind=='c':v=Q(n[1],n[2]);d={}
  elif kind=='actual_C0_scalar_function':v=atoms[n[1]];d={n[1]:Q(1)}
  elif kind=='sum':
   v=sum((Q(a,b)*val[j]for j,a,b in n[1]),Q(0));d={}
   for j,a,b in n[1]:d=add(d,der[j],Q(a,b))
  elif kind=='product':
   factors=n[3];pref=[Q(n[1],n[2])]
   for j in factors:pref.append(pref[-1]*val[j])
   suff=[Q(1)]*(len(factors)+1)
   for k in range(len(factors)-1,-1,-1):suff[k]=suff[k+1]*val[factors[k]]
   v=pref[-1];d={}
   for k,j in enumerate(factors):d=add(d,der[j],pref[k]*suff[k+1])
  else:raise AssertionError(n)
  val[i]=v;der[i]=d
 raw=[]
 for col in columns:
  d={}
  for row,v in col:
   if row%4==0 and weights[row//4]:d=add(d,der[v],weights[row//4])
  raw.append(d)
 raw += [{}for _ in range(140)]
 residual=[]
 for j in original:
  v=raw[j]
  for c,k in zip(coeff[j],xp):
   if c:v=add(v,raw[k],-c)
  residual.append(v)
 sensitivity=[]
 for j in finalids:
  v=residual[j]
  for c,k in zip(Wfactor[j],wp):
   if c:v=add(v,residual[k],-c)
  sensitivity.append([[i,c]for i,c in sorted(v.items())])
 retained=sorted({i for row in sensitivity for i,c in row});source_nodes=len(D['source']['DAG_scalar_nodes'])
 pure=[i for i in retained if i<source_nodes];nonpure=[i for i in retained if i>=source_nodes]
 assert source_nodes==41315
 assert not retained,'retained independent atom second derivatives require a separate source gate'
 # Whole cache-worker input gate. These are schema/origin/support checks, not
 # second-coefficient evaluation or a source/gauge solve.
 for name in ['source2_peer','adapter_peer']:assert D[name]['all_passed']
 J=D['source2'];fams=D['source']['families'];nodes=D['source']['DAG_scalar_nodes']
 ids={f:set(J['both57_source2'][f]['unknown57'])for f in ['CA','CB']}
 assert len(ids['CA'])==len(ids['CB'])==57 and len(ids['CA']&ids['CB'])==9 and len(ids['CA']|ids['CB'])==105
 for f in ['CA','CB']:
  sf=J['both57_source2'][f]
  assert sf['all57_equations_all55_coefficients_verified']
  assert all(len(v)==55 for i,v in sf['all57_source2jets'])
 assert len(D['physical_adapter2']['YB_in_YA_quadratic'])==9
 physicalrows=B['physical_target_rows12'];assert physicalrows==[18,19,28,29,60,67,72,63,64,68,73,74]
 targetroots=[]
 for f in ['CA','CB']:
  family=fams[f]
  targetroots+=family['all48_zero_IFT_quadratic_coefficient_equations']+J['both57_source2'][f]['F9_nodes']
  targetroots+=[family['actual_tracefree_nu']]
  targetroots += [i for row in physicalrows for i in family['complete_normalized_physical_w_coefficients'][0][row]]
 sr=set();todo=list(targetroots)
 while todo:
  i=todo.pop()
  if i in sr:continue
  sr.add(i);n=nodes[i]
  if n[0]in['a','m','n','i']:todo.extend(n[1:])
  else:assert n[0]in['v','c']
 leaves={i for i in sr if nodes[i][0]=='v'};assert leaves<=ids['CA']|ids['CB']
 def maxlive(order,dependencies,roots):
  uses=collections.Counter(roots)
  for i in order:uses.update(dependencies(i))
  live=set();peak=0
  for i in order:
   live.add(i);peak=max(peak,len(live))
   for j in dependencies(i):
    uses[j]-=1
    if not uses[j]:live.remove(j)
  return peak
 sdep=lambda i:nodes[i][1:]if nodes[i][0]in['a','m','n','i']else[]
 smax=maxlive(sorted(sr),sdep,targetroots)
 def bdep(i):
  n=ops[i]
  return [j for j,a,b in n[1]]if n[0]=='sum'else n[3]if n[0]=='product'else[]
 bmax=maxlive(sorted(needed),bdep,roots)
 N=D['native1'];P=D['proper1']
 assert all(list(map(Q,v[:4]))==[0,0,0,Q(v[3])]for v in N['G8'])
 assert list(map(Q,P['kappa9'][1][:4]))==[0,-2,Q(-1,2),0]
 assert all(not Q(v[0])and not Q(v[3])for v in [P['kappa9'][i]for i in [0,1,2,6]])
 assert all(not Q(c)for v in P['physicalF12']for c in v[:4])
 for col in D['radical']['native_tensor_Koszul6_HB2_resolution'][1]['columns']:
  for row,pol in col:
   if row==5:assert all(not mon[5]and sum(mon)<=2 for mon,c in pol)
 # C4 coefficient proof of the unit recovery, retaining every unknown trace slot.
 def c4mul(a,b):
  v=[Q(0)]*4
  for i,x in enumerate(a):
   for j,y in enumerate(b):
    if not i&j:v[i|j]+=x*y
  return v
 d0=[Q(0),Q(-2),Q(-1,2),Q(0)]
 for k in range(4):
  u=[Q(int(j==k))for j in range(4)]
  assert c4mul(d0,u)[1]==Q(-2 if k==0 else 0)
  assert c4mul(u,[Q(0),Q(0),Q(0),Q(1)])[1]==0
 budget={
  'all105_composed_source2_cache_slots':105*45,
  'all164_requested_scalar_source_values_upper':164*45,
  'source_AST_maximum_live45_vectors':smax*45,
  'old_scalar_AST_maximum_live45_vectors':bmax*45,
  'all79_existing_atom_firsts':79*9,
  'all640_weighted_old_row20_quadratic15_coeffs':640*15,
  'all140_native_row20_quadratic15_coeffs':140*15,
  'cached_N12_C36_and_row20_adapter_first':144*36+48*9,
  'F12_C180_and_u5_even2_scratch':12*180+5*36,
  'all780_raw_row20_quadratic15_coeffs_optional':780*15,
  'all849_minimal_first_relation40_coeffs_upper':849*40,
  'streaming_constant_elimination_and_fixed_rank40_minors':12000}
 assert sum(budget.values())<150000,(sum(budget.values()),budget)
 results={
  'same_S_frame':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA'],
  'm0_row7_single_actual_source':[20,Q(5,2)],'final_wedge_projection_is_exact_coordinate_selection':[6,24,26],
  'needed_m1_row_only':7,'canonical_m1_row7':'only r derivative of Vgamma_A(z=0,r)/r; coefficient is actual Vgamma_A pure-r² in saved C44',
  'm1_row7_z_derivatives_zero_by_canonical_definition':True,
  'no_unknown_m1_third_sector_or_pure_r3_needed_for_final_three_wedges':True,
  'm1_term_retained_as_actual_cached_callback_not_set_zero':True,
  'all15_needed_row20_X2_monomials':monomials,'existing_B11_monomials3':cached,'outside_B11_monomials12':missing,
  'ordered_Hadamard_rule':'target i=0,4,5 reads Xrow20 coefficient zi*zj only for j<=i or j=r; absent prefix coefficients are zero by exact restriction, not by data omission',
  'all681_actual_final_constant_source_combinations_have_Bold0_zero':True,
  'intermediate_nine_Bold0_residuals_still_nonzero':True,
  'entire_unknown_ResAT2_Bold0_term_zero_AFTER_both_constant_source_cancellations':True,
  'constant_source_definition':'c_j=(e_originalj−Xp*X0factor_j)−sum_k Wfactor_jk*(e_originalWp_k−Xp*X0factor_Wp_k)',
  'final681_original_X_column_labels':[original[j]for j in finalids],
  'retained_old_atom_second_derivative_sensitivities681':sensitivity,
  'retained_old_atom_second_nodes':retained,'retained_pure_source_second_nodes':pure,'retained_nonpure_second_nodes':nonpure,
  'old_atom_second_value_not_computed':True,
  'proper_native_row20_second_availability':{
   'native_source_row':'h generator coefficient, R9=(u5,h,cr,cl,rl)',
   'native_S35_rowh_coefficients_only':'u5 and cr,cl,rl; never h itself',
   'quadratic_core_scalar_second':'origin core variables have only tauA/tauB components; unknown core base2 changes have zero even component in their quadratic products, since nuA(0)=nuB(0)=0',
   'u5_even_second_unit_recovery':'Use d*G8=N12*F12, F12_0=0, G8_0=g*tauA*tauB, d0=−2*tauA−tauB/2. In pure-tauA component at base degree2, the only unknown G2 term is −2*G_even2; unknown N2F0 and d2G0 give zero in this component. Thus solve with known unit −2, never divide by d.',
   'required_unit_recovery_inputs':'all pure physical F12 seconds and monics2 from two cached57 solutions plus exact physical-F9 adapter2; N0/N1,d0/d1,G0/G1 already saved',
   'new_gauge2_solve_required_by_this_recovery':False},
  'actual_coefficient_worker_run':False,'new_source_gauge_jet_or_second_scalar_value':False,
  'new_visible_source_sensitivity_slots':sum(len(x)for x in sensitivity)+sum(len(x)for x in der.values()),
  'whole_cache_input_gate':{
   'joint105_source2_leaf_union_verified':True,'source_AST_closure_count':len(sr),'source_AST_maximum_live':smax,
   'old_scalar_AST_needed_count':len(needed),'old_scalar_AST_maximum_live':bmax,
   'literal_source_evaluation_root_order':targetroots,'literal_source_AST_closure_order':sorted(sr),
   'source2_same_physical_frame':'compose CB own57 source2 with frozen full9 YB(YA)2, then ell sA=rA; shared nine raw physical roots must agree exactly with CA57',
   'no_unknown79_second_callbacks_retained':True,'all_C4_origin_trace_unit_recovery_checks':True,
   'native_rowh_no_h_variable_and_degree_at_most2_checked':True},
  'whole_worker_budget_saved_plus_live_QQ':budget,'whole_worker_QQ_upper':sum(budget.values()),
  'worker_hard_bounds':{'seconds':60,'address_bytes':2<<30,'output_bytes':32<<20,'visible_QQ_slots':150000},
  'remaining_gate_before_any_cache_value_worker':'all source/cache/unit/sensitivity gates now closed; one cache-only value worker may be run under root conditional authorization; no source/gauge equation may be solved'}
 assert results['new_visible_source_sensitivity_slots']<150000
 ev=out/'actual-final681-source-sensitivities-and-cached-directions.json';blob=json.dumps(enc(results),separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;ev.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'SOURCE_SENSITIVITY_PREFLIGHT_no_new_coefficient_values',
  'all681_adapter_second_terms_cancel_after_actual_constant_eliminations':True,
  'needed_X2_row_only':20,'needed_X2_monomials':15,'already_B11':3,'outside_B11':12,
  'remaining_old_atom_second_node_count':len(retained),'remaining_nonpure_second_nodes':nonpure,
  'all_cache_only_worker_source_and_budget_gates_closed':True,'cache_worker_QQ_upper':sum(budget.values()),
  'old_second_atoms_or_full_first_relation_values_computed':False,
  'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ev.name:sha(ev)},
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()

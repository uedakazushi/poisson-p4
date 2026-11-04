#!/usr/bin/env python3
"""119 true cubic relation classes from actual X2 only; no X2 evaluation.

Checks all actual first-order hypotheses and constructs the exact quadratic
pivot correction a2.  Cubic-dual/adjoint extraction is a separate peer task.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63'
PATHS={
 'schur':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/actual-minimal5-by825-Schur-operator.json',
 'schur_peer':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW/report.json',
 'X0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
 'V':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
 'first640':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
 'native128':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json',
 'm1':'reports/AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/actual-projected-three-class-X2-operator-preflight.json',
 'linear':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
 'syzygies':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_PURE_R3_20261104_v63/actual-linear-syzygies-and-pure-r3-source-values.json',
 'cubic_rewrite':'reports/AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63/linear-module-rewrite-certificate.json',
 'cubic_closure':'reports/AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63/actual-known-quadratic-cubic-closure-certificate.json',
 'cubic_report':'reports/AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rec(x):return {int(i):Q(c) for i,c in x if Q(c)}
def enc(x):return [[i,str(c)] for i,c in sorted(x.items()) if c]
def add(a,b,c=Q(1)):
 a=dict(a)
 for i,v in b.items():
  z=a.get(i,Q(0))+c*v
  if z:a[i]=z
  else:a.pop(i,None)
 return a
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 (out/'phase.json').write_text(json.dumps({'phase':'cubic_source_homotopy_inputs_saved_before_math','inputs':inputs},indent=2)+'\n')
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 S,X,V,M,L,Y,C=[D[k] for k in ['schur','X0','V','m1','linear','syzygies','cubic_closure']]
 assert D['schur_peer']['all_passed'] and D['cubic_report']['all_passed']
 assert S['actual_minimal_generator_presentation_shape']==[5,825]
 assert C['certified_cubic_rank']==3 and C['reduced_target_count']==5
 selected=Y['selected_actual_minimal_cycle_indices37'];assert len(selected)==37
 first=D['first640']['source_columns640_rows56_sparse_dual9']+D['native128']['source_columns128_rows48_sparse_dual9']
 X1=[{8*r+k-1:Q(c) for r,p in col for k,c in p if k and Q(c)} for col in first]
 old741=list(map(rec,S['source741_zero_constant_combinations']));old657=list(map(rec,S['remaining657_source741_combinations']))
 combos=[];xx=[]
 for j in selected:
  combo={}
  for k,c in old657[j].items():combo=add(combo,old741[k],c)
  x={}
  for k,c in combo.items():x=add(x,X1[k],c)
  combos.append(combo);xx.append(x)
 m0=[list(map(Q,row)) for row in V['all24_by56_m0']]
 U=X['free_kernel_rows32'];pi=list(map(rec,S['constant_free32_quotient_projection5_rows']))
 for j,x in zip(selected,xx):
  assert all(sum(m0[i][r]*x.get(8*r+k,Q(0)) for r in range(56))==0 for i in range(24) for k in range(8))
  p={}
  for i,row in enumerate(pi):
   for k in range(8):
    a=sum(c*x.get(8*U[r]+k,Q(0)) for r,c in row.items())
    if a:p[8*i+k]=a
  assert p==rec(L['actual40_by657_linear_relation_columns'][j])
 mons2=list(itertools.combinations_with_replacement(range(8),2));mi={m:i for i,m in enumerate(mons2)}
 piv=[r for r in X['all768_fixed_elimination_receipts_and_constant_combinations'] if 'pivot' in r]
 pb={r['pivot']:rec(r['normalized']) for r in piv};pindex={r['pivot']:j for j,r in enumerate(piv)}
 xp0=list(map(rec,M['actual_pivot_X0_columns27_sparse56']));assert len(xp0)==27
 a2=[];known_x2=[]
 for r in Y['all119_linear_coefficient_syzygies']:
  syz=rec(r['coefficients296_sparse']);x2={}
  for index,c in syz.items():
   j,var=divmod(index,8)
   for rk,a in xx[j].items():
    row,k=divmod(rk,8);m=mi[tuple(sorted((var,k)))];key=36*row+m
    z=x2.get(key,Q(0))+c*a
    if z:x2[key]=z
    else:x2.pop(key,None)
  free={36*i+k:x2.get(36*row+k,Q(0)) for i,row in enumerate(U) for k in range(36) if x2.get(36*row+k,Q(0))};coeff={}
  for p in sorted(pb):
   for k in range(36):
    a=free.get(36*p+k,Q(0))
    if a:
     coeff[36*pindex[p]+k]=a
     free=add(free,{36*i+k:c for i,c in pb[p].items()},-a)
  assert not free,'linear syzygy x2 not in actual27-pivot origin span'
  check={}
  for jk,c in coeff.items():
   j,k=divmod(jk,36);check=add(check,{36*i+k:a for i,a in xp0[j].items()},c)
  assert check==x2,'full56 second-leading source span must hold'
  a2.append(enc(coeff));known_x2.append(enc(x2))
 assert len(a2)==119
 support=sum(map(len,combos))+sum(map(len,xx))+sum(map(len,a2))+sum(map(len,known_x2))+sum(map(len,M['actual_pivot_X1_by27_sparse56_times8']))+24*56
 assert support<150000
 ev={
  'actual_base_ring':'QQ[[zA0,...,zA6,rA]]=R0/(ell)',
  'actual_source_manifest':inputs,'actual_selected37_minimal_relation_indices':selected,
  'actual_selected37_original768_constant_combinations':list(map(enc,combos)),
  'actual_selected37_X1_sparse56_by8':list(map(enc,xx)),
  'all37_m0_Xj1_zero_and_pi_Xj1_equal_actual_Lj':True,
  'all119_linear_syzygies_source296':Y['all119_linear_coefficient_syzygies'],
  'all119_known_x2=sum_linear_l_times_Xj1_sparse56_times36':known_x2,
  'all119_actual_a2_sparse27_times36':a2,
  'all119_full56_x2_equals_Xp0_a2':True,
  'quadratic_monomials36':mons2,'cubic_monomials120':list(itertools.combinations_with_replacement(range(8),3)),
  'actual_Xp1_sparse27_columns':M['actual_pivot_X1_by27_sparse56_times8'],
  'typed_original_X2_callback_shape':[56,768,36],'actual_X2_callback_evaluated':False,
  'typed119_C3_formula':'sum_{j,var} l[j,var]*variable[var]*actual_X2(selected37_original_source_combo[j]) - sum_pivot actual_Xp1[pivot]*actual_a2[pivot]',
  'typed_cubic_two_channel_formula':'Pi_cubic*(pi*C3-Ystar1*U0^-1*kappa_degree3(m0*C3))',
  'Pi_cubic_source':'Compose all600 frozen cubic linear-rewrite normal forms with quotient by the frozen rank3 known-quadratic cubic-closure span.',
  'Pi_cubic_adjoint_extraction_owner':'a_zero_proof; no duplicate extraction or coefficient evaluation here',
  'Ystar1_U0_m0_pi':'same actual matrices pinned in minimal5 Schur and quadratic m1-cancellation packets',
  'degree3_kappa_operator':{'wi8':'(E+1)^-1 partial_ui(v8)','wij':'(E+2)^-1(partial_ui(vj)-partial_uj(vi))','E':'acts on OUTPUT polynomial after differentiation','applied_to':'original cubic leading coefficient, never an I²-normalform'},
  'unknowns_removed':['X3','m1','m2','Xp0*a3'],
  'exact_unknown_removal_reason':'z=sum lXj-Xp*a2 beginsorder3; (mz)3=m0*C3; pi*z3=pi*C3. m0Xp0=piXp0=0 removes the actual third pivot coefficient.',
  'Euler_commutator_scope':'leading wedge degree2 is in ker d2_linear, and H2 K(u7,0)=0 gives d3_linear times linear coefficients; its minimal cubic image is a known d3-quadratic times variable.',
  'literal_unprojected_cubic_equality_claimed':False,
  'only119_linear_syzygy_family_covered':True,
  'all_cubic_initial_relations_claimed_computed':False,
  'other_constant_quadratic_kernel_relations_may_need_X3':True,
  'native_effective_S2_scope':'Each selected37 source is already a zero-origin constant combination. The proved L2*S0 cancellation therefore applies before multiplying by l; a monomial-restricted B11 representative can only be used after its support is proved sufficient.',
  'no_new_source_second_or_third_coefficients_jet_solve_rank_GB_or_348_matrix':True,
  'expanded_QQ_support':support,'actual_HF3_total_T_length_or_nil_cube_computed':False}
 ep=out/'actual119-cubic-X2-only-source-operator-preflight.json';blob=json.dumps(ev,separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;ep.write_text(blob)
 report={'id':STEM,'all_passed':True,'actual_selected_source_relations':37,'actual_linear_syzygies':119,
  'all119_full56_quadratic_pivot_corrections_verified':True,'projected_cubic_X2_only_theorem_proved':True,
  'literal_unprojected_equality_or_complete_cubic_space_claimed':False,
  'new_source_coefficients_or_adjoints_evaluated':False,'T_length_HF3_depth_or_nil_cube_decided':False,
  'expanded_QQ_support':support,'bytes':len(blob.encode()),'input_hashes':inputs,
  'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 (out/'phase.json').write_text(json.dumps({'phase':'FINAL','report_sha256':sha(out/'report.json'),'new_X2_or_X3_values':False})+'\n')
 print(json.dumps(report),flush=True)
if __name__=='__main__':main()

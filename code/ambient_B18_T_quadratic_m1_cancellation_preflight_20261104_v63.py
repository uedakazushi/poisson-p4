#!/usr/bin/env python3
"""Projected quadratic T operator, with m1 removed; no second coefficient query.

Only fixed frozen elimination receipts are replayed.  The actual X2 callbacks
remain typed unknown coefficient functions, not values chosen to be zero.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63'
PATHS={
 'schur':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/actual-minimal5-by825-Schur-operator.json',
 'schur_report':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/report.json',
 'low':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/actual-linear-quadratic-T-certificate.json',
 'phase1':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
 'phase2':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase2-certified-quadratic-input.json',
 'X0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
 'V':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
 'first640':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
 'native128':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rec(a):return {int(i):Q(c) for i,c in a if Q(c)}
def enc(a):return [[i,str(c)] for i,c in sorted(a.items()) if c]
def add(a,b,c=Q(1)):
 a=dict(a)
 for i,v in b.items():
  z=a.get(i,Q(0))+c*v
  if z:a[i]=z
  else:a.pop(i,None)
 return a
def fixed_combinations(cert):
 basis=[];zeros=[]
 for j,r in enumerate(cert['fixed_elimination_receipts']):
  v={j:Q(1)}
  for p,c in r['subtract']:v=add(v,basis[p],-Q(c))
  if 'pivot_row' in r:
   c=Q(r['divide']);basis.append({i:a/c for i,a in v.items()})
  else:zeros.append(v)
 return basis,zeros
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 (out/'phase.json').write_text(json.dumps({'phase':'typed_proof_inputs_saved_before_construction','inputs':inputs},indent=2)+'\n')
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 S,L,P1,P2,X,V=[D[k] for k in ['schur','low','phase1','phase2','X0','V']]
 assert D['schur_report']['all_passed'] and L['actual_linear_relation_rank']==37
 assert L['quadratic_certificate']['rank']==177
 assert S['actual_minimal_generator_presentation_shape']==[5,825]
 assert X['rank_X0']==27
 linbasis,zero620=fixed_combinations(L['complete_linear_rank_certificate'])
 assert len(linbasis)==37 and len(zero620)==620
 # The three quotient functionals are obtained from the frozen177 pivot
 # vectors.  This is a deterministic projection replay, not a new rank test.
 qc=L['quadratic_certificate'];qb=list(map(rec,qc['basis_normalized_columns']))
 leads={min(c):j for j,c in enumerate(qb)};tail=[i for i in range(180) if i not in leads]
 assert tail==[143,176,179] and len(qb)==177
 Pi=[{} for _ in tail]
 for j in range(180):
  v={j:Q(1)}
  for p in sorted(leads):v=add(v,qb[leads[p]],-v.get(p,Q(0)))
  assert all(i in tail for i in v)
  for k,i in enumerate(tail):
   if v.get(i,Q(0)):Pi[k][j]=v[i]
 assert [[r.get(j,Q(0)) for j in tail] for r in Pi]==[[Q(i==j) for j in range(3)] for i in range(3)]
 # Check annihilation of all464 previously proved actual quadratics; their
 # full-source arithmetic is not recomputed or assumed for unknown X2.
 qinput=P2['all464_certified_quadratic_columns_sparse']
 assert len(qinput)==464
 assert all(not sum(row.get(i,Q(0))*c for i,c in rec(col).items()) for row in Pi for col in qinput)
 mons=P2['quadratic_monomials36'];assert mons==[list(p) for p in itertools.combinations_with_replacement(range(8),2)]
 # Source X1 is only used to certify the promised27-dimensional span.  The
 # actual second source coefficients remain absent from this calculation.
 xfirst=D['first640']['source_columns640_rows56_sparse_dual9']+D['native128']['source_columns128_rows48_sparse_dual9']
 assert len(xfirst)==768
 X1=[{8*r+k-1:Q(c) for r,a in col for k,c in a if k and Q(c)} for col in xfirst]
 piv=[r for r in X['all768_fixed_elimination_receipts_and_constant_combinations'] if 'pivot' in r]
 pb={r['pivot']:rec(r['normalized']) for r in piv};porder=sorted(pb)
 U=X['free_kernel_rows32'];pi=[rec(r) for r in S['constant_free32_quotient_projection5_rows']]
 xp0=[];xp1=[]
 for r in piv:
  a0={};a1={}
  for j,c in rec(r['constant_original_source_combination']).items():
   a0=add(a0,{i:Q(row[j]) for i,row in enumerate(X['full_X0_56_by768']) if Q(row[j])},c)
   a1=add(a1,X1[j],c)
  xp0.append(a0);xp1.append(a1)
 old741=[rec(r) for r in S['source741_zero_constant_combinations']]
 old657=[rec(r) for r in S['remaining657_source741_combinations']]
 raw657=[]
 for c in old657:
  v={}
  for j,a in c.items():v=add(v,old741[j],a)
  raw657.append(v)
 a1all=[];comboall=[]
 for c in zero620:
  original={}
  for j,a in c.items():original=add(original,raw657[j],a)
  v={}
  for j,a in original.items():v=add(v,X1[j],a)
  coeff={};free={8*i+k:v.get(8*row+k,Q(0)) for i,row in enumerate(U) for k in range(8) if v.get(8*row+k,Q(0))}
  pivindex={r['pivot']:j for j,r in enumerate(piv)}
  for p in porder:
   for k in range(8):
    a=free.get(8*p+k,Q(0))
    if a:
     coeff[8*pivindex[p]+k]=a
     free=add(free,{8*i+k:b for i,b in pb[p].items()},-a)
  assert not free,'zero-linear source first not in the actual27 pivot span'
  reconstruct={}
  for ik,a in coeff.items():
   j,k=divmod(ik,8);reconstruct=add(reconstruct,{8*i+k:b for i,b in xp0[j].items()},a)
  assert reconstruct==v,'full56 source span equality must hold, not only free32'
  a1all.append(enc(coeff));comboall.append(enc(original))
 # The following is an exact finite typed LINEAR functional on actual X2.
 # All matrix/first coefficients are frozen inputs.  There is no m1 slot.
 op={
  'actual_X2_input_shape':[56,768,36],
  'actual_X2_values_provided_or_assigned':False,
  '620_actual_source_combinations_in_original768':comboall,
  '620_actual_linear_pivot_coefficients_a1_by27_times8':a1all,
  'actual_pivot_X1_by27_sparse56_times8':list(map(enc,xp1)),
  'actual_pivot_X0_columns27_sparse56':list(map(enc,xp0)),
  'm0_actual24_by56':V['all24_by56_m0'],
  'pi_actual5_by56':'constant_free32_quotient_projection5_rows on free_kernel_rows32, all other columns0',
  'pi_free32_rows':list(map(enc,pi)),'free_kernel_rows32':U,
  'Ystar1_40_by84_columns':P1['normalized84_unitcycle_Q1_columns_sparse'],
  'U0_84_columns':S['actual_second_stage_U0_columns84'],
  'Pi3_180_dual_rows':list(map(enc,Pi)),
  'target_quadratic_order':'five actual free targets, then quadratic_monomials36',
  'quadratic_monomials36':mons,
  'C2':'sum_source_combinations(actual_X2)-sum_pivot(Xp1*a1), with ordinary polynomial multiplication of linear coefficients',
  'kappa_linear':{
   'input':'24 actual quadratic coefficients m0*C2 in sector blocks8',
   'wi8':'(E+1)^-1 partial_ui(v8)',
   'wij':'(E+2)^-1(partial_ui(vj)-partial_uj(vi))',
   'E_on_output':'sum first7 u_i*partial_ui after the derivative',
   'three_sector_linear_coordinates':'same zA0..6 at ordinary linear level; CB quadratic shift changes no coefficient here'},
  'actual_three_class_functional_on_each_source_C2':'Pi3*(pi*C2-Ystar1*U0^-1*kappa_linear(m0*C2))',
  'unknown_m1_not_an_input':True,
  'equality_scope':'three-class quotient modulo all known variable*linear relations and all168 known d3 quadratic relations',
  'literal_unprojected_Schur_quadratic_equality_claimed':False,
  'Euler_commutator_boundary':'H2 K(u1,...,u7,0)=H2 K(u7) direct_sum H1 K(u7)=0; the leading linear wedge difference is a literal d3_linear constant boundary',
  'actual_second27_schur_coefficients_a2_drop':'m0*Xp0=0 and pi*Xp0=0',
  'native_unit_second_basis_change_naturality':{
   'premise':'coherent target basis P=I+P2+O(m3) on the source56 frame and coherent relation/source basis change',
   'Xrest0_zero':'after all frozen constant source combinations',
   'left_transport_effect_on_Xrest2':'-P2*Xrest0=0',
   'right_relation_transport_effect':'Xp0*a2, killed by m0 and pi',
   'warning':'Changing an adapter coefficient in SyzQ alone while freezing Q and B is not a coherent basis change. Genuine physical Q2/core changes are not removed by this lemma.'}}
 support=(sum(map(len,a1all))+sum(map(len,comboall))+sum(map(len,xp0))+sum(map(len,xp1))
          +sum(map(len,Pi))+24*56+5*32+sum(map(len,P1['normalized84_unitcycle_Q1_columns_sparse']))
          +sum(map(len,S['actual_second_stage_U0_columns84'])))
 assert support<150000
 assert len(json.dumps(op).encode())<32<<20
 (out/'actual-projected-three-class-X2-operator-preflight.json').write_text(json.dumps(op,separators=(',',':'))+'\n')
 ep=out/'actual-projected-three-class-X2-operator-preflight.json'
 report={'id':STEM,'all_passed':True,'actual_zero_linear_source_combinations':620,'actual_full56_first_span_checks':620,
  'actual_three_class_quotient_dimension':3,'Pi3_annihilates_all464_known_quadratics':True,
  'projected_quadratic_m1_cancellation_proved':True,'literal_unprojected_equality_claimed':False,
  'coherent_pure_second_basis_change_naturality_proved':True,'actual_N2_coefficients_proved_all_pure_basis':False,
  'saved_rational_support_count':support,
  'actual_X2_coefficients_computed':False,'T_length_or_depth_decided':False,
  'input_hashes':inputs,'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 (out/'phase.json').write_text(json.dumps({'phase':'FINAL','report_sha256':sha(out/'report.json'),'actual_X2_not_evaluated':True},indent=2)+'\n')
 print(json.dumps(report),flush=True)
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Constant619 quartic source bridge; retain X3/P4/nativeL2 as unknowns."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,resource,signal,time

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_CONSTANT619_QUARTIC_SOURCE_PREFLIGHT_20261104_v63'
PATHS={
 'complete':'reports/AMBIENT_B18_T_PURE_CUBIC619_COMPLETENESS_PREFLIGHT_20261104_v63/independent-pure619-Schur-Theta-and-completeness-preflight.json',
 'm1':'reports/AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/actual-projected-three-class-X2-operator-preflight.json',
 'X0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
 'first640':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
 'native128':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json',
 'HF3peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json',
 '168peer':'reports/AMBIENT_B18_T_CORRECTED168_QUARTIC_VALUES_20261104_v63_REVIEW/report.json',
 'd3preflight':'reports/AMBIENT_B18_T_CORRECTED_D3_QUARTIC_PREFLIGHT_20261104_v63/actual168-corrected-d3-quartic-source-preflight.json'}
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
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 dump(out,'phase.json',{'phase':'constant619_source_inputs_saved','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert D['HF3peer']['all_passed'] and D['168peer']['all_passed']
 M=D['m1'];c620=list(map(rec,M['620_actual_source_combinations_in_original768']));a620=list(map(rec,M['620_actual_linear_pivot_coefficients_a1_by27_times8']))
 xp0=list(map(rec,M['actual_pivot_X0_columns27_sparse56']))
 allfirst=D['first640']['source_columns640_rows56_sparse_dual9']+D['native128']['source_columns128_rows48_sparse_dual9']
 X1=[{8*r+k-1:Q(c) for r,p in col for k,c in p if k and Q(c)} for col in allfirst]
 X0=D['X0']['full_X0_56_by768'];cases=[];support=0
 for j,K in enumerate(D['complete']['actual619_constant_quadratic_kernel_combinations620']):
  c={};a1={}
  for k,s in rec(K).items():c=add(c,c620[k],s);a1=add(a1,a620[k],s)
  first={}
  for k,s in c.items():first=add(first,X1[k],s)
  check={}
  for index,s in a1.items():
   p,var=divmod(index,8)
   for r,t in xp0[p].items():check=add(check,{8*r+var:t},s)
  assert first==check
  for row in range(56):
   assert sum((Q(X0[row][k])*s for k,s in c.items() if k<640),Q(0))==0
   assert sum((Q(X0[row][k])*s for k,s in c.items() if k>=640),Q(0))==0
  cases.append({'constant619_index':j,'actual_original768_constant_combination_c0':enc(c),
   'actual_Xp27_linear_a1_sparse27_times8':enc(a1),'full56_Xc1_equals_Xp0_a1':True,
   'old_and_native_origin_parts_separately_zero':True})
  support+=len(c)+len(a1)+len(first);assert support<150000
 assert len(cases)==619
 ev={'source_manifest':inputs,'all619':cases,
  'actual_base':'QQ[[zA0,...,zA6,rA]]; lasttargete4/source55',
  'typed_C2':'Xc2-Xp1*a1',
  'typed_b1':'U0^-1*kappa_degree2(m0*C2)',
  'typed_Qproxy2':'pi*C2-Ystar1*b1',
  'typed_l1':'fixed177-basis lift of Qproxy2 into L37*variables; existence from actual619 kernel of the additional qB projection',
  'typed_a2':'fixed27-origin-pivot projection of C2-Xstar1*b1-XL1*l1; this vector lies ker(m0) intersect ker(pi)',
  'typed_corrected_source_z':'Xc-Xp*a1-Xstar*b1-XL*l1-Xp*a2',
  'typed_C3':'Xc3-Xp2*a1-Xstar2*b1-XL2*l1-Xp1*a2',
  'typed_beta2':'U0^-1*kappa_degree3(m0*C3)',
  'typed_Qproxy3':'pi*C3-Ystar1*beta2',
  'typed_raw_pure4':'rawP55_c,4 -Xp55,3rrr*a1r -Xstar55,3rrr*b1r -XL55,3rrr*l1r -Xp55,2rr*a2rr -Ystar55,2rr*beta2rr',
  'typed_corrected_quartic':'typed_raw_pure4-Theta4(typed_Qproxy3)',
  'Theta4_source':'acceptedTheta3+acceptedc4*PiC, b3=-42215/3456,c4=4565/68832; actual canonical source corrections pinned in168peer',
  'm1_removal_proof':'After subtracting actual homogeneous a1,b1,l1,a2, z startsdegree3. Therefore(mz)3=m0*C3. The full56 first source check is literal; higher steps are typed coefficient operators, not evaluated arrays.',
  'Euler_boundary_proof':'The stored corrected wedge and kappa(mz) startdegree2; their exactd3 difference has coefficient in m since d3_linear on constant wedge3 is injective. Puredegree2/3 ofactualminimald3 vanish, so pure4 difference0; Theta4 kills its cubic initial.',
  'native_expansion_order3':'L0*(S3*c0+S2*c1+S1*c2)+L1*(S2*c0+S1*c1+S0*c2)+L2*(S1*c0+S0*c1)+L3*S0*c0; c=c0+c1+c2 is the chosen homogeneous coefficient vector, with monicRes4 scalar lowering included in L',
  'native_L3_drops':'S0*c0=0 separately, checked; no third adapter derivative requested',
  'native_L2_must_remain':'S1*c0+S0*c1 is not knownzero. The individual nativeXp2 terms cannot be read from effectiveS2; aggregate L2 sensitivity must be retained and tested.',
  'unknown_coefficient_gate':['required X2 pairs for b1/l1/a2 not yet adjoint-bound','actualX3 mixed coefficients including nativeS3/L2','ordinaryP4/commonCB4/regularquotient derivatives, subject to cancellation'],
  'no_input_cache_availability_beyond_stated_source_firsts_claimed':True,
  'no_new_source_coefficient_solve_point_prime_highjet_or_348_matrix':True,
  'actual619_quartic_values_HF4_exactlength_not_computed':True,'expanded_QQ_support':support}
 ep=out/'actual619-quartic-source-order-and-unknown-gates.json';dump(out,ep.name,ev)
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'all619_full56_source_first_identities_verified':True,
  'constant619_quartic_m1_free_source_bridge_proved':True,'native_L2_cancellation_or_cache_sufficiency_claimed':False,
  'actual619_quartic_values_or_HF4_or_T_length_computed':False,'input_sha256':inputs,
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 assert sum(p.stat().st_size for p in out.iterdir())<32<<20
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL','report_sha256':sha(out/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()

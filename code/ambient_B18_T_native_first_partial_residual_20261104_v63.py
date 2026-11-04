#!/usr/bin/env python3
"""Actual619 old/native first residuals: full equality need not split."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_NATIVE_FIRST_PARTIAL_RESIDUAL_20261104_v63'
PATHS={
 '619':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_SOURCE_PREFLIGHT_20261104_v63/actual619-quartic-source-order-and-unknown-gates.json',
 '619report':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_SOURCE_PREFLIGHT_20261104_v63/report.json',
 'S':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/actual-minimal5-by825-Schur-operator.json',
 'X0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
 'oldX1':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
 'nativeX1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rec(a):return {int(i):Q(c) for i,c in a if Q(c)}
def enc(a):return [[i,str(c)] for i,c in sorted(a.items()) if c]
def put(a,i,c):
 v=a.get(i,Q(0))+c
 if v:a[i]=v
 else:a.pop(i,None)
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 dump(out,'phase.json',{'phase':'actual_partial_source_array_inputs_saved','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()};assert D['619report']['all_passed']
 X0=D['X0']['full_X0_56_by768'];xp=[[],[]]
 for c in D['S']['first27_source_combinations']:
  cp=rec(c)
  for sector in range(2):
   q={}
   for r in range(56):
    a=sum((Q(X0[r][j])*s for j,s in cp.items() if (j>=640)==bool(sector)),Q(0))
    if a:q[r]=a
   xp[sector].append(q)
 first=D['oldX1']['source_columns640_rows56_sparse_dual9']+D['nativeX1']['source_columns128_rows48_sparse_dual9']
 X1=[{8*r+k-1:Q(c) for r,p in col for k,c in p if k and Q(c)} for col in first]
 allres=[];nonzero=[];support=0
 for row in D['619']['all619']:
  c=rec(row['actual_original768_constant_combination_c0']);a1=rec(row['actual_Xp27_linear_a1_sparse27_times8']);res=[{},{}]
  for j,s in c.items():
   sector=int(j>=640)
   for i,a in X1[j].items():put(res[sector],i,s*a)
  for i,s in a1.items():
   p,v=divmod(i,8)
   for sector in range(2):
    for r,a in xp[sector][p].items():put(res[sector],8*r+v,-s*a)
  assert res[0]=={i:-c for i,c in res[1].items()}
  j=row['constant619_index']
  if res[1]:nonzero.append(j)
  support+=len(res[1]);allres.append({'constant619_index':j,'native_partial_first_residual56_times8':enc(res[1]),'old_partial_is_exact_negative':True})
 assert len(allres)==619 and len(nonzero)==82 and support==634
 assert allres[288]['native_partial_first_residual56_times8']==[[41,'-1'],[49,'1/4'],[257,'-1']]
 ev={'source_manifest':inputs,'source_order':'original768 = old0..639 and native640..767',
  'external_variables':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA'],
  'slot_definition':'8*actual_fullV56_row+variable_index; native last8 rows48..55 are zero',
  'Xp_frame':'actual first27_source_combinations from frozen minimal5 Schur, applied to same fullX0 original768 matrix',
  'definition':'sum native X1_j*c0_j - sum native Xp0_p*a1_p; old residual formed separately on old0..639',
  'all619':allres,'nonzero_native_partial_indices':nonzero,'native_partial_QQ_support':support,
  'full56_identity_verified_by_old_plus_native_zero':True,
  'Omega1_source_meaning':'Omega1=L0^-1*native_partial_residual in the actual scalarRes4C4 frame; L0 is a unit, hence these82 are nonzero proper residuals as well',
  'N2_scope':'This only disproves automatic zero from the fullfirstidentity. It does not assert Gamma4*(L2*Omega1) nonzero; that exact adjoint fold remains necessary.',
  'no_new_source_coefficients_N2_values_jet_or_solve':True}
 ep=out/'actual619-native-old-first-partial-residuals.json';dump(out,ep.name,ev)
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'actual_nonzero_native_partial_sources':82,'actual_native_partial_QQ_support':634,
  'all_old_partial_residuals_exact_negative':True,'Gamma4_L2_cancellation_or_new_source_values_claimed':False,
  'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL','report_sha256':sha(out/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()

#!/usr/bin/env python3
"""ONE guarded known-QQ adjoint test; no physical or unknown-source values."""
import argparse,json,os,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_MIXED619_LOWER2_ADJOINT_20261105_v63'
class GuardStop(Exception):pass
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 inputs={};rawbytes=0;inputnum=0
 def read(p):
  nonlocal rawbytes,inputnum
  f=ROOT/p;rawbytes+=f.stat().st_size;inputs[p]={'path':p,'sha256':sha(f)};v=json.loads(f.read_text());inputnum+=count(v);return v
 gatep='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER2_MIXED_ADJOINT_INPUT_20261105_v63/report.json';gate=read(gatep);assert gate['all_passed']
 dp='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER2_MIXED_ADJOINT_INPUT_20261105_v63/actual136-mixed-output210-slot-adjoint-source-input.json';d=read(dp);assert sha(ROOT/dp)==gate['evidence_sha256'][Path(dp).name]
 gp='reports/AMBIENT_B18_ELL2_T2_MIXED619_ADJOINT_DENSE_GUARD_20261105_v63/report.json';g=read(gp);assert g['all_passed'] and g['proven_whole_numeric']==51562 and g['proven_whole_RAW']==16842088
 cp='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_COVECTORS_20261105_v63/actual-four-full918-source-covectors-and52-plane-adjoint-input.json';full=read(cp)
 assert inputnum<=g['source_metadata_numeric_cap_before_arithmetic']==6000 and rawbytes<=g['source_RAW_input_cap_before_arithmetic']==1000000
 for p,h in gate['input_sha256'].items():assert sha(ROOT/h['path'])==h['sha256']
 A=list(map(sp,d['matrix_columns_on210_domain']));targets=list(map(sp,d['target_covectors2_on210_domain']));domain=d['full918_domain_slots210'];n,m=len(A),len(domain);assert(n,m)==(91,210)
 assert d['all_other708_domain_slots_identically_zero'] and len(d['nonzero_output_indices'])==91
 norm={};reps={};selected=[];history=[];solutions=[];current_v=[];current_rep=[];pending={};phase='source_ready_before_elimination';column=-1;targetindex=-1;maxcommitted=0;maxtemporary=0
 hp=out/'pivot-history.jsonl';hp.write_text('')
 def dump_frontier(reason=None):
  temp=out/'frontier.tmp';final=out/'latest-frontier.json'
  with temp.open('w') as f:
   f.write('{"phase":'+json.dumps(phase)+',"source_column":'+str(column)+',"target_index":'+str(targetindex)+',"selected_original_columns":'+json.dumps(selected)+',"basis":[')
   for k,p in enumerate(sorted(norm)):
    if k:f.write(',')
    f.write(json.dumps({'pivot':p,'vector':[str(x) for x in norm[p]],'original_source_representation':[str(x) for x in reps[p]]},separators=(',',':')))
   f.write('],"current_vector":'+json.dumps([str(x) for x in current_v])+',"current_representation":'+json.dumps([str(x) for x in current_rep])+',"pending_operation":'+json.dumps(pending)+',"reason":'+json.dumps(reason)+',"partial_identity_semantics":'+json.dumps('basis[p]=A*representation[p]; current input residual=A*current_representation, or current target residual=target-A*current_representation')+'}\n')
  os.replace(temp,final)
  assert sum(p.stat().st_size for p in out.iterdir())<g['unchanged_RAW_cap']==16940016
 def alarm(_s,_f):raise GuardStop('60s time guard')
 signal.signal(signal.SIGALRM,alarm)
 def tempcheck(x):
  nonlocal maxtemporary
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length());maxtemporary=max(maxtemporary,b)
  if b>1153:raise GuardStop('uncommitted scalar guard')
  return x
 def commit(x):
  nonlocal maxcommitted
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length());maxcommitted=max(maxcommitted,b)
  if b>384:raise GuardStop('committed coefficient guard')
  return x
 def axpy(v,w,c):
  # Construct a separate vector; old complete residual remains intact if a
  # single new coefficient fails its guard.
  return [commit(tempcheck(x+tempcheck(c*y))) for x,y in zip(v,w,strict=True)]
 def scale(v,c):return [commit(tempcheck(x/c)) for x in v]
 def dot(v,w):
  a=Q(0)
  for x,y in zip(v,w,strict=True):a=commit(tempcheck(a+tempcheck(x*y)))
  return a
 def first(v):return next((i for i,x in enumerate(v) if x),None)
 def original_image(r):
  v=[Q(0)]*m
  for j,c in enumerate(r):
   if c:
    for i,x in A[j].items():v[i]=commit(tempcheck(v[i]+tempcheck(c*x)))
  return v
 try:
  phase='knownQQ_column_elimination';dump_frontier()
  for column,col in enumerate(A):
   current_v=[col.get(i,Q(0)) for i in range(m)];current_rep=[Q(int(j==column)) for j in range(n)];steps=[]
   while (p:=first(current_v)) in norm:
    c=current_v[p];pending={'kind':'subtract_basis','pivot':p,'coefficient':str(c)}
    v=axpy(current_v,norm[p],-c);a=axpy(current_rep,reps[p],-c);current_v,current_rep=v,a;steps.append([p,str(c)])
   p=first(current_v);divide=Q(1)
   if p is not None:
    divide=current_v[p];pending={'kind':'normalize_new_basis','pivot':p,'divide_by':str(divide)}
    v=scale(current_v,divide);a=scale(current_rep,divide);assert v[p]==1
    norm[p]=v;reps[p]=a;selected.append(column)
   receipt={'source_column':column,'pivot':p,'divide_by':str(divide),'subtract':steps};history.append(receipt)
   with hp.open('a') as f:f.write(json.dumps(receipt,separators=(',',':'))+'\n')
   assert original_image(current_rep)==current_v
   pending={};dump_frontier()
  # All input columns are now in the fixed basis. This is not an independent
  # mathematical rank assertion; only the following adjoint identities matter.
  phase='reduce_two_target_covectors'
  for targetindex,target in enumerate(targets):
   current_v=[target.get(i,Q(0)) for i in range(m)];current_rep=[Q(0)]*n
   for p in sorted(norm):
    c=current_v[p]
    if not c:continue
    pending={'kind':'target_subtract_basis','pivot':p,'coefficient':str(c)}
    v=axpy(current_v,norm[p],-c);a=axpy(current_rep,reps[p],c);current_v,current_rep=v,a
   if (k:=first(current_v)) is not None:
    # Exact dual of the fixed reduction, in reverse operator order.
    theta=[Q(int(i==k)) for i in range(m)]
    for p in reversed(sorted(norm)):
     c=dot(theta,norm[p])
     theta[p]=commit(theta[p]-c)
    assert all(dot([x for i,x in col.items()],[theta[i] for i,x in col.items()])==0 for col in A)
    pairing=dot([x for i,x in target.items()],[theta[i] for i,x in target.items()]);assert pairing
    phase='restricted_output_family_inconsistent';pending={};dump_frontier()
    result={'both_adjoint_identities_exist':False,'failed_target':targetindex,'exact_left_null_on210':[[i,str(x)] for i,x in enumerate(theta) if x],'nonzero_target_pairing':str(pairing),'scope':'Only136chosenquadraticoutputs. No assertionaboutall619orfullsourcecoefficientnecessity.'}
    break
   assert original_image(current_rep)==[target.get(i,Q(0)) for i in range(m)]
   solutions.append(current_rep);pending={};dump_frontier()
  else:
   alpha=list(map(sp,full['full918_covectors4']))[2:];slot={v:i for i,v in enumerate(domain)}
   for target,rho in zip(alpha,solutions,strict=True):
    for i in range(918):
     v=Q(0)
     if i in slot:
      for j,c in enumerate(rho):
       if c:v=commit(tempcheck(v+tempcheck(c*A[j].get(slot[i],Q(0)))))
     assert v==target.get(i,Q(0))
   result={'both_adjoint_identities_exist':True,'solutions_on91_nonzero_output_columns':[[[j,str(x)] for j,x in enumerate(v) if x] for v in solutions],
    'solutions_on136_literal_outputs':[[[d['nonzero_output_indices'][j],str(x)] for j,x in enumerate(v) if x] for v in solutions],
    'all2_times918_literal_identities_verified':True,'output_frame17_by8':d['monomial_output_order8'],'no_unknown_X2_or_619_values':True}
   phase='complete_known_adjoint_test';dump_frontier()
  ep=out/'actual-two-mixed-output-adjoint-certificates-or-obstruction.json';ep.write_text(json.dumps(result,separators=(',',':'))+'\n')
  report={'id':STEM,'all_passed':True,'both_adjoint_identities_exist':result['both_adjoint_identities_exist'],'no_physical_gauge_or_619_values':True,'dynamic_guards_completed':True,'actual_maximum_committed_bits':maxcommitted,'actual_maximum_uncommitted_bits':maxtemporary,
   'whole_numeric_bound':g['proven_whole_numeric'],'whole_RAW_bound':g['proven_whole_RAW'],'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
   'evidence_sha256':{ep.name:sha(ep),hp.name:sha(hp),'latest-frontier.json':sha(out/'latest-frontier.json')},'seconds':time.monotonic()-start}
  (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'both_adjoint_identities_exist':result['both_adjoint_identities_exist'],'actual_bits':maxcommitted,'report_sha256':sha(out/'report.json')}))
 except GuardStop as e:
  signal.alarm(0);dump_frontier(str(e));(out/'guard-stop.json').write_text(json.dumps({'stopped':str(e),'no_autoretry':True,'code_sha256':sha(Path(__file__)),'latest_frontier_sha256':sha(out/'latest-frontier.json'),'pivot_history_sha256':sha(hp)},separators=(',',':'))+'\n');print(json.dumps({'guard_stop':str(e),'frontier_saved':True}));raise SystemExit(2)
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);a=parser.parse_args();main(a.output.resolve())

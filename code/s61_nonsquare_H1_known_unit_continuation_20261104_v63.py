#!/usr/bin/env python3
"""ONE bounded continuation: reuse96 exact H1 rows and seven unit pivots; query only38 new rows."""
import argparse,gzip,hashlib,json,os,resource,signal,time
from collections import Counter
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_H1_KNOWN_UNIT_CONTINUATION_20261104_v63'
class Stop(Exception):pass
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
 signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(Stop('57s_internal_boundary')));signal.alarm(57)
 parent=ROOT/'reports/S61_NONSQUARE_H1_KNOWN_UNIT_ELIMINATION_20261104_v63';pr=read(parent/'report.json');pp=parent/'checkpoint.json.gz';assert sha(pp)==pr['evidence_sha256'][pp.name]
 saved=json.loads(gzip.decompress(pp.read_bytes()));assert saved['completed_pivots']==7 and len(saved['materialized_raw_rows'])==96
 base=ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63';ip=base/'certificate.json';assert sha(ip)==saved['source_input_sha256'];source=read(ip)
 nodes=source['scalar_circuit']['nodes'];atoms=source['scalar_circuit']['atoms_num_den'];roots=source['complete_H1_image_block_scalar_roots']
 R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();known=[a,d,500*a+90*e**2+3*e-24,(50*a+18*e**2-9*e-2)*d+60*e-40,100*a+18*e**2-9*e-2,5*a*d-6*e-4]
 decp=lambda p:R.from_dict({tuple(m):s.QQ(c)for m,c in p});dec=lambda c:K(decp(c[0]))/K(decp(c[1]))
 encp=lambda p:[[list(m),str(c)]for m,c in p.items()];enc=lambda c:[encp(c.numer),encp(c.denom)]
 materialized={i:[dec(c)for c in row]for i,row in saved['materialized_raw_rows']};initial_ids=set(materialized);unexpanded=[i for i in range(134)if i not in initial_ids];assert len(unexpanded)==38
 assert sum(i<113 for i in unexpanded)==17 and sum(i>=113 for i in unexpanded)==21
 steps=[]
 for st in saved['unit_pivots']:
  st=dict(st);st['_row']=[dec(c)for c in st['normalized23']];steps.append(st)
 initial_pivots=len(steps);available=saved['remaining_rows'][:];cols=saved['remaining_columns'][:];tested=saved['tested_row_pivotstage_pairs'][:];rowstatus=dict(saved['row_status'])
 new_completed=[];new_residual={};new_initial_coeffs={};evaluated_node_count=0;max_live_terms=0;cache={0:K.zero};atomcache={};queries=[];stop=None;pending=None
 def guard(c):
  if max(len(c.numer),len(c.denom))>100000:raise Stop('100000_output_terms_per_entry')
  if time.monotonic()-began>53:raise Stop('53s_arithmetic_boundary')
  return c
 def op(tag,x,y):
  if tag=='m':
   if not x or not y:return K.zero
   bound=max(len(x.numer)*len(y.numer),len(x.denom)*len(y.denom))
  else:
   if not y:return x
   if not x:return y if tag=='a'else-y
   bound=len(x.numer)+len(y.numer)if x.denom==y.denom else max(len(x.numer)*len(y.denom)+len(y.numer)*len(x.denom),len(x.denom)*len(y.denom))
  if bound>150000:raise Stop('150000_intermediate_term_bound')
  return guard(x*y if tag=='m'else x+y if tag=='a'else x-y)
 def evaluate(i):
  nonlocal evaluated_node_count,max_live_terms
  assert i in unexpanded and i not in initial_ids and i not in materialized;queries.append(i);wanted=roots[i]
  need=set();todo=wanted[:]
  while todo:
   n=todo.pop()
   if n in cache or n in need:continue
   need.add(n);nd=nodes[n]
   if nd[0]!='c':assert nd[0]in['a','s','m'];todo.extend(nd[1:])
  refs=Counter(wanted)
  for n in need:
   if nodes[n][0]!='c':refs.update(nodes[n][1:])
  live={};live_terms=0
  for n in sorted(need):
   nd=nodes[n]
   if nd[0]=='c':
    ai=nd[1]
    if ai not in atomcache:atomcache[ai]=dec(atoms[ai])
    v=atomcache[ai]
   else:
    x=cache[nd[1]]if nd[1]in cache else live[nd[1]];y=cache[nd[2]]if nd[2]in cache else live[nd[2]];v=op(nd[0],x,y)
   live[n]=v;evaluated_node_count+=1;live_terms+=len(v.numer)+len(v.denom);max_live_terms=max(max_live_terms,live_terms)
   if nd[0]!='c':
    for p in nd[1:]:
     refs[p]-=1
     if not refs[p]and p in live:
      z=live.pop(p);live_terms-=len(z.numer)+len(z.denom)
  ans=[cache[n]if n in cache else live[n]for n in wanted]
  for z in ans:guard(z)
  return ans
 def reduce_row(raw,start=0,cofs=None):
  row=raw[:];cofs=[]if cofs is None else cofs[:]
  for st in steps[start:]:
   c=row[st['pivot_column']];cofs.append(c)
   if c:row=[op('s',z,op('m',c,v))for z,v in zip(row,st['_row'])]
   assert not row[st['pivot_column']]
  return row,cofs
 def checkpoint():
  state={'source_input_sha256':sha(ip),'previous_checkpoint_sha256':sha(pp),'coefficient_domain':saved['coefficient_domain'],'shape_redundant_H1':[134,23],
   'completed_pivots':len(steps),'initial_reused_pivots':initial_pivots,'remaining_rows':available,'remaining_columns':cols,
   'materialized_raw_rows':[[i,[enc(c)for c in row]]for i,row in materialized.items()],
   'unit_pivots':[{k:v for k,v in st.items()if k!='_row'}for st in steps],
   'raw_row_preconditioner':saved['raw_row_preconditioner'],'completed_evaluated_rows':list(materialized),'row_status':rowstatus,'tested_row_pivotstage_pairs':tested,
   'initial96_row_ids_reused_without_scalar_queries':sorted(initial_ids),'new38_required_row_ids':unexpanded,'new_completed_raw_and_initial7_residual_row_ids':new_completed,
   'new_row_residuals_after_initial7':[[i,[enc(c)for c in new_residual[i]],[enc(c)for c in new_initial_coeffs[i]]]for i in new_completed],
   'exact_scalar_query_row_ids':queries,'no_old96_scalar_re_evaluation':not any(i in initial_ids for i in queries),
   'all_unexpanded_rows_retained_as_original_scalar_circuit_plus_saved_source_preconditioner':True,
   'redundant_full_source_target_after_pivots':[190-len(steps),121-len(steps)],'no_literal_free163_minus_pivots_claim':True,
   'lazy_full_source_pivot_lift':saved['lazy_full_source_pivot_lift'],'pure98_coefficient_queries':0,
   'evaluated_scalar_nodes_with_recomputation':evaluated_node_count,'max_live_terms':max_live_terms,'seconds':time.monotonic()-began}
  data=gzip.compress((json.dumps(state,separators=(',',':'))+'\n').encode(),compresslevel=3,mtime=0)
  if len(data)>30*1024*1024:raise Stop('30MiB_checkpoint_guard')
  tmp=out/'checkpoint.tmp';tmp.write_bytes(data);os.replace(tmp,out/'checkpoint.json.gz')
  (out/'progress.json').write_text(json.dumps({'completed_pivots':len(steps),'new_completed_rows':new_completed,'all134_completed':len(materialized)==134,'seconds':time.monotonic()-began},indent=2)+'\n')
 unitcache={}
 def unit(c):
  if not c:return None
  key=(tuple(c.numer.items()),tuple(c.denom.items()))
  if key in unitcache:return unitcache[key]
  ex=[];ground=[]
  for p in[c.numer,c.denom]:
   powers=[]
   for f in known:
    n=0
    while p and not p.rem(f):p=p.exquo(f);n+=1
    powers.append(n)
   if not p.is_ground:unitcache[key]=None;return None
   ex.append(powers);ground.append(p)
  answer={'constant':str(K(ground[0])/K(ground[1])),'a_d_F_B_C_J_powers':[x-y for x,y in zip(*ex)]};unitcache[key]=answer;return answer
 checkpoint()
 try:
  # No completed old row is queried: first finish the38 fresh exact coefficients, saving every completed row.
  for i in unexpanded:
   if time.monotonic()-began>48:raise Stop('48s_next_fresh_row_boundary')
   pending={'row':i,'operation':'fresh_raw_coefficient_evaluation'};raw=evaluate(i);pending['operation']='transport_initial7_knownunit_pivots';row,cofs=reduce_row(raw)
   materialized[i]=raw;new_residual[i]=row;new_initial_coeffs[i]=cofs;new_completed.append(i);pending=None;checkpoint()
   print(json.dumps({'new_completed_row':i,'new_count':len(new_completed),'quartic':i>=113,'seconds':time.monotonic()-began}),flush=True)
  # Only after all134 are complete, test as-yet-unscanned row/stage pairs; old scalar data remain cached.
  while available and cols:
   if time.monotonic()-began>48:raise Stop('48s_next_gaussian_row_boundary')
   unseen=[i for i in available if[i,len(steps)]not in tested]
   if not unseen:stop='all134_exactly_scanned_no_individual_FULLU_unit_entry';break
   i=min(unseen,key=lambda i:sum(len(c.numer)for c in materialized[i]));pending={'row':i,'operation':'knownunit_Gaussian_on_cached_exact_row'}
   if i in new_residual:row,cofs=reduce_row(new_residual[i],initial_pivots,new_initial_coeffs[i])
   else:row,cofs=reduce_row(materialized[i])
   tested.append([i,len(steps)]);rowstatus[str(i)]={'last_test_stage':len(steps),'nonzero_entries':sum(bool(row[j])for j in cols),'numerator_terms':sum(len(row[j].numer)for j in cols)}
   hits=[]
   for j in cols:
    u=unit(row[j])
    if u is not None:hits.append((sum(len(c.numer)for c in row),len(row[j].numer),j,u))
   if not hits:pending=None;continue
   _,_,j,u=min(hits,key=lambda z:z[:3]);p=row[j];normalized=[guard(c/p)for c in row];assert normalized[j]==K.one and not any(normalized[st['pivot_column']]for st in steps)
   proof=K.from_expr(s.sympify(u['constant']))
   for f,n in zip(known,u['a_d_F_B_C_J_powers']):proof*=K(f)**n
   assert proof==p
   steps.append({'raw134_source_row':i,'pivot_column':j,'pivot':enc(p),'unit_factorization':u,'subtract_previous_pivot_coefficients':[enc(c)for c in cofs],'normalized23':[enc(c)for c in normalized],'_row':normalized})
   available.remove(i);cols.remove(j);pending=None;checkpoint();print(json.dumps({'completed_pivot':len(steps),'row':i,'column':j,'seconds':time.monotonic()-began}),flush=True)
 except(Stop,MemoryError)as exc:stop=str(exc)or type(exc).__name__
 signal.alarm(0);checkpoint()
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),parent/'report.json',pp,base/'report.json',ip]
 result={'id':ID,'initial96_rows_and7_unit_pivots_reused':True,'old96_scalar_re_evaluation_count':0,'new_row_ids_only':queries,'new_fully_completed_rows':len(new_completed),'all134_rows_complete':len(materialized)==134,
  'new_completed_global_H1_pivots':len(steps)-initial_pivots,'total_global_H1_pivots':len(steps),'new_target_unit_lower_from_original239':118+len(steps),'lossless_redundant_source_target_shape':[190-len(steps),121-len(steps)],
  'stop_reason':stop,'unfinished_row_not_promoted':pending,'last_completed_row_boundary_saved':True,'all134_source_rows_and_all21_quartics_retained':True,
  'no_automatic_free163_minus_pivots_claim':True,'pure98_coefficient_queries':0,'known38_missing42_unchanged':True,'no_newpoint_prime_GB_targetcompletion_genericrank_or_full121_claim':True,
  'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-began,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256','new_row_ids_only']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())

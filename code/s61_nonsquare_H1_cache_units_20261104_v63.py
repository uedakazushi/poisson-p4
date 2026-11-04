#!/usr/bin/env python3
"""ONE cache-only unit Gaussian; zero scalar-DAG queries and zero quartic expansion."""
import argparse,gzip,hashlib,json,os,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_H1_CACHE_UNITS_20261104_v63'
class Stop(Exception):pass
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
 signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(Stop('57s_internal_boundary')));signal.alarm(57)
 parent=ROOT/'reports/S61_NONSQUARE_H1_KNOWN_UNIT_CONTINUATION_20261104_v63';pp=parent/'checkpoint.json.gz';pr=read(parent/'report.json');assert sha(pp)==pr['evidence_sha256'][pp.name]
 saved=json.loads(gzip.decompress(pp.read_bytes()));raw_encoded=dict(saved['materialized_raw_rows']);assert len(raw_encoded)==113 and set(raw_encoded)==set(range(113))
 assert saved['completed_pivots']==7;old96=set(saved['initial96_row_ids_reused_without_scalar_queries']);fresh17=set(saved['new_completed_raw_and_initial7_residual_row_ids']);assert len(old96)==96 and len(fresh17)==17 and not old96&fresh17
 assert read(ROOT/'reports/S61_NONSQUARE_H1_KNOWN_UNIT_ELIMINATION_20261104_v63_REVIEW/report.json')['all_passed']
 R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();known=[a,d,500*a+90*e**2+3*e-24,(50*a+18*e**2-9*e-2)*d+60*e-40,100*a+18*e**2-9*e-2,5*a*d-6*e-4]
 decp=lambda p:R.from_dict({tuple(m):s.QQ(c)for m,c in p});dec=lambda c:K(decp(c[0]))/K(decp(c[1]));encp=lambda p:[[list(m),str(c)]for m,c in p.items()];enc=lambda c:[encp(c.numer),encp(c.denom)]
 steps=[]
 for st in saved['unit_pivots']:
  st=dict(st);st['_row']=[dec(c)for c in st['normalized23']];steps.append(st)
 initial=len(steps);cols=saved['remaining_columns'][:];available=[i for i in range(113)if i not in[st['raw134_source_row']for st in steps]]
 residual={i:[7,[dec(c)for c in row],[dec(c)for c in cofs]]for i,row,cofs in saved['new_row_residuals_after_initial7']};assert set(residual)==fresh17
 tested=[];pending=None;stop=None;raw_decoded={};old96_arithmetic_rows=[];opcount=0
 def guard(c):
  if max(len(c.numer),len(c.denom))>100000:raise Stop('100000_output_terms_per_entry')
  if time.monotonic()-began>53:raise Stop('53s_arithmetic_boundary')
  return c
 def op(tag,x,y):
  nonlocal opcount
  if tag=='m':
   if not x or not y:return K.zero
   bound=max(len(x.numer)*len(y.numer),len(x.denom)*len(y.denom))
  else:
   if not y:return x
   if not x:return-y
   bound=len(x.numer)+len(y.numer)if x.denom==y.denom else max(len(x.numer)*len(y.denom)+len(y.numer)*len(x.denom),len(x.denom)*len(y.denom))
  if bound>150000:raise Stop('150000_intermediate_term_bound')
  opcount+=1;return guard(x*y if tag=='m'else x-y)
 def current(i):
  if i not in residual:
   # The96 old rows are decoded/reduced only AFTER a genuinely new pivot exists.
   assert len(steps)>initial or i not in old96
   if i in old96:old96_arithmetic_rows.append(i)
   raw_decoded[i]=[dec(c)for c in raw_encoded[i]];start=0;row=raw_decoded[i][:];cofs=[]
  else:start,row,cofs=residual[i];row=row[:];cofs=cofs[:]
  for st in steps[start:]:
   c=row[st['pivot_column']];cofs.append(c)
   if c:row=[op('s',z,op('m',c,v))for z,v in zip(row,st['_row'])]
   assert not row[st['pivot_column']]
  residual[i]=[len(steps),row,cofs];return row,cofs
 unitcache={}
 def unit(c):
  if not c:return None
  key=(tuple(c.numer.items()),tuple(c.denom.items()))
  if key in unitcache:return unitcache[key]
  powers=[];ground=[]
  for p in[c.numer,c.denom]:
   ex=[]
   for f in known:
    n=0
    while p and not p.rem(f):p=p.exquo(f);n+=1
    ex.append(n)
   if not p.is_ground:unitcache[key]=None;return None
   powers.append(ex);ground.append(p)
  z={'constant':str(K(ground[0])/K(ground[1])),'a_d_F_B_C_J_powers':[x-y for x,y in zip(*powers)]};unitcache[key]=z;return z
 def checkpoint():
  state={'parent_checkpoint_path':str(pp.relative_to(ROOT)),'parent_checkpoint_sha256':sha(pp),'parent113_cached_raw_rows_are_immutable_and_reused':True,
   'initial7_pivots_reused':initial,'completed_pivots':len(steps),'remaining_rows_among_cached113':available,'remaining_columns':cols,
   'unit_pivots':[{k:v for k,v in st.items()if k!='_row'}for st in steps],
   'residual_cache_by_row_stage':[[i,stage,[enc(c)for c in row],[enc(c)for c in cofs]]for i,(stage,row,cofs)in residual.items()],
   'row_stage_pairs_tested':tested,'old96_first_arithmetic_rows_after_newpivot_only':old96_arithmetic_rows,'old96_arithmetic_before_newpivot_count':0,
   'all21_quartics_retained_lazy':True,'all_unvisited_cached_rows_retained_with_exact_unit_schur_recipe':True,
   'source_preimage_recipe':saved['lazy_full_source_pivot_lift'],'raw134_source_preconditioner':saved['raw_row_preconditioner'],
   'safe_redundant_source_target_shape':[190-len(steps),121-len(steps)],'no_free163_minus_pivots_claim':True,
   'scalar_DAG_query_count':0,'quartic_expansion_count':0,'pure98_coefficient_query_count':0,
   'cached_rows_at_current_stage':[i for i,(n,row,cofs)in residual.items()if n==len(steps)],'arithmetic_operations':opcount,'seconds':time.monotonic()-began}
  data=gzip.compress((json.dumps(state,separators=(',',':'))+'\n').encode(),compresslevel=3,mtime=0)
  if len(data)>30*1024*1024:raise Stop('30MiB_checkpoint_guard')
  tmp=out/'checkpoint.tmp';tmp.write_bytes(data);os.replace(tmp,out/'checkpoint.json.gz')
  (out/'progress.json').write_text(json.dumps({'new_pivots':len(steps)-initial,'total_H1_pivots':len(steps),'remaining_H1_targets':len(cols),'projected_rows_at_current_stage':len(state['cached_rows_at_current_stage']),'seconds':time.monotonic()-began},indent=2)+'\n')
 checkpoint()
 try:
  while available and cols:
   if time.monotonic()-began>48:raise Stop('48s_next_row_boundary')
   allowed=[i for i in available if i in fresh17 or len(steps)>initial]
   unseen=[i for i in allowed if[i,len(steps)]not in tested]
   if not unseen:
    stop='fresh17_have_no_individual_unit_no_old96_arithmetic'if len(steps)==initial else'all113_cached_rows_scanned_no_individual_FULLU_unit_entry';break
   # Costs are read from already stored sparse coefficients; no old-row arithmetic is needed.
   def cost(i):
    return sum(len(c.numer)for c in residual[i][1])if i in residual else sum(len(c[0])for c in raw_encoded[i])
   i=min(unseen,key=lambda i:(cost(i),i));pending={'row':i,'stage':len(steps),'operation':'project_cached_row_in_new_schur_state'};row,cofs=current(i);tested.append([i,len(steps)])
   hits=[]
   for j in cols:
    u=unit(row[j])
    if u is not None:hits.append((len(row[j].numer),j,u))
   if not hits:pending=None;continue
   _,j,u=min(hits,key=lambda z:z[:2]);p=row[j];norm=[guard(c/p)for c in row]
   assert norm[j]==K.one and not any(norm[st['pivot_column']]for st in steps)
   proof=K.from_expr(s.sympify(u['constant']))
   for f,n in zip(known,u['a_d_F_B_C_J_powers']):proof*=K(f)**n
   assert proof==p
   steps.append({'raw134_source_row':i,'pivot_column':j,'pivot':enc(p),'unit_factorization':u,'subtract_previous_pivot_coefficients':[enc(c)for c in cofs],'normalized23':[enc(c)for c in norm],'_row':norm})
   available.remove(i);cols.remove(j);pending=None;checkpoint();print(json.dumps({'new_completed_pivot':len(steps)-initial,'total_pivots':len(steps),'row':i,'column':j,'seconds':time.monotonic()-began}),flush=True)
 except(Stop,MemoryError)as exc:stop=str(exc)or type(exc).__name__
 signal.alarm(0);checkpoint()
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),parent/'report.json',pp,ROOT/'reports/S61_NONSQUARE_H1_KNOWN_UNIT_ELIMINATION_20261104_v63_REVIEW/report.json']
 result={'id':ID,'new_completed_global_H1_unit_pivots':len(steps)-initial,'total_global_H1_pivots':len(steps),'global_original239_target_unit_lower':118+len(steps),
  'safe_redundant_source_target_shape':[190-len(steps),121-len(steps)],'not_automatic_free163_minus_n':True,
  'scalar_DAG_queries':0,'quartic_expansions':0,'pure98_queries':0,'old96_arithmetic_before_newpivot_count':0,
  'stop_reason':stop,'interrupted_row_not_promoted':pending,'all113_cached_rows_and_all21_lazy_quartics_retained':True,
  'known38_missing42_unchanged':True,'no_allpoint_H23_full239_GB_point_prime_or_targetcompletion_claim':True,
  'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-began,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())

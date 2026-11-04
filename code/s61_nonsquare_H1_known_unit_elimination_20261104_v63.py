#!/usr/bin/env python3
"""ONE live exact-coefficient/known-unit search on the complete redundant H1 image."""
import argparse, gzip, hashlib, json, os, resource, signal, time
from collections import Counter
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_H1_KNOWN_UNIT_ELIMINATION_20261104_v63'
class Stop(Exception): pass
def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 began=time.monotonic(); resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(Stop('57s_internal_boundary')))
 signal.alarm(57); out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63'; rp=base/'report.json'; ip=base/'certificate.json'; report=read(rp)
 assert sha(ip)==report['certificate_sha256']; source=read(ip)
 assert source['complete_H1_image_block_shape']==[134,23]
 nodes=source['scalar_circuit']['nodes']; atoms=source['scalar_circuit']['atoms_num_den']; roots=source['complete_H1_image_block_scalar_roots']
 R,a,e,d=ring('a,e,d',s.QQ); K=R.to_field(); F=500*a+90*e**2+3*e-24; B=(50*a+18*e**2-9*e-2)*d+60*e-40; C=100*a+18*e**2-9*e-2; J=5*a*d-6*e-4
 known=[a,d,F,B,C,J]
 dec=lambda p:R.from_dict({tuple(m):s.QQ(c)for m,c in p})
 encp=lambda p:[[list(m),str(c)]for m,c in p.items()]
 enc=lambda c:[encp(c.numer),encp(c.denom)]
 cache={0:K.zero}; atomcache={}; evaluated_node_count=0; max_live_terms=0
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
   if x.denom==y.denom:bound=len(x.numer)+len(y.numer)
   else:bound=max(len(x.numer)*len(y.denom)+len(y.numer)*len(x.denom),len(x.denom)*len(y.denom))
  if bound>150000:raise Stop('150000_intermediate_term_bound')
  return guard(x*y if tag=='m'else x+y if tag=='a'else x-y)
 def evaluate(wanted):
  nonlocal evaluated_node_count,max_live_terms
  need=set(); todo=list(wanted)
  while todo:
   i=todo.pop()
   if i in cache or i in need:continue
   need.add(i);n=nodes[i]
   if n[0]!='c':assert n[0]in['a','s','m'];todo.extend(n[1:])
  refs=Counter(wanted)
  for i in need:
   n=nodes[i]
   if n[0]!='c':refs.update(n[1:])
  live={};live_terms=0
  for i in sorted(need):
   n=nodes[i]
   if n[0]=='c':
    ai=n[1]
    if ai not in atomcache:atomcache[ai]=K(dec(atoms[ai][0]))/K(dec(atoms[ai][1]))
    val=atomcache[ai]
   else:
    x=cache[n[1]]if n[1]in cache else live[n[1]];y=cache[n[2]]if n[2]in cache else live[n[2]];val=op(n[0],x,y)
   live[i]=val;live_terms+=len(val.numer)+len(val.denom);evaluated_node_count+=1;max_live_terms=max(max_live_terms,live_terms)
   if n[0]!='c':
    for dep in n[1:]:
     refs[dep]-=1
     if not refs[dep]and dep in live:
      z=live.pop(dep);live_terms-=len(z.numer)+len(z.denom)
  result=[cache[i]if i in cache else live[i]for i in wanted]
  for z in result:guard(z)
  return result
 def native_rows(path,header=True):
  lines=path.read_text().splitlines();ans=[]
  for line in lines[1:]if header else lines:
   t=line.split();n=int(t[0]);assert len(t)==1+4*n
   ans.append(R.from_dict({tuple(map(int,t[k+1:k+4])):s.QQ(t[k])for k in range(1,len(t),4)}))
  return ans
 # Exact invertible RAW134 source preconditioning; no target/pure coefficient queries.
 compact=ROOT/'reports/S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63'; cp=compact/'compact-polynomials.txt'; cr=read(compact/'report.json'); assert sha(cp)==cr['artifact_sha256'][cp.name]
 assert read(ROOT/'reports/S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63_REVIEW/report.json')['all_passed']
 pol=native_rows(cp); oldm2=[list(map(int,z.split()))for z in(compact/'source-metadata.txt').read_text().splitlines()[1:30]]
 hi=oldm2.index([0,0,0,0,0,1,0,0]);kZ=K(pol[8+hi]);kW=K(pol[8+29+hi]);delta=K(-64000*C**2);assert K(pol[74])==delta
 pre=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63'; inverse_path=pre/'coordinate-inverse-numerators7x7.txt';N=native_rows(inverse_path);Delta=K(3200*a*C*d**2)/81
 NN=[[K(N[7*i+j])/Delta for j in range(7)]for i in range(7)]
 oldB=[NN[0],NN[1],NN[2],NN[6],NN[3],[NN[5][j]+NN[1][j]-3*NN[0][j]for j in range(7)],NN[4]]
 hcols=source['target122_Hdegree_classes']['1'];mons=source['parent_joint_target239_monomials'];kept=source['parent_target122_monomial_indices'];hmons=[tuple(mons[kept[j]][:7])for j in hcols]
 assert len(hmons)==23 and all(sum(m)==2 for m in hmons)
 removed=set(tuple(int(j==i)+int(j==k)for j in range(7))for i,k in[(0,0),(0,2),(1,5),(1,6),(1,1)])
 assert set(hmons).isdisjoint(removed)
 assert len(set(hmons)|removed)==28
 special={}
 for j in range(7):
  full={}
  for k in range(7):
   for v,c in[(5,kZ),(6,kW)]:
    m=tuple(int(t==k)+int(t==v)for t in range(7));full[m]=op('a',full.get(m,K.zero),op('m',c,oldB[j][k]))
  # All five pivot coefficients except ellZ/ellW vanish literally, so NF117 only drops these two.
  assert not any(full.get(m,K.zero)for m in removed if m not in[(0,1,0,0,0,1,0),(0,1,0,0,0,0,1)])
  special[43+j]=[full.get(m,K.zero)for m in hmons]
 available=list(range(134));cols=list(range(23));materialized={};steps=[];tested=[];rowstatus={};stop=None;pending=None
 unitcache={}
 def unit(c):
  if not c:return None
  key=(tuple(c.numer.items()),tuple(c.denom.items()))
  if key in unitcache:return unitcache[key]
  powers=[];const=[]
  for p in[c.numer,c.denom]:
   ex=[]
   for f in known:
    k=0
    while p and not p.rem(f):p=p.exquo(f);k+=1
    ex.append(k)
   if not p.is_ground:unitcache[key]=None;return None
   powers.append(ex);const.append(p)
  z={'constant':str(K(const[0])/K(const[1])),'a_d_F_B_C_J_powers':[x-y for x,y in zip(*powers)]};unitcache[key]=z;return z
 def reduce_row(raw):
  cur=raw[:];cofs=[]
  for st in steps:
   c=cur[st['pivot_column']];cofs.append(c)
   if c:cur=[op('s',z,op('m',c,v))for z,v in zip(cur,st['_row'])]
   assert not cur[st['pivot_column']]
  return cur,cofs
 def checkpoint():
  state={'source_input_sha256':sha(ip),'coefficient_domain':'FULLU=QQ[a,e,d,1/(adFBCJ)]','shape_redundant_H1':[134,23],
   'completed_pivots':len(steps),'remaining_rows':available,'remaining_columns':cols,'materialized_raw_rows':[[i,[enc(c)for c in row]]for i,row in materialized.items()],
   'unit_pivots':[{k:v for k,v in st.items()if k!='_row'}for st in steps],
   'raw_row_preconditioner':{'replace_rows':list(range(43,50)),'formula':'row(43+j)=delta*(oldrow(43+j)-oldrow(8+j))','delta':enc(delta),'source_determinant':'delta^7','oldlight_inverse7x7':[[enc(c)for c in row]for row in oldB],'compact_H_coefficients':[enc(kZ),enc(kW)],'literal_H_projection_drops_only_ellZ_ellW_for_these_rows':True},
   'completed_evaluated_rows':list(materialized),'row_status':rowstatus,'tested_row_pivotstage_pairs':tested,'all_unexpanded_rows_retained_as_original_scalar_circuit_plus_saved_source_preconditioner':True,
   'redundant_full_source_target_after_pivots':[190-len(steps),121-len(steps)],'no_literal_free163_minus_pivots_claim':True,
   'lazy_full_source_pivot_lift':'For each selected rawHrow take its actual full row (including lazy pure98), subtract saved coefficients times earlier normalized full pivot rows, divide by saved knownunit. Apply unit target Schur using that full pivot row to EVERY remaining generator; pure98 remains lazy.',
   'pure98_coefficient_queries':0,'evaluated_scalar_nodes_with_recomputation':evaluated_node_count,'max_live_terms':max_live_terms,'seconds':time.monotonic()-began}
  data=gzip.compress((json.dumps(state,separators=(',',':'))+'\n').encode(),compresslevel=3,mtime=0)
  if len(data)>30*1024*1024:raise Stop('30MiB_checkpoint_artifact_guard')
  tmp=out/'checkpoint.tmp';tmp.write_bytes(data);os.replace(tmp,out/'checkpoint.json.gz')
  (out/'progress.json').write_text(json.dumps({k:state[k]for k in['completed_pivots','remaining_rows','remaining_columns','completed_evaluated_rows','row_status','seconds']},indent=2)+'\n')
 # A source-bound starting checkpoint is recoverable even if the first expansion hits a cap.
 checkpoint()
 def cost(i):
  if i in special:return sum(len(z.numer)for z in special[i])
  return sum(len(atoms[nodes[r][1]][0])if nodes[r][0]=='c'else 500 for r in roots[i]if r)
 todo=sorted(range(134),key=lambda i:(i not in special,i>=113,cost(i),i))
 try:
  while available and cols:
   if time.monotonic()-began>48:raise Stop('48s_next_row_boundary')
   unseen=[i for i in available if[i,len(steps)]not in tested]
   if not unseen:stop='all134_rows_exactly_scanned_no_individual_FULLU_unit_entry';break
   candidates=[i for i in unseen if i in materialized]
   i=min(candidates,key=lambda i:sum(len(z.numer)for z in materialized[i]))if candidates else next(i for i in todo if i in unseen)
   pending={'row':i,'stage':len(steps),'operation':'expansion'}
   if i not in materialized:
    materialized[i]=special[i]if i in special else evaluate(roots[i])
   row,coefficients=reduce_row(materialized[i]);tested.append([i,len(steps)]);rowstatus[str(i)]={'last_test_stage':len(steps),'nonzero_entries':sum(bool(row[j])for j in cols),'numerator_terms':sum(len(row[j].numer)for j in cols)}
   hits=[]
   for j in cols:
    u=unit(row[j])
    if u is not None:hits.append((sum(len(z.numer)for z in row),len(row[j].numer),j,u))
   if not hits:
    pending=None
    if len(materialized)%8==0:checkpoint()
    continue
   _,_,j,u=min(hits,key=lambda z:z[:3]);p=row[j];normalized=[guard(z/p)for z in row];assert normalized[j]==K.one
   assert not any(normalized[st['pivot_column']]for st in steps)
   proof=K.from_expr(s.sympify(u['constant']))
   for f,n in zip(known,u['a_d_F_B_C_J_powers']):proof*=K(f)**n
   assert proof==p
   st={'raw134_source_row':i,'pivot_column':j,'pivot':enc(p),'unit_factorization':u,'subtract_previous_pivot_coefficients':[enc(c)for c in coefficients],'normalized23':[enc(c)for c in normalized],'_row':normalized}
   steps.append(st);available.remove(i);cols.remove(j);pending=None;checkpoint();print(json.dumps({'completed_pivot':len(steps),'row':i,'column':j,'remaining_H_targets':len(cols),'seconds':time.monotonic()-began}),flush=True)
 except(Stop,MemoryError)as exc:stop=str(exc)or type(exc).__name__
 signal.alarm(0)
 # Save partial searches too, but an interrupted row is never promoted to a pivot.
 checkpoint();pending_saved=pending
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),rp,ip,compact/'report.json',cp,compact/'source-metadata.txt',inverse_path,
  ROOT/'reports/S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63_REVIEW/report.json',ROOT/'reports/S61_NONSQUARE_JOINT_LELL_QPAIR_UNITS_20261104_v63/checkpoint.json']
 result={'id':ID,'completed_global_H1_known_unit_pivots':len(steps),'new_target_unit_lower_from_original239':118+len(steps),
  'previous117_quadratic_plus_one_H2_unit_pivots':118,'remaining_H1_target_columns':len(cols),'lossless_redundant_source_target_shape':[190-len(steps),121-len(steps)],
  'no_automatic_free163_minus_pivots_claim':True,'unit_pivot_source_preimages_recursive_and_full_pure98_lazy':True,
  'stop_reason':stop,'interrupted_row_not_promoted':pending_saved,'rows_exactly_expanded':len(materialized),'all134_source_rows_retained':True,'all21_remaining_genuine_quartics_retained':True,
  'known38_missing42_unchanged':True,'no_pure98_queries_point_prime_GB_genericrank_or_full121_claim':True,
  'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},
  'artifact_bytes':sum(p.stat().st_size for p in out.iterdir()if p.is_file()),'seconds':time.monotonic()-began,'frozen':True,'independent_review_pending':True}
 assert result['artifact_bytes']<32*1024*1024
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())

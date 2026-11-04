#!/usr/bin/env python3
"""ONE bounded named-FULLU-unit Gaussian on actual augmented25x32 input."""
import argparse,hashlib,json,os,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
class Stop(Exception):pass
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(Stop('57s_internal_boundary')));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 inp=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63';rp=inp/'report.json';r=read(rp);data=read(inp/'input.json')
 assert r['frozen']and r['allowed_input100k_and32MiB_guard_passed']and sha(inp/'input.json')==r['input_json_sha256']and r['input_profile']['core_terms']<=100000
 R,a,e,d=ring('a,e,d',s.QQ);field=R.to_field();decode=lambda ts:R.from_dict({tuple(m):s.QQ(v)for m,v in ts});encodep=lambda p:[[list(m),str(v)]for m,v in p.items()];encode=lambda q:[encodep(q.numer),encodep(q.denom)]
 A=[[field(decode(v))for v in row]for row in data['matrix25x32']];assert len(A)==25 and all(len(row)==32 for row in A)
 known=[a,d,500*a+90*e**2+3*e-24,(50*a+18*e**2-9*e-2)*d+60*e-40,100*a+18*e**2-9*e-2,5*a*d-6*e-4]
 rows=list(range(25));cols=list(range(32));steps=[];constant=Q(1);powers=[0]*6;cache={};pending=None;stop=None
 def unit(q,fast=False):
  if not q:return None
  if len(q.numer)==len(q.denom)==1:
   (nm,nc),=q.numer.items();(dm,dc),=q.denom.items()
   if nm[1]==dm[1]==0:return{'constant':str(nc/dc),'a_d_F_B_C_J_powers':[nm[0]-dm[0],nm[2]-dm[2],0,0,0,0]}
  if fast:return None
  key=str(q)
  if key in cache:return cache[key]
  factors=[];ground=[]
  for p in[q.numer,q.denom]:
   ex=[]
   for f in known:
    n=0
    while not p.rem(f):p=p.exquo(f);n+=1
    ex.append(n)
   if not p.is_ground:cache[key]=None;return None
   factors.append(ex);ground.append(p)
  answer={'constant':str(field(ground[0])/field(ground[1])),'a_d_F_B_C_J_powers':[x-y for x,y in zip(*factors)]};cache[key]=answer;return answer
 def guard(q):
  if max(len(q.numer),len(q.denom))>100000:raise Stop('entry100000_output_terms')
  if time.monotonic()-began>53:raise Stop('53s_new_arithmetic_boundary')
 def profile():
  values=[A[i][j]for i in rows for j in cols if A[i][j]]
  return{'rows':len(rows),'columns':len(cols),'nonzero_entries':len(values),'numerator_terms':sum(len(v.numer)for v in values),'denominator_terms':sum(len(v.denom)for v in values),'max_entry_numerator_terms':max([len(v.numer)for v in values]+[0])}
 def checkpoint():
  record={'completed_pivots':len(steps),'source_input_sha256':sha(inp/'input.json'),'matrix25x32':[[encode(v)for v in row]for row in A],'remaining_rows':rows,'remaining_columns':cols,'all_completed_unit_steps':steps,
   'selected_minor_rows':[x['pivot_row']for x in steps],'selected_minor_columns':[x['pivot_column']for x in steps],'selected_minor_unit_determinant':{'constant':str(constant),'a_d_F_B_C_J_powers':powers},'remaining_profile':profile(),'seconds':time.monotonic()-began}
  text=json.dumps(record,separators=(',',':'))+'\n';old=(out/'checkpoint.json').stat().st_size if(out/'checkpoint.json').exists()else 0
  if len(text)+old>32*1024*1024:raise Stop('atomic_checkpoint32MiB_guard')
  tmp=out/'checkpoint.tmp';tmp.write_text(text);os.replace(tmp,out/'checkpoint.json')
  (out/'progress.json').write_text(json.dumps({k:record[k]for k in['completed_pivots','remaining_rows','remaining_columns','remaining_profile','seconds']},indent=2)+'\n')
 checkpoint()
 try:
  while rows and cols:
   if time.monotonic()-began>48:raise Stop('48s_next_pivot_boundary')
   rowcost={i:sum(len(A[i][j].numer)for j in cols if A[i][j])for i in rows};colnz={j:sum(bool(A[i][j])for i in rows)for j in cols};choices=[]
   for fast in[True,False]:
    for i in rows:
     for j in cols:
      q=A[i][j]
      if not q:continue
      u=unit(q,fast)
      if u is not None:choices.append((colnz[j]*rowcost[i],rowcost[i],i,j,u))
     if time.monotonic()-began>48:raise Stop('48s_unit_scan_boundary')
    if choices:break
   if not choices:stop='no_individual_named_FULLU_unit_entry';break
   _,_,i,j,u=min(choices,key=lambda z:z[:4]);p=A[i][j];check=field.from_expr(s.sympify(u['constant']))
   for f,n in zip(known,u['a_d_F_B_C_J_powers']):check*=field(f)**n
   assert p==check
   pending={'pivot_row':i,'pivot_column':j,'pivot':encode(p),'known_unit':u,'eliminations':[],'candidate_rule':'a,d monomial entries first; column_nonzeros*row_numerator_terms'}
   normalized=[v/p for v in A[i]]
   for v in normalized:guard(v)
   A[i]=normalized
   for t in rows:
    if t==i:continue
    c=A[t][j]
    if not c:continue
    new=A[t][:]
    for k in cols:
     if k==j:new[k]=field.zero;continue
     if not normalized[k]:continue
     bound=len(c.numer)*len(normalized[k].numer)*len(A[t][k].denom)+len(A[t][k].numer)*len(c.denom)*len(normalized[k].denom)
     if bound>150000:raise Stop('entry150000_intermediate_term_bound')
     new[k]=A[t][k]-c*normalized[k];guard(new[k])
    A[t]=new;pending['eliminations'].append([t,encode(c)])
   rows.remove(i);cols.remove(j);steps.append(pending);pending=None;constant*=Q(u['constant']);powers=[x+y for x,y in zip(powers,u['a_d_F_B_C_J_powers'])];checkpoint()
 except(Stop,MemoryError)as err:stop=str(err)or type(err).__name__
 signal.alarm(0);state=read(out/'checkpoint.json');n=state['completed_pivots'];full=n==25
 assert full or state['remaining_profile']['nonzero_entries']>0 # generic rank25 already certified, not re-computed
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),rp,inp/'input.json',ROOT/'reports/S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63_REVIEW/report.json']
 report={'id':ID,'completed_global_known_unit_pivots':n,'global_augmented25x32_rank25_proved':full,'full_cubic57_all_FULLU_proved_primary_pending_review':full,
  'selected_minor_unit_determinant':state['selected_minor_unit_determinant'],'remaining_shape':[len(state['remaining_rows']),len(state['remaining_columns'])],'remaining_profile':state['remaining_profile'],
  'stop_reason':stop,'unfinished_step_not_in_checkpoint':pending is not None,'unfinished_step_coordinates':[pending['pivot_row'],pending['pivot_column']]if pending else None,
  'checkpoint_last_fully_completed_only':True,'unknown_core_factors_never_inverted':True,'no_allpoint_negative_from_stall':True,
  'quartic239_genericbeta_and_knownkernel38_missing42_unchanged':True,'no_newpoint_prime_genericrank_GB_or_Delta10_determinant':True,
  'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-began,'frozen':True,'independent_review_pending':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())

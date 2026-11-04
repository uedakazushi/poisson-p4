#!/usr/bin/env python3
"""ONE bounded known-FULLU-unit elimination on the frozen58x136 quadratic block."""
import argparse,hashlib,json,os,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_JOINT_LELL_QPAIR_UNITS_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
class Stop(Exception):pass
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(Stop('57s_internal_boundary')));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_OVERLAP_UNIT_20261104_v63';rp=base/'report.json';r=read(rp);assert r['all_exact_fixed_monic_reductions_and_unit_operations_passed']and r['extra_global_known_unit_pivots']==5
 ip=base/'remaining-pair-block.json';assert sha(ip)==r['evidence_sha256'][ip.name];source=read(ip);assert source['complete58_rows']and len(source['entries'])==58 and len(source['entries'][0])==136
 R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();known=[a,d,500*a+90*e**2+3*e-24,(50*a+18*e**2-9*e-2)*d+60*e-40,100*a+18*e**2-9*e-2,5*a*d-6*e-4]
 encp=lambda p:[[list(m),str(c)]for m,c in p.items()];enc=lambda c:[encp(c.numer),encp(c.denom)];dec=lambda ts:R.from_dict({tuple(m):s.QQ(c)for m,c in ts})
 A=[[K(dec(n))/K(dec(de))for n,de in row]for row in source['entries']];available=list(range(58));cols=list(range(136));steps=[];detconst=Q(1);detpowers=[0]*6;ucache={};pending=None;stop=None
 def factor_unit(c,fast=False):
  if not c:return None
  if len(c.numer)==len(c.denom)==1:
   (n,cn),=c.numer.items();(de,cd),=c.denom.items()
   if n[1]==de[1]==0:return{'constant':str(cn/cd),'a_d_F_B_C_J_powers':[n[0]-de[0],n[2]-de[2],0,0,0,0]}
  if fast:return None
  key=str(c)
  if key in ucache:return ucache[key]
  ex=[];ground=[]
  for p in[c.numer,c.denom]:
   powers=[]
   for f in known:
    n=0
    while not p.rem(f):p=p.exquo(f);n+=1
    powers.append(n)
   if not p.is_ground:ucache[key]=None;return None
   ex.append(powers);ground.append(p)
  answer={'constant':str(K(ground[0])/K(ground[1])),'a_d_F_B_C_J_powers':[x-y for x,y in zip(*ex)]};ucache[key]=answer;return answer
 def guard(c):
  if max(len(c.numer),len(c.denom))>100000:raise Stop('entry100000output_terms')
  if time.monotonic()-began>53:raise Stop('53s_new_arithmetic_boundary')
 def profile():
  return{'nonzero_entries':sum(bool(A[i][j])for i in available for j in cols),'numerator_terms':sum(len(A[i][j].numer)for i in available for j in cols if A[i][j]),'denominator_terms':sum(len(A[i][j].denom)for i in available for j in cols if A[i][j]),'max_entry_numerator_terms':max([len(A[i][j].numer)for i in available for j in cols]+[0])}
 def checkpoint():
  state={'completed_pivots':len(steps),'source_input_sha256':sha(ip),'matrix58x136':[[enc(c)for c in row]for row in A],'remaining_rows':available,'remaining_columns':cols,'all_completed_unit_steps':steps,'selected_minor_rows':[p['pivot_row']for p in steps],'selected_minor_columns':[p['pivot_column']for p in steps],
   'selected_minor_unit_determinant':{'constant':str(detconst),'a_d_F_B_C_J_powers':detpowers},'remaining_profile':profile(),'seconds':time.monotonic()-began}
  data=json.dumps(state,separators=(',',':'))+'\n';old=(out/'checkpoint.json').stat().st_size if(out/'checkpoint.json').exists()else 0
  if len(data)+old>32*1024*1024:raise Stop('checkpoint_atomic32MiB_guard')
  tmp=out/'checkpoint.tmp';tmp.write_text(data);os.replace(tmp,out/'checkpoint.json')
  (out/'progress.json').write_text(json.dumps({k:state[k]for k in['completed_pivots','remaining_rows','remaining_columns','remaining_profile','seconds']},indent=2)+'\n')
 checkpoint()
 try:
  while available and cols:
   if time.monotonic()-began>48:raise Stop('48s_next_pivot_boundary')
   rowcost={i:sum(len(A[i][j].numer)for j in cols if A[i][j])for i in available};colnz={j:sum(bool(A[i][j])for i in available)for j in cols};candidates=[]
   for fast in[True,False]:
    for i in available:
     for j in cols:
      p=A[i][j]
      if not p:continue
      u=factor_unit(p,fast)
      if u is not None:candidates.append((colnz[j]*rowcost[i],rowcost[i],i,j,u))
     if time.monotonic()-began>48:raise Stop('48s_unit_scan_boundary')
    if candidates:break
   if not candidates:stop='no_individual_FULLU_unit_entry';break
   _,_,i,j,u=min(candidates,key=lambda x:x[:4]);p=A[i][j];proof=K.from_expr(s.sympify(u['constant']))
   for f,n in zip(known,u['a_d_F_B_C_J_powers']):proof*=K(f)**n
   assert proof==p
   pending={'pivot_row':i,'pivot_column':j,'pivot':enc(p),'known_unit':u,'cost_rule':'column_nonzeros*row_numerator_terms; a,d-monomial candidates first','eliminations':[]}
   pivotrow=[c/p for c in A[i]]
   for c in pivotrow:guard(c)
   A[i]=pivotrow
   for t in available:
    if t==i:continue
    c=A[t][j]
    if not c:continue
    newrow=A[t][:]
    for k in cols:
     if k==j:newrow[k]=K.zero;continue
     if not pivotrow[k]:continue
     bound=len(c.numer)*len(pivotrow[k].numer)*len(A[t][k].denom)+len(A[t][k].numer)*len(c.denom)*len(pivotrow[k].denom)
     if bound>150000:raise Stop('entry150000_intermediate_term_bound')
     newrow[k]=A[t][k]-c*pivotrow[k];guard(newrow[k])
    A[t]=newrow;pending['eliminations'].append([t,enc(c)])
   available.remove(i);cols.remove(j);steps.append(pending);pending=None;detconst*=Q(u['constant']);detpowers=[x+y for x,y in zip(detpowers,u['a_d_F_B_C_J_powers'])]
   checkpoint()
   assert len(steps)<=49 # the accepted globalupper152 minus the already exact103.
 except(Stop,MemoryError)as error:stop=str(error)or type(error).__name__
 signal.alarm(0)
 state=read(out/'checkpoint.json');complete=state['completed_pivots'];lastrows=state['remaining_rows'];lastcols=state['remaining_columns'];literalzero=state['remaining_profile']['nonzero_entries']==0
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),rp,ip,base/'overlap-unit-certificate.json',ROOT/'reports/S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63/report.json']
 report={'id':ID,'completed_global_known_unit_pivots':complete,'previous_exact_quadratic_image_unit_rank':103,'new_proved_quadratic_image_unit_lower':103+complete,
 'remaining_quadratic_block_shape':[len(lastrows),len(lastcols)],'remaining_block_literal_zero':literalzero,'allpoint_weight4_quad_rank152_proved':complete==49 and literalzero,
 'stop_reason':stop,'unfinished_step_not_in_checkpoint':pending is not None,'unfinished_step_coordinates':[pending['pivot_row'],pending['pivot_column']]if pending else None,
 'checkpoint_is_last_fully_completed_pivot':True,'all_pivots_have_exact_FULLU_unit_factorization':True,'selected_minor_determinant_is_product_of_saved_unit_pivots':True,
 'known38_missing42_unchanged':True,'no_full183matrix_33quartic_expansion_newpoint_prime_GB_or_generic_rank':True,
 'all33_genuine_quartics_remain_in_lossless_source':True,'remaining_profile':state['remaining_profile'],
 'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-began,'artifact_bytes':sum(p.stat().st_size for p in out.iterdir()if p.is_file()),'frozen':True,'independent_review_pending':True}
 assert report['artifact_bytes']<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())

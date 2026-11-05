#!/usr/bin/env python3
"""Prepared two transpose-R covectors and all1836 identities only."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_RHO_20261105_v63'
GATE='reports/AMBIENT_B18_ELL2_T2_UPPER617_RHO_INPUT_GATE_20261105_v63/'
GATE_SHA='f8d649d0050e81b6671119ef9825cb14e30391bfea1693fc2aacd653d08d3845'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(v):return {int(i):F(x) for i,x in v if F(x)}
def rec(v):return [[i,str(x)] for i,x in sorted(v.items()) if x]
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except (ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def main(out,peer,custody):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True)
 pins={};cursor={};pending={};accumulator=F(0);maximum=0;rhos=[];unitrecords=[];identities=[];adj={};w=[];rhs=[];sol={};deferred=[]
 def read(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def dump(name,v):
  (out/name).write_text(json.dumps(v,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in out.iterdir())<2_000_000
 g=read(GATE+'report.json');assert g['all_passed'] and pins[GATE+'report.json']==GATE_SHA
 ip=GATE+'pruned-transpose-R-unit28-source-and-positive-rho-input.json';d=read(ip);assert pins[ip]==g['evidence_sha256'][Path(ip).name]
 p=read(peer);c=read(custody);assert p['all_passed'] and GATE_SHA in json.dumps(p) and c['all_passed'] and c['safe_help_literal_PASS'] and c['worker_sha256']==sha(Path(__file__))
 for path,h in c['literal_input_sha256'].items():assert sha(ROOT/path)==h
 assert count([g,p,c,pins])<1000 and g['ready_with_unchanged384'] and g['all1836_F_identity_positive_prior_bits']==201 and g['temporary_prior_bits']==604
 assert g['future_value_numeric_bound']==26867 and count(d)==5499
 alpha=list(map(parse,d['alpha_two_by918_sparse_rows']));program=dict(d['pruned_scalar_linear_nodes']);outputs=d['full_source_outputs918'];unit=d['active_unit28'];U=unit['matrix_integer28_by28'];uin=unit['input_nodes28'];us=d['unit_component_node_ids28'];steps=d['all28_transpose_singleton_steps'];A=[{j:F(U[j][i]) for j in range(28) if U[j][i]} for i in range(28)];lo=min(us);hi=max(us)
 assert us==list(range(lo,hi+1)) and len(us)==len(steps)==28 and all(program[i]==['unit28_component',0,j] for j,i in enumerate(us)) and max(uin)<lo
 G2=list(map(parse,d['actual_full_G2linear_columns']));pairs=d['quadratic_pair_order44'];pindex={tuple(v):i for i,v in enumerate(pairs)};env=[{i:(int(D),int(V)) for i,D,V in row} for row in d['positive_rho_envelopes_by_channel']]
 initial=[x for a in alpha for x in a.values()]+[F(n[1]) for n in program.values() if n[0]=='scale']+[x for row in A for x in row.values()]+[x for col in G2 for x in col.values()]
 assert max(max(abs(x.numerator).bit_length(),x.denominator.bit_length()) for x in initial)<=201
 del initial
 parsed=count([alpha,program,A,G2,pindex,env]);cap=26867
 # one332-node sparse adjoint map <=664; allunitRHS/solutions/deferred <=224;
 # two142-entry rho rows plus serialization <=568;1836 zeros twice3672;
 # auxiliary serialized unit/proof fields512; metadata/frontier reserve6k.
 before=2*count(d)+parsed+664+224+568+3672+512+6000;assert before<cap
 dump('phase.json',{'phase':'literal_SOURCE_custody_and_allstage_bounds_PASS_before_signedrho','source_gate_sha256':GATE_SHA,'positive_committed_prior201':201,'temporary_prior604':604,'hard_guards':[384,1153],'whole_numeric_cap':cap,'prearith_live_upper':before,'no_formalRHS_sourceZ2_b1_or_upper_response_evaluation':True})
 def stop(reason,x):
  state={'reason':reason,'cursor':dict(cursor),'pending_operation':dict(pending),'partial_accumulator_before_operation':str(accumulator),'uncommitted_fraction':str(x),'current_adjoint_map':rec(adj),'current_unit_rhs28':list(map(str,w)),'current_mutated_rhs28':list(map(str,rhs)),'current_unit_solutions':rec(sol),'current_deferred_snapshots':[[i,str(v)] for i,v in deferred],'completed_rho_rows':[rec(v) for v in rhos],'completed_unit_records':unitrecords,'completed_original_identities':identities,'input_sha256':pins,'no_retry_or_cap_raise':True}
  assert 2*count(d)+parsed+3*count(state)+1000<cap;dump('guard-frontier.json',state);raise RuntimeError(reason)
 def ck(x):
  nonlocal maximum
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length());maximum=max(maximum,b)
  if b>201:stop('source201_before_commit',x)
  return x
 def mul(a,b):pending.update(operation='multiply',operand_a=str(a),operand_b=str(b));return ck(a*b)
 def add(a,b):pending.update(operation='add',operand_a=str(a),operand_b=str(b));return ck(a+b)
 def div(a,b):assert b;pending.update(operation='divide_known_literal_unit',operand_a=str(a),operand_b=str(b));return ck(a/b)
 def envelope(x,D,V):return D>0 and D%x.denominator==0 and abs(x.numerator)*D<=V*x.denominator
 def put(i,q):adj[i]=add(adj.get(i,F(0)),q)
 for t,a in enumerate(alpha):
  adj={};w=[];rhs=[];sol={};deferred=[]
  for slot,x in a.items():
   cursor={'phase':'seed_actual_alpha_into_R_outputs','channel':t,'source_slot918':slot};i=outputs[slot]
   if i:put(i,x)
  def reverse(i):
   nonlocal accumulator
   if i not in adj:return
   value=adj.pop(i);n=program[i];op=n[0];cursor.update(node_id=i,node_kind=op)
   if op=='scale':put(n[2],mul(value,F(n[1])))
   elif op=='sum':
    for offset,j in enumerate(n[1]):cursor['child_offset']=offset;accumulator=adj.get(j,F(0));put(j,value)
   elif op=='input':put(i,value)
   else:assert op=='constant' and F(n[1])==0
  for i in sorted((i for i in program if i>hi),reverse=True):cursor={'phase':'reverse_post_unit_DAG','channel':t};reverse(i)
  w=[adj.pop(i,F(0)) for i in us];rhs=list(w)
  for k,st in enumerate(steps):
   r=st['row28'];column=st['column28'];pv=F(st['pivot_integer']);assert A[r][column]==pv
   cursor={'phase':'fixed_transpose28_forward','channel':t,'step28':k,'row28':r,'column28':column}
   if st['kind']=='row_singleton':
    assert [j for j in st['remaining_cols_before'] if j in A[r]]==[column];value=div(rhs[r],pv);sol[column]=value
    for i in st['remaining_rows_before']:
     if i!=r and A[i].get(column):cursor['updated_rhs_row28']=i;accumulator=rhs[i];term=mul(A[i][column],value);rhs[i]=add(rhs[i],-term)
   else:
    assert st['kind']=='column_singleton';deferred.append([k,rhs[r]])
  for k,b in reversed(deferred):
   st=steps[k];r=st['row28'];column=st['column28'];q=b
   for j,cof in st['literal_row_remaining_coefficients']:
    if j!=column:
     cursor={'phase':'fixed_transpose28_deferred_backsolve','channel':t,'step28':k,'row28':r,'column28':column,'remaining_column28':j};accumulator=q;term=mul(F(cof),sol[j]);q=add(q,-term)
   sol[column]=div(q,F(st['pivot_integer']))
  assert len(sol)==28;y=[sol[j] for j in range(28)];zeros=[]
  for i in range(28):
   q=F(0)
   for j,cof in A[i].items():cursor={'phase':'all56_literal_Utranspose_equations','channel':t,'row28':i,'column28':j};accumulator=q;term=mul(cof,y[j]);q=add(q,term)
   if q!=w[i]:stop('fixed_unit_transpose_identity_failed',q)
   zeros.append('0')
  for k,x in enumerate(w):D,V=map(int,d['positive_unit28_rhs_solution_by_channel'][t]['rhs28'][k]);assert envelope(x,D,V)
  for k,x in enumerate(y):D,V=map(int,d['positive_unit28_rhs_solution_by_channel'][t]['solution28'][k]);assert envelope(x,D,V)
  unitrecords.append({'rhs28':list(map(str,w)),'solution28':list(map(str,y)),'literal_Utranspose_residuals28':zeros})
  for k,(i,value) in enumerate(zip(uin,y)):
   cursor={'phase':'ADD_unit_pullback_to_existing_prefix_adjoints','channel':t,'component28':k,'input_node':i};accumulator=adj.get(i,F(0));put(i,value)
  for i in sorted((i for i in program if i<lo),reverse=True):cursor={'phase':'reverse_prefix_DAG','channel':t};reverse(i)
  rho={program[i][1]:x for i,x in adj.items() if program[i][0]=='input' and x};assert set(rho)<=set(d['possible_rho_input_slots_by_channel'][t])
  for slot in range(1144):D,V=env[t].get(slot,(1,0));assert envelope(rho.get(slot,F(0)),D,V)
  rhos.append(rho);checks=[];identities.append(checks)
  for j,col in enumerate(G2):
   for b in range(9):
    q=F(0)
    for slot,cof in col.items():
     target,a0=divmod(slot,9);pair=tuple(sorted([a0,b]))
     if pair==(8,8):continue
     outputslot=44*target+pindex[pair];cursor={'phase':'all1836_original_rhoF_equals_alpha','channel':t,'sourceF2_102':j,'source_var9':b,'F1quad_slot1144':outputslot};accumulator=q;term=mul(cof,rho.get(outputslot,F(0)));q=add(q,term)
    cursor={'phase':'all1836_original_rhoF_equals_alpha','channel':t,'sourceF2_102':j,'source_var9':b,'pending':'subtract_literal_alpha'};q=add(q,-alpha[t].get(9*j+b,F(0)))
    if q:stop('original_full918_rhoF_identity_failed',q)
    checks.append('0')
 ep='actual-two-rho1144-covectors-and-full918-identities.json';ev={'rho_two_by1144_sparse_rows':list(map(rec,rhos)),'actual_transpose_unit28_two_records':unitrecords,'all1836_original_rhoF_minus_alpha_literal_zeros':identities,'all2288_rho_envelopes_including_structural0_checked':True,'prefix_existing_adjoints_ADDED_not_overwritten':True,'full_source_slot_convention':'input1144=44*F1target26+pair44;source918=9*F2source102+var9','scope':'ONLYknowntransposeRcoefficientsandunit/sourceidentities.NORHSsourceZ2/b1/upper617/physical/gauge/HF/length.'}
 owned=2*count(d)+parsed+count([adj,w,rhs,sol,deferred,rhos,unitrecords,identities,ev,g,p,c,pins,cursor,pending])+1000;assert owned<cap;dump(ep,ev)
 report={'id':STEM,'all_passed':True,'actual_rho_coefficients_only_no_RHS_or_sourceZ2':True,'rho_supports':list(map(len,rhos)),'all56_unit_transpose_checks_exact0':True,'all1836_original_source_checks_exact0':True,'source_prior201':201,'temporary_prior604':604,'hard_guards':[384,1153],'actual_maximum_bits':maximum,'actual_owned_upper':owned,'whole_numeric_cap':cap,'RAW_output_cap':2000000,'no_b1_upper_HF_or_new_functionvalues':True,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};dump('report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','rho_supports','all1836_original_source_checks_exact0','actual_maximum_bits','actual_owned_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output,a.source_peer,a.custody)

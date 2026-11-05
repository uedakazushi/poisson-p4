#!/usr/bin/env python3
"""Fresh transpose-R SOURCE/positive/control/custody peer, no signed rho."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,math,resource,signal,subprocess,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,count,fields,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_RHO_INPUT_GATE_20261105_v63_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,v):(o/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 b='reports/AMBIENT_B18_ELL2_T2_UPPER617_RHO_INPUT_GATE_20261105_v63/';r=rd(b+'report.json');assert r['all_passed'] and pins[b+'report.json']=='f8d649d0050e81b6671119ef9825cb14e30391bfea1693fc2aacd653d08d3845'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in r['evidence_sha256'].items():assert sha(ROOT/b/p)==h
 d=rd(b+'pruned-transpose-R-unit28-source-and-positive-rho-input.json');assert count(d)==5499
 rp='reports/AMBIENT_B18_ELL2_T2_FULL918_RIGHTSOLVE_SOURCE_20261105_v63/actual-factored918-rightsolve-source-circuit-and-all-generator-identities.json';z=fld(rp,['scalar_linear_nodes','active_unit_block28','full_source_outputs918','quadratic_pair_order44','actual_full_G2linear_columns']);N=z['scalar_linear_nodes'];unit=z['active_unit_block28'];outputs=z['full_source_outputs918'];U=unit['matrix_integer28_by28'];uin=unit['input_nodes28']
 assert d['active_unit28']==unit and d['full_source_outputs918']==outputs and d['quadratic_pair_order44']==z['quadratic_pair_order44'] and d['actual_full_G2linear_columns']==z['actual_full_G2linear_columns']
 ap='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_CLEARING_20261105_v63/actual-two-full918-alpha-covectors-and836-G3-identities.json';alpha=fld(ap,['alpha_two_by918_sparse_rows'])['alpha_two_by918_sparse_rows'];assert alpha==d['alpha_two_by918_sparse_rows'];alpha=list(map(sp,alpha))
 for p,h in [('reports/AMBIENT_B18_ELL2_T2_FULL918_RIGHTSOLVE_SOURCE_20261105_v63_REVIEW/report.json','33d2adb5b0f9715125adc2d1ce48d8d86a8483566c3be9d2829f064660bdb9fa'),('reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_CLEARING_20261105_v63_REVIEW/report.json','2b9a5ee58eda108ba19df265442cc226194afde86d84012af2a242d081184b5d')]:assert rd(p)['all_passed'] and pins[p]==h
 program=dict(d['pruned_scalar_linear_nodes']);assert len(program)==332 and all(n==N[i] for i,n in program.items())
 routing='reports/AMBIENT_B18_ELL2_T2_UPPER617_RHO_INPUT_ROUTING_SOURCE_20261105_v63/report.json';assert rd(routing)['all_passed'] and pins[routing]=='fd1ec079113b96ea3968c443285e5604480de0c55b695fffed934945088b849d'
 original_inputs=[n[1] for n in N if n[0]=='input'];assert len(original_inputs)==1144 and sorted(original_inputs)==list(range(1144))
 us=[i for i,n in enumerate(N) if n[0]=='unit28_component'];assert us==d['unit_component_node_ids28']==list(range(us[0],us[0]+28));lo,hi=us[0],us[-1];assert max(uin)<lo
 inputmap={i:n[1] for i,n in program.items() if n[0]=='input'};assert len(set(inputmap.values()))==len(inputmap) and all(0<=j<1144 for j in inputmap.values())
 possible=[];reach=[]
 for a in alpha:
  todo=[outputs[i] for i in a];seen=set()
  while todo:
   i=todo.pop()
   if i in seen:continue
   seen.add(i);n=N[i]
   kids=n[1] if n[0]=='sum' else [n[2]] if n[0]=='scale' else uin if n[0]=='unit28_component' else []
   for k in kids:assert k<i;todo.append(k)
  reach.append(sorted(seen));possible.append(sorted(N[i][1] for i in seen if N[i][0]=='input'))
 assert reach==d['reachable_node_ids_by_channel'] and possible==d['possible_rho_input_slots_by_channel'] and list(map(len,possible))==[62,80]
 assert set(program)==set().union(*(set(v) for v in reach),set(us))
 A=[{j:Q(U[j][i]) for j in range(28) if U[j][i]} for i in range(28)];steps=d['all28_transpose_singleton_steps'];rows=set(range(28));cols=set(range(28))
 for s in steps:
  i,j=s['row28'],s['column28'];assert s['remaining_rows_before']==sorted(rows) and s['remaining_cols_before']==sorted(cols) and Q(s['pivot_integer'])==A[i][j]
  if s['kind']=='row_singleton':assert set(A[i])&cols=={j}
  else:assert s['kind']=='column_singleton' and {k for k in rows if j in A[k]}=={i}
  assert s['literal_row_remaining_coefficients']==[[k,int(A[i][k])] for k in sorted(cols) if k in A[i]];rows.remove(i);cols.remove(j)
 assert not rows and not cols and len(steps)==28
 maximum=0
 def mark(a):
  nonlocal maximum
  maximum=max(maximum,a[0].bit_length(),a[1].bit_length());assert maximum<=201;return a
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def plus(a,b0):
  if not a[1]:return b0
  if not b0[1]:return a
  D=math.lcm(a[0],b0[0]);return mark((D,a[1]*(D//a[0])+b0[1]*(D//b0[0])))
 def mul(a,b0):return mark((a[0]*b0[0],a[1]*b0[1])) if a[1] and b0[1] else (1,0)
 def divide(a,pv):return mark((a[0]*abs(pv),a[1])) if a[1] else (1,0)
 Z=(1,0);B=[];unitbounds=[]
 for a in alpha:
  adj={}
  def acc(i,q):adj[i]=plus(adj.get(i,Z),q)
  for slot,x in a.items():
   i=outputs[slot]
   if i:acc(i,ex(x))
  def rev(i):
   if i not in adj:return
   v=adj.pop(i);n=program[i]
   if n[0]=='scale':acc(n[2],mul(v,ex(n[1])))
   elif n[0]=='sum':
    for j in n[1]:acc(j,v)
   elif n[0]=='input':acc(i,v)
   else:assert n[0]=='constant' and Q(n[1])==0
  for i in sorted((i for i in program if i>hi),reverse=True):rev(i)
  w=[adj.pop(i,Z) for i in us];rhs=list(w);sol={};defer=[]
  for s in steps:
   i,j=s['row28'],s['column28'];pv=int(s['pivot_integer'])
   if s['kind']=='row_singleton':
    sol[j]=divide(rhs[i],pv)
    for k in s['remaining_rows_before']:
     if k!=i and j in A[k]:rhs[k]=plus(rhs[k],mul(ex(A[k][j]),sol[j]))
   else:defer.append((s,rhs[i]))
  for s,q in reversed(defer):
   i,j=s['row28'],s['column28']
   for k,c in s['literal_row_remaining_coefficients']:
    if k!=j:q=plus(q,mul(ex(c),sol[k]))
   sol[j]=divide(q,int(s['pivot_integer']))
  y=[sol[j] for j in range(28)]
  for i in range(28):
   q=Z
   for j,c in A[i].items():q=plus(q,mul(ex(c),y[j]))
   plus(q,w[i])
  unitbounds.append({'rhs28':[[str(D),str(V)] for D,V in w],'solution28':[[str(D),str(V)] for D,V in y]})
  for i,q in zip(uin,y):acc(i,q)
  for i in sorted((i for i in program if i<lo),reverse=True):rev(i)
  B.append({program[i][1]:v for i,v in adj.items() if program[i][0]=='input' and v[1]})
 assert maximum==r['reverse_program_positive_prior_bits']==200
 pairs=d['quadratic_pair_order44'];assert pairs==[list(v) for v in __import__('itertools').combinations_with_replacement(range(9),2) if v!=(8,8)];pix={tuple(v):i for i,v in enumerate(pairs)};G=list(map(sp,d['actual_full_G2linear_columns']))
 for t0,row in enumerate(B):
  for j,c in enumerate(G):
   for v in range(9):
    q=Z
    for k,x in c.items():
     target,u=divmod(k,9);pair=tuple(sorted([u,v]))
     if pair!=(8,8):q=plus(q,mul(ex(x),row.get(44*target+pix[pair],Z)))
    plus(q,ex(alpha[t0].get(9*j+v,Q(0))))
 assert maximum==r['all1836_F_identity_positive_prior_bits']==201 and 3*201+1==r['temporary_prior_bits']==604
 assert unitbounds==d['positive_unit28_rhs_solution_by_channel'] and [[[i,str(D),str(V)] for i,(D,V) in sorted(v.items())] for v in B]==d['positive_rho_envelopes_by_channel']
 cb='reports/AMBIENT_B18_ELL2_T2_UPPER617_RHO_PREPARED_CUSTODY_20261105_v63/';cu=rd(cb+'report.json');assert cu['all_passed'] and pins[cb+'report.json']=='f522d744b963bf03eb1bd88a09823ee2a4bc0e9ea82cfaa5eb241c3f51d16fc0'
 for p,h in cu['literal_input_sha256'].items():assert sha(ROOT/p)==h
 worker=ROOT/cu['worker'];assert sha(worker)==cu['worker_sha256']=='4f13266478d03b2d5d7a3294b583748e07371bf5ebf114b7d34df93a8a4725d1'
 helper=subprocess.run([str(ROOT/'.venv/bin/python'),str(worker),'--help'],capture_output=True,timeout=10);assert helper.returncode==0 and hashlib.sha256(helper.stdout).hexdigest()==cu['safe_help_stdout_sha256']
 code=worker.read_text();assert "for k,(i,value) in enumerate(zip(uin,y))" in code and "accumulator=adj.get(i,F(0));put(i,value)" in code and "def put(i,q):adj[i]=add(adj.get(i,F(0)),q)" in code
 assert "w=[adj.pop(i,F(0)) for i in us];rhs=list(w)" in code and "for k,b in reversed(deferred)" in code and "if q:stop('original_full918_rhoF_identity_failed',q)" in code
 assert all(s in code for s in ["'current_adjoint_map':rec(adj)","'current_unit_rhs28'","'current_mutated_rhs28'","'current_unit_solutions'","'current_deferred_snapshots'","'completed_rho_rows'","'completed_unit_records'","'completed_original_identities'","'pending_operation':dict(pending)","'no_retry_or_cap_raise':True"])
 env=[{i:(int(D),int(V)) for i,D,V in row} for row in d['positive_rho_envelopes_by_channel']];parsed=count([alpha,program,A,G,pix,env]);live=2*count(d)+parsed+664+224+568+3672+512+6000;assert live==cu['success_live_upper']==24679<26867
 ledger=rd(cb+'allstage-live-RAW-and-preoperation-ledger.json');assert pins[cb+'allstage-live-RAW-and-preoperation-ledger.json']==cu['ownership_ledger_sha256']
 assert ledger['guard_live_upper']==2*5499+parsed+3*3100+1000==23339==cu['guard_live_upper']<26867
 # At most one completed918-checkrow before the currentchannel reverse;
 # afterprefix the AD map contains at most80inputentries, not332nodes.
 reverse_guard=664+56+56+56+124+168+918+100
 checking_guard=160+56+56+56+284+168+1836+100
 assert max(reverse_guard,checking_guard)<3100
 rawset=set(cu['literal_input_sha256'])|{cu['worker'],cb+'report.json',cb+'allstage-live-RAW-and-preoperation-ledger.json'};rawbytes=sum((ROOT/p).stat().st_size for p in rawset);assert rawbytes+3100000<33554432
 proof={'literal_original_alpha_R332_unit_inputs_and_all102F_columns':True,'INPUT_node_to1144_slot_injective':True,'all28_Utranspose_literal_nofill_singletons_and_reverse_deferred_source':True,'all_source_reachable_inputs62_80_and_structuralzeros':True,'prefix_existing_ADJOINTS_ADDED_never_overwritten':True,'all2288_positive_rho_envelopes_and56_unit_envelopes_replayed':True,'all1836_originalF_identity_bounds201':True,'reverse200_G3identity201_temp604':True,'worker_literal_hash_safehelp_and_nohelpers':True,'success_live24679_guard23339_cap26867':True,'guard_state_stage_bounds':[reverse_guard,checking_guard],'whole_RAW_independent_upper':rawbytes+3100000,'custody_whole_RAW3161559':cu['whole_RAW_upper'],'scope':'rhoand56unit/1836sourceidentitiesONLY. NOformalRHSsourceZ2/b1/upper/HF/physical/gauge/highercallback.','no_signed_rho_or_unitadjoint_values_evaluated':True}
 ep='independent-rho-transpose-source-positive-private-state-and-custody.json';put(ep,proof)
 owned=count([r,d,z,N,unit,outputs,U,alpha,program,inputmap,possible,reach,A,steps,B,unitbounds,G,cu,ledger,proof,pins])+14000;assert owned<150000
 report={'id':STEM,'all_passed':True,'SOURCE_positive_code_custody_only_no_signedrho':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(o/ep)},'owned_numeric_bound':owned,'seconds':time.monotonic()-t};put('report.json',report);print(json.dumps({'all_passed':True,'sha256':sha(o/'report.json'),'owned':owned}))
if __name__=='__main__':main()

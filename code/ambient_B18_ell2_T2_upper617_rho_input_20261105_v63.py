#!/usr/bin/env python3
"""Pruned transpose-R/source prior only. No signed rho or RHS values."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER617_RHO_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={}
 def field(p,k):pins[p]=H.sha(ROOT/p);return H.selected(p,k)
 def read(p):pins[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 rp='reports/AMBIENT_B18_ELL2_T2_FULL918_RIGHTSOLVE_SOURCE_20261105_v63/actual-factored918-rightsolve-source-circuit-and-all-generator-identities.json'
 nodes=field(rp,'scalar_linear_nodes');unit=field(rp,'active_unit_block28');outputs=field(rp,'full_source_outputs918');pairs=field(rp,'quadratic_pair_order44');G2=field(rp,'actual_full_G2linear_columns');assert len(nodes)==2349 and len(outputs)==918 and len(pairs)==44 and len(G2)==102
 ap='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_CLEARING_20261105_v63/actual-two-full918-alpha-covectors-and836-G3-identities.json';alpha=list(map(H.parse,field(ap,'alpha_two_by918_sparse_rows')))
 assert list(map(len,alpha))==[34,52] and all(i//9<42 and i%9 in [0,4,5,7,8] for a in alpha for i in a)
 for p in ['reports/AMBIENT_B18_ELL2_T2_FULL918_RIGHTSOLVE_SOURCE_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_CLEARING_20261105_v63_REVIEW/report.json']:assert read(p)['all_passed']
 us=[i for i,n in enumerate(nodes) if n[0]=='unit28_component'];assert len(us)==28 and us==list(range(us[0],us[0]+28)) and all(nodes[i]==['unit28_component',0,j] for j,i in enumerate(us));lo=min(us);hi=max(us);uin=unit['input_nodes28'];U=unit['matrix_integer28_by28'];assert len(uin)==len(U)==28 and all(len(v)==28 for v in U) and max(uin)<lo
 # Original U^T has a complete no-fill singleton peeling. Its coefficients
 # never change; only a FUTURE right-hand side would be updated numerically.
 A=[{j:U[j][i] for j in range(28) if U[j][i]} for i in range(28)];rows=set(range(28));cols=set(range(28));steps=[]
 while rows:
  found=next((('row_singleton',r,[c for c in A[r] if c in cols][0]) for r in sorted(rows) if sum(c in cols for c in A[r])==1),None)
  if found is None:found=next((('column_singleton',[r for r in rows if c in A[r]][0],c) for c in sorted(cols) if sum(c in A[r] for r in rows)==1),None)
  assert found is not None,'fixedunittransposegraphfrontier_not_numerical_rank'
  kind,r,c=found;pv=A[r][c];assert pv
  steps.append({'kind':kind,'row28':r,'column28':c,'pivot_integer':pv,'remaining_rows_before':sorted(rows),'remaining_cols_before':sorted(cols),'literal_row_remaining_coefficients':[[j,A[r][j]] for j in sorted(cols) if j in A[r]]});rows.remove(r);cols.remove(c)
 support=read('reports/AMBIENT_B18_ELL2_T2_ALPHA_ACTIVE_SUPPORT_SOURCE_20261105_v63/literal-alpha-active-support-and-pruned-rightsolve-dependency.json')
 reach=[];possible=[]
 for t,a in enumerate(alpha):
  todo=[outputs[i] for i in a];seen=set()
  while todo:
   i=todo.pop()
   if i in seen:continue
   seen.add(i);n=nodes[i];op=n[0]
   if op=='sum':kids=n[1]
   elif op=='scale':kids=[n[2]]
   elif op=='unit28_component':kids=uin
   else:assert op in ['input','constant'];kids=[]
   assert all(k<i for k in kids);todo.extend(k for k in kids if k not in seen)
  assert seen==set(support['all_reachable_scalar_nodes'][t])
  inp=sorted(nodes[i][1] for i in seen if nodes[i][0]=='input');assert inp==support['all_possible_rho_input_slots_by_literal_DAG_dependency'][t]
  assert all(k//44<17 and all(v in [0,1,4,5,7,8] for v in pairs[k%44]) for k in inp);reach.append(seen);possible.append(inp)
 keep=set().union(*reach)|set(us);program={i:nodes[i] for i in keep}
 sourcebase=count([nodes,unit,outputs,pairs,G2,alpha,U,A,steps,support,reach,possible,program,pins]);assert sourcebase+20000<150000
 maximum=0
 def mark(q):
  nonlocal maximum
  maximum=max(maximum,q[0].bit_length(),q[1].bit_length());assert maximum<2048;return q
 def ex(x):x=F(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def add(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def divide_integer(a,pv):return mark((a[0]*abs(pv),a[1])) if a[1] else (1,0)
 zero=(1,0)
 def solve_bound(w):
  rhs=list(w);sol={};deferred=[]
  for st in steps:
   r=st['row28'];c=st['column28'];pv=st['pivot_integer']
   if st['kind']=='row_singleton':
    v=divide_integer(rhs[r],pv);sol[c]=v
    for i in st['remaining_rows_before']:
     if i!=r and A[i].get(c):rhs[i]=add(rhs[i],mul(ex(A[i][c]),v))
   else:deferred.append((st,rhs[r]))
  for st,b in reversed(deferred):
   r=st['row28'];c=st['column28'];q=b
   for j,a in st['literal_row_remaining_coefficients']:
    if j!=c:q=add(q,mul(ex(a),sol[j]))
   sol[c]=divide_integer(q,st['pivot_integer'])
  assert len(sol)==28
  for i in range(28):
   q=zero
   for j,a in A[i].items():q=add(q,mul(ex(a),sol[j]))
   add(q,w[i])
  return [sol[j] for j in range(28)]
 bounds=[];unit_bounds=[];stagemax=[]
 for t,a in enumerate(alpha):
  adj={}
  def put(i,q):adj[i]=add(adj.get(i,zero),q)
  for slot,x in a.items():
   i=outputs[slot]
   if i:put(i,ex(x))
  def reverse(i):
   if i not in adj:return
   w=adj.pop(i);n=program[i];op=n[0]
   if op=='scale':put(n[2],mul(w,ex(n[1])))
   elif op=='sum':
    for j in n[1]:put(j,w)
   elif op=='input':put(i,w)
   else:assert op=='constant' and F(n[1])==0
  for i in sorted((i for i in keep if i>hi),reverse=True):reverse(i)
  beforeunit=maximum;w=[adj.pop(i,zero) for i in us];y=solve_bound(w);unit_bounds.append({'rhs28':[[str(D),str(V)] for D,V in w],'solution28':[[str(D),str(V)] for D,V in y]});afterunit=maximum
  # Existing prefix adjoints from direct paths are ADDED, never replaced.
  for i,q in zip(uin,y):put(i,q)
  for i in sorted((i for i in keep if i<lo),reverse=True):reverse(i)
  rho={nodes[i][1]:v for i,v in adj.items() if nodes[i][0]=='input' and v[1]};assert set(rho)<=set(possible[t]);bounds.append(rho);stagemax.append([beforeunit,afterunit,maximum])
 reverse_prior=maximum
 # The complete102*9 source-column identities are tested after values, not
 # just selected pivots or source supports. Include their operation bounds.
 pindex={tuple(v):i for i,v in enumerate(pairs)};G=list(map(H.parse,G2))
 for t,rho in enumerate(bounds):
  for j,col in enumerate(G):
   for b in range(9):
    q=zero
    for slot,x in col.items():
     target,a=divmod(slot,9);pair=tuple(sorted([a,b]))
     if pair==(8,8):continue
     q=add(q,mul(ex(x),rho.get(44*target+pindex[pair],zero)))
    add(q,ex(alpha[t].get(9*j+b,F(0))))
 allprior=maximum
 def rec(v):return [[i,str(x)] for i,x in sorted(v.items())]
 payload={'alpha_two_by918_sparse_rows':list(map(rec,alpha)),'pruned_scalar_linear_nodes':sorted(program.items()),'full_source_outputs918':outputs,'active_unit28':unit,'all28_transpose_singleton_steps':steps,'reachable_node_ids_by_channel':list(map(sorted,reach)),'possible_rho_input_slots_by_channel':possible,'positive_rho_envelopes_by_channel':[[[i,str(D),str(V)] for i,(D,V) in sorted(v.items())] for v in bounds],'positive_unit28_rhs_solution_by_channel':unit_bounds,'quadratic_pair_order44':pairs,'actual_full_G2linear_columns':G2,'unit_component_node_ids28':us,'transpose_recipe':'reversepostblock;collectALL28componentadjoints;solveU^Tbyfixed28singletoncircuit;ADDtoearlierinputnodeadjoints;reverseprefix;rho1144.Otherinputslotsstructural0.','all1836_original_source_identities_required':'rhoF=alphaonALL9*102slots;ell²term0ONLYbecausecoefficientring.','stage_positive_maxima':stagemax,'no_actual_rho_unitRHS_b1_Z2_or_upper_values':True,'no_source_callback_descent_claim':True}
 n=count(payload);future=2*n+count([program,U,A,steps,alpha,G,bounds,unit_bounds])+12000;assert future<60000
 ep=out/'pruned-transpose-R-unit28-source-and-positive-rho-input.json';H.dump(ep,payload)
 report={'id':STEM,'all_passed':True,'source_positive_only_no_signed_rho':True,'no_fill_unit_transpose_steps28':28,'pruned_scalar_nodes':len(program),'possible_rho_supports':list(map(len,possible)),'reverse_program_positive_prior_bits':reverse_prior,'all1836_F_identity_positive_prior_bits':allprior,'temporary_prior_bits':3*allprior+1,'ready_with_unchanged384':allprior<=384,'source_owned_bound':sourcebase+20000,'future_value_numeric_bound':future,'pruned_numeric_records':n,'future_RAW_output_cap':2000000,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','no_fill_unit_transpose_steps28','pruned_scalar_nodes','possible_rho_supports','reverse_program_positive_prior_bits','all1836_F_identity_positive_prior_bits','temporary_prior_bits','ready_with_unchanged384','source_owned_bound','future_value_numeric_bound','pruned_numeric_records','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

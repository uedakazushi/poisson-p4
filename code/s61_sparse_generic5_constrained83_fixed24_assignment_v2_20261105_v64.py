#!/usr/bin/env python3
"""One fixed24 exact leading assignment certificate, no field rank/search."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,time,resource,signal
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_ASSIGNMENT_20261105_v64';GATE='S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_ASSIGNMENT_SOURCE_GATE_V2_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,list):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0
def leading(blob,weights):
 at=0
 def u():
  nonlocal at
  n=s=0
  while True:
   assert at<len(blob);v=blob[at];at+=1;n|=(v&127)<<s
   if not v&128:return n
   s+=7
 def atom():
  nonlocal at
  tag=u();n=tag//2;start=at;at+=n;assert at<=len(blob);return start,n,tag&1
 def value(t):
  start,n,neg=t;v=int.from_bytes(blob[start:start+n],'big');return-v if neg else v
 content=Q(value(atom()),value(atom()));assert content
 shift=[]
 for _ in range(5):
  z=u();shift.append(z//2 if not z&1 else-(z//2)-1)
 n=u();best=None
 for j in range(n):
  c=atom();e=[shift[k]+u()for k in range(5)];cost=sum(a*b for a,b in zip(weights,e))
  if best is None or cost<best['weight']:
   coef=content*value(c);assert coef;best={'weight':cost,'exponent5':e,'QQ_coefficient':str(coef),'source_term_index':j}
  elif cost==best['weight']:assert e==best['exponent5'],'Kronecker collision'
 assert at==len(blob)and n>0;return n,best

def hungarian(C):
 n=len(C);largest=max(abs(x)for row in C for x in row if x is not None);INF=(8*n+1)*(2*largest+1);A=[[INF if x is None else x for x in row]for row in C];u=[0]*(n+1);v=[0]*(n+1);p=[0]*(n+1);way=[0]*(n+1);steps=[]
 for i in range(1,n+1):
  p[0]=i;j0=0;dist=[INF*4]*(n+1);used=[False]*(n+1)
  while True:
   used[j0]=True;i0=p[j0];delta=INF*4;j1=0
   for j in range(1,n+1):
    if used[j]:continue
    cur=A[i0-1][j-1]-u[i0]-v[j]
    if cur<dist[j]:dist[j]=cur;way[j]=j0
    if dist[j]<delta:delta=dist[j];j1=j
   assert j1
   for j in range(n+1):
    if used[j]:u[p[j]]+=delta;v[j]-=delta
    else:dist[j]-=delta
   j0=j1
   if p[j0]==0:break
  while True:
   j1=way[j0];p[j0]=p[j1];j0=j1
   if j0==0:break
  steps.append({'inserted_row':i-1,'u':u[1:],'v':v[1:],'column_to_row':[x-1 for x in p[1:]]})
 match=[None]*n
 for j in range(1,n+1):match[p[j]-1]=j-1
 assert all(C[i][match[i]] is not None for i in range(n)),'fixed24 structural block has no finite perfect matching'
 reduced=[[None if C[i][j]is None else C[i][j]-u[i+1]-v[j+1]for j in range(n)]for i in range(n)]
 assert all(x is None or x>=0 for row in reduced for x in row)and all(reduced[i][match[i]]==0 for i in range(n))
 assert sum(C[i][match[i]]for i in range(n))==sum(u[1:])+sum(v[1:])
 return match,u[1:],v[1:],reduced,steps

def main():
 begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 gp=ROOT/f'reports/{GATE}/report.json';g=read(gp);assert g['all_passed']and g['ONE_fixed_assignment_job_safe'];bp=gp.parent/'whole-worker-source-gate.json';assert sha(bp)==g['evidence_sha256'][bp.name];b=read(bp)
 prepared_gate_bytes=gp.stat().st_size+bp.stat().st_size;assert prepared_gate_bytes<=131072
 actual_read_bytes=b['actual_READ_input_RAW_bytes']+prepared_gate_bytes
 for p,h in b['actual_READ_input_sha256'].items():assert sha(ROOT/p)==h
 assert sha(Path(__file__))==b['immutable_worker_sha256']
 sp=ROOT/b['source_plan'];frame=read(sp);rows=frame['rows24_fixed_before_any_rank_or_matching_test'];cols=frame['columns24'];weights=frame['injective_Kronecker_weights5'];n=24;assert len(rows)==len(cols)==n
 width=frame['determinant_box_width5'];want=[1]*5
 for j in range(3,-1,-1):want[j]=want[j+1]*width[j+1]
 assert weights==want and all(width[j]==24*(frame['global_entry_exponent_upper5'][j]-frame['global_entry_exponent_lower5'][j])+1 for j in range(5))
 entries=frame['actual_entry_sources576'];C=[[None]*n for _ in range(n)];leaders=[];by={};peak=numeric(frame)+25000;payload_peak=0
 for z in entries:
  i=rows.index(z['row']);j=cols.index(z['column']);ref=z['source'];blob=b''
  if z['terms']:
   p=ROOT/ref['file']
   if ref['kind']=='unchanged20_column_literal':
    with p.open('rb')as f:f.seek(ref['offset']);blob=f.read(ref['bytes'])
    assert hashlib.sha256(blob).hexdigest()==ref['payload_sha256']
   else:
    blob=p.read_bytes();assert sha(p)==ref['sha256']and blob[:9]==b'Q5CELL1\n\x01';blob=blob[9:]
   count,lead=leading(blob,weights);assert count==z['terms']
   assert all(z['lower5'][k]<=lead['exponent5'][k]<=z['upper5'][k]for k in range(5));C[i][j]=lead['weight'];q={'i':i,'j':j,'row':z['row'],'column':z['column'],'source':ref,**lead};leaders.append(q);by[i,j]=q
  payload_peak=max(payload_peak,z['terms']);peak=max(peak,numeric(frame)+numeric(leaders)+numeric(C)+2*payload_peak+25000);assert peak<=b['whole_owned_numeric_upper']
 assert sum(z['terms']for z in entries)==frame['logical_terms_not_resident']
 match,u,v,reduced,steps=hungarian(C)
 selected=[]
 # Independently rescan the24 chosen actual monomials, not just cost metadata.
 for i,j in enumerate(match):
  q=by[i,j];ref=q['source'];p=ROOT/ref['file']
  if ref['kind']=='unchanged20_column_literal':
   with p.open('rb')as f:f.seek(ref['offset']);blob=f.read(ref['bytes'])
   assert hashlib.sha256(blob).hexdigest()==ref['payload_sha256']
  else:blob=p.read_bytes()[9:]
  _,again=leading(blob,weights);assert all(again[k]==q[k]for k in again);selected.append(q)
 owner={j:i for i,j in enumerate(match)};edges=[sorted({owner[j]for j in range(n)if reduced[i][j]==0 and j!=match[i]})for i in range(n)];colors=[0]*n;stack=[];cycle=None;post=[]
 def dfs(i):
  nonlocal cycle
  colors[i]=1;stack.append(i)
  for j in edges[i]:
   if cycle is not None:break
   if colors[j]==0:dfs(j)
   elif colors[j]==1:cycle=stack[stack.index(j):]+[j];break
  stack.pop();colors[i]=2;post.append(i)
 for i in range(n):
  if colors[i]==0 and cycle is None:dfs(i)
 unique=cycle is None;certificate={'all_passed':True,'fixed_rows24':rows,'fixed_columns24':cols,'weights5':weights,'all_actual_leading_entries':leaders,'matching24':match,'selected24_monomials_rechecked':selected,'row_dual24':u,'column_dual24':v,'all576_reduced_costs':reduced,'all_reduced_costs_nonnegative':True,'all24_matched_reduced_costs_zero':True,'cost_primal_dual_equal':sum(C[i][match[i]]for i in range(n)),'deterministic_augmenting_receipts24':steps,'zero_reduced_cost_row_graph':edges,'unique_minimum_matching':unique,'no_other_rowset_or_weight_search':True,'not_field_rank_or_determinant_coefficient_calculation':True}
 if not unique:
  alt=match[:]
  for i,j in zip(cycle[:-1],cycle[1:]):assert j in edges[i];alt[i]=match[j]
  assert sorted(alt)==list(range(n))and alt!=match and sum(C[i][alt[i]]for i in range(n))==sum(C[i][match[i]]for i in range(n))
  certificate.update({'alternating_cycle_rows_positions':cycle,'alternating_cycle_actual_rows':[rows[i]for i in cycle],'second_distinct_matching24':alt,'both_matching_costs_identical':True,'determinant_or_actual_pivots_nonzero_not_proved':True,'no_inverse_program_claim':True})
 else:
  topo=post[::-1];pos={i:k for k,i in enumerate(topo)};assert len(topo)==n and all(pos[i]<pos[j]for i in range(n)for j in edges[i]);certificate['zero_cost_DAG_topological_row_order']=topo
  # Actual K5 RF Gaussian program. No field coefficient is expanded, no entry
  # is discarded by a leading valuation; every initial nonzero is literal.
  nodes=[]
  def node(op,*args):nodes.append([op,*args]);return len(nodes)-1
  zero=node('QQ','0');one=node('QQ','1');A=[]
  for ii in topo:
   row=[]
   for jj in topo:
    j=match[jj];row.append(node('actual_entry_ref',by[ii,j]['source'])if(ii,j)in by else zero)
   A.append(row)
  rhs=[node('selected_transformed_RHS',rows[i])for i in topo];piv=[]
  for k in range(n):
   p=A[k][k];piv.append({'step':k,'actual_pivot_RF_node':p,'nonzero_proof':'Unique leading matching on each matched principal prefix; actual pivot=prefix minor ratio.'})
   for i in range(k+1,n):
    q=node('div',A[i][k],p)
    for j in range(k+1,n):A[i][j]=node('sub',A[i][j],node('mul',q,A[k][j]))
    rhs[i]=node('sub',rhs[i],node('mul',q,rhs[k]));A[i][k]=zero
  solution=[None]*n
  for k in range(n-1,-1,-1):
   vR=rhs[k]
   for j in range(k+1,n):vR=node('sub',vR,node('mul',A[k][j],solution[j]))
   solution[k]=node('div',vR,A[k][k])
  instruction={'all_passed':True,'actual_RF_nodes':nodes,'all24_nonzero_actual_pivots':piv,'solution_in_topological_matching_column_order':solution,'column_order':[cols[match[i]]for i in topo],'selected_actual_row_order':[rows[i]for i in topo],'all_actual_entries_retained_no_series_zero_drop':True,'no_scalar_inverse_expansion':True,'normal42_RHS_transport':'Apply frozen49 unit reconstruction; apply ten actual constraint particular solutions and remaining RHS subtraction. Solve selected24. Undo free-column scaling by old_x_free=P290*y_free exactly once, restore ten pivot variables by h_part+N_alias*y_free, then frozen49 back substitution. Return original83 columns41..82.','consistency_transport':'For arbitrary236 RHS retain all87 unselected top-row equations and all66 transformed zero-row RHS. This is a left inverse on the actual column image, not surjectivity onto all row RHS.','full177_cokernel':'coker(full177)=coker(top111) direct_sumK5^66; no zero-row summand deleted.','chart':'Old base monomial units*P290*detA*detB plus this fixed nonzero24minor. Complement remains open; no generic fiber exhaustion claim.'}
  certificate['lazy_inverse_operator']='actual24-lazy-inverse-and-RHS-transport.json';(out/certificate['lazy_inverse_operator']).write_bytes(enc(instruction));peak=max(peak,numeric(instruction)+numeric(certificate)+numeric(entries)+25000);assert peak<=b['whole_owned_numeric_upper']
 (out/'exact-fixed24-assignment-and-uniqueness.json').write_bytes(enc(certificate));raw=actual_read_bytes+sum(p.stat().st_size for p in out.rglob('*')if p.is_file());assert raw<=b['fresh_RAW_upper']
 report={'id':ID,'all_passed':True,'fixed24_selected':True,'actual_leading_entries_checked':len(leaders),'selected_actual_monomials_rechecked':24,'exact_dual_and_all_reduced_costs':True,'unique_minimum_matching':unique,'whole_owned_numeric_upper':b['whole_owned_numeric_upper'],'checked_owned_numeric_peak':peak,'actual_fresh_RAW_before_small_report':raw,'prepared_gate_READ_bytes':prepared_gate_bytes,'all_actual_READ_bytes':actual_read_bytes,'no_alternative_row_search':True,'no_actual_pivot_nonzero_or_inverse_claim':not unique,'input_sha256':{**b['actual_READ_input_sha256'],str(gp.relative_to(ROOT)):sha(gp),str(bp.relative_to(ROOT)):sha(bp)},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));assert actual_read_bytes+sum(p.stat().st_size for p in out.rglob('*')if p.is_file())<=b['fresh_RAW_upper'];signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':main()

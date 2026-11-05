#!/usr/bin/env python3
"""SOURCE ONLY: actual11x11, injective radix valuations, finite prefix cost gates."""
import argparse,hashlib,itertools,json,math,resource,signal,struct,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_TROPICAL_SOURCE_PREFLIGHT_20261105_v63'
FULLBOX=((5,66),(0,102),(0,32))
NEG=-10**15
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def selected_inputs(path,roots):
 wanted={x for r in roots for x in r}|set(range(575,593));by_source={};diff={};kept=skipped=0;inventory=[]
 with path.open('rb')as f:
  def get(fmt):
   n=struct.calcsize('<'+fmt);b=f.read(n);assert len(b)==n;return struct.unpack('<'+fmt,b)[0]
  def z(limit,keep):
   n=get('H');sgn=get('B');assert sgn in[0,1]and n<=(limit+7)//8
   if not keep:f.seek(n,1);return None
   b=f.read(n);assert len(b)==n and(not n or b[0]);v=int.from_bytes(b,'big');assert v.bit_length()<=limit;return -v if sgn else v
  assert f.read(8)==b'K92SP001';n=get('I');assert n<65536;header=json.loads(f.read(n));assert header['source_completed']==572;np=get('I');assert np<=1024
  for _ in range(np):
   pid=get('I');role=get('B');src=get('i');term=get('i');row=get('i');mask=get('i');done=get('B');retain=src in wanted or(role==2 and term==-20);num=z(4096,retain);den=z(4096,retain);nt=get('I');co=[]if retain else None
   for j in range(nt):
    m=(get('i'),get('i'),get('i'));v=z(1024,retain)
    if retain:assert v and(not co or co[-1][0]<m);co.append((m,v))
   if retain:
    assert done;p=dict(state_id=pid,source=src,role=role,term=term,row=row,mask=mask,c=Q(num,den),co=co);assert p['c']>0;kept+=nt;inventory.append({k:v for k,v in p.items()if k not in['co','c']})
    if src>=0:by_source[src]=p
    else:assert row==1 and mask in[0,1,2];diff[mask]=p
   else:skipped+=nt
  pending=get('B')
  if pending:
   for _ in range(3):get('i')
   for _ in range(2):get('Q')
   for _ in range(3):get('q')
   z(16384,False);get('B');z(16384,False)
  assert not f.read(1)
 assert len(by_source)==106 and len(diff)==3 and kept==21690
 zero=dict(state_id=-1,c=Q(1),co=[])
 M=[[by_source[x]for x in r]+[zero]*3 for r in roots]
 for r in range(6):
  for j in range(3):M[5+r][8+j]=diff[j]if r==1 else by_source[575+3*r+j]
 return M,inventory,kept,skipped
def assignment(V):
 dp={0:(0,1)}
 for r in range(11):
  nex={}
  for mask,(v,count)in dp.items():
   for j,vj in enumerate(V[r]):
    if mask>>j&1 or vj==NEG:continue
    nm=mask|1<<j;value=v+vj
    if nm not in nex or value>nex[nm][0]:nex[nm]=(value,count)
    elif value==nex[nm][0]:nex[nm]=(value,count+nex[nm][1])
  dp=nex
 return dp[2047]
def dual_potentials(V):
 n=11;INF=10**15;u=[0]*(n+1);v=[0]*(n+1);p=[0]*(n+1);way=[0]*(n+1)
 for i in range(1,n+1):
  p[0]=i;j0=0;mins=[INF]*(n+1);used=[False]*(n+1)
  while True:
   used[j0]=True;i0=p[j0];delta=INF;j1=0
   for j in range(1,n+1):
    if used[j]:continue
    cost=INF//2 if V[i0-1][j-1]==NEG else -V[i0-1][j-1];cur=cost-u[i0]-v[j]
    if cur<mins[j]:mins[j]=cur;way[j]=j0
    if mins[j]<delta:delta=mins[j];j1=j
   for j in range(n+1):
    if used[j]:u[p[j]]+=delta;v[j]-=delta
    else:mins[j]-=delta
   j0=j1
   if not p[j0]:break
  while True:
   j1=way[j0];p[j0]=p[j1];j0=j1
   if not j0:break
 r=[-u[i]for i in range(1,n+1)];c=[-v[j]for j in range(1,n+1)];matching=[None]*n
 for j in range(1,n+1):matching[p[j]-1]=j-1
 assert all(V[i][matching[i]]!=NEG and r[i]+c[matching[i]]==V[i][matching[i]]for i in range(n))
 assert all(V[i][j]==NEG or r[i]+c[j]>=V[i][j]for i in range(n)for j in range(n))
 return r,c,matching
def radix(order):
 w=[0]*3;stride=1
 for j in order:w[j]=stride;stride*=FULLBOX[j][1]-FULLBOX[j][0]+1
 return tuple(w)
def tight_Newton_max(w,meta):
 best=None;point=None
 for a in range(FULLBOX[0][0],FULLBOX[0][1]+1):
  for e in range(FULLBOX[1][0],FULLBOX[1][1]+1):
   lo,hi=FULLBOX[2]
   for weight,(lower,upper)in zip(meta['weights'],meta['output_metadata']['linear']):
    base=weight[0]*a+weight[1]*e
    if weight[2]:lo=max(lo,(lower-base+weight[2]-1)//weight[2]);hi=min(hi,(upper-base)//weight[2])
    elif not lower<=base<=upper:hi=lo-1;break
   if lo<=hi:
    q=(a,e,hi);value=sum(x*y for x,y in zip(q,w))
    if best is None or value>best:best=value;point=q
 return best,point
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 b=src('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json');roots=b['actual_pure_matrix_roots'];path=ROOT/'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin';pins[str(path.relative_to(ROOT))]=sha(path);assert pins[str(path.relative_to(ROOT))]=='cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02'
 M,inventory,kept,skipped=selected_inputs(path,roots)
 wb=src('reports/S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63/weighted-output.json');norm=src('reports/S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63/normalization.json');assert tuple(map(tuple,wb['output_metadata']['box']))==FULLBOX and wb['row_known_a_factors']==[0,1,1,0,0,0,0,0,0,0,0]
 # Exactsource references rather than new determinant values.
 src('reports/S61_NONSQUARE_PREPOWER_K92_H3_SIGN_SCHEDULE_20261104_v63/report.json');src('reports/S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63/report.json')
 denom_punit=all(p['c'].denominator%32003 for row in M for p in row);assert denom_punit
 dump(out/'source-binding.json',dict(state_sha256=pins[str(path.relative_to(ROOT))],matrix11x11_state_ids=[[p['state_id']for p in row]for row in M],matrix11x8_old_DAG_roots=roots,
  physical_targets=[22,11,13,25,15,19,14,18,28,30,29],row_order=['Q9','G18','G39','G68','G20','A=G115/c10','B=(G120-G99)/c11','D=G92/c8','C=G85/c6','E=aG129/c12','F=G99/c9'],
  H3_first5rows_literal_zero=True,H3_lower6=['N10','N11-N9 savedliteral','N8','N6','N12','N9'],H3_column_normalization='divide known FULLU aFC²',P4_rows1and2_a_normalization='already in saved pure roots; H3 here is literalzero',det_normalization='original fixed det11=a² * determinant(this11x11); existingglobalQQ/a5 normalization separate',
  selected_polynomials=109,selected_terms=kept,skipped_terms_without_allocation=skipped,all_literal_QQdenominators_p32003units=denom_punit,inventory=inventory,NO_source572replay=True,NO_new_polynomial_coefficients=True))
 profiles=[]
 for order in itertools.permutations(range(3)):
  w=radix(order);V=[];leading=[]
  for row in M:
   vals=[];lead=[]
   for p in row:
    if not p['co']:vals.append(NEG);lead.append(None);continue
    values=[sum(x*y for x,y in zip(mon,w))for mon,v in p['co']];mx=max(values);idx=values.index(mx);assert values.count(mx)==1;mon,coef=p['co'][idx];vals.append(mx);lead.append(dict(state_id=p['state_id'],source_coefficient_index=idx,monomial=mon,integer_coefficient=str(coef),common_QQfactor=str(p['c'])))
   V.append(vals);leading.append(lead)
  degree,count=assignment(V);r,c,match=dual_potentials(V);assert sum(r)+sum(c)==degree;maximum,point=tight_Newton_max(w,wb);boxmax=sum(w[i]*FULLBOX[i][1]for i in range(3));Tmat=[]
  for i in range(11):Tmat.append([leading[i][j]if V[i][j]!=NEG and r[i]+c[j]==V[i][j]else None for j in range(11)])
  layer_cost=[]
  for K in[0,8,16,32]:
   entry_terms=0;rowbits=[];rows=[]
   for i,row in enumerate(M):
    selected=[];D=1
    for j,p in enumerate(row):
     for idx,(mon,v)in enumerate(p['co']):
      gap=r[i]+c[j]-sum(x*y for x,y in zip(mon,w));assert gap>=0
      if gap<=K:selected.append((j,idx,gap,v,p['c']));D=math.lcm(D,p['c'].denominator)
    L=sum(abs(v)*factor.numerator*(D//factor.denominator)for j,idx,gap,v,factor in selected);rowbits.append(L.bit_length());entry_terms+=len(selected);rows.append(dict(row=i,terms=len(selected),row_denominator_bits=D.bit_length(),row_integer_L1_bits=L.bit_length()))
   max_DP_cells=924*(K+1);matrix_dense=121*(K+1);integer_bits=sum(rowbits);raw=max_DP_cells*((integer_bits+7)//8+16)+entry_terms*((max(rowbits)+7)//8+40)+2*1024*1024
   layer_cost.append(dict(K=K,selected_existing_source_terms=entry_terms,zero_layers_FORCED_by_knownNewton_bound=max(0,degree-maximum),known_doomed_for_this_K=degree-K>maximum,
    QQinteger_coefficient_bit_upper=integer_bits,QQDP_coefficient_cells_upper=max_DP_cells,matrix_dense_slots_upper=matrix_dense,total_with21690source_and15000scratch=max_DP_cells+matrix_dense+21690+15000,QQprefix_RAW_upper=raw,
    Fp32003_DP_bytes_upper=(max_DP_cells+matrix_dense)*8+2*1024*1024,rows=rows))
  profiles.append(dict(variable_low_to_high_order=order,weight=w,injective_on_full_known_box=True,valuation_matrix=V,tropical_degree=degree,optimal_assignment_count=count,row_potentials=r,column_potentials=c,matching=match,
   dual_objective=sum(r)+sum(c),leading_tight_matrix_existing_coefficients=Tmat,known_full_Newton_max=maximum,attaining_upper_lattice_point=point,known_full_box_max=boxmax,forced_initial_zero_layers=max(0,degree-maximum),prefix_costs=layer_cost,NO_QQ_or_Fp_determinant_coefficient_evaluated=True))
  dump(out/'phase.json',dict(phase='SOURCE ONLY injective-radix support/assignment metadata',orders_profiled=len(profiles)))
 dump(out/'six-support-orders.json',dict(profiles=profiles,six_variable_orders_are_metadata_ONLY=True,no_weight_family_coefficient_probe=True))
 best=min(profiles,key=lambda p:(p['forced_initial_zero_layers'],p['weight']));feasible=[p for p in profiles if p['forced_initial_zero_layers']<=32]
 decision=dict(fixed_prefix_K_max_preflight=32,best_order=best['variable_low_to_high_order'],best_weight=best['weight'],best_tropical_degree=best['tropical_degree'],best_knownNewton_max=best['known_full_Newton_max'],best_forced_zero_layers=best['forced_initial_zero_layers'],
  any_K32_not_a_priori_doomed=bool(feasible),proposed_scalar_job='NONE if everyorderforcedgap>32; otherwise ROOTmayselectONEfixedweight/K andQQorFp method',
  QQ_method='row-clear denominators; truncated integer subset-DP at fixedK. No rationalfunctionGaussian or parameterpoint.',Fp_method='existing32003 only; all literalQQdenominators units. Same exactformalcoefficient extraction. Any NONZERO Fp coefficient proveschar0NONZERO, not fullQQdetcoefficient identity.',
  fullfield_and_chart_scope='Only fixed det11 nonzero if positive; DDelta407*Saux*T*N6*N7bar and all original complement gaps remain. No global57 or fullresidualrank2.',NO_scalar_job=True)
 dump(out/'decision.json',decision);note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));total=sum(p.stat().st_size for p in out.iterdir());assert kept+15000<150000 and total<12*1024*1024
 dump(out/'phase.json',dict(phase='FINAL SOURCE ONLY tropical gate; no determinant coefficient work'))
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,selected_source_terms=kept,no_source572_replay=True,no_new_polynomial_or_determinant_coefficients=True,all_literal_QQdenominators_p32003units=denom_punit,
  support_order_profiles=6,best_fixed_weight=best['weight'],best_forced_zero_layers=best['forced_initial_zero_layers'],any_K32_not_a_priori_doomed=bool(feasible),no_coeff_job_authorized=True,
  total_RAW_before_report=total,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

#!/usr/bin/env python3
"""SOURCE ONLY: support-certified Schur gap closure and first finite leading block."""
import argparse,hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_SCHUR_GAP_SOURCE_PREFLIGHT_20261105_v63';W=(1,2244,68);INF=10**12
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def min_assignment(V):
 n=len(V);dp={0:(0,1)}
 for row in V:
  nex={}
  for mask,(v,ct)in dp.items():
   for j,c in enumerate(row):
    if mask>>j&1 or c>=INF:continue
    nm=mask|1<<j;x=v+c
    if nm not in nex or x<nex[nm][0]:nex[nm]=(x,ct)
    elif x==nex[nm][0]:nex[nm]=(x,ct+nex[nm][1])
  dp=nex
 return dp[(1<<n)-1]
def dual_min(V):
 n=len(V);u=[0]*(n+1);v=[0]*(n+1);p=[0]*(n+1);way=[0]*(n+1)
 for i in range(1,n+1):
  p[0]=i;j0=0;mins=[INF]*(n+1);used=[False]*(n+1)
  while True:
   used[j0]=True;i0=p[j0];delta=INF;j1=0
   for j in range(1,n+1):
    if used[j]:continue
    cur=V[i0-1][j-1]-u[i0]-v[j]
    if cur<mins[j]:mins[j]=cur;way[j]=j0
    if mins[j]<delta:delta=mins[j];j1=j
   assert delta<INF
   for j in range(n+1):
    if used[j]:u[p[j]]+=delta;v[j]-=delta
    else:mins[j]-=delta
   j0=j1
   if not p[j0]:break
  while True:
   j1=way[j0];p[j0]=p[j1];j0=j1
   if not j0:break
 a=u[1:];b=v[1:];match=[None]*n
 for j in range(1,n+1):match[p[j]-1]=j-1
 assert all(a[i]+b[j]<=V[i][j]for i in range(n)for j in range(n))and all(a[i]+b[match[i]]==V[i][match[i]]for i in range(n))
 return a,b,match
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 base='reports/S61_NONSQUARE_PREPOWER_K92_ONE_ROW_PRECONDITION_20261105_v63';rpt=src(base+'/report.json');assert rpt['all_passed']and pins[base+'/report.json']=='d9d474f281598585255e91984dbf0e0a1f9f4c208d809d4bd7c3472dcc817a68';x=src(base+'/retropicalized-source.json');row=src(base+'/actual-row11.json');deltagate=src(base+'/next-delta-block-source-gate.json')
 r=x['new_leading_constant_rank'];assert r==4 and x['smallest_positive_gap']==1;R=[[Q(v)for v in rr]for rr in x['row_transform']];C=[[Q(v)for v in rr]for rr in x['column_transform']];assert all(C[1][j]==0 for j in range(r,11));assert len(deltagate['physical_delta_source_records'])==1 and deltagate['physical_delta_source_records'][0]['column']==1
 dump(out/'gap1-support-zero.json',dict(physical_gap1_entries=[[1,1]],Q_row1_residual_columns=[str(C[1][j])for j in range(r,11)],B1_and_D1_identically_zero_by_support=True,Schur_coefficient1_identically_zero=True,NO_coefficient_products_or_7det=True,targetweight232220_outside_knownNewton=deltagate['monomial_in_known_Newton_bounds']is False))
 helper=ROOT/'code/s61_nonsquare_prePOWER_K92_tropical_source_preflight_v2_20261105_v63.py';pins[str(helper.relative_to(ROOT))]=sha(helper);sp=importlib.util.spec_from_file_location('sourcegate',helper);h=importlib.util.module_from_spec(sp);sp.loader.exec_module(h);bind=src('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json');state=ROOT/'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin';pins[str(state.relative_to(ROOT))]=sha(state);assert pins[str(state.relative_to(ROOT))]=='cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02';M,inventory,kept,skip=h.selected_inputs(state,bind['actual_pure_matrix_roots'])
 for p in row['columns']:M[1][p['column']]=dict(state_id=-1001-p['column'],c=Q(p['QQ_common_factor']),co=[(tuple(m),int(v))for m,v in p['primitive_integer_polynomial']])
 rp=x['row_potentials'];cp=x['column_potentials'];Ds=list(map(int,x['constant_clearing_rowLCMs']));F=[[INF]*11 for i in range(11)];physical_min=[[INF]*11 for i in range(11)];gaprecords={}
 for i,rr in enumerate(M):
  for j,p in enumerate(rr):
   by_gap={}
   for idx,(mon,v)in enumerate(p['co']):
    gap=rp[i]+cp[j]-sum(a*b for a,b in zip(mon,W));assert gap>=0 and gap not in by_gap;by_gap[gap]=dict(source_index=idx,monomial=mon,integer=v)
   gaprecords[i,j]=by_gap;g=min((g for g in by_gap if g>0),default=INF);physical_min[i][j]=g
   for u in range(11):
    if R[u][i]:
     for v in range(11):
      if C[j][v]:F[u][v]=min(F[u][v],g)
 A=[rr[:r]for rr in F[:r]];B=[rr[r:]for rr in F[:r]];CC=[rr[:r]for rr in F[r:]];D=[rr[r:]for rr in F[r:]];star=[[0 if i==j else A[i][j]for j in range(r)]for i in range(r)]
 for k in range(r):
  for i in range(r):
   for j in range(r):star[i][j]=min(star[i][j],star[i][k]+star[k][j])
 n=11-r;route=[[min(CC[i][k]+star[k][l]+B[l][j]for k in range(r)for l in range(r))for j in range(n)]for i in range(n)];S=[[min(D[i][j],route[i][j])for j in range(n)]for i in range(n)];order,count=min_assignment(S);alpha,beta,matching=dual_min(S);assert sum(alpha)+sum(beta)==order
 tight=[(i,j)for i in range(n)for j in range(n)if S[i][j]<INF and alpha[i]+beta[j]==S[i][j]];assert all(D[i][j]==S[i][j]and route[i][j]>S[i][j]for i,j in tight)
 DR=math.lcm(*(v.denominator for rr in R for v in rr));DQ=math.lcm(*(v.denominator for rr in C for v in rr));RZ=[[int(v*DR)for v in rr]for rr in R];QZ=[[int(v*DQ)for v in rr]for rr in C];physical={};requests=[];scalars=0;maxnum=0
 for i,j in tight:
  u=r+i;v=r+j;g=S[i][j];contributors=[]
  for a in range(11):
   if not RZ[u][a]:continue
   for b in range(11):
    if not QZ[b][v]or g not in gaprecords[a,b]:continue
    atom=gaprecords[a,b][g];p=M[a][b];z=p['c']*atom['integer']*Ds[a];assert z.denominator==1;key=(a,b,g)
    if key not in physical:physical[key]=dict(id=len(physical),row=a,column=b,gap=g,source_state_id=p['state_id'],source_index=atom['source_index'],monomial=atom['monomial'],source_integer=str(atom['integer']),source_QQfactor=str(p['c']),row_LCM=str(Ds[a]),cleared_integer=str(z.numerator))
    contributors.append(dict(physical_id=physical[key]['id'],R_cleared_coefficient=str(RZ[u][a]),Q_cleared_coefficient=str(QZ[b][v])));scalars+=2
  assert contributors;bound=sum(abs(int(c['R_cleared_coefficient']))*abs(int(physical[next(k for k,p in physical.items()if p['id']==c['physical_id'])]['cleared_integer']))*abs(int(c['Q_cleared_coefficient']))for c in contributors);maxnum=max(maxnum,bound)
  requests.append(dict(residual_row=i,residual_column=j,global_row=u,global_column=v,gap=g,contributors=contributors,integer_numerator_abs_bound=str(bound)))
 degree=x['new_tropical_degree'];target=degree-order;z=target-5;mon=(5+z%68,(z//68)//33,(z//68)%33);wb=src('reports/S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63/weighted-output.json');inside=5<=mon[0]<=72 and 0<=mon[1]<=107 and 0<=mon[2]<=32 and all(lo<=sum(a*b for a,b in zip(mon,w))<=hi for w,(lo,hi)in zip(wb['weights'],wb['output_metadata']['linear']));assert inside
 dump(out/'minplus-proof.json',dict(INF=INF,physical_min_positive_gaps=physical_min,transformed_min_positive_gaps=F,Aplus_min=A,B_min=B,C_min=CC,D_min=D,Ainverse_minplus_closure=star,correction_C_Ainv_B_min=route,Schur_min=S,
  first_possible_determinant_gap=order,minimum_assignment_count=count,row_dual_potentials=alpha,column_dual_potentials=beta,matching=matching,all_tight_entries=tight,all_tight_corrections_STRICTLY_later=True,NO_scalar_Schur_coefficients_computed=True))
 numerator_bits=maxnum.bit_length();det_bits=n*numerator_bits+math.factorial(n).bit_length();raw_upper=(len(physical)+len(requests)+121+49)*(det_bits//3+300)+500000;live=kept+row['actual_total_terms']+len(physical)+len(requests)+121+49+15000
 dump(out/'finite-leading7-source.json',dict(weight=W,rank4_normalization=True,formal_matrix_B='rowLCM*t^rowpot*M(t^-w)*t^colpot',tropical_degree=degree,first_possible_determinant_gap=order,
  first_possible_original_weight=target,unique_original_monomial=mon,monomial_in_knownNewton=True,leading7_rows='residualindices0..6 corresponding transformedsource rows4..10',leading7_columns='residualindices0..6 corresponding transformedtargetcols4..10',
  leading7_entry='(R*B_gap*Q)[4+i,4+j] only on tight minassignment edges; all other entriesZERO',correction_C_Ainv_B_STRICTLY_later_on_all_tight_edges=True,common_denominator=str(DR*DQ),R_common_denominator=str(DR),Q_common_denominator=str(DQ),
  physical_atoms=list(physical.values()),entry_recipes=requests,entry_recipe_is_exact_finite_integer_sum=True,all_other49_entries_zero_by_dual_strict_inequality=True,
  determinant_recipe='detB coefficient atgap13255 = det(leading7)/(detR*detQ); detA(0)=1. OriginaldetM weightedcoefficient divides rowLCMproduct aswell.',
  NO_leading7_entries_OR_determinant_evaluated=True,positive_only_scope='If det7 nonzero then fixedactualdet11 nonzero; ifzero onlyleadinglayerzero, no fullzero/global57 inference',
  max_integer_entry_numerator_bits=numerator_bits,max_integer_det_bits=det_bits,scalar_products_for_all_entries=scalars,artifact_RAW_upper=raw_upper,whole_source_and_live_upper=live,
  next_proposed_JOB='ONE allfinite integer entry sums+exact7x7det and rank certificate; no seriesexpansion, no source572, no otherweight/K/parameterpoint/prime. RequiresROOTbudget.'))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);total=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert live<150000 and raw_upper+external<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,gap1_Schur_zero_proved_by_support=True,NO_gap1_coefficient_job=True,rank4_fixed_constant_normalization=True,residual_size=7,Schur_determinant_order_lower_bound=order,
  first_possible_original_monomial=mon,all_tight_corrections_strictly_later=True,physical_coefficient_atoms=len(physical),finite_integer_entry_requests=len(requests),finite_scalar_products=scalars,entry_numerator_bit_bound=numerator_bits,determinant_bit_bound=det_bits,
  whole_live_upper=live,next_JOB_RAW_upper=raw_upper,no_leading7_or_determinant_coefficients_evaluated=True,no_det11_nonzero_or_global57_claim=True,total_RAW_before_report=total,external_read_RAW=external,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

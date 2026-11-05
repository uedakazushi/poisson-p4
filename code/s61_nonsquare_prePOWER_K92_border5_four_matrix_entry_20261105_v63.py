#!/usr/bin/env python3
"""ONE integer four-matrix coefficient; run only under explicit ROOT allocation."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_BORDER5_FOUR_MATRIX_ENTRY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={};peak=0
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 base='reports/S61_NONSQUARE_PREPOWER_K92_BORDER5_FOUR_MATRIX_SOURCE_20261105_v63';r=src(base+'/report.json');assert r['all_passed']and pins[base+'/report.json']=='67be63d38c53ce2b79b3fd11f161933d222649fa82cd0921c263857d020de4ad';gate=src(base+'/resource-gate.json');plan=src(base+'/polynomial-border5-source.json');fixed=src(base+'/four-matrix-coefficient-proof.json');receipt=src('reports/S61_NONSQUARE_PREPOWER_K92_BORDER5_FOUR_MATRIX_WORKER_GATE_20261105_v63/report.json');assert receipt['all_passed']
 assert plan['whole_coefficient_order']==26794 and plan['required_border5_coefficient_degree']==141 and plan['maximum_Hbar_gap']==138;cap=138
 L=int(plan['source_uniform_denominator']);Dm=int(plan['DA_denominator']);HD=int(plan['Hbar_denominator']);TD=int(plan['Tbar_denominator']);Rd=int(plan['R2_denominator']);Qd=int(plan['Q2_denominator']);assert HD==L**5*Dm and TD==HD*Rd*Qd
 R=[[int(Q(x)*Rd)for x in rr]for rr in plan['R2_normalized']];QQ=[[int(Q(x)*Qd)for x in rr]for rr in plan['Q2']];DA={int(g):int(Q(x)*Dm)for g,x in plan['DA'].items()};I={};source_receipts=[];inputslots=0
 typ=src('reports/S61_NONSQUARE_PREPOWER_K92_FACTOR_RANK_PROFILE_20261105_v63/typed-next-residual-Schur.json');alpha=typ['alpha'];beta=typ['beta'];closed=src('reports/S61_NONSQUARE_PREPOWER_K92_RESIDUAL6_FINITE_SOURCE_20261105_v63/closed-inverse-source-proof.json')
 for rel in['reports/S61_NONSQUARE_PREPOWER_K92_EXACT_ABC72_20261105_v63/ABC72.jsonl','reports/S61_NONSQUARE_PREPOWER_K92_EXACT_D49_20261105_v63/D49.jsonl']:
  p=ROOT/rel;pins[rel]=sha(p)
  with p.open()as f:
   for ln,line in enumerate(f):
    e=json.loads(line);b=e.get('block','D');i=e['block_row']if b!='D'else e['residual_row'];j=e['block_column']if b!='D'else e['residual_column'];c=Q(e['QQ_content']);assert L%c.denominator==0;factor=c.numerator*(L//c.denominator);pp={int(g):int(v)*factor for g,v in e['sparse_integer_polynomial']};I[b,i,j]=pp;inputslots+=len(pp);source_receipts.append(dict(block=b,row=i,column=j,path=rel,line=ln,entry_index=e['entry_index'],QQ_content=str(c),integer_scaling_factor=str(factor),coefficient_terms=len(pp)))
 assert inputslots==23932
 HB={};TB={};Abar={};DP={};matrices={};cofs={};counters={'path_pairs':0,'H_contributions':0,'cofactor_products':0,'requested_border_products':0,'DP_products':0}
 def guard(scratch=0):
  nonlocal peak
  actual=inputslots+sum(map(len,Abar.values()))+sum(map(len,HB.values()))+sum(map(len,TB.values()))+sum(map(len,cofs.values()))+sum(map(len,DP.values()))+sum(len(m)*len(m[0])for m in matrices.values())+scratch+15000;peak=max(peak,actual)
  assert actual<=receipt['whole_owned_scalar_upper']<=150000
  if time.monotonic()-start>45:raise RuntimeError('45s_guard')
 def mul(p,q,limit,kind):
  outp={}
  for a,x in p.items():
   for b,y in q.items():
    if a+b<=limit:outp[a+b]=outp.get(a+b,0)+x*y;counters[kind]+=1
   guard(len(outp)+1)
  return{k:v for k,v in outp.items()if v}
 for rec in plan['finite_simple_path_proof']:
  k=rec['start'];l=rec['end'];pp={};limit=closed['Ainverse_minplus'][k][l]+cap
  for path in rec['paths']:
   v={0:1}
   for a,b in path['edges']:v=mul(v,I['A',a,b],limit,'path_pairs')
   v=mul(v,{0:Dm}if path['through1']else DA,limit,'path_pairs');scale=path['sign']*L**(3-len(path['edges']))
   for g,x in v.items():pp[g]=pp.get(g,0)+scale*x
  Abar[k,l]={g:x for g,x in pp.items()if x};guard(len(pp)*2)
 for i in range(7):
  for j in range(7):
   shift=alpha[i]+beta[j];pp={}
   for d,dv in I['D',i,j].items():
    for z,zv in DA.items():
     g=d+z-shift
     if 0<=g<=cap:pp[g]=pp.get(g,0)+L**4*dv*zv;counters['H_contributions']+=1
   for k in range(4):
    for l in range(4):
     A=Abar[k,l];C=I['C',i,k];B=I['B',l,j]
     if not A or not C or not B or min(C)+min(A)+min(B)-shift>cap:continue
     for c,cv in C.items():
      if c+min(A)+min(B)-shift>cap:continue
      for b,bv in B.items():
       if c+b+min(A)-shift>cap:continue
       for a,av in A.items():
        g=c+a+b-shift
        if 0<=g<=cap:pp[g]=pp.get(g,0)-cv*av*bv;counters['H_contributions']+=1
   HB[i,j]={g:x for g,x in pp.items()if x};guard(len(pp)*2)
 for u in range(7):
  for v in range(7):
   pp={}
   for i in range(7):
    for j in range(7):
     if R[u][i]and QQ[j][v]:
      for g,x in HB[i,j].items():pp[g]=pp.get(g,0)+R[u][i]*QQ[j][v]*x
   TB[u,v]={g:x for g,x in pp.items()if x};assert TB[u,v].get(0,0)==(TD if u==v==0 else 0);guard(len(pp)*2)
 dump(out/'literal-source-and-all-finite-path-identities.json',dict(source=source_receipts,common_source_denominator=str(L),D_A_integer=DA,Hbar_common_denominator=str(HD),Tbar_common_denominator=str(TD),all_finite_Abar=[dict(row=i,column=j,integer_polynomial=sorted(pp.items()))for(i,j),pp in Abar.items()],
  all_Hbar=[dict(row=i,column=j,integer_polynomial=sorted(pp.items()))for(i,j),pp in HB.items()],all_Tbar=[dict(row=i,column=j,integer_polynomial=sorted(pp.items()))for(i,j),pp in TB.items()],source_formula_and_signed_paths_bound_to_plan=True))
 # Release source/frame caches before the fixed four-matrix extraction.
 I.clear();Abar.clear();HB.clear();inputslots=0
 def diff(p,q):
  z=dict(p)
  for k,v in q.items():z[k]=z.get(k,0)-v
  return{k:v for k,v in z.items()if v}
 cofs['pivot2']=diff(mul(TB[0,0],TB[6,6],cap,'cofactor_products'),mul(TB[0,6],TB[6,0],cap,'cofactor_products'))
 for q in plan['freeQ']:
  j=1+q;cofs['A'+str(q)]=diff(mul(TB[0,6],TB[6,j],cap,'cofactor_products'),mul(TB[6,6],TB[0,j],cap,'cofactor_products'));cofs['B'+str(q)]=diff(mul(TB[6,0],TB[0,j],cap,'cofactor_products'),mul(TB[0,0],TB[6,j],cap,'cofactor_products'))
 c=Q(plan['source_secondpivot_leading']);assert cofs['pivot2'].get(1,0)==int(c*TD**2)
 def coefficient(p,q,k):
  z=0
  for a,x in p.items():
   if k-a in q:z+=x*q[k-a];counters['requested_border_products']+=1
  return z
 roots=[]
 for d in[2,3,67,68]:matrices[d]=[[0]*5 for i in range(5)]
 for rec in plan['entries']:
  u=rec['row'];v=rec['column'];i=rec['Tbar_row'];j=rec['Tbar_column'];q=plan['freeQ'][v];shift=rec['removed_t_power'];bound={int(g):int(n)for g,n in rec['needed_coefficient_upper']}
  for d in[0,1,2,3,67,68]:
   k=shift+d;terms=[coefficient(cofs['pivot2'],TB[i,j],k),coefficient(cofs['A'+str(q)],TB[i,0],k),coefficient(cofs['B'+str(q)],TB[i,6],k)];z=sum(terms)
   if d<2:assert z==0
   else:assert abs(z)<=bound.get(d,0);matrices[d][u][v]=z
   roots.append(dict(row=u,column=v,degree=d,Tbar_absolute_degree=k,three_exact_signed_integer_terms=list(map(str,terms)),integer_coefficient=str(z),common_denominator=str(TD**3),previous_ZERO_check=d<2))
 dump(out/'all100-four-matrix-and50-knownzero-identities.json',dict(matrices={str(k):v for k,v in matrices.items()},common_denominator=str(TD**3),all150_roots=roots,exact25_previous_coefficients0and1_ZERO=True,
  all11_cofactor_integer_polynomials={k:sorted(p.items())for k,p in cofs.items()}))
 DP={0:{(0,0):1}};stages=[]
 for i in range(5):
  nxt={}
  for mask,pp in DP.items():
   for j in range(5):
    if mask>>j&1:continue
    sign=(-1)**((mask>>(j+1)).bit_count());qq=nxt.setdefault(mask|1<<j,{})
    for (u,t),x in pp.items():
     for (du,dt),d in[((0,0),2),((0,1),3),((1,0),67),((1,1),68)]:
      z=matrices[d][i][j]
      if u+du<=2 and t+dt<=1 and z:qq[u+du,t+dt]=qq.get((u+du,t+dt),0)+sign*x*z;counters['DP_products']+=1
   guard(sum(map(len,nxt.values()))+100)
  DP={m:{k:v for k,v in p.items()if v}for m,p in nxt.items()};stages.append(dict(row_completed=i+1,all_states=[dict(mask=m,coefficients=[[u,t,str(v)]for(u,t),v in sorted(p.items())])for m,p in DP.items()]))
 target=DP.get(31,{}).get((2,1),0);assert abs(target)<=int(fixed['coefficient_abs_upper']);factor=Q(plan['original_detM_target_from_integer_det5_target_factor']);actual=target*factor
 dump(out/'all-fixed-bivariate-DP-identities.json',dict(variable_order=['u','t'],bounds=[2,1],all5_stages=stages,integer_coefficient_u2t1=str(target),all_signs='(-1)^count(previousselectedcolumns>j)',common_denominator=str(TD**15),
  exact_original_detM_a5e91d18_coefficient=str(actual),exact_coefficient_normalization_factor=str(factor),original_candidate_det11='a²*detM'))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert new+external<33554432 and new<gate['new_RAW_upper']
 report=dict(id=ID,all_passed=True,frozen=True,ONE_finite_four_matrix_job=True,only_entry_degrees=[2,3,67,68],all_source23932_and_paths_bound=True,all100_integer_values_verified=True,all50_previous_ZERO_checks=True,
  exact_target_bidegree=[2,1],exact_target_whole_order=26794,original_detM_monomial=[5,91,18],integer_target_coefficient=str(target),original_QQ_target_coefficient=str(actual),target_NONZERO=bool(target),
  source_operations=counters,whole_owned_scalar_peak=peak,no_full_det11_ZERO_or_global57_claim=True,no_new_weight_window_point_prime_or_source572=True,
  total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256','integer_target_coefficient','original_QQ_target_coefficient']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

#!/usr/bin/env python3
"""ONE integer four-matrix coefficient; run only under explicit ROOT allocation."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_BORDER4_DIRECT_ENTRY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={};peak=0
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 base='reports/S61_NONSQUARE_PREPOWER_K92_BORDER4_DIRECT_SOURCE_20261105_v63';r=src(base+'/report.json');assert r['all_passed']and pins[base+'/report.json']=='442bfad7b7c0f6cd7a8ab7f9ea07ec2739f5fa7e3b27d08fc1ab913690a17c61';gate=src(base+'/resource-gate.json');plan=src(base+'/direct-border4-source.json');receipt=src('reports/S61_NONSQUARE_PREPOWER_K92_BORDER4_DIRECT_WORKER_GATE_20261105_v63/report.json');assert receipt['all_passed']
 assert plan['whole_coefficient_order']==27338 and plan['required_border4_coefficient_degree']==345 and plan['maximum_Hbar_gap']==550;cap=550
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
  actual=inputslots+sum(map(len,Abar.values()))+sum(map(len,HB.values()))+sum(map(len,TB.values()))+sum(map(len,cofs.values()))+sum(map(len,DP.values()))+sum(map(len,matrices.values()))+scratch+15000;peak=max(peak,actual)
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
 oldprefix=src('reports/S61_NONSQUARE_PREPOWER_K92_BORDER5_FOUR_MATRIX_ENTRY_20261105_v63/literal-source-and-all-finite-path-identities.json');overlap=[]
 for k,entries,table in [('Hbar',oldprefix['all_Hbar'],HB),('Tbar',oldprefix['all_Tbar'],TB)]:
  for e in entries:
   pair=e['row'],e['column'];oldp={int(g):int(v)for g,v in e['integer_polynomial']};newp={g:v for g,v in table[pair].items()if g<=138};assert oldp==newp;overlap.append(dict(block=k,row=pair[0],column=pair[1],all_saved_prefix_coefficients_match=True))
 dump(out/'all98-old-prefix-overlap.json',dict(overlap=overlap,all_complete_polynomial_prefixes_through138_equal=True))
 dump(out/'literal-source-and-all-finite-path-identities.json',dict(source=source_receipts,common_source_denominator=str(L),D_A_integer=DA,Hbar_common_denominator=str(HD),Tbar_common_denominator=str(TD),all_finite_Abar=[dict(row=i,column=j,integer_polynomial=sorted(pp.items()))for(i,j),pp in Abar.items()],
  all_Hbar=[dict(row=i,column=j,integer_polynomial=sorted(pp.items()))for(i,j),pp in HB.items()],all_Tbar=[dict(row=i,column=j,integer_polynomial=sorted(pp.items()))for(i,j),pp in TB.items()],source_formula_and_signed_paths_bound_to_plan=True))
 # Release source/frame caches before the fixed four-matrix extraction.
 I.clear();Abar.clear();HB.clear();inputslots=0
 def diff(p,q):
  z=dict(p)
  for k,v in q.items():z[k]=z.get(k,0)-v
  return{k:v for k,v in z.items()if v}
 P=plan['pivot_rows'];assert P==plan['pivot_columns']==[0,6,1]
 for i in range(3):
  for j in range(3):
   rr=[P[k]for k in range(3)if k!=i];cc=[P[k]for k in range(3)if k!=j]
   q=diff(mul(TB[rr[0],cc[0]],TB[rr[1],cc[1]],cap,'cofactor_products'),mul(TB[rr[0],cc[1]],TB[rr[1],cc[0]],cap,'cofactor_products'));cofs[i,j]={g:(-1)**(i+j)*v for g,v in q.items()}
 detp={}
 for j in range(3):
  for g,v in mul(TB[P[0],P[j]],cofs[0,j],cap,'cofactor_products').items():detp[g]=detp.get(g,0)+v
 detp={g:v for g,v in detp.items()if v};assert min(detp)==71 and detp[71]==int(Q(plan['pivot_determinant_leading'])*TD**3);cofs['pivot3']=detp
 roots=[]
 for rec in plan['entries']:
  u=rec['row'];v=rec['column'];i=rec['Tbar_row'];j=rec['Tbar_column'];shift=rec['removed_t_power'];bound={int(g):int(n)for g,n in rec['needed_coefficient_upper']}
  q=mul(detp,TB[i,j],cap,'requested_border_products');termhashes=[]
  for k in range(3):
   for l in range(3):
    z=mul(mul(TB[i,P[k]],cofs[l,k],cap,'requested_border_products'),TB[P[l],j],cap,'requested_border_products');q=diff(q,z);termhashes.append(dict(left_pivot_index=k,right_pivot_index=l,cofactor_index=[l,k],term_coefficient_count=len(z)))
  assert all(not x for g,x in q.items()if g<=shift);tr={g-shift:x for g,x in q.items()if 1<=g-shift<=342};assert all(abs(x)<=bound.get(g,0)for g,x in tr.items());matrices[u,v]=tr;roots.append(dict(row=u,column=v,Tbar_row=i,Tbar_column=j,removed_t_power=shift,integer_polynomial=sorted(tr.items()),constant_and_negative_powers_ZERO=True,nine_reverse_signed_cofactor_terms=termhashes));guard(len(q)*2)
 dump(out/'all16-direct-four-border-identities.json',dict(common_denominator=str(TD**4),all16_roots=roots,all9_cofactor_integer_polynomials=[dict(row=i,column=j,integer_polynomial=sorted(cofs[i,j].items()))for i in range(3)for j in range(3)],
  pivot3_integer_polynomial=sorted(detp.items()),pivot_leading_QQ=plan['pivot_determinant_leading'],pivot_order71_verified=True,all_cofactor_reverse_signs='row_i,P[k] times cofactor[l,k] times Tbar[P[l],j], subtracted',all16_negative_and_constant_coefficients_ZERO=True))
 DP={0:{0:1}};stages=[]
 for i in range(4):
  nxt={}
  for mask,pp in DP.items():
   for j in range(4):
    if mask>>j&1:continue
    sign=(-1)**((mask>>(j+1)).bit_count());q=mul(pp,matrices[i,j],345-(3-i),'DP_products');dst=nxt.setdefault(mask|1<<j,{})
    for g,x in q.items():dst[g]=dst.get(g,0)+sign*x
   guard(sum(map(len,nxt.values()))+100)
  DP={m:{k:v for k,v in p.items()if v}for m,p in nxt.items()};stages.append(dict(row_completed=i+1,all_states=[dict(mask=m,coefficients=[[g,str(v)]for g,v in sorted(p.items())])for m,p in DP.items()]))
 target=DP.get(15,{}).get(345,0);assert abs(target)<=int(gate['target_integer_bound']);factor=Q(plan['original_detM_target_from_integer_det4_target_factor']);actual=target*factor
 dump(out/'all-fixed-univariate-DP-identities.json',dict(target_degree=345,all4_stages=stages,integer_target_coefficient=str(target),all_signs='(-1)^count(previousselectedcolumns>j)',common_denominator=str(TD**16),
  exact_original_detM_a5e91d10_coefficient=str(actual),exact_coefficient_normalization_factor=str(factor),original_candidate_det11='a²*detM'))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert new+external<33554432 and new<gate['new_RAW_upper']
 report=dict(id=ID,all_passed=True,frozen=True,ONE_finite_direct4border_job=True,all_source23932_and_paths_bound=True,all98_complete_old_prefix_polynomials_match=True,all16_finite_border_polynomials_verified=True,
  exact_target_degree=345,exact_target_whole_order=27338,original_detM_monomial=[5,91,10],integer_target_coefficient=str(target),original_QQ_target_coefficient=str(actual),target_NONZERO=bool(target),
  source_operations=counters,whole_owned_scalar_peak=peak,no_full_det11_ZERO_or_global57_claim=True,no_new_weight_window_point_prime_or_source572=True,
  total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256','integer_target_coefficient','original_QQ_target_coefficient']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

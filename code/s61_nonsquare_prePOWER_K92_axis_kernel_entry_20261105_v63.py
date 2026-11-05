#!/usr/bin/env python3
"""ONE complete110 finitepolynomial kernel propagation checks, no ranksearch."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_KERNEL_ENTRY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 base='reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_SOURCE_V2_20261105_v63';r=src(base+'/report.json');assert r['all_passed']and pins[base+'/report.json']=='00d66bf77479cd46700142f0d6f1c2fc77bf11dc04659347bd9078ae0c997f51';axis=src(base+'/whole-axis-literal-source.json');plan=src(base+'/fixed-kernel-propagation-source.json');g=src(base+'/resource-gate.json');receipt=src('reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_KERNEL_WORKER_GATE_20261105_v63/report.json');assert receipt['all_passed']
 A={(e['row'],e['column']):{int(k):int(v)for k,v in e['axis_integer_polynomial']}for e in axis['axis_coefficients']};assert len(A)==121 and sum(map(len,A.values()))==579;FC=[[Q(x)for x in row]for row in plan['right_factorization_FC']];R=[[Q(x)for x in row]for row in plan['left_bottom5_R']];RD=int(plan['right_clear']);LD=int(plan['left_clear']);assert all((RD*x).denominator==1 for row in FC for x in row)and all((LD*x).denominator==1 for row in R for x in row);pc=axis['fixed6minor_columns'];free=[j for j in range(11)if j not in pc];results=[];seen=set();adds=0;peak=0
 for job in plan['all110_finite_polynomial_linear_checks']:
  side=job['side'];i=job['row'];j=job['column'];key=(side,i,j);assert key not in seen;seen.add(key);expected=[]
  if side=='right_kernel':expected=[(i,free[j],RD)]+[(i,c,-int(RD*FC[z][free[j]]))for z,c in enumerate(pc)if FC[z][free[j]]]
  else:expected=[(z,j,int(LD*R[i][z]))for z in range(11)if R[i][z]]
  actual=[(q['row'],q['column'],int(q['integer_scale']))for q in job['source_terms']];assert expected==actual;pp={};N=int(job['integer_abs_upper'])
  for u,v,s in actual:
   for k,x in A[u,v].items():
    z=pp.get(k,0)+s*x;assert abs(z)<=N
    if z:pp[k]=z
    else:pp.pop(k,None)
    adds+=1
  assert len(pp)<=job['support_upper']and pp.get(0,0)==0
  results.append(dict(side=side,row=i,column=j,exact_integer_polynomial=sorted(pp.items()),known_denominator=str(RD if side=='right_kernel'else LD),is_ZERO=not pp,all_literal_terms_checked=actual));peak=max(peak,2*579+sum(len(e['exact_integer_polynomial'])for e in results)+121+2*len(pp)+15000);assert peak<=g['whole_owned_slots_upper']
 assert seen=={('right_kernel',i,j)for i in range(11)for j in range(5)}|{('left_kernel',i,j)for i in range(5)for j in range(11)};assert adds==g['scalar_additions_upper']
 right=all(e['is_ZERO']for e in results if e['side']=='right_kernel');left=all(e['is_ZERO']for e in results if e['side']=='left_kernel');rank6=right or left
 dump(out/'all110-polynomial-identities.json',dict(all110=results,all_constant_terms_ZERO=True,right55_ALL_ZERO=right,left55_ALL_ZERO=left,no_unknown_polynomial_inverse_or_ranksearch=True))
 dump(out/'exact-axis-rank-conclusion.json',dict(right_kernel_globally_extends=right,left_kernel_globally_extends=left,wholeaxis_rank_over_QQ_a_exact6=rank6,
  proof_if_rank6='complete5independentkernel vectors implyrank<=6; sourcefixed6minor hasnonzeroconstant so isnonzeropolynomial andrank>=6 overQQ(a)',
  NO_all_a_pointwise_rank6_claim=True,speciala_complement_retained='zerosoffixed6minor, plusnamedoriginalunitimages; d0isnotFULLU',
  both_constant_kernels_extend=right and left,next_source_if_both='Raxis*Caxis*Q0=diag(A6(a),0_5),A6(0)=I; extract(e,d) degree1 lower5block fromliteralMsource beforeanynewjob',
  if_obstruction='nonzerokernelpropagation onlyrejectsconstantkernel; doesNOTproveaxisrank>6 orfullminorNONZERO',no_full_detM_or_global57_claim=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert new<g['new_RAW_upper']and new+external<33554432
 report=dict(id=ID,all_passed=True,frozen=True,ONE_complete110_polynomial_linear_job=True,all579_axis_atoms_and_constkernel_terms_bound=True,all110_constant_terms_ZERO=True,right55_ALL_ZERO=right,left55_ALL_ZERO=left,
  axis_rank_over_QQ_a_exact6=rank6,independent_peer_pending=True,no_new_ranksearch_or_inverse=True,no_all_a_pointwise_or_fullU_cubic57_claim=True,scalar_additions=adds,owned_scalar_peak=peak,
  total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

#!/usr/bin/env python3
"""SOURCE ONLY: structural rank<=3 at Z1 and exact whole-Newton gap skip."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_RESIDUAL5_SUPPORT_ZERO_20261105_v63';W=(1,2244,68)
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 base='reports/S61_NONSQUARE_PREPOWER_K92_RESIDUAL5_FINITE_SOURCE_20261105_v63';r=src(base+'/report.json');assert r['all_passed']and pins[base+'/report.json']=='3d2cd38f949f5bc29507a542409638c17031e90bddac59e4e3b7280af96d3dbb';dag=src(base+'/finite-residual5-DAG.json');nodes=dag['nodes']
 for n in nodes:
  args=[nodes[i]for i in n['args']];assert all(i<n['id']for i in n['args'])
  if n['op']=='literal':v=Q(n['value']);D=v.denominator;N=abs(v.numerator)
  elif n['op']=='sum':D=math.lcm(*(int(a['denominator_upper'])for a in args));N=sum(abs(w)*int(a['integer_abs_upper'])*(D//int(a['denominator_upper']))for w,a in zip(n['weights'],args))
  elif n['op']=='product':D=math.prod(int(a['denominator_upper'])for a in args);N=math.prod(int(a['integer_abs_upper'])for a in args)
  elif n['op']=='power':D=int(args[0]['denominator_upper'])**n['exponent'];N=int(args[0]['integer_abs_upper'])**n['exponent']
  else:raise AssertionError(n['op'])
  assert D==int(n['denominator_upper'])and N==int(n['integer_abs_upper'])
 possible=[(e['row'],e['column'])for e in dag['residual5_roots']if int(nodes[e['node']]['integer_abs_upper'])];assert len(possible)==11 and all(i in[1,3]or j==1 for i,j in possible)
 constant_before=dag['known_total_order_beforeZ'];assert constant_before==26653;formal_bound=constant_before+5+(5-3);assert formal_bound==26660
 wb=src('reports/S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63/weighted-output.json');box=wb['output_metadata']['box'];weights=wb['weights'];linear=wb['output_metadata']['linear'];tropical=232227;limit=tropical-formal_bound;best=None;bestmon=None;checked=0
 for a in range(box[0][0],box[0][1]+1):
  for e in range(box[1][0],box[1][1]+1):
   lo,hi=box[2]
   for w,(lower,upper)in zip(weights,linear):
    base=w[0]*a+w[1]*e
    if w[2]:lo=max(lo,(lower-base+w[2]-1)//w[2]);hi=min(hi,(upper-base)//w[2])
    elif not lower<=base<=upper:hi=lo-1;break
   hi=min(hi,(limit-W[0]*a-W[1]*e)//W[2]);checked+=1
   if lo<=hi:
    value=W[0]*a+W[1]*e+W[2]*hi
    if best is None or value>best:best=value;bestmon=(a,e,hi)
 assert best==205522 and bestmon==(26,91,19);newbound=tropical-best;assert newbound==26705
 dump(out/'support-zero-proof.json',dict(possible_Z1_nonzero_cells=possible,all_other14_cells_zero_by_verified_DAG_bounds=True,rank_upper3_cover=dict(rows=[1,3],columns=[1]),
  matrix_decomposition='Z1=e1*row1 + e3*row3 + column1_with_rows1and3_removed * e1^T',
  rank_Z1_atmost3=True,Z_typed_constant_zero=True,formal_determinant_bound='Z=t*(Z1+tZ2+...); rankZ1<=3 =>orddetZ>=5+(5-3)=7',outer_order=26653,whole_order_atleast26660=True,
  no_Z1_values_or_determinant_computed=True,NO_full_detzero_claim=True))
 dump(out/'whole-Newton-gap-proof.json',dict(fixed_weight=W,original_tropical_degree=tropical,whole_polynomial_support_box=box,whole_polynomial_support_weights=weights,whole_polynomial_support_intervals=linear,
  previous_formal_order_lower_bound=formal_bound,maximum_weight_allowed_by_formal_order=limit,highest_remaining_possible_Newton_weight=best,highest_possible_monomial=bestmon,
  new_total_order_lower_bound=newbound,exponent_box_a_e_pairs_checked=checked,all_d_intervals_checked_without_materializinglattice=True,NO_parameter_values_or_numeric_points=True,
  proof='all determinant weights greaterthanlimit alreadyzero; the known exactNewton support has NO latticepoint between205523 and205567, so next possible weight is205522',
  next_possible_coefficient_NOT_evaluated=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);total=sum(p.stat().st_size for p in out.iterdir())
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,all_DAG_absolute_bounds_independently_rechecked=True,rank_firstpositive5_matrix_atmost3=True,formal_order_lower_bound=26660,Newton_sharpened_order_lower_bound=newbound,
  next_possible_monomial=bestmon,NO_new_Z_coefficients_or_determinant=True,no_full_det11_zero_nonzero_or_global57_claim=True,total_RAW_before_report=total,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

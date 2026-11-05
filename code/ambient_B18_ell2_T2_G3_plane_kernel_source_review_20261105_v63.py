#!/usr/bin/env python3
"""Independent source/initial-module review of the plane-kernel gate."""
import hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_T2_G3_PLANE_KERNEL_INPUT_GATE_20261105_v63';STEM=A+'_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sp(a):return {int(i):Q(c) for i,c in a}
def enc(d):return [[i,str(c)] for i,c in sorted(d.items())]
def add(a,b,c=Q(1)):
 d=dict(a)
 for i,v in b.items():
  t=d.get(i,Q(0))+c*v
  if t:d[i]=t
  else:d.pop(i,None)
 return d
def main():
 t=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 pp=ROOT/'reports'/A/'report.json';pr=json.loads(pp.read_text());assert pr['all_passed']
 inputs={'primary':{'path':str(pp.relative_to(ROOT)),'sha256':sha(pp)}};D={}
 for p,h in pr['input_sha256'].items():
  path=ROOT/p;assert sha(path)==h;inputs[p]={'path':p,'sha256':h}
  if path.name!='actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json':D[path.name]=json.loads(path.read_text())
 ep=pp.parent/'actual378-annihilation-checks-and105by238-condition.json';assert sha(ep)==pr['evidence_sha256'][ep.name]
 E=json.loads(ep.read_text());inputs['evidence']={'path':str(ep.relative_to(ROOT)),'sha256':sha(ep)}
 assert sha(ROOT/'code/ambient_B18_ell2_T2_G3_plane_kernel_source_20261105_v63.py')==pr['code_sha256']
 assert sha(ROOT/'research/notes'/f'{A}.md')==pr['note_sha256']
 lin=D['actual-complete-L48-and-degree2-rewrite-source-input.json'];nodes=lin['degree2_monomial_rewrite_nodes'];memo={}
 def nf(k):
  if k in memo:return memo[k]
  n=nodes[k]
  if n[0]=='standard':v={n[1]:Q(1)}
  elif n[0]=='zero_ell2':v={}
  else:
   v={}
   for j,c in n[2]:v=add(v,nf(j),Q(c))
  memo[k]=v;return v
 monidx={tuple(m):i for i,m in enumerate(lin['degree2_monomial_order44'])};checks=[]
 for target in range(6):
  for z in range(7):
   for x in range(9):
    m=[0]*9;m[z]+=1;m[x]+=1
    root=lin['all264_NF_root_nodes'][target][monidx[tuple(m)]]
    assert not nf(root);checks.append([target,z,x,root])
 assert len(checks)==378 and checks==E['all378_z_times_degree1_NF_zero_checks']
 prefix=next(p for p in pr['input_sha256'] if p.endswith('actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json'))
 F=fields(ROOT/prefix,['actual_active_G3linear42_by238_columns','PiY3_from102_columns'])
 active=list(map(sp,F['actual_active_G3linear42_by238_columns']));PiY=list(map(sp,F['PiY3_from102_columns']))
 assert len(active)==238 and len(PiY)==102 and all(not c for c in PiY[42:])
 conditions=[];projected=[]
 for col in active:
  assert all(0<=i//9<42 for i in col)
  fullplane={2*(i//9)+(i%9)-7:c for i,c in col.items() if i%9 in (7,8)};p={}
  for i,c in col.items():
   row,var=divmod(i,9)
   p=add(p,{9*a+var:c*q for a,q in PiY[row].items()})
  cond=add(fullplane,{84+7*(i//9)+(i%9):c for i,c in p.items() if i%9<7})
  conditions.append(cond);projected.append(p)
  # The projected plane rows are a literal constant linear combination of
  # the first84 full-plane rows; no rank search is needed for this implication.
  for a in range(3):
   for var in (7,8):
    assert p.get(9*a+var,0)==sum((PiY[row].get(a,Q(0))*fullplane.get(2*row+var-7,Q(0)) for row in range(42)),Q(0))
 assert list(map(enc,conditions))==E['condition105_by238_columns']
 assert list(map(enc,projected))==E['full_projected_active_G3linear_columns']
 B={'condition_input':sum(map(len,conditions)),'known_rank62_basis105_by62':105*62,
    'basis_reps_plus_kernel_partition238_by238':238*238,'pivot_metadata':4*238,
    'one_working_target_source_product_and_serialization':4*(105+238),'source_metadata':10000}
 assert B==E['next_rank_bounds']
 assert sum(B.values())==75684==E['next_rank_whole_QQ_bound']==pr['next_rank_whole_QQ_bound']<150000
 raw=(105*62+238*238)*(2*len(str(2**384))+24)+2_000_000
 assert raw==18167424==E['next_rank_raw_bound']==pr['next_rank_raw_bound']<32*1024*1024
 assert E['dynamic_commit_bits_stop']==384 and E['rank_not_computed'] and not E['new_quadratic_or_physical_coefficients_evaluated']
 proof={
 'sufficient_condition':'For each admitted constant c, full g1=G3_1 c has only z-linear coefficients and PiY0*g1=0. Hence its minimal6 column has zero first term.',
 'unknown_quotient_derivative':'Write the exact final target map on pureG3 columns as F(R):R102->R6. The first27 cancellation has bottom0 and is irrelevant. Later99 elimination gives F0=(PiY0;0), while F1 may contain both unknown quotient and bottom-unit derivatives. The quadratic column is F0*g2+F1*g1. Each term of F1*g1 has factor z_i*x_j and all378 target normal forms vanish.',
 'constant_graph_frame':'The accepted graph normalization is coherent and constant. Its action on a bottomzero column is identity; any subsequent formal unit derivative also multiplies its zero projected linear term. No unknown derivative is set zero.',
 'rank_bound':'The condition matrix is a fixed linear image of the accepted rank62 active full matrix, hence rank<=62. This does not choose a new basis or assert a rank.',
 'future_quadratic_source_gate':'For general q=D31I3c, retain its z projection. K0^-1 gives A126(q)=A126(q-Efree42*Pi42*q), since K0^-1*Efree42 has zero selected126 component. Only the corrected vector is asserted to be in the fixed182 span. No full-q span assumption.',
 'scope':'This is only an active sufficient family. It is neither the unrestricted projected393 kernel nor the full1012-dimensional initial kernel. Actual HF2 or total length is not determined.'}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 data={'all378_NF_annihilation_checks_independently_rebuilt':True,'all238_105row_columns_and_projected_plane_implications':True,
       'condition_shape':[105,238],'rank_upper_bound':62,'rank_or_kernel_not_computed':True,'source_proof':proof,
       'next_numeric_slots':75684,'next_raw_bytes':18167424,'dynamic384_guard_is_fail_closed_not_height_prediction':True}
 path=out/'independent-plane-kernel-source-and-bound-proof.json';path.write_text(json.dumps(data,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'sufficient_plane_kernel_and_unknown_derivative_NF_cancellation_accepted':True,
         'rank_or_quadratic_values_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),
         'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{path.name:sha(path)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
 print(json.dumps({'all_passed':True,'source_shape':[105,238],'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()

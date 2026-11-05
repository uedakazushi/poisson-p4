#!/usr/bin/env python3
"""Source-only enlargement using z*m zero in the known linear cover."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_G3_PLANE_KERNEL_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 for p in ['reports/AMBIENT_B18_ELL2_T2_FULL_G3_LINEAR_KERNEL_20261105_v63_REVIEW/report.json',
  'reports/AMBIENT_B18_ELL2_T2_G3_QUADRATIC_INPUT_GATE_20261105_v63_REVIEW/report.json',
  'reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_GB_COVER_20261105_v63_REVIEW/report.json']:
  assert read(p)['all_passed']
 lin=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json')
 nf=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_DEGREE2_20261105_v63/actual432-overlap-rank-and-known-linear-cover-degree2.json')
 memo={i:H.parse(v) for i,v in nf['all264_linear_NF_columns']};mi={tuple(m):i for i,m in enumerate(lin['degree2_monomial_order44'])};checks=[]
 for target in range(6):
  for z in range(7):
   for var in range(9):
    mon=[0]*9;mon[z]+=1;mon[var]+=1;idx=mi[tuple(mon)];node=lin['all264_NF_root_nodes'][target][idx]
    assert not memo[node];checks.append([target,z,var,node])
 assert len(checks)==378
 source='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json';hs[source]=H.sha(ROOT/source)
 active=list(map(H.parse,H.selected(source,'actual_active_G3linear42_by238_columns')));PiY=list(map(H.parse,H.selected(source,'PiY3_from102_columns')))
 assert all(not q for q in PiY[42:])
 condition=[];projected=[]
 for col in active:
  v={};p={}
  for i,c in col.items():
   row,var=divmod(i,9)
   if var in (7,8):v[row*2+var-7]=c
   for t,a in PiY[row].items():p=add(p,{t*9+var:a*c})
  for i,c in p.items():
   target,var=divmod(i,9)
   if var<7:v[84+target*7+var]=c
   else:assert not c or any(col.get(row*9+var) for row in range(42))
  condition.append(H.guard(v));projected.append(p)
 # The84 full-plane rows imply the omitted6 projected-plane rows.
 assert len(condition)==238
 bounds={'condition_input':sum(map(len,condition)),'known_rank62_basis105_by62':105*62,
  'basis_reps_plus_kernel_partition238_by238':238*238,'pivot_metadata':4*238,
  'one_working_target_source_product_and_serialization':4*(105+238),'source_metadata':10000}
 total=sum(bounds.values());assert total<150000
 raw=(105*62+238*238)*(2*len(str(2**384))+24)+2_000_000;assert raw<32*1024*1024
 proof={'condition':'full active G3linear plane(r,ell)=0 and PiY0G3linear z-part=0; together imply complete projectedG3linear=0',
  'not_fullzero_condition':'nonzero full z-linear terms are retained and may belong the99-unit span',
  'unit_projection_derivative':'G3firstc has onlyz terms; every unknown homogeneouslinear quotient derivative times it is a sum e_target*z_i*x_j, whose378 fixedNF checks vanish. Same for the known bottomunit derivative. This proves only its projected quadratic term is independent of the unknown derivatives.',
  'current_G3_quadratic_callback':'only PiY0 G3_2c remains; the prior finiteK(t) formula is valid. Allmixedz quadratic monomials reduce0, so onlyr²/r*ell innerchannels must be bound before values.',
  'normalized_frame':'existing coherent constant triangular target normalization; G3bottom0. No canonical unknownM1 constant is set0.',
  'nonactive_not_queried':'This gate is only an active238 source family; no automatic claim for remaining180 columns or completeHF2.',
  'rankupper62':'condition is a fixedlinear projection of the fullactiveG3linear rank62 matrix; no newrank search in sourcegate',
  'negative_or_nonzero_scope':'A later nonzero image supplies true quadratic initialrelations and an uppercover only; allzero stops onlythisfamily.'}
 ep=out/'actual378-annihilation-checks-and105by238-condition.json';H.dump(ep,{'all378_z_times_degree1_NF_zero_checks':checks,
  'condition105_by238_columns':[srec(v) for v in condition],'full_projected_active_G3linear_columns':[srec(v) for v in projected],
  'condition_row_order':'first84=42targets*(r,ell); final21=3projectedtargets*zA0..6',
  'rank_not_computed':True,'new_quadratic_or_physical_coefficients_evaluated':False,'source_proof':proof,
  'next_rank_bounds':bounds,'next_rank_whole_QQ_bound':total,'next_rank_raw_bound':raw,'dynamic_commit_bits_stop':384})
 report={'id':STEM,'all_passed':True,'all378_annihilation_checks':True,'rank_not_computed':True,'new_coefficients_evaluated':False,
  'shape':[105,238],'next_rank_whole_QQ_bound':total,'next_rank_raw_bound':raw,'input_sha256':hs,
  'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','shape','next_rank_whole_QQ_bound','next_rank_raw_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())

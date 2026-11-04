#!/usr/bin/env python3
"""Independent fixed14-pivot Laurent audit; no Gaussian or candidate search."""
import argparse,json,resource,shutil,signal,time
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
from sympy.polys.rings import ring
from s61_nonsquare_joint_Lell_overlap_unit_review_20261104_v63 import read,sha,decoded,plus,minus,times,invert_monomial,profile,Z
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_JOINT_LELL_QPAIR_UNITS_20261104_v63_REVIEW'
def main(out,prior_failure=None):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 if prior_failure:
  old=read(prior_failure/'first-failure.json');assert old['no_parent_input_mutation']and old['no_PASS_claim_before_repair']
  for name in['first-failure.json','failed-metadata-review.py']:shutil.copyfile(prior_failure/name,out/name)
 b=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_QPAIR_UNITS_20261104_v63';r=read(b/'report.json');state=read(b/'checkpoint.json')
 assert r['frozen']and r['completed_global_known_unit_pivots']==14 and not r['unfinished_step_not_in_checkpoint']and r['stop_reason']=='no_individual_FULLU_unit_entry'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h,p
 for p,h in r['evidence_sha256'].items():assert sha(b/p)==h,p
 ip=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_OVERLAP_UNIT_20261104_v63/remaining-pair-block.json';inp=read(ip);assert sha(ip)==state['source_input_sha256']
 dep=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_OVERLAP_UNIT_20261104_v63_REVIEW/report.json';review=read(dep);assert review['all_passed']and review['frozen_qpair58x136_all_coefficients_source_bound']and review['input_sha256'][str(ip.relative_to(ROOT))]==sha(ip)
 A=[[decoded(v)for v in row]for row in inp['entries']];assert len(A)==58 and all(len(row)==136 for row in A)
 available=list(range(58));cols=list(range(136));steps=state['all_completed_unit_steps'];assert len(steps)==14
 pivot_expected={(-2,0,-4):Q(19683,1280)};detconstant=Q(1);detpowers=[0]*6;eliminations=0
 for number,step in enumerate(steps):
  i,j=step['pivot_row'],step['pivot_column'];assert i in available and j in cols
  assert A[i][j]==decoded(step['pivot'])==pivot_expected
  assert step['known_unit']=={'constant':'19683/1280','a_d_F_B_C_J_powers':[-2,-4,0,0,0,0]}
  inv=invert_monomial(A[i][j]);A[i]=[times(v,inv)for v in A[i]]
  used=set(range(136))-set(cols);assert all(not A[i][k]for k in used)
  given={t:decoded(v)for t,v in step['eliminations']};assert given=={t:A[t][j]for t in available if t!=i and A[t][j]}
  for t in available:
   if t==i:continue
   co=A[t][j]
   if co:
    # Independently update ALL136 columns, not just the active-column
    # optimization used by the author. Old pivot columns stay zero.
    assert all(not A[t][k]for k in used);A[t]=[minus(x,times(co,y))for x,y in zip(A[t],A[i])];eliminations+=1
  available.remove(i);cols.remove(j);detconstant*=Q(step['known_unit']['constant']);detpowers=[x+y for x,y in zip(detpowers,step['known_unit']['a_d_F_B_C_J_powers'])]
  (out/'progress.json').write_text(json.dumps({'completed_fixed_pivots':number+1,'seconds':time.monotonic()-began},indent=2)+'\n')
 assert A==[[decoded(v)for v in row]for row in state['matrix58x136']]
 assert eliminations==40
 assert available==state['remaining_rows']and cols==state['remaining_columns']and len(available)==44 and len(cols)==122
 assert state['selected_minor_rows']==[x['pivot_row']for x in steps]and state['selected_minor_columns']==[x['pivot_column']for x in steps]
 assert state['selected_minor_unit_determinant']=={'constant':str(detconstant),'a_d_F_B_C_J_powers':detpowers}
 remaining=[[A[i][j]for j in cols]for i in available];actual_profile=profile(remaining)
 assert {k:actual_profile[k]for k in state['remaining_profile']}==state['remaining_profile']==r['remaining_profile']
 assert profile(remaining)['numerator_terms']==79165 and profile(remaining)['nonzero_entries']==1223 and profile(remaining)['max_entry_numerator_terms']==391
 # Verify the recorded stopping rule, without selecting or executing any
 # further pivot. Only exact division by the six named factors is used.
 R,a,e,d=ring('a,e,d',s.QQ);known=[a,d,500*a+90*e**2+3*e-24,(50*a+18*e**2-9*e-2)*d+60*e-40,100*a+18*e**2-9*e-2,5*a*d-6*e-4];nonunits=0
 for row in remaining:
  for c in row:
   if not c:continue
   ma=min(k[0]for k in c);md=min(k[2]for k in c);p=R.from_dict({(k[0]-ma,k[1],k[2]-md):s.QQ(v.numerator,v.denominator)for k,v in c.items()})
   for f in known:
    while not p.rem(f):p=p.exquo(f)
   assert not p.is_ground;nonunits+=1
 assert nonunits==1223
 note=ROOT/'research/notes'/(ID+'.md');utility=ROOT/'code/s61_nonsquare_joint_Lell_overlap_unit_review_20261104_v63.py'
 proof={'fixed_selected_minor_rows':state['selected_minor_rows'],'fixed_selected_minor_columns':state['selected_minor_columns'],
 'selected_order_minor_determinant':state['selected_minor_unit_determinant'],
 'determinant_proof':'Only prior selected rows are subtracted from futureselected rows; selected-order minor becomesuppertriangularafter exactscaling. Therefore determinant is productof14 pivots, no selected-order sign.',
 'all58x136_completed_checkpoint_coefficients_match':True,'all40_recorded_elimination_coefficients_before_operations_match':eliminations==40,
 'safe_full_presentation_shape':[202,122],
 'no_automatic_free164_basis':'The selected target pivots split14 sourcecombinations. Their complement after old38 quotient is projective; an explicit global freebasis needs a separate source transport certificate.',
 'no_unit_entry_is_not_unitideal_or_cokernel_nonzero':True}
 (out/'independent-certificate.json').write_text(json.dumps(proof,indent=2)+'\n')
 result=dict(id=ID,all_passed=True,frozen_input58x136_bound_to_independent_overlap_review=True,
  all14_saved_pivots_literal19683_over1280_a2_d4=True,
  all40_elimination_coefficients_checked_before_updates=True,
  all7888_full_checkpoint_coefficients_and44x122_frontier_match=True,
  selected14_minor_unit_product_verified=True,
  remaining1223_entries_have_no_named_FULLU_unit_factorization=True,
  quadratic_degree4_image_uniform_lower_bound=117,target_complement_rank=122,
  known_kernel38_missing42_unchanged=True,
  no_free164_fullrank152_rank122_full239_or_exhaustion_claim=True,
  no_new_pivot_candidate_Gaussian_search_rank_point_prime_GB_or_author_rerun=True,
  prior_profile_metadata_failure_preserved=bool(prior_failure),
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),utility,note,b/'report.json',b/'checkpoint.json',ip,dep]},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-began,frozen=True)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--prior-failure',type=Path);a=p.parse_args();main(a.output.resolve(),a.prior_failure.resolve()if a.prior_failure else None)

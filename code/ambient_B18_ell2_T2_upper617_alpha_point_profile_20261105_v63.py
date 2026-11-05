#!/usr/bin/env python3
"""Point99 transpose graph units, metadata only; no dual solve or rank."""
import argparse,json,resource,signal,time
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_POINT_PROFILE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True)
 mu='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json';raw=H.selected(mu,'actual_selected_Y0_columns102');cols=list(map(H.parse,raw));free={1,36,38};target=sorted(set(range(102))-free);assert len(cols)==len(target)==99
 # A=Ystar^T restricted to the99 nonfree target axes. Every entry is a
 # literal previously certified point coefficient. Graph elimination leaves
 # all surviving entries unchanged; no rational response is calculated.
 A=[{i:x for i,x in c.items() if i not in free} for c in cols];R=set(range(99));C=set(target);steps=[]
 while True:
  found=None
  for r in sorted(R):
   supp=[c for c in A[r] if c in C]
   if len(supp)==1:found=('row_singleton',r,supp[0]);break
  if found is None:
   for c in sorted(C):
    supp=[r for r in R if c in A[r]]
    if len(supp)==1:found=('column_singleton',supp[0],c);break
  if found is None:break
  kind,r,c=found;assert A[r][c]
  if kind=='row_singleton':
   recipe={'solve':'s_c=rhs_r/literalpivot','update_rhs_other_equations':'rhs_i -= literalA_i,c *s_c','point_compatible_no_fill':True}
  else:recipe={'deferred_backsolve':'s_c=(rhs_r-sum_remaining literalA_r,j*s_j)/literalpivot','other_rhs_unchanged':True,'point_compatible_no_fill':True}
  steps.append({'kind':kind,'equation_selected99':r,'unknown_target102':c,'literal_pivot':str(A[r][c]),'recipe':recipe});R.remove(r);C.remove(c)
 assert len(R)==len(C)
 ev={'source_selected99_Y0_columns102':raw,'point_free_axes':[1,36,38],'transpose_equation':'s*Ystar0=K => A*s=K withA=Ystar0^T restrictedtononfreeaxes','all_graph_singleton_unit_steps':steps,'remaining_equations_selected99':sorted(R),'remaining_unknown_target102':sorted(C),'residual_matrix_literal_rows':[[r,[[c,str(A[r][c])] for c in sorted(C) if c in A[r]]] for r in sorted(R)],'residual_invertibility_source':'actualpointY0rank99withEfreecomplementpreviouslypeered;nonzerographunitsimplyremainingmatrixinvertible.Noresidualrankordeterminantqueried.','no_signed_K_rhs_dual_alpha_or_new_inverse_coefficient_computed':True}
 ep=out/'point99-transpose-no-fill-unit-profile-and-reduced-source.json';H.dump(ep,ev);owned=count([raw,cols,A,R,C,steps,ev])+10000;assert owned<60000
 report={'id':STEM,'all_passed':True,'metadata_only_no_pointsolve_or_rank':True,'singleton_units':len(steps),'remaining_shape':[len(R),len(C)],'residual_nnz':sum(sum(c in A[r] for c in C) for r in R),'source_owned_bound':owned,'input_sha256':{mu:H.sha(ROOT/mu)},'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','singleton_units','remaining_shape','residual_nnz','source_owned_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

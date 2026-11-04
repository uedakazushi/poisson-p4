#!/usr/bin/env python3
"""Extract the nonzero10minor from the already-computed frozen QQ RREF columns."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_ELL_CUBIC_RANK10_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def determinant(a):
 a=[r[:]for r in a];d=Q(1)
 for j in range(len(a)):
  k=next((i for i in range(j,len(a))if a[i][j]),None);assert k is not None
  if k!=j:a[k],a[j]=a[j],a[k];d=-d
  p=a[j][j];d*=p
  for i in range(j+1,len(a)):
   t=a[i][j]/p
   for k in range(j+1,len(a)):a[i][k]-=t*a[j][k]
   a[i][j]=Q()
 return d
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 folder=ROOT/'reports/S61_NONSQUARE_ELL_COLON_CUBIC11_20261104_v63';rp=folder/'report.json';r=json.loads(rp.read_text());assert r['all_saved_source_and_exact_QQ_rank_checks_passed']and r['cubic11_rank_at_saved_point']==10
 for rel,h in r['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 for rel,h in r['evidence_sha256'].items():assert sha(folder/rel)==h,rel
 np=folder/'native.json';n=json.loads(np.read_text());S=[[Q(c)for c in row]for row in n['actual11x16']];cols=n['actual11_pivot_columns'];assert len(cols)==10
 # The columns are fixed by the previous RREF. Select its independent rows
 # by one exact echelon extraction, without trying any alternate minor/point.
 chosen=[];basis={}
 for i,row in enumerate(S):
  v=[row[j]for j in cols]
  for j in sorted(basis):
   c=v[j]
   if c:v=[a-c*b for a,b in zip(v,basis[j])]
  pivot=next((j for j,c in enumerate(v)if c),None)
  if pivot is not None:
   c=v[pivot];basis[pivot]=[a/c for a in v];chosen.append(i)
  if len(chosen)==10:break
 assert len(chosen)==10;minor=[[S[i][j]for j in cols]for i in chosen];delta=determinant(minor);assert delta
 d24=Q(n['quad19_minor'])*Q(n['heavy5_minor']);d34=d24*delta;assert d34
 fullrows=n['quad_pivot_rows']+[20+i for i in n['restricted_heavy5_pivot_rows']]+[20+n['remaining11_genuine_rows'][i]for i in chosen]
 fullcols=n['quad_pivot_columns']+list(range(35,40))+[n['quotient16_columns'][j]for j in cols]
 assert len(fullrows)==len(set(fullrows))==34 and len(fullcols)==len(set(fullcols))==34
 evidence={'same_already_computed_QQ_point':[1,2,3],'fixed_columns_from_completed_RREF':cols,'extracted10_rows':chosen,'exact10x10_minor':[[str(c)for c in row]for row in minor],'exact_nonzero10_minor_determinant':str(delta),
 'parametric_nonzero34minor':{'source':'the exact restricted36x40 polynomial circuit in parent','rows':fullrows,'columns':fullcols,'oldpoint_value':str(d34),'Schur_identity':'det34=detQuad19*detHeavy5*detR11minor10'},
 'generic_rank_interval_for_actual_R11':[10,11],'rank11_not_proved':True,'rank10_generic_upper_not_proved':True,'extra_chart_factor':'Delta10 nonzero, not a known FULLU unit',
 'additional_exact_reduction':'On U78,heavy7,projection122->V48 onto and thisDelta10 chart: split15knownkernel,split10canonicalsource/targetimages,then48projection; actualsamecoker map E49->V6; genericremainingkernel43',
 'Hq_lift_clarification':'Use originalsource KZ=H*(ellZ)-ell*fZ and KW=H*(ellW)-ell*fW, with fZ=HZ-pZ,fW=HW-pW; images exactlyell*pZ,ell*pW. Their C122 coordinates include fixedcommutator adjustments after NF_H, not a claim that unadjusted NF_H(HellZ)=ell*pZ literally.',
 'source_r9_lift':'ell*r_j; their images areell*r_j before dividingell through the exactcolon injection',
 'no_new_point_rank_job_field_GB_or_quartic_coefficient_evaluation':True}
 (out/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),rp,np,folder/'circuit.json',folder/'oldpoint-native-input.txt']
 report={'id':ID,'exact_saved_RREF_minor_extraction_passed':True,'generic_actual_cubic11_rank_atleast10_proved':True,'generic_rank11_or_upper10_not_claimed':True,'nonempty_Delta10_chart_actual49x6_samecoker_recipe':True,'no_allpoint_or_full239_claim':True,'no_new_point_rank_GB_or_functionfield_solver':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-started,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);a=p.parse_args();main(a.output.resolve())

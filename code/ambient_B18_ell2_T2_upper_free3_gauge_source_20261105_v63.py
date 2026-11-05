#!/usr/bin/env python3
"""Fixed-source free3 comparison-frame existence; no M1 or619 values."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER_FREE3_GAUGE_SOURCE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 up=read('reports/AMBIENT_B18_ELL2_T2_UPPER2_CUBIC_OUTPUT_SOURCE_20261105_v63/fixed-upper2-one-cubic-output-identity-and-source-frontier.json');s=[H.parse(up[k]) for k in ['s0_on_F2_source102','s1_on_F2_source102']]
 lin=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json');L=list(map(H.parse,lin['all48_true_initial_linear_columns54']));rewrite=[(p,H.parse(v)) for p,v in lin['normalized_linear_rewrite_rows']]
 W=[F(-3375,11),F(-50),F(25,2)];R=[{38:F(1)},{1:F(3652,1095)}];maximum=0
 def mark(q):
  nonlocal maximum
  maximum=max(maximum,q[0].bit_length(),q[1].bit_length());return q
 def ex(x):x=F(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  d=math.lcm(a[0],b[0]);return mark((d,a[1]*(d//a[0])+b[1]*(d//b[0])))
 for row in s:
  for v in R:
   b=(1,0)
   for i,x in v.items():b=plus(b,mul(ex(x),ex(row.get(i,F(0)))))
 for col in L:
  for var in [7,8]:
   q=[col.get(9*(3+i)+var,F(0)) for i in range(3)]
   for i in range(3):plus(mul(ex(q[i]),ex(W[2])),mul(ex(q[2]),ex(W[i])))
 for target in range(3):
  for z in range(7):
   v={9*target+z:(1,1)}
   for p,row in rewrite:
    if p not in v:continue
    c=v[p]
    for i,x in row.items():v[i]=plus(v.get(i,(1,0)),mul(c,ex(x)))
    v.pop(p,None)
 assert maximum<384
 H.dump(out/'phase.json',{'phase':'positive_fixed_source_bound_before_signed_checks','positive_prior_bits':maximum,'no_M1_or_619_values':True})
 actual=0
 def ck(x):
  nonlocal actual
  actual=max(actual,abs(x.numerator).bit_length(),x.denominator.bit_length());assert actual<=maximum;return x
 assert [[sum((ck(x*row.get(i,F(0))) for i,x in v.items()),F(0)) for v in R] for row in s]==[[F(1),F(0)],[F(0),F(1)]]
 plane=[]
 for j,col in enumerate(L):
  for var in [7,8]:
   q=[col.get(9*(3+i)+var,F(0)) for i in range(3)]
   for i in range(3):assert ck(q[i]*W[2])==ck(q[2]*W[i])
   if any(q):plane.append([j,var,[str(x) for x in q]])
 for target in range(3):
  for z in range(7):
   v={9*target+z:F(1)}
   for p,row in rewrite:
    c=v.get(p,F(0))
    if c:
     for i,x in row.items():
      q=ck(v.get(i,F(0))-ck(c*x))
      if q:v[i]=q
      else:v.pop(i,None)
   assert not v
 premises={
 'coherent_gauge':'reports/AMBIENT_B18_ELL2_T2_FIRST_GRAPH_COMPARISON_GAUGE_PREFLIGHT_20261105_v63/report.json',
 'known_actual_graph':'reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_VALUES_20261105_v63_REVIEW/report.json',
 'e23_invariants':'reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_20261105_v63_REVIEW/report.json',
 'e29_e30_invariants':'reports/AMBIENT_B18_ELL2_T2_ALPHA29_BETA29_20261105_v63_REVIEW/report.json',
 'exact_upper_source':'reports/AMBIENT_B18_ELL2_T2_UPPER619_SOURCE_CLEARING_20261105_v63/report.json'}
 for p in premises.values():assert read(p)['all_passed']
 ep=out/'fixed-free3-axis-gauge-right-inverse-and-full-L48-preservation.json';H.dump(ep,{'right_inverse_R_columns102':[[[i,str(x)] for i,x in sorted(v.items())] for v in R],'all_s0_s1_right_inverse_identities':True,'all48_lower_r_ell_coefficients_proportional_to_Wv':True,'nonzero_lower_plane_source_columns':plane,'all21_upper_z_vectors_in_L48':True,
 'actual_normalizer_typed_not_evaluated':'A0=-R*[lambda1 MrH0;lambda1 MellH0]*Pi3 in the existing Wv-normalized frame; A0H0Wv=0. Transport M,Y,Yp and actual27/99units coherently.',
 'coherent_new_frame_M1H0_output_r2ell':['sum_i alpha_i*B2_i_rell+beta_i*B2_i_rr','0'],'canonical_complementary_M1_entries_not_computed':True,'b1_unchanged_reason':'Ybarc1 changesbyA0Xbarc1=0;starY0andfixedG3linearunchanged','no_remaining619_values_or_new_source_coefficients':True,'source_premises':premises})
 report={'id':STEM,'all_passed':True,'existence_frame_source_lemma_only':True,'actual_comparison_values_queried':False,'full_L48_preservation_fixed_sources_checked':True,'positive_prior_bits':maximum,'actual_fixed_source_bits':actual,'owned_numeric_bound':count([up,s,lin,L,rewrite,R,plane,premises])+12000,'input_sha256':hs,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start,'independent_source_review_pending':True};assert report['owned_numeric_bound']<60000;H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','positive_prior_bits','actual_fixed_source_bits','owned_numeric_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

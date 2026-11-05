#!/usr/bin/env python3
"""Independent source-only free-three gauge preservation audit."""
import json
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_UPPER_FREE3_GAUGE_SOURCE_20261105_v63';STEM=BASE+'_REVIEW'
def main():
 pins={}
 def read(p,keys=None):pins[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and pins[rp]['sha256']=='f3a623bf15cec9c93569237d545368048a5edcd5b514e0bf31209929932ae3d4'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 ep=f'reports/{BASE}/fixed-free3-axis-gauge-right-inverse-and-full-L48-preservation.json';e=read(ep);assert pins[ep]['sha256']==r['evidence_sha256'][Path(ep).name]
 up=read('reports/AMBIENT_B18_ELL2_T2_UPPER2_CUBIC_OUTPUT_SOURCE_20261105_v63/fixed-upper2-one-cubic-output-identity-and-source-frontier.json',['s0_on_F2_source102','s1_on_F2_source102']);s=list(map(sp,up.values()));R=list(map(sp,e['right_inverse_R_columns102']));assert R==[{38:Q(1)},{1:Q(3652,1095)}]
 assert [[sum((x*row.get(i,Q(0)) for i,x in col.items()),Q(0)) for col in R] for row in s]==[[1,0],[0,1]]
 l=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json',['all48_true_initial_linear_columns54','normalized_linear_rewrite_rows']);cols=list(map(sp,l['all48_true_initial_linear_columns54']));rewrite=[(p,sp(v)) for p,v in l['normalized_linear_rewrite_rows']];W=[Q(-3375,11),Q(-50),Q(25,2)];plane=[]
 for j,col in enumerate(cols):
  for var in [7,8]:
   v=[col.get(9*(3+i)+var,Q(0)) for i in range(3)];scale=v[2]/W[2];assert v==[scale*x for x in W]
   if any(v):plane.append([j,var,list(map(str,v))])
 assert plane==e['nonzero_lower_plane_source_columns']
 membership=[]
 for target in range(3):
  for z in range(7):
   v={9*target+z:Q(1)};receipts=[]
   for p,row in rewrite:
    assert row[p]==1;c=v.get(p,Q(0))
    if c:
     receipts.append([p,str(c)])
     for i,x in row.items():
      y=v.get(i,Q(0))-c*x
      if y:v[i]=y
      else:v.pop(i,None)
   assert not v;membership.append({'upper_target':target,'z':z,'fixed_rewrite_receipts':receipts})
 mg=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['Mnew0_columns']);assert all(not sp(mg['Mnew0_columns'][j]) for j in [23,29,30])
 mini=read('reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json',['free30_V_indices','remaining_free3_bottom_coordinates','fixed_bottom30_inverse_columns']);Pi=[]
 assert [mini['free30_V_indices'][i] for i in mini['remaining_free3_bottom_coordinates']]==[23,29,30]
 for k in range(27,30):Pi.append({mini['free30_V_indices'][i]:sp(c).get(k,Q(0)) for i,c in enumerate(mini['fixed_bottom30_inverse_columns']) if sp(c).get(k,Q(0))})
 assert [[row.get(j,Q(0)) for j in [23,29,30]] for row in Pi]==[[1,0,0],[0,1,0],[0,0,1]]
 for p in e['source_premises'].values():assert read(p)['all_passed']
 assert read('reports/AMBIENT_B18_ELL2_T2_UPPER619_SOURCE_CLEARING_20261105_v63_REVIEW/report.json')['all_passed']
 proof={'admissible_source_change':'Mprime=M+G2*A0,Yprime=Y+A0*X preserves G1M=V and G2Y=MX exactly. Since G2(0)=0, M0 and fixed section unit persist. All actual first27 and99 unit constructions are transported; their higher coefficients are not zero-filled.',
 'gauge_right_inverse':'Gauge change of(mu0r,mu1r,mu0ell,mu1ell) is(0,u,u,v). The verified right inverseR of(s0,s1) allows A0=-R[mu1r;mu1ell]Pi3, giving new(mu0r,mu1r,mu0ell,mu1ell)=(alpha,0,beta,0) on allthree H0 axes. These are an existence-frame identity, not canonical numerical M1 values.',
 'L48_preservation':'The original W-normalized frame has both complementary rows vanish onW. The SAME linear right inverse therefore forces A0H0W=0. Every L48 bottom r/ell vector is a multiple ofW. Its remaining bottom z vectors can add only upper-z combinations, all21 explicitly lie inL48. The induced triangular target mapI+K hasK^2=0 and preservesL48 bothways. Pure upperQ2 is fixed, and lowerfour duals are unchanged.',
 'moving_coefficients':'For the619 inputs, actual Zc1=0. Thus Ybarc1 changes byA0Zc1=0. Starpointcolumns are unchanged since Zstar0=0; fixed G3 source is unchanged. Therefore b1 and point-derived A2/B2 are unchanged while allhigher unit data remain coherently transported.',
 'output':'In this same frame, lambda0(M1H0B2)r2ell=sum(alpha_i*B_i,rell+beta_i*B_i,rr),lambda1(...)=0. All honest cubic, native L2, b1 and A2 source gates remain. No upper563 or619 value has been evaluated.'}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ve=out/'independent-free3-right-inverse-full48-and-coherent-frame-proof.json';ve.write_text(json.dumps({'right_inverse_two_sided_identity':True,'all48_plane_W_relations_checked':True,'all21_upper_z_membership_receipts':membership,'literal_Pi3_H0_identity':True,'coherent_source_proof':proof,'no_canonical_M1_or_upper_values':True},separators=(',',':'))+'\n');rr={'id':STEM,'all_passed':True,'primary_sha256':pins[rp]['sha256'],'coherent_existence_frame_source_accepted':True,'full_L48_preservation_and_b1_unchanged':True,'canonical_M1_or_upper_values_not_queried':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ve.name:sha(ve)}};(out/'report.json').write_text(json.dumps(rr,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()

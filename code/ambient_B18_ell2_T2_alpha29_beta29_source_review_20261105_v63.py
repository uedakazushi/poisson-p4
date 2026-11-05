#!/usr/bin/env python3
"""Fresh universal-two-dot source/custody/bound check, no dot evaluation."""
import json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_ALPHA29_BETA29_INPUT_GATE_20261105_v63';STEM=BASE+'_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p,keys=None):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and r['new_pairing_count_only']==2 and r['no_pairing_or_new_function_value_evaluated']
 assert r['code_sha256']==sha(ROOT/'code/ambient_B18_ell2_T2_alpha29_beta29_gate_20261105_v63.py')
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 ep=f'reports/{BASE}/two-pairing-source-binding-and-positive-height-input.json';d=read(ep);assert sha(ROOT/ep)==r['evidence_sha256'][Path(ep).name]
 graph=read('reports/AMBIENT_B18_ELL2_T2_LINEAR_GRAPH4_SOURCE_20261105_v63/actual23-linear-graph-four-common-contractions-and-six-adjoint-input.json')
 assert list(map(sp,graph['H0_free3_columns']))==[{23:Q(1)},{29:Q(1)},{30:Q(1)}]
 M=list(map(sp,read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['Mnew0_columns'])['Mnew0_columns']))
 assert all(not M[i] for i in [23,29,30])
 rho=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63/actual3960-rho-factor-identities-or-obstructions.json',['actual_two_three_component_bar_rho33'])['actual_two_three_component_bar_rho33']
 assert rho==d['two_universal_bar_rho33']
 for stem in ['AMBIENT_B18_ELL2_T2_GRAPH_MIXED_ADJOINT_20261105_v63_REVIEW','AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63_REVIEW']:
  assert read('reports/'+stem+'/report.json')['all_passed']
 phy=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json')['physical12_coefficients'][7]
 assert phy['physical_row']==phy['roots']['physical_row']==63 and phy['roots']==d['source_raw63_roots']
 assert phy['roots']['delta_operator']==['subtract',20144,40922] and phy['roots']['bB_node']==38648
 V29=[sp(phy['delta']),{}, {i:-v for i,v in sp(phy['bB']).items()}];V30=[{},sp(phy['delta']),sp(phy['bA'])]
 assert list(map(sp,d['primitive_V29_raw63_mask1']))==V29 and list(map(sp,d['primitive_V30_raw63_mask2_source_check_only']))==V30
 W=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_VALUES_20261105_v63/phase1-actual-two-invariants-from-raw33-primitive.json')
 E23=read('reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_20261105_v63/phase1-actual-alpha23-beta23-and-primitive23.json')
 VW=list(map(sp,W['primitive_V_Wv_full33']));V23=list(map(sp,E23['actual_primitive_V23_raw33']))
 assert d['known_alpha_beta_on_Wv']==W['actual_two_invariants_alpha_beta'] and d['known_alpha23_beta23']==E23['actual_alpha23_beta23']
 assert sp(d['source_Wv_coefficients'])==sp(graph['actual_H0v_source56'])=={23:Q(-3375,11),29:Q(-50),30:Q(25,2)}
 maximum=0
 def mark(v):
  nonlocal maximum
  maximum=max(maximum,v[0].bit_length(),v[1].bit_length());return v
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator)))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 dots=[];recover=[]
 for ch,row in enumerate(map(lambda v:list(map(sp,v)),rho)):
  b=(1,0)
  for target,rr in enumerate(row):
   for i,c in rr.items():b=plus(b,mul(ex(c),ex(V29[target].get(i,Q(0)))))
  dots.append(list(map(str,b)));a=plus(ex(Q(d['known_alpha_beta_on_Wv'][ch])),mul(ex(Q(3375,11)),ex(Q(d['known_alpha23_beta23'][ch]))));a=plus(a,mul(ex(Q(50)),b));recover.append(list(map(str,mul(ex(Q(2,25)),a))))
 for target in range(3):
  for i in set(VW[target])|set(V23[target])|set(V29[target])|set(V30[target]):
   b=ex(VW[target].get(i,Q(0)))
   for c,v in [(Q(-3375,11),V23),(Q(-50),V29),(Q(25,2),V30)]:b=plus(b,mul(ex(c),ex(v[target].get(i,Q(0)))))
 assert dots==d['two_new_dot_positive_bounds'] and recover==d['two_recovered30_positive_bounds']
 assert maximum==r['source_committed_bit_bound']==76 and 3*maximum+1==r['source_temporary_bit_bound']==229
 for target in range(3):
  for i in range(33):assert Q(-3375,11)*V23[target].get(i,Q(0))-50*V29[target].get(i,Q(0))+Q(25,2)*V30[target].get(i,Q(0))==VW[target].get(i,Q(0))
 assert 2*count(d)+5000<=r['whole_numeric_bound']==5222 and r['whole_RAW_bound']==1090666
 cp='reports/AMBIENT_B18_ELL2_T2_ALPHA29_BETA29_PREPARED_CUSTODY_20261105_v63/report.json';cust=read(cp);wp=ROOT/cust['worker'];assert sha(wp)==cust['worker_sha256']=='f87cc6b354f0b4946ce8bd6fe38e7aa33c24a121f8f541ef99adcf17dc165d96'
 for p,h in cust['worker_and_helper_sha256'].items():assert sha(ROOT/p)==h
 assert cust['source_gate_sha256']==sha(ROOT/rp) and cust['safe_help_literal_PASS'] and cust['safe_help_exit']==0 and cust['exact_two_new_pairings_only']
 text=wp.read_text();assert text.count("d['primitive_V30_raw63_mask2_source_check_only']")==0 and "for ch,row in enumerate(rho)" in text
 inputs['literal_two_dot_worker']={'path':cust['worker'],'sha256':sha(wp)}
 proof={'universal_domain':'ForanyconstantvwithM0v=0,F=MvhasF0=0andG1F=Vv. Thepeered286-variablemixedadjointidentityisuniversalinF1/F2,independentoftheoriginalchoiceWv. Thecomponentwise3960raw33factorproofislikewiseindependentofv. Thusrho(Ve29)givesthesamegauge-invariantalpha/betaforv=e29.',
 'literal_source_vector':'Column29isphysicalrow63masktauAinstandardfree4basis:F63*tauAprojectsfrom(delta,bA,-bB,0)to(delta,0,-bB);theconstantnuA*bAcomponentisremovedwiththecommonunitbasegenerator. NoF9coordinateisusedastracevariable.',
 'recovery':'All99cachedfunction-vectorcoefficientsverifyVW=(-3375/11)V23-50V29+(25/2)V30. Rholinearitythenprovesu30=(2/25)(uW+(3375/11)u23+50u29). NoV30dotmustbeevaluated.',
 'scope':'OnlytwoM1invariantpairingsandlinearrecovery.Canonicalothermuvaluesandallupper619clearing/pivot/G2/sourcecoefficientsremainseparate.Nopsi/HF2/lengthclaim.'}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);e=out/'independent-universal29-source-vector-76-prior-and-two-dot-custody.json';e.write_text(json.dumps({'M0_all3_point_kernel_directions_exact':True,'literal_V29_raw63_mask1_verified':True,'all99_Wv_vector_identities_verified':True,'both76bit_positive_dot_and_recovery_bounds_rebuilt':True,'no_pairing_evaluated':True,'proof':proof},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'source_gate_sha256':sha(ROOT/rp),'two_dot_source_and_worker_condition_closed':True,'no_pairing_or_new_response_value':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{e.name:sha(e)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'whole_source_numeric':count(d),'gate_peer_custody_numeric':count([r,report,cust])}))
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Independent fixed-point and formal-order audit of the upper619 clearing."""
import hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_UPPER619_SOURCE_CLEARING_20261105_v63';STEM=BASE+'_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sp(v):return {i:Q(x) for i,x in v if Q(x)}
def dot(a,b):return sum((x*b.get(i,Q(0)) for i,x in a.items()),Q(0))
def image(A,b):
 z={}
 for j,c in b.items():
  for i,x in A[j].items():
   v=z.get(i,Q(0))+c*x
   if v:z[i]=v
   else:z.pop(i,None)
 return z
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));pins={}
 def read(p,keys=None):pins[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and pins[rp]['sha256']=='f6236b432b2db7aca9177e0ee4cf7391890cb2ec6ecae16c7773b6c6077dd920'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 ep=f'reports/{BASE}/exact619-upper-source-clearing-and-order-frontier.json';e=read(ep);assert pins[ep]['sha256']==r['evidence_sha256'][Path(ep).name]
 assert sha(ROOT/'code/ambient_B18_ell2_T2_upper619_source_clearing_20261105_v63.py')==r['code_sha256'] and sha(ROOT/'research/notes'/f'{BASE}.md')==r['note_sha256']
 for k,b in e['source_bindings'].items():
  assert sha(ROOT/b['path'])==b['sha256'];v=read(b['path'],[b['field']])[b['field']]
  if k in ['original_full_R0_X768','actual_Y768']:assert len(v)==768
  elif k in ['actual_comparison_M','complete_minimal6','first27_unit','next99_unit']:assert isinstance(v,int)
  else:assert len(v)==b['shape'][0]
  del v
 assert read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63_REVIEW/report.json')['all_passed']
 assert read('reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_RANK_20261105_v63_REVIEW/report.json')['all_passed']
 mu='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json'
 mg=read(mu,['full_G2_linear234_by102_columns','Mnew0_columns']);G=list(map(sp,mg['full_G2_linear234_by102_columns']));M=list(map(sp,mg['Mnew0_columns']));assert len(G)==102 and all(not M[i] for i in [23,29,30])
 up='reports/AMBIENT_B18_ELL2_T2_UPPER2_CUBIC_OUTPUT_SOURCE_20261105_v63/fixed-upper2-one-cubic-output-identity-and-source-frontier.json';u=read(up);lam=list(map(sp,u['lambda_two']));s0=sp(u['s0_on_F2_source102']);s1=sp(u['s1_on_F2_source102'])
 GR=[{i//9:x for i,x in c.items() if i%9==7} for c in G];GE=[{i//9:x for i,x in c.items() if i%9==8} for c in G]
 for j in range(102):assert dot(lam[0],GR[j])==0 and dot(lam[0],GE[j])==dot(lam[1],GR[j])==s0.get(j,Q(0)) and dot(lam[1],GE[j])==s1.get(j,Q(0))
 sy=list(map(sp,read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json',['actual_selected_Y0_columns102'])['actual_selected_Y0_columns102']));assert len(sy)==99
 # Independent one-denominator integer replay of all198 point annihilations.
 receipts=[]
 for channel,s in enumerate([s0,s1]):
  for j,col in enumerate(sy):
   terms=[a*col.get(i,Q(0)) for i,a in s.items()];D=math.lcm(*(x.denominator for x in terms));N=sum(x.numerator*(D//x.denominator) for x in terms);assert N==0
   receipts.append([channel,j,str(D),str(N)])
 assert e['fixed_point_fixture']['lambda_two']==u['lambda_two'] and list(map(sp,e['fixed_point_fixture']['s0_s1_two_rows']))==[s0,s1]
 # Reconstruct the upper extraction in the accepted known-linear NF, rather
 # than infer it from the preceding s0/s1 annihilations.
 pp='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json';Pi=list(map(sp,read(pp,['PiY3_from102_columns'])['PiY3_from102_columns']))
 lp='reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json';l=read(lp,['degree2_monomial_order44','all264_NF_root_nodes'])
 np='reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_DEGREE2_20261105_v63/actual432-overlap-rank-and-known-linear-cover-degree2.json';memo={i:sp(v) for i,v in read(np,['all264_linear_NF_columns'])['all264_linear_NF_columns']};ix={tuple(m):i for i,m in enumerate(l['degree2_monomial_order44'])}
 rr=[0]*9;rr[7]=2;re=[0]*9;re[7]=re[8]=1;top=[[memo[l['all264_NF_root_nodes'][i][ix[tuple(m)]]] for m in [rr,re]] for i in range(3)];Q2=list(map(sp,u['Q_two6']));assert Q2==top[2]
 for j in range(102):
  for mon in range(2):
   for k in range(6):
    actual=sum((c*top[i][mon].get(k,Q(0)) for i,c in Pi[j].items()),Q(0));expected=Q2[1].get(k,Q(0))*s0.get(j,Q(0)) if mon else Q2[0].get(k,Q(0))*s0.get(j,Q(0))+Q2[1].get(k,Q(0))*s1.get(j,Q(0));assert actual==expected
 for target in range(6):
  for z in range(7):
   for v in range(9):
    m=[0]*9;m[z]+=1;m[v]+=1;assert not memo[l['all264_NF_root_nodes'][target][ix[tuple(m)]]]
 mini=read('reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json',['free30_V_indices','fixed_bottom30_basis_Xp27_then_Efree3','fixed_bottom30_inverse_columns','remaining_free3_bottom_coordinates','remaining_free3_originalV_rows'])
 BB=list(map(sp,mini['fixed_bottom30_basis_Xp27_then_Efree3']));Bi=list(map(sp,mini['fixed_bottom30_inverse_columns']));assert len(BB)==len(Bi)==30
 for j in range(30):assert image(BB,Bi[j])==image(Bi,BB[j])=={j:Q(1)}
 assert mini['remaining_free3_originalV_rows']==[23,29,30] and [mini['free30_V_indices'][i] for i in mini['remaining_free3_bottom_coordinates']]==[23,29,30]
 retained=e['retained_unknowns'];assert any('b1' in x for x in retained) and any('M1' in x for x in retained) and any('native L2' in x for x in retained)
 assert e['all619_values_uncomputed'] and not r['remaining619_values_computed']
 proof={'moving_clearing':'Use the actual27 moving bottom chart and actual99 upper unit, with a single coherent M,Y frame. G2G3=0 gives exactly G2U=ME. The accepted full first identities and G3rank25 clearing imply U0=U1=B0=B1=0. Moving pivot coordinates vanish on E. Point splitting then forces E0=E1=0; degree2 of the identity gives M0E2=0, so E2=H0B2.',
 'cubic_identity':'Taking degree3 yields G2linU2=M0E3+M1E2; G2quadU1 and M2E1 vanish for the stated exact leading-zero reason. The independent NF identities extract only r^2 ell, with ordinary coefficients.',
 'coefficient_expansion':'A=ac-astar*b has A1=a1; E2=P2-Xp0A2 with P2=Xc2-Xp1a1-Xstar1b1. Thus A2=P27P2,B2=Pi3P2. Degree3 is Xc3-Xp2a1-Xstar2b1-Xp1A2-Xstar1b2-Xp0A3. The A3 term is killed by M0Xp0=0. Each b2 contraction is killed by the198 explicitly verified s0/s1Ystar0 identities.',
 'frame_and_counterconditions':'No canonical M1 is substituted into the normalized L48 frame; M+G2A,Y+AX must be transported together. b1 uses the full918-domain first vector. Unknown Yp0,A2,M1 and native L2 contributions remain. Individual effective native representatives cannot be substituted without their own S0 condition. No associated-graded exactness or higher coefficient availability is asserted.'}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ve=out/'independent198-fixtures-upperNF-and-exact-clearing-proof.json';ve.write_text(json.dumps({'all198_integer_lattice_annihilations':receipts,'all408_lambda_G2plane_identities':True,'all1224_upper_known_NF_identities':True,'all378_z_quadratic_NF_zero':True,'point30_chart_both_inverse_products':True,'H0_axes_M0_kernel_literal':True,'exact_formal_proof':proof,'no_remaining619_or_new_source_coefficients_evaluated':True},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'primary_sha256':pins[rp]['sha256'],'exact_upper_source_clearing_accepted':True,'all_unknown_Yp_b1_A2_M1_nativeL2_retained':True,'no_upper619_values_or_cache_sufficiency_claim':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ve.name:sha(ve)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'seconds':report['seconds']}))
if __name__=='__main__':main()

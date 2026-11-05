#!/usr/bin/env python3
"""Independent pruned-first/source and positive M1-operator input review."""
import json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count,image
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_NONACTIVE9_M1_FIRST_INPUT_GATE_V2_20261105_v63';STEM=BASE+'_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p,keys=None):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and not r['signed_first_M1_or_invariant_values_evaluated']
 assert sha(ROOT/'code/ambient_B18_ell2_T2_nonactive9_M1_first_gate_v2_20261105_v63.py')==r['code_sha256']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 ep=f'reports/{BASE}/pruned-known-firsts-and-positive-M1-four-invariant-test-input.json';d=read(ep);assert sha(ROOT/ep)==r['evidence_sha256'][Path(ep).name]
 known=read('reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_KNOWN_INPUT_GATE_20261105_v63/pruned-known619-lower4-input-and-prior-resource-gate.json')
 mini=read('reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json')
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json'
 mm=read(mug,['Mnew0_columns','full_G2_linear234_by102_columns'])
 C={j:sp(v) for j,v in d['selected99_original768']};Y={j:sp(v) for j,v in d['selected99_Y0']};R={j:sp(v) for j,v in d['original_r_first_columns']};E={j:sp(v) for j,v in d['original640_ell_first_columns']}
 used=d['used_selected99_indices'];assert len(used)==len(C)==34 and sorted(used)==sorted(C)
 assert [d['free30_V_indices'],d['fixed_B30_inverse_columns'],d['point_Xp27_original_indices']]==[mini['free30_V_indices'],mini['fixed_bottom30_inverse_columns'],mini['original27_X_source_columns']]
 assert d['M0_columns']==mm['Mnew0_columns'] and d['G2linear_columns']==mm['full_G2_linear234_by102_columns']
 for j in used:assert C[j]==sp(known['selected99_original768_combinations'][j]) and Y[j]==sp(known['selected99_Y0_columns102'][j])
 assert all(v==sp(known['original768_r_first_columns'][j]) for j,v in R.items())
 assert all(v==sp(known['original640_ell_first_columns'][j]) for j,v in E.items())
 assert d['native_L0_columns']==known['native_L0'] and d['native_S0_columns']==known['native_S0']
 assert all(sp(v)==sp(known['native_Se'][j]) for j,v in d['native_Se_columns'])
 Xp={}
 for off,path,key in [(0,'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json','source_columns640_rows56_sparse_dual9'),(640,'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json','source_columns128_rows48_sparse_dual9')]:
  raw=read(path,[key])[key]
  for j in d['point_Xp27_original_indices']:
   if off<=j<off+len(raw):Xp[j]={row:sp(v)[0] for row,v in raw[j-off] if sp(v).get(0)}
 assert list(map(sp,d['point_Xp27_columns']))==[Xp[j] for j in d['point_Xp27_original_indices']]
 support=read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_UPPER_SUPPORT_PREFLIGHT_20261105_v63/eighteen-formal-weights-and-coefficient-free-source-ancestor-inventory.json')
 assert d['formal_eighteen_weights']==[{'original418_source':x['original418_source'],'contraction':x['contraction'],'weights':x['formal_weights_on99psi2']} for x in support['eighteen_formal_contractions']]
 aux=read('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_INPUT_GATE_20261105_v63/actual-eight-source-root-folds-M1-operators-and-positive-AD-budget.json',['lambda26_two_channels'])
 assert d['lambda_two']==aux['lambda26_two_channels']
 graph=read('reports/AMBIENT_B18_ELL2_T2_LINEAR_GRAPH4_SOURCE_20261105_v63/actual23-linear-graph-four-common-contractions-and-six-adjoint-input.json')
 assert d['common_Wv56']==graph['actual_H0v_source56'] and sp(d['common_Wv56'])=={23:Q(-3375,11),29:Q(-50),30:Q(25,2)}
 maximum=0
 def mark(a):
  nonlocal maximum
  maximum=max(maximum,a[0].bit_length(),a[1].bit_length());return a
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator)))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def accum(v,k,b):
  if b[1]:v[k]=plus(v.get(k,(1,0)),b)
 def vec(v):return {i:ex(x) for i,x in v.items()}
 def padd(v,w,c=(1,1)):
  q=dict(v)
  for i,b in w.items():accum(q,i,mul(c,b))
  return q
 def apply(A,v):
  q={}
  for j,b in v.items():
   for i,x in A[j].items():accum(q,i,mul(ex(x),b))
  return q
 L=list(map(sp,d['native_L0_columns']));S0=list(map(sp,d['native_S0_columns']));Se={j:sp(v) for j,v in d['native_Se_columns']};inv=list(map(sp,d['fixed_B30_inverse_columns']));piv=list(map(sp,d['point_Xp27_columns']));M=list(map(sp,d['M0_columns']));G2=list(map(sp,d['G2linear_columns']));lam=list(map(sp,d['lambda_two']));W=sp(d['common_Wv56']);free=d['free30_V_indices']
 EB={j:vec(v) for j,v in E.items()}
 for j,v in Se.items():EB[640+j]=apply(L,vec(v))
 Z={}
 for j,c in C.items():
  cn={k-640:x for k,x in c.items() if k>=640};assert not image(S0,cn);apply(S0,vec(cn))
  zz=[]
  for direction,raw in enumerate([R,EB]):
   v={}
   for k,x in c.items():v=padd(v,vec(raw[k]) if direction==0 else raw[k],ex(x))
   a=apply(inv,{i:v[row] for i,row in enumerate(free) if row in v});a={i:b for i,b in a.items() if i<27};z=padd(v,apply(piv,a));zz.append(z)
   apply(M,z);rhs=[{i//9:x for i,x in v.items() if i%9==(7+direction)} for v in G2];apply(rhs,vec(Y[j]))
  Z[j]=zz
 p,s=d['fixed_extraction_lambda_row'],d['fixed_extraction_source_row'];assert (p,s)==(12,29) and lam[0][p]*W[s]
 for n,rec in enumerate(d['formal_eighteen_weights']):
  v=[[{},{}],[{},{}]]
  for flat,x in rec['weights']:
   j,ch=divmod(flat,2)
   for direction in range(2):v[ch][direction]=padd(v[ch][direction],Z[j][direction],ex(x))
  assert [[sorted(q) for q in row] for row in v]==d['positive_weighted_vector_supports'][n]
  op={}
  for ch,row in enumerate(lam):
   for target,x in row.items():
    for direction,sourcev in [(0,v[ch][1]),(1,v[ch][0])]:
     for source,b in sourcev.items():accum(op,(direction,target,source),mul(ex(x),b))
  lp,lp1=lam[0][p],lam[1].get(p,Q(0))
  bw=mul(op.get((1,p,s),(1,0)),ex(1/(lp*W[s])))
  be=plus(mul(op.get((1,p,23),(1,0)),ex(1/lp)),mul(bw,ex(W[23])))
  aw=mul(plus(op.get((0,p,s),(1,0)),mul(ex(lp1),mul(bw,ex(W[s])))),ex(1/(lp*W[s])))
  ae=plus(mul(plus(op.get((0,p,23),(1,0)),mul(ex(lp1),plus(mul(bw,ex(W[23])),be))),ex(1/lp)),mul(aw,ex(W[23])))
  assert [[str(D),str(U)] for D,U in [aw,bw,ae,be]]==d['positive_four_coefficient_bounds'][n]
  for target,x in lam[0].items():
   for source,w in W.items():accum(op,(0,target,source),mul(aw,mul(ex(x),ex(w))));accum(op,(1,target,source),mul(bw,mul(ex(x),ex(w))))
   accum(op,(0,target,23),mul(ae,ex(x)));accum(op,(1,target,23),mul(be,ex(x)))
  for target,x in lam[1].items():
   for source,w in W.items():accum(op,(0,target,source),mul(bw,mul(ex(x),ex(w))))
   accum(op,(0,target,23),mul(be,ex(x)))
 assert maximum==r['positive_prior_bits']==197 and r['future_temporary_prior_bits']==592
 assert count(d)==3307 and 2*count(d)+8000==14614<r['future_whole_numeric_bound']==16454 and r['future_RAW_bound']==4632640
 cp='reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_M1_FIRST_PREPARED_CUSTODY_20261105_v63/report.json';cust=read(cp);wp=ROOT/cust['worker'];assert sha(wp)==cust['worker_sha256']=='121fdb2b2c75df24e7d468d896a75c956448dc3ae7379f4418435b1eea09102a'
 assert cust['source_gate_sha256']==sha(ROOT/rp) and cust['prepared_not_executed_needs_ROOT_GO_and_FINAL_sourcepeer'] and cust['safe_help_exit']==0
 for h in cust['helpers_with_unchanged_hashes_and_previously_verified_safe_help']:assert sha(ROOT/h['path'])==h['sha256']
 inputs['prepared_known_first_worker']={'path':cust['worker'],'sha256':sha(wp)}
 err='reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_M1_FIRST_PRIOR_ERRATUM_20261105_v63/report.json';read(err)
 proof={'known_first_scope':'Only34cachedselected99sourcefirstsarefolded. IndividualnativeellLeS0isnotdiscarded: SeiseffectiveonlyaftereachfullselectedsourceS0cNzero. Actual34fullchainidentitiesand2912-slotrepresentationtestsarefutureworkerchecks.',
 'M1_operator':'Forch0/1,lambda_ch(Mr*Xeff_ell+Mell*Xeff_r). OriginalM1valuesremainindependentvariables.',
 'four_basis':'alphaW=lambda0MrW;betaW=lambda0MellW-lambda1MrW;alpha23=lambda0Mre23;beta23=lambda0Melle23-lambda1Mre23.',
 'correct_extraction':'UseMell(p,29)torecoverbetaW;Mell(p,23)forbeta23;THENMr+lambda1beta-carryforalphaW/alpha23. No r/ellbetaslotinterchange. All2912slotsmustdecideeachproposedrepresentation,notfourprobevaluesalone.',
 'W2_scope':'Samea1chartandPi3mapgiveM0W2=0,Pi3W2=0bylower9zero,pivot27W2=0;thefrozen56coordinatesplitisinjective. ThisprovesplaneW2zerojointlyfor9,notindividualpsi/M1zero.',
 'no_invariant_values':'Thetestaskscoefficientrowspanonly.Evenasuccessdoesnotevaluatealpha/beta,alph23/beta23,psi,HForT2length.'}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);e=out/'independent34-first-source197-prior-and-M1-four-basis-method.json';e.write_text(json.dumps({'all_pruned_original_firsts_and_sameB30_source_bound':True,'all34_nativeS0_combinations_zero':True,'positive197_and_all18_four_coefficient_bounds_rebuilt':True,'serialized_numeric_with_metadata':3307,'two_input_copies_plus8000':14614,'no_signed_weighted_M1_operator_or_invariant_value':True,'proof':proof},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'source_gate_sha256':sha(ROOT/rp),'known_first_worker_method_and_custody_accepted':True,'no_signed_M1_operator_or_values':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{e.name:sha(e)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'prior':197,'report_sha256':sha(out/'report.json'),'gate_peer_numeric':count([r,report])}))
if __name__=='__main__':main()

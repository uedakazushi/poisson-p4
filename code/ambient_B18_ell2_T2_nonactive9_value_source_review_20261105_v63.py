#!/usr/bin/env python3
"""Fresh whole9 cache/frame/SUM/positive/rank-method source review."""
import json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count,image
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_NONACTIVE9_VALUE_INPUT_GATE_20261105_v63';STEM=BASE+'_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p,keys=None):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and not r['value_job_started']
 assert r['code_sha256']==sha(ROOT/'code/ambient_B18_ell2_T2_nonactive9_value_gate_20261105_v63.py')
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 ep=f'reports/{BASE}/whole-old-native-pivot-nine-classes-positive-source-input.json';d=read(ep);assert sha(ROOT/ep)==r['evidence_sha256'][Path(ep).name]
 for stem in ['NONACTIVE9_LINEAR_AD','NONACTIVE9_M1_FIRST','NONACTIVE9_LOWER_VALUES']:
  assert read('reports/AMBIENT_B18_ELL2_T2_'+stem+'_20261105_v63_REVIEW/report.json')['all_passed']
 first=read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_M1_FIRST_INPUT_GATE_V2_20261105_v63/pruned-known-firsts-and-positive-M1-four-invariant-test-input.json')
 ad=read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_LINEAR_AD_20261105_v63/actual-eighteen-formal-old-second-sensitivities-and-native-source-combinations.json')['all18_formal_source_results'];assert d['old18_formal_results']==ad
 old=read('reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_20261105_v63/actual104-source-rell-values-and-two-by619-old-knownparts.json')
 CV={i:Q(v[3]) for i,v,grad in old['all104_original_C0_values_and_opaque_gradients']};needed=sorted({i for x in ad for i,v in x['exact_independent_original_C0_second_sensitivity']})
 assert d['literal_old104_second_values']==[[i,str(CV[i])] for i in needed] and len(needed)==10
 original='reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63/whole-cache-eight-psi-four-quadratic-classes-source-input.json'
 keys=['native_L0','native_Lr','native_Le','native_S0','native_Sr','native_Se','native_S2rell','pivot27_source_indices','pivot27_r_columns','pivot27_ell_recipes','Gamma_two','lambda_two','Q_two6_columns','Q_fixed_minor_rows','Q_fixed_minor_determinant']
 src=read(original,keys);assert all(d[k]==src[k] for k in keys)
 assert d['native_L0']==first['native_L0_columns'] and d['lambda_two']==first['lambda_two'] and d['pivot27_source_indices']==first['point_Xp27_original_indices']
 assert all(d['carry_first_inputs'][k]==first[k] for k in d['carry_first_inputs'])
 # Exact source of the finite-resolution cancellation, not extrapolation
 # from first-degree fixtures.
 cancel=read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_G2_MIXED_ZERO_SOURCE_20261105_v63/report.json');assert cancel['all_passed']
 for p,h in cancel['input_sha256'].items():assert sha(ROOT/p)==h
 ff=read('reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json',['actual_full_G3linear102_by418_columns'])
 G2=list(map(sp,first['G2linear_columns']));G3=list(map(sp,ff['actual_full_G3linear102_by418_columns']));lam=list(map(sp,d['lambda_two']))
 assert all(set(v)<=set(range(12,16)) for v in lam)
 for iz,z in enumerate([2,3,6]):
  for j in range(17):
   active={k:x for k,x in G2[42+17*iz+j].items() if k//9<17};assert active=={9*j+z:Q(-1)}
 for col in G3[238:]:assert all(42<=k//9<93 for k in col if k%9 in (7,8))
 assert cancel['allorder_G2_quadratic_zero_follows_exact_tensor_source_not_tangent_inference']
 maximum=0
 def mark(v):
  nonlocal maximum
  maximum=max(maximum,v[0].bit_length(),v[1].bit_length());return v
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator)))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  DD=math.lcm(a[0],b[0]);return mark((DD,a[1]*(DD//a[0])+b[1]*(DD//b[0])))
 def accum(v,k,b):
  if b[1]:v[k]=plus(v.get(k,(1,0)),b)
 def vec(v):return {i:ex(x) for i,x in v.items()}
 def padd(a,b):
  q=dict(a)
  for i,x in b.items():accum(q,i,x)
  return q
 def apply(A,v):
  q={}
  for j,b in v.items():
   for i,x in A[j].items():accum(q,i,mul(ex(x),b))
  return q
 def dot(row,v):
  q=(1,0)
  for i,x in row.items():q=plus(q,mul(ex(x),v.get(i,(1,0))))
  return q
 L0,Lr,Le,S0,Sr,Se,S2=[list(map(sp,d[k])) for k in keys[:7]];Gam=list(map(sp,d['Gamma_two']));QQ=list(map(sp,d['Q_two6_columns']))
 rows=d['Q_fixed_minor_rows'];det=QQ[0].get(rows[0],Q(0))*QQ[1].get(rows[1],Q(0))-QQ[0].get(rows[1],Q(0))*QQ[1].get(rows[0],Q(0));assert det==Q(d['Q_fixed_minor_determinant']) and det
 pivotE=[]
 for mode,v in d['pivot27_ell_recipes']:
  if mode=='literal':pivotE.append(vec(sp(v)))
  else:
   assert mode=='native'
   # Independent semantic SUM, including both sets of partial sums.
   pivotE.append(padd(apply(L0,vec(Se[v])),apply(Le,vec(S0[v]))))
 C={j:sp(v) for j,v in first['selected99_original768']};R={j:sp(v) for j,v in first['original_r_first_columns']};E={j:vec(sp(v)) for j,v in first['original640_ell_first_columns']};inv=list(map(sp,first['fixed_B30_inverse_columns']));free=first['free30_V_indices'];pivR=list(map(sp,d['pivot27_r_columns']));A={}
 for j in {k for c in C.values() for k in c if k>=640}:E[j]=apply(L0,vec(Se[j-640]))
 for j,c in C.items():
  cn={k-640:x for k,x in c.items() if k>=640};assert not image(S0,cn);apply(S0,vec(cn));pair=[]
  for direction in (0,1):
   v={}
   for k,x in c.items():
    for row,b in (vec(R[k]) if direction==0 else E[k]).items():accum(v,row,mul(ex(x),b))
   a=apply(inv,{i:v[row] for i,row in enumerate(free) if row in v});pair.append({i:b for i,b in a.items() if i<27})
  A[j]=pair
 psi=[]
 for n,(rec,wrec,bound) in enumerate(zip(ad,first['formal_eighteen_weights'],d['all18_component_positive_bounds'],strict=True)):
  assert (rec['original418_source'],rec['contraction'])==(wrec['original418_source'],wrec['contraction'])==(bound['source418'],bound['contraction'])
  aa=[[{},{}],[{},{}]]
  for flat,x in wrec['weights']:
   j,ch=divmod(flat,2)
   for direction in range(2):
    for i,b in A[j][direction].items():accum(aa[ch][direction],i,mul(ex(x),b))
  oldb=(1,0)
  for i,x in rec['exact_independent_original_C0_second_sensitivity']:oldb=plus(oldb,mul(ex(x),ex(CV[i])))
  nativebounds=[];pivotbounds=[]
  for ch,records in enumerate(rec['exact_original_native128_combinations_two_channels']):
   cn=sp(records);assert not image(S0,cn);cb=vec(cn);apply(S0,cb)
   # All THREE genuine second coefficient products, each added as a sum.
   native=padd(padd(apply(L0,apply(S2,cb)),apply(Lr,apply(Se,cb))),apply(Le,apply(Sr,cb)))
   nativebounds.append(dot(Gam[ch],native));ar,ae=aa[ch];pv=apply(pivR,ae)
   for i,b in ar.items():
    for row,x in pivotE[i].items():accum(pv,row,mul(b,x))
   pivotbounds.append(dot(Gam[ch],pv))
  total=oldb
  for b in nativebounds+pivotbounds:total=plus(total,b)
  expected={'source418':rec['original418_source'],'contraction':rec['contraction'],'old':list(map(str,oldb)),'native_channels':[list(map(str,b)) for b in nativebounds],'pivot_channels':[list(map(str,b)) for b in pivotbounds],'psi':list(map(str,total))}
  assert expected==bound;psi.append(total)
 for j in range(9):
  for row in range(6):plus(mul(ex(QQ[0].get(row,Q(0))),psi[2*j]),mul(ex(QQ[1].get(row,Q(0))),psi[2*j+1]))
 assert maximum==d['nonrank_positive_bound_bits']==r['nonrank_positive_bound_bits']==204
 DD=math.lcm(*(D for D,U in psi));UU=max(U*(DD//D) for D,U in psi)
 rankbits=max(DD.bit_length(),3*UU.bit_length()+8)
 assert str(DD)==d['coordinate_common_D'] and str(UU)==d['coordinate_abs_clear_U']
 assert rankbits==d['integer_rank_and_factor_bound_bits']==r['whole_value_height_bound_bits']==362<384
 assert max(3*maximum+1,rankbits)==r['temporary_value_height_bound_bits']==613<1153
 assert 2*count(d)+5000<r['whole_value_numeric_bound']==21856 and r['whole_RAW_bound']==3886640
 cust=read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_VALUES_PREPARED_CUSTODY_20261105_v63/report.json');wp=ROOT/cust['worker'];assert sha(wp)==cust['worker_sha256']=='9234ede81240f2992b20fdeef7453517b8b4bb466aecc72524b889c2899f3669'
 for p,h in cust['worker_and_helper_sha256'].items():assert sha(ROOT/p)==h
 assert cust['source_gate_sha256']==sha(ROOT/rp) and cust['safe_help_literal_PASS'] and cust['safe_help_exit']==0
 inputs['prepared_complete9_worker']={'path':cust['worker'],'sha256':sha(wp)}
 proof={'formula':'psi=old104+Gamma(L0S2c+LrSec+LeSrc)-Gamma(Xpr*aell+Xp_ell*ar). Xp_ell=L0Se+LeS0withSUM,notmaxunion. S0c=0eligibleonlyafterwholeweightedsourcecheck.',
 'source_cancellations':'M1operatorisliteral0(df568). ExacttensorG2mixedactive17block=-ziIhasnoquadandlambdas12..15annihilateremainingtargets17..25. NonactiveG3rawhigherPiY0projectionandrawbottomare0. W2plane0andz*mNF0removecoherentnormalizerjointly,nocanonicalunknownmu0assignment.',
 'rank_bound':'Ranktestusesinteger2-vectorsboundedU. Minor/Nnumeratorsorder<=2:<=2Ubits+1. Clearedfactor/kernelverificationusesatmost3Ubits+log9;3Ubits+8=362coversALLintegerpartials. RationalF(N,minor)isconstructedonlyafterintegeridentities,notgenericfractionelimination. Thusranktemporary362iscoveredbyoverall613.',
 'scope':'actualclass=-Qpsiwithfixednonzerominor.Full9factors/kernelandsixcoordinatechecksarefutureworkerrequirements.RankfromthisfamilyboundsHF2aboveonly;remaining619andtotalstayopen.'}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);e=out/'independent-whole9-literal-cache-SUM-prior-and-frame-proof.json';e.write_text(json.dumps({'all_literal_cache_payload_fields_exact':True,'all18_component_positive_envelopes_rebuilt':True,'individual_pivot_ell_is_actual_two_product_SUM':True,'all3_native_coefficient_products_retained':True,'nonrank_prior':204,'rank_committed_prior':362,'whole_temporary_prior':613,'prepared_worker_sha256':sha(wp),'no_actual_psi_or_class_values':True,'proof':proof},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'source_gate_sha256':sha(ROOT/rp),'whole_complete9_source_and_worker_condition_closed':True,'no_actual_psi_class_or_rank_values':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{e.name:sha(e)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'gate_peer_custody_numeric':count([r,report,cust]),'input_numeric':count(d)}))
if __name__=='__main__':main()

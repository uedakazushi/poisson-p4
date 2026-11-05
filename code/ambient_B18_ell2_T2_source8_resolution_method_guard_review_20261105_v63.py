#!/usr/bin/env python3
"""Independent positive prior-bound and finite-resolution method review."""
import itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_G3_quadratic_source_review_20261105_v63 import sha,sp
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_SOURCE8_RESOLUTION_METHOD_GUARD_20261105_v63_REVIEW'
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 inputs={};reports={}
 for stem,code in [('AMBIENT_B18_ELL2_T2_SOURCE8_RESOLUTION_PRIOR_HEIGHT_20261105_v63','ambient_B18_ell2_T2_source8_resolution_prior_height_20261105_v63.py'),('AMBIENT_B18_ELL2_T2_SOURCE8_RESOLUTION_VERIFY_HEIGHT_20261105_v63','ambient_B18_ell2_T2_source8_resolution_prior_verify_bound_20261105_v63.py')]:
  p=ROOT/'reports'/stem/'report.json';r=json.loads(p.read_text());assert r['all_passed'];reports[stem]=r;inputs[stem]={'path':str(p.relative_to(ROOT)),'sha256':sha(p)}
  assert sha(ROOT/'code'/code)==r['code_sha256'] and sha(ROOT/'research/notes'/f'{stem}.md')==r['note_sha256']
  for q,h in r['input_sha256'].items():assert sha(ROOT/q)==h;inputs[q]={'path':q,'sha256':h}
 t=json.loads((ROOT/'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63/actual-active-unit-bases-and26-frame-transport.json').read_text())
 oldpath=ROOT/'reports/AMBIENT_B18_ELL2_T2_SOURCE8_RESOLUTION_PRIOR_HEIGHT_20261105_v63/prior-input-receipt-product-and-partialsum-height-envelopes.json';old=json.loads(oldpath.read_text());maximum=0
 def mark(a):
  nonlocal maximum
  maximum=max(maximum,a[0].bit_length(),a[1].bit_length());return a
 def ex(x):q=Q(x);return mark((q.denominator,abs(q.numerator)))
 def plus(a,b):
  d=math.lcm(a[0],b[0]);return mark((d,a[1]*(d//a[0])+b[1]*(d//b[0])))
 def times(a,b):return mark((a[0]*b[0],a[1]*b[1]))
 def apply(cols,v):
  out={}
  for j,b in v.items():
   for i,a in cols[j].items():out[i]=plus(out.get(i,(1,0)),times(ex(a),b))
  return out
 def receipt_bounds(records):
  norms={};reps={};bound=0
  for r in records:
   v=(1,1)
   for p,c in r['subtract']:v=plus(v,times(ex(c),reps[p]))
   v=times(v,ex(1/Q(r['divide_by'])));norms[r['pivot']]=sp(r['normalized']);reps[r['pivot']]=v;bound=max(bound,v[0].bit_length(),v[1].bit_length())
  return norms,reps,bound
 def extract(v,norm,reps):
  q=dict(v);ans=(1,0)
  for p in sorted(norm):
   x=q.pop(p,(1,0))
   if not x[1]:continue
   ans=plus(ans,times(x,reps[p]))
   for row,a in norm[p].items():
    if row!=p:q[row]=plus(q.get(row,(1,0)),times(x,ex(a)))
  return ans
 def corrected(v,Pi,free):
  q=dict(v)
  for i,b in apply(Pi,v).items():q[free[i]]=plus(q.get(free[i],(1,0)),b)
  return q
 def raw(degree):
  active=[0,1,4,5,7];target=[w for q in range(degree) for w in itertools.combinations(active,q)];look={w:i for i,w in enumerate(target)};out=[]
  for q in range(degree+1):
   for w in itertools.combinations(active,q):
    p=degree-q
    for j in range(14):
     v={}
     if p:v[(14*look[w]+j)*6+5]=Q(1)
     for s,z in enumerate(w):v[(14*look[w[:s]+w[s+1:]]+j)*6+active.index(z)]=Q((-1)**(p+s))
     out.append(v)
  return out
 n1,r1,b1=receipt_bounds(t['active_F1_fixed67_receipts']);n2,r2,b2=receipt_bounds(t['active_F2_fixed182_receipts']);Pi17=list(map(sp,t['active_Pi17_columns']));Pi42=list(map(sp,t['active_Pi42_columns']))
 g2=json.loads((ROOT/'reports/AMBIENT_B18_ELL2_C2_G2_LINEAR_SOURCE_PROFILE_20261105_v63/actual-G2-linear-fixed-source-profile.json').read_text());I2=list(map(sp,g2['actual_I2_source_origin42_columns']))
 p=ROOT/'reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json';I3=list(map(sp,fields(p,['actual_active_I3_origin_columns'])['actual_active_I3_origin_columns']))
 fold=json.loads((ROOT/'reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_INPUT_GATE_20261105_v63/actual-eight-source-root-folds-M1-operators-and-positive-AD-budget.json').read_text());four=json.loads((ROOT/'reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE_PROFILE_20261105_v63/actual-four-active-source-complements-and-unknown-derivative-profile.json').read_text());raw2=raw(2);raw3=raw(3);records=[]
 def direction(cols,var):return [{i//6:a for i,a in c.items() if i%6==var} for c in cols]
 for r in fold['eight_source_folds']:
  for channel,ys in enumerate(r['pure_known_G2quadratic_Y0_two_channels']):
   w=apply(I2,{i:ex(a) for i,a in ys if i<42})
   for var in (4,5):
    a=extract(corrected(apply(direction(raw2,var),w),Pi17,t['active_F1_free17']),n1,r1)
    records.append({'kind':'G2_first67_coordinate_L1','extra':r['extra_active_source'],'contraction':r['contraction'],'channel':channel,'direction':['r','ell'][var-4],'D_U':list(map(str,a))})
 for r in four['remaining_four_source_directions']:
  w=apply(I3,{i:ex(a) for i,a in r['actual_projected_zero_linear_source238_combination']})
  for var in range(6):
   a=extract(corrected(apply(direction(raw3,var),w),Pi42,t['active_F2_free42']),n2,r2)
   records.append({'kind':'G3_first182_coordinate_L1','extra':r['source_extra_index'],'direction':var,'D_U':list(map(str,a))})
 assert records==old['all56_span_coordinate_positive_L1_bounds'] and b1==old['repr67_L1_height_bound_bits'] and b2==old['repr182_L1_height_bound_bits']
 for var,exp in [(4,[1,0]),(5,[0,1])]:
  for col in t['active11_source_columns_polynomial_r_ell']:apply(Pi17,{row:ex(x) for row,terms in col for m,x in terms if m==exp})
  for j in t['active56_source_columns']:apply(Pi17,{i//6:ex(x) for i,x in raw2[j].items() if i%6==var})
 assert maximum==243==old['all_factor_construction_height_bound_bits']
 classes={};membership=[];committed=maximum
 for kind,key in [('G2','active_F1_fixed67_columns'),('G3','active_F2_fixed182_columns')]:
  values=[Q(x) for c in t[key] for i,x in c];den=math.lcm(*(x.denominator for x in values));up=max(abs(x.numerator)*(den//x.denominator) for x in values);classes[kind]=[den,up]
 for r in records:
  kind=r['kind'][:2];d,u=map(int,r['D_U']);dm,um=classes[kind];height=max((d*dm).bit_length(),(u*um).bit_length());membership.append(height);committed=max(committed,height)
 for kind,b in [('G2',b1),('G3',b2)]:committed=max(committed,b+max(v.bit_length() for v in classes[kind]))
 verify=reports['AMBIENT_B18_ELL2_T2_SOURCE8_RESOLUTION_VERIFY_HEIGHT_20261105_v63'];assert classes=={k:list(map(int,v)) for k,v in verify['fixed_matrix_common_D_U'].items()}
 assert max(membership)==committed==280==verify['all_factor_and_membership_commit_height_bits'] and 3*committed+1==841==verify['uncommitted_scalar_product_sum_height_bits']
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'independent-prior-height-and-source-method-proof.json'
 proof={'activeG2':'P17 last inverse[a11,D2sel56,E17]; free17 correction is a11 multiple, killed byP. Source clearing D2sel killed byP; cachedI2 hasD20I2=0, hence G2=P17D21I2 exactly. Quadrell inverse derivative -Pi17(K1r*aell+K1ell*ar).',
        'tensor':'mixed51 differential to zero1 target is +G1active and to active1 target has z factor, so its pureplane quad sign is+. zero9 gives onlyz terms.',
        'G3':'Use cachedI3 and q-E42Pi42q before fixed182 extraction. Quotient derivative stays in separateeightpsi, not discarded.',
        'bound':'Positive common-denominator/L1 envelopes include all fixed elimination and scalar divisions; separate fixed-matrix Dm/Um covers all membership partial sums.280 commits/841 temporaries fit unchanged384/1153 guards.',
        'scope':'Only source factor construction may proceed within combined<60k inputs/scratch,32MiB/60s/2GB. No signed response or Phi dot authorization.'}
 ep.write_text(json.dumps({'all56_positive_envelopes_independently_replayed':True,'prior_factor_bits':243,'all_membership_commit_bits':280,'uncommitted_bits':841,'method_proof':proof},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'finite_resolution_factor_source_method_and_prior_guard_accepted':True,'new_factor_or_signed_response_values_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'commit_bits':280,'temporary_bits':841}))
if __name__=='__main__':main()

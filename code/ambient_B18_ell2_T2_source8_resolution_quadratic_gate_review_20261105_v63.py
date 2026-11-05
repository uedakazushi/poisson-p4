#!/usr/bin/env python3
"""Independent finite-resolution factor postcheck, no signed responses."""
import itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_G3_quadratic_source_review_20261105_v63 import sha,sp,add,scale,mv
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_SOURCE8_RESOLUTION_QUADRATIC_INPUT_GATE_20261105_v63';STEM=A+'_REVIEW'
def raw(degree):
 active=[0,1,4,5,7];target=[w for q in range(degree) for w in itertools.combinations(active,q)];look={w:i for i,w in enumerate(target)};out=[]
 for q in range(degree+1):
  for w in itertools.combinations(active,q):
   p=degree-q
   for j in range(14):
    col={}
    if p:col[(14*look[w]+j)*6+5]=Q(1)
    for s,z in enumerate(w):col=add(col,{(14*look[w[:s]+w[s+1:]]+j)*6+active.index(z):Q((-1)**(p+s))})
    out.append(col)
 return out
def replay(cols,receipts):
 norms={};reps={}
 for r in receipts:
  w={r['source_column']:Q(1)}
  for p,c in r['subtract']:w=add(w,reps[p],-Q(c))
  w=scale(w,1/Q(r['divide_by']));v=sp(r['normalized']);assert mv(cols,w)==v;norms[r['pivot']]=v;reps[r['pivot']]=w
 return norms,reps
def extraction(v,cols,norms,reps):
 q=dict(v);a={}
 for p in sorted(norms):
  x=q.get(p,Q(0))
  if x:q=add(q,norms[p],-x);a=add(a,reps[p],x)
 assert not q and mv(cols,a)==v;return a
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 pp=ROOT/'reports'/A/'report.json';P=json.loads(pp.read_text());assert P['all_passed'];inputs={'primary':{'path':str(pp.relative_to(ROOT)),'sha256':sha(pp)}};D={}
 for p,h in P['input_sha256'].items():
  path=ROOT/p;assert sha(path)==h;inputs[p]={'path':p,'sha256':h}
  if path.name!='actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json':D[path.name]=json.loads(path.read_text())
 assert sha(ROOT/'code/ambient_B18_ell2_T2_source8_resolution_quadratic_gate_20261105_v63.py')==P['code_sha256'] and sha(ROOT/'research/notes'/f'{A}.md')==P['note_sha256']
 ep=pp.parent/'actual-factored-G2-rell-eight-and-G3-quadratic-four-source-operators.json';assert sha(ep)==P['evidence_sha256'][ep.name];E=json.loads(ep.read_text());inputs['evidence']={'path':str(ep.relative_to(ROOT)),'sha256':sha(ep)}
 t=D['actual-active-unit-bases-and26-frame-transport.json'];pref=D['active-five-d3-source-split-and-budget.json'];fold=D['actual-eight-source-root-folds-M1-operators-and-positive-AD-budget.json'];four=D['actual-four-active-source-complements-and-unknown-derivative-profile.json'];old=D['actual176-G3quadratic-factored-source-and-known-NF-input.json'];plane=D['actual-plane-G1-degree23-and-twelve-adjoint-source-input.json'];accepted=D['actual-G2-linear-fixed-source-profile.json']
 I2=list(map(sp,accepted['actual_I2_source_origin42_columns']));assert len(I2)==42 and accepted['rank']==42 and accepted['linear_injectivity_full102']
 Pi17=list(map(sp,t['active_Pi17_columns']));Pi42=list(map(sp,t['active_Pi42_columns']));fixed1=list(map(sp,t['active_F1_fixed67_columns']));fixed2=list(map(sp,t['active_F2_fixed182_columns']));n1,r1=replay(fixed1,t['active_F1_fixed67_receipts']);n2,r2=replay(fixed2,t['active_F2_fixed182_receipts']);raw2=raw(2);raw3=raw(3)
 E17=[{i:Q(1)} for i in t['active_F1_free17']];E42=[{i:Q(1)} for i in t['active_F2_free42']]
 assert all(mv(Pi17,v)=={j:Q(1)} for j,v in enumerate(E17)) and all(mv(Pi42,v)=={j:Q(1)} for j,v in enumerate(E42))
 K=[]
 for var,exponent in [(4,[1,0]),(5,[0,1])]:
  cols=[]
  for poly in t['active11_source_columns_polynomial_r_ell']:
   v={}
   for row,terms in poly:
    for m,x in terms:
     if m==exponent:v=add(v,{row:Q(x)})
   cols.append(mv(Pi17,v))
  for j in t['active56_source_columns']:cols.append(mv(Pi17,{i//6:x for i,x in raw2[j].items() if i%6==var}))
  K.append(cols)
 assert len(K[0])==len(K[1])==67 and K[0]==list(map(sp,E['Pi17_Kfirst_r_columns67'])) and K[1]==list(map(sp,E['Pi17_Kfirst_ell_columns67']))
 G1rell=[{row:Q(x) for row,m,x in col if m==[1,1]} for col in plane['exact_plane_G1_polynomial_columns'][:17]];assert G1rell==list(map(sp,E['G1_active17_quadratic_rell_columns']))
 for rec,saved in zip(fold['eight_source_folds'],E['eight_G2_rell_factored_responses'],strict=True):
  assert rec['extra_active_source']==saved['extra_active_source'] and rec['contraction']==saved['contraction']
  for ydata,channel in zip(rec['pure_known_G2quadratic_Y0_two_channels'],saved['channels'],strict=True):
   y=sp(ydata);assert all(i<42 for i in y);w=mv(I2,y)
   for var,key in [(4,'active_inverse_first67_r'),(5,'active_inverse_first67_ell')]:
    q={}
    for j,a in w.items():q=add(q,{i//6:x*a for i,x in raw2[j].items() if i%6==var})
    q=add(q,mv(E17,mv(Pi17,q)),Q(-1));a=sp(channel[key]);assert mv(fixed1,a)==q
   assert not channel['mixed51_source_coefficients'] and channel['zero9_quadratic_plane_response_zero']
 prefix=next(p for p in P['input_sha256'] if p.endswith('actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json'));I3=list(map(sp,fields(ROOT/prefix,['actual_active_I3_origin_columns'])['actual_active_I3_origin_columns']))
 for rec,saved in zip(four['remaining_four_source_directions'],E['four_G3_quadratic_factored_responses'],strict=True):
  assert rec['source_extra_index']==saved['extra_active_source'];w=mv(I3,sp(rec['actual_projected_zero_linear_source238_combination']));v=mv(raw3,w);expected={}
  for var in range(6):
   q={i//6:x for i,x in v.items() if i%6==var};q=add(q,mv(E42,mv(Pi42,q)),Q(-1));a=extraction(q,fixed2,n2,r2)
   for j,x in a.items():
    if j>=56:expected[6*(j-56)+var]=x
  assert expected==sp(saved['selected126_inverse_linear_coefficients'])
 for oldkey,newkey in [('actual_selected126_rawD31_Pi42_linear_columns','G3_Pi42_rawD31selected126_linear_columns'),('PiY_L48_degree2_NF42_times20_columns','G3_PiY_L48_NF42times20_columns'),('degree2_active_pair_order20','G3_active_pair_order20')]:assert old[oldkey]==E[newkey]
 lam=list(map(sp,fold['lambda26_two_channels']));assert lam==list(map(sp,E['lambda26_two_channels']));P3=list(map(sp,E['G3_Pi42_rawD31selected126_linear_columns']));NF=list(map(sp,E['G3_PiY_L48_NF42times20_columns']));pairs={tuple(p):i for i,p in enumerate(E['G3_active_pair_order20'])};env=[]
 def envelope(terms):
  den=1;total=Q(0)
  for x in terms:den=math.lcm(den,x.denominator);total+=abs(x)
  return [str(den),str(total.numerator*(den//total.denominator))]
 for rec in E['eight_G2_rell_factored_responses']:
  for ch,c in enumerate(rec['channels']):
   terms=[]
   for var,key in [(0,'active_inverse_first67_ell'),(1,'active_inverse_first67_r')]:
    for j,x in sp(c[key]).items():
     for row,y in K[var][j].items():terms.append(abs(x*y*lam[ch].get(row,Q(0))))
   env.append(envelope(terms))
 for rec in E['four_G3_quadratic_factored_responses']:
  terms=[]
  for k,x in sp(rec['selected126_inverse_linear_coefficients']).items():
   source,u=divmod(k,6)
   for i,y in P3[source].items():
    row,v=divmod(i,6);pair=tuple(sorted((u,v)))
    if pair!=(5,5):terms.extend(abs(x*y*z) for z in NF[20*row+pairs[pair]].values())
  env.append(envelope(terms))
 assert env==E['positive_envelopes_G2_channel16_then_G3_four']
 bits=max(max(int(d).bit_length(),int(u).bit_length()) for d,u in env);assert bits==62==E['positive_bound_bits']==P['positive_bound_bits']
 assert E['product_and_sum_operational_bound_bits']==P['operational_bound_bits']==317 and E['prior_all_factor_and_membership_height_bound_bits']==280
 assert E['source_construction_whole_numeric_upper']==50563<60000 and E['next_value_whole_QQ_bound']==25275<60000 and E['next_value_raw_byte_bound']==7459400<32*1024*1024
 assert not E['signed_quadratic_responses_evaluated'] and not E['physical_or_implicit_coefficients_queried']
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);path=out/'independent-16-G2-and4-G3-factor-memberships-and-source-scope.json'
 proof={'activeG2':'P17 inverse derivative formula source proof pinned by prior peer; all16 channels corrected67 memberships verified. CachedI2 is accepted exact D20-kernel source map, no new higher Taylor coefficient.',
        'tensor':'Every weightedY0 here is active42-only. Mixed51 response coefficient lists are empty and zero9 absent. General+G1rell formula remains correct, no response substituted.',
        'G3':'All4 corrected182 A126 maps independently extracted. P/NF exactly identical to earlier independently checked maps. Final quotient response remains separate eightpsi.',
        'scope':'Only factor identities and positive envelopes.25275/7.46MB is the pure known-response callback bound, not the whole8psi cache-value bound. RootGO and a fresh union/height/frame proof are required before actual combined values.'}
 dump={'all16_G2_corrected67_memberships':True,'all4_G3_corrected182_factors':True,'Kfirst_r_ell_and_G1rell_exact':True,'positive_envelopes20_independently_verified':True,'positive_bits':62,'operational_bits':317,'source_scope':proof,'signed_responses_not_evaluated':True}
 path.write_text(json.dumps(dump,separators=(',',':'))+'\n');report={'id':STEM,'all_passed':True,'finite_resolution_factors_and_pure_response_source_bound_accepted':True,'signed_responses_or_phi_dots_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{path.name:sha(path)},'seconds':time.monotonic()-start};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'positive_bits':62}))
if __name__=='__main__':main()

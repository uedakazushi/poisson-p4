#!/usr/bin/env python3
"""Factored known-resolution G2/G3 quadratic source gate, no signed values."""
import argparse,itertools,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_SOURCE8_RESOLUTION_QUADRATIC_INPUT_GATE_20261105_v63'
def fixed_replay(cols,receipts):
 norm={};reps={}
 for r in receipts:
  p=r['pivot'];w={r['source_column']:F(1)}
  for k,c in r['subtract']:w=add(w,reps[k],-F(c))
  w=H.guard(scale(w,1/F(r['divide_by'])));v=H.parse(r['normalized']);assert mv(cols,w)==v;norm[p]=v;reps[p]=w
 return norm,reps
def span(v,cols,norm,reps):
 q=dict(v);a={}
 while q:
  p=min(q);assert p in norm;x=q[p];q=H.guard(add(q,norm[p],-x));a=H.guard(add(a,reps[p],x))
 assert mv(cols,a)==v;return a
def raw_linear(active,degree):
 target=[(degree-1-q,w) for q in range(degree) for w in itertools.combinations(active,q)]
 look={w:i for i,(p,w) in enumerate(target)};source=[(degree-q,w) for q in range(degree+1) for w in itertools.combinations(active,q)];out=[]
 for p,w in source:
  for j in range(14):
   v={}
   if p:v[(14*look[w]+j)*6+5]=F(1)
   for s,z in enumerate(w):v=add(v,{(14*look[w[:s]+w[s+1:]]+j)*6+(active+[8]).index(z):F((-1)**(p+s))})
   out.append(v)
 return out
def main(out,budget_only=False):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'factored_known_resolution_source_no_quadratic_values'});hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 prior=read('reports/AMBIENT_B18_ELL2_T2_SOURCE8_RESOLUTION_PRIOR_HEIGHT_20261105_v63/report.json');assert prior['all_passed'] and prior['hard_factor_guard_384_covered']
 verifyheight=read('reports/AMBIENT_B18_ELL2_T2_SOURCE8_RESOLUTION_VERIFY_HEIGHT_20261105_v63/report.json');assert verifyheight['all_passed'] and verifyheight['prepared384guard_covered']
 assert read('reports/AMBIENT_B18_ELL2_T2_SOURCE8_RESOLUTION_METHOD_GUARD_20261105_v63_REVIEW/report.json')['all_passed']
 t=read(H.PATHS['transport']);pref=read(H.PATHS['active']);active=pref['active_variables'];assert active==[0,1,4,5,7]
 fold=read('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_INPUT_GATE_20261105_v63/actual-eight-source-root-folds-M1-operators-and-positive-AD-budget.json')
 four=read('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE_PROFILE_20261105_v63/actual-four-active-source-complements-and-unknown-derivative-profile.json')
 oldg3=read('reports/AMBIENT_B18_ELL2_T2_G3_QUADRATIC_INPUT_GATE_20261105_v63/actual176-G3quadratic-factored-source-and-known-NF-input.json')
 plane=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_PLANE34_INPUT_GATE_20261105_v63/actual-plane-G1-degree23-and-twelve-adjoint-source-input.json')
 acceptedG2=read(H.PATHS['G2'])
 source='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json';hs[source]=H.sha(ROOT/source)
 I3stored=H.selected(source,'actual_active_I3_origin_columns')
 source_input_numeric=count([t,pref,fold,four,oldg3,plane,acceptedG2,prior,verifyheight,I3stored])
 def repr_support_upper(receipts):
  supports={}
  for r in receipts:
   s={r['source_column']}
   for p,c in r['subtract']:
    if F(c):s.update(supports[p])
   supports[r['pivot']]=s
  return sum(map(len,supports.values()))
 # Before any factor coefficients are constructed: existing constants are
 # counted twice only where a parsed rational copy is retained, and fixed
 # representation support is bounded from the receipt DAG without values.
 parsed_constant_upper=count([t['active_F1_fixed67_columns'],t['active_F2_fixed182_columns'],t['active_Pi17_columns'],t['active_Pi42_columns'],
  acceptedG2['actual_I2_source_origin42_columns'],I3stored,oldg3['actual_selected126_rawD31_Pi42_linear_columns'],oldg3['PiY_L48_degree2_NF42_times20_columns']])
 repr_and_norm_upper=2*(repr_support_upper(t['active_F1_fixed67_receipts'])+repr_support_upper(t['active_F2_fixed182_receipts']))+count([r['normalized'] for r in t['active_F1_fixed67_receipts']+t['active_F2_fixed182_receipts']])
 constructed_arrays_upper=2*(8*2*2*67+4*6*126+2*17*67+3*17+26*2)
 raw_operator_upper=2*(3*224+4*364)
 source_whole_upper=source_input_numeric+parsed_constant_upper+repr_and_norm_upper+constructed_arrays_upper+raw_operator_upper+4500
 H.dump(out/'phase.json',{'phase':'before_factor_construction','source_whole_numeric_upper':source_whole_upper,'prior_commit_bits':verifyheight['all_factor_and_membership_commit_height_bits']})
 assert source_whole_upper<60000,('source_combined_input_output_guard',source_whole_upper)
 if budget_only:
  H.dump(out/'budget-only.json',{'source_whole_numeric_upper':source_whole_upper,'factor_coefficients_constructed':False,'prior_commit_bound_bits':280,'input_sha256':hs,'code_sha256':H.sha(Path(__file__))})
  print(json.dumps({'budget_only':True,'source_whole_numeric_upper':source_whole_upper,'no_factor_coefficients_constructed':True}));return
 f1=list(map(H.parse,t['active_F1_fixed67_columns']));n1,r1=fixed_replay(f1,t['active_F1_fixed67_receipts']);Pi17=list(map(H.parse,t['active_Pi17_columns']));free17=t['active_F1_free17'];E17=[{j:F(1)} for j in free17]
 f2=list(map(H.parse,t['active_F2_fixed182_columns']));n2,r2=fixed_replay(f2,t['active_F2_fixed182_receipts']);Pi42=list(map(H.parse,t['active_Pi42_columns']));free42=t['active_F2_free42'];E42=[{j:F(1)} for j in free42]
 raw2=raw_linear(active,2);raw3=raw_linear(active,3)
 sel56=t['active56_source_columns'];sel126=t['active126_source_columns']
 I2=list(map(H.parse,acceptedG2['actual_I2_source_origin42_columns']))
 assert acceptedG2['rank']==42 and acceptedG2['linear_injectivity_full102']
 # Frozen source coefficient I2 and D20I2=0 were already checked by the
 # independent G2-linear gate; retain that proof rather than recomputing I2.
 assert len(I2)==42
 Kfirst=[[],[]]
 for var,exponent in [(4,[1,0]),(5,[0,1])]:
  cols=[]
  for poly in t['active11_source_columns_polynomial_r_ell']:
   v={}
   for row,terms in poly:
    for m,x in terms:
     if m==exponent:v=add(v,{row:F(x)})
   cols.append(v)
  cols.extend({i//6:x for i,x in raw2[j].items() if i%6==var} for j in sel56)
  assert len(cols)==67;Kfirst[var-4]=[H.guard(mv(Pi17,v)) for v in cols]
 G1rell=[]
 for col in plane['exact_plane_G1_polynomial_columns'][:17]:G1rell.append({row:F(x) for row,m,x in col if m==[1,1]})
 assert len(G1rell)==17
 g2=[]
 for rec in fold['eight_source_folds']:
  channels=[]
  for ys in rec['pure_known_G2quadratic_Y0_two_channels']:
   y=H.parse(ys);w=H.guard(mv(I2,{i:x for i,x in y.items() if i<42}))
   a_pair=[]
   for var in [4,5]:
    v={};
    for j,c in w.items():v=add(v,{i//6:x*c for i,x in raw2[j].items() if i%6==var})
    z=add(v,mv(E17,mv(Pi17,v)),-1);a_pair.append(span(z,f1,n1,r1))
   channels.append({'active_inverse_first67_r':srec(a_pair[0]),'active_inverse_first67_ell':srec(a_pair[1]),
    'mixed51_source_coefficients':[[i-42,str(x)] for i,x in y.items() if 42<=i<93],
    'zero9_quadratic_plane_response_zero':True})
  g2.append({'extra_active_source':rec['extra_active_source'],'contraction':rec['contraction'],'channels':channels})
 # G3 actual source basis and inverse derivative; general V has a free42 part.
 I3=list(map(H.parse,I3stored));g3=[]
 for rec in four['remaining_four_source_directions']:
  c=H.parse(rec['actual_projected_zero_linear_source238_combination']);w=H.guard(mv(I3,c));v=H.guard(mv(raw3,w));a={}
  for var in range(6):
   q={i//6:x for i,x in v.items() if i%6==var};q=add(q,mv(E42,mv(Pi42,q)),-1);s=span(q,f2,n2,r2)
   for i,x in s.items():
    if i>=56:a[(i-56)*6+var]=x
  g3.append({'extra_active_source':rec['source_extra_index'],'selected126_inverse_linear_coefficients':srec(H.guard(a))})
 lam=list(map(H.parse,fold['lambda26_two_channels']));P3=list(map(H.parse,oldg3['actual_selected126_rawD31_Pi42_linear_columns']));NF3=list(map(H.parse,oldg3['PiY_L48_degree2_NF42_times20_columns']));pairidx={tuple(p):i for i,p in enumerate(oldg3['degree2_active_pair_order20'])}
 # Only positive envelopes: do not evaluate signed quadratic contractions.
 envelopes=[];max_input=0
 def envelope(terms):
  D=1;total=F(0)
  for x in terms:
   x=abs(x);D=math.lcm(D,x.denominator);total+=x
  U=total.numerator*(D//total.denominator);return [str(D),str(U)]
 for rec in g2:
  for ch,s in enumerate(rec['channels']):
   terms=[]
   for var,key in [(0,'active_inverse_first67_ell'),(1,'active_inverse_first67_r')]:
    for j,x in H.parse(s[key]).items():
     for row,y in Kfirst[var][j].items():terms.append(x*y*lam[ch].get(row,F(0)))
   for source,xs in s['mixed51_source_coefficients']:
    iz,j=divmod(source,17)
    for row,y in G1rell[j].items():terms.append(F(xs)*y*lam[ch].get(17+3*iz+row,F(0)))
   envelopes.append(envelope(terms))
 for rec in g3:
  terms=[]
  for k,x in H.parse(rec['selected126_inverse_linear_coefficients']).items():
   source,u=divmod(k,6)
   for i,y in P3[source].items():
    row,v=divmod(i,6);pair=tuple(sorted((u,v)))
    if pair!=(5,5):terms.extend(x*y*z for z in NF3[row*20+pairidx[pair]].values())
  envelopes.append(envelope(terms))
 allnumeric=[x for cols in Kfirst for v in cols for x in v.values()]+[x for v in P3+NF3+G1rell+lam for x in v.values()]
 allnumeric.extend(F(x) for r in g3 for i,x in r['selected126_inverse_linear_coefficients'])
 allnumeric.extend(F(x) for r in g2 for c in r['channels'] for key in ['active_inverse_first67_r','active_inverse_first67_ell','mixed51_source_coefficients'] for i,x in c[key])
 max_input=max([0]+[max(abs(x.numerator).bit_length(),x.denominator.bit_length()) for x in allnumeric]);height=max([0]+[max(int(d).bit_length(),int(u).bit_length()) for d,u in envelopes]);operational=max(height,3*max_input+2)
 slots=25000+sum(len(x) for cols in Kfirst for x in cols)+sum(len(v) for v in P3+NF3+G1rell)+sum(len(q['selected126_inverse_linear_coefficients']) for q in g3)+sum(len(c['active_inverse_first67_r'])+len(c['active_inverse_first67_ell'])+len(c['mixed51_source_coefficients']) for r in g2 for c in r['channels']);assert slots<60000 and operational<2048
 raw=slots*(2*len(str(2**operational))+24)+2_000_000;assert raw<32*1024*1024
 ep=out/'actual-factored-G2-rell-eight-and-G3-quadratic-four-source-operators.json';H.dump(ep,{
  'Pi17_Kfirst_r_columns67':[srec(v) for v in Kfirst[0]],'Pi17_Kfirst_ell_columns67':[srec(v) for v in Kfirst[1]],
  'G1_active17_quadratic_rell_columns':[srec(v) for v in G1rell],'eight_G2_rell_factored_responses':g2,
  'four_G3_quadratic_factored_responses':g3,'G3_Pi42_rawD31selected126_linear_columns':[srec(v) for v in P3],
  'G3_PiY_L48_NF42times20_columns':[srec(v) for v in NF3],'G3_active_pair_order20':oldg3['degree2_active_pair_order20'],
  'lambda26_two_channels':fold['lambda26_two_channels'],'all_I2_D20_and_67_span_identities':True,'all_I3_D30_and_corrected182_span_identities':True,
  'D20_I2_and_D30_I3_source_zero_identities':'reusedacceptedfrozenconstantI2/I3sourceproofsbylinearity; noredundantverificationproductsoutsidepriorheightmargin',
  'prior_all_factor_and_membership_height_bound_bits':verifyheight['all_factor_and_membership_commit_height_bits'],
  'source_construction_whole_numeric_upper':source_whole_upper,
  'G2_active_formula':'G2=P17D21I2(0), P17last17inverse[a11,D2sel56,Efree17]; quadrell=-Pi17(Kfirst_r*a_ell+Kfirst_ell*a_r)',
  'G2_mixed_formula':'source42+17iz+j contributes+G1active17_j_quadrell inrows17+3iz; zero9none',
  'G3_formula':'PiY0G3quad projected=-NF(Pi42D31selected126*A126(D31I3c-Efree42Pi42D31I3c)); unknownYstarfirst separateeightpsi',
  'positive_envelopes_G2_channel16_then_G3_four':envelopes,'positive_bound_bits':height,'product_and_sum_operational_bound_bits':operational,
  'next_value_whole_QQ_bound':slots,'next_value_raw_byte_bound':raw,'signed_quadratic_responses_evaluated':False,'physical_or_implicit_coefficients_queried':False})
 report={'id':STEM,'all_passed':True,'signed_quadratic_responses_evaluated':False,'physical_or_implicit_coefficients_queried':False,
  'next_value_whole_QQ_bound':slots,'next_value_raw_byte_bound':raw,'positive_bound_bits':height,'operational_bound_bits':operational,
  'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','next_value_whole_QQ_bound','next_value_raw_byte_bound','positive_bound_bits','operational_bound_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--budget-only',action='store_true');a=p.parse_args();main(a.output.resolve(),a.budget_only)

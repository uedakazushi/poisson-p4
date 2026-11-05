#!/usr/bin/env python3
"""PRIOR positive bound for finite-resolution factor construction only."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_source8_resolution_quadratic_gate_20261105_v63 import raw_linear
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_SOURCE8_RESOLUTION_PRIOR_HEIGHT_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={};maximum=0;operations=0
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def check(a):
  nonlocal maximum,operations
  operations+=1;maximum=max(maximum,a[0].bit_length(),a[1].bit_length());return a
 def exact(x):x=F(x);return check((x.denominator,abs(x.numerator)))
 def plus(a,b):
  d=math.lcm(a[0],b[0]);return check((d,a[1]*(d//a[0])+b[1]*(d//b[0])))
 def times(a,b):return check((a[0]*b[0],a[1]*b[1]))
 def apply(columns,v):
  q={}
  for j,b in v.items():
   for row,x in columns[j].items():q[row]=plus(q.get(row,(1,0)),times(exact(x),b))
  return q
 def repr_bounds(receipts):
  norm={};reps={};bound=0
  for r in receipts:
   a=(1,1)
   for p,c in r['subtract']:a=plus(a,times(exact(c),reps[p]))
   a=times(a,exact(1/F(r['divide_by'])));p=r['pivot'];reps[p]=a;norm[p]=H.parse(r['normalized']);bound=max(bound,a[0].bit_length(),a[1].bit_length())
  return norm,reps,bound
 def extraction(v,norm,reps):
  q=dict(v);a=(1,0)
  for p in sorted(norm):
   x=q.pop(p,(1,0))
   if not x[1]:continue
   a=plus(a,times(x,reps[p]))
   for row,c in norm[p].items():
    if row!=p:q[row]=plus(q.get(row,(1,0)),times(x,exact(c)))
  return a
 def correct(v,Pi,free):
  p=apply(Pi,v);q=dict(v)
  for k,b in p.items():q[free[k]]=plus(q.get(free[k],(1,0)),b)
  return q
 t=read(H.PATHS['transport']);active=[0,1,4,5,7];Pi17=list(map(H.parse,t['active_Pi17_columns']));Pi42=list(map(H.parse,t['active_Pi42_columns']))
 n1,r1,b1=repr_bounds(t['active_F1_fixed67_receipts']);n2,r2,b2=repr_bounds(t['active_F2_fixed182_receipts'])
 g=read(H.PATHS['G2']);I2=list(map(H.parse,g['actual_I2_source_origin42_columns']))
 prefix='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json';hs[prefix]=H.sha(ROOT/prefix);I3=list(map(H.parse,H.selected(prefix,'actual_active_I3_origin_columns')))
 fold=read('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_INPUT_GATE_20261105_v63/actual-eight-source-root-folds-M1-operators-and-positive-AD-budget.json')
 four=read('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE_PROFILE_20261105_v63/actual-four-active-source-complements-and-unknown-derivative-profile.json')
 raw2=raw_linear(active,2);raw3=raw_linear(active,3)
 def rawvar(raw,var):return [{k//6:x for k,x in c.items() if k%6==var} for c in raw]
 raw2pair=[rawvar(raw2,v) for v in [4,5]];raw3six=[rawvar(raw3,v) for v in range(6)]
 records=[]
 for rec in fold['eight_source_folds']:
  for channel,ys in enumerate(rec['pure_known_G2quadratic_Y0_two_channels']):
   y={i:exact(x) for i,x in ys if i<42};w=apply(I2,y)
   for var in range(2):
    v=apply(raw2pair[var],w);q=correct(v,Pi17,t['active_F1_free17']);a=extraction(q,n1,r1)
    records.append({'kind':'G2_first67_coordinate_L1','extra':rec['extra_active_source'],'contraction':rec['contraction'],'channel':channel,'direction':['r','ell'][var],'D_U':list(map(str,a))})
 for rec in four['remaining_four_source_directions']:
  c={i:exact(x) for i,x in rec['actual_projected_zero_linear_source238_combination']};w=apply(I3,c)
  for var in range(6):
   v=apply(raw3six[var],w);q=correct(v,Pi42,t['active_F2_free42']);a=extraction(q,n2,r2)
   records.append({'kind':'G3_first182_coordinate_L1','extra':rec['source_extra_index'],'direction':var,'D_U':list(map(str,a))})
 # Pi17*Kfirst coefficients are sums/products of literal frozen inputs.
 prior=maximum
 for which,exponent in [(4,[1,0]),(5,[0,1])]:
  for col in t['active11_source_columns_polynomial_r_ell']:
   v={row:exact(x) for row,terms in col for m,x in terms if m==exponent};apply(Pi17,v)
  for j in t['active56_source_columns']:apply(Pi17,{i//6:exact(x) for i,x in raw2[j].items() if i%6==which})
 assert maximum<=2048,('source_prior_envelope_exceeds_declared_ceiling',maximum)
 ep=out/'prior-input-receipt-product-and-partialsum-height-envelopes.json';H.dump(ep,{
  'repr67_L1_height_bound_bits':b1,'repr182_L1_height_bound_bits':b2,'span_construction_bound_before_PiKfirst_bits':prior,
  'all_factor_construction_height_bound_bits':maximum,'uncommitted_scalar_product_sum_bound_bits':3*maximum+1,
  'all56_span_coordinate_positive_L1_bounds':records,
  'method':'Everybound(D,U)ensuresallcoeffdenominatorsdivideD andsumabsoluteD-clearednumerators<=U. ExactfrozenreceiptoperationspropagatecommonLCM/positivepartialsumanddivisionbounds. Coordinatewisepositiveeliminationignoresonlyprovedzeroresidualfreecoordinates, notunknownsourcevalues. No signednewfactororquadraticresponseconstructed.',
  'cached_I2_I3_used_instead_of_recomputation':True,'new_G2_G3_factor_vectors_constructed':False,
  'hard_factor_guard_384_covered':maximum<=384,'no_cap_increase_or_value_authorization':True})
 report={'id':STEM,'all_passed':True,'factor_construction_prior_height_bound_bits':maximum,'hard_factor_guard_384_covered':maximum<=384,
  'scalar_temporary_prior_bound_bits':3*maximum+1,'new_factor_or_quadratic_values_constructed':False,
  'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'positive_operations':operations,'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','factor_construction_prior_height_bound_bits','hard_factor_guard_384_covered','scalar_temporary_prior_bound_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())

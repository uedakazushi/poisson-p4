#!/usr/bin/env python3
"""Formal SOURCE8 old rell sensitivities; no second source values."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_AD_20261105_v63'
def main(out,gate):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'authorized_SOURCE8_formal_AD_no_second_values'})
 gp=ROOT/gate;gr=json.loads(gp.read_text());assert gr['all_passed'] and not gr['AD_started']
 ip=gp.parent/'actual-eight-source-root-folds-M1-operators-and-positive-AD-budget.json';assert H.sha(ip)==gr['evidence_sha256'][ip.name];d=json.loads(ip.read_text())
 atoms={i:list(map(F,v)) for i,v in d['actual_atom_origin_r_ell_firsts']};orig={i:F(v) for i,v in d['typed_constants']};deps={i:set(v) for i,v in d['typed_structural_second_dependency_labels']}
 vals={};grad={};maxbits=0;temporary_peak=0
 def guard(v,commit=True):
  nonlocal maxbits,temporary_peak
  for c in v:
   b=max(abs(c.numerator).bit_length(),c.denominator.bit_length());temporary_peak=max(temporary_peak,b)
   assert b<=2048,('formal_AD_bit_guard',b)
   if commit:maxbits=max(maxbits,b)
  return v
 def plus(a,b):return guard([x+y for x,y in zip(a,b)])
 def times(a,b):
  return guard([a[0]*b[0],a[0]*b[1]+b[0]*a[1],a[0]*b[2]+b[0]*a[2],a[0]*b[3]+b[0]*a[3]+a[1]*b[2]+a[2]*b[1]])
 for i,n in d['needed_typed_nodes']:
  if n[0]=='c':v=[F(n[1],n[2]),F(0),F(0),F(0)];g={}
  elif n[0]=='actual_C0_scalar_function':v=atoms[n[1]];g={n[1]:F(1)}
  elif n[0]=='sum':
   v=[F(0)]*4;g={}
   for j,a,b in n[1]:c=F(a,b);v=plus(v,[c*x for x in vals[j]]);g=add(g,grad[j],c);guard(g.values())
  else:
   assert n[0]=='product';v=[F(n[1],n[2]),F(0),F(0),F(0)];g={}
   for j in n[3]:g=add(scale(g,vals[j][0]),grad[j],v[0]);guard(g.values());v=times(v,vals[j])
  guard(v);guard(g.values());assert v[0]==orig[i] and set(g)<=deps[i];vals[i]=v;grad[i]=g
 results=[]
 for r,b in zip(d['eight_source_folds'],d['positive_output_bound_rows']):
  known=F(0);unknown={}
  for (kind,i),cs in r['old_typed_or_literal_C0_root_functional']:
   c=F(cs)
   if kind:unknown=add(unknown,{i:c})
   else:known+=c*vals[i][3];unknown=add(unknown,grad[i],c)
  guard([known]);guard(unknown.values())
  D,U=map(int,b['old_known_cross_bound']);assert D%known.denominator==0 and abs(known.numerator)*(D//known.denominator)<=U
  D,U=map(int,b['old_gradient_L1_bound']);assert all(D%c.denominator==0 for c in unknown.values())
  assert sum(abs(c.numerator)*(D//c.denominator) for c in unknown.values())<=U
  results.append({'extra_active_source':r['extra_active_source'],'contraction':r['contraction'],
   'known_first_cross_response_only':str(known),'literal_original_C0_second_rell_sensitivity':srec(unknown),
   'actual_second_or_phi_values_not_substituted':True})
 assert maxbits<=d['positive_AD_height_bound_bits']
 slots=sum(4+len(grad[i]) for i in vals)+sum(1+len(r['literal_original_C0_second_rell_sensitivity']) for r in results)+d['next_AD_slot_breakdown']['folded_source_and_M1_factored_records']+14000
 assert slots<=d['next_AD_whole_QQ_bound']<60000
 ep=out/'actual-eight-old-formal-rell-sensitivities-with-M1-native-retained.json';H.dump(ep,{
  'all_eight_formal_old_results':results,'typed_known4_lower_first_coefficients':[[i,list(map(str,v))] for i,v in sorted(vals.items())],
  'typed_literal_second_gradients':[[i,srec(v)] for i,v in sorted(grad.items())],
  'unknown_original_C0_second_label_list':sorted({i for r in results for i,c in r['literal_original_C0_second_rell_sensitivity']}),
  'M1_A0_native_Xp1_G2quadratic_factored_responses_unchanged':{'source_path':str(ip.relative_to(ROOT)),'sha256':H.sha(ip)},
  'all_opaque_seconds_independent':True,'no_actual_second_phi_comparison_or_gauge_values':True,
  'whole_combined_eight_sensitivity_not_yet_closed':True,'actual_T2_HF2_or_length_not_claimed':True})
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report={'id':STEM,'all_passed':True,'formal_old_sensitivity_nonzeros':sum(len(r['literal_original_C0_second_rell_sensitivity']) for r in results),
  'unknown_original_second_label_count':len({i for r in results for i,c in r['literal_original_C0_second_rell_sensitivity']}),
  'actual_second_phi_or_comparison_values_computed':False,'new_physical_or_implicit_coefficients':False,
  'input_sha256':{gate:H.sha(gp),str(ip.relative_to(ROOT)):H.sha(ip)},'code_sha256':H.sha(Path(__file__)),
  'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:H.sha(ep)},
  'owned_QQ_slot_bound':slots,'actual_committed_height_bits':maxbits,'source_height_bound_bits':d['positive_AD_height_bound_bits'],
  'actual_temporary_peak_bits':temporary_peak,'source_temporary_upper_bits':3*d['positive_AD_height_bound_bits']+1,'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','formal_old_sensitivity_nonzeros','unknown_original_second_label_count','owned_QQ_slot_bound','actual_committed_height_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);a=p.parse_args();main(a.output.resolve(),a.gate)

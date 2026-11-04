#!/usr/bin/env python3
"""Exact polynomial inputs on two factoredcurves; monic a division only."""
import hashlib,importlib.util,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_CURVE_COMPONENTS_INPUT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def ser(p):return[[list(m),str(c)]for m,c in sorted(p.items(),reverse=True)]
def main():
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);helper=ROOT/'code/s61_nonsquare_raw_H_B_pair_norm_input_20261104_v63.py';spec=importlib.util.spec_from_file_location('arith',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
 fb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_COMMON_CURVE_FACTOR_20261104_v63';fp=fb/'certificate.json';assert sha(fp)==load(fb/'report.json')['evidence_sha256'][fp.name];factors=load(fp)['factors'];Qf=ar.poly(factors[0]['factor']);Lf=ar.poly(factors[1]['factor']);assert max(m[0]for m in Qf)==2 and max(m[0]for m in Lf)==1
 hb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_MONIC_COMPLETION_20261104_v63';hp=hb/'certificate.json';assert sha(hp)==load(hb/'report.json')['evidence_sha256'][hp.name];H=ar.poly(load(hp)['monic_H_laurent']);Hn=ar.shift(H,(5,0,0));assert all(m[0]>=0 for m in Hn)
 pb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_COPRIME_PAIR_PREFLIGHT_20261104_v63';pp=pb/'certificate.json';assert sha(pp)==load(pb/'report.json')['evidence_sha256'][pp.name];g=ar.poly(load(pp)['selected_pair_source']['polynomials'][1]);kb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_ADDITIONAL_MINOR_NORM_INPUT_20261104_v63';kp=kb/'certificate.json';assert sha(kp)==load(kb/'report.json')['evidence_sha256'][kp.name];k=ar.poly(load(kp)['source_k'])
 def divide(p,F):
  deg=max(m[0]for m in F);assert {m:c for m,c in F.items()if m[0]==deg}=={(deg,0,0):Q(1)};r=dict(p);q={};tr=[]
  while r and max(m[0]for m in r)>=deg:
   n=max(m[0]for m in r);t={(m[0]-deg,m[1],m[2]):c for m,c in r.items()if m[0]==n};q=ar.add(q,t);r=ar.add(r,ar.scale(ar.mul(t,F),-1));tr.append([n,len(t),len(r)])
  assert ar.add(ar.mul(q,F),r)==p
  return q,r,tr
 rec=[]
 for name,F in [('Q',Qf),('L',Lf)]:
  pol=[]
  for label,P in [('a5H',Hn),('g',g),('k',k)]:
   q,r,tr=divide(P,F);pol.append({'source_label':label,'quotient':ser(q),'remainder':ser(r),'division_steps':tr,'remainder_stats':ar.stats(r),'min_d_degree':min([m[2]for m in r],default=-1),'dleading_coefficient':ser({(m[0],m[1],0):c for m,c in r.items()if m[2]==max([m[2]for m in r],default=-1)})})
  rec.append({'factor_name':name,'factor':ser(F),'factor_degree_a':max(m[0]for m in F),'polynomials':pol,'total_remainder_terms':sum(len(x['remainder'])for x in pol),'all_three_d0_terms_zero':all(x['min_d_degree']>=1 for x in pol)})
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);cert={'variable_order':['a','e','d'],'H_source_cleared_by_known_a5':True,'component_records':rec,'all6_literal_monic_a_division_identities_checked':True,'unknown_factor_inverted':False,'field_GCD_subresultant_NOT_computed':True,'generic_field_form_L':'QQ(e)[d] with a=-(9e^2+3e-6)/50','generic_field_form_Q':'QQ(e)[a]/Q is a degree2 field by exactirreducibility+Gauss; coefficientsrepresentedmodmonicQ','denominator_exception_requirement':'Any future d-pseudo sequence LC must be nonzero in componentfunctionfield and its zero divisor retained as finite exceptions; do not invert itglobally.','no_new_resultant_GB_point_prime':True};(out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');paths=[Path(__file__),helper,ROOT/'research/notes'/(ID+'.md'),fb/'report.json',fp,hb/'report.json',hp,pb/'report.json',pp,kb/'report.json',kp];report={'id':ID,'component_profiles':[{'name':x['factor_name'],'terms':[len(p['remainder'])for p in x['polynomials']],'total':x['total_remainder_terms'],'d_min_degrees':[p['min_d_degree']for p in x['polynomials']],'all_d0_terms_zero':x['all_three_d0_terms_zero']}for x in rec],'stage_peak_polynomial_terms':ar.peak,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{'certificate.json':sha(out/'certificate.json')},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Sharper full918 clearing prior from accepted four actual point rows."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_CLEARING_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={}
 def read(p):pins[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):pins[p]=H.sha(ROOT/p);return H.selected(p,k)
 rows='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_ROWS_GUARDED_20261105_v63/';r=read(rows+'report.json');assert r['all_passed'];actual=field(rows+'actual-four-specific-B99-point-dual-rows-STOP-only384.json','actual_four_rows_channel_pair')
 assert read('reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_ROWS_GUARDED_20261105_v63_REVIEW/report.json')['all_passed']
 old='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_INPUT_GATE_20261105_v63/Kfull396-to-full918-alpha-factor-input-and-positive-stage-bounds.json'
 raw={k:field(old,k) for k in ['PiY3_from102_columns','selected25_G3_source_indices','J25_full918_columns','F25_columns27','all418_G3_full918_columns','formal_Kfull_to_alpha_factor_program']}
 Pi=list(map(H.parse,raw['PiY3_from102_columns']));J=list(map(H.parse,raw['J25_full918_columns']));FF=list(map(H.parse,raw['F25_columns27']));G3=list(map(H.parse,raw['all418_G3_full918_columns']));assert [len(Pi),len(J),len(FF),len(G3)]==[102,25,27,418]
 rr={};re={}
 for v in actual:
  assert v['all99_source_and_normalized_memberships'];(rr if v['pair_rr_rell']==0 else re)[v['channel']]=H.parse(v['point_row102'])
 assert set(rr)==set(re)=={0,1};pref=[{**{9*i+7:x for i,x in re[t].items()},**{9*i+8:x for i,x in rr[t].items()}} for t in range(2)]
 sourcebase=count([r,actual,raw,Pi,J,FF,G3,rr,re,pref,pins]);assert sourcebase+20000<150000
 maximum=0
 def mark(q):
  nonlocal maximum
  maximum=max(maximum,q[0].bit_length(),q[1].bit_length());assert maximum<2048;return q
 def ex(x):x=F(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def add(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 zero=(1,0);ljbounds=[];alpha=[]
 for t in range(2):
  lj=[]
  for col in J:
   q=zero
   for i,x in col.items():q=add(q,mul(ex(pref[t].get(i,F(0))),ex(x)))
   lj.append(q)
  row={}
  for i in range(102):
   for var in range(9):
    q=ex(pref[t].get(9*i+var,F(0)))
    for target,pi in Pi[i].items():
     for b,x in FF[9*target+var].items():q=add(q,mul(lj[b],mul(ex(pi),ex(x))))
    if q[1]:row[9*i+var]=q
  ljbounds.append(lj);alpha.append(row)
 clearing=maximum
 for row in alpha:
  for col in G3:
   q=zero
   for i,x in col.items():q=add(q,mul(row.get(i,zero),ex(x)))
 allprior=maximum
 payload={**raw,'actual_four_point_rows':actual,'pre_two_full918_rows':[[[i,str(x)] for i,x in sorted(v.items())] for v in pref],'positive_two_full918_alpha_envelopes':[[[i,str(D),str(U)] for i,(D,U) in sorted(v.items())] for v in alpha],'positive25_LJ_envelopes':[[[str(D),str(U)] for D,U in v] for v in ljbounds],'full_clearing_formula':'alpha_t=pre_t-(pre_tJ25)F25PiY3;slot9*F2source102+var9,ordinarycoefficients','full836_G3_annihilation_checks_required':True,'old424_forecast_preserved':'99f34321b075836ace94919945af0ebc6ea6e9d15707d8706a4f0d2af0187435','actual_alpha_or_G3dot_values_not_evaluated':True}
 n=count(payload);future=2*n+count([Pi,J,FF,G3,pref])+16000;assert future<60000
 ep=out/'actual-four-row-to-full918-alpha-clearing-pruned-input.json';H.dump(ep,payload)
 report={'id':STEM,'all_passed':True,'source_positive_only_no_signed_clearing_or_G3dots':True,'clearing_positive_prior_bits':clearing,'all836_G3_check_positive_prior_bits':allprior,'temporary_prior_bits':3*allprior+1,'ready_with_unchanged384':allprior<=384,'source_owned_bound':sourcebase+20000,'future_value_numeric_bound':future,'pruned_numeric_records':n,'future_RAW_output_cap':2000000,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','clearing_positive_prior_bits','all836_G3_check_positive_prior_bits','temporary_prior_bits','ready_with_unchanged384','source_owned_bound','future_value_numeric_bound','pruned_numeric_records','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

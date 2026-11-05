#!/usr/bin/env python3
"""Prepared TWO raw33 pairings; no V30 pairing or source-jet query."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_ALPHA29_BETA29_20261105_v63'
def main(out,gate,peer,custody):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 gp=ROOT/gate;g=json.loads(gp.read_text());pp=ROOT/peer;p=json.loads(pp.read_text());cp=ROOT/custody;c=json.loads(cp.read_text());assert g['all_passed'] and p['all_passed'] and H.sha(gp) in json.dumps(p);assert c['all_passed'] and c['safe_help_literal_PASS'] and c['worker_sha256']==H.sha(Path(__file__)) and c['source_gate_sha256']==H.sha(gp)
 for path,h in c['worker_and_helper_sha256'].items():assert H.sha(ROOT/path)==h
 for path,h in g['input_sha256'].items():assert H.sha(ROOT/path)==h
 ip=gp.parent/'two-pairing-source-binding-and-positive-height-input.json';assert H.sha(ip)==g['evidence_sha256'][ip.name];d=json.loads(ip.read_text());assert d['no_pairing_or_new_function_values'] and g['new_pairing_count_only']==2;assert count([d,g,p,c])+1000<g['whole_numeric_bound']==5222
 cap=g['source_committed_bit_bound'];tmp=g['source_temporary_bit_bound'];assert cap==76 and tmp==229
 H.dump(out/'phase.json',{'phase':'authorized_TWO_cached_pairings_only','gate_sha256':H.sha(gp),'source_peer_sha256':H.sha(pp),'custody_sha256':H.sha(cp),'no_physical_gauge_or_jet_query':True})
 maximum=0;cursor={}
 def ck(x):
  nonlocal maximum
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length());maximum=max(maximum,b)
  if b>cap:
   H.dump(out/'guard-frontier.json',{'cursor':cursor,'completed_new_pairings':list(map(str,a29)),'completed_recovered30':list(map(str,a30)),'offending_uncommitted_value':str(x),'bits':b,'no_retry_or_cap_raise':True});raise RuntimeError('source76guard')
  return x
 V=[{i:F(x) for i,x in v} for v in d['primitive_V29_raw63_mask1']];rho=[[{i:F(x) for i,x in v} for v in row] for row in d['two_universal_bar_rho33']];aW=list(map(F,d['known_alpha_beta_on_Wv']));a23=list(map(F,d['known_alpha23_beta23']));a29=[];a30=[];pieces=[]
 for ch,row in enumerate(rho):
  cursor={'stage':'two_new_pairings','channel':ch};q=F(0);component=[]
  for t,rr in enumerate(row):
   z=F(0)
   for i,x in rr.items():cursor.update({'component':t,'basis33':i,'partial_component':str(z),'partial_pairing':str(q)});z=ck(z+ck(x*V[t].get(i,F(0))))
   component.append(str(z));q=ck(q+z)
  D,U=map(int,d['two_new_dot_positive_bounds'][ch]);assert D%q.denominator==0 and abs(q.numerator)*(D//q.denominator)<=U;a29.append(q);pieces.append(component)
  cursor={'stage':'recover30_by_proved_vector_linearity','channel':ch};z=ck(F(2,25)*ck(ck(aW[ch]+ck(F(3375,11)*a23[ch]))+ck(F(50)*q)));D,U=map(int,d['two_recovered30_positive_bounds'][ch]);assert D%z.denominator==0 and abs(z.numerator)*(D//z.denominator)<=U;a30.append(z)
 ep=out/'actual-two29-pairings-and-recovered30-kernel-invariants.json';H.dump(ep,{'actual_alpha29_beta29':list(map(str,a29)),'actual_alpha30_beta30_recovered_without_new_pairing':list(map(str,a30)),'known_alpha23_beta23':list(map(str,a23)),'known_alphaW_betaW':list(map(str,aW)),'two_pairing_component_sums':pieces,'actual_all_three_point_kernel_direction_invariants':{'23':list(map(str,a23)),'29':list(map(str,a29)),'30':list(map(str,a30))},'Wv_recovery_follows_verified99_source_vector_identity':True,'only_TWO_new_pairings_evaluated':True,'no_new_Phi_physical_gauge_implicit_or_jet_coefficients':True,'scope':'Onlygauge-invariantM1contractionsatfixedpointkerneldirections. Fullupper619typedclearing,b1/A2/Yp0/M1carriesandactualquadraticimagesremainseparate;noHF2ortotalclaim.'})
 assert sum(v.stat().st_size for v in out.iterdir())<g['whole_RAW_bound']
 report={'id':STEM,'all_passed':True,'new_pairing_count':2,'actual_alpha29_beta29':list(map(str,a29)),'actual_alpha30_beta30_recovered':list(map(str,a30)),'actual_upper619_values_or_HF2_not_computed':True,'no_new_physical_gauge_or_jet_coefficients':True,'input_sha256':{gate:H.sha(gp),peer:H.sha(pp),custody:H.sha(cp),str(ip.relative_to(ROOT)):H.sha(ip)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:H.sha(ep)},'source_committed_bits':cap,'temporary_bits':tmp,'actual_maximum_bits':maximum,'whole_numeric_bound':g['whole_numeric_bound'],'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','actual_alpha29_beta29','actual_alpha30_beta30_recovered','actual_maximum_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output.resolve(),a.gate,a.source_peer,a.custody)

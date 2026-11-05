#!/usr/bin/env python3
"""Known coefficient source/prior for D=Cfull-rhoM0; no D or RHS values."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER617_D_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={}
 def read(p):pins[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):pins[p]=H.sha(ROOT/p);return H.selected(p,k)
 cp='reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_20261105_v63/actual224-Cfull-coefficients-with-native-pivot-SUM-and-Psi.json';C=field(cp,'Cfull_two_by2_by56');cc=[[[F(x) for x in row] for row in ch] for ch in C]
 rp='reports/AMBIENT_B18_ELL2_T2_UPPER617_RHO_20261105_v63/actual-two-rho1144-covectors-and-full918-identities.json';rr=field(rp,'rho_two_by1144_sparse_rows');rho=list(map(H.parse,rr))
 qp='reports/AMBIENT_B18_ELL2_T2_UPPER617_RHO_INPUT_GATE_20261105_v63/pruned-transpose-R-unit28-source-and-positive-rho-input.json';pairs=field(qp,'quadratic_pair_order44');pindex={tuple(v):i for i,v in enumerate(pairs)}
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';rawM=field(mug,'Mnew0_columns');M=list(map(H.parse,rawM));assert len(M)==56
 active_pairs=sorted({(7,7),(7,8)}|{tuple(pairs[i%44]) for row in rho for i in row});assert active_pairs==[(0,7),(0,8),(1,7),(1,8),(4,7),(4,8),(5,7),(5,8),(7,7),(7,8)]
 for p in ['reports/AMBIENT_B18_ELL2_T2_UPPER617_RHO_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/report.json']:assert read(p)['all_passed']
 maximum=0
 def mark(q):
  nonlocal maximum
  maximum=max(maximum,q[0].bit_length(),q[1].bit_length());assert maximum<=384;return q
 def ex(x):x=F(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def add(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 env=[]
 for t in range(2):
  ch=[]
  for pair in active_pairs:
   row=[];mi=pindex[pair]
   for j,col in enumerate(M):
    q=ex(cc[t][1][j]) if pair==(7,7) else ex(cc[t][0][j]) if pair==(7,8) else (1,0)
    for target,x in col.items():q=add(q,mul(ex(rho[t].get(44*target+mi,F(0))),ex(x)))
    row.append([str(q[0]),str(q[1])])
   ch.append(row)
  env.append(ch)
 payload={'Cfull_two_by2_by56':C,'rho_two_by1144_sparse_rows':rr,'Mnew0_columns56':rawM,'quadratic_pair_order44':pairs,'active_pairs10':[list(v) for v in active_pairs],'positive_D_two_by10_by56':env,'formula':'D_t,pair,row56=Cfull_ell[row]ifrr+Cfull_r[row]ifrell-sum_targetrho_t[44target+pair44]*Mnew0[row][target].','full_upper_source_formula':'upper_t=Gamma_t(Xc3-Xp2*a1)_(r²ell)+D_t(Xc2-Xp1*a1);all617lower-kernelrecipes.','no_OLDpivot_honest2_cache_assumption':True,'native_individualXp2_full_L2S0_retained':True,'all_actual_D_or_RHS_functionvalues_not_computed':True,'scope':'Onlycoefficientfunctionalinput/prior.SourceOLD/nativecallbacks,617weights,honestADanduppervaluesremainseparate.'}
 n=count(payload);future=2*n+count([cc,rho,M,pindex])+8000;assert future<30000
 ep=out/'pruned-D-coefficient-functional-and-ten-pair-source-input.json';H.dump(ep,payload)
 sourceowned=count([C,cc,rr,rho,pairs,pindex,rawM,M,env,payload,pins])+10000;assert sourceowned<60000
 report={'id':STEM,'all_passed':True,'positive_source_only_no_signed_D':True,'coefficient_count':1120,'source_positive_prior_bits':maximum,'temporary_prior_bits':3*maximum+1,'pruned_numeric_records':n,'future_value_numeric_bound':future,'future_RAW_output_cap':1000000,'source_owned_bound':sourceowned,'no_RHS_Z2_or_upper_response_evaluated':True,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','source_positive_prior_bits','temporary_prior_bits','pruned_numeric_records','future_value_numeric_bound','source_owned_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

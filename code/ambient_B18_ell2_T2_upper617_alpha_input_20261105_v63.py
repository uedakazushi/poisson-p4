#!/usr/bin/env python3
"""Positive/source-only Kfull-to-alpha918 gate. No signed dual/alpha."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={}
 def read(p):pins[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):pins[p]=H.sha(ROOT/p);return H.selected(p,k)
 kp='reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_20261105_v63/actual396-Kfull-P5-and-native23-existing-cache-values.json';K=field(kp,'K99_two_channels_rr_rell');assert len(K)==99 and all(len(v)==2 and all(len(q)==2 for q in v) for v in K)
 assert read('reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_20261105_v63_REVIEW/report.json')['all_passed']
 mu='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json';normraw=field(mu,'normalized_basis_rows');repraw=field(mu,'normalized_basis_original_selected_source_representations');Yraw=field(mu,'actual_selected_Y0_columns102');norm={i:H.parse(v) for i,v in normraw};rep={i:H.parse(v) for i,v in repraw};assert len(norm)==len(rep)==99 and set(norm)==set(rep)
 assert set(range(102))-set(norm)=={1,36,38}
 for i,v in norm.items():assert v.get(i)==1 and all(j>=i for j in v)
 pre='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json';Pi=list(map(H.parse,field(pre,'PiY3_from102_columns')));G3raw=field(pre,'actual_full_G3linear102_by418_columns');G3=list(map(H.parse,G3raw));assert len(Pi)==102 and len(G3)==418
 r25='reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json';indices=field(r25,'selected_columns');assert len(indices)==25;J=[G3[i] for i in indices]
 fixed='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_FIXED_ROWS_20261105_v63/actual-fixed-specific-dual-rows-and-F25-with-next-clearing-prior.json';Fraw=field(fixed,'actual_F25_columns27');FF=list(map(H.parse,Fraw));assert len(FF)==27
 for p in ['reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_FIXED_ROWS_20261105_v63/report.json','reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER2_ADJOINT_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_INPUT_GATE_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_20261105_v63_REVIEW/report.json']:assert read(p)['all_passed']
 sourcecount=count([K,normraw,repraw,Yraw,norm,rep,Pi,G3raw,G3,J,FF,Fraw,indices,pins]);assert sourcecount+30000<150000
 maximum=0
 def mark(v):
  nonlocal maximum
  maximum=max(maximum,v[0].bit_length(),v[1].bit_length());assert maximum<2048,'positiveforecastonly:stopwithoutanysignedalpha';return v
 def ex(x):x=F(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def add(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 zero=(1,0);stage={};dual=[];rhsall=[]
 # No signed K*representation or back substitution takes place. These are
 # absolute envelopes, using source coefficient indices and cached magnitudes.
 for t in range(2):
  ch=[];rhs_ch=[]
  for pair in range(2):
   rhs={}
   for p,v in rep.items():
    q=zero
    for j,x in v.items():q=add(q,mul(ex(K[j][t][pair]),ex(x)))
    rhs[p]=q
   row={}
   for p in reversed(sorted(norm)):
    q=rhs[p]
    for j,x in norm[p].items():
     if j!=p:q=add(q,mul(ex(x),row.get(j,zero)))
    if q[1]:row[p]=q
   ch.append(row);rhs_ch.append(rhs)
  dual.append(ch);rhsall.append(rhs_ch)
 stage['four_specific_B99_dual_rows']=maximum;raw_prior=maximum
 # Include all source membership products, not just the row-construction rhs.
 for t in range(2):
  for pair in range(2):
   for col in list(map(H.parse,Yraw))+list(norm.values()):
    q=zero
    for i,x in col.items():q=add(q,mul(dual[t][pair].get(i,zero),ex(x)))
 stage['four_rows_and_all99_source_and_normalized_memberships']=maximum;row_membership_prior=maximum
 pref=[]
 for t in range(2):
  q={9*i+8:x for i,x in dual[t][0].items()};q.update({9*i+7:x for i,x in dual[t][1].items()});pref.append(q)
 alpha=[];ljs=[]
 for t in range(2):
  lj=[]
  for col in J:
   q=zero
   for i,x in col.items():q=add(q,mul(pref[t].get(i,zero),ex(x)))
   lj.append(q)
  row={}
  for i in range(102):
   for var in range(9):
    q=pref[t].get(9*i+var,zero)
    for target,pi in Pi[i].items():
     for b,x in FF[9*target+var].items():q=add(q,mul(lj[b],mul(ex(pi),ex(x))))
    if q[1]:row[9*i+var]=q
  alpha.append(row);ljs.append(lj)
 stage['two_full918_fixed25_clearing_covectors']=maximum;full_prior=maximum
 # Full418 equations are a distinct known-array verification after values;
 # their operation bounds are nevertheless included before any permission.
 for row in alpha:
  for col in G3:
   q=zero
   for i,x in col.items():q=add(q,mul(row.get(i,zero),ex(x)))
 stage['all836_fullG3linear_annihilation_checks']=maximum;all_prior=maximum
 formal={'point_B99_row_definition':'s_t,pair*Ystar_j0=K[j,t,pair],s_t,pair[E1,E36,E38]=0. Obtainrhs_p=sum_j K[j,t,pair]*repr[p,j];backsubdescendingpivotnorm99.','raw_full918_orientation':'pre_t[9*i+7]=s_t,rell[i];pre_t[9*i+8]=s_t,rr[i];allotherpreinputs0.','fixed25_clearing':'alpha_t=pre_t-(pre_t*J25)*F25*PiY3;PiY3actsseparatelyoneachvar9. ALL918 slotsretained.','truefunctional':'alpha_Kfull(Ybar1)=Kfull_rr*b1_ell+Kfull_rell*b1_r,whereb1=B99*(I-J25F25PiY3)Ybar1.','full_K_not_pureU':'AllKinputsincludeCfull/Psi. Sourcefull393zero+rank190full918strandexactgivesactualkernelannihilation;workerMUSTstillverifyall836literalG3dots.','no_arbitrary_full_X2_or_raw33_descent':'Adjointoutputsupport/rightsolveandactualsourceZ2cachegateSTILLSEPARATE.'}
 def boundrec(v):return [[i,str(D),str(U)] for i,(D,U) in sorted(v.items())]
 payload={'K99_two_channels_rr_rell':K,'selected99_Y0_columns102':Yraw,'normalized99_basis':normraw,'normalized99_selected_representations':repraw,'free_Y3_point_axes':[1,36,38],'PiY3_from102_columns':field(pre,'PiY3_from102_columns'),'selected25_G3_source_indices':indices,'J25_full918_columns':[G3raw[i] for i in indices],'F25_columns27':Fraw,'all418_G3_full918_columns':G3raw,'formal_Kfull_to_alpha_factor_program':formal,'positive_four_dual_rows':[[boundrec(v) for v in ch] for ch in dual],'positive_two_full918_alpha_rows':list(map(boundrec,alpha)),'positive_stage_maxima':stage,'signed_dual_alpha_or_G3dot_values_computed':False}
 n=count(payload);row_future=2*n+4*102+count([norm,rep,Pi,J,FF])+5000;full_future=2*n+2*918*2+4*102+count([norm,rep,Pi,G3,J,FF])+8000
 assert row_future<150000 and full_future<150000
 ep=out/'Kfull396-to-full918-alpha-factor-input-and-positive-stage-bounds.json';H.dump(ep,payload)
 report={'id':STEM,'all_passed':True,'positive_source_only_no_signed_alpha':True,'raw_row_construction_prior_bits':raw_prior,'row_and_source_membership_prior_bits':row_membership_prior,'full918_clearing_prior_bits':full_prior,'all_G3_annihilation_prior_bits':all_prior,'row_stage_ready_with_unchanged384':row_membership_prior<=384,'full_stage_ready_with_unchanged384':all_prior<=384,'no_cap_raise_or_signed_operation_if_notready':True,'source_numeric_bound':sourcecount+30000,'future_row_numeric_bound':row_future,'future_full_numeric_bound':full_future,'pruned_numeric_records':n,'future_RAW_upper':4000000,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','raw_row_construction_prior_bits','row_and_source_membership_prior_bits','full918_clearing_prior_bits','all_G3_annihilation_prior_bits','row_stage_ready_with_unchanged384','full_stage_ready_with_unchanged384','future_row_numeric_bound','future_full_numeric_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

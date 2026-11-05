#!/usr/bin/env python3
"""Pruned existing-first inputs and positive prior only; no signed M1 fold."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
from ambient_B18_ell2_raw33_quotient_20261105_v63 import srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_NONACTIVE9_M1_FIRST_INPUT_GATE_V2_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 support=read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_UPPER_SUPPORT_PREFLIGHT_20261105_v63/eighteen-formal-weights-and-coefficient-free-source-ancestor-inventory.json');used=support['selected99_source_indices_used']
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';mu='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json'
 allc=field(mug,'all741_zero_origin_source_combinations');selected=field(mu,'actual_selected_Y0_source_columns');C={j:H.parse(allc[selected[j]]['original768']) for j in used};del allc
 allY=field(mu,'actual_selected_Y0_columns102');Y={j:H.parse(allY[j]) for j in used};del allY
 mini=read('reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json');piv=mini['original27_X_source_columns'];free=mini['free30_V_indices'];inv=list(map(H.parse,mini['fixed_bottom30_inverse_columns']))
 S=set().union(*(set(v) for v in C.values()));R={};X0={};peakread=0
 for offset,path,key in [(0,H.PATHS['old640'],'source_columns640_rows56_sparse_dual9'),(640,H.PATHS['old128'],'source_columns128_rows48_sparse_dual9')]:
  raw=field(path,key);peakread=max(peakread,count(raw)+count([support,C,Y,mini])+8000);assert peakread<60000
  for k,col in enumerate(raw):
   j=k+offset
   if j not in S and j not in piv:continue
   v={row:H.parse(a) for row,a in col}
   if j in S:R[j]={row:a[8] for row,a in v.items() if a.get(8)}
   if j in piv:X0[j]={row:a[0] for row,a in v.items() if a.get(0)}
  del raw
 olde=field('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json','actual_first640_ell_columns');E={j:H.parse(olde[j]) for j in S if j<640};del olde
 L=list(map(H.parse,field(mug,'L0_Res4_Ntranspose_columns')));Seall=field(mug,'native_S1ell_columns');Se={j-640:H.parse(Seall[j-640]) for j in S if j>=640};del Seall
 S0all=field(mug,'native_S0_columns');S0=list(map(H.parse,S0all));del S0all
 M=list(map(H.parse,field(mug,'Mnew0_columns')));G2=list(map(H.parse,field(mug,'full_G2_linear234_by102_columns')))
 lambdas=list(map(H.parse,field('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_INPUT_GATE_20261105_v63/actual-eight-source-root-folds-M1-operators-and-positive-AD-budget.json','lambda26_two_channels')))
 graph=read('reports/AMBIENT_B18_ELL2_T2_LINEAR_GRAPH4_SOURCE_20261105_v63/actual23-linear-graph-four-common-contractions-and-six-adjoint-input.json');W=H.parse(graph['actual_H0v_source56']);assert graph['two_gauge_invariants']=='alpha=mu[0,r], beta=mu[0,ell]-mu[1,r]' and W[29] and W[23]
 p=min(lambdas[0]);s=29;assert lambdas[0][p]
 initial=count([support,C,Y,R,X0,E,L,Se,S0,M,G2,inv,free,piv,lambdas,W])+8000;assert initial<60000
 H.dump(out/'phase.json',{'phase':'positive_bound_only_before_signed_first_or_M1_fold','owned_pruned_inputs_and_scratch':initial,'maximum_input_read_ownership':peakread})
 maximum=0
 def ck(b):
  nonlocal maximum
  maximum=max(maximum,b[0].bit_length(),b[1].bit_length());return b
 def ex(x):x=F(x);return ck((x.denominator,abs(x.numerator)))
 def mul(a,b):return ck((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  d=math.lcm(a[0],b[0]);return ck((d,a[1]*(d//a[0])+b[1]*(d//b[0])))
 def acc(v,k,b):
  if b[1]:v[k]=plus(v.get(k,(1,0)),b)
 def vector(v):return {k:ex(x) for k,x in v.items() if x}
 def padd(v,w,c=(1,1)):
  q=dict(v)
  for k,b in w.items():acc(q,k,mul(c,b))
  return q
 def apply(A,v):
  q={}
  for j,b in v.items():
   for k,x in A[j].items():acc(q,k,mul(ex(x),b))
  return q
 EB={j:vector(v) for j,v in E.items()}
 for j,v in Se.items():EB[640+j]=apply(L,vector(v))
 Z={};AB={}
 for j,c in C.items():
  apply(S0,{k-640:ex(x) for k,x in c.items() if k>=640})
  z=[];ab=[]
  for raw in [R,EB]:
   q={}
   for k,x in c.items():q=padd(q,vector(raw[k]) if raw is R else raw[k],ex(x))
   a=apply(inv,{i:q[row] for i,row in enumerate(free) if row in q});a={i:b for i,b in a.items() if i<27};ab.append(a);z.append(padd(q,apply([X0[k] for k in piv],a)))
  Z[j]=z;AB[j]=ab
  # Prior for every later source-first chain membership, not just the M1 fold.
  for var,z in zip([7,8],z):
   apply(M,z);rhs=[{i//9:x for i,x in v.items() if i%9==var} for v in G2];apply(rhs,vector(Y[j]))
 weights=[];weighted=[];cofbounds=[];opbounds=[]
 for r in support['eighteen_formal_contractions']:
  ww=H.parse(r['formal_weights_on99psi2']);weights.append({'original418_source':r['original418_source'],'contraction':r['contraction'],'weights':r['formal_weights_on99psi2']});v=[[{},{}],[{},{}]]
  for k,x in ww.items():
   j,ch=divmod(k,2)
   for direction in range(2):v[ch][direction]=padd(v[ch][direction],Z[j][direction],ex(x))
  weighted.append(v);op={}
  for ch,lam in enumerate(lambdas):
   for target,x in lam.items():
    for direction,sourcev in [(0,v[ch][1]),(1,v[ch][0])]:
     for source,b in sourcev.items():acc(op,(direction,target,source),mul(ex(x),b))
  # Fixed extraction of alpha_W,beta_W,alpha23,beta23; these are coefficient
  # bounds only. The eventual full2912-slot test decides representation.
  ellps=op.get((1,p,s),(1,0));ellp23=op.get((1,p,23),(1,0));lp=lambdas[0][p];l1=lambdas[1].get(p,F(0));bW=mul(ellps,ex(1/(lp*W[s])));be=plus(mul(ellp23,ex(1/lp)),mul(bW,ex(W[23])))
  ar_s=plus(op.get((0,p,s),(1,0)),mul(ex(l1),mul(bW,ex(W[s]))));aW=mul(ar_s,ex(1/(lp*W[s])));ar_23=plus(op.get((0,p,23),(1,0)),mul(ex(l1),plus(mul(bW,ex(W[23])),be)));ae=plus(mul(ar_23,ex(1/lp)),mul(aW,ex(W[23])))
  cofbounds.append([aW,bW,ae,be]);opbounds.append(op)
  # Include every fixed-basis verification product/partial sum before values.
  for t,x in lambdas[0].items():
   for source,w in W.items():
    acc(op,(0,t,source),mul(aW,mul(ex(x),ex(w))))
    acc(op,(1,t,source),mul(bW,mul(ex(x),ex(w))))
   acc(op,(0,t,23),mul(ae,ex(x)))
   acc(op,(1,t,23),mul(be,ex(x)))
  for t,x in lambdas[1].items():
   for source,w in W.items():acc(op,(0,t,source),mul(bW,mul(ex(x),ex(w))))
   acc(op,(0,t,23),mul(be,ex(x)))
  assert initial+count([EB,Z,AB,weighted,cofbounds,opbounds])+5000<60000
 # Pruned source input; broad phase dictionaries are not needed by the worker.
 data={'used_selected99_indices':used,'selected99_original768':[[j,srec(C[j])] for j in used],'selected99_Y0':[[j,srec(Y[j])] for j in used],
  'original_r_first_columns':[[j,srec(v)] for j,v in sorted(R.items())],'original640_ell_first_columns':[[j,srec(v)] for j,v in sorted(E.items())],
  'point_Xp27_columns':[srec(X0[j]) for j in piv],'point_Xp27_original_indices':piv,'free30_V_indices':free,'fixed_B30_inverse_columns':[srec(v) for v in inv],
  'native_L0_columns':[srec(v) for v in L],'native_Se_columns':[[j,srec(v)] for j,v in sorted(Se.items())],'native_S0_columns':[srec(v) for v in S0],
  'M0_columns':[srec(v) for v in M],'G2linear_columns':[srec(v) for v in G2],'lambda_two':[srec(v) for v in lambdas],'common_Wv56':srec(W),'fixed_extraction_lambda_row':p,'fixed_extraction_source_row':s,'formal_eighteen_weights':weights,
  'positive_weighted_vector_supports':[[[sorted(q) for q in vv] for vv in v] for v in weighted],
  'positive_four_coefficient_bounds':[[list(map(str,b)) for b in v] for v in cofbounds],
  'no_signed_weighted_first_or_M1_operator_or_invariant_value_evaluated':True,'four_basis_definitions':['lambda0 Mr Wv','lambda0 Mell Wv-lambda1 Mr Wv','lambda0 Mr e23','lambda0 Mell e23-lambda1 Mr e23']}
 rawinput=count(data);zrecords=sum(len(q) for v in Z.values() for q in v);wrecords=sum(len(q) for v in weighted for vv in v for q in vv);maxop=max(map(len,opbounds));future=2*(rawinput+5)+4*zrecords+4*wrecords+4*maxop+8000;source=initial+count([EB,Z,AB,weighted,cofbounds,opbounds])+count(data)+5000
 rawbound=2_000_000+future*(2*math.ceil(maximum*math.log10(2))+40);ready=maximum<=384 and future<60000 and source<60000 and rawbound<32*1024*1024
 data.update({'source_positive_prior_bits':maximum,'future_whole_numeric_bound':future,'future_RAW_bound':rawbound,'source_positive_profile_ready':ready,'source_ownership_bound':source})
 ep=out/'pruned-known-firsts-and-positive-M1-four-invariant-test-input.json';H.dump(ep,data)
 report={'id':STEM,'all_passed':ready,'signed_first_M1_or_invariant_values_evaluated':False,'positive_prior_bits':maximum,'future_temporary_prior_bits':3*maximum+1,'future_whole_numeric_bound':future,'future_RAW_bound':rawbound,'source_ownership_bound':source,'input_sha256':hs,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','positive_prior_bits','future_whole_numeric_bound','future_RAW_bound','source_ownership_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())

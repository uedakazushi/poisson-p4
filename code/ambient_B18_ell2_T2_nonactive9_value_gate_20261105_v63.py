#!/usr/bin/env python3
"""Whole source/cache and positive-height gate; no actual psi values."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
from ambient_B18_ell2_raw33_quotient_20261105_v63 import srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_NONACTIVE9_VALUE_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 for stem in ['NONACTIVE9_LINEAR_AD','NONACTIVE9_M1_FIRST','NONACTIVE9_LOWER_VALUES']:
  assert read('reports/AMBIENT_B18_ELL2_T2_'+stem+'_20261105_v63_REVIEW/report.json')['all_passed']
 cancel=read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_G2_MIXED_ZERO_SOURCE_20261105_v63/report.json');assert cancel['all_passed']
 source='reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_M1_FIRST_INPUT_GATE_V2_20261105_v63/pruned-known-firsts-and-positive-M1-four-invariant-test-input.json';first=read(source)
 ad=read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_LINEAR_AD_20261105_v63/actual-eighteen-formal-old-second-sensitivities-and-native-source-combinations.json')['all18_formal_source_results']
 proof=read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_M1_FIRST_20261105_v63/actual-known-first-weighted-M1-operators-and-four-invariant-tests.json');assert all(r['represented_by_four_known_invariant_functionals'] and all(F(v)==0 for v in r['proposed_four_coefficients']) for r in proof['all18_results'])
 oldpath='reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_20261105_v63/actual104-source-rell-values-and-two-by619-old-knownparts.json';old=read(oldpath);CV={i:F(v[3]) for i,v,grad in old['all104_original_C0_values_and_opaque_gradients']};needed=sorted({i for r in ad for i,x in r['exact_independent_original_C0_second_sensitivity']});assert set(needed)<=set(CV)
 assert read('reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_20261105_v63_REVIEW/report.json')['all_passed']
 base='reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63/whole-cache-eight-psi-four-quadratic-classes-source-input.json'
 keys=['native_L0','native_Lr','native_Le','native_S0','native_Sr','native_Se','native_S2rell','pivot27_source_indices','pivot27_r_columns','pivot27_ell_recipes','Gamma_two','lambda_two','Q_two6_columns','Q_fixed_minor_rows','Q_fixed_minor_determinant']
 cached={k:field(base,k) for k in keys};assert cached['lambda_two']==first['lambda_two'];assert cached['native_L0']==first['native_L0_columns'];assert cached['pivot27_source_indices']==first['point_Xp27_original_indices']
 L0,Lr,Le,S0,Sr,Se,S2=[list(map(H.parse,cached[k])) for k in keys[:7]];Gamma=list(map(H.parse,cached['Gamma_two']));Q=list(map(H.parse,cached['Q_two6_columns']));rows=cached['Q_fixed_minor_rows'];det=Q[0].get(rows[0],F(0))*Q[1].get(rows[1],F(0))-Q[0].get(rows[1],F(0))*Q[1].get(rows[0],F(0));assert det==F(cached['Q_fixed_minor_determinant']) and det
 C={j:H.parse(v) for j,v in first['selected99_original768']};R={j:H.parse(v) for j,v in first['original_r_first_columns']};E={j:H.parse(v) for j,v in first['original640_ell_first_columns']};inv=list(map(H.parse,first['fixed_B30_inverse_columns']));free=first['free30_V_indices'];pivotR=list(map(H.parse,cached['pivot27_r_columns']))
 fixed=count([first,ad,proof,old,cached,C,R,E,inv,Gamma,Q])+6000;assert fixed<60000
 H.dump(out/'phase.json',{'phase':'whole_positive_only_before_any_psi_native_or_carry_value','fixed_numeric_and_scratch':fixed,'all_source_values_existing':True})
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
 def vector(v):return {i:ex(x) for i,x in v.items() if x}
 def padd(a,b):
  q=dict(a)
  for i,x in b.items():acc(q,i,x)
  return q
 def apply(A,v):
  q={}
  for j,b in v.items():
   for i,x in A[j].items():acc(q,i,mul(ex(x),b))
  return q
 def dot(row,v):
  q=(1,0)
  for i,x in row.items():q=plus(q,mul(ex(x),v.get(i,(1,0))))
  return q
 # The individual native pivot ell derivative is a SUM of two matrix
 # products. Neither term is discarded or bounded by a max-union.
 epiv=[]
 for mode,v in cached['pivot27_ell_recipes']:
  if mode=='literal':epiv.append(vector(H.parse(v)))
  else:assert mode=='native';epiv.append(padd(apply(L0,vector(Se[v])),apply(Le,vector(S0[v]))))
 # Existing sourcefirsts for34 zero-origin columns; effectiveL0Se is used
 # only here, where their separate S0-zero source certificate applies.
 epos={j:vector(v) for j,v in E.items()}
 for j in {k for c in C.values() for k in c if k>=640}:epos[j]=apply(L0,vector(Se[j-640]))
 A={}
 for j,c in C.items():
  apply(S0,{k-640:ex(x) for k,x in c.items() if k>=640});pair=[]
  for direction in [0,1]:
   q={}
   for k,x in c.items():
    for row,b in (vector(R[k]) if direction==0 else epos[k]).items():acc(q,row,mul(ex(x),b))
   aa=apply(inv,{i:q[row] for i,row in enumerate(free) if row in q});pair.append({i:b for i,b in aa.items() if i<27})
  A[j]=pair
 bounds=[];psi_bounds=[];weighted_a_support=0;maximum_native_support=0
 for n,(r,wrec) in enumerate(zip(ad,first['formal_eighteen_weights'])):
  assert (r['original418_source'],r['contraction'])==(wrec['original418_source'],wrec['contraction']);ww=H.parse(wrec['weights']);aa=[[{},{}],[{},{}]]
  for k,x in ww.items():
   j,ch=divmod(k,2)
   for direction in range(2):
    for i,b in A[j][direction].items():acc(aa[ch][direction],i,mul(ex(x),b))
  weighted_a_support+=sum(len(q) for ch in aa for q in ch);oldb=(1,0)
  for i,x in r['exact_independent_original_C0_second_sensitivity']:oldb=plus(oldb,mul(ex(x),ex(CV[i])))
  nbs=[];pbs=[]
  for ch,cnrec in enumerate(r['exact_original_native128_combinations_two_channels']):
   cn=vector(H.parse(cnrec));apply(S0,cn)
   native=padd(padd(apply(L0,apply(S2,cn)),apply(Lr,apply(Se,cn))),apply(Le,apply(Sr,cn)))
   maximum_native_support=max(maximum_native_support,len(native));nbs.append(dot(Gamma[ch],native))
   ar,ae=aa[ch];pivot=apply(pivotR,ae)
   for i,b in ar.items():
    for row,x in epiv[i].items():acc(pivot,row,mul(b,x))
   pbs.append(dot(Gamma[ch],pivot))
  total=oldb
  for b in nbs+pbs:total=plus(total,b)
  psi_bounds.append(total);bounds.append({'source418':r['original418_source'],'contraction':r['contraction'],'old':list(map(str,oldb)),'native_channels':[list(map(str,b)) for b in nbs],'pivot_channels':[list(map(str,b)) for b in pbs],'psi':list(map(str,total))})
  assert fixed+count([A,epiv,epos,bounds,aa,native])+5000<60000
 # The complete actual class is -Q*(psi_A,psi_B); rawG3/M1/G2/bottom are
 # separately proved zero. Include every six-coordinate output bound.
 for j in range(9):
  for row in range(6):plus(mul(ex(Q[0].get(row,F(0))),psi_bounds[2*j]),mul(ex(Q[1].get(row,F(0))),psi_bounds[2*j+1]))
 nonrank=maximum;D=math.lcm(*(b[0] for b in psi_bounds));U=max(b[1]*(D//b[0]) for b in psi_bounds);integer_bound=max(D.bit_length(),3*U.bit_length()+8);guard=max(nonrank,integer_bound);temporary=max(3*nonrank+1,integer_bound)
 carrykeys=['selected99_original768','original_r_first_columns','original640_ell_first_columns','free30_V_indices','fixed_B30_inverse_columns','formal_eighteen_weights']
 data={'old18_formal_results':ad,'literal_old104_second_values':[[i,str(CV[i])] for i in needed],'carry_first_inputs':{k:first[k] for k in carrykeys},**cached,
  'all18_component_positive_bounds':bounds,'nonrank_positive_bound_bits':nonrank,'whole_value_height_bound_bits':guard,'temporary_value_height_bound_bits':temporary,
  'coordinate_common_D':str(D),'coordinate_abs_clear_U':str(U),'integer_rank_and_factor_bound_bits':integer_bound,
  'actual_M1_G2_rawG3_and_bottom_responses_identically_zero':True,'G2_zero_source_path':'reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_G2_MIXED_ZERO_SOURCE_20261105_v63/report.json','value_job_started':False,
  'exact_value_formula':'psi=oldliteral104+Gamma(L0S2c+LrSec+LeSrc)-Gamma(Xp_r*aell+Xp_ell*ar); nativeS0zeroeligibleONLYweightedsourcecombinations; individualnativepivotell=L0Se+LeS0. actualclass=-Q*(psiA,psiB).',
  'source_scope':'Complete9nonactivecandidatesonly;162othernonactivesourcezero andactive222closed. Remaining619mixedcyclefamilyOPEN; actualHF2equality/totalnotclaimed.'}
 anew=sum(len(q) for pair in A.values() for q in pair);future=2*(count(data)+4)+6*anew+6*weighted_a_support+4*maximum_native_support+10000;owned=fixed+count([A,epiv,epos,bounds])+count(data)+5000
 raw=2*len(json.dumps(data).encode())+3000*(2*math.ceil(guard*math.log10(2))+40)+3_000_000;ready=guard<=2048 and temporary<=2048 and future<60000 and owned<60000 and raw<32*1024*1024
 data.update({'whole_value_numeric_bound':future,'whole_RAW_bound':raw,'source_builder_owned_bound':owned,'all_source_and_resource_premises_ready':ready});ep=out/'whole-old-native-pivot-nine-classes-positive-source-input.json';H.dump(ep,data)
 report={'id':STEM,'all_passed':ready,'value_job_started':False,'whole_value_numeric_bound':future,'whole_RAW_bound':raw,'whole_value_height_bound_bits':guard,'nonrank_positive_bound_bits':nonrank,'temporary_value_height_bound_bits':temporary,'source_builder_owned_bound':owned,'input_sha256':hs,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','whole_value_numeric_bound','whole_RAW_bound','whole_value_height_bound_bits','nonrank_positive_bound_bits','temporary_value_height_bound_bits','source_builder_owned_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())

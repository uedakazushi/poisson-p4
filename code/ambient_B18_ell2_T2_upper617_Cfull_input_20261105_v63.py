#!/usr/bin/env python3
"""Prior-only source gate for the224 known Cfull coefficients; no signed C."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER617_CFULL_INPUT_GATE_20261105_v63'

def main(out):
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);pins={};maximum=0
 def field(p,k):pins[p]=H.sha(ROOT/p);return H.selected(p,k)
 def read(p):pins[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 cache='reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63/whole-cache-eight-psi-four-quadratic-classes-source-input.json'
 raw={k:field(cache,k) for k in ['Gamma_two','pivot27_source_indices','pivot27_r_columns','pivot27_ell_recipes','native_L0','native_Le','native_S0','native_Se']}
 chart='reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json'
 free=field(chart,'free30_V_indices');B=list(map(H.parse,field(chart,'fixed_bottom30_inverse_columns')))
 assert len(free)==len(B)==30 and len(set(free))==30
 P=[{} for _ in range(56)];Pi=[{} for _ in range(56)]
 for v,col in zip(free,B):
  P[v]={i:x for i,x in col.items() if i<27}
  Pi[v]={i-27:x for i,x in col.items() if i>=27}
 pre='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json'
 assert Pi==list(map(H.parse,field(pre,'PiBottom3_from56_columns')))
 corr='reports/AMBIENT_B18_ELL2_T2_UPPER617_FULL_COORDINATE_SOURCE_20261105_v63/'
 cr=read(corr+'report.json');assert cr['all_passed'] and cr['actual_K_or_upper_values_queried'] is False
 delta=field(corr+'known-Psi-point-correction-and-full-b1-source-annihilation.json','DeltaC_known_point_rows')
 invp='reports/AMBIENT_B18_ELL2_T2_ALPHA29_BETA29_20261105_v63/actual-two29-pairings-and-recovered30-kernel-invariants.json'
 inv=field(invp,'actual_all_three_point_kernel_direction_invariants');alpha=[F(inv[str(j)][0]) for j in (23,29,30)];beta=[F(inv[str(j)][1]) for j in (23,29,30)]
 assert read('reports/AMBIENT_B18_ELL2_T2_ALPHA29_BETA29_20261105_v63_REVIEW/report.json')['all_passed']
 assert read('reports/AMBIENT_B18_ELL2_T2_UPPER_FREE3_GAUGE_SOURCE_20261105_v63_REVIEW/report.json')['all_passed']
 # These are the only retained numerical source fields. The20k scratch reserve
 # covers parsed duplicates, per-pivot positive vectors and224 output envelopes.
 source_count=count([raw,free,B,P,Pi,cr,delta,inv,alpha,beta,pins])
 assert source_count+20000<60000
 def mark(b):
  nonlocal maximum
  maximum=max(maximum,b[0].bit_length(),b[1].bit_length());assert maximum<=384;return b
 def exact(x):x=F(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def add(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def boundcol(col):return {i:exact(x) for i,x in H.parse(col).items()}
 def pluscols(a,b):
  r=dict(a)
  for i,x in b.items():r[i]=add(r.get(i,(1,0)),x)
  return r
 def mv(cols,x):
  y={}
  for j,c in x.items():
   for i,v in cols[j].items():y[i]=add(y.get(i,(1,0)),mul(exact(c),exact(v)))
  return y
 def dot(row,col):
  a=(1,0)
  for i,c in row.items():a=add(a,mul(exact(c),col.get(i,(1,0))))
  return a
 G=list(map(H.parse,raw['Gamma_two']));L0=list(map(H.parse,raw['native_L0']));Le=list(map(H.parse,raw['native_Le']));S0=list(map(H.parse,raw['native_S0']));Se=list(map(H.parse,raw['native_Se']))
 assert len(G)==2 and len(raw['pivot27_r_columns'])==len(raw['pivot27_ell_recipes'])==27
 assert len(L0)==len(Le)==48 and len(S0)==len(Se)==128
 Xr=[boundcol(v) for v in raw['pivot27_r_columns']];Xe=[]
 for kind,arg in raw['pivot27_ell_recipes']:
  if kind=='literal':Xe.append(boundcol(arg))
  else:
   assert kind=='native'
   # Both terms are an actual SUM, never a max/union of alternatives.
   Xe.append(pluscols(mv(L0,Se[arg]),mv(Le,S0[arg])))
 gamma_bound=[[[dot(g,v) for v in Xr],[dot(g,v) for v in Xe]] for g in G]
 delta_rows=[(H.parse(v['r_row56']),H.parse(v['ell_row56'])) for v in delta]
 bounds=[]
 for t in range(2):
  ch=[]
  for k in range(2):
   vals=[]
   for j in range(56):
    q=exact(delta_rows[t][k].get(j,F(0)))
    if t==0:
     weights=alpha if k==0 else beta
     for i,c in Pi[j].items():q=add(q,mul(exact(weights[i]),exact(c)))
    for p,c in P[j].items():q=add(q,mul(gamma_bound[t][k][p],exact(c)))
    vals.append([str(q[0]),str(q[1])])
   ch.append(vals)
  bounds.append(ch)
 pruned={'Gamma_two':raw['Gamma_two'],'pivot27_source_indices':raw['pivot27_source_indices'],'pivot27_r_columns':raw['pivot27_r_columns'],'pivot27_ell_recipes':raw['pivot27_ell_recipes'],'native_L0':raw['native_L0'],'native_Le':raw['native_Le'],'native_S0':raw['native_S0'],'native_Se':raw['native_Se'],'P27_columns56':[[[i,str(x)] for i,x in sorted(c.items())] for c in P],'Pi3_columns56':[[[i,str(x)] for i,x in sorted(c.items())] for c in Pi],'DeltaC_known_point_rows':delta,'alpha_three':list(map(str,alpha)),'beta_three':list(map(str,beta)),'positive_Cfull_two_by2_by56':bounds,'source_formula':'Cfull_r=alpha*Pi3(forchannel0)-GammaXp_rP27+DeltaC_r;Cfull_ell=beta*Pi3(forchannel0)-GammaXp_ellP27+DeltaC_ell;channel1M1term0ONLYinthepeeredexistenceframe.','native_pivot_ell_is_SUM_L0Se_LeS0':True,'full_coordinate_not_pureU':True,'no_signed_C_or_K_or_upper_image_computed':True}
 numeric=count(pruned);future=2*numeric+count([G,L0,Le,S0,Se,P,Pi,Xr,Xe,gamma_bound])+4000
 assert future<60000
 ep=out/'pruned-known-Cfull224-input-and-positive-source-envelopes.json';H.dump(ep,pruned)
 report={'id':STEM,'all_passed':True,'known_source_prior_only_no_signed_C_K_or_relation_values':True,'source_positive_prior_bits':maximum,'source_builder_prior_owned_bound':source_count+20000,'future_Cfull_numeric_bound':future,'pruned_numeric_records':numeric,'native_pivot_ell_SUM_semantics_checked':True,'full_bottom_point_correction_retained':True,'actual_C_K_or_upper_values_not_computed':True,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','source_positive_prior_bits','pruned_numeric_records','future_Cfull_numeric_bound','seconds']}))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

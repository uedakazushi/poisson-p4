#!/usr/bin/env python3
"""Pruned known-input/positive-prior gate, no new covectors or values."""
import json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_G3_quadratic_source_review_20261105_v63 import sha,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_PRIOR_DIAGNOSTIC_20261105_v63'
def count(x):
 if isinstance(x,int):return 1
 if isinstance(x,str):
  try:Q(x);return 1
  except(ValueError,ZeroDivisionError):return 0
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 return 0
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={};maximum=0;stage="input";stage_maxima={}
 def select(p,keys):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys)
 def read(p):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return json.loads((ROOT/p).read_text())
 assert read('reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_SOURCE_PREFLIGHT_20261105_v63/report.json')['all_passed']
 mini=read('reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json');C=mini['free30_V_indices'];bottomInv=list(map(sp,mini['fixed_bottom30_inverse_columns']))
 Pi3=[{k-27:v for k,v in col.items() if k>=27} if i in C else {} for i in range(56) for col in [bottomInv[C.index(i)] if i in C else {}]]
 mup='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json';mu=select(mup,['actual_selected_Y0_source_columns','actual_selected_Y0_columns102','normalized_basis_rows','normalized_basis_original_selected_source_representations'])
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';m=select(mug,['all741_zero_origin_source_combinations','native_S0_columns','native_S1ell_columns','L0_Res4_Ntranspose_columns','full_G2_linear234_by102_columns'])
 selected=mu['actual_selected_Y0_source_columns'];combs=[m['all741_zero_origin_source_combinations'][i]['original768'] for i in selected]
 XR=[]
 for path,key in [('reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json','source_columns640_rows56_sparse_dual9'),('reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json','source_columns128_rows48_sparse_dual9')]:
  F=select(path,[key])
  for column in F[key]:XR.append([[row,str(v)] for row,p in column if (v:=sp(p).get(8,Q(0)))])
 XE=select('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json',['actual_first640_ell_columns'])['actual_first640_ell_columns']
 gpath='reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json';g=read(gpath);Jindices=g['selected_columns'];assert len(Jindices)==25
 p='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json';F=select(p,['PiY3_from102_columns','actual_full_G3linear102_by418_columns']);J=[F['actual_full_G3linear102_by418_columns'][i] for i in Jindices]
 d={'Pi3_source56_columns':[[[k,str(v)] for k,v in sorted(col.items())] for col in Pi3],
 'selected99_original768_combinations':combs,'selected99_Y0_columns102':mu['actual_selected_Y0_columns102'],'B99_normalized_basis99':mu['normalized_basis_rows'],'B99_normalized_selected99_representations':mu['normalized_basis_original_selected_source_representations'],
 'original768_r_first_columns':XR,'original640_ell_first_columns':XE,'native_L0':m['L0_Res4_Ntranspose_columns'],'native_Se':m['native_S1ell_columns'],'native_S0':m['native_S0_columns'],
 'PiY0_source102_columns':F['PiY3_from102_columns'],'J25_full918_columns':J,'J25_projected27_columns':[g['actual_projected_G3linear_columns'][i] for i in Jindices],
 'F25_normalized_basis':g['normalized_basis'],'F25_fixed_receipts':g['fixed_rank_receipts'],'G2linear_source102_columns':m['full_G2_linear234_by102_columns'],
 'lower_duals2_rows3':[[[1,'1'],[2,'4']],[[0,'11'],[2,'270']]],'target_order4':['A_rr','A_rell','B_rr','B_rell'],
 'source_theorem':'Pi3Z2-Hstar1B99(I-J25F25PiY0)Ybar1;G2linYbar1=M0Z2; noM1/G2quadY0. Retainzslots andall25clearings.',
 'no_actual_covectors_adjoint_solve_or_X2_values_yet':True}
 input_numeric=count(d);assert input_numeric<50000
 def mark(a):
  nonlocal maximum
  maximum=max(maximum,a[0].bit_length(),a[1].bit_length());stage_maxima[stage]=max(stage_maxima.get(stage,0),a[0].bit_length(),a[1].bit_length());return a
 def exact(x):x=Q(x);return mark((x.denominator,abs(x.numerator)))
 def plus(a,b):
  if not a[1]:return mark(b)
  if not b[1]:return mark(a)
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def times(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def accum(q,k,b):
  if b[1]:q[k]=plus(q.get(k,(1,0)),b)
 def apply(cols,v):
  out={}
  for j,b in v.items():
   for i,a in cols[j].items():accum(out,i,times(exact(a),b))
  return out
 def positive_extract(norm,reps,n):
  out=[]
  for i in range(n):
   q={i:(1,1)};a={}
   for p in sorted(norm):
    x=q.pop(p,(1,0))
    if not x[1]:continue
    for j,c in norm[p].items():
     if j!=p:accum(q,j,times(x,exact(c)))
    for j,c in reps[p].items():accum(a,j,times(x,c))
   out.append(a)
  return out
 stage='B99_inverse'
 norms={p:sp(v) for p,v in mu['normalized_basis_rows']};reps={p:{j:exact(c) for j,c in v} for p,v in mu['normalized_basis_original_selected_source_representations']};B=positive_extract(norms,reps,102)
 stage='F25_inverse'
 gn={p:sp(v) for p,v in g['normalized_basis']};gr={}
 for r in g['fixed_rank_receipts']:
  a={r['basis_column']:(1,1)}
  for p,c in r['subtract']:
   for j,v in gr[p].items():accum(a,j,times(exact(c),v))
  gr[r['pivot']]={j:times(v,exact(1/Q(r['divide_by']))) for j,v in a.items()}
 FF=positive_extract(gn,gr,27);R=list(map(sp,XR));Eold=list(map(sp,XE));L0=list(map(sp,m['L0_Res4_Ntranspose_columns']));Se=list(map(sp,m['native_S1ell_columns']));H=[[],[]]
 stage='Hstar_firsts'
 for cdata in combs:
  c={j:exact(v) for j,v in cdata};r=apply(Pi3,apply(R,c));e=apply(Eold,{j:b for j,b in c.items() if j<640});ne=apply(L0,apply(Se,{j-640:b for j,b in c.items() if j>=640}))
  for row,b in ne.items():accum(e,row,b)
  H[0].append(r);H[1].append(apply(Pi3,e))
 stage='raw_covectors'
 Lraw=[]
 for ell in d['lower_duals2_rows3']:
  sv=[]
  for var in range(2):
   h={j:(1,0) for j in range(99)}
   for j in range(99):
    for row,c in ell:accum(h,j,times(exact(c),H[var][j].get(row,(1,0))))
   s={}
   for i,column in enumerate(B):
    for j,b in column.items():accum(s,i,times(h[j],b))
   sv.append(s)
  rr={9*i+7:b for i,b in sv[0].items()};re={9*i+7:b for i,b in sv[1].items()}
  for i,b in sv[0].items():accum(re,9*i+8,b)
  Lraw.extend([rr,re])
 stage='cleared_covectors'
 PiY=list(map(sp,d['PiY0_source102_columns']));JJ=list(map(sp,J));Alpha=[]
 for raw in Lraw:
  w=[]
  for col in JJ:
   b=(1,0)
   for i,c in col.items():b=plus(b,times(exact(c),raw.get(i,(1,0))))
   w.append(b)
  projected={}
  for i,col in enumerate(FF):
   for j,b in col.items():accum(projected,i,times(w[j],b))
  a=dict(raw)
  for row,pcol in enumerate(PiY):
   for target,c in pcol.items():
    for var in range(9):accum(a,9*row+var,times(exact(c),projected.get(9*target+var,(1,0))))
  Alpha.append(a)
 # Prior bounds for the future fixed identity verifications too.
 stage='B99_membership'
 for y in map(sp,mu['actual_selected_Y0_columns102']):
  for out in range(99):
   b=(1,0)
   for i,c in y.items():b=plus(b,times(exact(c),B[i].get(out,(1,0))))
 stage='F25_membership'
 for y in map(sp,d['J25_projected27_columns']):
  for out in range(25):
   b=(1,0)
   for i,c in y.items():b=plus(b,times(exact(c),FF[i].get(out,(1,0))))
 stage='C_J_zero_check'
 for a in Alpha:
  for col in JJ:
   b=(1,0)
   for i,c in col.items():b=plus(b,times(exact(c),a.get(i,(1,0))))
 positive_owned=count([B,FF,H,Lraw,Alpha])+input_numeric+12000;assert positive_owned<150000
 future=input_numeric+2*(sum(map(len,B+FF))+sum(len(v) for hv in H for v in hv)+4*918+52*204)+20000
 # Sparse52matrix has at most3times the plane coefficient support, not dense.
 plane_support=sum(sum(i%9 in(7,8) for i,c in col) for col in d['G2linear_source102_columns']);future-=2*52*204;future+=2*3*plane_support
 assert future<150000
 output_coefficients=4*918+sum(len(v) for hv in H for v in hv)+3*plane_support+600
 raw_bound=output_coefficients*(2*len(str(2**maximum))+40)+len(json.dumps(d).encode())+2_000_000
 ready=maximum<=384 and raw_bound<32*1024*1024
 d.update({'prior_stage_maxima':stage_maxima,'positive_prior_bound_bits':maximum,'proposed_commit_guard_bits':384,'proposed_scalar_temporary_bound_bits':3*maximum+1,
 'prior_guard_covered':ready,'prior_positive_owned_bound':positive_owned,'next_known_source_numeric_bound':future,'next_known_source_RAW_bound':raw_bound,
 'Alpha4_positive_sparse_support_upper':[len(a) for a in Alpha],'B99_positive_support_upper':sum(map(len,B)),'F25_positive_support_upper':sum(map(len,FF)),
 'next_52_adjoint_row_domain':'918linearY1slots;52planeoutputs. Mustfirstcheck714zRHSslots0, no solveauthorizedbythisinputgate.',
 'no_signed_covectors_or_physical_values_evaluated':True})
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'pruned-known619-lower4-input-and-prior-resource-gate.json';ep.write_text(json.dumps(d,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'source_only':True,'prior_guard_covered':ready,'positive_prior_bits':maximum,'next_known_source_numeric_bound':future,'next_known_source_RAW_bound':raw_bound,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'no_signed_covectors_or_values_evaluated':True,'seconds':time.monotonic()-start};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'stage_maxima':stage_maxima,'ready_under384':ready,'bits':maximum,'numeric':future,'RAW':raw_bound,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
